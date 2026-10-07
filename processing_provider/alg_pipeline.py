# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""One-click hydrology pipeline (F7) with its run report (F8), v0.17.

Chains the QEHT tools as child algorithms - the pipeline adds no hydrology
of its own, so every step is the same code, with the same tests, as the
individual tool. Writes one output folder:

    layers/   crossings, catchments, flow paths, streams, candidates, QA layers (GeoPackage)
    rasters/  filled DEM, flow direction, accumulation, streams, Strahler, erosion/
    tables/   catchment_characteristics.csv + one CSV per layer
    report/   run_report.html
    package/  design_hydrology.gpkg (only when asked)
    settings.json   every parameter, for "Re-run from settings"
"""

import json
import os
import shutil
import time

from qgis.core import (
    QgsProcessing, QgsProcessingException, QgsProcessingMultiStepFeedback,
    QgsProcessingParameterRasterLayer, QgsProcessingParameterFeatureSource,
    QgsProcessingParameterNumber, QgsProcessingParameterBoolean, QgsProcessingParameterEnum,
    QgsProcessingParameterString, QgsProcessingParameterFile,
    QgsProcessingParameterFolderDestination, QgsProcessingOutputFile,
    QgsProcessingOutputString, QgsProcessingContext, QgsProject,
)

from .base import QehtAlgorithm
from ..core.raster import read_dem, audit_nodata, audit_resampling

FLAT_LABELS = ["toward lower terrain", "Barnes 2014"]
STAGES = ["Full run", "Stop after crossing candidates (review them, then run again with "
          "the reviewed layer as Crossings)"]
PASS_PREFIXES = ("SOIL", "LANDCOVER", "CN_", "RC_CSV", "RAIN_")


class HydrologyPipelineAlgorithm(QehtAlgorithm):

    def name(self): return "hydrologypipeline"
    def displayName(self): return "Run hydrology pipeline (one click)"
    def group(self): return "Workflow"
    def groupId(self): return "workflow"

    def shortHelpString(self):
        return (
            "Runs the whole QEHT chain in one go and writes one output folder with a "
            "fixed layout: <b>layers/</b> (crossings, catchments, flow paths, streams, "
            "candidates and the QA layers), <b>rasters/</b> (filled DEM, flow direction, "
            "accumulation, streams, Strahler order, erosion), <b>tables/</b> (the catchment "
            "characteristics table - one row per crossing, ready for a spreadsheet - and one "
            "CSV per layer), <b>report/run_report.html</b> and <b>settings.json</b>.\n\n"
            "<b>Steps:</b> DEM checks (NoData, nearest-neighbour resampling) → fill → D8 "
            "(Barnes by default) → accumulation → streams → with a road: crossing candidates "
            "(optionally burnt through the embankment, then routed again) → erosion "
            "indices / RUSLE → crossings, catchments and flow paths with soils, curve number, "
            "rainfall and erosion blocks → drainage coverage check → flat-method check → "
            "tables and report. Each step is the stand-alone tool of the same name.\n\n"
            "<b>Crossings:</b> with a road and no crossing layer, the recommended candidates "
            "are used. To review them first, choose 'Stop after crossing candidates': edit "
            "<i>layers/crossing_candidates.gpkg</i> (set status to accepted / rejected, move "
            "or add points), then run again with that layer as <i>Crossings</i> (quickest: "
            "Re-run from settings and set Crossings). Hand-placed pour points work too; they "
            "are snapped to the nearest stream within 5 cells.\n\n"
            "<b>Re-run from settings:</b> give a settings.json from an earlier run and every "
            "parameter is taken from it (a new output folder, and Crossings if you set it, "
            "override).\n\n"
            "<b>Design hydrology package:</b> off by default. Tick it to keep the "
            "self-describing GeoPackage (package/design_hydrology.gpkg) used by design "
            "software such as HEAS; the stand-alone outputs are the same either way.\n\n"
            "All values come from the individual tools; see their help for methods and "
            "limits. Curve numbers from WorldCover are a PROXY, rainfall-derived R is an "
            "ESTIMATE, and both are labelled so in every output.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile(
            "SETTINGS", "Re-run from settings (settings.json of an earlier run; optional)",
            optional=True, fileFilter="JSON (*.json)"))
        self.addParameter(QgsProcessingParameterRasterLayer("DEM", "DEM (raw, projected, metres)",
                                                            optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "ROAD", "Road alignment (centreline; optional)",
            [QgsProcessing.SourceType.TypeVectorLine], optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "CROSSINGS", "Crossings: reviewed candidates or pour points (optional with a road)",
            [QgsProcessing.SourceType.TypeVectorPoint], optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            "START", "Start chainage (m)", QgsProcessingParameterNumber.Type.Double, defaultValue=0.0))
        self.addParameter(QgsProcessingParameterBoolean(
            "REVERSE", "Reverse chainage direction", defaultValue=False))
        self.addParameter(QgsProcessingParameterEnum("STAGE", "Run", options=STAGES, defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber(
            "STREAM_KM2", "Stream threshold (km2)", QgsProcessingParameterNumber.Type.Double,
            defaultValue=1.0, minValue=0.0001))
        self.addParameter(QgsProcessingParameterEnum(
            "FLAT_METHOD", "Flat resolution", options=["Toward lower terrain",
                                                       "Barnes 2014 (recommended)"], defaultValue=1))
        self.addParameter(QgsProcessingParameterBoolean(
            "BURN", "Burn crossings through embankments, then route again (needs a road)",
            defaultValue=False))
        self.addParameter(QgsProcessingParameterNumber(
            "NODATA_OVERRIDE", "Treat this DEM value as NoData (blank = trust the file header)",
            QgsProcessingParameterNumber.Type.Double, optional=True))
        self.addParameter(QgsProcessingParameterString(
            "DEM_SOURCE", "DEM source for the record (e.g. FABDEM 30 m, bilinear to UTM 37S)",
            optional=True))
        self.addParameter(QgsProcessingParameterString(
            "RUN_NAME", "Run name (report title)", optional=True))
        self.add_soil_parameters("SOIL", optional=True)
        self.add_runoff_parameters()
        self.add_rainfall_parameters()
        self.addParameter(QgsProcessingParameterBoolean(
            "EROSION", "Erosion indices and RUSLE (LS-only when R, K or C is missing)",
            defaultValue=True))
        self.addParameter(QgsProcessingParameterRasterLayer(
            "R_RASTER", "R erosivity raster (e.g. GloREDa; optional)", optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            "R_VALUE", "R single value (used when no raster)",
            QgsProcessingParameterNumber.Type.Double, optional=True, minValue=0.0))
        self.addParameter(QgsProcessingParameterBoolean(
            "COVERAGE", "Drainage coverage check along the road", defaultValue=True))
        self.addParameter(QgsProcessingParameterBoolean(
            "PACKAGE", "Keep the design hydrology package (GeoPackage for design software)",
            defaultValue=False))
        self.addParameter(QgsProcessingParameterBoolean(
            "LOAD", "Load the main layers when finished", defaultValue=True))
        self.addParameter(QgsProcessingParameterFolderDestination("OUTPUT", "Output folder"))
        self.addOutput(QgsProcessingOutputFile("REPORT", "Run report"))
        self.addOutput(QgsProcessingOutputFile("CHARACTERISTICS", "Catchment characteristics CSV"))
        self.addOutput(QgsProcessingOutputFile("PACKAGE_FILE", "Design hydrology package"))
        self.addOutput(QgsProcessingOutputString("SUMMARY", "Summary"))

    # -- helpers ---------------------------------------------------------------
    def _settings(self, parameters, context, feedback):
        path = self.parameterAsFile(parameters, "SETTINGS", context) if parameters.get("SETTINGS") else ""
        if not path:
            return dict(parameters)
        try:
            with open(path, encoding="utf-8") as f:
                saved = json.load(f)
        except (OSError, ValueError) as e:
            raise QgsProcessingException(f"Could not read the settings file: {e}")
        p = dict(saved.get("parameters", saved))
        for key in ("OUTPUT", "CROSSINGS", "STAGE"):
            if parameters.get(key) not in (None, "") and key in parameters:
                p[key] = parameters[key]
        p["SETTINGS"] = ""
        feedback.pushInfo(f"Parameters from {path} (saved by QEHT {saved.get('qeht_version', '?')}).")
        return {k: v for k, v in p.items() if self.parameterDefinition(k) is not None}

    def _run(self, alg, params, context, feedback):
        import processing
        try:
            return processing.run(f"qeht:{alg}", params, context=context, feedback=feedback,
                                  is_child_algorithm=True)
        except Exception as e:                  # noqa: BLE001 - re-raised with the step name
            raise QgsProcessingException(f"Step '{alg}' failed: {e}")

    @staticmethod
    def _copy_layers(package, folder, names):
        """Copy package layers to one GeoPackage each (in-process GDAL)."""
        from osgeo import gdal
        gdal.UseExceptions()
        out = {}
        for name in names:
            dst = os.path.join(folder, f"{name}.gpkg")
            if os.path.exists(dst):
                os.remove(dst)
            gdal.VectorTranslate(dst, package, options=gdal.VectorTranslateOptions(
                format="GPKG", layers=[name], layerName=name))
            out[name] = dst
        return out

    # -- run -------------------------------------------------------------------
    def processAlgorithm(self, parameters, context, feedback):
        t0 = time.time()
        p = self._settings(parameters, context, feedback)
        if not p.get("DEM"):
            raise QgsProcessingException("Give a DEM (or a settings file that names one).")
        folder = self.parameterAsFileOutput(
            parameters if parameters.get("OUTPUT") not in (None, "") else p, "OUTPUT", context)
        if not folder:
            raise QgsProcessingException("Choose an output folder.")
        d = {k: os.path.join(folder, k) for k in ("layers", "rasters", "tables", "report")}
        for v in d.values():
            os.makedirs(v, exist_ok=True)
        R = lambda name: os.path.join(d["rasters"], name)  # noqa: E731
        L = lambda name: os.path.join(d["layers"], name)   # noqa: E731
        raw = self.raster_path(p, "DEM", context)
        road = p.get("ROAD") or None
        crossings_in = p.get("CROSSINGS") or None
        stage = self.parameterAsEnum(p, "STAGE", context) if "STAGE" in p else 0
        burn = road is not None and self.parameterAsBool(p, "BURN", context)
        flat_idx = self.parameterAsEnum(p, "FLAT_METHOD", context) if "FLAT_METHOD" in p else 1
        warnings, outputs = [], {}

        def warn(msg):
            warnings.append(msg)
            feedback.pushWarning(msg)
        if road is None and crossings_in is None:
            raise QgsProcessingException("Give a road alignment, crossings / pour points, or both.")

        steps = 12
        fb = QgsProcessingMultiStepFeedback(steps, feedback)
        # 0. DEM checks -----------------------------------------------------------
        fb.setCurrentStep(0)
        fb.setProgressText("DEM checks")
        if p.get("NODATA_OVERRIDE") not in (None, ""):
            from osgeo import gdal
            gdal.UseExceptions()
            nd = self.parameterAsDouble(p, "NODATA_OVERRIDE", context)
            gdal.Translate(R("dem_input.tif"), raw, options=gdal.TranslateOptions(
                noData=nd, creationOptions=["COMPRESS=DEFLATE", "TILED=YES"]))
            feedback.pushInfo(f"DEM value {nd:g} declared NoData in rasters/dem_input.tif; "
                              "every step uses that copy.")
            raw = R("dem_input.tif")
        dem, valid, info = read_dem(raw)
        if not info.projection_wkt:
            raise QgsProcessingException("The DEM has no CRS; QEHT needs a projected, metric CRS.")
        self.check_size(feedback, info.rows, info.cols)
        audit = audit_nodata(dem, valid, info.nodata) + audit_resampling(dem, valid)
        for a in audit:
            warn("DEM check: " + a)
        cell_area = info.cell_width * info.cell_height
        km2 = self.parameterAsDouble(p, "STREAM_KM2", context) if p.get("STREAM_KM2") not in (None, "") else 1.0
        thr = max(1.0, round(km2 * 1e6 / cell_area))
        feedback.pushInfo(f"Stream threshold {km2:g} km2 = {thr:,.0f} cells of "
                          f"{info.cell_width:g} x {info.cell_height:g} m")
        del dem, valid

        def route(dem_path, tag):
            fb.setProgressText(f"Fill, flow direction, accumulation, streams ({tag})")
            self._run("filldepressions", {"DEM": dem_path, "MIN_SLOPE": 0.0,
                                          "OUTPUT": R("filled.tif")}, context, fb)
            self._run("flowdirection", {"DEM": R("filled.tif"), "RESOLVE_FLATS": True,
                                        "FLAT_METHOD": flat_idx,
                                        "OUTPUT": R("flow_direction.tif")}, context, fb)
            self._run("flowaccumulation", {"FDR": R("flow_direction.tif"), "QUANTITY": 0,
                                           "OUTPUT": R("flow_accumulation.tif")}, context, fb)
            self._run("streamnetwork", {"FDR": R("flow_direction.tif"),
                                        "FAC": R("flow_accumulation.tif"), "MODE": 0,
                                        "THRESHOLD": thr, "STREAMS": R("streams.tif"),
                                        "ORDER": R("strahler.tif")}, context, fb)

        # 1. routing ----------------------------------------------------------------
        fb.setCurrentStep(1)
        route(raw, "raw DEM")
        # 2. candidates -------------------------------------------------------------
        fb.setCurrentStep(2)
        cand_path = None
        if road is not None and crossings_in is None:
            fb.setProgressText("Road crossing candidates")
            self._run("crossingcandidates", {
                "FDR": R("flow_direction.tif"), "FAC": R("flow_accumulation.tif"),
                "STREAMS": R("streams.tif"), "THRESHOLD": thr, "ORDER": R("strahler.tif"),
                "ROAD": road, "START": p.get("START", 0.0), "REVERSE": p.get("REVERSE", False),
                "CANDIDATES": L("crossing_candidates.gpkg"),
                "PARALLEL": L("parallel_reaches.gpkg")}, context, fb)
            cand_path = L("crossing_candidates.gpkg")
            outputs["Crossing candidates"] = cand_path
        crossing_layer = crossings_in or cand_path
        # 3. burn and route again ---------------------------------------------------
        fb.setCurrentStep(3)
        burn_log = None
        route_dem = raw
        if burn:
            fb.setProgressText("Burn crossings through embankments")
            self._run("burncrossings", {"DEM": raw, "CROSSINGS": crossing_layer, "ROAD": road,
                                        "ALL": False, "OUTPUT": R("dem_burned.tif"),
                                        "LOG": L("burn_log.gpkg")}, context, fb)
            burn_log = L("burn_log.gpkg")
            route_dem = R("dem_burned.tif")
            route(route_dem, "burned DEM")
        fb.setCurrentStep(4)
        self._run("streamlines", {"FDR": R("flow_direction.tif"), "FAC": R("flow_accumulation.tif"),
                                  "DEM": R("filled.tif"), "MODE": 0, "THRESHOLD": thr,
                                  "OUTPUT": L("streams.gpkg")}, context, fb)
        outputs["Streams"] = L("streams.gpkg")
        settings_path = os.path.join(folder, "settings.json")
        self._write_settings(settings_path, p, context, folder)

        if stage == 1:
            if cand_path is None:
                warn("'Stop after crossing candidates' needs a road and no crossing layer; "
                     "the full run continues.")
            else:
                msg = ("Stopped for review. Edit layers/crossing_candidates.gpkg (status = "
                       "accepted / rejected; move or add points), then run again: Re-run from "
                       f"settings = {settings_path}, Crossings = the edited layer, Run = Full run.")
                feedback.pushInfo(msg)
                self._load(context, p, [("Crossing candidates", cand_path),
                                        ("Streams", L("streams.gpkg"))])
                return {"OUTPUT": folder, "REPORT": "", "CHARACTERISTICS": "", "PACKAGE_FILE": "",
                        "SUMMARY": msg}

        # 5. erosion ---------------------------------------------------------------
        fb.setCurrentStep(5)
        ero_folder = ""
        if self.parameterAsBool(p, "EROSION", context) if "EROSION" in p else True:
            fb.setProgressText("Erosion indices and RUSLE")
            ep = {"RAW_DEM": raw, "FAC": R("flow_accumulation.tif"), "CHANNEL_CELLS": thr,
                  "OUTPUT": os.path.join(d["rasters"], "erosion")}
            if p.get("R_RASTER"):
                ep["R_RASTER"] = p["R_RASTER"]
            elif p.get("R_VALUE") not in (None, ""):
                ep["R_VALUE"] = p["R_VALUE"]
            elif p.get("RAIN_MAP") and int(p.get("RAIN_R_RELATION") or 0) > 0:
                ep["R_MAP"] = p["RAIN_MAP"]
                ep["R_RELATION"] = int(p["RAIN_R_RELATION"]) - 1
            if p.get("LANDCOVER"):
                ep["WORLDCOVER"] = p["LANDCOVER"]
            ep.update({k: v for k, v in p.items() if k.startswith("SOIL") and v not in (None, "")})
            self._run("erosionindices", ep, context, fb)
            ero_folder = ep["OUTPUT"]
            outputs["Erosion rasters"] = ero_folder

        # 6. package (always built; kept only when asked) ----------------------------
        fb.setCurrentStep(6)
        fb.setProgressText("Crossings, catchments and flow paths")
        keep = self.parameterAsBool(p, "PACKAGE", context) if "PACKAGE" in p else False
        pdir = os.path.join(folder, "package" if keep else "_work")
        os.makedirs(pdir, exist_ok=True)
        package = os.path.join(pdir, "design_hydrology.gpkg")
        bp = {"FDR": R("flow_direction.tif"), "FAC": R("flow_accumulation.tif"), "RAW_DEM": raw,
              "ORDER": R("strahler.tif"), "POINTS": crossing_layer,
              "SNAP": 5, "SNAP_THRESHOLD": thr, "LOCAL": False,
              "DEM_SOURCE": p.get("DEM_SOURCE") or "",
              "FILLED": R("filled.tif"), "CSV": False, "OUTPUT": package,
              "COVERAGE": road is not None and (self.parameterAsBool(p, "COVERAGE", context)
                                                if "COVERAGE" in p else True)}
        if road is not None:
            bp.update({"ROAD": road, "START": p.get("START", 0.0), "REVERSE": p.get("REVERSE", False)})
        if ero_folder:
            bp["EROSION"] = ero_folder
        if burn_log:
            bp["BURN_LOG"] = burn_log
        bp.update({k: v for k, v in p.items() if k.startswith(PASS_PREFIXES) and v not in (None, "")})
        res = self._run("buildheasexchange", bp, context, fb)
        if res.get("SUMMARY"):
            feedback.pushInfo(res["SUMMARY"])

        # 7. flat-method check ------------------------------------------------------
        fb.setCurrentStep(7)
        fb.setProgressText("Flat-method check")
        from ..core.interop import gpkg
        from ..core.flow.sensitivity import flat_method_sensitivity
        filled, fvalid, finfo = read_dem(R("filled.tif"))
        crs_rows = gpkg.read_table(package, "crossings", with_geometry=False)
        gt = finfo.geotransform
        rc = [(int((c["outlet_y"] - gt[3]) // gt[5]), int((c["outlet_x"] - gt[0]) // gt[1]))
              for c in crs_rows]
        flat = {}
        if rc:
            res_f = flat_method_sensitivity(filled, fvalid, finfo.cell_width, finfo.cell_height, rc)
            flat = {c["outlet_uid"]: r for c, r in zip(crs_rows, res_f)}
            n_f = sum(1 for r in flat.values() if r.get("flat_sensitive"))
            if n_f:
                warn(f"{n_f} crossing(s) change area by more than 10 % between the two flat "
                     "methods: " + ", ".join(u for u, r in flat.items() if r.get("flat_sensitive"))
                     + " - verify against mapped drainage or on site.")
        del filled, fvalid

        # 8. stand-alone layers and tables ------------------------------------------
        fb.setCurrentStep(8)
        fb.setProgressText("Layers and tables")
        tabs = [t for t, _ in gpkg.list_tables(package)]
        names = [n for n in ("crossings", "catchments", "flowpaths", "coverage_check",
                             "sag_points", "flat_stretches", "alignment_profile",
                             "road_alignment") if n in tabs]
        outputs.update({f"Layer {k}": v for k, v in self._copy_layers(package, d["layers"], names).items()})
        from ..core.interop.heas_exchange import write_csvs
        for csv_path in write_csvs(package, d["tables"]):
            base = os.path.basename(csv_path).replace("design_hydrology_", "")
            os.replace(csv_path, os.path.join(d["tables"], base))
        from ..core.report.characteristics import write_characteristics_csv
        chars, n_rows = write_characteristics_csv(
            package, os.path.join(d["tables"], "catchment_characteristics.csv"), flat)
        outputs["Catchment characteristics"] = chars

        # 9. report -------------------------------------------------------------------
        fb.setCurrentStep(9)
        md = {r["key"]: r["value"] for r in gpkg.read_table(package, "qeht_run_metadata")}
        for w in (md.get("n_proposed") and int(md["n_proposed"]) and
                  [f"{md['n_proposed']} proposed crossing(s) from the coverage check - adopt or "
                   "delete them, then use Renumber and relink."] or []):
            warn(w)
        from ..core.report.run_report import write_report
        report = os.path.join(d["report"], "run_report.html")
        outputs["Settings"] = settings_path
        if keep:
            outputs["Design hydrology package"] = package
        write_report(package, report, title=p.get("RUN_NAME") or "QEHT hydrology run",
                     extra={"warnings": warnings, "dem_audit": audit, "flat_check": flat,
                            "settings_path": settings_path,
                            "outputs": {k: os.path.relpath(v, folder) for k, v in outputs.items()}})
        if not keep:
            shutil.rmtree(pdir, ignore_errors=True)

        # 10. load --------------------------------------------------------------------
        fb.setCurrentStep(10)
        self._load(context, p, [("Catchments", L("catchments.gpkg")),
                                ("Flow paths", L("flowpaths.gpkg")),
                                ("Streams", L("streams.gpkg")),
                                ("Crossings", L("crossings.gpkg"))]
                   + [(n.replace("_", " ").capitalize(), L(f"{n}.gpkg"))
                      for n in ("coverage_check", "sag_points", "flat_stretches") if n in names])
        summary = (f"{n_rows} crossing(s); {len(warnings)} warning(s); "
                   f"{time.time() - t0:,.0f} s. Report: {report}")
        feedback.pushInfo(summary)
        return {"OUTPUT": folder, "REPORT": report, "CHARACTERISTICS": chars,
                "PACKAGE_FILE": package if keep else "", "SUMMARY": summary}

    def _write_settings(self, path, p, context, folder):
        from ..core.interop.heas_exchange import plugin_version
        rec = {}
        for k, v in p.items():
            d = self.parameterDefinition(k)
            if d is None or k in ("SETTINGS",):
                continue
            text = "" if v is None else v
            kind = d.type()
            if text not in ("", None) and kind == "raster":
                lyr = self.parameterAsRasterLayer(p, k, context)
                text = lyr.source() if lyr is not None else str(text)
            elif text not in ("", None) and kind in ("source", "vector"):
                lyr = self.parameterAsVectorLayer(p, k, context)
                text = lyr.source() if lyr is not None else str(text)
            elif not isinstance(text, (str, int, float, bool)):
                text = str(text)
            rec[k] = text
        rec["OUTPUT"] = folder
        with open(path, "w", encoding="utf-8") as f:
            json.dump({"qeht_version": plugin_version(), "algorithm": "qeht:hydrologypipeline",
                       "parameters": rec}, f, indent=1, sort_keys=True)

    def _load(self, context, p, items):
        if "LOAD" in p and not self.parameterAsBool(p, "LOAD", context):
            return
        group = f"QEHT {p.get('RUN_NAME') or 'run'}"
        for label, path in items:
            if not os.path.exists(path):
                continue
            details = QgsProcessingContext.LayerDetails(label, QgsProject.instance(), label)
            if hasattr(details, "groupName"):
                details.groupName = group
            context.addLayerToLoadOnCompletion(path, details)
