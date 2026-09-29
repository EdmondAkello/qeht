# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Build HEAS exchange package - one self-describing GeoPackage per run."""

import os

from qgis.core import (
    QgsProcessing, QgsProcessingException, QgsProcessingContext,
    QgsProcessingParameterRasterLayer, QgsProcessingParameterFeatureSource,
    QgsProcessingParameterField, QgsProcessingParameterString,
    QgsProcessingParameterEnum, QgsProcessingParameterNumber,
    QgsProcessingParameterBoolean, QgsProcessingParameterFileDestination,
    QgsProcessingOutputString,
    QgsCoordinateReferenceSystem,
)
import numpy as np

from .base import QehtAlgorithm
from ..core.raster import read_dem
from ..core.grid import decode_d8
from ..core.watershed.delineate import extract_streams
from ..core.network.crossings import select_crossings, CANDIDATE_FIELDS
from ..core.interop.heas_exchange import (build_exchange_records, write_exchange,
                                          validate_exchange, crs_check,
                                          file_fingerprint, plugin_version,
                                          ExchangeError, SCHEMA_VERSION)

FDR = "FDR"; DEM = "DEM"; RAW_DEM = "RAW_DEM"; FAC = "FAC"; ORDER = "ORDER"
ROAD = "ROAD"; START = "START"; REVERSE = "REVERSE"
POINTS = "POINTS"; ID_FIELD = "ID_FIELD"; ID_PREFIX = "ID_PREFIX"; ID_ORDER = "ID_ORDER"
SNAP = "SNAP"; SNAP_THRESHOLD = "SNAP_THRESHOLD"; LOCAL = "LOCAL"
DEM_SOURCE = "DEM_SOURCE"; CONDITIONING = "CONDITIONING"; FLAT_METHOD = "FLAT_METHOD"
CSV = "CSV"; OUTPUT = "OUTPUT"; SUMMARY = "SUMMARY"

ORDER_OPTIONS = ["Automatic: along the chainage when known, else downstream first",
                 "Downstream first (largest contributing area = 001)",
                 "Pour-point layer order"]
ORDER_KEYS = ["auto", "downstream", "input"]


class BuildHeasExchangeAlgorithm(QehtAlgorithm):

    def name(self): return "buildheasexchange"
    def displayName(self): return "Build HEAS exchange package"
    def group(self): return "Interoperability"
    def groupId(self): return "interop"

    def shortHelpString(self):
        return (
            "Runs snapping, catchment delineation, longest flow path and "
            "catchment characteristics for every pour point and writes ONE "
            f"GeoPackage (schema {SCHEMA_VERSION}) that HEAS imports with no field "
            "mapping.\n\n"
            "<b>Layers:</b> crossings (snapped outlets), catchments, flowpaths, "
            "plus qeht_run_metadata (DEM, method, thresholds, QEHT version) and "
            "qeht_field_dictionary (meaning, unit and method of every field).\n\n"
            "<b>outlet_uid</b> links the three layers. It is assigned once per "
            "run: taken from the ID attribute you choose (e.g. existing culvert "
            "numbers), or sequential with your prefix (X001, X002 ...), numbered "
            "downstream-first by default. Duplicate or empty IDs stop the run "
            "with a list. <i>outlet_id</i> (the pour-point feature id) is still "
            "written but is NOT stable - do not link on it.\n\n"
            "<b>A projected, metric CRS is required.</b> A DEM in degrees is "
            "refused; reproject it (e.g. to UTM) first.\n\n"
            "<b>Crossing candidates:</b> give the layer from 'Road crossing "
            "candidates' as the crossings. The accepted candidates are used (or the "
            "recommended ones if none is accepted); each keeps its outlet cell "
            "(no snapping) and its chainage, IDs are numbered along the chainage, "
            "and the full candidate layer is stored in the package as "
            "crossing_candidates. With hand-placed points, give the road alignment "
            "to get chainages.\n\n"
            "<b>Catchments:</b> full upstream area per crossing by default "
            "(overlapping - what a culvert design flow needs). 'Local' gives "
            "non-overlapping areas, computed upstream-first so the result does "
            "not depend on the order of the pour-point layer.\n\n"
            "<b>Slopes (four domains, never merged):</b> catch_slope_horn (Horn "
            "mean), catch_relief_ratio (relief / LFP length), lfp_slope (drop / "
            "length) and lfp_slope_1085 (10-85 along the LFP, measured from the "
            "OUTLET; lfp_L10_m, lfp_L85_m, lfp_z10_m, lfp_z85_m let you check it "
            "by hand).\n\n"
            "Supply the RAW DEM for reported elevations and slopes; routing uses "
            "the flow-direction grid. Everything runs in-process.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(FDR, "Flow direction (D8-coded)"))
        self.addParameter(QgsProcessingParameterRasterLayer(FAC, "Flow accumulation"))
        self.addParameter(QgsProcessingParameterRasterLayer(
            RAW_DEM, "Raw DEM (for reported elevations and slopes)"))
        self.addParameter(QgsProcessingParameterRasterLayer(
            ORDER, "Stream order raster (optional, for crossing attributes)", optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource(
            POINTS, "Crossings / pour points", [QgsProcessing.SourceType.TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterFeatureSource(
            ROAD, "Road alignment (optional; chainage for hand-placed points)",
            [QgsProcessing.SourceType.TypeVectorLine], optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            START, "Start chainage (m)", QgsProcessingParameterNumber.Type.Double,
            defaultValue=0.0))
        self.addParameter(QgsProcessingParameterBoolean(
            REVERSE, "Reverse chainage direction", defaultValue=False))
        self.addParameter(QgsProcessingParameterField(
            ID_FIELD, "ID attribute (optional; blank = sequential IDs)",
            parentLayerParameterName=POINTS, optional=True))
        self.addParameter(QgsProcessingParameterString(
            ID_PREFIX, "ID prefix", defaultValue="X", optional=True))
        self.addParameter(QgsProcessingParameterEnum(
            ID_ORDER, "Sequential numbering order", options=ORDER_OPTIONS, defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber(
            SNAP, "Snap radius (cells)", QgsProcessingParameterNumber.Type.Integer,
            defaultValue=5, minValue=0, maxValue=100))
        self.addParameter(QgsProcessingParameterNumber(
            SNAP_THRESHOLD, "Snap-to-stream threshold (cells; 0 = max accumulation)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=200.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterBoolean(
            LOCAL, "Non-overlapping (local) catchments", defaultValue=False))
        self.addParameter(QgsProcessingParameterString(
            DEM_SOURCE, "DEM source (recorded in metadata, e.g. 'ALOS 30 m')", optional=True))
        self.addParameter(QgsProcessingParameterString(
            CONDITIONING, "Conditioning used (recorded, e.g. 'fill, min_slope 0')",
            optional=True))
        self.addParameter(QgsProcessingParameterString(
            FLAT_METHOD, "Flat resolution used for flow direction (recorded, e.g. 'toward')",
            optional=True))
        self.addParameter(QgsProcessingParameterBoolean(
            CSV, "Also write one CSV per layer (for spreadsheets)", defaultValue=False))
        self.addParameter(QgsProcessingParameterFileDestination(
            OUTPUT, "Exchange GeoPackage", fileFilter="GeoPackage (*.gpkg)"))
        self.addOutput(QgsProcessingOutputString(SUMMARY, "Summary"))

    def processAlgorithm(self, parameters, context, feedback):
        fdr_path = self.raster_path(parameters, FDR, context)
        fac_path = self.raster_path(parameters, FAC, context)
        raw_path = self.raster_path(parameters, RAW_DEM, context)
        out_path = self.parameterAsFileOutput(parameters, OUTPUT, context)
        if not out_path.lower().endswith(".gpkg"):
            out_path += ".gpkg"

        d8, valid, info = read_dem(fdr_path)

        # CRS rule first - fail before any heavy work. The WKT check decides
        # (geographic refused, metres required); QGIS only supplies a clean
        # EPSG code and name when it recognises the CRS.
        ok, epsg, msg = crs_check(info.projection_wkt)
        dem_crs = QgsCoordinateReferenceSystem()
        dem_crs.createFromWkt(info.projection_wkt)
        if dem_crs.isValid() and dem_crs.isGeographic():
            ok, msg = False, ("The DEM is in a geographic CRS (degrees). Reproject it to "
                              "a projected, metric CRS (e.g. UTM) and re-run.")
        if not ok:
            raise QgsProcessingException(msg)
        if dem_crs.isValid() and dem_crs.authid().upper().startswith("EPSG:"):
            epsg = int(dem_crs.authid().split(":")[1])

        self.check_size(feedback, info.rows, info.cols)
        direction = decode_d8(d8.astype(np.int32))
        accum, _, _ = read_dem(fac_path)
        elevation, raw_valid, raw_info = read_dem(raw_path)
        if (raw_info.rows, raw_info.cols) != (info.rows, info.cols):
            raise QgsProcessingException(
                "The raw DEM and the flow-direction grid have different sizes. "
                "Use the raw DEM the flow direction was derived from.")
        # NoData in the raw DEM must not become an elevation (-9999 would
        # pass as the catchment minimum): mark it NaN, statistics skip NaN.
        elevation = np.where(raw_valid, elevation, np.nan)
        stream_order = None
        if self.parameterAsRasterLayer(parameters, ORDER, context) is not None:
            so, so_valid, _ = read_dem(self.raster_path(parameters, ORDER, context))
            stream_order = np.where(so_valid, so, np.nan)

        id_field = self.field_parameter(parameters, ID_FIELD, context)
        prefix = (self.parameterAsString(parameters, ID_PREFIX, context) or "").strip()
        order = ORDER_KEYS[self.parameterAsEnum(parameters, ID_ORDER, context)]
        snap_radius = self.parameterAsInt(parameters, SNAP, context)
        snap_threshold = self.parameterAsDouble(parameters, SNAP_THRESHOLD, context)
        local = self.parameterAsBool(parameters, LOCAL, context)

        points = self.read_pour_points(parameters, POINTS, context, info, feedback,
                                       id_field=id_field,
                                       extra_fields=[f for f, _ in CANDIDATE_FIELDS])
        extra_layers = []
        candidate_mode = any("attr_status" in p and "attr_outlet_x" in p for p in points)
        if candidate_mode:
            extra_layers.append(("crossing_candidates", "POINT", CANDIDATE_FIELDS,
                                 [((p["x"], p["y"]), {f: p.get("attr_" + f)
                                                      for f, _ in CANDIDATE_FIELDS})
                                  for p in points],
                                 "All crossing candidates (audit trail)"))
            idx, rule = select_crossings([{"status": p.get("attr_status"),
                                           "recommended": p.get("attr_recommended")}
                                          for p in points])
            points = [points[k] for k in idx]
            if not points:
                raise QgsProcessingException(
                    "No candidate is accepted or recommended - nothing to export.")
            for p in points:
                p["outlet_x"], p["outlet_y"] = p.get("attr_outlet_x"), p.get("attr_outlet_y")
                p["chainage"] = p.get("attr_chainage_m")
            feedback.pushInfo(f"Candidate layer: {len(points)} crossing(s) selected ({rule}).")

        alignment = self.read_alignment(
            parameters, ROAD, context, info, feedback,
            start_chainage=self.parameterAsDouble(parameters, START, context),
            reverse=self.parameterAsBool(parameters, REVERSE, context))
        if alignment is not None:
            missing = [p for p in points if p.get("chainage") is None]
            if missing:
                ch, off, _ = alignment.locate([p["x"] for p in missing], [p["y"] for p in missing])
                for p, c in zip(missing, ch):
                    p["chainage"] = float(c)
            extra_layers.append((
                "road_alignment", "LINESTRING",
                [("part", "int"), ("chainage_from_m", "real"), ("chainage_to_m", "real")],
                [([tuple(v) for v in part], {"part": k + 1,
                  "chainage_from_m": float(alignment.ch0[alignment.part == k].min()),
                  "chainage_to_m": float((alignment.ch0 + alignment.seg_len)[alignment.part == k].max())})
                 for k, part in enumerate(alignment.parts)],
                "Road alignment used for chainage"))
        if order == "auto":
            order = ("chainage" if all(p.get("chainage") is not None for p in points)
                     else "downstream")
        elif order == "chainage" and not all(p.get("chainage") is not None for p in points):
            raise QgsProcessingException("Chainage numbering needs a candidate layer or a road alignment.")
        stream_mask = None
        if snap_radius > 0 and snap_threshold > 0:
            stream_mask = extract_streams(accum, valid, threshold_cells=snap_threshold)

        feedback.pushInfo(f"{len(points)} pour points; "
                          f"{'local' if local else 'full upstream'} catchments; "
                          f"IDs {'from ' + id_field if id_field else 'sequential ' + (prefix or '') + '001...'}")
        try:
            crossings, catchments, flowpaths, issues, id_info = build_exchange_records(
                direction, valid, accum, elevation, info.geotransform, points,
                snap_radius_cells=snap_radius, stream_mask=stream_mask, local=local,
                stream_order=stream_order,
                id_scheme="attribute" if id_field else "sequential",
                id_prefix=prefix, id_order=order,
                progress=self.make_progress(feedback, weight=0.9))
        except ExchangeError as e:
            raise QgsProcessingException(str(e))
        for msg in issues:
            feedback.pushWarning(msg)

        cell_area = info.cell_width * info.cell_height
        md = dict(id_info)
        md.update({
            "qeht_version": plugin_version(),
            "cell_size_m": f"{info.cell_width:g} x {info.cell_height:g}",
            "dem_path": fdr_path, "raw_dem_path": raw_path,
            "dem_sha256": file_fingerprint(raw_path),
            "dem_source": self.parameterAsString(parameters, DEM_SOURCE, context) or "",
            "conditioning": self.parameterAsString(parameters, CONDITIONING, context) or "not recorded",
            "flat_method": self.parameterAsString(parameters, FLAT_METHOD, context) or "not recorded",
            "tie_rule": "QEHT D8: steepest drop/distance; ties by lowest internal index",
            "stream_threshold_cells": f"{snap_threshold:g}",
            "stream_threshold_km2": f"{snap_threshold * cell_area / 1e6:g}",
            "snap_radius_cells": str(snap_radius),
            "snap_strategy": ("none" if snap_radius <= 0 else
                              "nearest_stream" if stream_mask is not None else "max_accumulation"),
            "catchment_mode": "local" if local else "full",
            "id_attribute": id_field or "",
            "crossing_source": ("crossing candidates" if candidate_mode else "pour points"),
            "chainage_start_m": (f"{alignment.start_chainage:g}" if alignment is not None else ""),
            "parameters_json": {k: str(v) for k, v in parameters.items()},
        })
        try:
            write_exchange(out_path, crossings, catchments, flowpaths, md,
                           info.projection_wkt, crs_epsg=epsg,
                           crs_name=dem_crs.description() if dem_crs.isValid() else "",
                           csv_dir=(os.path.splitext(out_path)[0] + "_csv")
                           if self.parameterAsBool(parameters, CSV, context) else None,
                           extra_layers=extra_layers)
        except ExchangeError as e:
            raise QgsProcessingException(str(e))

        errors, warnings = validate_exchange(out_path)
        for w in warnings:
            feedback.pushWarning(w)
        if errors:
            raise QgsProcessingException("The package failed validation: " + "; ".join(errors))

        for _, a in catchments:
            feedback.pushInfo(
                f"  {a['outlet_uid']}: A={a['area_km2']:.4f} km2  "
                f"Shorn={a.get('catch_slope_horn', float('nan')):.4f}  "
                f"RR={a.get('catch_relief_ratio', float('nan')):.4f}  "
                f"L={a.get('lfp_length_km', float('nan')):.3f} km")

        stem = os.path.splitext(os.path.basename(out_path))[0]
        for lname in ("catchments", "flowpaths", "crossings"):
            try:
                context.addLayerToLoadOnCompletion(
                    f"{out_path}|layername={lname}",
                    QgsProcessingContext.LayerDetails(f"{stem} - {lname}",
                                                      context.project(), lname))
            except Exception as e:      # loading is a convenience only
                feedback.pushInfo(f"(Could not auto-load {lname}: {e})")

        summary = (f"{len(crossings)} crossings, {len(catchments)} catchments, "
                   f"{len(flowpaths)} flow paths -> {out_path} ({SCHEMA_VERSION})")
        feedback.pushInfo(summary)
        return {OUTPUT: out_path, SUMMARY: summary}
