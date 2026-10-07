# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Shared base for QEHT Processing algorithms.

This is the ONLY place where QGIS and `qeht.core` meet. The core never
imports QGIS; the algorithms below never implement hydrology. If you
find yourself writing a loop over cells in this package, it belongs in
core instead.
"""

import os

from qgis.core import (
    QgsProcessingAlgorithm,
    QgsProcessingException,
)


class QehtAlgorithm(QgsProcessingAlgorithm):
    """Common metadata and helpers."""

    def group(self):
        return "Terrain and drainage"

    def groupId(self):
        return "terrain"

    def createInstance(self):
        return type(self)()

    # -- helpers ----------------------------------------------------

    def raster_path(self, parameters, name, context):
        """Resolve a raster parameter to a filesystem path.

        We read through GDAL directly rather than through the QGIS raster
        provider, so we need a real path. In-memory and remote layers are
        rejected with a clear message rather than failing obscurely later.
        """
        layer = self.parameterAsRasterLayer(parameters, name, context)
        if layer is None:
            raise QgsProcessingException(f"Could not resolve raster parameter '{name}'.")
        path = layer.source().split("|")[0]
        import os
        if not os.path.exists(path):
            raise QgsProcessingException(
                f"Layer '{layer.name()}' is not backed by a file GDAL can open "
                f"directly ({path}). Save it to GeoTIFF first.")
        return path

    def make_progress(self, feedback, weight=1.0, offset=0.0):
        """Adapt core's progress(fraction, message) to QGIS feedback."""
        def _progress(fraction, message=""):
            if feedback.isCanceled():
                raise QgsProcessingException("Cancelled by user.")
            feedback.setProgress(int((offset + fraction * weight) * 100))
            if message:
                feedback.setProgressText(message)
        return _progress

    def check_size(self, feedback, rows, cols, steps=None):
        """Report the estimated peak memory before a run; warn when it is
        close to the free RAM (or large and free RAM is unknown).

        QEHT processes the grid whole - there is no tiling. Estimates come
        from measured per-step peaks (core/memory.py). Clipping to the
        catchment plus a buffer is standard practice and does not change
        results inside the clip.
        """
        from ..core.memory import estimate_peak_memory, available_memory_gb
        est = estimate_peak_memory(rows, cols, steps)
        free = available_memory_gb()
        feedback.pushInfo(f"Grid {rows:,} x {cols:,} = {est['megapixels']:,.1f} megapixels; "
                          f"estimated peak memory {est['peak_gb']:,.2f} GB"
                          + (f" (free: {free:,.1f} GB)" if free else ""))
        if (free and est["peak_gb"] > 0.8 * free) or (free is None and est["peak_gb"] > 8.0):
            feedback.pushWarning(
                f"This run needs roughly {est['peak_gb']:,.1f} GB at peak"
                + (f" and only {free:,.1f} GB is free" if free else "")
                + ". If it runs out of memory, clip the DEM to your catchment plus a buffer "
                "and run again - results inside the clip are unchanged.")
        return est["megapixels"]

    def read_alignment(self, parameters, name, context, info, feedback,
                       start_chainage=0.0, reverse=False):
        """Road line layer -> core Alignment in the DEM CRS (None if absent).

        Features are chained in layer order; a multi-part feature contributes
        its parts in order. Reverse flips the direction of chainage.
        """
        from qgis.core import (QgsCoordinateReferenceSystem,
                               QgsCoordinateTransform, QgsProject)
        from ..core.network.alignment import Alignment
        source = self.parameterAsSource(parameters, name, context)
        if source is None:
            return None
        dem_crs = QgsCoordinateReferenceSystem()
        dem_crs.createFromWkt(info.projection_wkt)
        transform = None
        if dem_crs.isValid() and source.sourceCrs() != dem_crs:
            transform = QgsCoordinateTransform(source.sourceCrs(), dem_crs,
                                               QgsProject.instance())
        parts = []
        for feature in source.getFeatures():
            geom = feature.geometry()
            if geom is None or geom.isEmpty():
                continue
            if transform is not None:
                geom.transform(transform)
            lines = geom.asMultiPolyline() if geom.isMultipart() else [geom.asPolyline()]
            for line in lines:
                if len(line) >= 2:
                    parts.append([(p.x(), p.y()) for p in line])
        if not parts:
            raise QgsProcessingException("The road layer has no line geometry.")
        al = Alignment(parts, start_chainage=start_chainage, reverse=reverse)
        feedback.pushInfo(f"Alignment: {len(al.parts)} part(s), {al.length:,.1f} m, "
                          f"chainage {al.start_chainage:,.1f} -> {al.end_chainage:,.1f}")
        return al

    def read_pour_points(self, parameters, name, context, info, feedback,
                         id_field=None, extra_fields=()):
        """Pour points -> list of dicts {x, y, fid, source_id} in the DEM CRS.

        Shared by every tool that takes pour points, so reprojection,
        multipoint handling and out-of-grid checks cannot drift between
        tools. `source_id` is the value of `id_field` (None if not given).
        """
        from qgis.core import (QgsCoordinateReferenceSystem,
                               QgsCoordinateTransform, QgsProject)
        source = self.parameterAsSource(parameters, name, context)
        if source is None:
            raise QgsProcessingException("Could not read the pour-point layer.")
        dem_crs = QgsCoordinateReferenceSystem()
        dem_crs.createFromWkt(info.projection_wkt)
        transform = None
        if dem_crs.isValid() and source.sourceCrs() != dem_crs:
            transform = QgsCoordinateTransform(source.sourceCrs(), dem_crs,
                                               QgsProject.instance())
            feedback.pushInfo(f"Reprojecting pour points "
                              f"{source.sourceCrs().authid()} -> {dem_crs.authid()}")
        points = []
        for feature in source.getFeatures():
            geom = feature.geometry()
            if geom is None or geom.isEmpty():
                feedback.pushWarning(f"Pour point fid {feature.id()} has no geometry. Skipped.")
                continue
            if transform is not None:
                geom.transform(transform)
            pt = geom.asMultiPoint()[0] if geom.isMultipart() else geom.asPoint()
            row, col = info.xy_to_rowcol(pt.x(), pt.y())
            if not (0 <= row < info.rows and 0 <= col < info.cols):
                feedback.pushWarning(f"Pour point fid {feature.id()} "
                                     f"({pt.x():.1f}, {pt.y():.1f}) falls outside the DEM. Skipped.")
                continue
            sid = None
            if id_field:
                v = feature[id_field]
                sid = None if v is None or str(v) == "NULL" else v
            rec = {"x": pt.x(), "y": pt.y(), "fid": int(feature.id()), "source_id": sid}
            names = feature.fields().names()
            for f in extra_fields:
                # a shapefile truncates field names to 10 characters
                src_name = f if f in names else (f[:10] if f[:10] in names else None)
                if src_name is not None:
                    v = feature[src_name]
                    rec["attr_" + f] = None if v is None or str(v) == "NULL" else v
            points.append(rec)
        if not points:
            raise QgsProcessingException("No usable pour points.")
        return points

    def field_parameter(self, parameters, name, context):
        """First field chosen in a field parameter, or None.

        parameterAsStrings() replaces parameterAsFields() from QGIS 3.32;
        the plugin still supports 3.22, so use whichever exists.
        """
        getter = getattr(self, "parameterAsStrings", None) or self.parameterAsFields
        values = getter(parameters, name, context)
        return values[0] if values else None

    def outlet_uids(self, points, outlet_rc, accumulation, id_field, prefix):
        """outlet_uid per outlet for the single-purpose tools.

        Attribute IDs when `id_field` is set; otherwise sequential with
        `prefix`, downstream-first when an accumulation grid is available,
        else in pour-point order. Duplicates raise with a list.
        """
        from ..core.linking.ids import assign_uids, OutletIdError
        try:
            if id_field:
                uids, _ = assign_uids(len(points), scheme="attribute", prefix="",
                                      attribute_values=[p["source_id"] for p in points])
            elif accumulation is not None:
                uids, _ = assign_uids(len(points), prefix=prefix, order="downstream",
                                      accumulation=[float(accumulation[r, c])
                                                    for r, c in outlet_rc])
            else:
                uids, _ = assign_uids(len(points), prefix=prefix, order="input")
        except OutletIdError as e:
            raise QgsProcessingException(str(e))
        return uids

    # -- soils (v0.11) ------------------------------------------------------

    def parameters_record(self, parameters, context):
        """Parameters for the run metadata, with layers recorded by their file
        source (not the project layer id) and empty inputs as ''."""
        out = {}
        for k, v in parameters.items():
            text = "" if v is None else str(v)
            if text in ("NULL", "None"):
                text = ""
            d = self.parameterDefinition(k)
            kind = d.type() if d is not None else ""
            if text and kind == "raster":
                lyr = self.parameterAsRasterLayer(parameters, k, context)
                text = lyr.source() if lyr is not None else text
            elif text and kind in ("source", "vector"):
                lyr = self.parameterAsVectorLayer(parameters, k, context)
                text = lyr.source() if lyr is not None else text
            out[k] = text
        return out

    def candidate_selection(self, points, feedback):
        """Pour points read from a 'Road crossing candidates' layer: keep the
        accepted candidates (or the recommended ones), each at its own outlet
        cell (no snapping) with its chainage. Returns (points, is_candidate_layer,
        all_points). Other layers pass through unchanged."""
        from ..core.network.crossings import select_crossings
        is_cand = any("attr_status" in p and "attr_outlet_x" in p for p in points)
        if not is_cand:
            return points, False, points
        idx, rule = select_crossings([{"status": p.get("attr_status"),
                                       "recommended": p.get("attr_recommended")} for p in points])
        chosen = [points[k] for k in idx]
        if not chosen:
            raise QgsProcessingException(
                "No candidate is accepted or recommended - nothing to export.")
        for p in chosen:
            p["outlet_x"], p["outlet_y"] = p.get("attr_outlet_x"), p.get("attr_outlet_y")
            p["chainage"] = p.get("attr_chainage_m")
        feedback.pushInfo(f"Candidate layer: {len(chosen)} of {len(points)} crossing(s) used ({rule}).")
        return chosen, True, points

    def _advanced(self, param):
        """Mark a parameter as advanced (QGIS 3.22 - 4)."""
        try:
            from qgis.core import Qgis
            flag = Qgis.ProcessingParameterFlag.Advanced
        except AttributeError:
            from qgis.core import QgsProcessingParameterDefinition
            flag = QgsProcessingParameterDefinition.FlagAdvanced
        param.setFlags(param.flags() | flag)
        return param

    def add_soil_parameters(self, prefix="SOIL", optional=True):
        """Soil inputs (v0.15: any source).

        Units: soil polygons (+ SOTWIS database or CSV table, or texture fields
        mapped by name), or a soil-unit raster + CSV table. Texture rasters:
        sand / clay / (silt) / organic carbon (+ bulk density, coarse
        fragments), or a SoilGrids folder with depth layers. Overrides: a
        hydrologic soil group raster or field, a K raster or field.
        """
        from qgis.core import (QgsProcessingParameterFeatureSource, QgsProcessing,
                               QgsProcessingParameterFile, QgsProcessingParameterString,
                               QgsProcessingParameterNumber, QgsProcessingParameterField,
                               QgsProcessingParameterRasterLayer, QgsProcessingParameterEnum)
        P = prefix
        self.addParameter(QgsProcessingParameterFeatureSource(
            P + "_POLYGONS", "Soil map polygons (SOTWIS / SOTER, national map, any polygons)",
            [QgsProcessing.SourceType.TypeVectorPolygon], optional=optional))
        self.addParameter(QgsProcessingParameterFile(
            P + "_DB", "SOTWIS / SOTER SQLite database (optional; full unit composition)",
            optional=True, fileFilter="SQLite (*.db *.sqlite)"))
        self.addParameter(QgsProcessingParameterFile(
            P + "_CSV", "Soil table CSV (optional; unit field + sand, silt, clay, oc %, "
            "[bulk, cfrag, drain])", optional=True, fileFilter="CSV (*.csv *.txt)"))
        self.addParameter(QgsProcessingParameterString(
            P + "_UNIT_FIELD", "Soil unit field linking polygons / unit raster to the database or table",
            defaultValue="NEWSUID", optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            P + "_TOP", "Soil depth from (cm)", QgsProcessingParameterNumber.Type.Double,
            defaultValue=0.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber(
            P + "_BOTTOM", "Soil depth to (cm)", QgsProcessingParameterNumber.Type.Double,
            defaultValue=20.0, minValue=1.0))
        adv = self._advanced
        # G2 - polygon texture fields by name
        for key, label in (("SAND", "sand %"), ("SILT", "silt %"), ("CLAY", "clay %"),
                           ("OC", "organic carbon"), ("BULK", "bulk density g/cm3"),
                           ("CFRAG", "coarse fragments vol %"), ("DRAIN", "FAO drainage class"),
                           ("HSG_FIELD", "hydrologic soil group (A, B, C, D, A/D ...)"),
                           ("K_FIELD", "USLE K, SI units")):
            self.addParameter(adv(QgsProcessingParameterField(
                f"{P}_F_{key}", f"Soil polygons: field with {label} (optional)",
                parentLayerParameterName=P + "_POLYGONS", optional=True)))
        self.addParameter(adv(QgsProcessingParameterEnum(
            P + "_OC_UNITS", "Soil polygons: organic carbon field units",
            options=["percent", "g/kg"], defaultValue=0)))
        # G3 - rasters
        self.addParameter(adv(QgsProcessingParameterFile(
            P + "_SOILGRIDS", "SoilGrids folder (sand/silt/clay/soc/bdod/cfvo_<depth>cm_mean.tif; "
            "depth-weighted)", behavior=QgsProcessingParameterFile.Behavior.Folder, optional=True)))
        for key, label in (("SAND_R", "Sand raster"), ("CLAY_R", "Clay raster"),
                           ("SILT_R", "Silt raster (optional; else 100 - sand - clay)"),
                           ("OC_R", "Organic carbon raster"), ("BULK_R", "Bulk density raster"),
                           ("CFRAG_R", "Coarse fragments raster")):
            self.addParameter(adv(QgsProcessingParameterRasterLayer(
                f"{P}_{key}", label + " (optional)", optional=True)))
        self.addParameter(adv(QgsProcessingParameterEnum(
            P + "_RASTER_UNITS", "Units of the soil rasters",
            options=["SoilGrids 2.0 (g/kg; soc dg/kg; bdod cg/cm3; cfvo per mille)",
                     "Percent (sand/silt/clay/oc %), g/cm3, vol %"], defaultValue=0)))
        self.addParameter(adv(QgsProcessingParameterRasterLayer(
            P + "_UNIT_RASTER", "Soil unit raster (codes joined to the CSV table, e.g. HWSD v2)",
            optional=True)))
        # G4 - direct overrides
        self.addParameter(adv(QgsProcessingParameterRasterLayer(
            P + "_HSG_R", "Hydrologic soil group raster (HYSOGs250m codes 1-4, 11-14)",
            optional=True)))
        self.addParameter(adv(QgsProcessingParameterRasterLayer(
            P + "_K_R", "USLE K raster, SI (overrides K computed from texture)", optional=True)))

    def load_soil(self, parameters, context, info, feedback, prefix="SOIL"):
        """Soil inputs -> core.soils.sources.SoilGrid on the DEM grid, or None.

        Units (polygons with the SOTWIS database, a CSV table, SOTWIS-named,
        plain or user-mapped texture fields; or a unit raster + CSV) are
        painted onto their cells; texture rasters give per-cell values; an
        HSG raster/field and a K raster/field override the proxies.
        """
        from ..core.soils.sources import (SoilGrid, hsg_from_codes, soilgrids_layers,
                                          depth_weighted, SOILGRIDS, SOILGRIDS_FACTORS,
                                          SOILGRIDS_PROP, HYSOGS)
        from ..core.raster import warp_to_grid
        import numpy as np
        P = prefix

        def has(name):
            return name in parameters and parameters[name] not in (None, "")

        def rpath(name):
            if not has(name) or self.parameterAsRasterLayer(parameters, name, context) is None:
                return None
            return self.raster_path(parameters, name, context)

        top = self.parameterAsDouble(parameters, P + "_TOP", context)
        bottom = self.parameterAsDouble(parameters, P + "_BOTTOM", context)
        if bottom <= top:
            raise QgsProcessingException("Soil depth 'to' must be greater than 'from'.")
        unit_field = (self.parameterAsString(parameters, P + "_UNIT_FIELD", context)
                      or "NEWSUID").strip()
        csv_path = self.parameterAsFile(parameters, P + "_CSV", context) if has(P + "_CSV") else ""
        sg = None
        extra = {}

        src = self.parameterAsSource(parameters, P + "_POLYGONS", context) \
            if has(P + "_POLYGONS") else None
        unit_raster = rpath(P + "_UNIT_RASTER")
        sg_dir = self.parameterAsFile(parameters, P + "_SOILGRIDS", context) \
            if has(P + "_SOILGRIDS") else ""
        sand_r = rpath(P + "_SAND_R")
        if src is not None:
            sg, extra = self._soil_from_polygons(parameters, context, info, feedback, P, src,
                                                 unit_field, csv_path, top, bottom)
        elif unit_raster:
            if not csv_path:
                raise QgsProcessingException("A soil unit raster needs the soil table CSV.")
            from ..core.soils.sotwis import load_csv_units
            codes, desc = warp_to_grid(unit_raster, info, resampling="near")
            try:
                units, sinfo = load_csv_units(csv_path, unit_field)
            except ValueError as e:
                raise QgsProcessingException(str(e))
            keys = {}
            grid = np.zeros(codes.shape, dtype=np.int32)
            fin = np.isfinite(codes)
            for v in np.unique(codes[fin]):
                key = (str(int(v)) if float(v).is_integer() else str(v))
                keys[key] = len(keys) + 1
                grid[fin & (codes == v)] = keys[key]
            idx = {i: k for k, i in keys.items()}
            missing = [k for k in keys if k not in units]
            sinfo = dict(sinfo, soil_dataset=f"{sinfo['soil_dataset']}; units from {desc}",
                         soil_depth_cm="as tabulated")
            feedback.pushInfo(f"Soil: {len(keys)} unit code(s) on the grid from {desc}"
                              + (f"; {len(missing)} without a table row (e.g. {missing[:5]})"
                                 if missing else ""))
            sg = SoilGrid.from_units(grid, idx, units, sinfo)
        elif sg_dir or sand_r:
            soilgrids_units = self.parameterAsEnum(parameters, P + "_RASTER_UNITS", context) == 0 \
                if has(P + "_RASTER_UNITS") else True
            read = lambda p, rs="bilinear": warp_to_grid(p, info, resampling=rs)[0]
            arrays = {}
            used = []
            if sg_dir:
                layers = soilgrids_layers(sg_dir)
                if not layers:
                    raise QgsProcessingException(
                        f"No SoilGrids files (e.g. sand_0-5cm_mean.tif) found in {sg_dir}.")
                for prop, lay in layers.items():
                    a, ints = depth_weighted(lay, top, bottom, read)
                    if a is not None:
                        arrays[SOILGRIDS_PROP[prop]] = a / SOILGRIDS_FACTORS[prop]
                        used.append(f"{prop} {','.join(ints)} cm")
                label = f"{SOILGRIDS} folder {os.path.basename(os.path.normpath(sg_dir))}"
                depth = f"{top:g}-{bottom:g}"
            else:
                for key, prop, sgk in (("SAND_R", "sand", "sand"), ("CLAY_R", "clay", "clay"),
                                       ("SILT_R", "silt", "silt"), ("OC_R", "oc_pct", "soc"),
                                       ("BULK_R", "bulk", "bdod"), ("CFRAG_R", "cfrag", "cfvo")):
                    path = rpath(f"{P}_{key}")
                    if path:
                        a = read(path)
                        arrays[prop] = a / SOILGRIDS_FACTORS[sgk] if soilgrids_units else a
                        used.append(os.path.basename(path))
                label = ((SOILGRIDS + " rasters: ") if soilgrids_units else "soil rasters: ") \
                    + ", ".join(used)
                depth = "as supplied"
            for need in ("sand", "clay"):
                if need not in arrays:
                    raise QgsProcessingException(f"Soil rasters: {need} is required.")
            if "oc_pct" not in arrays:
                feedback.pushWarning("Soil rasters: no organic carbon - K cannot be computed "
                                     "(texture and HSG proxy only).")
            sg = SoilGrid.from_texture(arrays["sand"], arrays["clay"], arrays.get("oc_pct"),
                                       silt=arrays.get("silt"), bulk=arrays.get("bulk"),
                                       cfrag=arrays.get("cfrag"),
                                       info={"soil_dataset": label, "soil_depth_cm": depth})
            feedback.pushInfo(f"Soil: {label} (depth {depth} cm)")

        hsg_r = rpath(P + "_HSG_R")
        k_r = rpath(P + "_K_R")
        if sg is None and (hsg_r or k_r):
            sg = SoilGrid((info.rows, info.cols), {"soil_dataset": "none (HSG / K only)",
                                                   "soil_depth_cm": ""})
        if sg is None:
            return None
        if "hsg" in extra:
            sg.override_hsg(extra["hsg"], extra["hsg_source"])
        if "k" in extra:
            sg.override_k(extra["k"], extra["k_source"])
        if hsg_r:
            codes, desc = warp_to_grid(hsg_r, info, resampling="near")
            h = hsg_from_codes(codes, "hysogs")
            sg.override_hsg(h, f"{HYSOGS if 'hysog' in desc.lower() else 'HSG raster'}: {desc}")
            feedback.pushInfo(f"Soil: hydrologic soil group from {desc} "
                              f"({100.0 * (h > 0).mean():.0f}% of the grid known)")
        if k_r:
            kk, desc = warp_to_grid(k_r, info)
            sg.override_k(kk, f"K raster {desc}")
            feedback.pushInfo(f"Soil: K from {desc}")
        return sg

    def _soil_from_polygons(self, parameters, context, info, feedback, P, src, unit_field,
                            csv_path, top, bottom):
        """Polygon soil sources (v0.11-0.14 paths + G2 field mapping). -> (SoilGrid, extra)."""
        from qgis.core import (QgsCoordinateReferenceSystem, QgsCoordinateTransform,
                               QgsProject, QgsRectangle, QgsGeometry, QgsFeatureRequest)
        from ..core.soils.sotwis import load_sotwis, load_attributes
        from ..core.soils.sources import SoilGrid, hsg_from_codes
        from ..core.geometry.rasterize import rasterize_polygons
        import numpy as np
        db = self.parameterAsFile(parameters, P + "_DB", context) if P + "_DB" in parameters else ""
        names = src.fields().names()
        lower = {n.lower(): n for n in names}

        def fld(key):
            name = P + "_F_" + key
            if name not in parameters or not parameters[name]:
                return None
            if hasattr(self, "parameterAsStrings"):          # QGIS >= 3.32
                v = self.parameterAsStrings(parameters, name, context)
            else:
                v = self.parameterAsFields(parameters, name, context)
            return v[0] if v else None
        mapped = {k: fld(k) for k in ("SAND", "SILT", "CLAY", "OC", "BULK", "CFRAG", "DRAIN",
                                      "HSG_FIELD", "K_FIELD")}
        user_map = all(mapped[k] for k in ("SAND", "CLAY", "OC"))

        dem_crs = QgsCoordinateReferenceSystem()
        dem_crs.createFromWkt(info.projection_wkt)
        tr = None
        if dem_crs.isValid() and src.sourceCrs() != dem_crs:
            tr = QgsCoordinateTransform(src.sourceCrs(), dem_crs, QgsProject.instance())
        gt = info.geotransform
        ext = QgsRectangle(gt[0], gt[3] + info.rows * gt[5], gt[0] + info.cols * gt[1], gt[3])

        records, polys, index_of = [], [], {}
        sotwis_fields = all(f in names for f in ("SDTO", "STPC", "CLPC", "TOTC"))
        plain = all(f in lower for f in ("sand", "silt", "clay", "oc"))
        if csv_path and not db:
            if unit_field not in names:
                raise QgsProcessingException(
                    f"The soil polygons have no field '{unit_field}' to join to the soil table.")
        elif db:
            if unit_field not in names:
                raise QgsProcessingException(
                    f"The soil polygons have no field '{unit_field}' to join to the database.")
        elif not (user_map or sotwis_fields or plain or mapped["HSG_FIELD"] or mapped["K_FIELD"]):
            raise QgsProcessingException(
                "Soil polygons need the SOTWIS database, a soil table CSV, SOTWIS fields "
                "(SDTO, STPC, CLPC, TOTC), fields sand, silt, clay, oc (percent), or the "
                "sand / clay / organic carbon fields chosen under the advanced parameters.")
        # Soil maps often carry invalid polygons (the SOTWIS Kenya shapefile
        # does). Processing's default check would abort; the cell-centre
        # even-odd rasteriser handles self-touching rings sensibly, so read
        # them unchecked and report how many there were.
        req = QgsFeatureRequest()
        try:
            from qgis.core import Qgis
            nocheck = Qgis.InvalidGeometryCheck.NoCheck
        except AttributeError:
            nocheck = QgsFeatureRequest.InvalidGeometryCheck.GeometryNoCheck
        req.setInvalidGeometryCheck(nocheck)
        try:                                   # QGIS >= 3.36
            from qgis.core import Qgis
            skip = Qgis.ProcessingFeatureSourceFlag.SkipGeometryValidityChecks
        except AttributeError:
            from qgis.core import QgsProcessingFeatureSource
            skip = QgsProcessingFeatureSource.Flag.FlagSkipGeometryValidityChecks
        n_invalid = 0
        hsg_val, k_val = {}, {}
        for feat in src.getFeatures(req, skip):
            g = QgsGeometry(feat.geometry())
            if g is None or g.isEmpty():
                continue
            if tr is not None:
                g.transform(tr)
            if not g.boundingBox().intersects(ext):
                continue
            if not g.isGeosValid():
                n_invalid += 1
            if db:
                key = feat[unit_field]
            elif csv_path:
                key = str(feat[unit_field]).strip()
            else:
                key = int(feat.id())
                rec = {"id": key}
                for n in names:
                    rec[n] = feat[n]
                records.append(rec)
            if key not in index_of:
                index_of[key] = len(index_of) + 1
                if mapped["HSG_FIELD"]:
                    v = feat[mapped["HSG_FIELD"]]
                    hsg_val[index_of[key]] = None if v is None or str(v) == "NULL" else str(v)
                if mapped["K_FIELD"]:
                    v = feat[mapped["K_FIELD"]]
                    try:
                        k_val[index_of[key]] = float(v)
                    except (TypeError, ValueError):
                        pass
            parts = g.asMultiPolygon() if g.isMultipart() else [g.asPolygon()]
            for part in parts:
                polys.append(([[(p.x(), p.y()) for p in ring] for ring in part], index_of[key]))
        if not polys:
            feedback.pushWarning("No soil polygon overlaps the DEM - soil fields will be empty.")
        if n_invalid:
            feedback.pushInfo(f"Soil: {n_invalid} polygon(s) over the DEM have invalid geometry; "
                              "rasterised by cell centre (even-odd), not repaired.")
        if db:
            units, sinfo = load_sotwis(db, top, bottom)
        elif csv_path:
            from ..core.soils.sotwis import load_csv_units
            try:
                units, sinfo = load_csv_units(csv_path, unit_field)
            except ValueError as e:
                raise QgsProcessingException(str(e))
        elif user_map:
            gkg = (self.parameterAsEnum(parameters, P + "_OC_UNITS", context) == 1
                   if P + "_OC_UNITS" in parameters else False)
            if not mapped["SILT"]:                # silt = 100 - sand - clay
                for r in records:
                    try:
                        r["__silt"] = max(100.0 - float(r[mapped["SAND"]]) - float(r[mapped["CLAY"]]), 0.0)
                    except (TypeError, ValueError):
                        r["__silt"] = None
            units, sinfo = load_attributes(
                records, sand=mapped["SAND"], silt=mapped["SILT"] or "__silt",
                clay=mapped["CLAY"], oc_pct=mapped["OC"], oc_is_gkg=gkg,
                bulk=mapped["BULK"], cfrag=mapped["CFRAG"], drain=mapped["DRAIN"],
                label=f"{src.sourceName()} polygon fields "
                      f"(sand {mapped['SAND']}, clay {mapped['CLAY']}, oc {mapped['OC']}"
                      f"{' g/kg' if gkg else ' %'})")
        elif sotwis_fields:
            units, sinfo = load_attributes(records, sand="SDTO", silt="STPC", clay="CLPC",
                                           oc_pct="TOTC", oc_is_gkg=True,
                                           bulk="BULK" if "BULK" in names else None,
                                           cfrag="CFRAG" if "CFRAG" in names else None,
                                           drain="DRAIN" if "DRAIN" in names else None,
                                           label="SOTWIS polygon attributes (dominant soil only)")
        elif plain:
            units, sinfo = load_attributes(records, sand=lower["sand"], silt=lower["silt"],
                                           clay=lower["clay"], oc_pct=lower["oc"],
                                           bulk=lower.get("bulk"), cfrag=lower.get("cfrag"),
                                           drain=lower.get("drain"))
        else:
            units, sinfo = {}, {"soil_dataset": f"{src.sourceName()} (HSG / K fields only)",
                                "soil_depth_cm": ""}
        feedback.setProgressText("Rasterising soil polygons")
        grid = rasterize_polygons(polys, gt, (info.rows, info.cols))
        idx = {v: k for k, v in index_of.items()}
        missing = sum(1 for k in index_of if k not in units)
        feedback.pushInfo(f"Soil: {len(index_of)} unit(s) over the DEM from {sinfo['soil_dataset']}, "
                          f"depth {sinfo['soil_depth_cm']} cm"
                          + (f"; {missing} unit(s) without data (e.g. water, towns)"
                             if missing and units else ""))
        sg = SoilGrid.from_units(grid, idx, units, sinfo)
        extra = {}
        n = max(index_of.values(), default=0) + 1
        if hsg_val:
            lut = np.zeros(n, dtype=np.int8)
            for i, v in hsg_val.items():
                lut[i] = hsg_from_codes(np.array([v], dtype=object), "letters")[0]
            extra["hsg"] = np.where(grid > 0, lut[np.clip(grid, 0, n - 1)], 0)
            extra["hsg_source"] = f"{src.sourceName()} field {mapped['HSG_FIELD']}"
        if k_val:
            lut = np.full(n, np.nan)
            for i, v in k_val.items():
                lut[i] = v
            extra["k"] = np.where(grid > 0, lut[np.clip(grid, 0, n - 1)], np.nan)
            extra["k_source"] = f"{src.sourceName()} field {mapped['K_FIELD']}"
        return sg, extra

    def add_runoff_parameters(self):
        """Land cover + lookups for the curve number / Rational C block (F1)."""
        from qgis.core import (QgsProcessingParameterRasterLayer, QgsProcessingParameterEnum,
                               QgsProcessingParameterFile)
        self.addParameter(QgsProcessingParameterRasterLayer(
            "LANDCOVER", "Land cover (ESA WorldCover classes; optional, for CN and land-cover "
            "shares)", optional=True))
        self.addParameter(QgsProcessingParameterEnum(
            "CN_CONDITION", "Hydrologic condition for the TR-55 curve numbers",
            options=["fair", "good", "poor"], defaultValue=0))
        self.addParameter(QgsProcessingParameterEnum(
            "CN_AMC", "Antecedent moisture condition of the exported CN",
            options=["II (average)", "III (wet)", "I (dry)"], defaultValue=0))
        self.addParameter(self._advanced(QgsProcessingParameterFile(
            "CN_CSV", "Curve number lookup CSV (class, A, B, C, D; replaces the TR-55 default)",
            optional=True, fileFilter="CSV (*.csv *.txt)")))
        self.addParameter(self._advanced(QgsProcessingParameterFile(
            "RC_CSV", "Rational C lookup CSV (class, A, B, C, D or class, C; none ships)",
            optional=True, fileFilter="CSV (*.csv *.txt)")))

    def load_runoff(self, parameters, context, info, soil, feedback):
        """-> core.runoff.curve_number.RunoffInputs or None (no land cover)."""
        from ..core.runoff.curve_number import RunoffInputs, read_lookup_csv, CONDITIONS, AMC
        from ..core.raster import warp_to_grid
        if "LANDCOVER" not in parameters or not parameters["LANDCOVER"] or \
                self.parameterAsRasterLayer(parameters, "LANDCOVER", context) is None:
            return None
        lc, desc = warp_to_grid(self.raster_path(parameters, "LANDCOVER", context), info,
                                resampling="near")
        cond = CONDITIONS[self.parameterAsEnum(parameters, "CN_CONDITION", context)] \
            if "CN_CONDITION" in parameters else "fair"
        amc = AMC[self.parameterAsEnum(parameters, "CN_AMC", context)] \
            if "CN_AMC" in parameters else "II"
        cn_path = self.parameterAsFile(parameters, "CN_CSV", context) if parameters.get("CN_CSV") else ""
        rc_path = self.parameterAsFile(parameters, "RC_CSV", context) if parameters.get("RC_CSV") else ""
        try:
            cn_lut = read_lookup_csv(cn_path, "CN") if cn_path else None
            rc_lut = read_lookup_csv(rc_path, "C") if rc_path else None
        except ValueError as e:
            raise QgsProcessingException(str(e))
        hsg = soil.hsg if soil is not None else None
        ro = RunoffInputs(lc, hsg=hsg, condition=cond, amc=amc, cn_lookup=cn_lut,
                          cn_lookup_path=cn_path or None, rc_lookup=rc_lut,
                          rc_lookup_path=rc_path or None, lc_dataset=desc)
        if hsg is None:
            feedback.pushWarning("Curve numbers need hydrologic soil groups: give a soil source "
                                 "(or an HSG raster). Land-cover shares only.")
        else:
            feedback.pushInfo(f"Curve numbers: {ro.cn_lookup_name}; exported at AMC {amc}"
                              + ("" if cn_path else " - PROXY lookup, check before design use"))
        if ro.note:
            feedback.pushWarning(ro.note)
        return ro

    def add_rainfall_parameters(self, with_relation=True):
        """Rainfall zones, mean annual rainfall and R from rainfall (F3)."""
        from qgis.core import (QgsProcessingParameterFeatureSource, QgsProcessingParameterField,
                               QgsProcessingParameterRasterLayer, QgsProcessingParameterString,
                               QgsProcessingParameterEnum, QgsProcessingParameterNumber,
                               QgsProcessing)
        self.addParameter(QgsProcessingParameterFeatureSource(
            "RAIN_ZONES", "Rainfall zones (polygons; optional)",
            [QgsProcessing.SourceType.TypeVectorPolygon], optional=True))
        self.addParameter(QgsProcessingParameterField(
            "RAIN_ZONE_FIELD", "Zone name field", parentLayerParameterName="RAIN_ZONES",
            optional=True))
        self.addParameter(QgsProcessingParameterRasterLayer(
            "RAIN_MAP", "Mean annual rainfall raster, mm/yr (optional, e.g. a CHIRPS climatology)",
            optional=True))
        self.addParameter(self._advanced(QgsProcessingParameterString(
            "RAIN_MAP_DATASET", "Rainfall dataset name for the record (e.g. CHIRPS v2.0 1991-2020)",
            optional=True)))
        if with_relation:
            self.addParameter(self._advanced(QgsProcessingParameterEnum(
                "RAIN_R_RELATION", "RUSLE R from rainfall (an ESTIMATE; prefer an erosivity raster)",
                options=["none", "Renard & Freimund 1994", "Lo et al. 1985"], defaultValue=0)))
        self.addParameter(self._advanced(QgsProcessingParameterNumber(
            "RAIN_ZONE_QA", "Warn when the dominant zone covers less than (%)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=80.0, minValue=0.0,
            maxValue=100.0)))

    def load_rainfall(self, parameters, context, info, feedback):
        """-> core.runoff.rainfall.RainfallInputs or None (no zones and no rainfall raster)."""
        from qgis.core import (QgsCoordinateReferenceSystem, QgsCoordinateTransform,
                               QgsProject, QgsRectangle, QgsGeometry)
        from ..core.runoff.rainfall import RainfallInputs, R_RELATION_KEYS
        from ..core.raster import warp_to_grid
        from ..core.geometry.rasterize import rasterize_polygons
        zone_grid, zone_names, zone_src = None, {}, ""
        src = self.parameterAsSource(parameters, "RAIN_ZONES", context) \
            if parameters.get("RAIN_ZONES") else None
        if src is not None:
            field = self.field_parameter(parameters, "RAIN_ZONE_FIELD", context) \
                if parameters.get("RAIN_ZONE_FIELD") else None
            if not field:
                raise QgsProcessingException("Rainfall zones: choose the zone name field.")
            dem_crs = QgsCoordinateReferenceSystem()
            dem_crs.createFromWkt(info.projection_wkt)
            tr = None
            if dem_crs.isValid() and src.sourceCrs() != dem_crs:
                tr = QgsCoordinateTransform(src.sourceCrs(), dem_crs, QgsProject.instance())
            gt = info.geotransform
            ext = QgsRectangle(gt[0], gt[3] + info.rows * gt[5], gt[0] + info.cols * gt[1], gt[3])
            index_of, polys, unnamed = {}, [], 0
            for feat in src.getFeatures():
                g = QgsGeometry(feat.geometry())
                if g is None or g.isEmpty():
                    continue
                if tr is not None:
                    g.transform(tr)
                if not g.boundingBox().intersects(ext):
                    continue
                v = feat[field]
                name = "" if v is None or str(v) == "NULL" else str(v).strip()
                if not name:
                    unnamed += 1
                    continue
                if name not in index_of:
                    index_of[name] = len(index_of) + 1
                parts = g.asMultiPolygon() if g.isMultipart() else [g.asPolygon()]
                for part in parts:
                    polys.append(([[(p.x(), p.y()) for p in ring] for ring in part],
                                  index_of[name]))
            if unnamed:
                feedback.pushWarning(f"Rainfall zones: {unnamed} polygon(s) without a zone name "
                                     "ignored.")
            if not polys:
                feedback.pushWarning("No rainfall zone polygon overlaps the DEM - zone fields "
                                     "will be empty.")
            # polygons are burnt in layer order: where zones overlap, the later one wins
            zone_grid = rasterize_polygons(polys, gt, (info.rows, info.cols))
            zone_names = {i: n for n, i in index_of.items()}
            zone_src = f"{src.sourceName()} field {field}"
            feedback.pushInfo(f"Rainfall zones: {len(index_of)} zone(s) over the DEM from "
                              f"{zone_src}")
        map_grid, map_ds = None, ""
        if parameters.get("RAIN_MAP") and \
                self.parameterAsRasterLayer(parameters, "RAIN_MAP", context) is not None:
            path = self.raster_path(parameters, "RAIN_MAP", context)
            map_grid, desc = warp_to_grid(path, info, resampling="bilinear")
            label = self.parameterAsString(parameters, "RAIN_MAP_DATASET", context) \
                if parameters.get("RAIN_MAP_DATASET") else ""
            map_ds = f"{label} - {desc}" if label else desc
        rel = R_RELATION_KEYS[self.parameterAsEnum(parameters, "RAIN_R_RELATION", context)] \
            if "RAIN_R_RELATION" in parameters else "none"
        qa = self.parameterAsDouble(parameters, "RAIN_ZONE_QA", context) \
            if parameters.get("RAIN_ZONE_QA") not in (None, "") else 80.0
        if zone_grid is None and map_grid is None:
            if rel != "none":
                feedback.pushWarning("An R-P relation was chosen but no rainfall raster was "
                                     "given; no rainfall block.")
            return None
        try:
            rain = RainfallInputs(zone_grid, zone_names, zone_src, map_grid, map_ds, rel, qa)
        except ValueError as e:
            raise QgsProcessingException(str(e))
        for w in rain.warnings:
            feedback.pushWarning(w)
        if rain.r_grid is not None:
            feedback.pushWarning(f"rusle_r is an {rain.r_method}")
        return rain

    @staticmethod
    def qgs_field(name, kind):
        """QgsField for text/int/real that works on QGIS 3.22 and 4."""
        from qgis.core import QgsField
        try:
            from qgis.PyQt.QtCore import QMetaType
            t = {"real": QMetaType.Type.Double, "int": QMetaType.Type.Int,
                 "text": QMetaType.Type.QString}[kind]
            return QgsField(name, t)
        except (TypeError, AttributeError):
            from qgis.PyQt.QtCore import QVariant
            t = {"real": QVariant.Type.Double, "int": QVariant.Type.Int,
                 "text": QVariant.Type.String}[kind]
            return QgsField(name, t)

    @staticmethod
    def write_vector(path, layer_name, projection_wkt, geom_kind, fields, rows):
        """Write (geometry, attrs) rows to a vector file with OGR (format from
        the extension: .shp, .geojson ... else GeoPackage).

        geom_kind: 'point' ((x, y)), 'line' ([(x, y), ...]) or
        'multipolygon' (core polygonize output). fields: [(name, type)] with
        type text/int/real (or float). None/NaN are written as NULL.
        """
        import math
        import os
        from osgeo import ogr, osr
        from ..core.interop.gpkg import wkb_point, wkb_linestring, wkb_multipolygon
        to_wkb = {"point": lambda g: wkb_point(*g), "line": wkb_linestring,
                  "multipolygon": wkb_multipolygon}[geom_kind]
        gtype = {"point": ogr.wkbPoint, "line": ogr.wkbLineString,
                 "multipolygon": ogr.wkbMultiPolygon}[geom_kind]
        otype = {"int": ogr.OFTInteger, "real": ogr.OFTReal, "float": ogr.OFTReal,
                 "text": ogr.OFTString}
        srs = None
        if projection_wkt:
            srs = osr.SpatialReference(); srs.ImportFromWkt(projection_wkt)
        from ..core.raster import ogr_driver_for
        drv = ogr_driver_for(path)
        if os.path.exists(path):
            drv.DeleteDataSource(path)
        ds = drv.CreateDataSource(path)
        layer = ds.CreateLayer(layer_name, srs=srs, geom_type=gtype)
        for fname, ftype in fields:
            layer.CreateField(ogr.FieldDefn(fname, otype[ftype]))
        defn = layer.GetLayerDefn()
        # by index: a shapefile truncates names to 10 characters
        for geom, attrs in rows:
            feat = ogr.Feature(defn)
            for i, (fname, ftype) in enumerate(fields):
                val = attrs.get(fname)
                if val is None or (ftype in ("int", "real", "float")
                                   and not math.isfinite(float(val))):
                    feat.SetFieldNull(i)
                elif ftype == "int":
                    feat.SetField(i, int(val))
                elif ftype in ("real", "float"):
                    feat.SetField(i, float(val))
                else:
                    feat.SetField(i, str(val))
            feat.SetGeometry(ogr.CreateGeometryFromWkb(to_wkb(geom)))
            layer.CreateFeature(feat)
            feat = None
        ds = None
        return path

    def report_stats(self, feedback, title, stats):
        """Print a stats block. Handles str/int/float without assuming type -
        a string value here previously raised
        'ValueError: Cannot specify comma with s' and killed the algorithm
        AFTER all the real work was done."""
        feedback.pushInfo(f"--- {title} ---")
        for key, value in stats.items():
            if isinstance(value, bool):
                text = str(value)
            elif isinstance(value, float):
                text = f"{value:,.6g}"
            elif isinstance(value, int):
                text = f"{value:,}"
            else:
                text = str(value)
            feedback.pushInfo(f"    {key}: {text}")
