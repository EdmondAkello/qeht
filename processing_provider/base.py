# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Shared base for QEHT Processing algorithms.

This is the ONLY place where QGIS and `qeht.core` meet. The core never
imports QGIS; the algorithms below never implement hydrology. If you
find yourself writing a loop over cells in this package, it belongs in
core instead.
"""

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
                if f in names:
                    v = feature[f]
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

    def add_soil_parameters(self, prefix="SOIL", optional=True):
        """Standard soil inputs: polygons, optional SOTWIS database, depth."""
        from qgis.core import (QgsProcessingParameterFeatureSource, QgsProcessing,
                               QgsProcessingParameterFile, QgsProcessingParameterString,
                               QgsProcessingParameterNumber)
        self.addParameter(QgsProcessingParameterFeatureSource(
            prefix + "_POLYGONS", "Soil map polygons (e.g. SOTWIS KEN_SOTWISv1_t1s1d1)",
            [QgsProcessing.SourceType.TypeVectorPolygon], optional=optional))
        self.addParameter(QgsProcessingParameterFile(
            prefix + "_DB", "SOTWIS SQLite database (optional; full unit composition)",
            optional=True, fileFilter="SQLite (*.db *.sqlite)"))
        self.addParameter(QgsProcessingParameterFile(
            prefix + "_CSV", "Soil table CSV (optional; unit field + sand, silt, clay, oc %, "
            "[bulk, cfrag, drain])", optional=True, fileFilter="CSV (*.csv *.txt)"))
        self.addParameter(QgsProcessingParameterString(
            prefix + "_UNIT_FIELD", "Soil unit field linking polygons to the database or table",
            defaultValue="NEWSUID", optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            prefix + "_TOP", "Soil depth from (cm)", QgsProcessingParameterNumber.Type.Double,
            defaultValue=0.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber(
            prefix + "_BOTTOM", "Soil depth to (cm)", QgsProcessingParameterNumber.Type.Double,
            defaultValue=20.0, minValue=1.0))

    def load_soil(self, parameters, context, info, feedback, prefix="SOIL"):
        """Soil inputs -> (unit_grid, index_to_unit, units, info) on the DEM grid,
        or None when no soil polygons were given.

        With the SOTWIS database: full component composition per unit
        (SOTERunitComposition x SOTERparameterEstimates), depth-weighted.
        Without it: polygon attributes - SOTWIS field names (SDTO, STPC,
        CLPC, TOTC g/kg, BULK, CFRAG, DRAIN; dominant soil only) or plain
        names (sand, silt, clay, oc [%], bulk, cfrag, drain).
        """
        from qgis.core import (QgsCoordinateReferenceSystem, QgsCoordinateTransform,
                               QgsProject, QgsRectangle, QgsGeometry)
        from ..core.soils.sotwis import load_sotwis, load_attributes
        from ..core.geometry.rasterize import rasterize_polygons
        import numpy as np
        src = self.parameterAsSource(parameters, prefix + "_POLYGONS", context)
        if src is None:
            return None
        db = self.parameterAsFile(parameters, prefix + "_DB", context)
        csv_path = self.parameterAsFile(parameters, prefix + "_CSV", context) \
            if (prefix + "_CSV") in parameters else ""
        unit_field = (self.parameterAsString(parameters, prefix + "_UNIT_FIELD", context)
                      or "NEWSUID").strip()
        top = self.parameterAsDouble(parameters, prefix + "_TOP", context)
        bottom = self.parameterAsDouble(parameters, prefix + "_BOTTOM", context)
        names = src.fields().names()
        lower = {n.lower(): n for n in names}

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
        elif not (sotwis_fields or plain):
            raise QgsProcessingException(
                "Soil polygons need either the SOTWIS database, SOTWIS fields "
                "(SDTO, STPC, CLPC, TOTC) or fields sand, silt, clay, oc (percent).")
        # Soil maps often carry invalid polygons (the SOTWIS Kenya shapefile
        # does). Processing's default check would abort; the cell-centre
        # even-odd rasteriser handles self-touching rings sensibly, so read
        # them unchecked and report how many there were.
        from qgis.core import QgsFeatureRequest
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
        elif sotwis_fields:
            units, sinfo = load_attributes(records, sand="SDTO", silt="STPC", clay="CLPC",
                                           oc_pct="TOTC", oc_is_gkg=True,
                                           bulk="BULK" if "BULK" in names else None,
                                           cfrag="CFRAG" if "CFRAG" in names else None,
                                           drain="DRAIN" if "DRAIN" in names else None,
                                           label="SOTWIS polygon attributes (dominant soil only)")
        else:
            units, sinfo = load_attributes(records, sand=lower["sand"], silt=lower["silt"],
                                           clay=lower["clay"], oc_pct=lower["oc"],
                                           bulk=lower.get("bulk"), cfrag=lower.get("cfrag"),
                                           drain=lower.get("drain"))
        feedback.setProgressText("Rasterising soil polygons")
        grid = rasterize_polygons(polys, gt, (info.rows, info.cols))
        idx = {v: k for k, v in index_of.items()}
        missing = sum(1 for k in index_of if k not in units)
        feedback.pushInfo(f"Soil: {len(index_of)} unit(s) over the DEM from {sinfo['soil_dataset']}, "
                          f"depth {sinfo['soil_depth_cm']} cm"
                          + (f"; {missing} unit(s) without data (e.g. water, towns)" if missing else ""))
        return grid, idx, units, sinfo

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
        """Write (geometry, attrs) rows to a GeoPackage layer with OGR.

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
        drv = ogr.GetDriverByName("GPKG")
        if os.path.exists(path):
            drv.DeleteDataSource(path)
        ds = drv.CreateDataSource(path)
        layer = ds.CreateLayer(layer_name, srs=srs, geom_type=gtype)
        for fname, ftype in fields:
            layer.CreateField(ogr.FieldDefn(fname, otype[ftype]))
        defn = layer.GetLayerDefn()
        for geom, attrs in rows:
            feat = ogr.Feature(defn)
            for fname, ftype in fields:
                val = attrs.get(fname)
                if val is None or (ftype in ("int", "real", "float")
                                   and not math.isfinite(float(val))):
                    feat.SetFieldNull(fname)
                elif ftype == "int":
                    feat.SetField(fname, int(val))
                elif ftype in ("real", "float"):
                    feat.SetField(fname, float(val))
                else:
                    feat.SetField(fname, str(val))
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
