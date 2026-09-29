# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Renumber and relink an edited exchange package (decision D3)."""

import os

from qgis.core import (
    QgsProcessingParameterFile, QgsProcessingParameterRasterLayer,
    QgsProcessingParameterFeatureSource, QgsProcessingParameterNumber,
    QgsProcessingParameterBoolean, QgsProcessingParameterString,
    QgsProcessingParameterFileDestination, QgsProcessing, QgsProcessingException,
    QgsProcessingContext,
)
import numpy as np

from .base import QehtAlgorithm
from ..core.raster import read_dem
from ..core.grid import decode_d8
from ..core.watershed.delineate import extract_streams
from ..core.interop.heas_exchange import (build_exchange_records, write_exchange,
                                          validate_exchange, crs_check,
                                          file_fingerprint, plugin_version, ExchangeError)
from ..core.interop import gpkg
from ..core.linking.relink import renumber_log, LOG_FIELDS

PACKAGE = "PACKAGE"; FDR = "FDR"; FAC = "FAC"; RAW_DEM = "RAW_DEM"; ORDER = "ORDER"
ROAD = "ROAD"; START = "START"; REVERSE = "REVERSE"; PREFIX = "PREFIX"
SNAP = "SNAP"; SNAP_THRESHOLD = "SNAP_THRESHOLD"; LOCAL = "LOCAL"; OUTPUT = "OUTPUT"


class RenumberRelinkAlgorithm(QehtAlgorithm):

    def name(self): return "renumberrelink"
    def displayName(self): return "Renumber and relink exchange package"
    def group(self): return "Interoperability"
    def groupId(self): return "interop"

    def shortHelpString(self):
        return (
            "After editing the <i>crossings</i> layer of an exchange package in QGIS "
            "- deleting, moving or adding points - run this to bring the package "
            "back into a consistent state:\n"
            "• crossings are re-sorted (by chainage when every crossing has one, "
            "else downstream-first) and given gapless outlet_uid values;\n"
            "• catchments, flow paths and every attribute are recomputed;\n"
            "• a renumber_log table records old -> new IDs (unchanged, renumbered, "
            "new, deleted, and how far a moved point moved), so references "
            "elsewhere can be updated.\n\n"
            "Unmoved crossings keep their outlet cell exactly; moved and added "
            "points are snapped to the nearest stream cell. Give the road "
            "alignment to (re)compute chainages of moved and added points. IDs "
            "are only ever re-issued when you run this tool. A NEW package is "
            "written; the edited one is left as it is.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile(
            PACKAGE, "Edited exchange package", extension="gpkg"))
        self.addParameter(QgsProcessingParameterRasterLayer(FDR, "Flow direction (D8-coded)"))
        self.addParameter(QgsProcessingParameterRasterLayer(FAC, "Flow accumulation"))
        self.addParameter(QgsProcessingParameterRasterLayer(
            RAW_DEM, "Raw DEM (for reported elevations and slopes)"))
        self.addParameter(QgsProcessingParameterRasterLayer(
            ORDER, "Stream order raster (optional)", optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource(
            ROAD, "Road alignment (optional; chainage of moved/added points)",
            [QgsProcessing.SourceType.TypeVectorLine], optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            START, "Start chainage (m)", QgsProcessingParameterNumber.Type.Double, defaultValue=0.0))
        self.addParameter(QgsProcessingParameterBoolean(
            REVERSE, "Reverse chainage direction", defaultValue=False))
        self.addParameter(QgsProcessingParameterString(
            PREFIX, "ID prefix (blank = keep the package's prefix)", optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            SNAP, "Snap radius for moved/added points (cells)",
            QgsProcessingParameterNumber.Type.Integer, defaultValue=5, minValue=0, maxValue=100))
        self.addParameter(QgsProcessingParameterNumber(
            SNAP_THRESHOLD, "Snap-to-stream threshold (cells)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=200.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterBoolean(
            LOCAL, "Non-overlapping (local) catchments", defaultValue=False))
        self.addParameter(QgsProcessingParameterFileDestination(
            OUTPUT, "Relinked exchange package", fileFilter="GeoPackage (*.gpkg)"))

    def processAlgorithm(self, parameters, context, feedback):
        from osgeo import ogr
        pkg = self.parameterAsFile(parameters, PACKAGE, context)
        out_path = self.parameterAsFileOutput(parameters, OUTPUT, context)
        if not out_path.lower().endswith(".gpkg"):
            out_path += ".gpkg"
        if os.path.abspath(out_path) == os.path.abspath(pkg):
            raise QgsProcessingException("Write the relinked package to a new file.")

        old_md = {r["key"]: r["value"] for r in gpkg.read_table(pkg, "qeht_run_metadata")}
        old_uids = [r.get("outlet_uid") for r in gpkg.read_table(pkg, "catchments",
                                                                 with_geometry=False)]
        d8, valid, info = read_dem(self.raster_path(parameters, FDR, context))
        ok, epsg, msg = crs_check(info.projection_wkt)
        if not ok:
            raise QgsProcessingException(msg)

        ds = ogr.Open(pkg)
        layer = ds.GetLayerByName("crossings") if ds is not None else None
        if layer is None:
            raise QgsProcessingException("The package has no 'crossings' layer.")
        outlets, old_by_fid = [], {}
        for f in layer:
            g = f.GetGeometryRef()
            if g is None:
                continue
            x, y = g.GetX(), g.GetY()
            fid = int(f.GetFID())
            names = [f.GetFieldDefnRef(k).GetName() for k in range(f.GetFieldCount())]
            get = lambda n: (f.GetField(n) if n in names and f.IsFieldSetAndNotNull(n) else None)
            ox, oy, uid = get("outlet_x"), get("outlet_y"), get("outlet_uid")
            unmoved = (ox is not None and oy is not None
                       and abs(ox - x) < 1e-6 and abs(oy - y) < 1e-6)
            o = {"x": x, "y": y, "fid": fid, "source_id": None,
                 "chainage": get("chainage_m") if unmoved else None}
            if unmoved:
                o["outlet_x"], o["outlet_y"] = ox, oy
            outlets.append(o)
            old_by_fid[fid] = {"uid": uid, "x": ox, "y": oy}
        ds = None
        if not outlets:
            raise QgsProcessingException("The crossings layer is empty.")

        alignment = self.read_alignment(
            parameters, ROAD, context, info, feedback,
            start_chainage=self.parameterAsDouble(parameters, START, context),
            reverse=self.parameterAsBool(parameters, REVERSE, context))
        if alignment is not None:
            miss = [o for o in outlets if o["chainage"] is None]
            if miss:
                ch, _, _ = alignment.locate([o["x"] for o in miss], [o["y"] for o in miss])
                for o, c in zip(miss, ch):
                    o["chainage"] = float(c)
        order = "chainage" if all(o["chainage"] is not None for o in outlets) else "downstream"
        if order == "downstream" and any(o["chainage"] is not None for o in outlets):
            feedback.pushWarning("Some crossings have no chainage (moved/added without a road "
                                 "alignment): numbering downstream-first instead.")

        direction = decode_d8(d8.astype(np.int32))
        accum, _, _ = read_dem(self.raster_path(parameters, FAC, context))
        raw_path = self.raster_path(parameters, RAW_DEM, context)
        elevation, raw_valid, _ = read_dem(raw_path)
        elevation = np.where(raw_valid, elevation, np.nan)
        stream_order = None
        if self.parameterAsRasterLayer(parameters, ORDER, context) is not None:
            so, sov, _ = read_dem(self.raster_path(parameters, ORDER, context))
            stream_order = np.where(sov, so, np.nan)
        thr = self.parameterAsDouble(parameters, SNAP_THRESHOLD, context)
        snap = self.parameterAsInt(parameters, SNAP, context)
        stream_mask = extract_streams(accum, valid, threshold_cells=thr) if thr > 0 else None
        prefix = (self.parameterAsString(parameters, PREFIX, context) or "").strip() \
            or old_md.get("id_prefix", "X")
        local = self.parameterAsBool(parameters, LOCAL, context)
        try:
            cr, ca, fp, issues, id_info = build_exchange_records(
                direction, valid, accum, elevation, info.geotransform, outlets,
                snap_radius_cells=snap, stream_mask=stream_mask, local=local,
                stream_order=stream_order, id_prefix=prefix, id_order=order,
                progress=self.make_progress(feedback, weight=0.9))
        except ExchangeError as e:
            raise QgsProcessingException(str(e))
        for m in issues:
            feedback.pushWarning(m)

        log = renumber_log(cr, old_by_fid, old_uids)
        md = {k: v for k, v in old_md.items()
              if k not in ("schema_version", "qeht_version", "run_utc", "n_crossings")}
        md.update(id_info)
        md.update({"qeht_version": plugin_version(), "relinked_from": pkg,
                   "crossing_source": "relinked package",
                   "raw_dem_path": raw_path, "dem_sha256": file_fingerprint(raw_path),
                   "catchment_mode": "local" if local else "full",
                   "parameters_json": {k: str(v) for k, v in parameters.items()}})
        try:
            write_exchange(out_path, cr, ca, fp, md, info.projection_wkt, crs_epsg=epsg,
                           extra_tables=[("renumber_log", LOG_FIELDS, log,
                                          "old -> new outlet_uid (renumber and relink)")])
        except ExchangeError as e:
            raise QgsProcessingException(str(e))
        errors, warnings = validate_exchange(out_path)
        for w in warnings:
            feedback.pushWarning(w)
        if errors:
            raise QgsProcessingException("The package failed validation: " + "; ".join(errors))
        for r in log:
            feedback.pushInfo(f"  {r['old_uid'] or '-':>8} -> {r['new_uid'] or '-':<8} {r['change']}"
                              + (f"  moved {r['moved_m']:.1f} m" if r.get("moved_m") else ""))
        stem = os.path.splitext(os.path.basename(out_path))[0]
        for lname in ("catchments", "flowpaths", "crossings"):
            try:
                context.addLayerToLoadOnCompletion(
                    f"{out_path}|layername={lname}",
                    QgsProcessingContext.LayerDetails(f"{stem} - {lname}", context.project(), lname))
            except Exception as e:
                feedback.pushInfo(f"(Could not auto-load {lname}: {e})")
        return {OUTPUT: out_path}
