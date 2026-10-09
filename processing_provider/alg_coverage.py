# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Drainage coverage check along a road (A2/A3, v0.16)."""

import numpy as np
from qgis.core import (
    QgsProcessingParameterRasterLayer, QgsProcessingParameterFeatureSource,
    QgsProcessingParameterNumber, QgsProcessingParameterBoolean,
    QgsProcessingParameterVectorDestination, QgsProcessing, QgsProcessingException,
)

from .base import QehtAlgorithm
from ..core.raster import read_dem
from ..core.grid import decode_d8
from ..core.network.profile import profile_with_crossings
from ..core.network.coverage import (cluster_cover, missing_crossings, ponding_sags, flat_stretches,
                                     walled_accumulation, sag_areas, COVERAGE_FIELDS,
                                     SAG_FIELDS, FLAT_FIELDS, _nearest)
from ..core.network.crossings import CANDIDATE_FIELDS


class DrainageCoverageAlgorithm(QehtAlgorithm):

    def name(self): return "drainagecoverage"
    def displayName(self): return "Drainage coverage check along a road"
    def group(self): return "Road drainage"
    def groupId(self): return "roaddrainage"

    def shortHelpString(self):
        return (
            "Checks a road alignment for drainage that has no crossing, from the DEM alone.\n\n"
            "<b>Missing crossings:</b> places where a stream of at least the given area "
            "crosses the centreline (exact D8 link intersections) with no existing "
            "crossing within the search distance of chainage.\n"
            "<b>Sag points:</b> low points of the ground profile (smoothed; depth = how far "
            "the ground rises on both sides) where water ponds against the embankment, "
            "with the local area draining to each one measured with the alignment as a "
            "wall - the places for equalisers / relief culverts.\n"
            "<b>Flat stretches:</b> runs flatter than the slope limit along and across the "
            "road, longer than the minimum length, where relief culverts for sheet flow "
            "may be needed (common practice: single pipes at a maximum spacing, e.g. "
            "200 m). QEHT does not place or size them.\n\n"
            "Existing crossings: pour points or a layer from 'Road crossing candidates' "
            "(accepted / recommended ones). Without them every stream crossing is "
            "reported. 'Build design hydrology package' runs the same check and also "
            "delineates and characterises each missing crossing as a proposed one.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer("RAW_DEM", "Raw DEM"))
        self.addParameter(QgsProcessingParameterRasterLayer(
            "FILLED", "Filled DEM (optional; filled from the raw DEM if not given)", optional=True))
        self.addParameter(QgsProcessingParameterRasterLayer("FDR", "Flow direction (D8-coded)"))
        self.addParameter(QgsProcessingParameterRasterLayer("FAC", "Flow accumulation"))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "ROAD", "Road centreline", [QgsProcessing.SourceType.TypeVectorLine]))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "CROSSINGS", "Existing crossings (pour points or crossing candidates; optional)",
            [QgsProcessing.SourceType.TypeVectorPoint], optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            "START", "Start chainage (m)", QgsProcessingParameterNumber.Type.Double, defaultValue=0.0))
        self.addParameter(QgsProcessingParameterBoolean(
            "REVERSE", "Reverse chainage direction", defaultValue=False))
        for key, label, default in (
                ("THRESHOLD", "Stream threshold (cells)", 200.0),
                ("MIN_AREA", "Smallest stream area to need a crossing (km2; 0 = stream threshold)", 0.0),
                ("SEARCH", "An existing crossing within this chainage covers a stream (m)", 50.0),
                ("MERGE", "Merge stream stations within (m)", 30.0),
                ("STEP", "Profile station spacing (m)", 10.0),
                ("SAG_SMOOTH", "Sags: smoothing window (m)", 30.0),
                ("SAG_DEPTH", "Sags: minimum depth (m)", 0.3),
                ("FLAT_SLOPE", "Flat stretches: slope below (%)", 0.5),
                ("FLAT_CROSSFALL", "Flat stretches: cross-fall measured over (m)", 100.0),
                ("FLAT_MIN_LEN", "Flat stretches: minimum length (m)", 300.0)):
            self.addParameter(QgsProcessingParameterNumber(
                key, label, QgsProcessingParameterNumber.Type.Double, defaultValue=default,
                minValue=0.0))
        self.addParameter(QgsProcessingParameterVectorDestination(
            "COVERAGE_OUT", "Coverage findings", QgsProcessing.SourceType.TypeVectorPoint))
        self.addParameter(QgsProcessingParameterVectorDestination(
            "SAGS_OUT", "Sag points", QgsProcessing.SourceType.TypeVectorPoint))
        self.addParameter(QgsProcessingParameterVectorDestination(
            "FLATS_OUT", "Flat stretches", QgsProcessing.SourceType.TypeVectorLine))

    def processAlgorithm(self, parameters, context, feedback):
        P = lambda k: self.parameterAsDouble(parameters, k, context)
        raw, rv, info = read_dem(self.raster_path(parameters, "RAW_DEM", context))
        raw = np.where(rv, raw, np.nan)
        d8, dv, di = read_dem(self.raster_path(parameters, "FDR", context))
        acc, av, ai = read_dem(self.raster_path(parameters, "FAC", context))
        if (di.rows, di.cols) != (info.rows, info.cols) or (ai.rows, ai.cols) != (info.rows, info.cols):
            raise QgsProcessingException("Flow direction, accumulation and the raw DEM must be on one grid.")
        direction = decode_d8(d8.astype(np.int32))
        valid = dv & rv
        if self.parameterAsRasterLayer(parameters, "FILLED", context) is not None:
            fz, fv, _ = read_dem(self.raster_path(parameters, "FILLED", context))
            filled = np.where(fv, fz, np.nan)
        else:
            from ..core.conditioning.fill import fill_depressions
            fz, _, _ = fill_depressions(np.nan_to_num(raw, nan=0.0), valid,
                                        cell_width=info.cell_width, cell_height=info.cell_height)
            filled = np.where(rv, fz, np.nan)
        al = self.read_alignment(parameters, "ROAD", context, info, feedback,
                                 start_chainage=P("START"),
                                 reverse=self.parameterAsBool(parameters, "REVERSE", context))
        thr = P("THRESHOLD")
        prof, _ = profile_with_crossings(al, raw, info.geotransform, direction=direction,
                                         valid=valid, accumulation=acc, filled=filled,
                                         stream_threshold_cells=thr, step=P("STEP"))
        existing, pts, existing_cov = [], [], []
        if parameters.get("CROSSINGS"):
            pts = self.read_pour_points(parameters, "CROSSINGS", context, info, feedback,
                                        extra_fields=[f for f, _ in CANDIDATE_FIELDS])
            pts, cand_mode, all_pts = self.candidate_selection(pts, feedback, info.cell_width)
            ch, _, _ = al.locate([p.get("outlet_x") or p["x"] for p in pts],
                                 [p.get("outlet_y") or p["y"] for p in pts])
            existing = [(f"existing {k + 1}" if p.get("source_id") is None else str(p["source_id"]),
                         float(c)) for k, (p, c) in enumerate(zip(pts, ch))]
            if cand_mode:
                # the other intersections of a clustered stream are covered by its crossing
                existing_names = [u for u, _ in existing]
                existing_cov = cluster_cover(
                    [(p.get("attr_cluster_id"), p.get("attr_chainage_m")) for p in all_pts],
                    [(p.get("attr_cluster_id"), u) for p, u in zip(pts, existing_names)])
            else:
                existing_cov = []
        cell_km2 = info.cell_width * info.cell_height / 1e6
        min_a = P("MIN_AREA") or (thr + 1.0) * cell_km2
        miss = missing_crossings(prof, existing + existing_cov, min_a, P("MERGE"), P("SEARCH"))
        sags = ponding_sags(prof, existing + existing_cov, P("SEARCH"), P("SAG_SMOOTH"), P("SAG_DEPTH"))
        if sags:
            feedback.setProgressText("Routing with the alignment as a wall (sag areas)")
            gaps = [(p_.get("outlet_x") or p_["x"], p_.get("outlet_y") or p_["y"])
                    for p_ in (pts if parameters.get("CROSSINGS") else [])] + \
                [(r["x"], r["y"]) for r in prof if r.get("stream") == 1]
            wacc, wall, wdir = walled_accumulation(filled, valid, al, info.geotransform,
                                                   info.cell_width, info.cell_height, gaps=gaps)
            sag_areas(sags, wacc, wall, info.geotransform, al, direction=wdir)
        flats = flat_stretches(prof, raw, info.geotransform, al, P("FLAT_SLOPE"),
                               P("FLAT_CROSSFALL"), P("FLAT_MIN_LEN"), existing)
        rows = []
        for m in miss:
            if not m["covered"]:
                rows.append(((m["x"], m["y"]), {
                    "issue": "missing_crossing", "chainage_m": m["chainage_m"],
                    "area_km2": m["area_km2"], "nearest_uid": m["nearest_uid"],
                    "nearest_m": m["nearest_m"],
                    "note": f"stream of {m['area_km2']:.3g} km2 with no crossing within {P('SEARCH'):g} m"}))
        for s in sags:
            u, d = _nearest(s["chainage_m"], existing)
            rows.append(((s["x"], s["y"]), {
                "issue": "sag_point", "chainage_m": s["chainage_m"], "area_km2": s.get("sag_area_km2"),
                "nearest_uid": u, "nearest_m": d,
                "note": (f"low point {s['sag_depth_m']:.2f} m deep; ponds against the embankment "
                         f"({s['side']} side, {s['sag_area_km2']:.3g} km2)" if s.get("sag_area_km2")
                         else f"low point {s['sag_depth_m']:.2f} m deep; no water ponds here")}))
        for f in flats:
            rows.append((al.point_at(0.5 * (f["chainage_m"] + f["chainage_to_m"])), {
                "issue": "flat_stretch", "chainage_m": f["chainage_m"],
                "chainage_to_m": f["chainage_to_m"],
                "note": f"{f['length_m']:.0f} m flat, {f['n_crossings_within']} crossing(s) within"}))
        rows.sort(key=lambda t: t[1]["chainage_m"])
        out = {}
        out["COVERAGE_OUT"] = self.parameterAsOutputLayer(parameters, "COVERAGE_OUT", context)
        self.write_vector(out["COVERAGE_OUT"], "coverage_check", info.projection_wkt, "point",
                          COVERAGE_FIELDS, rows)
        out["SAGS_OUT"] = self.parameterAsOutputLayer(parameters, "SAGS_OUT", context)
        self.write_vector(out["SAGS_OUT"], "sag_points", info.projection_wkt, "point", SAG_FIELDS,
                          [((s["x"], s["y"]), s) for s in sags])
        out["FLATS_OUT"] = self.parameterAsOutputLayer(parameters, "FLATS_OUT", context)
        lines = []
        for f in flats:
            cs = [f["chainage_m"]] + [r["chainage_m"] for r in prof
                                      if f["chainage_m"] < r["chainage_m"] < f["chainage_to_m"]] \
                + [f["chainage_to_m"]]
            lines.append(([al.point_at(c) for c in cs], f))
        self.write_vector(out["FLATS_OUT"], "flat_stretches", info.projection_wkt, "line",
                          FLAT_FIELDS, lines)
        n_m = sum(1 for m in miss if not m["covered"])
        feedback.pushInfo(f"{n_m} stream crossing(s) without a culvert, {len(sags)} sag point(s), "
                          f"{len(flats)} flat stretch(es) along {al.length:,.0f} m.")
        for _, a in rows:
            feedback.pushInfo(f"  ch {a['chainage_m']:,.0f}: {a['issue']} - {a['note']}")
        return out
