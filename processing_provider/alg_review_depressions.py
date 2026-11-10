# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Review fill depressions (1.2): one polygon per filled depression, with
depth, area, volume and, given a road, distance, chainage and a flag for
deep fills at or near the road."""

from qgis.core import (
    QgsProcessingParameterRasterLayer, QgsProcessingParameterFeatureSource,
    QgsProcessingParameterNumber, QgsProcessingParameterBoolean, QgsProcessing,
    QgsProcessingParameterVectorDestination, QgsProcessingParameterRasterDestination,
    QgsProcessingException,
)

from .base import QehtAlgorithm

FIELD_TYPES = {"int": "int", "float": "real", "str": "text"}


def add_review_parameters(alg, prefix=""):
    """Thresholds shared by the depression and package review tools."""
    alg.addParameter(QgsProcessingParameterNumber(
        prefix + "MIN_DEPTH", "List depressions filled at least this deep (m)",
        QgsProcessingParameterNumber.Type.Double, defaultValue=0.5, minValue=0.0))
    alg.addParameter(QgsProcessingParameterNumber(
        prefix + "DEEP", "Flag fills at least this deep near the road (m)",
        QgsProcessingParameterNumber.Type.Double, defaultValue=1.0, minValue=0.0))
    alg.addParameter(QgsProcessingParameterNumber(
        prefix + "HALF_WIDTH", "Road half-width: a fill this close to the centreline touches "
        "the road (m)", QgsProcessingParameterNumber.Type.Double, defaultValue=15.0, minValue=0.0))
    alg.addParameter(QgsProcessingParameterNumber(
        prefix + "NEAR", "Near the road / near a crossing: within this distance (m)",
        QgsProcessingParameterNumber.Type.Double, defaultValue=100.0, minValue=0.0))


def review_settings(alg, parameters, context, prefix=""):
    d = lambda k: alg.parameterAsDouble(parameters, prefix + k, context)   # noqa: E731
    return {"min_depth": d("MIN_DEPTH"), "deep_m": d("DEEP"),
            "road_half_width_m": d("HALF_WIDTH"), "near_m": d("NEAR")}


def write_depressions(alg, path, rows, projection_wkt):
    from ..core.review.depressions import DEPRESSION_FIELDS
    fields = [(n, FIELD_TYPES[t]) for n, t in DEPRESSION_FIELDS]
    return alg.write_vector(path, "fill_depressions", projection_wkt, "multipolygon", fields,
                            [(r["_geom"], r) for r in rows if r.get("_geom")])


class ReviewDepressionsAlgorithm(QehtAlgorithm):

    def name(self): return "reviewdepressions"
    def displayName(self): return "Review fill depressions"
    def group(self): return "Review"
    def groupId(self): return "review"

    def shortHelpString(self):
        return (
            "Lists the depressions that filling raised, one polygon each, deepest first: "
            "cells, area, maximum and mean fill depth, volume and the deepest point.\n\n"
            "With a <b>road</b>, each depression also gets its distance from the centreline, "
            "the chainage of its deepest point, whether it touches the road (within the "
            "half-width) and, with <b>crossings</b>, the nearest crossing along the road. "
            "Deep fills against or near the road are flagged: the DEM may be treating the "
            "embankment as a dam, so check for a culvert the DEM does not show (and burn it "
            "in) or a missing one.\n\n"
            "The thresholds are review settings, not design criteria. On a DEM stored in whole "
            "metres many 1 m pits are noise; raise the depths to suit.\n\n"
            "Give the DEM the fill ran on (after any breaching) and the filled DEM; without a "
            "filled DEM the tool fills a copy itself (no minimum slope)."
        )

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer("DEM", "DEM before filling"))
        self.addParameter(QgsProcessingParameterRasterLayer(
            "FILLED", "Filled DEM (blank = fill here)", optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "ROAD", "Road alignment (optional)", [QgsProcessing.SourceType.TypeVectorLine],
            optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            "START", "Chainage at the start of the road (m)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=0.0))
        self.addParameter(QgsProcessingParameterBoolean(
            "REVERSE", "Chainage runs from the end of the line", defaultValue=False))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "CROSSINGS", "Crossings (optional; located on the road for chainage)",
            [QgsProcessing.SourceType.TypeVectorPoint], optional=True))
        self.addParameter(QgsProcessingParameterVectorDestination(
            "OUTPUT", "Fill depressions", type=QgsProcessing.SourceType.TypeVectorPolygon))
        add_review_parameters(self)
        self.addParameter(QgsProcessingParameterRasterDestination(
            "DEPTH", "Fill depth", optional=True, createByDefault=False))

    def processAlgorithm(self, parameters, context, feedback):
        from ..core.raster import read_dem, write_raster
        from ..core.conditioning.fill import fill_depressions
        from ..core.review.depressions import depression_table
        dem, valid, info = read_dem(self.raster_path(parameters, "DEM", context))
        self.check_size(feedback, info.rows, info.cols)
        if self.parameterAsRasterLayer(parameters, "FILLED", context) is not None:
            filled, fvalid, finfo = read_dem(self.raster_path(parameters, "FILLED", context))
            if (finfo.rows, finfo.cols) != (info.rows, info.cols) or any(
                    abs(a - b) > 1e-6 * max(1.0, abs(a))
                    for a, b in zip(finfo.geotransform, info.geotransform)):
                raise QgsProcessingException("The filled DEM is not on the DEM's grid.")
            valid = valid & fvalid
        else:
            feedback.pushInfo("Filling a copy of the DEM (priority flood, no minimum slope).")
            filled, _, _ = fill_depressions(dem, valid, cell_width=info.cell_width,
                                            cell_height=info.cell_height,
                                            progress=self.make_progress(feedback, 0.6))
        road = self.read_alignment(parameters, "ROAD", context, info, feedback,
                                   start_chainage=self.parameterAsDouble(parameters, "START", context),
                                   reverse=self.parameterAsBool(parameters, "REVERSE", context))
        crossings = None
        if road is not None and self.parameterAsSource(parameters, "CROSSINGS", context) is not None:
            pts = self.read_pour_points(parameters, "CROSSINGS", context, info, feedback,
                                        extra_fields=("outlet_uid", "crossing_id"))
            ch, _, _ = road.locate([p["x"] for p in pts], [p["y"] for p in pts])
            crossings = [(str(p.get("attr_outlet_uid") or p.get("attr_crossing_id")
                              or f"fid {p['fid']}"), float(c)) for p, c in zip(pts, ch)]
        s = review_settings(self, parameters, context)
        rows, depth, summ = depression_table(dem, filled, valid, info.geotransform,
                                             alignment=road, crossings=crossings, **s)
        out = self.parameterAsOutputLayer(parameters, "OUTPUT", context)
        write_depressions(self, out, rows, info.projection_wkt)
        self.report_stats(feedback, "Fill depressions", {
            "depressions": summ["depressions"], "listed": summ["reported"],
            "raised cells": summ["raised_cells"], "fill volume (m3)": summ["volume_m3"],
            "flagged against the road": summ["flag_at_road"],
            "flagged near the road": summ["flag_near_road"]})
        for r in rows:
            if r["flag"]:
                feedback.pushWarning(f"Depression {r['dep_id']} at chainage "
                                     f"{r['chainage_m']:,.0f}: {r['note']}")
        result = {"OUTPUT": out}
        dpath = self.parameterAsOutputLayer(parameters, "DEPTH", context)
        if dpath:
            write_raster(dpath, depth.astype("float32"), info, valid=valid,
                         metadata={"QEHT_TOOL": "Review fill depressions"})
            result["DEPTH"] = dpath
        return result
