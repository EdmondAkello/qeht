# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""D8 Flow Direction, standard D8-coded output."""

from qgis.core import (
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterRasterDestination,
)

from .base import QehtAlgorithm
from ..core.raster import read_dem, write_raster, audit_nodata
from ..core.flow.direction import d8_direction
from ..core.grid import encode_d8

DEM = "DEM"
RESOLVE_FLATS = "RESOLVE_FLATS"
FLAT_METHOD = "FLAT_METHOD"
OUTPUT = "OUTPUT"


class FlowDirectionAlgorithm(QehtAlgorithm):

    def name(self):
        return "flowdirection"

    def displayName(self):
        return "D8 flow direction"

    def shortHelpString(self):
        return (
            "Computes D8 flow direction by steepest descent, weighted by "
            "distance so that diagonal neighbours are compared over "
            "cellsize x sqrt(2).\n\n"
            "Output uses the standard D8 encoding (E=1, SE=2, S=4, SW=8, W=16, "
            "NW=32, N=64, NE=128; 0 for sinks and cells draining off the "
            "grid edge), so the raster is directly interchangeable with "
            "other GIS/hydrology toolsets that use the same convention.\n\n"
            "<b>Input should already be conditioned</b> - run Fill "
            "depressions first. Running this on a raw DEM leaves every "
            "closed depression as an unrouted sink.\n\n"
            "<b>Flat resolution</b> (cells with no downhill neighbour after "
            "filling):\n"
            "• <i>Toward lower terrain</i> (default): each flat cell drains along "
            "the shortest path to the flat's outlet. Validated against a reference "
            "hydrology toolset: stream IoU 0.60-0.67 on flat coastal terrain.\n"
            "• <i>Barnes 2014</i>: two gradients, toward the outlet and away from "
            "higher ground, flat_mask = 2*toward + (height - away); gives "
            "convergent drainage on broad flats. <b>Fixed in 0.12:</b> up to 0.11 "
            "its outlet rule could point two cells at each other, forming small "
            "flow loops that cut main channels off from their upstream area; the "
            "earlier finding that Barnes matched the reference poorly on coastal "
            "flats (IoU ~0.19) used that code and should be repeated.\n"
            "• <i>Hybrid</i>: flat_mask = w*toward + (height - away). w = 2 is "
            "Barnes; large w follows shortest paths to the outlet; 1 < w < 2 "
            "pushes flow off high edges more strongly. Any w > 1 is loop-free. "
            "On broad real flats w changes little: how ties between equally "
            "short paths are broken matters more. Optionally, flats smaller than "
            "a given size use the toward gradient only.\n\n"
            "Ties between equally steep neighbours follow the fixed priority S, W, "
            "N, E, SE, SW, NW, NE (recovered from a reference platform). Filling "
            "with a small minimum slope removes most flats before routing."
        )

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(
            DEM, "Conditioned DEM"))
        self.addParameter(QgsProcessingParameterBoolean(
            RESOLVE_FLATS, "Resolve flats", defaultValue=True))
        from qgis.core import QgsProcessingParameterEnum
        self.addParameter(QgsProcessingParameterEnum(
            FLAT_METHOD, "Flat resolution method",
            options=["Toward lower terrain (recommended - closest to reference platforms)",
                     "Barnes 2014 convergent (iterated)",
                     "Hybrid: weighted gradient w*toward + (height - away)"],
            defaultValue=0))
        from qgis.core import QgsProcessingParameterNumber
        self.addParameter(QgsProcessingParameterNumber(
            "FLAT_WEIGHT", "Hybrid weight w (> 1; 2 = Barnes, large = toward-lower)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=2.0, minValue=1.0001))
        self.addParameter(QgsProcessingParameterNumber(
            "SMALL_FLATS", "Hybrid: flats smaller than this many cells use toward only (0 = off)",
            QgsProcessingParameterNumber.Type.Integer, defaultValue=0, minValue=0))
        self.addParameter(QgsProcessingParameterRasterDestination(
            OUTPUT, "Flow direction (D8-coded)"))

    def processAlgorithm(self, parameters, context, feedback):
        dem_path = self.raster_path(parameters, DEM, context)
        resolve = self.parameterAsBool(parameters, RESOLVE_FLATS, context)
        flat_method = ["toward", "barnes", "hybrid"][self.parameterAsEnum(parameters, FLAT_METHOD, context)]
        flat_weight = self.parameterAsDouble(parameters, "FLAT_WEIGHT", context) \
            if "FLAT_WEIGHT" in parameters else 2.0
        small_flats = self.parameterAsInt(parameters, "SMALL_FLATS", context) \
            if "SMALL_FLATS" in parameters else 0
        out_path = self.parameterAsOutputLayer(parameters, OUTPUT, context)

        dem, valid, info = read_dem(dem_path)
        self.check_size(feedback, info.rows, info.cols,
                        ["flow_direction" if flat_method == "toward" else "flow_direction_barnes"])
        for warning in audit_nodata(dem, valid, info.nodata):
            feedback.pushWarning("DEM QA: " + warning)
        direction, stats = d8_direction(
            dem, valid, cell_width=info.cell_width,
            cell_height=info.cell_height, resolve_flats=resolve,
            flat_method=flat_method, flat_weight=flat_weight, small_flat_cells=small_flats,
            progress=self.make_progress(feedback, weight=0.9))

        self.report_stats(feedback, "Flow direction", stats)
        no_descent_share = 100.0 * stats["cells_without_descent"] / max(1, stats["cells_valid"])
        if no_descent_share > 10.0:
            feedback.pushWarning(
                f"{stats['cells_without_descent']:,} cells ({no_descent_share:.1f}%) had no "
                "downhill neighbour. Above ~10% this usually means the DEM was not "
                "conditioned, or that undeclared NoData (often 0) is being treated as "
                "real terrain - a flat plateau at 0 m has no descent anywhere. Check the "
                "DEM QA warnings and set the NoData override on Fill depressions.")

        if stats["cells_still_unrouted"] > stats["cells_valid"] * 0.01:
            feedback.pushWarning(
                f"{stats['cells_still_unrouted']:,} cells remain unrouted "
                f"({100.0*stats['cells_still_unrouted']/max(1,stats['cells_valid']):.1f}%). "
                "If this DEM was not conditioned, run Fill depressions first.")

        d8 = encode_d8(direction, valid)
        write_raster(out_path, d8, info, dtype="int32", nodata=-1, valid=valid)
        return {OUTPUT: out_path}
