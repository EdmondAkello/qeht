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
            "Flat resolution uses a gradient-toward-lower-terrain BFS "
            "(the only method wired in v0.7). It routes flat cells toward "
            "the nearest outlet. On broad flats this can produce somewhat "
            "parallel flow rather than the convergent pattern a reference "
            "GIS platform produces; validated against a reference hydrology "
            "toolset's accumulation it still "
            "reproduces the channel network to a stream IoU of 0.60-0.67. "
            "A convergent two-gradient Barnes 2014 resolver (iterated to "
            "convergence) is available as an option. On validation it "
            "produces correct convergent drainage on synthetic saddles, but "
            "on real flat coastal terrain it disperses flow MORE than the "
            "toward-lower method and matches the reference platform's "
            "accumulation less well "
            "(stream IoU ~0.19 vs ~0.63 at threshold 1000). The reference "
            "platform's actual "
            "flat behaviour, despite citing Garbrecht and Martz, empirically "
            "resembles the simpler method here. Toward-lower is therefore the "
            "recommended default. "
            "The 'gradient away from higher terrain' component of "
            "Garbrecht & Martz (1997) is not implemented, so drainage "
            "across wide flats converges into fewer channels than the "
            "reference platform "
            "produces. Filling with a small minimum slope avoids this."
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
                     "Barnes 2014 convergent (iterated)"],
            defaultValue=0))
        self.addParameter(QgsProcessingParameterRasterDestination(
            OUTPUT, "Flow direction (D8-coded)"))

    def processAlgorithm(self, parameters, context, feedback):
        dem_path = self.raster_path(parameters, DEM, context)
        resolve = self.parameterAsBool(parameters, RESOLVE_FLATS, context)
        flat_method = "toward" if self.parameterAsEnum(parameters, FLAT_METHOD, context) == 0 else "barnes"
        out_path = self.parameterAsOutputLayer(parameters, OUTPUT, context)

        dem, valid, info = read_dem(dem_path)
        self.check_size(feedback, info.rows, info.cols)
        for warning in audit_nodata(dem, valid, info.nodata):
            feedback.pushWarning("DEM QA: " + warning)
        direction, stats = d8_direction(
            dem, valid, cell_width=info.cell_width,
            cell_height=info.cell_height, resolve_flats=resolve,
            flat_method=flat_method,
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
