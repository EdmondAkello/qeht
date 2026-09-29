# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""D8 Flow Direction, standard D8-coded output."""

from qgis.core import (
    QgsProcessingException,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterRasterDestination,
)

from .base import QehtAlgorithm
from ..core.raster import read_dem, write_raster, audit_nodata, audit_resampling
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
            "• <i>Barnes 2014</i> (default since 0.13): two gradients, toward the "
            "flat's outlet and away from higher ground, give convergent drainage "
            "across flats. QEHT reproduces Barnes' own implementation (RichDEM) "
            "cell for cell, and it agreed best with TauDEM in the WP-G benchmark "
            "(80 Kenyan areas, flat to mountainous, ALOS and FABDEM).\n"
            "• <i>Toward lower terrain</i>: each flat cell drains along the "
            "shortest path to the flat's outlet. Gives parallel drainage lines on "
            "broad flats; use it to match outputs of platforms that route flats "
            "this way (it was the default up to 0.12).\n"
            "The choice matters on large flats - common on DEMs stored in whole "
            "metres (e.g. ALOS AW3D30), rare on floating-point DEMs (e.g. FABDEM). "
            "Flats that touch the grid edge or NoData drain to it.\n\n"
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
            # index order kept from 0.12 (0 = toward) so saved models keep their method
            options=["Toward lower terrain (shortest path to the outlet)",
                     "Barnes 2014 convergent (recommended)"],
            defaultValue=1))
        self.addParameter(QgsProcessingParameterRasterDestination(
            OUTPUT, "Flow direction (D8-coded)"))

    def processAlgorithm(self, parameters, context, feedback):
        dem_path = self.raster_path(parameters, DEM, context)
        resolve = self.parameterAsBool(parameters, RESOLVE_FLATS, context)
        choice = self.parameterAsEnum(parameters, FLAT_METHOD, context)
        if choice not in (0, 1):
            raise QgsProcessingException(
                "The hybrid flat method was removed in QEHT 0.13 (it reproduced Barnes). "
                "Choose Barnes 2014 or Toward lower terrain.")
        flat_method = ["toward", "barnes"][choice]
        out_path = self.parameterAsOutputLayer(parameters, OUTPUT, context)

        dem, valid, info = read_dem(dem_path)
        self.check_size(feedback, info.rows, info.cols,
                        ["flow_direction" if flat_method == "toward" else "flow_direction_barnes"])
        for warning in audit_nodata(dem, valid, info.nodata) + audit_resampling(dem, valid):
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
