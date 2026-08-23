# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Flow accumulation with selectable output quantity."""

from qgis.core import (
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterEnum,
    QgsProcessingParameterRasterDestination,
)
import numpy as np

from .base import QehtAlgorithm
from ..core.raster import read_dem, write_raster
from ..core.grid import decode_d8
from ..core.flow.accumulation import flow_accumulation

FDR = "FDR"
QUANTITY = "QUANTITY"
WEIGHT = "WEIGHT"
OUTPUT = "OUTPUT"

QUANTITIES = ["Cell count (standard D8 convention)",
              "Contributing area (map units squared)",
              "Weighted (use weight raster)"]


class FlowAccumulationAlgorithm(QehtAlgorithm):

    def name(self):
        return "flowaccumulation"

    def displayName(self):
        return "Flow accumulation"

    def shortHelpString(self):
        return (
            "Accumulates upstream contribution over a D8 flow direction "
            "raster by topological traversal - O(N), no recursion.\n\n"
            "Input is a standard D8-coded flow direction raster, so a grid "
            "produced by another GIS/hydrology toolset can be used directly. "
            "That is the "
            "cleanest way to cross-check this toolkit against a trusted "
            "production result: feed both the same direction grid and "
            "compare accumulation cell by cell.\n\n"
            "<b>Cell count</b> follows the standard D8 convention and excludes the "
            "cell itself, so a ridge cell is 0.\n"
            "<b>Contributing area</b> multiplies by cell area, which is the "
            "quantity you actually want for a design flow calculation and "
            "which survives a change in DEM resolution."
        )

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(
            FDR, "Flow direction (D8-coded)"))
        self.addParameter(QgsProcessingParameterEnum(
            QUANTITY, "Accumulated quantity", options=QUANTITIES,
            defaultValue=0))
        self.addParameter(QgsProcessingParameterRasterLayer(
            WEIGHT, "Weight raster", optional=True))
        self.addParameter(QgsProcessingParameterRasterDestination(
            OUTPUT, "Flow accumulation"))

    def processAlgorithm(self, parameters, context, feedback):
        fdr_path = self.raster_path(parameters, FDR, context)
        quantity = self.parameterAsEnum(parameters, QUANTITY, context)
        out_path = self.parameterAsOutputLayer(parameters, OUTPUT, context)

        d8, valid, info = read_dem(fdr_path)
        direction = decode_d8(d8.astype(np.int32))

        weights = None
        if quantity == 1:
            weights = np.full(d8.shape, info.cell_area, dtype=np.float64)
            feedback.pushInfo(f"Cell area: {info.cell_area:,.2f} map units squared")
        elif quantity == 2:
            wlayer = self.parameterAsRasterLayer(parameters, WEIGHT, context)
            if wlayer is None:
                from qgis.core import QgsProcessingException
                raise QgsProcessingException(
                    "Weighted accumulation requires a weight raster.")
            weights, wvalid, winfo = read_dem(self.raster_path(parameters, WEIGHT, context))
            if weights.shape != d8.shape:
                from qgis.core import QgsProcessingException
                raise QgsProcessingException(
                    "Weight raster must match the flow direction grid exactly "
                    f"({winfo.rows}x{winfo.cols} vs {info.rows}x{info.cols}).")
            weights = np.where(wvalid, weights, 0.0)

        accum, stats = flow_accumulation(
            direction, valid, weights=weights,
            progress=self.make_progress(feedback, weight=0.9))

        self.report_stats(feedback, "Accumulation", stats)
        if stats["cells_in_cycles"]:
            feedback.pushWarning(
                f"{stats['cells_in_cycles']:,} cells form circular flow paths "
                "and were not accumulated. The flow direction raster is not a "
                "valid drainage tree.")

        write_raster(out_path, accum, info, valid=valid)
        return {OUTPUT: out_path}
