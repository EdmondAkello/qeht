# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Stream network extraction with Strahler ordering."""

from qgis.core import (
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterEnum,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterDestination,
)
import numpy as np

from .base import QehtAlgorithm
from ..core.raster import read_dem, write_raster
from ..core.grid import decode_d8
from ..core.watershed.delineate import extract_streams
from ..core.flow.accumulation import strahler_order

FDR = "FDR"
FAC = "FAC"
MODE = "MODE"
THRESHOLD = "THRESHOLD"
STREAMS = "STREAMS"
ORDER = "ORDER"

MODES = ["Accumulation threshold (cells)", "Contributing area threshold (map units squared)"]


class StreamNetworkAlgorithm(QehtAlgorithm):

    def name(self):
        return "streamnetwork"

    def displayName(self):
        return "Stream network and Strahler order"

    def shortHelpString(self):
        return (
            "Extracts a channel network by thresholding flow accumulation, "
            "then assigns Strahler order to the extracted network.\n\n"
            "Strahler order is computed <b>on the stream network</b>, not "
            "over every cell of the DEM. These are different quantities and "
            "are routinely confused - whole-DEM ordering has no hydrological "
            "meaning off-channel.\n\n"
            "Prefer the area threshold in reports: it is independent of DEM "
            "resolution, so a reviewer can reproduce your network from a "
            "different source DEM."
        )

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(FDR, "Flow direction (D8-coded)"))
        self.addParameter(QgsProcessingParameterRasterLayer(FAC, "Flow accumulation (cell count)"))
        self.addParameter(QgsProcessingParameterEnum(MODE, "Threshold type", options=MODES, defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber(
            THRESHOLD, "Threshold value", QgsProcessingParameterNumber.Type.Double,
            defaultValue=1000.0, minValue=1.0))
        self.addParameter(QgsProcessingParameterRasterDestination(STREAMS, "Stream mask"))
        self.addParameter(QgsProcessingParameterRasterDestination(ORDER, "Strahler order"))

    def processAlgorithm(self, parameters, context, feedback):
        fdr_path = self.raster_path(parameters, FDR, context)
        fac_path = self.raster_path(parameters, FAC, context)
        mode = self.parameterAsEnum(parameters, MODE, context)
        threshold = self.parameterAsDouble(parameters, THRESHOLD, context)
        streams_path = self.parameterAsOutputLayer(parameters, STREAMS, context)
        order_path = self.parameterAsOutputLayer(parameters, ORDER, context)

        d8, valid, info = read_dem(fdr_path)
        accum, _, ainfo = read_dem(fac_path)
        if accum.shape != d8.shape:
            from qgis.core import QgsProcessingException
            raise QgsProcessingException("Accumulation and direction rasters must match.")

        direction = decode_d8(d8.astype(np.int32))

        if mode == 0:
            mask = extract_streams(accum, valid, threshold_cells=threshold)
            feedback.pushInfo(f"Threshold: {threshold:,.0f} cells "
                              f"= {threshold * info.cell_area:,.0f} map units squared")
        else:
            mask = extract_streams(accum, valid, threshold_area=threshold,
                                   cell_area=info.cell_area)
            feedback.pushInfo(f"Threshold: {threshold:,.0f} map units squared "
                              f"= {threshold / info.cell_area:,.1f} cells")

        order = strahler_order(direction, valid, mask,
                               progress=self.make_progress(feedback, weight=0.8))

        n_stream = int(mask.sum())
        self.report_stats(feedback, "Stream network", {
            "stream cells": n_stream,
            "percent of catchment": 100.0 * n_stream / max(1, int(valid.sum())),
            "maximum Strahler order": int(order.max()),
            "channel length (map units)": n_stream * info.cell_width,
        })
        if n_stream == 0:
            feedback.pushWarning("No cells exceed the threshold. Lower it.")

        write_raster(streams_path, mask.astype(np.int32), info, dtype="int32",
                     nodata=-1, valid=valid)
        write_raster(order_path, order, info, dtype="int32", nodata=-1, valid=valid)
        return {STREAMS: streams_path, ORDER: order_path}
