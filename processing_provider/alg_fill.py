# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Fill Depressions - DEM conditioning.

Deliberately a separate, explicit step. Conditioning is a modelling
decision that must be documented in a drainage report, not an invisible
side effect of computing flow direction.
"""

from qgis.core import (
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterDestination,
)

from .base import QehtAlgorithm
from ..core.raster import read_dem, write_raster, audit_nodata
from ..core.conditioning.fill import fill_depressions, depression_depth

DEM = "DEM"
MIN_SLOPE = "MIN_SLOPE"
NODATA_OVERRIDE = "NODATA_OVERRIDE"
OUTPUT = "OUTPUT"
DEPTH = "DEPTH"


class FillDepressionsAlgorithm(QehtAlgorithm):

    def name(self):
        return "filldepressions"

    def displayName(self):
        return "Fill depressions"

    def shortHelpString(self):
        return (
            "Fills closed depressions in a DEM using the priority-flood "
            "algorithm (Barnes et al. 2014).\n\n"
            "Equivalent to the 'Fill' tool in commercial raster hydrology toolsets.\n\n"
            "<b>Minimum slope</b> imposes a small gradient across filled "
            "surfaces so that flow routing across them is defined. Leave at "
            "0 for a classic flat fill; 0.0001 is a common choice that "
            "largely eliminates flats before routing.\n\n"
            "<b>Fill depth</b> is an optional QA raster showing how much "
            "each cell was raised. Inspect it before accepting the result: "
            "large contiguous depths usually mark a road embankment or dam "
            "that the DEM treats as a barrier, and may indicate a culvert "
            "that should be burned in instead of filled over."
        )

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(DEM, "Input DEM"))
        self.addParameter(QgsProcessingParameterNumber(
            MIN_SLOPE, "Minimum slope across filled surfaces (m/m)",
            QgsProcessingParameterNumber.Double, defaultValue=0.0,
            minValue=0.0, maxValue=0.1))
        self.addParameter(QgsProcessingParameterNumber(
            NODATA_OVERRIDE, "Treat this value as NoData (leave blank to trust the header)",
            QgsProcessingParameterNumber.Double, optional=True))
        self.addParameter(QgsProcessingParameterRasterDestination(
            OUTPUT, "Conditioned DEM"))
        self.addParameter(QgsProcessingParameterRasterDestination(
            DEPTH, "Fill depth (QA)", optional=True, createByDefault=False))

    def processAlgorithm(self, parameters, context, feedback):
        dem_path = self.raster_path(parameters, DEM, context)
        min_slope = self.parameterAsDouble(parameters, MIN_SLOPE, context)
        out_path = self.parameterAsOutputLayer(parameters, OUTPUT, context)
        depth_path = self.parameterAsOutputLayer(parameters, DEPTH, context)

        feedback.pushInfo(f"Reading DEM: {dem_path}")
        override = None
        if parameters.get(NODATA_OVERRIDE) is not None:
            override = self.parameterAsDouble(parameters, NODATA_OVERRIDE, context)

        dem, valid, info = read_dem(dem_path, nodata_override=override)
        feedback.pushInfo(f"{info}  valid cells: {int(valid.sum()):,}")

        self.check_size(feedback, info.rows, info.cols)
        for warning in audit_nodata(dem, valid, info.nodata):
            feedback.pushWarning("DEM QA: " + warning)

        filled, n_filled, seeds = fill_depressions(
            dem, valid, min_slope=min_slope,
            cell_width=info.cell_width, cell_height=info.cell_height,
            progress=self.make_progress(feedback, weight=0.9))

        self.report_stats(feedback, "Fill", {
            "boundary outlet cells seeded": seeds,
            "cells raised": n_filled,
            "percent of valid cells raised":
                100.0 * n_filled / max(1, int(valid.sum())),
        })

        write_raster(out_path, filled, info, valid=valid)
        results = {OUTPUT: out_path}

        if depth_path:
            depth = depression_depth(dem, filled, valid)
            write_raster(depth_path, depth, info, valid=valid)
            feedback.pushInfo(
                f"Maximum fill depth: {float(depth[valid].max()):.2f} m")
            results[DEPTH] = depth_path

        return results
