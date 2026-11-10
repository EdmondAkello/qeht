# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Per-crossing sensitivity table from a design hydrology package (1.2)."""

import os

from qgis.core import (QgsProcessingParameterFile, QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterFileDestination, QgsProcessingException)

from .base import QehtAlgorithm


def run_sensitivity(package, filled, valid, info, base_direction=None, progress=None):
    """-> (rows, summary) for every crossing of the package on this grid."""
    from ..core.review.sensitivity import crossing_sensitivity, package_inputs
    uids, rc, st, err, orc = package_inputs(package, info.geotransform, (info.rows, info.cols))
    if not uids:
        return [], []
    return crossing_sensitivity(filled, valid, info.cell_width, info.cell_height, rc, uids,
                                st["threshold_cells"], st["snap_cells"], st["flat_method"],
                                dem_error=err, base_direction=base_direction, outlets_rc=orc,
                                progress=progress)


class CrossingSensitivityAlgorithm(QehtAlgorithm):

    def name(self): return "crossingsensitivity"
    def displayName(self): return "Crossing sensitivity table"
    def group(self): return "Review"
    def groupId(self): return "review"

    def shortHelpString(self):
        return (
            "How much each crossing's catchment area and longest flow path move when one "
            "setting changes. Every crossing is re-measured with:\n"
            "- the other flat method (Barnes / toward lower);\n"
            "- half and twice the snap radius;\n"
            "- half and twice the stream threshold (the snapping target);\n"
            "- the DEM-error P10 and P90, copied from the package when DEM uncertainty was run.\n\n"
            "The base row is the run's own outlet, so it equals the package. The snap and "
            "threshold variants snap the crossing's input point again; the flat-method variant "
            "snaps the run's outlet onto the other method's streams. Outlet moved = the snapped outlet is more "
            "than one cell from the base outlet; lost = no stream within the snap radius.\n\n"
            "Writes the long table (one row per crossing and variant) and a summary (largest "
            "change in area and flow path per crossing and the factor behind it). Give the "
            "filled DEM the run routed on; the run's flow direction is optional and saves "
            "routing again. The pipeline writes both tables and adds them to the crossing "
            "review sheets.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile(
            "PACKAGE", "Design hydrology package", fileFilter="GeoPackage (*.gpkg)"))
        self.addParameter(QgsProcessingParameterRasterLayer("FILLED", "Filled DEM of the run"))
        self.addParameter(QgsProcessingParameterRasterLayer(
            "FDR", "Flow direction of the run (optional)", optional=True))
        self.addParameter(QgsProcessingParameterFileDestination(
            "OUTPUT", "Sensitivity table", fileFilter="CSV (*.csv)"))
        self.addParameter(QgsProcessingParameterFileDestination(
            "SUMMARY_CSV", "Sensitivity summary", fileFilter="CSV (*.csv)", optional=True,
            createByDefault=True))

    def processAlgorithm(self, parameters, context, feedback):
        from ..core.interop.heas_exchange import validate_exchange
        from ..core.raster import read_dem
        from ..core.grid import decode_d8
        from ..core.review.sensitivity import write_csv, SUMMARY_FIELDS
        pkg = self.parameterAsFile(parameters, "PACKAGE", context)
        errors, _ = validate_exchange(pkg)
        if errors:
            raise QgsProcessingException("Not a valid design hydrology package: " + "; ".join(errors))
        filled, valid, info = read_dem(self.raster_path(parameters, "FILLED", context))
        self.check_size(feedback, info.rows, info.cols, steps=("flow_direction_barnes",
                                                                "flow_direction"))
        base = None
        if self.parameterAsRasterLayer(parameters, "FDR", context) is not None:
            d8, dv, di = read_dem(self.raster_path(parameters, "FDR", context))
            if (di.rows, di.cols) != (info.rows, info.cols):
                raise QgsProcessingException("The flow direction is not on the filled DEM's grid.")
            base = decode_d8(d8.astype(int))
        rows, summ = run_sensitivity(pkg, filled, valid, info, base,
                                     self.make_progress(feedback))
        if not rows:
            raise QgsProcessingException("No crossing of the package falls on the filled DEM.")
        out = write_csv(rows, self.parameterAsFileOutput(parameters, "OUTPUT", context))
        result = {"OUTPUT": out}
        sp = self.parameterAsFileOutput(parameters, "SUMMARY_CSV", context)
        if sp:
            result["SUMMARY_CSV"] = write_csv(summ, sp, SUMMARY_FIELDS)
        top = sorted((s for s in summ if s["sens_max_area_pct"] is not None),
                     key=lambda s: -s["sens_max_area_pct"])[:5]
        for s in top:
            feedback.pushInfo(f"{s['outlet_uid']}: area up to {s['sens_max_area_pct']:.0f} % "
                              f"({s['sens_area_factor']})")
        feedback.pushInfo(f"{len(summ)} crossing(s) -> {out}")
        return result
