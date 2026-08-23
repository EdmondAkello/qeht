# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Shared base for QEHT Processing algorithms.

This is the ONLY place where QGIS and `qeht.core` meet. The core never
imports QGIS; the algorithms below never implement hydrology. If you
find yourself writing a loop over cells in this package, it belongs in
core instead.
"""

from qgis.core import (
    QgsProcessingAlgorithm,
    QgsProcessingException,
)


class QehtAlgorithm(QgsProcessingAlgorithm):
    """Common metadata and helpers."""

    def group(self):
        return "Terrain and drainage"

    def groupId(self):
        return "terrain"

    def createInstance(self):
        return type(self)()

    # -- helpers ----------------------------------------------------

    def raster_path(self, parameters, name, context):
        """Resolve a raster parameter to a filesystem path.

        We read through GDAL directly rather than through the QGIS raster
        provider, so we need a real path. In-memory and remote layers are
        rejected with a clear message rather than failing obscurely later.
        """
        layer = self.parameterAsRasterLayer(parameters, name, context)
        if layer is None:
            raise QgsProcessingException(f"Could not resolve raster parameter '{name}'.")
        path = layer.source().split("|")[0]
        import os
        if not os.path.exists(path):
            raise QgsProcessingException(
                f"Layer '{layer.name()}' is not backed by a file GDAL can open "
                f"directly ({path}). Save it to GeoTIFF first.")
        return path

    def make_progress(self, feedback, weight=1.0, offset=0.0):
        """Adapt core's progress(fraction, message) to QGIS feedback."""
        def _progress(fraction, message=""):
            if feedback.isCanceled():
                raise QgsProcessingException("Cancelled by user.")
            feedback.setProgress(int((offset + fraction * weight) * 100))
            if message:
                feedback.setProgressText(message)
        return _progress

    def check_size(self, feedback, rows, cols):
        """Warn before a very large DEM exhausts memory.

        QEHT holds the whole grid in memory - several float64 and int64
        copies during routing. Past validation, a 41-megapixel DEM
        (6400x6400) needed well over 8 GB and was killed on a 16 GB
        machine at the flow-direction step. There is no tiling in this
        version, so the honest guidance is to clip first.
        """
        megapixels = rows * cols / 1.0e6
        # rough: ~8 working copies of float64 = 64 bytes/cell peak
        est_gb = rows * cols * 64 / 1.0e9
        if megapixels > 25:
            feedback.pushWarning(
                f"Large DEM: {rows:,} x {cols:,} = {megapixels:,.0f} megapixels, "
                f"needing roughly {est_gb:,.1f} GB of RAM at peak. QEHT processes "
                "the grid whole - there is no tiling yet. If this runs out of "
                "memory, clip the DEM to your catchment plus a buffer and run "
                "again. Clipping is standard practice and does not affect "
                "results within the clipped area.")
        return megapixels

    def report_stats(self, feedback, title, stats):
        """Print a stats block. Handles str/int/float without assuming type -
        a string value here previously raised
        'ValueError: Cannot specify comma with s' and killed the algorithm
        AFTER all the real work was done."""
        feedback.pushInfo(f"--- {title} ---")
        for key, value in stats.items():
            if isinstance(value, bool):
                text = str(value)
            elif isinstance(value, float):
                text = f"{value:,.6g}"
            elif isinstance(value, int):
                text = f"{value:,}"
            else:
                text = str(value)
            feedback.pushInfo(f"    {key}: {text}")
