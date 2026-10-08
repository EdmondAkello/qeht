# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Prepare DEM for hydrology (F9, v0.21)."""

import os

from qgis.core import (
    QgsProcessing, QgsProcessingException, QgsProcessingParameterMultipleLayers,
    QgsProcessingParameterFeatureSource, QgsProcessingParameterExtent,
    QgsProcessingParameterNumber, QgsProcessingParameterEnum, QgsProcessingParameterCrs,
    QgsProcessingParameterString, QgsProcessingParameterRasterDestination,
    QgsProcessingParameterFileDestination, QgsProcessingOutputString,
)

from .base import QehtAlgorithm
from ..core.conditioning.prepare import prepare_dem, RESAMPLING

CRS_OPTIONS = ["Auto: WGS 84 / UTM zone of the extent centre (UPS near the poles)",
               "The CRS chosen below (any projected CRS, e.g. a national grid)"]


class PrepareDemAlgorithm(QehtAlgorithm):

    def name(self): return "preparedem"
    def displayName(self): return "Prepare DEM for hydrology"

    def shortHelpString(self):
        return (
            "Merges DEM tiles, clips them, reprojects them to a projected metric CRS and checks "
            "the result, all in-process with GDAL. Most bad runs start here: a DEM reprojected "
            "with nearest neighbour, an export finer than the native grid (repeated rows and "
            "columns), or NoData zeros the file does not declare.\n\n"
            "<b>Steps:</b> a virtual mosaic of the tiles → optional clip to a polygon layer or an "
            "extent, widened by the buffer (default 2,000 m) → warp to the target CRS and cell "
            "size, float32, NoData −9999, pixels aligned to the cell size → NoData and "
            "resampling audits on the input and the output → tags in the GeoTIFF.\n\n"
            "<b>Target CRS:</b> by default the WGS 84 / UTM zone of the extent centre, anywhere "
            "on Earth (e.g. 37S at 37.5° E, 1° S; 33N at 15° E, 50° N; UPS beyond 84° N or "
            "80° S; the Norway and Svalbard zone exceptions apply). The log notes an extent "
            "that spans several zones. Or choose any projected CRS, such as a national grid. "
            "<b>Cell size:</b> blank = the native cell converted to metres "
            "(the larger of x and y, rounded to 0.1 m).\n\n"
            "<b>Resampling:</b> bilinear (default) or cubic. Nearest neighbour is not offered: on "
            "a finer or rotated grid it repeats whole source rows and columns, which creates "
            "artificial flats and ties, and flow routing on them is arbitrary.\n\n"
            "<b>If the input already repeats rows or columns</b>, the damage is upstream and "
            "resampling cannot undo it: get the native product (the tiles at their own grid) "
            "and prepare it here.\n\n"
            "<b>NoData override:</b> a value to treat as NoData in the tiles (e.g. 0 where the "
            "header declares nothing).\n\n"
            "<b>DEM source</b> (e.g. 'FABDEM v1.2') is written to the GeoTIFF tags with the "
            "preparation record; the design hydrology package reads it when its own DEM source "
            "is left blank. The log (JSON) records the inputs, CRS, cell size, resampling, both "
            "audits and a native-grid estimate.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterMultipleLayers(
            "DEMS", "DEM raster(s) (tiles of one product)", QgsProcessing.SourceType.TypeRaster))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "CLIP", "Clip to polygons (optional)", [QgsProcessing.SourceType.TypeVectorPolygon],
            optional=True))
        self.addParameter(QgsProcessingParameterExtent(
            "EXTENT", "Clip to extent (optional; used when no clip polygons)", optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            "BUFFER", "Clip buffer (m)", QgsProcessingParameterNumber.Type.Double,
            defaultValue=2000.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterEnum("CRS_MODE", "Target CRS", options=CRS_OPTIONS,
                                                     defaultValue=0))
        self.addParameter(QgsProcessingParameterCrs("TARGET_CRS", "Target CRS (when chosen above)",
                                                    optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            "CELL", "Cell size (m; blank = native, converted to metres)",
            QgsProcessingParameterNumber.Type.Double, optional=True, minValue=0.0))
        self.addParameter(QgsProcessingParameterEnum(
            "RESAMPLING", "Resampling", options=["Bilinear (recommended)", "Cubic"], defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber(
            "NODATA", "Treat this input value as NoData (optional)",
            QgsProcessingParameterNumber.Type.Double, optional=True))
        self.addParameter(QgsProcessingParameterString(
            "SOURCE", "DEM source for the record (e.g. FABDEM v1.2)", optional=True))
        self.addParameter(QgsProcessingParameterRasterDestination("OUTPUT", "Prepared DEM"))
        self.addParameter(QgsProcessingParameterFileDestination(
            "LOG", "Preparation log", fileFilter="JSON (*.json)", optional=True))
        self.addOutput(QgsProcessingOutputString("SUMMARY", "Summary"))

    def processAlgorithm(self, parameters, context, feedback):
        layers = self.parameterAsLayerList(parameters, "DEMS", context)
        paths = []
        for lyr in layers:
            src = lyr.source().split("|")[0]
            if not os.path.exists(src):
                raise QgsProcessingException(f"{lyr.name()}: QEHT reads rasters from files; "
                                             f"'{src}' is not a file.")
            paths.append(src)
        if not paths:
            raise QgsProcessingException("Give at least one DEM raster.")
        mode = self.parameterAsEnum(parameters, "CRS_MODE", context)
        if mode == 0:
            target = "auto"
        else:
            crs = self.parameterAsCrs(parameters, "TARGET_CRS", context)
            if not crs.isValid():
                raise QgsProcessingException("Choose the target CRS.")
            if crs.isGeographic():
                raise QgsProcessingException("The target CRS is geographic; choose a projected, "
                                             "metric CRS (e.g. UTM).")
            target = crs.toWkt()
        clip_bounds, clip_crs = None, None
        clip = self.parameterAsSource(parameters, "CLIP", context)
        if clip is not None and clip.featureCount() > 0:
            ext = clip.sourceExtent()
            clip_bounds = (ext.xMinimum(), ext.yMinimum(), ext.xMaximum(), ext.yMaximum())
            clip_crs = clip.sourceCrs().toWkt()
        elif parameters.get("EXTENT") not in (None, ""):
            ecrs = self.parameterAsExtentCrs(parameters, "EXTENT", context)
            ext = self.parameterAsExtent(parameters, "EXTENT", context, ecrs)
            if not ext.isNull() and ext.width() > 0:
                clip_bounds = (ext.xMinimum(), ext.yMinimum(), ext.xMaximum(), ext.yMaximum())
                clip_crs = ecrs.toWkt() if ecrs.isValid() else None
        cell = self.parameterAsDouble(parameters, "CELL", context) \
            if parameters.get("CELL") not in (None, "") else None
        nodata = self.parameterAsDouble(parameters, "NODATA", context) \
            if parameters.get("NODATA") not in (None, "") else None
        out = self.parameterAsOutputLayer(parameters, "OUTPUT", context)
        log_path = self.parameterAsFileOutput(parameters, "LOG", context) \
            or os.path.splitext(out)[0] + "_prepare_log.json"
        feedback.pushInfo(f"{len(paths)} tile(s); target {CRS_OPTIONS[mode]}")
        try:
            log = prepare_dem(paths, out, target=target, cell_m=cell or None,
                              resampling=RESAMPLING[self.parameterAsEnum(parameters, "RESAMPLING", context)],
                              clip_bounds=clip_bounds, clip_crs=clip_crs,
                              buffer_m=self.parameterAsDouble(parameters, "BUFFER", context),
                              nodata_override=nodata,
                              source_text=(self.parameterAsString(parameters, "SOURCE", context) or "").strip(),
                              log_path=log_path,
                              progress=(lambda f: feedback.setProgress(100.0 * f)) if feedback else None)
        except (ValueError, RuntimeError) as e:
            raise QgsProcessingException(str(e))
        if log.get("zone_note"):
            feedback.pushWarning(log["zone_note"])
        for w in log["audit_input"]["warnings"]:
            feedback.pushWarning("Input: " + w)
        if log["input_damaged_upstream"]:
            feedback.pushWarning(log["advice"])
        for w in log["audit_output"]["warnings"]:
            feedback.pushWarning("Output: " + w)
        oi = log["output_info"]
        summary = (f"{oi['cols']} x {oi['rows']} cells of {log['cell_m']:g} m in {log['target_crs']}, "
                   f"{log['resampling']}; input {log['input_crs']} "
                   f"({log['input_cell_m'][0]:.1f} x {log['input_cell_m'][1]:.1f} m)")
        feedback.pushInfo(summary)
        feedback.pushInfo(f"Log: {log_path}")
        return {"OUTPUT": out, "LOG": log_path, "SUMMARY": summary}
