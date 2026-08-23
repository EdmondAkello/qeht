# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Catchment delineation from pour points, with snapping."""

from qgis.core import (
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterVectorDestination,
    QgsProcessing,
    QgsProcessingException,
    QgsCoordinateTransform,
    QgsProject,
    QgsCoordinateReferenceSystem,
)
import numpy as np

from .base import QehtAlgorithm
from ..core.raster import read_dem, write_raster, polygonize
from ..core.grid import decode_d8
from ..core.watershed.delineate import (delineate_catchment, snap_pour_point,
                                         extract_streams)

FDR = "FDR"
FAC = "FAC"
POINTS = "POINTS"
SNAP = "SNAP"
SNAP_THRESHOLD = "SNAP_THRESHOLD"
RASTER_OUT = "RASTER_OUT"
POLY_OUT = "POLY_OUT"


class DelineateCatchmentAlgorithm(QehtAlgorithm):

    def name(self):
        return "delineatecatchment"

    def displayName(self):
        return "Delineate catchments from pour points"

    def shortHelpString(self):
        return (
            "Delineates the upstream contributing area for each pour point "
            "and writes both a labelled raster and dissolved polygons.\n\n"
            "<b>Snap radius</b> moves each pour point to the highest-"
            "accumulation cell within the given number of cells. A surveyed "
            "culvert or outfall coordinate almost never lands exactly on the "
            "DEM-derived channel; without snapping you get a catchment of a "
            "few cells and an absurd design flow. Set to 0 to disable.\n\n"
            "Where catchments nest, the FIRST point in the input wins, so "
            "order upstream points before downstream ones to carve out "
            "nested subcatchments.\n\n"
            "<b>Snap-to-stream threshold</b> controls the snapping strategy. "
            "With a threshold set, the point moves to the NEAREST stream cell, "
            "which preserves which tributary it sits on. With 0, it moves to "
            "the highest accumulation in the radius - which validation against "
            "a reference hydrology toolset showed pulls road crossings off small tributaries onto "
            "the adjacent trunk stream, inflating catchments by up to 40x. "
            "Keep the threshold non-zero.\n\n"
            "Polygons are produced by gdal.Polygonize() in-process. No "
            "external executable is launched at any point."
        )

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(FDR, "Flow direction (D8-coded)"))
        self.addParameter(QgsProcessingParameterRasterLayer(
            FAC, "Flow accumulation (for snapping)", optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource(
            POINTS, "Pour points", [QgsProcessing.TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterNumber(
            SNAP, "Snap radius (cells)", QgsProcessingParameterNumber.Integer,
            defaultValue=5, minValue=0, maxValue=100))
        self.addParameter(QgsProcessingParameterNumber(
            SNAP_THRESHOLD, "Snap-to-stream threshold (cells; 0 = snap to max accumulation)",
            QgsProcessingParameterNumber.Double, defaultValue=200.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterRasterDestination(
            RASTER_OUT, "Catchment raster"))
        self.addParameter(QgsProcessingParameterVectorDestination(
            POLY_OUT, "Catchment polygons"))

    def processAlgorithm(self, parameters, context, feedback):
        fdr_path = self.raster_path(parameters, FDR, context)
        source = self.parameterAsSource(parameters, POINTS, context)
        snap_radius = self.parameterAsInt(parameters, SNAP, context)
        snap_threshold = self.parameterAsDouble(parameters, SNAP_THRESHOLD, context)
        raster_out = self.parameterAsOutputLayer(parameters, RASTER_OUT, context)
        poly_out = self.parameterAsOutputLayer(parameters, POLY_OUT, context)

        d8, valid, info = read_dem(fdr_path)
        direction = decode_d8(d8.astype(np.int32))

        accum = None
        if snap_radius > 0:
            fac_layer = self.parameterAsRasterLayer(parameters, FAC, context)
            if fac_layer is None:
                feedback.pushWarning(
                    "Snapping requested but no accumulation raster supplied. "
                    "Pour points will be used exactly as given.")
                snap_radius = 0
            else:
                accum, _, _ = read_dem(self.raster_path(parameters, FAC, context))

        # Transform pour points into the DEM CRS.
        dem_crs = QgsCoordinateReferenceSystem()
        dem_crs.createFromWkt(info.projection_wkt)
        transform = None
        if dem_crs.isValid() and source.sourceCrs() != dem_crs:
            transform = QgsCoordinateTransform(source.sourceCrs(), dem_crs,
                                               QgsProject.instance())
            feedback.pushInfo(f"Reprojecting pour points "
                              f"{source.sourceCrs().authid()} -> {dem_crs.authid()}")

        stream_mask = None
        if accum is not None and snap_threshold > 0:
            stream_mask = extract_streams(accum, valid, threshold_cells=snap_threshold)
            feedback.pushInfo(
                f"Snapping to the stream network ({int(stream_mask.sum()):,} cells "
                f"at >= {snap_threshold:,.0f} accumulation).")
            feedback.pushInfo(
                "  Lower the threshold if a pour point sits on a small tributary: "
                "too high a threshold snaps it onto the trunk stream and inflates "
                "the catchment.")

        outlets = []
        for feature in source.getFeatures():
            geom = feature.geometry()
            if transform is not None:
                geom.transform(transform)
            pt = geom.asPoint()
            row, col = info.xy_to_rowcol(pt.x(), pt.y())
            if not (0 <= row < info.rows and 0 <= col < info.cols):
                feedback.pushWarning(f"Pour point {pt.x():.1f}, {pt.y():.1f} "
                                     "falls outside the DEM. Skipped.")
                continue
            if snap_radius > 0 and accum is not None:
                row, col, moved, acc_at = snap_pour_point(
                    row, col, accum, valid, search_radius_cells=snap_radius,
                    stream_mask=stream_mask)
                feedback.pushInfo(
                    f"  point -> cell ({row}, {col}), snapped {moved} cells, "
                    f"accumulation {acc_at:,.0f}")
            outlets.append((row, col))

        if not outlets:
            raise QgsProcessingException("No usable pour points.")

        feedback.setProgressText("Delineating catchments")
        labels = delineate_catchment(direction, valid, outlets)

        areas = {}
        for idx in range(1, len(outlets) + 1):
            n_cells = int((labels == idx).sum())
            areas[f"catchment {idx} area (map units squared)"] = n_cells * info.cell_area
        self.report_stats(feedback, "Catchments", areas)

        write_raster(raster_out, labels, info, dtype="int32", nodata=-1)
        feedback.setProgress(85)

        result = polygonize(labels, info, poly_out, layer_name="catchments",
                            field_name="DN", ignore_value=0, dissolve=True)
        feedback.pushInfo(f"Polygonized in-process: {result['features_raw']} raw parts "
                          f"-> {result['features']} dissolved catchment(s)")
        return {RASTER_OUT: raster_out, POLY_OUT: poly_out}
