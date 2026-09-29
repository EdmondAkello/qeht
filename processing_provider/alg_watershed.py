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
ID_FIELD = "ID_FIELD"
ID_PREFIX = "ID_PREFIX"
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
            POINTS, "Pour points", [QgsProcessing.SourceType.TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterNumber(
            SNAP, "Snap radius (cells)", QgsProcessingParameterNumber.Type.Integer,
            defaultValue=5, minValue=0, maxValue=100))
        self.addParameter(QgsProcessingParameterNumber(
            SNAP_THRESHOLD, "Snap-to-stream threshold (cells; 0 = snap to max accumulation)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=200.0, minValue=0.0))
        from qgis.core import QgsProcessingParameterField, QgsProcessingParameterString
        self.addParameter(QgsProcessingParameterField(
            ID_FIELD, "ID attribute for outlet_uid (optional; blank = sequential)",
            parentLayerParameterName=POINTS, optional=True))
        self.addParameter(QgsProcessingParameterString(
            ID_PREFIX, "ID prefix (sequential IDs only; ignored with an ID attribute)",
            defaultValue="X", optional=True))
        self.addParameter(QgsProcessingParameterRasterDestination(
            RASTER_OUT, "Catchment raster"))
        self.addParameter(QgsProcessingParameterVectorDestination(
            POLY_OUT, "Catchment polygons"))

    def processAlgorithm(self, parameters, context, feedback):
        fdr_path = self.raster_path(parameters, FDR, context)
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

        id_field = self.field_parameter(parameters, ID_FIELD, context)
        prefix = (self.parameterAsString(parameters, ID_PREFIX, context) or "").strip()
        points = self.read_pour_points(parameters, POINTS, context, info, feedback,
                                       id_field=id_field)
        outlets = []
        for p in points:
            row, col = info.xy_to_rowcol(p["x"], p["y"])
            if snap_radius > 0 and accum is not None:
                row, col, moved, acc_at = snap_pour_point(
                    row, col, accum, valid, search_radius_cells=snap_radius,
                    stream_mask=stream_mask)
                feedback.pushInfo(
                    f"  point fid {p['fid']} -> cell ({row}, {col}), snapped {moved} cells, "
                    f"accumulation {acc_at:,.0f}")
            outlets.append((row, col))
        uids = self.outlet_uids(points, outlets, accum, id_field, prefix)

        if not outlets:
            raise QgsProcessingException("No usable pour points.")

        feedback.setProgressText("Delineating catchments")
        labels = delineate_catchment(direction, valid, outlets)

        areas = {}
        for idx in range(1, len(outlets) + 1):
            n_cells = int((labels == idx).sum())
            areas[f"{uids[idx - 1]} area (map units squared)"] = n_cells * info.cell_area
        self.report_stats(feedback, "Catchments", areas)

        write_raster(raster_out, labels, info, dtype="int32", nodata=-1)
        feedback.setProgress(85)

        result = polygonize(labels, info, poly_out, layer_name="catchments",
                            field_name="DN", ignore_value=0, dissolve=True)
        feedback.pushInfo(f"Polygonized in-process: {result['features_raw']} raw parts "
                          f"-> {result['features']} dissolved catchment(s)")

        # outlet_uid on the polygons (DN = position in the pour-point list)
        from osgeo import ogr
        ds = ogr.Open(poly_out, 1)
        if ds is not None:
            layer = ds.GetLayer(0)
            layer.CreateField(ogr.FieldDefn("outlet_uid", ogr.OFTString))
            layer.ResetReading()
            for feat in layer:
                dn = feat.GetField("DN")
                if dn is not None and 1 <= int(dn) <= len(uids):
                    feat.SetField("outlet_uid", uids[int(dn) - 1])
                    layer.SetFeature(feat)
            ds = None
        return {RASTER_OUT: raster_out, POLY_OUT: poly_out}
