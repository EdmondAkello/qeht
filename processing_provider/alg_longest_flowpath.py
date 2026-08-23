# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Longest flow path - the length term for time of concentration."""

from qgis.core import (
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterNumber,
    QgsProcessingParameterVectorDestination,
    QgsProcessing,
    QgsProcessingException,
    QgsCoordinateTransform,
    QgsProject,
    QgsCoordinateReferenceSystem,
)
import numpy as np

from .base import QehtAlgorithm
from ..core.raster import read_dem, polyline_from_cells
from ..core.grid import decode_d8
from ..core.watershed.delineate import (delineate_catchment, snap_pour_point,
                                        longest_flow_path, extract_streams,
                                        longest_flow_paths_by_catchment)

FDR = "FDR"
DEM = "DEM"
FAC = "FAC"
POINTS = "POINTS"
SNAP = "SNAP"
SNAP_THRESHOLD = "SNAP_THRESHOLD"
NESTED = "NESTED"
OUTPUT = "OUTPUT"


class LongestFlowPathAlgorithm(QehtAlgorithm):

    def name(self):
        return "longestflowpath"

    def displayName(self):
        return "Longest flow path"

    def shortHelpString(self):
        return (
            "Traces the longest flow path from the catchment divide to each "
            "pour point, and reports its length, elevation drop and average "
            "slope.\n\n"
            "This is one of the harder tools to find a free/open equivalent for, and "
            "the one that matters most for drainage design: length and slope "
            "along this path are the inputs to Kirpich, Bransby-Williams and "
            "the TRRL time-of-concentration methods.\n\n"
            "Length is accumulated per link using true diagonal weighting "
            "(cellsize x sqrt(2)), not a cell count, so the result is a real "
            "planimetric distance.\n\n"
            "One path is written PER POUR POINT, each bounded by that "
            "outlet's own catchment. Use 'non-overlapping' when you want "
            "local catchments rather than full upstream areas.\n\n"
            "Note this is the planimetric length. If you need the slope "
            "along the channel profile rather than the average divide-to-"
            "outlet gradient, extract the profile from the output line."
        )

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(FDR, "Flow direction (D8-coded)"))
        self.addParameter(QgsProcessingParameterRasterLayer(DEM, "Conditioned DEM (for slope)"))
        self.addParameter(QgsProcessingParameterRasterLayer(
            FAC, "Flow accumulation (for snapping)", optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource(
            POINTS, "Pour points", [QgsProcessing.TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterNumber(
            SNAP, "Snap radius (cells)", QgsProcessingParameterNumber.Integer,
            defaultValue=5, minValue=0, maxValue=100))
        self.addParameter(QgsProcessingParameterNumber(
            SNAP_THRESHOLD, "Snap-to-stream threshold (cells; 0 = max accumulation)",
            QgsProcessingParameterNumber.Double, defaultValue=200.0, minValue=0.0))
        from qgis.core import QgsProcessingParameterBoolean
        self.addParameter(QgsProcessingParameterBoolean(
            NESTED, "Non-overlapping (local) catchments", defaultValue=False))
        self.addParameter(QgsProcessingParameterVectorDestination(
            OUTPUT, "Longest flow paths"))

    def processAlgorithm(self, parameters, context, feedback):
        fdr_path = self.raster_path(parameters, FDR, context)
        dem_path = self.raster_path(parameters, DEM, context)
        source = self.parameterAsSource(parameters, POINTS, context)
        snap_radius = self.parameterAsInt(parameters, SNAP, context)
        out_path = self.parameterAsOutputLayer(parameters, OUTPUT, context)

        d8, valid, info = read_dem(fdr_path)
        elevation, _, _ = read_dem(dem_path)
        direction = decode_d8(d8.astype(np.int32))

        snap_threshold = self.parameterAsDouble(parameters, SNAP_THRESHOLD, context)
        nested = self.parameterAsBool(parameters, NESTED, context)

        accum = None
        if snap_radius > 0 and self.parameterAsRasterLayer(parameters, FAC, context):
            accum, _, _ = read_dem(self.raster_path(parameters, FAC, context))

        stream_mask = None
        if accum is not None and snap_threshold > 0:
            stream_mask = extract_streams(accum, valid, threshold_cells=snap_threshold)

        dem_crs = QgsCoordinateReferenceSystem()
        dem_crs.createFromWkt(info.projection_wkt)
        transform = None
        if dem_crs.isValid() and source.sourceCrs() != dem_crs:
            transform = QgsCoordinateTransform(source.sourceCrs(), dem_crs,
                                               QgsProject.instance())

        outlets, attrs = [], []
        for feature in source.getFeatures():
            geom = feature.geometry()
            if transform is not None:
                geom.transform(transform)
            if geom.isMultipart():
                pt = geom.asMultiPoint()[0]
            else:
                pt = geom.asPoint()
            row, col = info.xy_to_rowcol(pt.x(), pt.y())
            if not (0 <= row < info.rows and 0 <= col < info.cols):
                feedback.pushWarning("A pour point falls outside the DEM. Skipped.")
                continue
            if accum is not None:
                row, col, moved, _ = snap_pour_point(
                    row, col, accum, valid, search_radius_cells=snap_radius,
                    stream_mask=stream_mask)
            outlets.append((row, col))
            attrs.append(feature.id())

        if not outlets:
            raise QgsProcessingException("No usable pour points.")

        feedback.pushInfo(f"Tracing {len(outlets)} longest flow paths"
                          f"{' (local catchments)' if nested else ' (full upstream areas)'}")
        paths = longest_flow_paths_by_catchment(
            direction, valid, outlets, elevation=elevation,
            cell_width=info.cell_width, cell_height=info.cell_height,
            nested=nested, progress=self.make_progress(feedback, weight=0.9))

        from osgeo import ogr, osr
        import os
        srs = None
        if info.projection_wkt:
            srs = osr.SpatialReference(); srs.ImportFromWkt(info.projection_wkt)
        drv = ogr.GetDriverByName("GPKG")
        if os.path.exists(out_path):
            drv.DeleteDataSource(out_path)
        vds = drv.CreateDataSource(out_path)
        layer = vds.CreateLayer("longest_flow_paths", srs=srs,
                                geom_type=ogr.wkbLineString)
        for fname, ftype in [("outlet_id", ogr.OFTInteger), ("length", ogr.OFTReal),
                             ("drop", ogr.OFTReal), ("slope", ogr.OFTReal),
                             ("area", ogr.OFTReal)]:
            layer.CreateField(ogr.FieldDefn(fname, ftype))
        defn = layer.GetLayerDefn()

        written = 0
        for lfp, fid in zip(paths, attrs):
            if len(lfp["cells"]) < 2:
                feedback.pushWarning(f"Outlet {fid}: no upstream path. Skipped.")
                continue
            line = ogr.Geometry(ogr.wkbLineString)
            for row, col in lfp["cells"]:
                x, y = info.rowcol_to_xy(row, col)
                line.AddPoint_2D(x, y)
            feat = ogr.Feature(defn)
            feat.SetField("outlet_id", int(fid))
            feat.SetField("length", float(lfp["length"]))
            feat.SetField("drop", float(lfp.get("drop", 0.0)))
            feat.SetField("slope", float(lfp.get("slope", 0.0)))
            feat.SetField("area", float(lfp.get("catchment_area", 0.0)))
            feat.SetGeometry(line)
            layer.CreateFeature(feat); feat = None
            written += 1
            feedback.pushInfo(
                f"  outlet {fid}: L={lfp['length']:,.1f} m  "
                f"S={lfp.get('slope', 0):.5f}  A={lfp.get('catchment_area', 0)/1e6:,.3f} km2")
        vds = None
        feedback.pushInfo(f"Wrote {written} longest flow paths.")
        return {OUTPUT: out_path}
