# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Vectorise the stream network into polyline reaches."""

from qgis.core import (
    QgsProcessingParameterRasterLayer, QgsProcessingParameterEnum,
    QgsProcessingParameterNumber, QgsProcessingParameterVectorDestination,
    QgsProcessingException,
)
import numpy as np

from .base import QehtAlgorithm
from ..core.raster import read_dem
from ..core.grid import decode_d8
from ..core.watershed.delineate import extract_streams
from ..core.flow.accumulation import strahler_order
from ..core.flow.streamlines import vectorize_streams, network_summary

FDR="FDR"; FAC="FAC"; DEM="DEM"; MODE="MODE"; THRESHOLD="THRESHOLD"; OUTPUT="OUTPUT"
MODES=["Accumulation threshold (cells)","Contributing area threshold (map units squared)"]


class StreamLinesAlgorithm(QehtAlgorithm):

    def name(self): return "streamlines"
    def displayName(self): return "Stream network to polylines"

    def shortHelpString(self):
        return ("Extracts the channel network at a flow accumulation threshold "
                "and traces it into polyline REACHES - one feature per channel "
                "segment between sources, junctions and the outlet.\n\n"
                "Equivalent to the 'Drainage Line Processing' tool in commercial "
                "hydrology toolsets.\n\n"
                "Each reach carries its Strahler order, length, elevation drop "
                "and slope, so the output is directly usable for drainage "
                "density, channel slope and network morphometry without "
                "further processing.\n\n"
                "Lower the threshold to extend the network into smaller "
                "tributaries. If you are snapping pour points to this network, "
                "use the same threshold in both tools.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(FDR,"Flow direction (D8-coded)"))
        self.addParameter(QgsProcessingParameterRasterLayer(FAC,"Flow accumulation (cell count)"))
        self.addParameter(QgsProcessingParameterRasterLayer(DEM,"Conditioned DEM (for slope)",optional=True))
        self.addParameter(QgsProcessingParameterEnum(MODE,"Threshold type",options=MODES,defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber(
            THRESHOLD,"Threshold value",QgsProcessingParameterNumber.Type.Double,
            defaultValue=1000.0,minValue=1.0))
        self.addParameter(QgsProcessingParameterVectorDestination(OUTPUT,"Stream reaches"))

    def processAlgorithm(self, parameters, context, feedback):
        from osgeo import ogr, osr
        import os

        fdr_path=self.raster_path(parameters,FDR,context)
        fac_path=self.raster_path(parameters,FAC,context)
        mode=self.parameterAsEnum(parameters,MODE,context)
        threshold=self.parameterAsDouble(parameters,THRESHOLD,context)
        out_path=self.parameterAsOutputLayer(parameters,OUTPUT,context)

        d8,valid,info=read_dem(fdr_path)
        accum,_,_=read_dem(fac_path)
        if accum.shape!=d8.shape:
            raise QgsProcessingException("Accumulation and direction rasters must match.")
        direction=decode_d8(d8.astype(np.int32))

        elevation=None
        if self.parameterAsRasterLayer(parameters,DEM,context) is not None:
            elevation,_,_=read_dem(self.raster_path(parameters,DEM,context))

        if mode==0:
            mask=extract_streams(accum,valid,threshold_cells=threshold)
        else:
            mask=extract_streams(accum,valid,threshold_area=threshold,cell_area=info.cell_area)
        if not mask.any():
            raise QgsProcessingException("No cells exceed the threshold. Lower it.")

        feedback.setProgressText("Strahler ordering")
        order=strahler_order(direction,valid,mask)

        reaches=vectorize_streams(direction,valid,mask,
                                  cell_width=info.cell_width,cell_height=info.cell_height,
                                  elevation=elevation,order=order,
                                  progress=self.make_progress(feedback,weight=0.8))

        area=float(valid.sum())*info.cell_area
        self.report_stats(feedback,"Stream network",network_summary(reaches,catchment_area=area))

        srs=None
        if info.projection_wkt:
            srs=osr.SpatialReference(); srs.ImportFromWkt(info.projection_wkt)
        drv=ogr.GetDriverByName("GPKG")
        if os.path.exists(out_path): drv.DeleteDataSource(out_path)
        vds=drv.CreateDataSource(out_path)
        layer=vds.CreateLayer("stream_reaches",srs=srs,geom_type=ogr.wkbLineString)
        for fname,ftype in [("reach_id",ogr.OFTInteger),("order",ogr.OFTInteger),
                            ("length",ogr.OFTReal),("drop",ogr.OFTReal),
                            ("slope",ogr.OFTReal),("head_type",ogr.OFTString)]:
            layer.CreateField(ogr.FieldDefn(fname,ftype))
        defn=layer.GetLayerDefn()
        for r in reaches:
            line=ogr.Geometry(ogr.wkbLineString)
            for row,col in r["cells"]:
                x,y=info.rowcol_to_xy(row,col); line.AddPoint_2D(x,y)
            feat=ogr.Feature(defn)
            feat.SetField("reach_id",int(r["reach_id"]))
            feat.SetField("order",int(r.get("order",0)))
            feat.SetField("length",float(r["length"]))
            feat.SetField("drop",float(r.get("drop",0.0)))
            feat.SetField("slope",float(r.get("slope",0.0)))
            feat.SetField("head_type",str(r.get("head_type","")))
            feat.SetGeometry(line); layer.CreateFeature(feat); feat=None
        vds=None
        feedback.pushInfo(f"Wrote {len(reaches):,} reaches.")
        return {OUTPUT:out_path}
