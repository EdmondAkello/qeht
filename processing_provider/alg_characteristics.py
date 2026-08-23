# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Catchment and longest-flow-path characteristics in one pass."""

from qgis.core import (
    QgsProcessingParameterRasterLayer, QgsProcessingParameterFeatureSource,
    QgsProcessingParameterNumber, QgsProcessingParameterBoolean,
    QgsProcessingParameterVectorDestination, QgsProcessing,
    QgsProcessingException, QgsCoordinateTransform, QgsProject,
    QgsCoordinateReferenceSystem,
)
import numpy as np

from .base import QehtAlgorithm
from ..core.raster import read_dem, polygonize
from ..core.grid import decode_d8
from ..core.watershed.delineate import (delineate_catchment, snap_pour_point,
                                        longest_flow_path, extract_streams)
from ..core.watershed.statistics import (catchment_characteristics, horn_slope,
                                         CATCHMENT_FIELDS, FLOWPATH_FIELDS)

FDR="FDR"; DEM="DEM"; RAW_DEM="RAW_DEM"; FAC="FAC"; POINTS="POINTS"
SNAP="SNAP"; SNAP_THRESHOLD="SNAP_THRESHOLD"; NESTED="NESTED"
CATCH_OUT="CATCH_OUT"; PATH_OUT="PATH_OUT"


class CatchmentCharacteristicsAlgorithm(QehtAlgorithm):

    def name(self): return "catchmentcharacteristics"
    def displayName(self): return "Catchment and flow path characteristics"

    def shortHelpString(self):
        return (
            "Delineates a catchment and its longest flow path for each pour "
            "point, and writes the morphometry needed for a design flood "
            "calculation.\n\n"
            "<b>Catchment polygons carry:</b> area (km2), highest and lowest "
            "elevation (m), relief, mean terrain slope, relief ratio and the "
            "longest flow path length.\n\n"
            "<b>Flow path lines carry:</b> length (km), highest and lowest "
            "elevation (m), drop, slope and the 10-85 slope.\n\n"
            "<b>Two different catchment slopes are reported, deliberately.</b> "
            "<i>slope_mean</i> is the average ground gradient over every cell, "
            "by Horn's 3x3 method - a standard slope algorithm used by most "
            "GIS raster toolsets - and is what runoff coefficient and curve number tables "
            "assume. <i>slope_relief_ratio</i> is relief divided by longest "
            "flow path length, which is the 'catchment slope' of most road "
            "drainage manuals and the term that goes into Kirpich. They are "
            "not interchangeable; state which you used.\n\n"
            "<b>Supply the RAW DEM for elevations.</b> Reported heights should "
            "be real ground, not the fill surface. Routing still uses the "
            "conditioned DEM. If you leave the raw DEM blank the conditioned "
            "one is used for both, and elevations inside filled depressions "
            "will read high.\n\n"
            "<i>lfp_slope_1085</i> excludes the top 10% and bottom 15% of the "
            "path, removing the steep headwater and flat outlet reach that "
            "distort a whole-path average. It is required by several UK and "
            "TRRL methods and is usually the more defensible design figure."
        )

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(FDR,"Flow direction (D8-coded)"))
        self.addParameter(QgsProcessingParameterRasterLayer(DEM,"Conditioned DEM"))
        self.addParameter(QgsProcessingParameterRasterLayer(
            RAW_DEM,"Raw DEM (for reported elevations)",optional=True))
        self.addParameter(QgsProcessingParameterRasterLayer(FAC,"Flow accumulation"))
        self.addParameter(QgsProcessingParameterFeatureSource(
            POINTS,"Pour points",[QgsProcessing.SourceType.TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterNumber(
            SNAP,"Snap radius (cells)",QgsProcessingParameterNumber.Type.Integer,
            defaultValue=5,minValue=0,maxValue=100))
        self.addParameter(QgsProcessingParameterNumber(
            SNAP_THRESHOLD,"Snap-to-stream threshold (cells; 0 = max accumulation)",
            QgsProcessingParameterNumber.Type.Double,defaultValue=200.0,minValue=0.0))
        self.addParameter(QgsProcessingParameterBoolean(
            NESTED,"Non-overlapping (local) catchments",defaultValue=False))
        self.addParameter(QgsProcessingParameterVectorDestination(
            CATCH_OUT,"Catchments with characteristics"))
        self.addParameter(QgsProcessingParameterVectorDestination(
            PATH_OUT,"Longest flow paths with characteristics"))

    def processAlgorithm(self, parameters, context, feedback):
        from osgeo import ogr, osr
        import os

        fdr_path=self.raster_path(parameters,FDR,context)
        dem_path=self.raster_path(parameters,DEM,context)
        fac_path=self.raster_path(parameters,FAC,context)
        source=self.parameterAsSource(parameters,POINTS,context)
        snap_radius=self.parameterAsInt(parameters,SNAP,context)
        snap_threshold=self.parameterAsDouble(parameters,SNAP_THRESHOLD,context)
        nested=self.parameterAsBool(parameters,NESTED,context)
        catch_out=self.parameterAsOutputLayer(parameters,CATCH_OUT,context)
        path_out=self.parameterAsOutputLayer(parameters,PATH_OUT,context)

        d8,valid,info=read_dem(fdr_path)
        conditioned,_,_=read_dem(dem_path)
        accum,_,_=read_dem(fac_path)
        direction=decode_d8(d8.astype(np.int32))

        if self.parameterAsRasterLayer(parameters,RAW_DEM,context) is not None:
            elevation,_,_=read_dem(self.raster_path(parameters,RAW_DEM,context))
            feedback.pushInfo("Reporting elevations from the raw DEM.")
        else:
            elevation=conditioned
            feedback.pushWarning(
                "No raw DEM supplied - elevations are taken from the conditioned "
                "DEM and will read high inside filled depressions.")

        feedback.setProgressText("Computing terrain slope (Horn 3x3)")
        slope_raster=horn_slope(elevation,valid,info.cell_width,info.cell_height)

        stream_mask=None
        if snap_threshold>0:
            stream_mask=extract_streams(accum,valid,threshold_cells=snap_threshold)

        dem_crs=QgsCoordinateReferenceSystem(); dem_crs.createFromWkt(info.projection_wkt)
        transform=None
        if dem_crs.isValid() and source.sourceCrs()!=dem_crs:
            transform=QgsCoordinateTransform(source.sourceCrs(),dem_crs,QgsProject.instance())

        outlets,fids=[],[]
        for feature in source.getFeatures():
            geom=feature.geometry()
            if transform is not None: geom.transform(transform)
            pt=geom.asMultiPoint()[0] if geom.isMultipart() else geom.asPoint()
            row,col=info.xy_to_rowcol(pt.x(),pt.y())
            if not (0<=row<info.rows and 0<=col<info.cols):
                feedback.pushWarning("A pour point falls outside the DEM. Skipped."); continue
            if snap_radius>0:
                row,col,_,_=snap_pour_point(row,col,accum,valid,
                                            search_radius_cells=snap_radius,
                                            stream_mask=stream_mask)
            outlets.append((row,col)); fids.append(feature.id())
        if not outlets:
            raise QgsProcessingException("No usable pour points.")

        labels=delineate_catchment(direction,valid,outlets) if nested else None

        srs=None
        if info.projection_wkt:
            srs=osr.SpatialReference(); srs.ImportFromWkt(info.projection_wkt)
        drv=ogr.GetDriverByName("GPKG")
        for p in (catch_out,path_out):
            if os.path.exists(p): drv.DeleteDataSource(p)
        cds=drv.CreateDataSource(catch_out)
        clayer=cds.CreateLayer("catchments",srs=srs,geom_type=ogr.wkbPolygon)
        pds=drv.CreateDataSource(path_out)
        player=pds.CreateLayer("longest_flow_paths",srs=srs,geom_type=ogr.wkbLineString)
        for layer,fields in ((clayer,CATCHMENT_FIELDS),(player,FLOWPATH_FIELDS)):
            for fname,ftype in fields:
                layer.CreateField(ogr.FieldDefn(
                    fname, ogr.OFTInteger if ftype=="int" else ogr.OFTReal))

        label_grid=np.zeros(direction.shape,dtype=np.int32)
        stats_rows=[]
        for i,((row,col),fid) in enumerate(zip(outlets,fids)):
            feedback.setProgress(int(100.0*i/len(outlets)))
            mask=(labels==(i+1)) if nested else (delineate_catchment(direction,valid,[(row,col)])>0)
            if not mask.any():
                feedback.pushWarning(f"Outlet {fid}: empty catchment. Skipped."); continue

            lfp=longest_flow_path(direction,valid,(row,col),elevation=elevation,
                                  cell_width=info.cell_width,cell_height=info.cell_height,
                                  catchment_mask=mask)
            ch=catchment_characteristics(mask,elevation,valid,info.cell_width,
                                         info.cell_height,flow_path=lfp,
                                         slope_raster=slope_raster)
            if not ch: continue
            ch["outlet_id"]=int(fid)
            stats_rows.append(ch)
            label_grid[mask & (label_grid==0)]=i+1

            if len(lfp["cells"])>=2:
                line=ogr.Geometry(ogr.wkbLineString)
                for r,c in lfp["cells"]:
                    x,y=info.rowcol_to_xy(r,c); line.AddPoint_2D(x,y)
                feat=ogr.Feature(player.GetLayerDefn())
                for fname,ftype in FLOWPATH_FIELDS:
                    val=ch.get(fname,0)
                    feat.SetField(fname, int(val) if ftype=="int"
                                  else (float(val) if np.isfinite(val) else 0.0))
                feat.SetGeometry(line); player.CreateFeature(feat); feat=None

            feedback.pushInfo(
                f"  outlet {fid}: A={ch['area_km2']:.4f} km2  "
                f"Hmax={ch['elev_max_m']:.1f}  Hmin={ch['elev_min_m']:.1f}  "
                f"Smean={ch['slope_mean']:.5f}  Srelief={ch.get('slope_relief_ratio',0):.5f}  "
                f"L={ch.get('lfp_length_km',0):.3f} km  Slfp={ch.get('lfp_slope',0):.5f}")
        pds=None

        # Polygonize the label grid, then join the statistics back on.
        tmp_poly=catch_out+".tmp.gpkg"
        polygonize(label_grid,info,tmp_poly,layer_name="catchments",
                   field_name="DN",ignore_value=0,dissolve=True)
        src_ds=ogr.Open(tmp_poly)
        src_layer=src_ds.GetLayer(0)
        by_index={i+1:s for i,s in enumerate(stats_rows)}
        defn=clayer.GetLayerDefn()
        for feat_in in src_layer:
            dn=feat_in.GetField("DN")
            ch=by_index.get(dn)
            if ch is None: continue
            feat=ogr.Feature(defn)
            for fname,ftype in CATCHMENT_FIELDS:
                val=ch.get(fname,0)
                feat.SetField(fname, int(val) if ftype=="int"
                              else (float(val) if np.isfinite(val) else 0.0))
            feat.SetGeometry(feat_in.GetGeometryRef().Clone())
            clayer.CreateFeature(feat); feat=None
        src_ds=None; cds=None
        try: os.remove(tmp_poly)
        except OSError: pass

        feedback.pushInfo(f"Wrote {len(stats_rows)} catchments and their flow paths.")
        return {CATCH_OUT:catch_out, PATH_OUT:path_out}
