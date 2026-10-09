# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Catchment and longest-flow-path characteristics in one pass."""

from qgis.core import (
    QgsProcessingParameterRasterLayer, QgsProcessingParameterFeatureSource,
    QgsProcessingParameterNumber, QgsProcessingParameterBoolean,
    QgsProcessingParameterVectorDestination, QgsProcessing,
    QgsProcessingException,
    QgsCoordinateReferenceSystem,
)
import numpy as np

from .base import QehtAlgorithm
from ..core.raster import read_dem
from ..core.grid import decode_d8
from ..core.watershed.delineate import extract_streams
from ..core.watershed.statistics import CATCHMENT_FIELDS, FLOWPATH_FIELDS
from ..core.interop.heas_exchange import build_exchange_records, ExchangeError
from ..core.interop.gpkg import wkb_multipolygon, wkb_linestring

FDR="FDR"; DEM="DEM"; RAW_DEM="RAW_DEM"; FAC="FAC"; POINTS="POINTS"
SNAP="SNAP"; SNAP_THRESHOLD="SNAP_THRESHOLD"; NESTED="NESTED"
CATCH_OUT="CATCH_OUT"; PATH_OUT="PATH_OUT"; ID_FIELD="ID_FIELD"; ID_PREFIX="ID_PREFIX"
FLAT_CHECK="FLAT_CHECK"; FLAT_TOL="FLAT_TOL"
SENSITIVITY_FIELDS=[("area_barnes_km2","float"),("area_toward_km2","float"),
                    ("flat_sensitivity_pct","float"),("flat_sensitive","int")]


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
            "<i>catch_slope_horn</i> is the average ground gradient over every cell, "
            "by Horn's 3x3 method - a standard slope algorithm used by most "
            "GIS raster toolsets - and is what runoff coefficient and curve number tables "
            "assume. <i>catch_relief_ratio</i> is relief divided by longest "
            "flow path length (the relief ratio), a basin-steepness index; Tc "
            "formulas use the flow-path slope (lfp_slope / lfp_slope_1085) from "
            "Longest flow path. They are not interchangeable; "
            "state which you used. The v0.8 names slope_mean and "
            "slope_relief_ratio are still written as aliases for this release.\n\n"
            "<b>outlet_uid</b> is the stable identifier shared by a catchment and "
            "its flow path: from the ID attribute you choose, or sequential with "
            "your prefix (X001 = most downstream). <i>outlet_id</i> is the "
            "pour-point feature id and is NOT stable.\n\n"
            "<b>Supply the RAW DEM for elevations.</b> Reported heights should "
            "be real ground, not the fill surface. Routing still uses the "
            "conditioned DEM. If you leave the raw DEM blank the conditioned "
            "one is used for both, and elevations inside filled depressions "
            "will read high.\n\n"
            "<i>lfp_slope_1085</i> is the slope between the points at 10% and 85% "
            "of the path length measured from the OUTLET (the conventional "
            "definition), so it excludes the flat bottom 10% and the steep top "
            "15%. Elevations are interpolated along the path; lfp_L10_m, "
            "lfp_L85_m, lfp_z10_m and lfp_z85_m are written so it can be checked "
            "by hand. QEHT 0.8.3 and earlier measured from the divide instead.\n\n"
            "Every polygon holds its catchment's full area, also where "
            "catchments overlap (fixed in 0.9: earlier versions could clip or "
            "drop the polygon of a crossing listed after a downstream one).\n\n"
            "<b>Flat-method check</b> (on by default, 0.13.1): the conditioned DEM "
            "is also routed with both flat methods (Barnes 2014 and toward lower "
            "terrain) and the contributing area at each outlet is reported under "
            "each (<i>area_barnes_km2</i>, <i>area_toward_km2</i>). Where they differ "
            "by more than the tolerance, <i>flat_sensitive</i> = 1: flats or exactly "
            "tied cells on the drainage line decide where the flow goes, not the "
            "terrain, so verify that catchment against mapped drainage or on site. "
            "Adds about two routing passes of run time.\n\n"
            "Flow paths are also split at the channel head (first cell reaching the "
            "snap-to-stream threshold): overland and channel lengths and slopes, the "
            "10-85 slope of the channel part, and sheet (up to the cap) and shallow "
            "overland lengths. Catchments also carry perimeter, form factor, "
            "elongation and circularity ratios, drainage density, stream frequency "
            "and highest Strahler order.\n\n"
            "A layer from 'Road crossing candidates' can be given as the pour points: "
            "the accepted (else recommended) candidates are used at their own outlet "
            "cells, without snapping.\n\n"
            "For one linked package, use 'Build design hydrology package', which writes the same "
            "values into one self-describing GeoPackage."
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
        from qgis.core import QgsProcessingParameterField, QgsProcessingParameterString
        self.addParameter(QgsProcessingParameterField(
            ID_FIELD,"ID attribute for outlet_uid (optional; blank = sequential)",
            parentLayerParameterName=POINTS,optional=True))
        self.addParameter(QgsProcessingParameterString(
            ID_PREFIX,"ID prefix (sequential IDs only; ignored with an ID attribute)",
            defaultValue="X",optional=True))
        self.addParameter(QgsProcessingParameterBoolean(
            FLAT_CHECK,"Flat-method check: area under both flat methods, flag differences",
            defaultValue=True))
        self.addParameter(QgsProcessingParameterNumber(
            FLAT_TOL,"Flat-method check tolerance (%)",QgsProcessingParameterNumber.Type.Double,
            defaultValue=10.0,minValue=0.0,maxValue=100.0))
        self.addParameter(QgsProcessingParameterNumber(
            "SHEET_CAP","Sheet-flow cap within the overland part of the flow path (m)",
            QgsProcessingParameterNumber.Type.Double,defaultValue=100.0,minValue=0.0))
        self.addParameter(QgsProcessingParameterVectorDestination(
            CATCH_OUT,"Catchments with characteristics"))
        self.addParameter(QgsProcessingParameterVectorDestination(
            PATH_OUT,"Longest flow paths with characteristics"))

    def processAlgorithm(self, parameters, context, feedback):
        from osgeo import ogr, osr
        import os

        fdr_path=self.raster_path(parameters,FDR,context)
        fac_path=self.raster_path(parameters,FAC,context)
        snap_radius=self.parameterAsInt(parameters,SNAP,context)
        snap_threshold=self.parameterAsDouble(parameters,SNAP_THRESHOLD,context)
        nested=self.parameterAsBool(parameters,NESTED,context)
        catch_out=self.parameterAsOutputLayer(parameters,CATCH_OUT,context)
        path_out=self.parameterAsOutputLayer(parameters,PATH_OUT,context)
        id_field = self.field_parameter(parameters, ID_FIELD, context)
        prefix=(self.parameterAsString(parameters,ID_PREFIX,context) or "").strip()
        if id_field:
            prefix=""   # the prefix applies to sequential IDs only (0.13.1)

        d8,valid,info=read_dem(fdr_path)
        accum,_,_=read_dem(fac_path)
        direction=decode_d8(d8.astype(np.int32))

        if self.parameterAsRasterLayer(parameters,RAW_DEM,context) is not None:
            elevation,ev,_=read_dem(self.raster_path(parameters,RAW_DEM,context))
            feedback.pushInfo("Reporting elevations from the raw DEM.")
        else:
            elevation,ev,_=read_dem(self.raster_path(parameters,DEM,context))
            feedback.pushWarning(
                "No raw DEM supplied - elevations are taken from the conditioned "
                "DEM and will read high inside filled depressions.")
        elevation=np.where(ev,elevation,np.nan)   # NoData never becomes a height

        dem_crs=QgsCoordinateReferenceSystem(); dem_crs.createFromWkt(info.projection_wkt)
        if dem_crs.isValid() and dem_crs.isGeographic():
            feedback.pushWarning(
                "The DEM is in a geographic CRS: areas, lengths and slopes will be "
                "in degree units and are not meaningful. Reproject to a projected CRS.")

        stream_mask=None
        if snap_threshold>0:
            stream_mask=extract_streams(accum,valid,threshold_cells=snap_threshold)
        from ..core.network.crossings import CANDIDATE_FIELDS
        points=self.read_pour_points(parameters,POINTS,context,info,feedback,id_field=id_field,
                                     extra_fields=[f for f,_ in CANDIDATE_FIELDS])
        points,_,_=self.candidate_selection(points, feedback, info.cell_width)

        try:
            crossings,catchments,flowpaths,issues,_=build_exchange_records(
                direction,valid,accum,elevation,info.geotransform,points,
                snap_radius_cells=snap_radius,stream_mask=stream_mask,local=nested,
                id_scheme="attribute" if id_field else "sequential",id_prefix=prefix,
                progress=self.make_progress(feedback,weight=0.9),
                channel_threshold_cells=snap_threshold if snap_threshold>0 else None,
                sheet_cap_m=self.parameterAsDouble(parameters,"SHEET_CAP",context)
                if "SHEET_CAP" in parameters else 100.0)
        except ExchangeError as e:
            raise QgsProcessingException(str(e))
        for msg in issues: feedback.pushWarning(msg)

        srs=None
        if info.projection_wkt:
            srs=osr.SpatialReference(); srs.ImportFromWkt(info.projection_wkt)
        from ..core.raster import ogr_driver_for
        otype={"int":ogr.OFTInteger,"float":ogr.OFTReal,"text":ogr.OFTString}

        def write(path,lname,gtype,fields,rows,to_wkb):
            drv=ogr_driver_for(path)
            if os.path.exists(path): drv.DeleteDataSource(path)
            ds=drv.CreateDataSource(path)
            layer=ds.CreateLayer(lname,srs=srs,geom_type=gtype)
            for fname,ftype in fields:
                layer.CreateField(ogr.FieldDefn(fname,otype[ftype]))
            defn=layer.GetLayerDefn()
            for geom,attrs in rows:
                feat=ogr.Feature(defn)
                for i,(fname,ftype) in enumerate(fields):   # by index (.shp truncates names)
                    val=attrs.get(fname)
                    if val is None or (ftype!="text" and not np.isfinite(float(val))):
                        feat.SetFieldNull(i)       # never write NaN as 0
                    elif ftype=="int": feat.SetField(i,int(val))
                    elif ftype=="float": feat.SetField(i,float(val))
                    else: feat.SetField(i,str(val))
                feat.SetGeometry(ogr.CreateGeometryFromWkb(to_wkb(geom)))
                layer.CreateFeature(feat); feat=None
            ds=None

        catch_fields=list(CATCHMENT_FIELDS)
        if FLAT_CHECK not in parameters or self.parameterAsBool(parameters,FLAT_CHECK,context):
            from ..core.flow.sensitivity import flat_method_sensitivity
            tol=self.parameterAsDouble(parameters,FLAT_TOL,context) if FLAT_TOL in parameters else 10.0
            cond,cv,_=read_dem(self.raster_path(parameters,DEM,context))
            outlets=[info.xy_to_rowcol(cr["outlet_x"],cr["outlet_y"]) for _,cr in crossings]
            by_uid={cr["outlet_uid"]:res for (_,cr),res in zip(crossings,flat_method_sensitivity(
                cond,cv,info.cell_width,info.cell_height,outlets,tolerance=tol/100.0))}
            for _,ch in catchments:
                ch.update(by_uid.get(ch["outlet_uid"],{}))
            catch_fields+=SENSITIVITY_FIELDS
            flagged=[u for u,r in by_uid.items() if r["flat_sensitive"]]
            if flagged:
                feedback.pushWarning(
                    f"Flat-method check: {len(flagged)} outlet(s) change area by more than "
                    f"{tol:g}% between Barnes and toward-lower routing: "
                    + ", ".join(f"{u} ({by_uid[u]['area_barnes_km2']:.3f} vs "
                                f"{by_uid[u]['area_toward_km2']:.3f} km2)" for u in flagged)
                    + ". Flats or tied cells decide these catchments - verify against "
                    "mapped drainage or on site (see also the DEM QA check in Fill).")
            else:
                feedback.pushInfo(f"Flat-method check: every outlet within {tol:g}% under both methods.")
        write(catch_out,"catchments",ogr.wkbMultiPolygon,catch_fields,catchments,wkb_multipolygon)
        write(path_out,"longest_flow_paths",ogr.wkbLineString,FLOWPATH_FIELDS,flowpaths,wkb_linestring)

        for _,ch in catchments:
            feedback.pushInfo(
                f"  {ch['outlet_uid']} (fid {ch['outlet_id']}): A={ch['area_km2']:.4f} km2  "
                f"Hmax={ch['elev_max_m']:.1f}  Hmin={ch['elev_min_m']:.1f}  "
                f"Shorn={ch['catch_slope_horn']:.5f}  RR={ch.get('catch_relief_ratio',float('nan')):.5f}  "
                f"L={ch.get('lfp_length_km',float('nan')):.3f} km  "
                f"S1085={ch.get('lfp_slope_1085',float('nan')):.5f}")
        from ..core.interop.heas_exchange import plugin_version
        feedback.pushInfo(f"Wrote {len(catchments)} catchments and {len(flowpaths)} flow paths.")
        feedback.pushInfo(f"QEHT {plugin_version()} \u00b7 flat method: as used for the flow "
                          "direction raster (see the D8 flow direction log; default Barnes "
                          "since 0.13) \u00b7 10-85 reference outlet")
        return {CATCH_OUT:catch_out, PATH_OUT:path_out}
