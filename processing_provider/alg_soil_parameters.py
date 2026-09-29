# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Soil parameters for catchment polygons (WP-C)."""

from qgis.core import (
    QgsProcessingParameterRasterLayer, QgsProcessingParameterFeatureSource,
    QgsProcessingParameterFeatureSink, QgsProcessingParameterRasterDestination,
    QgsProcessing, QgsProcessingException, QgsFeatureSink, QgsFields, QgsFeature,
    QgsCoordinateReferenceSystem, QgsCoordinateTransform, QgsProject, QgsGeometry,
)
import numpy as np

from .base import QehtAlgorithm
from ..core.raster import read_dem, write_raster
from ..core.geometry.rasterize import rasterize_polygons
from ..core.soils.catchment import soil_block, SOIL_FIELDS

REF = "REF"; CATCHMENTS = "CATCHMENTS"; OUTPUT = "OUTPUT"
K_RASTER = "K_RASTER"; UNIT_RASTER = "UNIT_RASTER"


class SoilParametersAlgorithm(QehtAlgorithm):

    def name(self): return "soilparameters"
    def displayName(self): return "Soil parameters for catchments"
    def group(self): return "Soils and erosion"
    def groupId(self): return "soils"

    def shortHelpString(self):
        return (
            "Adds a soil block to catchment polygons (any polygon layer - e.g. the "
            "catchments of an exchange package): area-weighted topsoil sand, silt, "
            "clay, organic carbon, coarse fragments and bulk density; USDA texture; "
            "dominant FAO drainage class; USLE K (Williams/EPIC, SI units, plus the "
            "Renard Dg-based alternative); the coarse-fragment factor CFRG; a "
            "texture-based hydrologic soil group PROXY; and the share of the "
            "catchment covered by soil data.\n\n"
            "<b>SOTWIS:</b> give the polygons (KEN_SOTWISv1_t1s1d1) and the SQLite "
            "database (KEN_SOTWISv1.db - never the .mdb). Each soil unit's "
            "components (up to 10 profiles with their shares) are depth-weighted "
            "over the chosen interval (default 0-20 cm) and K is computed per "
            "component, then weighted. Without the database the polygons' own "
            "attributes are used (dominant soil only).\n\n"
            "The soil map is rasterised on the reference grid (cell-centre rule), "
            "so every catchment is weighted by cells exactly like the hydrology. "
            "Optional rasters: soil unit index and USLE K per cell.\n\n"
            "<i>soil_hsg_proxy</i> is derived from texture and drainage, not measured "
            "infiltration - label it as such in a report.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFeatureSource(
            CATCHMENTS, "Catchment polygons", [QgsProcessing.SourceType.TypeVectorPolygon]))
        self.addParameter(QgsProcessingParameterRasterLayer(
            REF, "Reference grid (the DEM or flow-direction raster used for the catchments)"))
        self.add_soil_parameters("SOIL", optional=False)
        self.addParameter(QgsProcessingParameterFeatureSink(
            OUTPUT, "Catchments with soil parameters", QgsProcessing.SourceType.TypeVectorPolygon))
        self.addParameter(QgsProcessingParameterRasterDestination(
            K_RASTER, "USLE K raster (optional)", optional=True, createByDefault=False))
        self.addParameter(QgsProcessingParameterRasterDestination(
            UNIT_RASTER, "Soil unit index raster (optional)", optional=True, createByDefault=False))

    def processAlgorithm(self, parameters, context, feedback):
        _, _, info = read_dem(self.raster_path(parameters, REF, context))
        soil = self.load_soil(parameters, context, info, feedback, "SOIL")
        grid, idx, units, sinfo = soil
        src = self.parameterAsSource(parameters, CATCHMENTS, context)
        fields = QgsFields(src.fields())
        existing = set(fields.names())
        for name, kind in SOIL_FIELDS:
            if name not in existing:
                fields.append(self.qgs_field(name, kind))
        sink, dest = self.parameterAsSink(parameters, OUTPUT, context, fields,
                                          src.wkbType(), src.sourceCrs())
        dem_crs = QgsCoordinateReferenceSystem(); dem_crs.createFromWkt(info.projection_wkt)
        tr = None
        if dem_crs.isValid() and src.sourceCrs() != dem_crs:
            tr = QgsCoordinateTransform(src.sourceCrs(), dem_crs, QgsProject.instance())
        n = max(1, src.featureCount())
        for k, feat in enumerate(src.getFeatures()):
            if feedback.isCanceled():
                break
            g = QgsGeometry(feat.geometry())
            if tr is not None:
                g.transform(tr)
            parts = g.asMultiPolygon() if g.isMultipart() else [g.asPolygon()]
            mask = rasterize_polygons(
                [([[(p.x(), p.y()) for p in ring] for ring in part], 1) for part in parts],
                info.geotransform, (info.rows, info.cols)) > 0
            block = soil_block(mask, grid, units, idx, sinfo)
            out = QgsFeature(fields)
            out.setGeometry(feat.geometry())
            attrs = dict(zip(src.fields().names(), feat.attributes()))
            for name in fields.names():
                if name in block:
                    v = block[name]
                    attrs[name] = None if (v is None or (isinstance(v, float) and not np.isfinite(v))) else v
            out.setAttributes([attrs.get(name) for name in fields.names()])
            sink.addFeature(out, QgsFeatureSink.Flag.FastInsert)
            feedback.pushInfo(f"  feature {feat.id()}: K={block['usle_k'] if block['usle_k'] is None else round(block['usle_k'], 5)}  "
                              f"texture={block['soil_texture']}  coverage={block['soil_coverage_pct']}")
            feedback.setProgress(int(100.0 * (k + 1) / n))
        result = {OUTPUT: dest}
        k_out = self.parameterAsOutputLayer(parameters, K_RASTER, context)
        if k_out:
            kgrid = np.full(grid.shape, np.nan)
            for i, uid in idx.items():
                u = units.get(uid)
                if u is not None and u.k_si is not None and np.isfinite(u.k_si):
                    kgrid[grid == i] = u.k_si
            write_raster(k_out, kgrid, info, nodata=-9999.0)
            result[K_RASTER] = k_out
        u_out = self.parameterAsOutputLayer(parameters, UNIT_RASTER, context)
        if u_out:
            write_raster(u_out, grid, info, dtype="int32", nodata=0)
            result[UNIT_RASTER] = u_out
        return result
