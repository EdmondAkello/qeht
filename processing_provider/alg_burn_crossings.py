# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Burn crossings through road embankments (WP-D)."""

from qgis.core import (
    QgsProcessingParameterRasterLayer, QgsProcessingParameterFeatureSource,
    QgsProcessingParameterNumber, QgsProcessingParameterBoolean,
    QgsProcessingParameterRasterDestination, QgsProcessingParameterVectorDestination,
    QgsProcessing, QgsProcessingException,
)
import numpy as np

from .base import QehtAlgorithm
from ..core.raster import read_dem, write_raster
from ..core.conditioning.burn import burn_crossings
from ..core.network.crossings import select_crossings

DEM = "DEM"; CROSSINGS = "CROSSINGS"; ROAD = "ROAD"; ALL = "ALL"
HALF = "HALF"; SEARCH = "SEARCH"; OUTPUT = "OUTPUT"; LOG = "LOG"

LOG_FIELDS = [("crossing", "text"), ("cells", "int"), ("max_cut_m", "real"),
              ("total_cut_m3", "real"), ("z_up_m", "real"), ("z_down_m", "real"),
              ("length_m", "real"), ("note", "text")]


class BurnCrossingsAlgorithm(QehtAlgorithm):

    def name(self): return "burncrossings"
    def displayName(self): return "Burn crossings through embankments"
    def group(self): return "Road drainage"
    def groupId(self): return "roaddrainage"

    def shortHelpString(self):
        return (
            "A DSM or survey DEM shows a road embankment as a dam; the culvert under "
            "it is invisible. Filling then ponds the valley and sends the flow over "
            "the crest wherever it is lowest, which corrupts the catchment, the "
            "longest flow path and its slopes.\n\n"
            "At each crossing this tool cuts a short straight breach across the road: "
            "from the lowest cell near a point <i>half-length</i> upstream of the "
            "centreline to the lowest cell the same distance downstream, "
            "perpendicular to the road. Along it the ground is lowered to a straight "
            "grade between the two ends (never raised). Nothing else changes.\n\n"
            "<b>Crossings:</b> a crossing-candidate layer (uses accepted candidates, "
            "else recommended ones, and the road direction stored on each) or any "
            "point layer plus the road alignment (for the direction).\n\n"
            "Run it on the RAW DEM, then fill, flow direction and the rest on the "
            "burned DEM. Report the breaches (log layer) and state in the exchange "
            "package's 'conditioning' field that crossings were burned. Report "
            "elevations from the unburned raw DEM.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(DEM, "Raw DEM"))
        self.addParameter(QgsProcessingParameterFeatureSource(
            CROSSINGS, "Crossings (candidate layer or points)",
            [QgsProcessing.SourceType.TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterFeatureSource(
            ROAD, "Road alignment (needed if crossings carry no road_azimuth_deg)",
            [QgsProcessing.SourceType.TypeVectorLine], optional=True))
        self.addParameter(QgsProcessingParameterBoolean(
            ALL, "Burn every point (ignore candidate status)", defaultValue=False))
        self.addParameter(QgsProcessingParameterNumber(
            HALF, "Breach half-length either side of the centreline (m)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=30.0, minValue=1.0))
        self.addParameter(QgsProcessingParameterNumber(
            SEARCH, "Search radius for the breach ends (cells)",
            QgsProcessingParameterNumber.Type.Integer, defaultValue=2, minValue=0, maxValue=20))
        self.addParameter(QgsProcessingParameterRasterDestination(OUTPUT, "Burned DEM"))
        self.addParameter(QgsProcessingParameterVectorDestination(
            LOG, "Breach log", QgsProcessing.SourceType.TypeVectorLine))

    def processAlgorithm(self, parameters, context, feedback):
        dem_path = self.raster_path(parameters, DEM, context)
        dem, valid, info = read_dem(dem_path)
        pts = self.read_pour_points(parameters, CROSSINGS, context, info, feedback,
                                    extra_fields=("status", "recommended", "cand_id",
                                                  "outlet_uid", "road_azimuth_deg"))
        if not self.parameterAsBool(parameters, ALL, context) and \
                any("attr_status" in p for p in pts):
            idx, rule = select_crossings([{"status": p.get("attr_status"),
                                           "recommended": p.get("attr_recommended")}
                                          for p in pts])
            pts = [pts[k] for k in idx]
            feedback.pushInfo(f"Candidate layer: burning {len(pts)} crossing(s) ({rule}).")
        if not pts:
            raise QgsProcessingException("No crossings to burn.")
        alignment = None
        if any(p.get("attr_road_azimuth_deg") is None for p in pts):
            alignment = self.read_alignment(parameters, ROAD, context, info, feedback)
            if alignment is None:
                raise QgsProcessingException(
                    "The crossings carry no road direction (road_azimuth_deg); "
                    "supply the road alignment.")
        crossings = []
        for p in pts:
            cid = p.get("attr_outlet_uid") or p.get("attr_cand_id") or f"fid {p['fid']}"
            if p.get("attr_road_azimuth_deg") is not None:
                crossings.append({"x": p["x"], "y": p["y"], "id": cid,
                                  "azimuth_deg": float(p["attr_road_azimuth_deg"])})
            else:
                _, _, seg = alignment.locate([p["x"]], [p["y"]])
                tx, ty = alignment.tangent(int(seg[0]))
                crossings.append({"x": p["x"], "y": p["y"], "id": cid,
                                  "tx": float(tx), "ty": float(ty)})

        burned, log = burn_crossings(
            dem, valid, info.geotransform, crossings,
            half_length_m=self.parameterAsDouble(parameters, HALF, context),
            search_radius_cells=self.parameterAsInt(parameters, SEARCH, context))
        out = self.parameterAsOutputLayer(parameters, OUTPUT, context)
        write_raster(out, burned, info, valid=valid,
                     nodata=info.nodata if info.nodata is not None else -9999.0)

        rows = []
        for e in log:
            feedback.pushInfo(f"  {e['id']}: {e['cells']} cells, max cut {e['max_cut_m']:.2f} m, "
                              f"{e['total_cut_m3']:,.0f} m3 {e['note']}")
            path = e.get("path") if e.get("path") and len(e["path"]) >= 2 else None
            if path is None:            # not burned: keep it in the log at the crossing
                cxy = next((c["x"], c["y"]) for c in crossings if c["id"] == e["id"])
                path = [cxy, cxy]
            rows.append((path, {"crossing": e["id"], "cells": e["cells"],
                                         "max_cut_m": e["max_cut_m"],
                                         "total_cut_m3": e["total_cut_m3"],
                                         "z_up_m": e["z_up"], "z_down_m": e["z_down"],
                                         "length_m": e["length_m"], "note": e["note"]}))
        log_out = self.parameterAsOutputLayer(parameters, LOG, context)
        self.write_vector(log_out, "breaches", info.projection_wkt, "line", LOG_FIELDS, rows)
        feedback.pushInfo(f"Burned {sum(1 for e in log if e['cells'])} of {len(log)} crossings "
                          "(0 cells = no ground above the breach grade, i.e. no embankment "
                          "in the DEM there).")
        return {OUTPUT: out, LOG: log_out}
