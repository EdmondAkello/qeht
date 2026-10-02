# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Alignment ground profile (A1, v0.15)."""

import csv
import os

from qgis.core import (
    QgsProcessingParameterRasterLayer, QgsProcessingParameterFeatureSource,
    QgsProcessingParameterNumber, QgsProcessingParameterBoolean,
    QgsProcessingParameterString, QgsProcessingParameterVectorDestination,
    QgsProcessingParameterFileDestination, QgsProcessing, QgsProcessingException,
)
import numpy as np

from .base import QehtAlgorithm
from ..core.raster import read_dem
from ..core.grid import decode_d8
from ..core.network.profile import (profile_with_crossings, PROFILE_FIELDS,
                                    profile_summary, profile_chart)


class AlignmentProfileAlgorithm(QehtAlgorithm):

    def name(self): return "alignmentprofile"
    def displayName(self): return "Alignment ground profile"
    def group(self): return "Road drainage"
    def groupId(self): return "roaddrainage"

    def shortHelpString(self):
        return (
            "Samples the ground along a road centreline every step (default 10 m) "
            "of chainage: DEM ground (bilinear), filled-DEM level and the ponding "
            "depth the fill removed, contributing area, the stream crossings and "
            "their Strahler order, and the longitudinal ground slope.\n\n"
            "Use it for provisional culvert and relief levels before the geometric "
            "design exists, and later to check the DEM against the design "
            "existing-ground profile (DEM height bias). It is a terrain profile, not "
            "a design long section.\n\n"
            "<b>Stream crossings</b> are the exact places a D8 stream link crosses "
            "the centreline (the same intersection as 'Road crossing candidates'), "
            "at the stream threshold given. Each flags the station nearest to it. "
            "Contributing area at a station is the largest within one cell "
            "(search radius) so a stream one cell off the line is not missed.\n\n"
            "Only the raw DEM and the centreline are required; flow direction and "
            "accumulation add the stream fields, the filled DEM adds the ponding "
            "depth. Chainage runs in the digitised direction from the start value "
            "(tick Reverse to flip). A rough centreline is fine.\n\n"
            "'Build design hydrology package' writes the same profile into the "
            "package (layer alignment_profile) when a road alignment is given.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer("RAW_DEM", "Raw DEM"))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "ROAD", "Road centreline", [QgsProcessing.SourceType.TypeVectorLine]))
        self.addParameter(QgsProcessingParameterRasterLayer(
            "FILLED", "Filled DEM (optional; ponding depth)", optional=True))
        self.addParameter(QgsProcessingParameterRasterLayer(
            "FDR", "Flow direction (D8-coded; optional, for stream crossings)", optional=True))
        self.addParameter(QgsProcessingParameterRasterLayer(
            "FAC", "Flow accumulation (optional; contributing area)", optional=True))
        self.addParameter(QgsProcessingParameterRasterLayer(
            "ORDER", "Stream order raster (optional)", optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            "THRESHOLD", "Stream threshold (cells)", QgsProcessingParameterNumber.Type.Double,
            defaultValue=200.0, minValue=1.0))
        self.addParameter(QgsProcessingParameterNumber(
            "START", "Start chainage (m)", QgsProcessingParameterNumber.Type.Double,
            defaultValue=0.0))
        self.addParameter(QgsProcessingParameterBoolean(
            "REVERSE", "Reverse chainage direction", defaultValue=False))
        self.addParameter(QgsProcessingParameterNumber(
            "STEP", "Station spacing (m)", QgsProcessingParameterNumber.Type.Double,
            defaultValue=10.0, minValue=0.5))
        self.addParameter(QgsProcessingParameterString(
            "NAME", "Alignment name (written to every station)", optional=True))
        self.addParameter(QgsProcessingParameterVectorDestination(
            "PROFILE", "Alignment profile stations", QgsProcessing.SourceType.TypeVectorPoint))
        self.addParameter(QgsProcessingParameterFileDestination(
            "CSV", "Profile table (CSV)", fileFilter="CSV (*.csv)", optional=True,
            createByDefault=False))
        self.addParameter(QgsProcessingParameterFileDestination(
            "CHART", "Profile chart (PNG)", fileFilter="PNG (*.png)", optional=True,
            createByDefault=False))

    def processAlgorithm(self, parameters, context, feedback):
        raw, rv, info = read_dem(self.raster_path(parameters, "RAW_DEM", context))
        raw = np.where(rv, raw, np.nan)

        def _opt(name):
            if self.parameterAsRasterLayer(parameters, name, context) is None:
                return None, None
            a, v, i = read_dem(self.raster_path(parameters, name, context))
            if (i.rows, i.cols) != (info.rows, info.cols):
                raise QgsProcessingException(
                    f"{name}: the raster is not on the raw DEM grid ({i.rows}x{i.cols} vs "
                    f"{info.rows}x{info.cols}). Use rasters derived from this DEM.")
            return a, v
        filled, fv = _opt("FILLED")
        if filled is not None:
            filled = np.where(fv, filled, np.nan)
        d8, dv = _opt("FDR")
        acc, av = _opt("FAC")
        order, ov = _opt("ORDER")
        if order is not None:
            order = np.where(ov, order, np.nan)
        direction = decode_d8(d8.astype(np.int32)) if d8 is not None else None
        if direction is None and acc is not None:
            feedback.pushInfo("No flow direction given: contributing area only, no stream crossings.")

        name = (self.parameterAsString(parameters, "NAME", context) or "").strip()
        if not name:
            src = self.parameterAsSource(parameters, "ROAD", context)
            name = src.sourceName() if src is not None else ""
        alignment = self.read_alignment(
            parameters, "ROAD", context, info, feedback,
            start_chainage=self.parameterAsDouble(parameters, "START", context),
            reverse=self.parameterAsBool(parameters, "REVERSE", context))
        rows, crossings = profile_with_crossings(
            alignment, raw, info.geotransform, direction=direction,
            valid=(dv if dv is not None else av), accumulation=acc, filled=filled,
            stream_threshold_cells=self.parameterAsDouble(parameters, "THRESHOLD", context),
            stream_order=order, step=self.parameterAsDouble(parameters, "STEP", context),
            align_name=name)
        out = self.parameterAsOutputLayer(parameters, "PROFILE", context)
        self.write_vector(out, "alignment_profile", info.projection_wkt, "point",
                          PROFILE_FIELDS, [((r["x"], r["y"]), r) for r in rows])
        result = {"PROFILE": out}
        csv_path = self.parameterAsFileOutput(parameters, "CSV", context)
        if csv_path:
            with open(csv_path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow([k for k, _ in PROFILE_FIELDS])
                for r in rows:
                    w.writerow(["" if r[k] is None else r[k] for k, _ in PROFILE_FIELDS])
            result["CSV"] = csv_path
        chart = self.parameterAsFileOutput(parameters, "CHART", context)
        if chart:
            if profile_chart(chart, rows, f"Ground profile - {name}" if name else
                             "Alignment ground profile"):
                result["CHART"] = chart
            else:
                feedback.pushWarning("matplotlib is not available - chart skipped.")
        feedback.pushInfo(profile_summary(rows))
        for c in crossings:
            feedback.pushInfo(f"  stream at ch {c['chainage_m']:,.1f}  A={c['acc_km2']:.4f} km2"
                              + (f"  order {c['stream_order']}" if c.get("stream_order") else ""))
        return result
