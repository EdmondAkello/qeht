# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Clip DEM to the road's contributing area (F13, v0.25)."""

import numpy as np
from qgis.core import (QgsProcessing, QgsProcessingException, QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterFeatureSource, QgsProcessingParameterNumber,
                       QgsProcessingParameterEnum, QgsProcessingParameterRasterDestination,
                       QgsProcessingOutputString)

from .base import QehtAlgorithm
from ..core.raster import read_dem
from ..core.conditioning.autoclip import contributing_window, write_clip, memory_note

MODES = ["Bounding box (recommended)", "Mask (NoData outside the contributing area)"]


class AutoClipAlgorithm(QehtAlgorithm):

    def name(self): return "autoclip"
    def displayName(self): return "Clip DEM to the road's contributing area"
    def group(self): return "Road drainage"
    def groupId(self): return "roaddrainage"

    def shortHelpString(self):
        return (
            "Cuts a large DEM down to the area that drains to the road, so the rest of QEHT "
            "fits in memory, without cutting catchments by hand.\n\n"
            "<b>Method:</b> a coarse copy of the DEM (the minimum of each k x k block, so "
            "valleys are kept; k chosen so the copy has at most 4 Mcells) is filled and routed "
            "with D8. Coarse cells within the road buffer (200 m) are marked, every cell that "
            "drains to them is found, and that area is widened by the margin (1 km plus 2 "
            "coarse cells). The bounding box (or the mask) of it is read from the full-resolution "
            "DEM.\n\n"
            "<b>Check after the run:</b> the clipped DEM is tagged; the design hydrology package "
            "marks every catchment that touches the clip edge (clip_edge = 1). Enlarge the margin "
            "and run again if any does.\n\n"
            "Reports the memory estimate before and after.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer("DEM", "DEM (projected, metres)"))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "ROAD", "Road alignment", [QgsProcessing.SourceType.TypeVectorLine]))
        for key, label, default in (("BUFFER", "Road buffer (m)", 200.0),
                                    ("MARGIN", "Margin around the contributing area (m)", 1000.0),
                                    ("MAX_MCELLS", "Coarse grid at most (million cells)", 4.0)):
            prm = QgsProcessingParameterNumber(
                key, label, QgsProcessingParameterNumber.Type.Double, defaultValue=default,
                minValue=0.0 if key != "MAX_MCELLS" else 0.01)
            self.addParameter(self._advanced(prm) if key == "MAX_MCELLS" else prm)
        self.addParameter(QgsProcessingParameterEnum("MODE", "Clip to", options=MODES, defaultValue=0))
        self.addParameter(QgsProcessingParameterRasterDestination("OUTPUT", "Clipped DEM"))
        self.addOutput(QgsProcessingOutputString("SUMMARY", "Summary"))

    def processAlgorithm(self, parameters, context, feedback):
        path = self.raster_path(parameters, "DEM", context)
        z, valid, info = read_dem(path)
        road = self.read_alignment(parameters, "ROAD", context, info, feedback)
        if road is None:
            raise QgsProcessingException("Give the road alignment.")
        mode = ["bbox", "mask"][self.parameterAsEnum(parameters, "MODE", context)]
        buf = self.parameterAsDouble(parameters, "BUFFER", context)
        mar = self.parameterAsDouble(parameters, "MARGIN", context)
        try:
            w = contributing_window(np.where(valid, z, np.nan), valid, info.geotransform, road,
                                    road_buffer_m=buf, margin_m=mar,
                                    max_mcells=self.parameterAsDouble(parameters, "MAX_MCELLS", context),
                                    mode=mode)
        except ValueError as e:
            raise QgsProcessingException(str(e))
        del z
        if w["k"] == 1:
            feedback.pushWarning("The DEM is already under the coarse-grid limit, so the clip was "
                                 "found at full resolution (a full fill and D8 run). The clip pays "
                                 "off on DEMs larger than the limit.")
        mem = memory_note(w["full_cells"], w.get("clip_cells_valid", w["clip_cells"]))
        meta = {"k": w["k"], "road_buffer_m": buf, "margin_m": mar, "mode": mode,
                "full_shape": [info.rows, info.cols], "source": path, **mem}
        out = self.parameterAsOutputLayer(parameters, "OUTPUT", context)
        write_clip(path, out, info, w["window"], mask=w.get("mask"), meta=meta)
        r0, r1, c0, c1 = w["window"]
        summary = (f"Clip rows {r0}-{r1}, cols {c0}-{c1} of {info.rows} x {info.cols} "
                   f"({100.0 * w['clip_cells'] / w['full_cells']:.0f} % of the cells; coarse factor "
                   f"{w['k']}); peak memory about {mem['full_peak_gb']:.2f} GB -> "
                   f"{mem['clip_peak_gb']:.2f} GB")
        feedback.pushInfo(summary)
        return {"OUTPUT": out, "SUMMARY": summary}
