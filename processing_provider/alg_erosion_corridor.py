# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Sample erosion indices along a road alignment (WP-F Level 2, v0.14)."""

import os

import numpy as np
from qgis.core import (QgsProcessingParameterFile, QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterNumber, QgsProcessingParameterBoolean,
                       QgsProcessingParameterVectorDestination,
                       QgsProcessingParameterFileDestination, QgsProcessing,
                       QgsProcessingException)

from .base import QehtAlgorithm
from ..core.raster import read_dem
from ..core.erosion import io as eio
from ..core.erosion.corridor import sample_corridor, profile_chart, GRIDS
from ..core.erosion.classes import COMBINED


class ErosionCorridorAlgorithm(QehtAlgorithm):

    def name(self): return "erosioncorridor"
    def displayName(self): return "Sample erosion along alignment"
    def group(self): return "Soils and erosion"
    def groupId(self): return "soils"

    def shortHelpString(self):
        return (
            "Samples the erosion folder written by 'Erosion indices and RUSLE soil "
            "loss' along a road centreline.\n\n"
            "Stations every step (default 10 m) of chainage. At each station the "
            "grids are read on the centreline and on perpendicular offsets out to "
            "the half-width on the left and right. Stations carry max and mean "
            "ln(SPI), LS, RUSLE soil loss and TWI per side, the combined erosion "
            "class on the centreline and the worst class within the buffer.\n\n"
            "Consecutive stations with the same worst class form a reach. Reach "
            "boundaries sit half-way between stations, so reach lengths add up to "
            "the alignment length. Reaches show where concentrated flow and "
            "gully-prone ground meet the road, i.e. where washouts and silting are "
            "most likely.\n\n"
            "Chainage runs from the start value in the direction the line was "
            "digitised (tick Reverse to flip it). A rough centreline is fine; "
            "chainage is provisional until the geometric design exists. The chart "
            "needs matplotlib, which ships with QGIS on most installs; it is skipped "
            "otherwise."
        )

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile(
            "EROSION", "Erosion output folder",
            behavior=QgsProcessingParameterFile.Behavior.Folder))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "ROAD", "Road centreline", [QgsProcessing.SourceType.TypeVectorLine]))
        self.addParameter(QgsProcessingParameterNumber(
            "START", "Chainage at the start of the line (m)", QgsProcessingParameterNumber.Type.Double,
            defaultValue=0.0))
        self.addParameter(QgsProcessingParameterBoolean("REVERSE", "Reverse the chainage direction",
                                                        defaultValue=False))
        self.addParameter(QgsProcessingParameterNumber(
            "STEP", "Station spacing (m)", QgsProcessingParameterNumber.Type.Double,
            defaultValue=10.0, minValue=1.0))
        self.addParameter(QgsProcessingParameterNumber(
            "HALF_WIDTH", "Buffer half-width each side (m)", QgsProcessingParameterNumber.Type.Double,
            defaultValue=50.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber(
            "OFFSET_STEP", "Sampling interval across the road (m)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=10.0, minValue=1.0))
        self.addParameter(QgsProcessingParameterVectorDestination(
            "STATIONS", "Erosion stations", QgsProcessing.SourceType.TypeVectorPoint))
        self.addParameter(QgsProcessingParameterVectorDestination(
            "REACHES", "Erosion reaches", QgsProcessing.SourceType.TypeVectorLine))
        self.addParameter(QgsProcessingParameterFileDestination(
            "CHART", "Chainage profile chart", fileFilter="PNG (*.png)", optional=True,
            createByDefault=False))

    def processAlgorithm(self, parameters, context, feedback):
        folder = self.parameterAsFile(parameters, "EROSION", context)
        run = eio.read_run(folder)
        comb, cv, info = read_dem(eio.path_of(folder, "combined_class"))
        grids = {}
        for g in GRIDS:
            p = eio.path_of(folder, g)
            if os.path.exists(p):
                a, v, _ = read_dem(p)
                grids[g] = np.where(v, a, np.nan)
        road = self.read_alignment(parameters, "ROAD", context, info, feedback,
                                   start_chainage=self.parameterAsDouble(parameters, "START", context),
                                   reverse=self.parameterAsBool(parameters, "REVERSE", context))
        if road is None:
            raise QgsProcessingException("No road centreline.")
        stations, reaches = sample_corridor(
            road, grids, np.where(cv, comb, 0).astype(np.uint8), info.geotransform,
            step=self.parameterAsDouble(parameters, "STEP", context),
            half_width=self.parameterAsDouble(parameters, "HALF_WIDTH", context),
            offset_step=self.parameterAsDouble(parameters, "OFFSET_STEP", context))
        names = COMBINED.names
        def cname(s):
            return names[s - 1] if s else None
        keys = [k for k in stations[0] if k not in ("x", "y")] if stations else []
        fields = [(k, "text" if k == "worst_side" else ("int" if k.endswith("score") else "real"))
                  for k in keys] + [("centre_class", "text"), ("worst_class", "text")]
        rows = [((s["x"], s["y"]), dict(s, centre_class=cname(s["centre_score"]),
                                        worst_class=cname(s["worst_score"]))) for s in stations]
        st_out = self.parameterAsOutputLayer(parameters, "STATIONS", context)
        self.write_vector(st_out, "erosion_stations", info.projection_wkt, "point", fields, rows)
        rrows = []
        for r in reaches:
            chs = [r["ch_start"]] + [s["chainage"] for s in stations
                                     if r["ch_start"] < s["chainage"] < r["ch_end"]] + [r["ch_end"]]
            line = [road.point_at(c) for c in chs]
            if len(line) >= 2:
                rrows.append((line, dict(r, worst_class=cname(r["worst_score"]))))
        rfields = [("ch_start", "real"), ("ch_end", "real"), ("length_m", "real"),
                   ("worst_score", "int"), ("worst_class", "text"), ("n_stations", "int"),
                   ("ln_spi_max", "real"), ("soil_loss_max", "real")]
        re_out = self.parameterAsOutputLayer(parameters, "REACHES", context)
        self.write_vector(re_out, "erosion_reaches", info.projection_wkt, "line", rfields, rrows)
        result = {"STATIONS": st_out, "REACHES": re_out}
        chart = self.parameterAsFileOutput(parameters, "CHART", context)
        if chart:
            if profile_chart(chart, stations, f"Erosion along the alignment ({run.get('mode')})"):
                result["CHART"] = chart
            else:
                feedback.pushWarning("matplotlib is not available - chart skipped.")
        bad = [r for r in reaches if (r["worst_score"] or 0) >= 4]
        feedback.pushInfo(f"{len(stations)} stations, {len(reaches)} reaches; "
                          f"{len(bad)} reach(es) High or Very high, "
                          f"{sum(r['length_m'] for r in bad):,.0f} m in total.")
        for r in bad[:30]:
            feedback.pushInfo(f"  ch {r['ch_start']:,.0f} - {r['ch_end']:,.0f}: "
                              f"{cname(r['worst_score'])}, max ln SPI "
                              f"{r['ln_spi_max'] if r['ln_spi_max'] is None else round(r['ln_spi_max'], 2)}")
        return result
