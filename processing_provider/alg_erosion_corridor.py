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
            "otherwise.\n\n"
            "<b>Side-drain siltation (STI R3)</b>, with the flow direction: on each side, the "
            "corridor cells whose D8 path reaches the road strip (within the road half-width of "
            "the centreline) without entering a channel are the slope that feeds the side drain. "
            "Per station and side the 50th and 90th percentiles of their overland STI "
            "(sti_overland.tif) and the cell count; LHS / RHS are left / right looking "
            "up-chainage. Classes low / moderate / high / very high are RELATIVE to this "
            "corridor (its 50 / 75 / 90th percentiles): relative transport capacity, not "
            "severity. With the raw DEM, the longitudinal ground slope is added and a station "
            "is flagged (siltation_lhs / _rhs = 1) where the class is high or very high and the "
            "slope is below the drain slope limit (default 1 %): a low-gradient drain fed by a "
            "high transport-capacity slope. Screening only; not used in the erosion classes, "
            "RUSLE or the composite score."
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
        from qgis.core import QgsProcessingParameterRasterLayer as _RL
        self.addParameter(_RL("FDR", "Flow direction (D8-coded; optional, for the side-drain "
                              "siltation indicator)", optional=True))
        self.addParameter(_RL("RAW_DEM", "Raw DEM (optional; longitudinal slope for the "
                              "siltation flag)", optional=True))
        for key, label, default in (
                ("ROAD_HALF", "Side drains: road strip half-width (m)", 5.0),
                ("DRAIN_SLOPE", "Side drains: flag below this longitudinal slope (%)", 1.0)):
            self.addParameter(self._advanced(QgsProcessingParameterNumber(
                key, label, QgsProcessingParameterNumber.Type.Double, defaultValue=default,
                minValue=0.0)))
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
        sti_md = None
        if self.parameterAsRasterLayer(parameters, "FDR", context) is not None:
            sti_md = self._side_drains(parameters, context, feedback, folder, road, info,
                                       stations, reaches)
        elif os.path.exists(eio.path_of(folder, "sti")):
            feedback.pushInfo("Give the flow direction for the side-drain siltation indicator.")
        names = COMBINED.names
        def cname(s):
            return names[s - 1] if s else None
        keys = [k for k in stations[0] if k not in ("x", "y")] if stations else []
        from ..core.erosion.side_drain import STATION_FIELDS, REACH_FIELDS
        sd_keys = {f for f, _ in STATION_FIELDS}
        fields = [(k, "text" if k == "worst_side" else ("int" if k.endswith("score") else "real"))
                  for k in keys if k not in sd_keys] + [("centre_class", "text"),
                                                        ("worst_class", "text")] + STATION_FIELDS
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
                   ("ln_spi_max", "real"), ("soil_loss_max", "real")] + REACH_FIELDS
        re_out = self.parameterAsOutputLayer(parameters, "REACHES", context)
        self.write_vector(re_out, "erosion_reaches", info.projection_wkt, "line", rfields, rrows)
        result = {"STATIONS": st_out, "REACHES": re_out}
        chart = self.parameterAsFileOutput(parameters, "CHART", context)
        if chart:
            if profile_chart(chart, stations, f"Erosion along the alignment ({run.get('mode')})"):
                result["CHART"] = chart
            else:
                feedback.pushWarning("matplotlib is not available - chart skipped.")
        if sti_md is not None:
            for side, lbl in (("lhs", "left"), ("rhs", "right")):
                flen = sum(r.get(f"siltation_len_{side}_m") or 0.0 for r in reaches)
                feedback.pushInfo(f"Side drains, {lbl}: {flen:,.0f} m flagged for siltation "
                                  "(screening).")
        bad = [r for r in reaches if (r["worst_score"] or 0) >= 4]
        feedback.pushInfo(f"{len(stations)} stations, {len(reaches)} reaches; "
                          f"{len(bad)} reach(es) High or Very high, "
                          f"{sum(r['length_m'] for r in bad):,.0f} m in total.")
        for r in bad[:30]:
            feedback.pushInfo(f"  ch {r['ch_start']:,.0f} - {r['ch_end']:,.0f}: "
                              f"{cname(r['worst_score'])}, max ln SPI "
                              f"{r['ln_spi_max'] if r['ln_spi_max'] is None else round(r['ln_spi_max'], 2)}")
        return result

    def _side_drains(self, parameters, context, feedback, folder, road, info, stations, reaches):
        """STI R3: adds the side-drain fields to stations and reaches in place."""
        from ..core.grid import decode_d8
        from ..core.erosion.side_drain import corridor_sti, params_json
        sti_path = eio.path_of(folder, "sti")
        if not os.path.exists(sti_path):
            feedback.pushWarning("No sti_overland.tif in the erosion folder (STI switched off): "
                                 "side-drain siltation indicator skipped.")
            return None
        d8, dv, dinfo = read_dem(self.raster_path(parameters, "FDR", context))
        if (dinfo.rows, dinfo.cols) != (info.rows, info.cols):
            raise QgsProcessingException("The flow direction is not on the erosion folder's grid.")
        sti, sv, _ = read_dem(sti_path)
        chan = None
        if os.path.exists(eio.path_of(folder, "channel")):
            ch, cv, _ = read_dem(eio.path_of(folder, "channel"))
            chan = cv & (ch > 0)
        slope_at, src = None, "none (no raw DEM: siltation flags empty)"
        step = self.parameterAsDouble(parameters, "STEP", context)
        if self.parameterAsRasterLayer(parameters, "RAW_DEM", context) is not None:
            from ..core.network.profile import alignment_profile
            raw, rv, rinfo = read_dem(self.raster_path(parameters, "RAW_DEM", context))
            if (rinfo.rows, rinfo.cols) != (info.rows, info.cols):
                raise QgsProcessingException("The raw DEM is not on the erosion folder's grid.")
            prof = alignment_profile(road, np.where(rv, raw, np.nan), info.geotransform, step=step)
            pch = np.array([r["chainage_m"] for r in prof])
            psl = np.array([np.nan if r["slope_long_pct"] is None else r["slope_long_pct"]
                            for r in prof])

            def slope_at(ch):
                return psl[np.abs(pch[None, :] - np.asarray(ch)[:, None]).argmin(axis=1)]
            src = "alignment profile on the raw DEM (centred difference)"
        half = self.parameterAsDouble(parameters, "HALF_WIDTH", context)
        off = self.parameterAsDouble(parameters, "OFFSET_STEP", context)
        rh = self.parameterAsDouble(parameters, "ROAD_HALF", context)
        ds = self.parameterAsDouble(parameters, "DRAIN_SLOPE", context)
        st, rr, br = corridor_sti(road, np.where(sv, sti, np.nan), decode_d8(d8.astype(np.int32)),
                                  dv, info.geotransform, reaches, channel=chan, slope_at=slope_at,
                                  step=step, half_width=half, offset_step=off, road_half_m=rh,
                                  drain_slope_pct=ds)
        for a, b in zip(stations, st):
            a.update({k: v for k, v in b.items() if k not in ("chainage", "x", "y")})
        for a, b in zip(reaches, rr):
            a.update(b)
        feedback.pushInfo("Side-drain siltation indicator: class breaks (P50/P75/P90 of station "
                          f"STI p90) {br}; slope from {src}.")
        return params_json(step, half, off, rh, ds, br, src)
