# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Erosion indices (SPI, TWI, LS), RUSLE soil loss and severity classes (WP-F, v0.14)."""

import json
import os

import numpy as np
from qgis.core import (QgsProcessingParameterRasterLayer, QgsProcessingParameterNumber,
                       QgsProcessingParameterEnum, QgsProcessingParameterString,
                       QgsProcessingParameterBoolean, QgsProcessingParameterFolderDestination,
                       QgsProcessingOutputString, QgsProcessingException)

from .base import QehtAlgorithm
from ..core.raster import read_dem, write_raster, warp_to_grid
from ..core.erosion.terrain import erosion_indices, channel_mask
from ..core.erosion.rusle import Factor, rusle, c_from_worldcover, WORLDCOVER_C, WORLDCOVER_P
from ..core.erosion.classes import (get_scheme, combine_scores, class_extents,
                                    write_class_raster, COMBINED, DEFAULT_MATRIX)
from ..core.erosion.summary import MCDMA_WEIGHTS, DEFAULT_BULK_KGM3
from ..core.erosion import io as eio

LS_OPTIONS = [("moore_burch", 0.4, "Moore & Burch (1986), m = 0.4 (default)"),
              ("moore_burch", "slope_classes", "Moore & Burch, m by slope class (A14 2025)"),
              ("desmet_govers", None, "Desmet & Govers (1996) with McCool m and S")]
SPI_OPTIONS = [("spi_percentile", "Percentiles 50/75/90/97 of this extent (relative, D7)"),
               ("spi_a14", "Fixed ln(SPI) 0 / 5 / 10: Low, Moderate, High, Severe (A14 2025)")]
RUSLE_OPTIONS = [("rusle_d7", "t/ha/yr 5 / 12 / 25 / 50: Very low ... Very high (D7)"),
                 ("rusle_a14", "t/ha/yr 5 / 10 / 20 / 40: Slight ... Severe (A14 2025)")]


def _floats(text, n=None, what="values"):
    if not text or not text.strip():
        return None
    try:
        vals = [float(v) for v in text.replace(";", ",").split(",") if v.strip()]
    except ValueError:
        raise QgsProcessingException(f"Could not read the {what} '{text}' - use numbers "
                                     "separated by commas.")
    if n is not None and len(vals) != n:
        raise QgsProcessingException(f"{what}: expected {n} numbers, got {len(vals)}.")
    return vals


class ErosionIndicesAlgorithm(QehtAlgorithm):

    def name(self): return "erosionindices"
    def displayName(self): return "Erosion indices and RUSLE soil loss"
    def group(self): return "Soils and erosion"
    def groupId(self): return "soils"

    def shortHelpString(self):
        return (
            "Terrain erosion indices, optional RUSLE soil loss and erosion severity "
            "classes, written to one output folder that the corridor sampler and "
            "the HEAS exchange package read.\n\n"
            "<b>Indices</b> (raw DEM for slope; flow accumulation for area): slope, "
            "specific catchment area A_s = (upslope cells + 1) x cell area / cell "
            "size, SPI = A_s tan(beta), ln(SPI), TWI = ln(A_s / tan(beta)) and the "
            "LS factor (Moore & Burch 1986 with m = 0.4, or with m by slope class, "
            "or Desmet & Govers 1996). tan(beta) is floored at 0.001 inside the "
            "indices only.\n\n"
            "<b>RUSLE</b> A = R K LS C P (t/ha/yr). Each factor is a raster, a single "
            "value or a dataset, in that order of precedence:\n"
            "• R: raster (e.g. GloREDa 2023 annual erosivity, recommended) or value, "
            "MJ mm/(ha h yr).\n"
            "• K: raster, value, or the SOTWIS soil map (Williams/EPIC K, 0-20 cm by "
            "default), SI units.\n"
            "• C: raster, value, or ESA WorldCover 2021 through an editable class "
            "lookup (a land-cover PROXY - flagged as such in every output).\n"
            "• P: raster, value (default 1.0 = no support practice), or the WorldCover "
            "P lookup of the A14 study.\n"
            "Rasters on other grids are resampled in-process onto the DEM grid "
            "(bilinear; nearest for WorldCover). If R, K or C is missing, no soil "
            "loss is computed: QEHT never invents a factor and falls back to LS-only "
            "classes labelled 'terrain potential'.\n\n"
            "<b>Classes</b> (uint8 rasters with colours, names, attribute table and a "
            "QGIS style): SPI by percentiles of this extent or fixed ln(SPI) 0/5/10; "
            "RUSLE by 5/12/25/50 or 5/10/20/40 t/ha/yr (or your own breaks); the "
            "combined class is the higher of the two severity scores unless you give "
            "a 5 x 5 matrix. class_extents.csv gives area and share per class, split "
            "into channel and hillslope cells. erosion_run.json records every "
            "factor's source, the schemes and the settings."
        )

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer("RAW_DEM", "Raw (unfilled) DEM"))
        self.addParameter(QgsProcessingParameterRasterLayer(
            "FAC", "Flow accumulation (cells, from the conditioned DEM)"))
        self.addParameter(QgsProcessingParameterEnum(
            "LS_METHOD", "LS factor", options=[o[2] for o in LS_OPTIONS], defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber(
            "CHANNEL_CELLS", "Channel threshold (upslope cells) for the channel/hillslope split",
            QgsProcessingParameterNumber.Type.Double, defaultValue=1000.0, minValue=1.0))
        for f, label in (("R", "R erosivity"), ("K", "K erodibility (SI)"), ("C", "C cover")):
            self.addParameter(QgsProcessingParameterRasterLayer(
                f + "_RASTER", f"{label} raster", optional=True))
            self.addParameter(QgsProcessingParameterNumber(
                f + "_VALUE", f"{label} single value (used when no raster)",
                QgsProcessingParameterNumber.Type.Double, optional=True, minValue=0.0))
        self.add_soil_parameters("SOIL", optional=True)
        self.addParameter(QgsProcessingParameterRasterLayer(
            "WORLDCOVER", "ESA WorldCover 2021 classes (for C, and P if chosen)", optional=True))
        self.addParameter(QgsProcessingParameterRasterLayer("P_RASTER", "P support raster",
                                                            optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            "P_VALUE", "P single value", QgsProcessingParameterNumber.Type.Double,
            defaultValue=1.0, minValue=0.0, maxValue=1.0))
        self.addParameter(QgsProcessingParameterBoolean(
            "P_FROM_WORLDCOVER", "P from WorldCover (A14 2025 lookup) instead of the value",
            defaultValue=False))
        self.addParameter(QgsProcessingParameterEnum(
            "SPI_SCHEME", "SPI classes", options=[o[1] for o in SPI_OPTIONS], defaultValue=0))
        self.addParameter(QgsProcessingParameterString(
            "SPI_BREAKS", "Custom ln(SPI) breaks (comma-separated; overrides the scheme)",
            optional=True))
        self.addParameter(QgsProcessingParameterEnum(
            "RUSLE_SCHEME", "RUSLE classes", options=[o[1] for o in RUSLE_OPTIONS], defaultValue=0))
        self.addParameter(QgsProcessingParameterString(
            "RUSLE_BREAKS", "Custom RUSLE breaks t/ha/yr (4 values; overrides the scheme)",
            optional=True))
        self.addParameter(QgsProcessingParameterString(
            "MATRIX", "Combination matrix, 25 values row by row (SPI score x RUSLE/LS score; "
            "blank = higher of the two)", optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            "BULK_DENSITY", "Bulk density for sediment volume when no soil data (kg/m3)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=DEFAULT_BULK_KGM3, minValue=100.0))
        self.addParameter(QgsProcessingParameterString(
            "WEIGHTS", "Crossing impact weights SPI, RUSLE, sediment", defaultValue="0.4,0.3,0.3"))
        self.addParameter(QgsProcessingParameterFolderDestination("OUTPUT", "Erosion output folder"))
        self.addOutput(QgsProcessingOutputString("SUMMARY", "Summary"))

    def _factor(self, parameters, context, info, key, name, feedback):
        lyr = self.parameterAsRasterLayer(parameters, key + "_RASTER", context)
        if lyr is not None:
            a, desc = warp_to_grid(self.raster_path(parameters, key + "_RASTER", context), info)
            return Factor(name, grid=a, source=desc)
        if parameters.get(key + "_VALUE") not in (None, ""):
            return Factor(name, value=self.parameterAsDouble(parameters, key + "_VALUE", context))
        return Factor(name)

    def processAlgorithm(self, parameters, context, feedback):
        folder = self.parameterAsString(parameters, "OUTPUT", context)
        os.makedirs(folder, exist_ok=True)
        raw, rv, info = read_dem(self.raster_path(parameters, "RAW_DEM", context))
        acc, av, ainfo = read_dem(self.raster_path(parameters, "FAC", context))
        if (ainfo.rows, ainfo.cols) != (info.rows, info.cols):
            raise QgsProcessingException("The flow accumulation and the DEM must be on the same grid.")
        valid = rv & av
        self.check_size(feedback, info.rows, info.cols, ["accumulation"])
        method, m_exp, _ = LS_OPTIONS[self.parameterAsEnum(parameters, "LS_METHOD", context)]
        feedback.pushInfo("Computing slope, A_s, SPI, TWI and LS ...")
        ind = erosion_indices(raw, valid, acc, info.cell_width, info.cell_height,
                              ls_method=method, m_exponent=m_exp if m_exp is not None else 0.4)
        chan = channel_mask(acc, valid, self.parameterAsDouble(parameters, "CHANNEL_CELLS", context))
        written = {}
        for key in ("slope_deg", "tan_beta", "spec_catch_area", "spi", "ln_spi", "twi", "ls"):
            written[key] = write_raster(eio.path_of(folder, key), ind[key], info, valid=valid)
        write_raster(eio.path_of(folder, "channel"), chan.astype(np.int32), info, dtype="int32",
                     nodata=-1, valid=valid)

        # -- factors ---------------------------------------------------------
        R = self._factor(parameters, context, info, "R", "R", feedback)
        K = self._factor(parameters, context, info, "K", "K", feedback)
        if not K.present:
            soil = self.load_soil(parameters, context, info, feedback)
            if soil is not None:
                grid, idx, units, sinfo = soil
                kg = np.full(grid.shape, np.nan)
                for i, uid in idx.items():
                    u = units.get(uid)
                    if u is not None and u.k_si is not None and np.isfinite(u.k_si):
                        kg[grid == i] = u.k_si
                K = Factor("K", grid=kg, source=f"{sinfo.get('soil_dataset')} Williams/EPIC K, "
                                               f"{sinfo.get('soil_depth_cm')} cm")
        C = self._factor(parameters, context, info, "C", "C", feedback)
        wc = None
        if self.parameterAsRasterLayer(parameters, "WORLDCOVER", context) is not None:
            wc, wdesc = warp_to_grid(self.raster_path(parameters, "WORLDCOVER", context), info,
                                     resampling="near")
        if not C.present and wc is not None:
            cg, cinfo = c_from_worldcover(np.nan_to_num(wc, nan=0).astype(np.int64), valid)
            if cinfo["unknown_classes"]:
                feedback.pushWarning(f"WorldCover classes without a C value: {cinfo['unknown_classes']}")
            C = Factor("C", grid=cg, source="ESA WorldCover 2021 class lookup (" + wdesc + ")",
                       proxy=True, note="land-cover class proxy, not a measured C")
        if self.parameterAsRasterLayer(parameters, "P_RASTER", context) is not None:
            pa, pdesc = warp_to_grid(self.raster_path(parameters, "P_RASTER", context), info)
            P = Factor("P", grid=pa, source=pdesc)
        elif self.parameterAsBool(parameters, "P_FROM_WORLDCOVER", context):
            if wc is None:
                raise QgsProcessingException("P from WorldCover needs the WorldCover raster.")
            pg, _ = c_from_worldcover(np.nan_to_num(wc, nan=0).astype(np.int64), valid,
                                      lookup=WORLDCOVER_P)
            P = Factor("P", grid=pg, source="ESA WorldCover 2021, A14 2025 P lookup", proxy=True,
                       note="land-cover class proxy")
        else:
            P = Factor("P", value=self.parameterAsDouble(parameters, "P_VALUE", context),
                       source="single value")
        a, finfo = rusle(ind["ls"], valid, R, K, C, P)
        for key, f in (("r", R), ("k", K), ("c", C), ("p", P)):
            if f.present:
                write_raster(eio.path_of(folder, key), f.as_grid(raw.shape), info, valid=valid)
        if a is not None:
            write_raster(eio.path_of(folder, "soil_loss"), a, info, valid=valid)
            fin = a[valid & np.isfinite(a)]
            feedback.pushInfo(f"RUSLE soil loss: mean {fin.mean():.2f}, p90 {np.percentile(fin, 90):.2f}, "
                              f"max {fin.max():.1f} t/ha/yr")
        else:
            feedback.pushWarning(finfo["note"])
        for d in finfo["factors"]:
            feedback.pushInfo(f"  {d['factor']}: {d['kind']} - {d['source']} ({d['basis']})")

        # -- classes -----------------------------------------------------------
        spi_key = SPI_OPTIONS[self.parameterAsEnum(parameters, "SPI_SCHEME", context)][0]
        spi_b = _floats(self.parameterAsString(parameters, "SPI_BREAKS", context), what="SPI breaks")
        if spi_b is not None and len(spi_b) not in (3, 4):
            raise QgsProcessingException("SPI breaks: give 3 (4 classes) or 4 (5 classes) values.")
        if spi_b is not None:
            base = "spi_a14" if len(spi_b) == 3 else "spi_percentile"
            spi_s = get_scheme(base, breaks=spi_b)
        else:
            spi_s = get_scheme(spi_key, ind["ln_spi"], valid)
        spi_codes = spi_s.classify(ind["ln_spi"], valid)
        if a is not None:
            ru_key = RUSLE_OPTIONS[self.parameterAsEnum(parameters, "RUSLE_SCHEME", context)][0]
            rb = _floats(self.parameterAsString(parameters, "RUSLE_BREAKS", context), 4, "RUSLE breaks")
            ru_s = get_scheme(ru_key, breaks=rb)
            ru_codes = ru_s.classify(a, valid)
            ru_file = "rusle_class"
        else:
            ru_s = get_scheme("ls_percentile", ind["ls"], valid)
            ru_codes = ru_s.classify(ind["ls"], valid)
            ru_file = "ls_class"
        mvals = _floats(self.parameterAsString(parameters, "MATRIX", context), 25, "matrix")
        matrix = [mvals[i * 5:(i + 1) * 5] for i in range(5)] if mvals else None
        comb = combine_scores(spi_s.score_of(spi_codes), ru_s.score_of(ru_codes), matrix)
        gt, wkt = info.geotransform, info.projection_wkt
        write_class_raster(eio.path_of(folder, "spi_class"), spi_codes, spi_s, gt, wkt, "SPI class")
        write_class_raster(eio.path_of(folder, ru_file), ru_codes, ru_s, gt, wkt,
                           "RUSLE class" if a is not None else "LS class (terrain potential)")
        write_class_raster(eio.path_of(folder, "combined_class"), comb, COMBINED, gt, wkt,
                           "Combined erosion severity")
        with open(os.path.join(folder, eio.EXTENTS_CSV), "w", encoding="utf-8") as f:
            f.write("raster,class,name,score,cells,area_km2,pct,channel_km2,hillslope_km2\n")
            for label, codes, sch in (("spi_class", spi_codes, spi_s), (ru_file, ru_codes, ru_s),
                                      ("combined_class", comb, COMBINED)):
                for r in class_extents(codes, sch, info.cell_area, channel=chan):
                    f.write(f"{label},{r['class']},{r['name']},{r['score']},{r['cells']},"
                            f"{r['area_km2']:.6f},{r['pct']:.4f},{r['channel_km2']:.6f},"
                            f"{r['hillslope_km2']:.6f}\n")
        w = _floats(self.parameterAsString(parameters, "WEIGHTS", context), 3, "weights") or \
            [MCDMA_WEIGHTS["spi"], MCDMA_WEIGHTS["rusle"], MCDMA_WEIGHTS["sediment"]]
        record = {
            "qeht_erosion_version": 1, "mode": finfo["mode"], "factors": finfo["factors"],
            "indices": ind["meta"], "channel_threshold_cells":
                self.parameterAsDouble(parameters, "CHANNEL_CELLS", context),
            "schemes": {"spi": spi_s.describe(), "rusle": ru_s.describe() if a is not None else None,
                        "ls": ru_s.describe() if a is None else
                        get_scheme("ls_percentile", ind["ls"], valid).describe(),
                        "volume": get_scheme("volume_a14").describe(),
                        "combined": COMBINED.describe()},
            "combination_matrix": matrix or DEFAULT_MATRIX,
            "mcdma_weights": {"spi": w[0], "rusle": w[1], "sediment": w[2]},
            "bulk_density_kgm3": self.parameterAsDouble(parameters, "BULK_DENSITY", context),
            "sdr_model": "0.565 * A_km2^-0.125, capped at 1 (FAO, area-based)",
            "c_lookup": {str(k): v for k, v in WORLDCOVER_C.items()},
            "dem": self.raster_path(parameters, "RAW_DEM", context),
            "grid": {"rows": info.rows, "cols": info.cols, "geotransform": list(gt)},
        }
        eio.write_run(folder, record)
        summary = (f"{finfo['mode']}; SPI classes: {spi_s.basis}; "
                   f"{'RUSLE' if a is not None else 'LS'} classes: {ru_s.basis}; folder {folder}")
        feedback.pushInfo(summary)
        return {"OUTPUT": folder, "SUMMARY": summary}
