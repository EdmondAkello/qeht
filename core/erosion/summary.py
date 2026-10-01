# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Catchment and crossing erosion summaries + crossing MCDMA (WP-F Level 3).

Catchment block (per catchment mask):
  mean / p90 RUSLE A (t/ha/yr), mean / p90 LS, channel ln(SPI) p90,
  area-weighted K, C, P (MUSLE inputs), shares of the combined severity
  score, gross soil loss, sediment delivery and sediment volume.
Crossing block (at the snapped outlet cell):
  3x3 max ln(SPI), approach-channel ln(SPI) p50/p90, local slope/LS/TWI,
  SPI class, worst combined score in the 3x3 window, and the multi-criteria
  hydrodynamic impact score of Akello & Omosa (2025):
      composite = 0.4 * SPI score + 0.3 * RUSLE score + 0.3 * sediment score
  with scores on the common 1-5 severity scale (classes.py) and levels
  <= 2.0 Low, <= 3.0 Moderate, <= 4.0 High, > 4.0 Severe.
Sediment delivery (area-based, FAO; as used on the A14 corridor):
      SDR = 0.565 * A_km2^-0.125 (capped at 1)
      sediment yield (t/yr) = gross soil loss (t/yr) * SDR
      volume (m3/yr) = sediment yield * 1000 / bulk density (kg/m3)
"""

import numpy as np

from .classes import get_scheme, combine_scores, COMBINED

MCDMA_WEIGHTS = {"spi": 0.4, "rusle": 0.3, "sediment": 0.3}
IMPACT_LEVELS = ((2.0, "Low", "Routine inspection and debris removal"),
                 (3.0, "Moderate", "Regular maintenance; sediment traps / check dams"),
                 (4.0, "High", "Energy dissipators, sediment traps, upsized structures, frequent "
                               "monitoring"),
                 (float("inf"), "Severe", "Energy dissipators, upsized structures, robust sediment "
                                          "control, upstream erosion control, frequent monitoring"))
DEFAULT_BULK_KGM3 = 1400.0
APPROACH_CELLS = 10


def sdr_area(area_km2, a=0.565, b=-0.125):
    """Area-based sediment delivery ratio (FAO), capped at 1."""
    if not area_km2 or area_km2 <= 0:
        return float("nan")
    return float(min(1.0, a * area_km2 ** b))


def impact_level(score):
    for upper, name, action in IMPACT_LEVELS:
        if score <= upper + 1e-9:
            return name, action
    return IMPACT_LEVELS[-1][1], IMPACT_LEVELS[-1][2]


def composite_score(spi_score, rusle_score, sediment_score, weights=None):
    w = dict(MCDMA_WEIGHTS if weights is None else weights)
    tot = w["spi"] + w["rusle"] + w["sediment"]
    return (w["spi"] * spi_score + w["rusle"] * rusle_score + w["sediment"] * sediment_score) / tot


class ErosionInputs(object):
    """Grids (all on the DEM grid) and settings for the erosion blocks.

    ln_spi, ls, twi, tan_beta : from terrain.erosion_indices
    soil_loss : RUSLE A grid (t/ha/yr) or None (LS-only fallback)
    k, c, p : factor grids or None
    channel : bool channel mask
    spi_scheme, rusle_scheme : ClassScheme objects (rusle_scheme is the
        LS-percentile scheme in LS-only mode)
    """

    def __init__(self, ln_spi, ls, twi, tan_beta, channel, spi_scheme, rusle_scheme,
                 soil_loss=None, k=None, c=None, p=None, volume_scheme=None,
                 weights=None, bulk_kgm3=DEFAULT_BULK_KGM3, matrix=None):
        self.ln_spi, self.ls, self.twi, self.tan_beta = ln_spi, ls, twi, tan_beta
        self.channel = channel
        self.soil_loss, self.k, self.c, self.p = soil_loss, k, c, p
        self.spi_scheme, self.rusle_scheme = spi_scheme, rusle_scheme
        self.volume_scheme = volume_scheme or get_scheme("volume_a14")
        self.weights = dict(MCDMA_WEIGHTS if weights is None else weights)
        self.bulk_kgm3 = float(bulk_kgm3)
        valid = np.isfinite(ln_spi)
        self.spi_codes = spi_scheme.classify(ln_spi, valid)
        base = soil_loss if soil_loss is not None else ls
        self.rusle_codes = rusle_scheme.classify(base, np.isfinite(base))
        self.combined = combine_scores(spi_scheme.score_of(self.spi_codes),
                                       rusle_scheme.score_of(self.rusle_codes), matrix)

    @property
    def mode(self):
        return "RUSLE" if self.soil_loss is not None else "LS-only"


def _stat(x, mask, fn):
    v = x[mask]
    v = v[np.isfinite(v)]
    if v.size == 0:
        return None
    return float(fn(v))


def catchment_erosion(mask, inp, cell_area_m2, bulk_kgm3=None):
    """Erosion block for one catchment (dict of ero_* fields)."""
    n = int(mask.sum())
    area_km2 = n * cell_area_m2 / 1e6
    out = {"ero_mode": inp.mode}
    out["ero_ls_mean"] = _stat(inp.ls, mask, np.mean)
    out["ero_ls_p90"] = _stat(inp.ls, mask, lambda v: np.percentile(v, 90))
    chm = mask & inp.channel
    out["ero_lnspi_ch_p90"] = _stat(inp.ln_spi, chm, lambda v: np.percentile(v, 90))
    for key, grid in (("ero_k_mean", inp.k), ("ero_c_mean", inp.c), ("ero_p_mean", inp.p)):
        out[key] = _stat(grid, mask, np.mean) if grid is not None else None
    sc = inp.combined[mask]
    sc = sc[sc > 0]
    for s in range(1, 6):
        out[f"ero_pct_s{s}"] = 100.0 * float((sc == s).sum()) / sc.size if sc.size else None
    if inp.soil_loss is not None:
        a_mean = _stat(inp.soil_loss, mask, np.mean)
        out["ero_a_mean_tha"] = a_mean
        out["ero_a_p90_tha"] = _stat(inp.soil_loss, mask, lambda v: np.percentile(v, 90))
        code = int(inp.rusle_scheme.classify(np.array([a_mean if a_mean is not None else np.nan]),
                                             np.array([a_mean is not None]))[0])
        out["ero_rusle_class"] = inp.rusle_scheme.names[code - 1] if code else None
        out["ero_rusle_score"] = inp.rusle_scheme.scores[code - 1] if code else None
        bd = float(bulk_kgm3) if bulk_kgm3 else inp.bulk_kgm3
        gross = (a_mean or 0.0) * area_km2 * 100.0          # t/ha/yr x ha
        sdr = sdr_area(area_km2)
        sy = gross * sdr if np.isfinite(sdr) else None
        vol = sy * 1000.0 / bd if sy is not None else None
        out.update(ero_gross_t_yr=gross, ero_sdr=sdr, ero_sy_t_yr=sy, ero_sy_m3_yr=vol,
                   ero_bulk_kgm3=bd)
        if vol is not None:
            vc = int(inp.volume_scheme.classify(np.array([vol]), np.array([True]))[0])
            out["ero_sy_class"] = inp.volume_scheme.names[vc - 1]
            out["ero_sy_score"] = inp.volume_scheme.scores[vc - 1]
    else:
        for k in ("ero_a_mean_tha", "ero_a_p90_tha", "ero_rusle_class", "ero_rusle_score",
                  "ero_gross_t_yr", "ero_sdr", "ero_sy_t_yr", "ero_sy_m3_yr", "ero_bulk_kgm3",
                  "ero_sy_class", "ero_sy_score"):
            out[k] = None
    return out


def _window(grid, r, c, half=1):
    rows, cols = grid.shape
    return grid[max(0, r - half):min(rows, r + half + 1), max(0, c - half):min(cols, c + half + 1)]


def approach_cells(direction, mask, channel, r, c, n=APPROACH_CELLS):
    """Channel cells upstream of (r, c) within n D8 steps, inside `mask`."""
    from ..grid import DROW, DCOL
    rows, cols = direction.shape
    seen = {(r, c)}
    front = [(r, c)]
    out = []
    for _ in range(int(n)):
        nxt = []
        for (a, b) in front:
            for k in range(8):
                ua, ub = a - int(DROW[k]), b - int(DCOL[k])      # neighbour that would flow in
                if 0 <= ua < rows and 0 <= ub < cols and (ua, ub) not in seen \
                        and mask[ua, ub] and channel[ua, ub] and int(direction[ua, ub]) == k:
                    seen.add((ua, ub)); nxt.append((ua, ub)); out.append((ua, ub))
        front = nxt
        if not front:
            break
    return out


def crossing_erosion(r, c, mask, direction, inp, catchment_block=None):
    """Erosion block for one crossing at outlet cell (r, c)."""
    w = _window(inp.ln_spi, r, c)
    out = {"ero_lnspi_max3x3": float(np.nanmax(w)) if np.isfinite(w).any() else None}
    app = approach_cells(direction, mask, inp.channel, r, c)
    vals = np.array([inp.ln_spi[a, b] for a, b in app], dtype=np.float64)
    vals = vals[np.isfinite(vals)]
    out["ero_lnspi_app_p50"] = float(np.percentile(vals, 50)) if vals.size else None
    out["ero_lnspi_app_p90"] = float(np.percentile(vals, 90)) if vals.size else None
    tw = _window(inp.tan_beta, r, c)
    out["ero_slope_pct"] = float(100.0 * np.nanmean(tw)) if np.isfinite(tw).any() else None
    for key, g in (("ero_ls_local", inp.ls), ("ero_twi_local", inp.twi)):
        ww = _window(g, r, c)
        out[key] = float(np.nanmean(ww)) if np.isfinite(ww).any() else None
    spi_val = out["ero_lnspi_max3x3"]
    if spi_val is not None:
        code = int(inp.spi_scheme.classify(np.array([spi_val]), np.array([True]))[0])
        out["ero_spi_class"] = inp.spi_scheme.names[code - 1]
        out["ero_spi_score"] = inp.spi_scheme.scores[code - 1]
    else:
        out["ero_spi_class"] = out["ero_spi_score"] = None
    cw = _window(inp.combined, r, c)
    out["ero_worst_score_3x3"] = int(cw.max()) if cw.size and cw.max() > 0 else None
    cb = catchment_block or {}
    rs, ss, ps = cb.get("ero_rusle_score"), cb.get("ero_sy_score"), out["ero_spi_score"]
    if None not in (rs, ss, ps):
        comp = composite_score(ps, rs, ss, inp.weights)
        lvl, act = impact_level(comp)
        out.update(ero_composite=comp, ero_impact=lvl, ero_mitigation=act)
    else:
        out.update(ero_composite=None, ero_impact=None,
                   ero_mitigation=None if inp.soil_loss is not None else
                   "LS-only mode: no soil loss, so no sediment score or composite")
    return out
