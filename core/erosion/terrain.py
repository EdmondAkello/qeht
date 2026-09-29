# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Terrain erosion indices (WP-F, v0.14): slope, A_s, SPI, TWI, LS.

Numerical decisions (all recorded in the run metadata by the tools):

* Slope beta: Horn 3x3 on the RAW DEM (not the filled one - filling
  flattens depressions and would zero the slope there). tan(beta) is
  floored at BETA_MIN_TAN = 0.001 inside the indices so SPI, TWI and LS
  stay finite on flats; the reported slope raster is not floored.
* Specific catchment area A_s = (upslope cells + 1) * cell area / flow
  width, with flow width = cell size (sqrt(cell width * cell height)) -
  the usual D8 convention (Moore et al. 1991).
* SPI = A_s * tan(beta);  ln_spi = ln(SPI);  TWI = ln(A_s / tan(beta)).
* LS, default Moore & Burch (1986):
      LS = (A_s / 22.13)^m * (sin(beta) / 0.0896)^1.3,  m = 0.4
  or with m by slope class (0-1 % 0.2, 1-3 % 0.3, 3-5 % 0.4, 5-10 % 0.45,
  10-20 % 0.5, 20-30 % 0.55, > 30 % 0.6; after Renard et al. 1997 and
  McCool et al. 1987, as used for the A14 corridor, Akello & Omosa 2025).
  Note: that study used flow accumulation x cell size for A_s; QEHT adds
  the cell itself (acc + 1), so headwater cells are not zero.
  alternative Desmet & Govers (1996), per cell, with McCool (1987/1989)
  rill/interrill exponent m and slope factor S:
      L = ((A_in + D^2)^(m+1) - A_in^(m+1)) / (D^(m+2) * x^m * 22.13^m)
      m = F / (1 + F),  F = (sin b / 0.0896) / (3 sin(b)^0.8 + 0.56)
      S = 10.8 sin b + 0.03  (tan b < 0.09),  16.8 sin b - 0.50  otherwise
  where A_in is the area draining INTO the cell, D the cell size and
  x = |sin a| + |cos a| for aspect a.

References: Moore, Grayson & Ladson (1991) Hydrol. Process. 5, 3-30;
Moore & Burch (1986) SSSAJ 50, 1294-1298; Desmet & Govers (1996) JSWC 51,
427-433; McCool et al. (1987, 1989) Trans. ASAE.
"""

import numpy as np

BETA_MIN_TAN = 0.001
LS_METHODS = ("moore_burch", "desmet_govers")


def horn_gradients(elevation, valid, cell_width, cell_height):
    """dz/dx (east), dz/dy (south) by Horn's 3x3 weights; NaN outside `valid`."""
    z = np.asarray(elevation, dtype=np.float64)
    p = np.pad(np.where(valid, z, np.nan), 1, mode="edge")
    a = p[0:-2, 0:-2]; b = p[0:-2, 1:-1]; c = p[0:-2, 2:]
    d = p[1:-1, 0:-2];                    f = p[1:-1, 2:]
    g = p[2:, 0:-2];   h = p[2:, 1:-1];   i = p[2:, 2:]
    with np.errstate(invalid="ignore"):
        dzdx = ((c + 2.0 * f + i) - (a + 2.0 * d + g)) / (8.0 * cell_width)
        dzdy = ((g + 2.0 * h + i) - (a + 2.0 * b + c)) / (8.0 * cell_height)
    dzdx[~valid] = np.nan
    dzdy[~valid] = np.nan
    return dzdx, dzdy


M_SLOPE_CLASSES = ((1.0, 0.2), (3.0, 0.3), (5.0, 0.4), (10.0, 0.45), (20.0, 0.5),
                   (30.0, 0.55), (float("inf"), 0.6))   # (upper slope %, m)


def m_by_slope_class(tan_b):
    """m exponent by percent-slope class (A14 2025 Table 5)."""
    pct = 100.0 * np.asarray(tan_b, dtype=np.float64)
    uppers = np.array([u for u, _ in M_SLOPE_CLASSES[:-1]])
    ms = np.array([m for _, m in M_SLOPE_CLASSES])
    return ms[np.searchsorted(uppers, pct, side="left")]


def mccool_m(sin_b):
    """Rill/interrill slope-length exponent m (McCool et al. 1989)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        f = (sin_b / 0.0896) / (3.0 * np.power(sin_b, 0.8) + 0.56)
    return f / (1.0 + f)


def mccool_s(sin_b, tan_b):
    """RUSLE slope-steepness factor S (McCool et al. 1987)."""
    return np.where(tan_b < 0.09, 10.8 * sin_b + 0.03, 16.8 * sin_b - 0.50)


def erosion_indices(raw_dem, valid, accumulation, cell_width, cell_height,
                    ls_method="moore_burch", beta_min_tan=BETA_MIN_TAN, m_exponent=0.4):
    """Terrain erosion indices on the DEM grid.

    raw_dem : unfilled DEM (slope is taken from it)
    m_exponent : Moore & Burch slope-length exponent - a number (default
        0.4) or "slope_classes" (m by percent slope, A14 2025)
    accumulation : D8 upslope cell count (QEHT convention: excludes the cell)
    Returns dict of float64 grids (NaN outside `valid`):
      slope_deg, tan_beta (unfloored), spec_catch_area (m), spi, ln_spi,
      twi, ls, plus 'meta' with the numerical decisions.
    """
    if ls_method not in LS_METHODS:
        raise ValueError(f"ls_method must be one of {LS_METHODS}")
    cw, ch = float(abs(cell_width)), float(abs(cell_height))
    size = float(np.sqrt(cw * ch))
    dzdx, dzdy = horn_gradients(raw_dem, valid, cw, ch)
    tan_raw = np.hypot(dzdx, dzdy)
    tan_b = np.maximum(np.nan_to_num(tan_raw, nan=0.0), beta_min_tan)
    beta = np.arctan(tan_b)
    sin_b = np.sin(beta)
    acc = np.where(valid, np.asarray(accumulation, dtype=np.float64), np.nan)
    a_s = (acc + 1.0) * cw * ch / size
    with np.errstate(invalid="ignore", divide="ignore"):
        spi = a_s * tan_b
        ln_spi = np.log(spi)
        twi = np.log(a_s / tan_b)
        if ls_method == "moore_burch":
            m = m_by_slope_class(tan_b) if m_exponent == "slope_classes" else float(m_exponent)
            ls = np.power(a_s / 22.13, m) * np.power(sin_b / 0.0896, 1.3)
        else:
            m = mccool_m(sin_b)
            s = mccool_s(sin_b, tan_b)
            a_in = acc * cw * ch
            asp = np.arctan2(np.nan_to_num(dzdx), np.nan_to_num(dzdy))
            x = np.abs(np.sin(asp)) + np.abs(np.cos(asp))
            x = np.where(x > 0, x, 1.0)
            l_fac = ((np.power(a_in + size * size, m + 1.0) - np.power(a_in, m + 1.0))
                     / (np.power(size, m + 2.0) * np.power(x, m) * np.power(22.13, m)))
            ls = l_fac * s
    out = {"slope_deg": np.degrees(np.arctan(tan_raw)), "tan_beta": tan_raw,
           "spec_catch_area": a_s, "spi": spi, "ln_spi": ln_spi, "twi": twi, "ls": ls}
    for k, v in out.items():
        v[~valid] = np.nan
    out["meta"] = {"slope_method": "Horn 3x3 on raw DEM", "beta_min_tan": float(beta_min_tan),
                   "flow_width_m": size, "ls_method": ls_method,
                   "ls_m": (m_exponent if ls_method == "moore_burch" else "McCool (1989), per cell"),
                   "spec_catch_area": "(upslope cells + 1) * cell area / cell size"}
    return out


def channel_mask(accumulation, valid, threshold_cells):
    """Cells with at least `threshold_cells` upslope cells (concentrated flow)."""
    acc = np.asarray(accumulation, dtype=np.float64)
    return valid & (np.nan_to_num(acc, nan=-1.0) >= float(threshold_cells))
