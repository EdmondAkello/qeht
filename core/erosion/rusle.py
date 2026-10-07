# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""RUSLE A = R * K * LS * C * P with source-tagged factors (WP-F, v0.14).

Every factor is supplied as a raster already on the DEM grid, a single
value, or (C only) an ESA WorldCover 2021 class raster translated by an
editable lookup. Each factor carries provenance: source text and whether
it is measured/derived or a class PROXY. If R, K or C is missing, no soil
loss is computed and the tools fall back to LS-only "terrain potential" -
QEHT never invents a factor.

Units: R in MJ mm ha-1 h-1 yr-1, K in t ha h ha-1 MJ-1 mm-1 (SI), so A is in
t ha-1 yr-1. (US-customary K must be multiplied by 0.1317 first.)
"""

import numpy as np

# ESA WorldCover 2021 (v200) class -> C proxy: (C, low, high, basis).
# Classes 10-60: the scheme published for the A14 corridor (Akello & Omosa
# 2025, AJERI, Table 6), C = mean of the typical range, which cites Renard
# et al. (1997), Panagos et al. (2015), Wischmeier & Smith (1978) and Ganasri
# & Ramesh (2016). Classes 70-100 (absent from that study area) from
# Panagos et al. (2015) and Linard et al. (2014). A class-based C is a
# land-cover PROXY, not a measured cover-management factor, and is labelled
# so in every output. Edit the table (or pass a lookup) for other regions.
WORLDCOVER_C = {
    10: (0.0255, 0.001, 0.05, "Tree cover 0.001-0.05 (A14 2025; Renard 1997; Panagos 2015)"),
    20: (0.03, 0.01, 0.05, "Shrubland 0.01-0.05 (A14 2025; Ganasri & Ramesh 2016; Panagos 2015)"),
    30: (0.08, 0.01, 0.15, "Grassland 0.01-0.15 (A14 2025; Wischmeier & Smith 1978; Renard 1997)"),
    40: (0.25, 0.10, 0.40, "Cropland 0.1-0.4 (A14 2025; Wischmeier & Smith 1978; Renard 1997)"),
    50: (0.125, 0.05, 0.20, "Built-up 0.05-0.2 (A14 2025; Panagos 2015; Ganasri & Ramesh 2016)"),
    60: (0.50, 0.40, 0.60, "Bare / sparse vegetation 0.4-0.6 (A14 2025; Renard 1997; Panagos 2015)"),
    70: (0.0, 0.0, 0.0, "Snow and ice - no soil erosion"),
    80: (0.0, 0.0, 0.0, "Permanent water bodies - no soil erosion"),
    90: (0.001, 0.001, 0.003, "Herbaceous wetland 0.001 (Linard 2014)"),
    95: (0.001, 0.001, 0.003, "Mangroves (woody wetland) 0.001 (Linard 2014)"),
    100: (0.01, 0.01, 0.08, "Moss and lichen - natural grassland lower bound (Panagos 2015)"),
}
# Optional class -> P proxy (A14 2025, Table 7: natural cover assumed to act
# as a support practice in a pastoral area with no formal conservation
# works). Default P is 1.0 (no support practice) unless this is chosen.
WORLDCOVER_P = {10: 0.6, 20: 0.8, 30: 0.9, 40: 0.7, 50: 1.0, 60: 1.0,
                70: 1.0, 80: 1.0, 90: 1.0, 95: 1.0, 100: 1.0}
WORLDCOVER_NAMES = {10: "Tree cover", 20: "Shrubland", 30: "Grassland", 40: "Cropland",
                    50: "Built-up", 60: "Bare / sparse vegetation", 70: "Snow and ice",
                    80: "Permanent water bodies", 90: "Herbaceous wetland", 95: "Mangroves",
                    100: "Moss and lichen"}


class Factor(object):
    """One RUSLE factor on the DEM grid with its provenance."""

    def __init__(self, name, grid=None, value=None, source="", proxy=False, note="",
                 basis=None):
        self.name = name
        self.basis = basis
        self.grid = None if grid is None else np.asarray(grid, dtype=np.float64)
        self.value = None if value is None else float(value)
        self.source = source or ("single value" if value is not None else "")
        self.proxy = bool(proxy)
        self.note = note

    @property
    def present(self):
        return self.grid is not None or self.value is not None

    def as_grid(self, shape):
        if self.grid is not None:
            return self.grid
        return np.full(shape, self.value, dtype=np.float64)

    def describe(self):
        kind = "raster" if self.grid is not None else ("value" if self.value is not None else "missing")
        d = {"factor": self.name, "kind": kind, "source": self.source,
             "basis": self.basis or ("class proxy" if self.proxy else "user value" if kind == "value"
                       else "raster as supplied" if self.present else "missing")}
        if self.value is not None:
            d["value"] = self.value
        if self.note:
            d["note"] = self.note
        return d


WORLDCOVER_CODES = frozenset(WORLDCOVER_C)
# plausible ranges (SI units). Outside them the factor is refused or flagged.
FACTOR_RANGES = {"R": (0.0, 30000.0), "K": (0.0, 0.1), "C": (0.0, 1.0), "P": (0.0, 1.0)}


class FactorError(ValueError):
    """A RUSLE factor is outside its physical range (wrong units or wrong raster)."""


def check_factor(name, grid=None, value=None):
    """Range check of one RUSLE factor; returns a list of warnings, raises
    FactorError when the factor cannot be right.

    Catches the common mistakes: a land-cover CLASS raster (codes 10-100)
    given as C or P; K in US units (x 7.59 of SI) or in percent; R or K
    of the wrong order of magnitude.
    """
    lo, hi = FACTOR_RANGES[name]
    if grid is not None:
        g = np.asarray(grid, dtype=np.float64)
        v = g[np.isfinite(g)]
        if v.size == 0:
            raise FactorError(f"{name}: the raster has no data over the DEM.")
        vmin, vmax = float(v.min()), float(np.percentile(v, 99.5))
        what = "raster"
    elif value is not None:
        vmin = vmax = float(value)
        v = np.array([vmin])
        what = "value"
    else:
        return []
    warnings = []
    if vmin < lo - 1e-9:
        raise FactorError(f"{name}: negative values in the {what} ({vmin:g}).")
    if name in ("C", "P") and vmax > hi + 1e-6:
        sub = v[:: max(1, v.size // 200000)]
        codes = np.isin(sub, sorted(WORLDCOVER_CODES))       # exact class codes
        sample = np.unique(sub[codes])
        # a class raster resampled bilinearly still has mostly exact codes
        if grid is not None and codes.mean() >= 0.8:
            raise FactorError(
                f"{name}: the raster holds land-cover class codes "
                f"({', '.join(str(int(x)) for x in sample[:8])}), not {name} values (0-1). "
                f"Give it as the 'ESA WorldCover' input instead - QEHT then applies the "
                f"class lookup with nearest-neighbour resampling.")
        raise FactorError(f"{name}: values up to {vmax:g}, but {name} lies between 0 and 1. "
                          "Check the raster / value.")
    if name == "K" and vmax > hi:
        if vmax <= 0.8:
            raise FactorError(
                f"K: values up to {vmax:g} look like US customary units "
                "(t.ac.h/(100.ac.ft.tonf.in)). Multiply by 0.1317 for SI t.ha.h/(ha.MJ.mm).")
        raise FactorError(f"K: values up to {vmax:g} are far outside 0-0.1 (SI). Check the units.")
    if name == "R" and vmax > hi:
        raise FactorError(f"R: values up to {vmax:g} exceed any observed annual erosivity "
                          f"(< {hi:g} MJ.mm/(ha.h.yr)). Check the units.")
    if name == "R" and float(np.median(v)) < 200:
        warnings.append(f"R: median {float(np.median(v)):g} MJ.mm/(ha.h.yr) is very low for "
                        "most climates (Kenya about 1,000-6,000). US units? multiply by 17.02. "
                        "Soil loss scales directly with R.")
    if name == "K" and vmax < 0.001:
        warnings.append(f"K: values up to {vmax:g} are unusually small for SI units.")
    return warnings


def c_from_worldcover(classes, valid=None, lookup=None):
    """C (or, with lookup=WORLDCOVER_P, P) grid from a WorldCover class grid.

    Unknown classes -> NaN and are reported.
    """
    lut = dict(WORLDCOVER_C) if lookup is None else dict(lookup)
    cls = np.asarray(classes)
    c = np.full(cls.shape, np.nan, dtype=np.float64)
    for code, row in lut.items():
        val = row[0] if isinstance(row, (tuple, list)) else float(row)
        c[cls == code] = val
    if valid is not None:
        c[~valid] = np.nan
    present = np.unique(cls[np.isfinite(c)]) if c.size else []
    unknown = sorted(set(np.unique(cls[~np.isfinite(c) & (valid if valid is not None else True)]).tolist())
                     - {0})
    return c, {"classes_present": [int(x) for x in present], "unknown_classes": unknown}


def rusle(ls, valid, R, K, C, P=None):
    """Soil loss A = R K LS C P (t/ha/yr) or None when R, K or C is missing.

    Returns (A or None, info) where info lists every factor's provenance and
    whether the LS-only fallback applies.
    """
    P = P if P is not None and P.present else Factor("P", value=1.0,
                                                     source="default 1.0 (no support practice)")
    factors = [R, K, C, P]
    info = {"factors": [f.describe() for f in factors]}
    missing = [f.name for f in (R, K, C) if not f.present]
    if missing:
        info["mode"] = "LS-only terrain potential"
        info["note"] = ("No soil loss computed: factor(s) " + ", ".join(missing) + " missing. "
                        "Classes use LS percentiles ('terrain potential'), not t/ha/yr.")
        return None, info
    shape = np.shape(ls)
    with np.errstate(invalid="ignore"):
        a = (R.as_grid(shape) * K.as_grid(shape) * np.asarray(ls, dtype=np.float64)
             * C.as_grid(shape) * P.as_grid(shape))
    a[~valid] = np.nan
    info["mode"] = "RUSLE"
    info["proxy_factors"] = [f.name for f in factors if f.proxy]
    return a, info
