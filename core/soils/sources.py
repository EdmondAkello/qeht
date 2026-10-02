# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Soils from any source on one per-cell model (G1-G4, v0.15).

Every soil source becomes a SoilGrid: per-cell topsoil properties on the
DEM grid, plus provenance. The catchment soil block (SOIL_FIELDS) is then
computed the same way whatever the source, so the design hydrology package
fields mean the same thing for SOTWIS in Kenya, SoilGrids anywhere, HWSD v2
or a national soil map.

Sources
-------
* unit-based (SOTWIS database, polygon attributes, polygons + CSV table,
  unit-code raster + CSV): each unit's description (component-weighted K,
  texture, drainage, HSG proxy) is painted onto the cells of that unit;
  the catchment block is identical to the v0.11-0.14 unit-weighted block.
* texture rasters (e.g. SoilGrids 2.0): sand, clay, (silt), organic carbon,
  optional bulk density and coarse fragments per cell, depth-weighted over
  the requested interval when several depth layers are given. K, texture
  and the HSG proxy are computed per cell.
* direct overrides: a hydrologic soil group raster/field (e.g. HYSOGs250m)
  replaces the texture proxy for the HSG shares; a K raster/field replaces
  the computed K. Precedence: direct input, then computed from texture,
  then empty. Each catchment records which one it used.

SoilGrids 2.0 units (ISRIC; divide the stored value by the factor):
    sand, silt, clay  g/kg      /10  -> %
    soc               dg/kg     /10  -> g/kg  (/100 -> %)
    bdod              cg/cm3    /100 -> g/cm3
    cfvo              cm3/dm3   /10  -> vol %
Depth intervals 0-5, 5-15, 15-30, 30-60, 60-100, 100-200 cm.

HYSOGs250m (Ross et al. 2018, ORNL DAAC) codes: 1-4 = A-D; 11-14 = A/D,
B/D, C/D, D/D (high runoff unless drained); 255 = NoData. Dual groups are
counted as D (undrained, the NRCS convention for design) and their share is
reported separately (hsg_pct_dual).

No QGIS imports. GDAL only through the `read` callables passed in.
"""

import math
import os
import re

import numpy as np

from .usle_k import (williams_k_arr, dg_k_arr, usda_texture_arr, hsg_proxy_arr, TEXTURES,
                     HSG_CODES, usda_texture, hsg_proxy, cfrg)

NAN = float("nan")
PROPS = ("sand", "silt", "clay", "oc_pct", "bulk", "cfrag")

SOILGRIDS = "SoilGrids 2.0 (ISRIC)"
SOILGRIDS_FACTORS = {"sand": 10.0, "silt": 10.0, "clay": 10.0, "soc": 100.0,   # dg/kg -> %
                     "bdod": 100.0, "cfvo": 10.0}
SOILGRIDS_PROP = {"sand": "sand", "silt": "silt", "clay": "clay", "soc": "oc_pct",
                  "bdod": "bulk", "cfvo": "cfrag"}
_SG_NAME = re.compile(r"^(sand|silt|clay|soc|bdod|cfvo)_(\d+)-(\d+)cm_(mean|Q0\.5)\.(tif|tiff|vrt)$",
                      re.I)

HYSOGS = "HYSOGs250m (Ross et al. 2018, ORNL DAAC)"
HYSOGS_CODES = {1: 1, 2: 2, 3: 3, 4: 4, 11: 14, 12: 14, 13: 14, 14: 14}  # -> 1..4, 14 = dual (D)
_HSG_TEXT = {"A": 1, "B": 2, "C": 3, "D": 4, "A/D": 14, "B/D": 14, "C/D": 14, "D/D": 14}


def hsg_from_codes(codes, scheme="hysogs"):
    """HSG code grid: 1..4 = A..D, 14 = dual group (counted as D), 0 = unknown.

    scheme "hysogs": HYSOGs250m integer codes (also fits rasters coded 1-4).
    scheme "letters": strings A, B, C, D, A/D ... (polygon fields).
    """
    if scheme == "letters":
        a = np.asarray(codes, dtype=object)
        out = np.zeros(a.shape, dtype=np.int8)
        for k, v in _HSG_TEXT.items():
            out[np.vectorize(lambda x: str(x).strip().upper() == k if x is not None else False,
                             otypes=[bool])(a)] = v
        return out
    a = np.asarray(codes, dtype=float)
    out = np.zeros(a.shape, dtype=np.int8)
    fin = np.isfinite(a)
    for k, v in HYSOGS_CODES.items():
        out[fin & (a == k)] = v
    return out


def soilgrids_layers(folder):
    """{property: [(top, bottom, path), ...]} for SoilGrids-named files in `folder`
    (e.g. sand_0-5cm_mean.tif). 'mean' is preferred over 'Q0.5'."""
    found = {}
    for name in sorted(os.listdir(folder)):
        m = _SG_NAME.match(name)
        if not m:
            continue
        prop, top, bot, stat = m.group(1).lower(), int(m.group(2)), int(m.group(3)), m.group(4)
        found.setdefault((prop, top, bot), {})[stat.lower()] = os.path.join(folder, name)
    out = {}
    for (prop, top, bot), stats in found.items():
        out.setdefault(prop, []).append((top, bot, stats.get("mean") or stats.get("q0.5")))
    for v in out.values():
        v.sort()
    return out


def depth_weighted(layers, top, bottom, read):
    """Overlap-weighted mean of depth layers [(t, b, path)] over [top, bottom] cm.

    `read(path)` returns a float array (NaN = NoData) on the DEM grid. A cell
    missing in one layer is averaged over the layers it has. Returns
    (array, list of the intervals used) or (None, []) when nothing overlaps.
    """
    num = den = None
    used = []
    for t, b, path in layers:
        ov = min(bottom, b) - max(top, t)
        if ov <= 0:
            continue
        a = np.asarray(read(path), dtype=float)
        ok = np.isfinite(a)
        if num is None:
            num = np.zeros(a.shape); den = np.zeros(a.shape)
        num += np.where(ok, a, 0.0) * ov
        den += ok * ov
        used.append(f"{t}-{b}")
    if num is None:
        return None, []
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan), used


class SoilGrid(object):
    """Per-cell topsoil description on the DEM grid (see the module docstring)."""

    def __init__(self, shape, info=None):
        self.shape = tuple(shape)
        for p in PROPS + ("k_si", "k_dg"):
            setattr(self, p, np.full(self.shape, np.nan))
        self.drain = None                 # array of FAO letters ('' = unknown) or None
        self.hsg = np.zeros(self.shape, dtype=np.int8)   # 1..4, 14 dual(D), 0 unknown
        self.unit_grid = None             # int grid of unit indices (0 = none)
        self.index_to_unit = {}
        self.units = {}
        self.info = dict(info or {})
        self.info.setdefault("soil_dataset", "")
        self.info.setdefault("soil_depth_cm", "")
        self.info.setdefault("hsg_source", "")
        self.info.setdefault("k_source", "")

    # -- construction ---------------------------------------------------------
    @classmethod
    def from_units(cls, unit_grid, index_to_unit, units, info):
        """Paint unit descriptions (core.soils.sotwis.UnitSoil) onto their cells."""
        g = cls(unit_grid.shape, info)
        g.unit_grid = np.asarray(unit_grid)
        g.index_to_unit = dict(index_to_unit)
        g.units = units
        n = int(max(index_to_unit.keys(), default=0)) + 1
        luts = {p: np.full(n, np.nan) for p in PROPS + ("k_si", "k_dg")}
        drain = np.full(n, "", dtype="<U1")
        hsg = np.zeros(n, dtype=np.int8)
        for i, uid in index_to_unit.items():
            u = units.get(uid)
            if u is None or u.k_si is None or not math.isfinite(u.k_si):
                continue      # units without K are "no data", as in the v0.11 block
            for p in PROPS + ("k_si", "k_dg"):
                v = getattr(u, p)
                luts[p][i] = v if v is not None and math.isfinite(v) else np.nan
            drain[i] = (u.drain or "")[:1].upper()
            hsg[i] = {"A": 1, "B": 2, "C": 3, "D": 4}.get(u.hsg or "", 0)
        idx = np.clip(g.unit_grid.astype(np.int64), 0, n - 1)
        inside = (g.unit_grid > 0) & (g.unit_grid < n)
        for p in PROPS + ("k_si", "k_dg"):
            setattr(g, p, np.where(inside, luts[p][idx], np.nan))
        g.drain = np.where(inside, drain[idx], "")
        g.hsg = np.where(inside, hsg[idx], 0).astype(np.int8)
        g.info["k_source"] = g.info.get("k_source") or (
            f"Williams/EPIC from {info.get('soil_dataset', 'soil units')}")
        g.info["hsg_source"] = g.info.get("hsg_source") or "proxy: texture and drainage class"
        return g

    @classmethod
    def from_texture(cls, sand, clay, oc_pct, silt=None, bulk=None, cfrag=None, info=None,
                     drain=None):
        """Per-cell texture grids (percent; oc in percent). Silt = 100 - sand - clay
        when not given. K (Williams/EPIC and Dg), texture and HSG proxy per cell."""
        sand = np.asarray(sand, float)
        g = cls(sand.shape, info)
        clay = np.asarray(clay, float)
        if silt is None:
            silt = np.where(np.isfinite(sand) & np.isfinite(clay),
                            np.clip(100.0 - sand - clay, 0.0, None), np.nan)
        g.sand, g.silt, g.clay = sand, np.asarray(silt, float), clay
        g.oc_pct = np.asarray(oc_pct, float) if oc_pct is not None else np.full(g.shape, np.nan)
        if bulk is not None:
            g.bulk = np.asarray(bulk, float)
        if cfrag is not None:
            g.cfrag = np.asarray(cfrag, float)
        g.k_si = williams_k_arr(g.sand, g.silt, g.clay, g.oc_pct)
        g.k_dg = dg_k_arr(g.sand, g.silt, g.clay)
        g.drain = np.asarray(drain) if drain is not None else None
        g.hsg = hsg_proxy_arr(usda_texture_arr(g.sand, g.silt, g.clay), g.drain)
        g.info["k_source"] = g.info.get("k_source") or "Williams/EPIC per cell from texture rasters"
        g.info["hsg_source"] = g.info.get("hsg_source") or "proxy: texture (per cell)"
        return g

    # -- overrides --------------------------------------------------------------
    def override_k(self, k, source):
        """Direct K (SI) where finite; computed K elsewhere."""
        k = np.asarray(k, float)
        has = np.isfinite(k) & (k >= 0)
        self.k_si = np.where(has, k, self.k_si)
        self._k_direct = has
        self.info["k_source"] = source

    def override_hsg(self, hsg_codes, source):
        """Direct HSG codes (from hsg_from_codes) where known; proxy elsewhere."""
        h = np.asarray(hsg_codes, dtype=np.int8)
        self._hsg_direct = h > 0
        self.hsg = np.where(h > 0, h, self.hsg).astype(np.int8)
        self.info["hsg_source"] = source

    def k_grid(self):
        return self.k_si

    # -- catchment block ---------------------------------------------------------
    def block(self, mask):
        """SOIL_FIELDS for one catchment mask (core.soils.catchment)."""
        from .catchment import SOIL_FIELDS
        out = {f: None for f, _ in SOIL_FIELDS}
        out["soil_dataset"] = self.info.get("soil_dataset")
        out["soil_depth_cm"] = self.info.get("soil_depth_cm")
        mask = np.asarray(mask, bool)
        n = int(mask.sum())
        if n == 0:
            return out
        cov = mask & np.isfinite(self.k_si)
        covered = int(cov.sum())
        out["soil_coverage_pct"] = float(100.0 * covered / n)
        # HSG shares over the catchment cells with a known group (direct or proxy)
        hs = self.hsg[mask]
        known = hs > 0
        if known.any():
            nk = float(known.sum())
            for k, letter in enumerate(HSG_CODES, start=1):
                cnt = np.sum(hs == k) + (np.sum(hs == 14) if letter == "D" else 0)
                out["hsg_pct_" + letter.lower()] = float(100.0 * cnt / nk)
            out["hsg_pct_dual"] = float(100.0 * np.sum(hs == 14) / nk)
            shares = [out["hsg_pct_" + x.lower()] for x in HSG_CODES]
            out["soil_hsg"] = HSG_CODES[int(np.argmax(shares))]
            direct = getattr(self, "_hsg_direct", None)
            out["soil_hsg_source"] = (self.info["hsg_source"] if direct is not None
                                      and bool((direct & mask).any())
                                      else (self.info["hsg_source"] if direct is None
                                            else "proxy: texture and drainage class"))
        if covered == 0:
            return out
        if self.unit_grid is not None:
            ug = self.unit_grid[cov]
            idx, counts = np.unique(ug, return_counts=True)
            uids = [self.index_to_unit.get(int(i)) for i in idx]
            out["soil_units"] = int(len(set(uids)))
            by = {}
            for u, c in zip(uids, counts.tolist()):
                by[u] = by.get(u, 0) + c
            dom = max(by, key=by.get)
            out["soil_dominant_unit"] = str(dom)
            u = self.units.get(dom)
            out["soil_ttr_main"] = getattr(u, "ttr_main", None) if u is not None else None
        for f, attr in (("soil_sand_pct", "sand"), ("soil_silt_pct", "silt"),
                        ("soil_clay_pct", "clay"), ("soil_oc_pct", "oc_pct"),
                        ("soil_cfrag_pct", "cfrag"), ("soil_bulk_gcm3", "bulk"),
                        ("usle_k", "k_si"), ("usle_k_dg", "k_dg")):
            v = getattr(self, attr)[cov]
            v = v[np.isfinite(v)]
            out[f] = float(v.mean()) if v.size else None
        if self.drain is not None:
            d = np.asarray(self.drain)[cov]
            keep = d != ""
            if keep.any():
                if self.unit_grid is not None:
                    # ties resolved in unit order, exactly as the v0.11 block
                    order = self.unit_grid[cov][keep]
                    first = {}
                    for letter, u in zip(d[keep].tolist(), order.tolist()):
                        first.setdefault(letter, u)
                    vals, cnt = np.unique(d[keep], return_counts=True)
                    by = sorted(zip(vals.tolist(), cnt.tolist()), key=lambda t: first[t[0]])
                    out["soil_drain_class"] = max(by, key=lambda t: t[1])[0]
                else:
                    vals, cnt = np.unique(d[keep], return_counts=True)
                    out["soil_drain_class"] = str(vals[int(np.argmax(cnt))])
        if None not in (out["soil_sand_pct"], out["soil_silt_pct"], out["soil_clay_pct"]):
            out["soil_texture"] = usda_texture(out["soil_sand_pct"], out["soil_silt_pct"],
                                               out["soil_clay_pct"])
            out["soil_hsg_proxy"] = hsg_proxy(out["soil_texture"], out["soil_drain_class"])
        if out["soil_cfrag_pct"] is not None:
            out["usle_cfrg"] = cfrg(out["soil_cfrag_pct"])
        kd = getattr(self, "_k_direct", None)
        out["usle_k_source"] = (self.info["k_source"] if kd is None or bool((kd & cov).any())
                                else "Williams/EPIC from texture (no direct K in this catchment)")
        return out
