# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Curve number and Rational C per catchment (F1, v0.15).

Per cell: land-cover class (ESA WorldCover, warped to the DEM grid with
nearest neighbour) x hydrologic soil group (direct HSG input or the
texture/drainage proxy, core.soils.sources) -> CN (AMC II) from an editable
lookup. The catchment value is the area-weighted mean of the cell CNs
(cn_ii); cn_export is that composite converted to the chosen antecedent
moisture condition. Rational C works the same way from a user lookup.

Default CN lookup (a PROXY, flagged as such in the metadata): USDA-NRCS
TR-55 (1986) Table 2-2, each WorldCover class matched to one TR-55 cover by
judgement - check it against local practice before design use:

    10 Tree cover          Woods (condition)
    20 Shrubland           Brush - brush-weed-grass mixture (condition)
    30 Grassland           Pasture, grassland or range (condition)
    40 Cropland            Row crops, straight row (good; poor if 'poor')
    50 Built-up            Commercial and business, 85 % impervious
    60 Bare / sparse       Fallow, bare soil
    80 Permanent water     100 (convention: all rainfall runs off)
    70, 90, 95, 100        no TR-55 equivalent - no CN unless the user's
                           lookup gives one (cells reported, not invented)

Condition good / fair / poor is a run parameter (default fair). Dual HSGs
(A/D, B/D, C/D) are used as D (undrained).

AMC conversions (Chow, Maidment & Mays 1988, eq. 5.5.4-5.5.5):
    CN_I   = 4.2 CN_II / (10 - 0.058 CN_II)
    CN_III = 23 CN_II / (10 + 0.13 CN_II)

Rational C: no default ships (no invented values). Give a lookup CSV with
columns class, A, B, C, D (or class, C for one value per class).

No QGIS imports, no GDAL.
"""

import csv
import json
import math

import numpy as np

TR55 = {  # (cover, {condition: (A, B, C, D)})
    10: ("Woods", {"poor": (45, 66, 77, 83), "fair": (36, 60, 73, 79), "good": (30, 55, 70, 77)}),
    20: ("Brush", {"poor": (48, 67, 77, 83), "fair": (35, 56, 70, 77), "good": (30, 48, 65, 73)}),
    30: ("Pasture, grassland or range",
         {"poor": (68, 79, 86, 89), "fair": (49, 69, 79, 84), "good": (39, 61, 74, 80)}),
    40: ("Row crops, straight row",
         {"poor": (72, 81, 88, 91), "fair": (67, 78, 85, 89), "good": (67, 78, 85, 89)}),
    50: ("Commercial and business, 85% impervious", {"any": (89, 92, 94, 95)}),
    60: ("Fallow, bare soil", {"any": (77, 86, 91, 94)}),
    80: ("Open water (convention)", {"any": (100, 100, 100, 100)}),
}
LC_GROUPS = {"tree": (10,), "shrub": (20,), "grass": (30,), "crop": (40,), "built": (50,),
             "bare": (60,), "water": (80,), "wetland": (90, 95), "other": (70, 100)}
RUNOFF_FIELDS = (["cn_ii", "cn_amc", "cn_export", "cn_coverage_pct", "rational_c",
                  "rc_coverage_pct"] + [f"lc_pct_{g}" for g in LC_GROUPS]
                 + ["cn_lookup_id", "lc_dataset"])
CONDITIONS = ("fair", "good", "poor")
AMC = ("II", "III", "I")


def tr55_lookup(condition="fair"):
    """{class: (A, B, C, D)} from TR-55 for the condition."""
    if condition not in CONDITIONS:
        raise ValueError(f"condition must be one of {CONDITIONS}")
    out = {}
    for cls, (_, table) in TR55.items():
        out[cls] = table.get(condition) or table["any"]
    return out


def lookup_id(condition, user_path=None):
    if user_path:
        return f"user lookup {user_path}"
    return f"TR-55 (1986) Table 2-2, {condition} condition, matched to ESA WorldCover (QEHT proxy)"


def read_lookup_csv(path, value_name="CN"):
    """Lookup CSV -> {class: (A, B, C, D)}. Columns: class, A, B, C, D, or
    class and one value column (same value for every HSG)."""
    with open(path, newline="", encoding="utf-8-sig") as f:
        sample = f.read(4096)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        except csv.Error:
            dialect = csv.excel
        rows = list(csv.DictReader(f, dialect=dialect))
    if not rows:
        raise ValueError(f"{path}: no rows")
    cols = {k.strip().lower(): k for k in rows[0] if k}
    if "class" not in cols:
        raise ValueError(f"{path}: needs a 'class' column")
    groups = [g for g in ("a", "b", "c", "d") if g in cols]
    if len(groups) == 4:
        single = []
    else:                       # e.g. "class, C" for a Rational C table
        single = [k for k in cols if k not in ("class", "name", "description")]
    if len(groups) != 4 and not single:
        raise ValueError(f"{path}: give columns A, B, C, D or one {value_name} column")
    out = {}
    for r in rows:
        try:
            cls = int(float(r[cols["class"]]))
        except (TypeError, ValueError):
            continue
        def num(k):
            try:
                v = float(r[cols[k]])
                return v if math.isfinite(v) else None
            except (TypeError, ValueError, KeyError):
                return None
        if len(groups) == 4:
            out[cls] = tuple(num(g) for g in ("a", "b", "c", "d"))
        else:
            v = num(single[0])
            out[cls] = (v, v, v, v)
    return out


def value_grid(classes, hsg, lookup):
    """Per-cell value from {class: (A, B, C, D)}; NaN where class or HSG unknown.

    hsg codes: 1..4 = A..D, 14 = dual group (used as D), 0 = unknown.
    """
    cls = np.nan_to_num(np.asarray(classes, float), nan=-1).astype(np.int64)
    h = np.asarray(hsg).astype(np.int64)
    h = np.where(h == 14, 4, h)
    out = np.full(cls.shape, np.nan)
    for c, row in lookup.items():
        if row is None:
            continue
        for k in range(4):
            v = row[k]
            if v is None:
                continue
            out[(cls == c) & (h == k + 1)] = float(v)
    return out


def amc_convert(cn_ii, amc="II"):
    if cn_ii is None or not math.isfinite(cn_ii):
        return None
    if amc == "III":
        return 23.0 * cn_ii / (10.0 + 0.13 * cn_ii)
    if amc == "I":
        return 4.2 * cn_ii / (10.0 - 0.058 * cn_ii)
    return cn_ii


def runoff_block(mask, classes=None, cn=None, rc=None, amc="II", lookup_name="",
                 lc_dataset=""):
    """RUNOFF_FIELDS for one catchment mask."""
    out = {k: None for k in RUNOFF_FIELDS}
    mask = np.asarray(mask, bool)
    n = int(mask.sum())
    if n == 0:
        return out
    if classes is not None:
        cl = np.nan_to_num(np.asarray(classes, float)[mask], nan=-1)
        known = cl > 0
        nk = int(known.sum())
        if nk:
            for g, codes in LC_GROUPS.items():
                out[f"lc_pct_{g}"] = float(100.0 * np.isin(cl[known], codes).sum() / nk)
        out["lc_dataset"] = lc_dataset
    if cn is not None:
        v = np.asarray(cn, float)[mask]
        ok = np.isfinite(v)
        out["cn_coverage_pct"] = float(100.0 * ok.sum() / n)
        if ok.any():
            out["cn_ii"] = float(v[ok].mean())
            out["cn_amc"] = amc
            out["cn_export"] = amc_convert(out["cn_ii"], amc)
        out["cn_lookup_id"] = lookup_name
    if rc is not None:
        v = np.asarray(rc, float)[mask]
        ok = np.isfinite(v)
        out["rc_coverage_pct"] = float(100.0 * ok.sum() / n)
        if ok.any():
            out["rational_c"] = float(v[ok].mean())
    return out


class RunoffInputs(object):
    """Per-cell land cover, CN and Rational C on the DEM grid + provenance."""

    def __init__(self, classes, hsg=None, condition="fair", amc="II", cn_lookup=None,
                 cn_lookup_path=None, rc_lookup=None, rc_lookup_path=None, lc_dataset=""):
        self.classes = np.asarray(classes, float)
        self.amc = amc
        self.lc_dataset = lc_dataset
        self.cn_lookup = cn_lookup or tr55_lookup(condition)
        self.cn_lookup_name = lookup_id(condition, cn_lookup_path)
        self.cn = value_grid(self.classes, hsg, self.cn_lookup) if hsg is not None else None
        self.rc = None
        self.note = ""
        if rc_lookup:
            per_group = any(r is not None and len(set(r)) > 1 for r in rc_lookup.values())
            if hsg is None and per_group:
                self.note = ("The Rational C lookup varies by soil group but no soil data was "
                             "given - rational_c left empty.")
            else:
                self.rc = value_grid(self.classes, hsg if hsg is not None else
                                     np.ones(self.classes.shape, np.int8), rc_lookup)
        self.rc_lookup_name = f"user lookup {rc_lookup_path}" if rc_lookup else ""
        self.meta = {"cn_lookup_id": self.cn_lookup_name, "cn_condition": condition,
                     "cn_amc": amc, "cn_proxy": not bool(cn_lookup_path),
                     "cn_lookup": {str(k): v for k, v in self.cn_lookup.items()},
                     "rc_lookup_id": self.rc_lookup_name,
                     "rc_lookup": {str(k): v for k, v in (rc_lookup or {}).items()},
                     "lc_dataset": lc_dataset}

    def block(self, mask):
        return runoff_block(mask, self.classes, self.cn, self.rc, self.amc,
                            self.cn_lookup_name, self.lc_dataset)

    def meta_json(self):
        return json.dumps(self.meta, sort_keys=True)
