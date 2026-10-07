# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Standalone catchment characteristics table (§0 of the F1-F8 plan, v0.17).

One row per crossing, joined on outlet_uid from the crossings, catchments
and flow paths of a design hydrology package: the table engineers otherwise
assemble by hand in a spreadsheet. Read from the package with no
recalculation. Optional groups (soils, curve number, rainfall, erosion,
flat-method check) appear only when at least one crossing has a value.

No QGIS imports, no GDAL.
"""

import csv
import math
import os

from ..interop import gpkg

# (column, source layer) in table order; groups after the first are optional
CORE = [
    ("outlet_uid", "crossings"), ("status", "crossings"), ("chainage_m", "crossings"),
    ("outlet_x", "crossings"), ("outlet_y", "crossings"), ("stream_order", "crossings"),
    ("area_km2", "catchments"), ("elev_max_m", "catchments"), ("elev_min_m", "catchments"),
    ("relief_m", "catchments"), ("catch_slope_horn", "catchments"),
    ("catch_relief_ratio", "catchments"),
    ("lfp_length_m", "flowpaths"), ("lfp_drop_m", "flowpaths"), ("lfp_slope", "flowpaths"),
    ("lfp_slope_1085", "flowpaths"),
]
GROUPS = [
    ("flow path segments", [("lfp_overland_m", "flowpaths"), ("lfp_overland_slope", "flowpaths"),
                            ("lfp_sheet_m", "flowpaths"), ("lfp_shallow_m", "flowpaths"),
                            ("lfp_channel_m", "flowpaths"), ("lfp_channel_slope", "flowpaths"),
                            ("lfp_channel_slope_1085", "flowpaths")]),
    ("channel at the crossing", [("ch_slope_us", "crossings"), ("ch_slope_ds", "crossings")]),
    ("flat-method check", [("flat_sensitivity_pct", "extra"), ("flat_sensitive", "extra")]),
    ("soils", [("soil_texture", "catchments"), ("soil_hsg", "catchments"),
               ("usle_k", "catchments"), ("soil_coverage_pct", "catchments")]),
    ("curve number", [("cn_ii", "catchments"), ("cn_amc", "catchments"),
                      ("cn_export", "catchments"), ("rational_c", "catchments")]),
    ("rainfall", [("rain_zone", "catchments"), ("rain_zone_pct", "catchments"),
                  ("map_mm", "catchments"), ("rusle_r", "catchments")]),
    ("shape and network", [("perimeter_km", "catchments"), ("form_factor", "catchments"),
                           ("elongation_ratio", "catchments"), ("circularity_ratio", "catchments"),
                           ("drainage_density", "catchments"), ("stream_frequency", "catchments"),
                           ("max_strahler", "catchments")]),
    ("erosion", [("ero_a_mean_tha", "catchments"), ("ero_rusle_class", "catchments"),
                 ("ero_sy_m3_yr", "catchments"), ("ero_impact", "crossings"),
                 ("ero_dep_flag", "crossings")]),
]


def _has(v):
    return v is not None and not (isinstance(v, float) and not math.isfinite(v)) and v != ""


def characteristics_rows(package_path, extra=None):
    """-> (columns, rows). extra: {outlet_uid: {field: value}} (e.g. the
    flat-method check), joined like a layer."""
    tables = dict(gpkg.list_tables(package_path))
    layers = {}
    for name in ("crossings", "catchments", "flowpaths"):
        if name in tables:
            layers[name] = {r["outlet_uid"]: r for r in
                            gpkg.read_table(package_path, name, with_geometry=False)}
    layers["extra"] = dict(extra or {})
    order = list(layers.get("crossings", {}))
    if not order:
        return [c for c, _ in CORE], []

    def value(uid, col, src):
        return layers.get(src, {}).get(uid, {}).get(col)

    cols = list(CORE)
    for _, group in GROUPS:
        if any(_has(value(u, c, s)) for u in order for c, s in group):
            cols += group
    if all(not _has(value(u, "chainage_m", "crossings")) for u in order):
        cols = [cs for cs in cols if cs[0] != "chainage_m"]
    rows = []
    for u in order:
        rows.append({c: value(u, c, s) for c, s in cols})

    def key(r):                      # along the road when chainage exists, else as stored
        ch = r.get("chainage_m")
        return (0 if _has(ch) else 1, ch if _has(ch) else 0)
    rows.sort(key=key)
    return [c for c, _ in cols], rows


def write_characteristics_csv(package_path, out_csv, extra=None):
    cols, rows = characteristics_rows(package_path, extra)
    os.makedirs(os.path.dirname(os.path.abspath(out_csv)), exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(cols)
        for r in rows:
            wr.writerow(["" if not _has(r[c]) else
                         (round(r[c], 6) if isinstance(r[c], float) else r[c]) for c in cols])
    return out_csv, len(rows)
