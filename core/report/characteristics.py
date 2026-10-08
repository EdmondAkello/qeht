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
    ("time of concentration", [("tc_kirpich_min", "crossings"), ("tc_kerby_kirpich_min", "crossings"),
                               ("tc_scs_lag_min", "crossings"), ("tc_tr55_min", "crossings"),
                               ("tc_bransby_williams_min", "crossings")]),
    ("DEM uncertainty", [("unc_area_p10", "crossings"), ("unc_area_p90", "crossings"),
                         ("unc_area_cv", "crossings"), ("unc_switch_pct", "crossings"),
                         ("unc_lost_pct", "crossings")]),
    ("land-cover scenario", [("scenario_name", "catchments"), ("cn_ii_scn", "catchments"),
                             ("d_cn", "catchments"), ("d_rational_c", "catchments"),
                             ("d_ero_a_mean_tha", "catchments")]),
    ("mapped drainage", [("map_agrees", "crossings"), ("map_river_dist_m", "crossings"),
                         ("map_precision", "catchments"), ("map_recall", "catchments")]),
    ("channel section", [("xs_bankfull_w_m", "crossings"), ("xs_bankfull_d_m", "crossings"),
                         ("xs_w_1p0_m", "crossings"), ("xs_quality", "crossings")]),
    ("floodplain width", [("fp_w_0p5_m", "crossings"), ("fp_w_1p0_m", "crossings"),
                          ("fp_w_2p0_m", "crossings"), ("fp_hand_w_1p0_m", "crossings"),
                          ("fp_note", "crossings")]),
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


TC_COLUMNS = ["outlet_uid", "status", "chainage_m", "area_km2", "lfp_length_m", "lfp_slope_1085",
              "tc_kirpich_min", "tc_kirpich_flag", "tc_kerby_kirpich_min", "tc_kerby_kirpich_flag",
              "tc_scs_lag_min", "tc_scs_lag_flag", "tc_tr55_min", "tc_tr55_flag",
              "tc_bransby_williams_min", "tc_bransby_williams_flag", "tc_note"]


def tc_rows(package_path):
    """Time of concentration table (F10): one row per crossing, every method side by side."""
    tables = dict(gpkg.list_tables(package_path))
    if "crossings" not in tables:
        return TC_COLUMNS, []
    lay = {n: {r["outlet_uid"]: r for r in gpkg.read_table(package_path, n, with_geometry=False)}
           for n in ("crossings", "catchments", "flowpaths") if n in tables}
    rows = []
    for uid, x in lay["crossings"].items():
        r = {}
        for c in TC_COLUMNS:
            for n in ("crossings", "catchments", "flowpaths"):
                v = lay.get(n, {}).get(uid, {}).get(c)
                if _has(v):
                    r[c] = v
                    break
            else:
                r[c] = None
        rows.append(r)
    rows.sort(key=lambda r: (0 if _has(r.get("chainage_m")) else 1, r.get("chainage_m") or 0))
    return TC_COLUMNS, rows


def write_tc_csv(package_path, out_csv):
    cols, rows = tc_rows(package_path)
    os.makedirs(os.path.dirname(os.path.abspath(out_csv)), exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(cols)
        for r in rows:
            wr.writerow(["" if not _has(r[c]) else
                         (round(r[c], 6) if isinstance(r[c], float) else r[c]) for c in cols])
    return out_csv, len(rows)
