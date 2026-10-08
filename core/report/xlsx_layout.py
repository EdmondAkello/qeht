# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Column order of the XLSX workbook (F11, v0.23).

One config list per sheet, so the order can be matched to a house
workbook later without code changes. Each entry:
(field, source layer, header, unit, number format). Source layers are
crossings, catchments, flowpaths and extra (values joined by outlet_uid
from outside the package, such as the flat-method check). Every column is
always written, empty where the package has no value, so the layout stays
fixed from one run to the next.

No QGIS imports, no GDAL.
"""

import json
import math

from ..interop import gpkg

CHARACTERISTICS = [
    # 1. identity and location
    ("outlet_uid", "crossings", "ID", "", None),
    ("status", "crossings", "Status", "", None),
    ("chainage_m", "crossings", "Chainage", "m", "0.0"),
    ("outlet_x", "crossings", "Easting", "m", "0.00"),
    ("outlet_y", "crossings", "Northing", "m", "0.00"),
    # 2. area
    ("area_km2", "catchments", "Area", "km²", "0.0000"),
    # 3. flow path length and drop
    ("lfp_length_m", "flowpaths", "LFP length", "m", "0.0"),
    ("lfp_drop_m", "flowpaths", "LFP drop", "m", "0.00"),
    # 4. elevations and relief
    ("elev_max_m", "catchments", "Elevation max", "m", "0.0"),
    ("elev_min_m", "catchments", "Elevation min", "m", "0.0"),
    ("elev_mean_m", "catchments", "Elevation mean", "m", "0.0"),
    ("relief_m", "catchments", "Relief", "m", "0.0"),
    # 5. slopes
    ("lfp_slope_1085", "flowpaths", "LFP slope 10-85", "m/m", "0.00000"),
    ("lfp_slope", "flowpaths", "LFP slope drop/length", "m/m", "0.00000"),
    ("catch_slope_horn", "catchments", "Catchment slope (Horn)", "m/m", "0.0000"),
    ("catch_relief_ratio", "catchments", "Relief ratio", "m/m", "0.0000"),
    # 6. flow-path segments
    ("lfp_overland_m", "flowpaths", "Overland length", "m", "0.0"),
    ("lfp_overland_slope", "flowpaths", "Overland slope", "m/m", "0.0000"),
    ("lfp_sheet_m", "flowpaths", "Sheet-flow length", "m", "0.0"),
    ("lfp_shallow_m", "flowpaths", "Shallow-flow length", "m", "0.0"),
    ("lfp_channel_m", "flowpaths", "Channel length", "m", "0.0"),
    ("lfp_channel_slope", "flowpaths", "Channel slope", "m/m", "0.00000"),
    # 7. channel slopes at the crossing
    ("ch_slope_us", "crossings", "Channel slope upstream", "m/m", "0.00000"),
    ("ch_slope_ds", "crossings", "Channel slope downstream", "m/m", "0.00000"),
    # 8. time of concentration
    ("tc_kirpich_min", "crossings", "Tc Kirpich", "min", "0.0"),
    ("tc_kerby_kirpich_min", "crossings", "Tc Kerby + Kirpich", "min", "0.0"),
    ("tc_scs_lag_min", "crossings", "Tc SCS lag", "min", "0.0"),
    ("tc_tr55_min", "crossings", "Tc TR-55", "min", "0.0"),
    ("tc_bransby_williams_min", "crossings", "Tc Bransby-Williams", "min", "0.0"),
    # 9. runoff
    ("cn_ii", "catchments", "CN (AMC II)", "", "0.0"),
    ("cn_amc", "catchments", "AMC", "", None),
    ("cn_export", "catchments", "CN exported", "", "0.0"),
    ("rational_c", "catchments", "Rational C", "", "0.00"),
    # 10. rainfall
    ("rain_zone", "catchments", "Rain zone", "", None),
    ("map_mm", "catchments", "Mean annual rainfall", "mm/yr", "0"),
    # 11. soils
    ("soil_texture", "catchments", "Soil texture", "", None),
    ("soil_hsg", "catchments", "HSG", "", None),
    ("usle_k", "catchments", "USLE K", "t.ha.h/(ha.MJ.mm)", "0.0000"),
    # 12. widths
    ("fp_w_1p0_m", "crossings", "Floodplain width (bed + 1 m)", "m", "0.0"),
    ("xs_bankfull_w_m", "crossings", "Channel bank-full width", "m", "0.0"),
    ("xs_bankfull_d_m", "crossings", "Channel bank-full depth", "m", "0.00"),
    ("xs_w_1p0_m", "crossings", "Channel width (bed + 1 m)", "m", "0.0"),
    # 13. flags
    ("flat_sensitive", "extra", "Flat-sensitive", "0/1", "0"),
    ("xs_quality", "crossings", "Channel section quality", "", None),
    ("map_agrees", "crossings", "Agrees with mapped drainage", "0/1", "0"),
]

TC = [
    ("outlet_uid", "crossings", "ID", "", None),
    ("chainage_m", "crossings", "Chainage", "m", "0.0"),
    ("area_km2", "catchments", "Area", "km²", "0.0000"),
    ("lfp_length_m", "flowpaths", "LFP length", "m", "0.0"),
    ("lfp_slope_1085", "flowpaths", "LFP slope 10-85", "m/m", "0.00000"),
    ("tc_kirpich_min", "crossings", "Kirpich", "min", "0.0"),
    ("tc_kirpich_flag", "crossings", "Kirpich validity", "", None),
    ("tc_kerby_kirpich_min", "crossings", "Kerby + Kirpich", "min", "0.0"),
    ("tc_kerby_kirpich_flag", "crossings", "Kerby validity", "", None),
    ("tc_scs_lag_min", "crossings", "SCS lag", "min", "0.0"),
    ("tc_scs_lag_flag", "crossings", "SCS validity", "", None),
    ("tc_tr55_min", "crossings", "TR-55", "min", "0.0"),
    ("tc_tr55_flag", "crossings", "TR-55 validity", "", None),
    ("tc_bransby_williams_min", "crossings", "Bransby-Williams", "min", "0.0"),
    ("tc_bransby_williams_flag", "crossings", "Bransby-Williams validity", "", None),
    ("tc_note", "crossings", "Missing inputs", "", None),
]

COVERAGE = [
    ("issue", "Finding", "", None), ("chainage_m", "Chainage", "m", "0.0"),
    ("chainage_to_m", "To chainage", "m", "0.0"), ("area_km2", "Area", "km²", "0.0000"),
    ("uid", "Proposed / flagged crossing", "", None), ("nearest_uid", "Nearest existing", "", None),
    ("nearest_m", "Distance to it", "m", "0.0"), ("note", "Note", "", None),
]


def _ok(v):
    return v is not None and v != "" and not (isinstance(v, float) and not math.isfinite(v))


def joined(package_path, layout, extra=None):
    """Rows (lists) in layout order, sorted along the road when chainage exists."""
    tabs = dict(gpkg.list_tables(package_path))
    lay = {n: {r["outlet_uid"]: r for r in gpkg.read_table(package_path, n, with_geometry=False)}
           for n in ("crossings", "catchments", "flowpaths") if n in tabs}
    lay["extra"] = dict(extra or {})
    uids = list(lay.get("crossings", {}))
    rows = [[lay.get(src, {}).get(u, {}).get(f) for f, src, *_ in layout] for u in uids]
    chi = [f for f, *_ in layout].index("chainage_m") if any(f == "chainage_m" for f, *_ in layout) else None
    if chi is not None:
        rows.sort(key=lambda r: (0 if _ok(r[chi]) else 1, r[chi] if _ok(r[chi]) else 0))
    return [[v if _ok(v) else None for v in r] for r in rows]


def workbook_sheets(package_path, extra=None):
    """The five sheets of the QEHT workbook as core.report.xlsx.Sheet objects."""
    from .xlsx import Sheet
    tabs = dict(gpkg.list_tables(package_path))
    sheets = []
    for name, layout in (("Catchment characteristics", CHARACTERISTICS),
                         ("Time of concentration", TC)):
        sheets.append(Sheet(name, [h for _, _, h, _, _ in layout], joined(package_path, layout, extra),
                            units=[u for _, _, _, u, _ in layout],
                            formats=[f for *_, f in layout], freeze=(2, 1)))
    cov = gpkg.read_table(package_path, "coverage_check", with_geometry=False) \
        if "coverage_check" in tabs else []
    sheets.append(Sheet("Coverage findings", [h for _, h, _, _ in COVERAGE],
                        [[r.get(f) for f, *_ in COVERAGE] for r in cov],
                        units=[u for _, _, u, _ in COVERAGE], formats=[f for *_, f in COVERAGE],
                        freeze=(2, 0)))
    from ..interop.field_dictionary import METADATA_KEYS
    meaning = dict(METADATA_KEYS)
    md = gpkg.read_table(package_path, "qeht_run_metadata")
    rows = []
    for r in md:
        v = r["value"]
        if isinstance(v, str) and len(v) > 32000:          # Excel cell limit
            v = v[:32000] + " ..."
        rows.append([r["key"], v, meaning.get(r["key"], "")])
    sheets.append(Sheet("Run metadata", ["Key", "Value", "Meaning"], rows, widths=[28, 80, 60],
                        freeze=(1, 0)))
    fdict = gpkg.read_table(package_path, "qeht_field_dictionary") \
        if "qeht_field_dictionary" in tabs else []
    cols = ["layer", "field", "type", "unit", "meaning", "method", "downstream_use"]
    sheets.append(Sheet("Field dictionary", [c.replace("_", " ").capitalize() for c in cols],
                        [[r.get(c) for c in cols] for r in fdict],
                        widths=[16, 26, 6, 12, 50, 60, 40], freeze=(1, 2)))
    return sheets


def write_workbook(package_path, out_path, extra=None, title=""):
    from .xlsx import write_xlsx
    return write_xlsx(out_path, workbook_sheets(package_path, extra), title=title)


def layout_json():
    return json.dumps([f for f, *_ in CHARACTERISTICS])
