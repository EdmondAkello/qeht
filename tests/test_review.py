# -*- coding: utf-8 -*-
"""1.2 review tools: fill-depression table, review register, crossing sheets.
Bare Python + NumPy, no QGIS, no GDAL.

    python -m qeht.tests.test_review
"""

import csv
import os
import sys
import tempfile

import numpy as np

from .test_interop import synthetic_run, UTM37S_WKT
from ..core.interop.heas_exchange import write_exchange
from ..core.network.alignment import Alignment
from ..core.conditioning.fill import fill_depressions
from ..core.review.depressions import (depression_table, label_depressions, fill_depth,
                                       FLAG_AT_ROAD, FLAG_NEAR_ROAD)
from ..core.review.register import build_register, register_counts, write_register_csv
from ..core.review.sheets import build_sheets, alignment_from_package

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def _pits():
    """Tilted 80 x 100 plane (10 m cells) with three pits and a road along row 40."""
    rows, cols, cs = 80, 100, 10.0
    r, c = np.mgrid[0:rows, 0:cols]
    dem = 100.0 + 0.02 * cs * r + 0.005 * cs * c
    dem[38:41, 20:24] -= 3.0          # A: 3 m pit on the road (rows 38-40)
    dem[50:52, 60:62] -= 1.5          # B: 1.5 m pit 100-120 m off the road
    dem[10:12, 80:82] -= 0.3          # C: shallow pit, far away
    dem[70, 5] -= 2.0                 # D: deep single-cell pit, far from the road
    valid = np.ones(dem.shape, bool)
    gt = (500000.0, cs, 0.0, 9900000.0, 0.0, -cs)
    filled, _, _ = fill_depressions(dem, valid, cell_width=cs, cell_height=cs)
    y = gt[3] - 39.5 * cs
    road = Alignment([[(gt[0], y), (gt[0] + cols * cs, y)]], start_chainage=1000.0)
    return dem, filled, valid, gt, road


def test_depressions():
    print("\n1. Fill-depression table")
    dem, filled, valid, gt, road = _pits()
    depth = fill_depth(dem, filled, valid)
    labels, n = label_depressions(depth)
    check("four separate depressions labelled (8-connected)", n == 4, f"n={n}")
    rows, _, summ = depression_table(dem, filled, valid, gt, alignment=road,
                                     crossings=[("X001", 1220.0), ("X002", 1900.0)],
                                     min_depth=0.5, deep_m=1.0, near_m=150.0)
    check("min depth drops the 0.3 m pit", len(rows) == 3 and summ["reported"] == 3,
          str([round(r["max_depth_m"], 2) for r in rows]))
    a = rows[0]
    check("deepest first; A (12 cells) filled 2.5-3 m (spill on the tilted plane)",
          a["cells"] == 12 and 2.5 <= a["max_depth_m"] <= 3.0, f"{a['max_depth_m']:.2f}")
    check("area = cells x cell area", a["area_m2"] == a["cells"] * 100.0 and a["cells"] == 12)
    check("volume = sum of depth x cell area",
          abs(a["volume_m3"] - float(depth[38:41, 20:24].sum()) * 100.0) < 1e-6)
    check("A touches the road and is flagged deep_at_road",
          a["touches_road"] == 1 and a["flag"] == FLAG_AT_ROAD)
    check("A chainage from the deepest cell (start 1000)", 1195 <= a["chainage_m"] <= 1245,
          f"{a['chainage_m']:.0f}")
    check("A nearest crossing X001 within ~25 m of chainage",
          a["nearest_uid"] == "X001" and a["nearest_ch_dist_m"] < 30, str(a["nearest_ch_dist_m"]))
    b = next(r for r in rows if r["cells"] == 4)
    check("B 100+ m off the road is deep_near_road with near 150 m",
          b["flag"] == FLAG_NEAR_ROAD and 90 <= b["dist_road_m"] <= 125, f"{b['dist_road_m']:.1f}")
    d = next(r for r in rows if r["cells"] == 1)
    check("D far from the road is not flagged", d["flag"] == "" and d["dist_road_m"] > 150)
    check("polygons traced, area matches", abs(sum(
        abs(sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(ring, ring[1:]))) / 2
        for poly in a["_geom"] for ring in poly[:1]) - a["area_m2"]) < 1e-6)
    rows2, _, _ = depression_table(dem, filled, valid, gt, min_depth=0.5)
    check("no road: no distances or flags", all(r["dist_road_m"] is None and r["flag"] == ""
                                                for r in rows2))
    flat = np.where(valid, 100.0, np.nan)
    rows3, _, s3 = depression_table(flat, flat, valid, gt)
    check("no depressions: empty table", rows3 == [] and s3["depressions"] == 0)
    rows4, _, s4 = depression_table(dem, filled, valid, gt, alignment=road, min_depth=0.5,
                                    max_rows=1, with_polygons=False)
    check("max_rows keeps the deepest; polygons optional",
          len(rows4) == 1 and "_geom" not in rows4[0] and s4["reported"] == 1)
    return rows


def _package(path, with_road=True, truncate=False):
    (cr, ca, fp, issues, id_info), g = synthetic_run()
    for k, (_, a) in enumerate(cr):
        a["chainage_m"] = [300.0, 100.0, 200.0][k]
    for k, (_, a) in enumerate(ca):        # the synthetic surface reaches its edges
        a["truncated"] = int(truncate and k == 0)
        a["truncated_cells"] = 7 if truncate and k == 0 else 0
    cr[1][1]["unc_switch_pct"] = 20.0
    cr[1][1]["unc_lost_pct"] = 0.0
    cr[2][1]["map_agrees"] = 0
    cr[2][1]["map_river_dist_m"] = 140.0
    gt = g["gt"]
    extra = []
    if with_road:
        y = gt[3] - 40 * 10.0
        extra.append(("road_alignment", "LINESTRING",
                      [("part", "int"), ("chainage_from_m", "real"), ("chainage_to_m", "real")],
                      [([(gt[0], y), (gt[0] + 500.0, y)],
                        {"part": 1, "chainage_from_m": 50.0, "chainage_to_m": 550.0})],
                      "Road alignment used for chainage"))
    md = dict(id_info, dem_path="synthetic", raw_dem_path="synthetic", dem_sha256="",
              dem_source="synthetic", flat_method="barnes", stream_threshold_km2="0.01",
              cell_size_m="10 x 10", qeht_version="test")
    write_exchange(path, cr, ca, fp, md, UTM37S_WKT, crs_name="WGS 84 / UTM zone 37S",
                   extra_layers=extra)
    return path, {a["outlet_uid"]: a["chainage_m"] for _, a in cr}


def test_register(dep_rows):
    print("\n2. Review register")
    tmp = tempfile.mkdtemp()
    pkg, ch = _package(os.path.join(tmp, "p.gpkg"), truncate=True)
    uids = sorted(ch)
    flat = {uids[0]: {"flat_sensitive": True, "flat_sensitivity_pct": 14.0,
                      "area_barnes_km2": 0.2, "area_toward_km2": 0.23},
            uids[1]: {"flat_moved": True, "area_barnes_km2": 0.1, "area_toward_km2": 0.02}}
    extra = {"flat_check": flat, "dem_audit": ["Vertical step 1 m: whole-metre elevations."],
             "warnings": ["DEM check: Vertical step 1 m: whole-metre elevations.",
                          "2 crossing(s) change area by more than 10 % between the two flat "
                          "methods: X001", "Exports failed: test"]}
    reg = build_register(pkg, extra, depressions=dep_rows, near_m=100.0)
    ids = [r["finding_id"] for r in reg]
    check("IDs R-001.. in order, unique", ids == [f"R-{i:03d}" for i in range(1, len(reg) + 1)])
    sev = [r["severity"] for r in reg]
    check("sorted by severity", sev == sorted(sev, key=["high", "medium", "low", "info"].index), str(sev))
    cats = {(r["category"], r["outlet_uid"]) for r in reg}
    trunc_uid = next(u for u in uids if any(r["category"] == "catchment" and r["outlet_uid"] == u
                                            for r in reg))
    check("truncated catchment is high", any(r["category"] == "catchment" and r["severity"] == "high"
                                             for r in reg), trunc_uid)
    check("flat_moved high, flat_sensitive medium",
          ("flat routing", uids[1]) in cats and ("flat routing", uids[0]) in cats
          and next(r for r in reg if r["category"] == "flat routing"
                   and r["outlet_uid"] == uids[1])["severity"] == "high")
    check("uncertainty switching listed", any(r["category"] == "DEM uncertainty" for r in reg))
    check("mapped drainage disagreement listed",
          any(r["category"] == "mapped drainage" and "140 m" in r["evidence"] for r in reg))
    check("DEM audit listed once (log copy dropped)",
          sum(1 for r in reg if "whole-metre" in r["evidence"]) == 1)
    check("flat warning not duplicated as a run note",
          not any(r["category"] == "run" and "flat methods" in r["evidence"] for r in reg))
    check("other warnings kept as info", any(r["category"] == "run" and r["severity"] == "info"
                                             and "Exports failed" in r["evidence"] for r in reg))
    fills = [r for r in reg if r["category"] == "fill"]
    check("flagged depressions listed (2 of 3)", len(fills) == 2, str(len(fills)))
    check("status open, response empty", all(r["status"] == "open" and r["response"] == "" for r in reg))
    n = register_counts(reg)
    check("counts add up", sum(n.values()) == len(reg), str(n))
    out = write_register_csv(reg, os.path.join(tmp, "reg.csv"))
    with open(out, encoding="utf-8-sig") as f:
        back = list(csv.DictReader(f))
    check("CSV round trip", len(back) == len(reg) and back[0]["finding_id"] == "R-001")
    pkg2, _ = _package(os.path.join(tmp, "p2.gpkg"))
    reg2 = build_register(pkg2, {})
    check("clean package: only the planted unc / map findings",
          {r["category"] for r in reg2} == {"DEM uncertainty", "mapped drainage"}, str(reg2))
    return pkg, reg, flat


def test_sheets(pkg, reg, flat, dep_rows):
    print("\n3. Crossing sheets")
    al = alignment_from_package(pkg)
    check("alignment rebuilt from the package with its start chainage",
          al is not None and abs(al.start_chainage - 50.0) < 1e-9 and abs(al.length - 500.0) < 1e-9)
    uids = sorted({r["outlet_uid"] for r in reg if r["outlet_uid"]})
    sens = {uids[0]: [{"factor": "snap radius", "setting": "10 cells", "area_km2": 0.21,
                       "d_area_pct": 5.0, "lfp_length_m": 900.0, "d_lfp_pct": 1.0,
                       "outlet_moved": False}]}
    html = build_sheets(pkg, reg, {"flat_check": flat, "sensitivity": sens}, dep_rows, title="T & <x>")
    check("one sheet per crossing", html.count("<section class='sheet'>") == 3)
    check("one plan per crossing", html.count("<svg") == 3)
    check("register findings shown with their IDs", all(r["finding_id"] in html for r in reg
                                                        if r["outlet_uid"]))
    check("sensitivity table rendered", "snap radius" in html and "Sensitivity" in html)
    check("title escaped", "T &amp; &lt;x&gt;" in html and "<x>" not in html)
    check("chainage order (100, 200, 300)",
          html.index("chainage 100") < html.index("chainage 200") < html.index("chainage 300"))
    check("sign-off box on every sheet", html.count("Reviewer, date") == 3)
    tmp = tempfile.mkdtemp()
    pkg2, _ = _package(os.path.join(tmp, "nr.gpkg"), with_road=False)
    check("no road: alignment None, sheets still build",
          alignment_from_package(pkg2) is None and "<section class='sheet'>" in build_sheets(pkg2))


def main():
    rows = test_depressions()
    pkg, reg, flat = test_register(rows)
    test_sheets(pkg, reg, flat, rows)
    print(f"\n{'ALL PASSED' if not FAILURES else 'FAILED: ' + ', '.join(FAILURES)}")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
