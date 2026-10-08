# -*- coding: utf-8 -*-
"""v0.24 validation: check against mapped drainage (F12).
Bare Python + NumPy, no QGIS, no GDAL. A synthetic V valley on a 10 m grid
whose stream runs down column 20, and mapped lines placed on it, shifted
off it, or diverted from it.

    python -m qeht.tests.test_mapped
"""

import json
import math
import sys

import numpy as np

from ..core.network.mapped import (agreement, crossing_fields, uncovered_rivers,
                                   divergence_reaches, clip_lines, dilate, params_json,
                                   MAP_FIELDS_X, MAP_FIELDS_C)
from ..core.network.alignment import Alignment
from ..core.flow.direction import d8_direction
from ..core.flow.accumulation import flow_accumulation

FAILURES = []
ROWS, COLS, CS, C0 = 60, 41, 10.0, 20
GT = (500000.0, CS, 0.0, 9900600.0, 0.0, -CS)


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def xy(r, c):
    return GT[0] + (c + 0.5) * CS, GT[3] - (r + 0.5) * CS


def valley():
    r, c = np.mgrid[0:ROWS, 0:COLS].astype(float)
    dem = 200.0 - 0.01 * r * CS + 0.05 * np.abs(c - C0) * CS
    valid = np.ones(dem.shape, bool)
    d, _ = d8_direction(dem, valid, CS, CS)
    acc, _ = flow_accumulation(d, valid)
    stream = (acc >= 30) & valid
    return d, valid, acc, stream


def test_network():
    print("\n1. Network agreement")
    d, valid, acc, stream = valley()
    rs = np.flatnonzero(stream[:, C0])
    check("synthetic DEM stream is column 20 only", stream.sum() == rs.size and rs.size > 20,
          f"{stream.sum()} cells")
    same = [[xy(rs[0], C0), xy(ROWS - 1, C0)]]
    p, r, f, _ = agreement(stream, d, valid, same, GT, 20.0)
    check("mapped line on the DEM stream: precision = recall = F1 = 1",
          math.isclose(p, 1.0) and math.isclose(r, 1.0) and math.isclose(f, 1.0), f"{p}, {r}, {f}")
    for k, (pe, re_) in {2: (1.0, 1.0), 3: (0.0, 0.0)}.items():
        sh = [[(x + k * CS, y) for x, y in same[0]]]
        p, r, f, _ = agreement(stream, d, valid, sh, GT, 20.0)
        check(f"map shifted {k} cells against a 2-cell (20 m) tolerance: P {pe:g}, R {re_:g}",
              math.isclose(p, pe) and math.isclose(r, re_), f"{p:.3f}, {r:.3f}")
    half = [[xy(rs[0], C0), xy(30, C0), xy(30, C0 + 15), xy(ROWS - 1, C0 + 15)]]
    p, r, f, _ = agreement(stream, d, valid, half, GT, 15.0)
    check("map following the stream for half its length: precision about one half, recall lower "
          "(the off-stream part is unmatched)", 0.4 < p < 0.6 and r < p, f"P {p:.3f} R {r:.3f}")
    reg = np.zeros(stream.shape, bool); reg[:30] = True
    p2, r2, _, _ = agreement(stream, d, valid, half, GT, 15.0, region=reg)
    check("restricted to the upper half (a catchment): P = R = 1", math.isclose(p2, 1.0)
          and math.isclose(r2, 1.0), f"{p2}, {r2}")
    m = np.zeros((9, 9), bool); m[4, 4] = True
    dd = dilate(m, 2.0)
    check("dilation is a disk of radius 2 cells (13 cells)", dd.sum() == 13)


def test_crossings_and_road():
    print("\n2. Crossings and road")
    lines = [[xy(0, C0), xy(ROWS - 1, C0)]]
    x0, y0 = xy(40, C0)
    cr = [{"outlet_x": x0 + 10.0, "outlet_y": y0, "acc_at_outlet_km2": 0.5},
          {"outlet_x": x0 + 500.0, "outlet_y": y0, "acc_at_outlet_km2": 0.5},
          {"outlet_x": x0 + 1500.0, "outlet_y": y0, "acc_at_outlet_km2": 0.5},
          {"outlet_x": x0, "outlet_y": y0, "acc_at_outlet_km2": 0.001}]
    b = crossing_fields(cr, lines, ["Mto <A>"], 30.0, 0.01)
    check("10 m off: agrees, distance and name", b[0]["map_agrees"] == 1
          and math.isclose(b[0]["map_river_dist_m"], 10.0) and b[0]["map_river_name"] == "Mto <A>")
    check("500 m off: map_agrees 0 with a note", b[1]["map_agrees"] == 0 and "500 m" in b[1]["map_note"])
    check("beyond 1 km: NULL", b[2]["map_agrees"] is None and "1000 m" in b[2]["map_note"])
    check("below the comparison threshold: NULL with a note",
          b[3]["map_agrees"] is None and "below the comparison threshold" in b[3]["map_note"])
    check("keys", set(b[0]) == set(MAP_FIELDS_X))
    yr = xy(45, 0)[1]
    al = Alignment([[(GT[0], yr), (GT[0] + COLS * CS, yr)]], start_chainage=1000.0)
    f = uncovered_rivers(lines, ["Mto <A>"], al, [1050.0, 1300.0])
    check("mapped river across the road, nearest crossing 100 m away: one finding at its chainage",
          len(f) == 1 and f[0][1]["issue"] == "mapped_river_uncovered"
          and math.isclose(f[0][1]["chainage_m"], 1205.0) and math.isclose(f[0][1]["nearest_m"], 95.0)
          and "Mto <A>" in f[0][1]["note"], str(f[0][1] if f else None))
    check("a crossing within 50 m covers it: no finding",
          uncovered_rivers(lines, ["x"], al, [1180.0]) == [])
    cl = clip_lines([[(GT[0] - 500.0, xy(10, 0)[1]), (GT[0] + 200.0, xy(10, 0)[1])]], ["a"], GT, (ROWS, COLS))
    check("clip to the grid: the part inside, 200 m", len(cl) == 1 and abs(cl[0][1]["length_m"] - 200.0) < 5.0,
          str(cl[0][1] if cl else None))


def test_divergence():
    print("\n3. Divergence reaches")
    d, valid, acc, stream = valley()
    rs = np.flatnonzero(stream[:, C0])
    div = [[xy(rs[0], C0), xy(29, C0), xy(29, C0 + 15), xy(ROWS - 1, C0 + 15)]]
    _, _, _, agree = agreement(stream, d, valid, div, GT, 15.0)
    reaches = divergence_reaches(stream, agree, d, valid, acc, GT, min_cells=10)
    check("one divergence reach where the stream leaves the mapped course (rows 31-59)",
          len(reaches) == 1 and reaches[0][1]["cells"] == 29, str([r[1] for r in reaches]))
    if reaches:
        a = reaches[0][1]
        check("length 28 x 10 m, area = accumulation at its first cell, does not rejoin",
              math.isclose(a["length_m"], 280.0) and math.isclose(a["area_km2"], (acc[31, C0] + 1) * 1e-4)
              and a["rejoins"] == 0, str(a))
    check("min_cells above the run length: none", divergence_reaches(
        stream, agree, d, valid, acc, GT, min_cells=40) == [])
    _, _, _, ag2 = agreement(stream, d, valid, [[xy(0, C0), xy(ROWS - 1, C0)]], GT, 15.0)
    check("full agreement: no divergence", divergence_reaches(stream, ag2, d, valid, acc, GT) == [])


def test_contract():
    print("\n4. Contract")
    from ..core.interop.field_dictionary import field_names, METADATA_KEYS, OPTIONAL_LAYERS
    check("crossing and catchment fields in the contract",
          all(f in field_names("crossings") for f in MAP_FIELDS_X)
          and all(f in field_names("catchments") for f in MAP_FIELDS_C))
    check("optional layers drainage_divergence and mapped_rivers_used; mapped_drainage_json",
          "drainage_divergence" in OPTIONAL_LAYERS and "mapped_rivers_used" in OPTIONAL_LAYERS
          and "mapped_drainage_json" in [k for k, _ in METADATA_KEYS])
    pj = json.loads(params_json("OSM", 1.0, 60.0, 10, (0.8, 0.7, 0.75), 3, "name"))
    check("params JSON: source, threshold, tolerance, overall scores",
          pj["source"] == "OSM" and pj["tolerance_m"] == 60.0 and pj["f1"] == 0.75)


def main(argv=None):
    test_network()
    test_crossings_and_road()
    test_divergence()
    test_contract()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
