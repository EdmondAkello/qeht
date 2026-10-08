# -*- coding: utf-8 -*-
"""v0.19 validation: approach-channel cross-sections at crossings (F5).
Bare Python + NumPy, no QGIS, no GDAL. Synthetic channels with analytic
sections; the grids are 1 m where a 14 m channel must be resolved exactly.

    python -m qeht.tests.test_section
"""

import json
import math
import sys

import numpy as np

from ..core.watershed.section import (section_at_outlet, section_block, section_location,
                                      params_json, XS_FIELDS)
from ..core.flow.direction import d8_direction

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def trapezoid(u, bed=6.0, hv_l=2.0, h_l=2.0, hv_r=2.0, h_r=2.0, fp=0.05):
    """Cross profile above the bed at signed offset u (+ = side 'r').
    Flat bed of width `bed`, banks H:V hv over height h, then a 1:20 floodplain."""
    u = np.asarray(u, float)
    hv = np.where(u >= 0, hv_r, hv_l)
    h = np.where(u >= 0, h_r, h_l)
    a = np.abs(u) - bed / 2.0
    return np.where(a <= 0, 0.0, np.where(a <= hv * h, a / hv, h + (a - hv * h) * fp))


def grid_ns(rows, cols, cs, profile, c0, long_slope=0.001):
    """Channel along column c0 draining south; profile(u) with u = + towards the WEST
    (the right bank looking downstream), in metres."""
    r, c = np.mgrid[0:rows, 0:cols].astype(float)
    u = (c0 - c) * cs
    return 100.0 + long_slope * (rows - 1 - r) * cs + profile(u)


def gt_for(cs, rows):
    return (500000.0, cs, 0.0, 9900000.0 + rows * cs, 0.0, -cs)


def run(dem, cs, r, c, **kw):
    valid = np.isfinite(dem)
    d, _ = d8_direction(np.where(valid, dem, 0.0), valid, cs, cs)
    gt = gt_for(cs, dem.shape[0])
    return section_at_outlet(d, valid, np.where(valid, dem, np.nan), gt, r, c, **kw), d, gt


def test_trapezoid():
    print("\n1. Analytic trapezoidal channel, N-S, 1 m grid")
    rows, cols, cs, c0 = 90, 81, 1.0, 40
    dem = grid_ns(rows, cols, cs, trapezoid, c0)
    (f, t), _, _ = run(dem, cs, 20, c0, half_m=30.0)
    bed = 100.0 + 0.001 * (rows - 1 - 50)
    check("section 30 m down the receivers, at row 50", f["xs_dist_m"] == 30.0)
    check("bed level exact", math.isclose(f["xs_bed_m"], bed, abs_tol=1e-9),
          f"{f['xs_bed_m']:.6f} vs {bed:.6f}")
    check("bank-full width 6 + 2 x 2 x 2 = 14 m (within one 0.5 m step)",
          abs(f["xs_bankfull_w_m"] - 14.0) <= 0.5, f"{f['xs_bankfull_w_m']:.3f}")
    check("bank-full depth 2 m", abs(f["xs_bankfull_d_m"] - 2.0) <= 0.05, f"{f['xs_bankfull_d_m']:.3f}")
    check("side slopes 2.0 +/- 0.1 (H:V)",
          abs(f["xs_side_slope_l"] - 2.0) <= 0.1 and abs(f["xs_side_slope_r"] - 2.0) <= 0.1,
          f"{f['xs_side_slope_l']:.3f} / {f['xs_side_slope_r']:.3f}")
    check("fixed stages: 6 + 2 x 2 x dz below bank top (8 / 10 m), 14 + 2 x 20 x 0 at 2 m",
          math.isclose(f["xs_w_0p5_m"], 8.0) and math.isclose(f["xs_w_1p0_m"], 10.0)
          and math.isclose(f["xs_w_2p0_m"], 14.0),
          f"{f['xs_w_0p5_m']}, {f['xs_w_1p0_m']}, {f['xs_w_2p0_m']}")
    check("quality high, no note", f["xs_quality"] == "high" and f["xs_note"] is None, str(f["xs_note"]))
    st = json.loads(f["xs_station_elev_json"])
    check("station list: -30 .. +30 m every 0.5 m, rounded to 0.01 m",
          len(st) == 121 and st[0][0] == -30.0 and st[-1][0] == 30.0
          and all(round(z, 2) == z for _, z in st))
    check("right side (+ offset) is the WEST bank looking downstream (south)",
          t["line"][1][0] < t["line"][0][0] and math.isclose(t["line"][0][1], t["line"][1][1]))
    check("keys = XS_FIELDS", set(f) == set(XS_FIELDS))


def test_diagonal():
    print("\n2. Diagonal channel (NE-SW)")
    rows = cols = 130
    cs, r0, c0 = 1.0, 65, 65
    r, c = np.mgrid[0:rows, 0:cols].astype(float)
    u = ((c - c0) + (r - r0)) / math.sqrt(2.0)          # across (+ = SE side)
    t_ = ((c - c0) - (r - r0)) / math.sqrt(2.0)         # along, + towards NE
    dem = 100.0 + 0.001 * t_ + trapezoid(u)
    (f, t), d, _ = run(dem, cs, 35, 95, half_m=30.0)
    fx, fy = t["flow"]
    ax = np.array([-1.0, -1.0]) / math.sqrt(2.0)       # SW in map coordinates
    (lx0, ly0), (lx1, ly1) = t["line"]
    ang = math.degrees(math.acos(abs((lx1 - lx0) * ax[0] + (ly1 - ly0) * ax[1])
                                 / math.hypot(lx1 - lx0, ly1 - ly0)))
    check("transect perpendicular to the channel (within 5 degrees of 90)", abs(ang - 90.0) <= 5.0,
          f"{ang:.2f} deg; flow ({fx:.3f}, {fy:.3f})")
    check("section 22 diagonal steps down: 30 m requested -> 31.1 m reached",
          math.isclose(f["xs_dist_m"], 22 * math.sqrt(2.0)), f"{f['xs_dist_m']:.2f}")
    rows2, cols2 = 90, 81
    (g, _), _, _ = run(grid_ns(rows2, cols2, cs, trapezoid, 40), cs, 20, 40, half_m=30.0)
    # bed + 2 m is exactly the bank-top corner, which bilinear interpolation rounds off
    # across a diagonal: that one stage is held to one cell, the others to one step
    same = all(abs(f[k] - g[k]) <= 0.5 for k in ("xs_bankfull_w_m", "xs_w_0p5_m", "xs_w_1p0_m")) \
        and abs(f["xs_w_2p0_m"] - g["xs_w_2p0_m"]) <= 1.0
    check("widths equal the N-S case within one step (bank-top corner within one cell)",
          same, ", ".join(f"{k[3:]} {f[k]:.2f}/{g[k]:.2f}" for k in
                          ("xs_bankfull_w_m", "xs_w_0p5_m", "xs_w_1p0_m", "xs_w_2p0_m")))
    check("depth and side slopes as N-S", abs(f["xs_bankfull_d_m"] - 2.0) <= 0.1
          and abs(f["xs_side_slope_l"] - 2.0) <= 0.2 and abs(f["xs_side_slope_r"] - 2.0) <= 0.2,
          f"d {f['xs_bankfull_d_m']:.2f}, {f['xs_side_slope_l']:.2f} / {f['xs_side_slope_r']:.2f}")


def test_subgrid():
    print("\n3. Sub-grid channel (single-cell notch, 10 m grid)")
    rows, cols, cs, c0 = 40, 41, 10.0, 20

    def notch(u):
        return 0.01 * np.abs(u) + np.where(np.abs(u) < 1.0, 0.0, 2.0)
    dem = grid_ns(rows, cols, cs, notch, c0, long_slope=0.002)
    (f, _), _, _ = run(dem, cs, 5, c0)
    check("quality low with a sub-grid note",
          f["xs_quality"] == "low" and "under 3 cells" in (f["xs_note"] or ""),
          f"{f['xs_quality']}: {f['xs_note']} (w {f['xs_bankfull_w_m']})")
    check("bank-full width is the bilinear V between cell centres (20 m)",
          math.isclose(f["xs_bankfull_w_m"], 20.0, abs_tol=1e-6), f"{f['xs_bankfull_w_m']}")
    vp = lambda u: 0.1 * np.abs(u)
    (g, _), _, _ = run(grid_ns(90, 81, 1.0, vp, 40), 1.0, 20, 40, half_m=30.0)
    check("no bank on either side (V valley): low, depth and width NULL",
          g["xs_quality"] == "low" and g["xs_bankfull_w_m"] is None
          and "no bank" in g["xs_note"])


def test_location():
    print("\n4. Section location and truncation")
    rows, cols, cs, c0 = 30, 41, 10.0, 20
    big = lambda u: trapezoid(u, bed=60.0, h_l=20.0, h_r=20.0)
    dem = grid_ns(rows, cols, cs, big, c0)
    valid = np.ones(dem.shape, bool)
    d, _ = d8_direction(dem, valid, cs, cs)
    path, reached, cut = section_location(d, valid, 5, c0, 30.0, cs, cs)
    check("30 m on a 10 m grid: 3 receiver steps, row 8", reached == 30.0 and path[-1][:2] == (8, c0)
          and not cut)
    (f, _), _, _ = run(dem, cs, 27, c0, half_m=150.0)
    check("walk cut by the grid edge after 20 m: noted, quality lowered",
          f["xs_dist_m"] == 20.0 and "walk stopped after 20 m" in f["xs_note"]
          and f["xs_quality"] == "medium", f"{f['xs_quality']}: {f['xs_note']}")
    dem2 = dem.copy()
    dem2[:, c0 - 6:c0 - 4] = np.nan                       # NoData strip on the right (west) bank
    (g, _), _, _ = run(dem2, cs, 5, c0, half_m=150.0)
    check("NoData across the right bank: medium, noted",
          g["xs_quality"] == "medium" and "NoData" in g["xs_note"], f"{g['xs_quality']}: {g['xs_note']}")


def test_asymmetric():
    print("\n5. Asymmetric banks")
    prof = lambda u: trapezoid(u, hv_l=2.0, h_l=2.0, hv_r=1.0, h_r=3.0)
    (f, _), _, _ = run(grid_ns(90, 81, 1.0, prof, 40), 1.0, 20, 40, half_m=30.0)
    check("depth from the lower (left, 2 m) bank", abs(f["xs_bankfull_d_m"] - 2.0) <= 0.05,
          f"{f['xs_bankfull_d_m']:.3f}")
    check("width at 2 m: 6 + 4 (left 1:2) + 2 (right 1:1) = 12 m",
          abs(f["xs_bankfull_w_m"] - 12.0) <= 0.5, f"{f['xs_bankfull_w_m']:.3f}")
    check("side slopes differ: left 2.0, right 1.0",
          abs(f["xs_side_slope_l"] - 2.0) <= 0.1 and abs(f["xs_side_slope_r"] - 1.0) <= 0.1,
          f"{f['xs_side_slope_l']:.3f} / {f['xs_side_slope_r']:.3f}")
    check("quality high", f["xs_quality"] == "high", str(f["xs_note"]))


def test_v_widths():
    print("\n6. Fixed-stage widths on a V channel")
    for hv in (5.0, 10.0):
        vp = lambda u, hv=hv: np.abs(u) / hv
        (f, _), _, _ = run(grid_ns(90, 101, 1.0, vp, 50), 1.0, 20, 50, half_m=45.0)
        exp = [2 * dz * hv for dz in (0.5, 1.0, 2.0)]
        got = [f["xs_w_0p5_m"], f["xs_w_1p0_m"], f["xs_w_2p0_m"]]
        check(f"H:V {hv:g}: widths 2 dz H = {exp}", all(math.isclose(a, b, abs_tol=1e-6)
                                                       for a, b in zip(got, exp)), str(got))
        check(f"H:V {hv:g}: side slopes from the bed to bed + 1 m",
              math.isclose(f["xs_side_slope_l"], hv, rel_tol=1e-6)
              and math.isclose(f["xs_side_slope_r"], hv, rel_tol=1e-6))


def test_contract():
    print("\n7. Contract")
    from ..core.interop.field_dictionary import field_names, METADATA_KEYS, OPTIONAL_LAYERS
    check("all fields in the crossings contract", all(f in field_names("crossings") for f in XS_FIELDS))
    check("xs_params_json in the metadata keys", "xs_params_json" in [k for k, _ in METADATA_KEYS])
    check("xs_transects in OPTIONAL_LAYERS as LINESTRING",
          OPTIONAL_LAYERS.get("xs_transects", (None,))[0] == "LINESTRING"
          and {"outlet_uid", "xs_dist_m", "xs_station_elev_json"}
          <= {f[0] for f in OPTIONAL_LAYERS["xs_transects"][1]})
    pj = json.loads(params_json(30.0, 150.0, None, 0.05))
    check("params JSON: dist, half width, step, bank slope, stages",
          pj["dist_m"] == 30.0 and pj["half_m"] == 150.0 and pj["bank_slope"] == 0.05
          and pj["dz_m"] == [0.5, 1.0, 2.0] and "Indicative" in pj["note"])
    rows, cols, cs, c0 = 90, 81, 1.0, 40
    dem = grid_ns(rows, cols, cs, trapezoid, c0)
    valid = np.ones(dem.shape, bool)
    d, _ = d8_direction(dem, valid, cs, cs)
    gt = gt_for(cs, rows)
    cr = [{"outlet_uid": "X001", "outlet_x": gt[0] + c0 + 0.5, "outlet_y": gt[3] - 20.5},
          {"outlet_uid": "X002", "outlet_x": gt[0] - 500.0, "outlet_y": gt[3] - 20.5}]
    b, lines = section_block(cr, d, valid, dem, gt, half_m=30.0)
    check("block: one transect per crossing inside the grid, outside -> low with a note",
          len(lines) == 1 and lines[0][1]["outlet_uid"] == "X001" and b[0]["xs_quality"] == "high"
          and b[1]["xs_quality"] == "low" and "outside" in b[1]["xs_note"])


def main(argv=None):
    test_trapezoid()
    test_diagonal()
    test_subgrid()
    test_location()
    test_asymmetric()
    test_v_widths()
    test_contract()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
