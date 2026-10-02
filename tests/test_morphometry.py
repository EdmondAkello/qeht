# -*- coding: utf-8 -*-
"""v0.15 validation: flow-path segments (F2) and basin shape / network
indices (F6). Bare Python + NumPy, no QGIS, no GDAL.

    python -m qeht.tests.test_morphometry
"""

import math
import sys

import numpy as np

from ..core.watershed.morphometry import (lfp_split, contour_perimeter, shape_indices,
                                          channel_network_stats)
from ..core.watershed.delineate import longest_flow_path, delineate_catchment
from ..core.interop.heas_exchange import build_exchange_records
from .test_crossings import valley_dem, route, gt_for

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def test_split_straight_path():
    print("\n1. F2 overland / channel split on a straight path")
    cs = 10.0
    n = 50
    cells = [(0, c) for c in range(n)]                 # 49 links x 10 m = 490 m
    z = np.zeros((1, n)); z[0, :] = 200.0 - 0.05 * np.arange(n) * cs   # 5 % fall
    z[0, 20:] = z[0, 20] - 0.01 * (np.arange(n - 20)) * cs              # 1 % after the head
    acc = np.arange(n, dtype=float).reshape(1, n) * 10.0                # 10 cells per step
    out = lfp_split(cells, z, acc, threshold_cells=200.0, cell_width=cs, cell_height=cs,
                    sheet_cap_m=100.0)
    check("channel head at the first cell with acc >= threshold (cell 20)",
          abs(out["lfp_overland_m"] - 200.0) < 1e-9, f"overland {out['lfp_overland_m']} m")
    check("overland + channel = path length (490 m)",
          abs(out["lfp_overland_m"] + out["lfp_channel_m"] - 490.0) < 1e-9)
    check("overland slope 5 %, channel slope 1 %",
          abs(out["lfp_overland_slope"] - 0.05) < 1e-12 and abs(out["lfp_channel_slope"] - 0.01) < 1e-12,
          f"{out['lfp_overland_slope']:.4f}, {out['lfp_channel_slope']:.4f}")
    check("10-85 slope of a uniform channel = its slope",
          abs(out["lfp_channel_slope_1085"] - 0.01) < 1e-12)
    check("sheet = 100 m cap, shallow = 100 m remainder",
          out["lfp_sheet_m"] == 100.0 and abs(out["lfp_shallow_m"] - 100.0) < 1e-9)
    check("threshold recorded as km2 (200 cells x 100 m2 = 0.02 km2), no_channel 0",
          abs(out["lfp_threshold_km2"] - 0.02) < 1e-12 and out["lfp_no_channel"] == 0)
    hi = lfp_split(cells, z, acc, threshold_cells=300.0, cell_width=cs, cell_height=cs)
    check("raising the threshold moves the head downstream (300 m)",
          abs(hi["lfp_overland_m"] - 300.0) < 1e-9)
    none = lfp_split(cells, z, acc, threshold_cells=1e9, cell_width=cs, cell_height=cs)
    check("no channel cell: whole path overland, channel 0, flag 1",
          abs(none["lfp_overland_m"] - 490.0) < 1e-9 and none["lfp_channel_m"] == 0.0
          and none["lfp_no_channel"] == 1 and none["lfp_channel_slope"] is None)
    diag = [(k, k) for k in range(5)]
    d = lfp_split(diag, np.zeros((5, 5)), np.eye(5) * np.arange(5)[:, None] * 100,
                  threshold_cells=200.0, cell_width=cs, cell_height=cs)
    check("diagonal steps measured as cell diagonals",
          abs(d["lfp_overland_m"] - 2 * cs * math.sqrt(2)) < 1e-9, f"{d['lfp_overland_m']:.3f}")


def test_shape_ratios():
    print("\n2. F6 shape ratios on a disc and a square")
    cs = 1.0
    n = 401
    r, c = np.mgrid[0:n, 0:n]
    R = 150.0
    disc = (r - 200) ** 2 + (c - 200) ** 2 <= R ** 2
    p = contour_perimeter(disc, cs, cs)
    check("disc perimeter within 1 % of 2 pi R", abs(p / (2 * math.pi * R) - 1) < 0.01,
          f"{p:.1f} vs {2 * math.pi * R:.1f}")
    a_km2 = disc.sum() * cs * cs / 1e6
    out = shape_indices(disc, a_km2, 2 * R, cs, cs)
    check("disc circularity ratio ~ 1 (Miller)", abs(out["circularity_ratio"] - 1) < 0.02,
          f"{out['circularity_ratio']:.4f}")
    check("disc with L = diameter: form factor pi/4, elongation 1",
          abs(out["form_factor"] - math.pi / 4) < 0.01 and abs(out["elongation_ratio"] - 1) < 0.005,
          f"Rf {out['form_factor']:.4f}, Re {out['elongation_ratio']:.4f}")
    sq = np.zeros((100, 100), bool); sq[20:80, 20:80] = True       # 60 x 60 cells, 10 m
    ps = contour_perimeter(sq, 10.0, 10.0)
    check("square (600 m side) outline within 2 % of 2400 m (corners rounded)", abs(ps / 2400.0 - 1) < 0.02,
          f"{ps:.1f} m")
    stair = 4 * 600.0
    check("disc estimate far closer than the cell-edge staircase",
          abs(p - 2 * math.pi * R) < 0.2 * abs(8 * R - 2 * math.pi * R))
    sq_out = shape_indices(sq, 0.36, 600.0, 10.0, 10.0)
    check("square: circularity ~ pi/4 (+-0.03)", abs(sq_out["circularity_ratio"] - math.pi / 4) < 0.03,
          f"{sq_out['circularity_ratio']:.4f}")
    ring = np.zeros((50, 50), bool); ring[10:40, 10:40] = True; ring[20:30, 20:30] = False
    check("a hole adds its own outline",
          contour_perimeter(ring, 1, 1) > contour_perimeter(np.pad(np.ones((30, 30), bool), 10), 1, 1))


def test_network_indices():
    print("\n3. F6 drainage density, links and Strahler on a synthetic network")
    rows, cols, cs = 80, 80, 10.0
    gt = gt_for(rows, cs)
    # main valley N -> S along x = 405 and a tributary from the NE joining at y = 400
    dem = valley_dem([[(405.0, 800.0), (405.0, 0.0)], [(795.0, 795.0), (405.0, 400.0)]],
                     rows, cols, cs, slope=0.02, side=0.3)
    v, f, d, a, st = route(dem, cs, threshold=30)
    out_rc = (rows - 1, int(405 // cs))
    mask = delineate_catchment(d, v, [out_rc]) > 0
    length, links, mx = channel_network_stats(mask, d, v, st, cs, cs)
    # independent count: links = sources + confluences
    from ..core.grid import DROW, DCOL
    donors = np.zeros((rows, cols), int)
    rr, cc = np.nonzero(st & mask)
    for r0, c0 in zip(rr, cc):
        k = d[r0, c0]
        if k >= 0:
            r1, c1 = r0 + DROW[k], c0 + DCOL[k]
            if 0 <= r1 < rows and 0 <= c1 < cols and st[r1, c1] and mask[r1, c1]:
                donors[r1, c1] += 1
    nd = donors[rr, cc]
    check("links = sources + confluence cells (independent count)",
          links == int((nd == 0).sum() + (nd >= 2).sum()) and links >= 3, f"{links} links")
    check("channel length >= straight main valley (>= 780 m)", length >= 780.0, f"{length:.0f} m")
    from ..core.flow.accumulation import strahler_order
    so = strahler_order(d, v, st).astype(float)
    _, _, mx = channel_network_stats(mask, d, v, st, cs, cs, strahler=so)
    check("max Strahler >= 2 where two first-order streams meet", mx is not None and mx >= 2,
          f"order {mx}")
    pts = [{"x": (out_rc[1] + 0.5) * cs, "y": gt[3] - (out_rc[0] + 0.5) * cs, "fid": 1,
            "source_id": None}]
    cr, ca, fp, issues, _ = build_exchange_records(d, v, a, dem, gt, pts, snap_radius_cells=0,
                                                   channel_threshold_cells=30)
    c0 = ca[0][1]
    check("package catchment carries F6 fields consistent with the helper",
          abs(c0["drainage_density"] - (length / 1000.0) / c0["area_km2"]) < 1e-9
          and c0["max_strahler"] == mx and c0["perimeter_km"] > 0,
          f"Dd {c0['drainage_density']:.2f} km/km2, Fs {c0['stream_frequency']:.1f}/km2, "
          f"Rc {c0['circularity_ratio']:.3f}")
    f0 = fp[0][1]
    check("package flow path carries F2 fields; overland + channel = length",
          abs(f0["lfp_overland_m"] + f0["lfp_channel_m"] - f0["lfp_length_m"]) < 1e-6
          and f0["lfp_no_channel"] == 0, f"overland {f0['lfp_overland_m']:.0f} m, "
          f"channel {f0['lfp_channel_m']:.0f} m")
    cr2, ca2, fp2, _, _ = build_exchange_records(d, v, a, dem, gt, pts, snap_radius_cells=0)
    check("without a channel threshold: no F2 values, shape ratios still given",
          fp2[0][1].get("lfp_overland_m") is None and ca2[0][1]["circularity_ratio"] is not None
          and ca2[0][1]["drainage_density"] is None)


def main(argv=None):
    test_split_straight_path()
    test_shape_ratios()
    test_network_indices()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
