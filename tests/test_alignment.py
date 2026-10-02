# -*- coding: utf-8 -*-
"""v0.15 validation: alignment ground profile (A1) and later alignment work.
Bare Python + NumPy, no QGIS, no GDAL.

    python -m qeht.tests.test_alignment
"""

import math
import sys

import numpy as np

from ..core.network.alignment import Alignment
from ..core.network.profile import (alignment_profile, profile_with_crossings, bilinear,
                                    longitudinal_slope_pct, stations)
from .test_crossings import valley_dem, route, gt_for

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def test_plane_profile():
    print("\n1. Alignment profile on an analytic plane")
    rows, cols, cs = 50, 80, 10.0
    gt = gt_for(rows, cs)
    r, c = np.mgrid[0:rows, 0:cols].astype(float)
    x = (c + 0.5) * cs
    y = rows * cs - (r + 0.5) * cs
    z = 100.0 + 0.02 * x + 0.005 * y          # plane: 2 % in x, 0.5 % in y
    al = Alignment([[(55.0, 250.0), (705.0, 250.0)]], start_chainage=1000.0)
    prof = alignment_profile(al, z, gt, step=10.0)
    ch = np.array([p["chainage_m"] for p in prof])
    zz = np.array([p["z_dem_m"] for p in prof])
    xs = np.array([p["x"] for p in prof])
    exact = 100.0 + 0.02 * xs + 0.005 * 250.0
    check("stations every 10 m from the start chainage, end included",
          ch[0] == 1000.0 and abs(ch[1] - 1010.0) < 1e-9 and abs(ch[-1] - 1650.0) < 1e-9,
          f"{len(ch)} stations, {ch[0]:g} -> {ch[-1]:g}")
    check("bilinear z is exact on a plane", np.nanmax(np.abs(zz - exact)) < 1e-9,
          f"max error {np.nanmax(np.abs(zz - exact)):.2e} m")
    sl = np.array([p["slope_long_pct"] for p in prof])
    check("longitudinal slope = 2.000 % along x", np.nanmax(np.abs(sl - 2.0)) < 1e-9,
          f"{np.nanmin(sl):.6f} .. {np.nanmax(sl):.6f}")
    rev = alignment_profile(Alignment([[(55.0, 250.0), (705.0, 250.0)]], reverse=True),
                            z, gt, step=10.0)
    check("reversed alignment: chainage from the other end, slope -2 %",
          abs(rev[0]["x"] - 705.0) < 1e-9 and abs(rev[0]["slope_long_pct"] + 2.0) < 1e-9)
    diag = alignment_profile(Alignment([[(55.0, 55.0), (455.0, 355.0)]]), z, gt, step=10.0)
    # direction (0.8, 0.6): slope = 0.02*0.8 + 0.005*0.6 = 0.019
    dsl = np.array([p["slope_long_pct"] for p in diag])
    check("oblique alignment: slope = grad . tangent (1.9 %)",
          np.nanmax(np.abs(dsl - 1.9)) < 1e-9, f"{np.nanmean(dsl):.6f} %")
    check("no filled DEM / accumulation -> empty fields, stream 0",
          all(p["z_fill_m"] is None and p["acc_km2"] is None and p["stream"] == 0 for p in prof))
    s2 = stations(Alignment([[(0, 0), (25, 0)]]), 10.0)
    check("last station at the alignment end", list(s2) == [0.0, 10.0, 20.0, 25.0], str(list(s2)))
    check("slope helper one-sided at ends",
          np.allclose(longitudinal_slope_pct([0, 10, 20], [0, 1, 3]), [10, 15, 20]))
    g = np.arange(12, dtype=float).reshape(3, 4)
    g[1, 1] = np.nan
    v = bilinear(g, (0, 1, 0, 3, 0, -1), np.array([1.0, 3.5, 9.0]), np.array([1.5, 2.5, 1.0]))
    check("bilinear NaN next to NoData and outside the raster",
          math.isnan(v[0]) and abs(v[1] - 3.0) < 1e-12 and math.isnan(v[2]), str(v))


def test_valley_crossing():
    print("\n2. V-valley crossing: stream flag, area, ponding")
    rows, cols, cs = 60, 60, 10.0
    gt = gt_for(rows, cs)
    # valley along x = 305 flowing south (ends at the bottom edge)
    dem = valley_dem([[(305.0, 600.0), (305.0, 0.0)]], rows, cols, cs, slope=0.02, side=0.3)
    # road across it at y = 205 with an embankment-like bump: a pit upstream
    v, f, d, a, st = route(dem, cs, threshold=30)
    al = Alignment([[(5.0, 205.0), (595.0, 205.0)]])
    prof, cr = profile_with_crossings(al, dem, gt, direction=d, valid=v, accumulation=a,
                                      filled=f, stream_threshold_cells=30, step=10.0)
    flagged = [p for p in prof if p["stream"] == 1]
    check("exactly one stream crossing flagged", len(flagged) == 1 and len(cr) == 1,
          f"{len(flagged)} flagged, {len(cr)} crossings")
    if flagged:
        p = flagged[0]
        check("flagged station at the valley floor (x = 305 +- 5 m)", abs(p["x"] - 305.0) <= 5.0,
              f"x = {p['x']:.1f}")
        # cells upstream of y=205 along the valley drain to it; compare with the grid
        rr = int((gt[3] - 205.0) // cs)
        col = int(305.0 // cs)
        ref = (a[rr - 1:rr + 2, col - 1:col + 2].max() + 1) * cs * cs / 1e6
        check("acc_km2 = largest (acc+1) x cell area within 1 cell",
              abs(p["acc_km2"] - max(ref, cr[0]["acc_km2"])) < 1e-12,
              f"{p['acc_km2']:.4f} km2")
        check("valley floor is the profile minimum", p["z_dem_m"] <= min(
            q["z_dem_m"] for q in prof) + 1e-6)
    pond = [p["pond_depth_m"] for p in prof]
    check("pond depth >= 0 everywhere and defined with a filled DEM",
          all(x is not None and x >= 0 for x in pond))
    # a pit on the road line gets a positive pond depth
    dem2 = dem.copy()
    rr = int((gt[3] - 205.0) // cs)
    dem2[rr, 10] -= 10.0
    v2, f2, d2, a2, _ = route(dem2, cs, threshold=30)
    prof2, _ = profile_with_crossings(al, dem2, gt, direction=d2, valid=v2, accumulation=a2,
                                      filled=f2, stream_threshold_cells=30, step=10.0)
    mx = max(p["pond_depth_m"] for p in prof2)
    check("a 10 m pit on the line shows as ponding depth", mx > 5.0, f"max {mx:.2f} m")


def main(argv=None):
    test_plane_profile()
    test_valley_crossing()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
