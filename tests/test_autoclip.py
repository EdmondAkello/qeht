# -*- coding: utf-8 -*-
"""v0.25 validation: corridor auto-clip for large DEMs (F13).
NumPy (the window writer uses GDAL). A synthetic 600 x 500 grid of 10 m:
valleys draining south across an E-W road, a divide at row 50 (the north
strip drains north), an eastern block that drains east and the ground
south of the road draining away from it.

    python -m qeht.tests.test_autoclip
"""

import json
import math
import os
import sys
import tempfile

import numpy as np

from ..core.conditioning.autoclip import (contributing_window, coarse_factor, block_min, drains_to,
                                          dilate_cells, touches_edge, memory_note, write_clip, TAG)
from ..core.conditioning.fill import fill_depressions
from ..core.flow.direction import d8_direction
from ..core.flow.accumulation import flow_accumulation
from ..core.network.alignment import Alignment

FAILURES = []
ROWS, COLS, CS = 600, 500, 10.0
GT = (300000.0, CS, 0.0, 9800000.0, 0.0, -CS)
ROAD_ROW = 250


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def terrain():
    r, c = np.mgrid[0:ROWS, 0:COLS].astype(float)
    ridges = 8.0 * np.abs(np.sin(c * math.pi / 40.0))       # valleys every 40 columns
    z = 500.0 - 0.5 * r + ridges + 0.001 * c                 # drains south across the road
    z = np.where(r < 50, 475.0 + 0.5 * r + ridges, z)        # north strip drains north
    return np.where(c >= 220, 700.0 - 2.0 * (c - 220) - 0.1 * r, z)   # east block drains east


def route(z, valid=None):
    valid = np.ones(z.shape, bool) if valid is None else valid
    f, _, _ = fill_depressions(np.where(valid, z, 0.0), valid, cell_width=CS, cell_height=CS)
    d, _ = d8_direction(f, valid, CS, CS)
    a, _ = flow_accumulation(d, valid)
    return d, a, valid


def road():
    y = GT[3] - (ROAD_ROW + 0.5) * CS
    return Alignment([[(GT[0] + 55.0, y), (GT[0] + 2105.0, y)]])


def test_clip():
    print("\n1. Clip contains every catchment of the road crossings")
    z = terrain()
    d, a, valid = route(z)
    xs = [c for c in range(5, 211) if a[ROAD_ROW, c] >= 200]
    check("road crossings found on the full grid", len(xs) >= 4, str(xs))
    check("coarse factor for 4 Mcells (default) on 0.3 Mcells is 1; on 0.012 Mcells, 5",
          coarse_factor(ROWS, COLS) == 1 and coarse_factor(ROWS, COLS, 0.012) == 5)
    w = contributing_window(z, valid, GT, road(), max_mcells=0.012)
    r0, r1, c0, c1 = w["window"]
    check("window leaves out the east block and the ground beyond the margin south of the road",
          c1 <= 340 and r1 < ROWS and r1 >= ROAD_ROW + 100 and r0 == 0,
          f"rows {r0}-{r1}, cols {c0}-{c1}, k {w['k']}")
    zc = z[r0:r1, c0:c1]
    dc, ac, vc = route(zc)
    same = all(a[ROAD_ROW, c] == ac[ROAD_ROW - r0, c - c0] for c in xs)
    check("areas at every crossing equal the unclipped run exactly", same,
          ", ".join(f"{a[ROAD_ROW, c]:.0f}/{ac[ROAD_ROW - r0, c - c0]:.0f}" for c in xs))
    edges = []
    for c in xs:
        m = np.zeros(zc.shape, bool); m[ROAD_ROW - r0, c - c0] = True
        edges.append(touches_edge(drains_to(m, dc, vc), vc))
    check("no catchment touches the clip edge", not any(edges), str(edges))
    mem = memory_note(w["full_cells"], w["clip_cells"])
    check("memory estimate lower than unclipped", mem["clip_peak_gb"] < mem["full_peak_gb"]
          and w["clip_cells"] < w["full_cells"], str(mem))
    w0 = contributing_window(z, valid, GT, road(), road_buffer_m=0.0, margin_m=0.0, max_mcells=0.012,
                             pad_cells=0)
    r0, r1, c0, c1 = w0["window"]
    check("margin and padding 0: the tightest window still holds the catchments (starts at the "
          "divide, row 50, or above)", r0 <= 50 and r1 >= ROAD_ROW, str(w0["window"]))
    r0 += 20                                              # cut the catchments on purpose
    zc0 = z[r0:r1, c0:c1]
    dc0, ac0, vc0 = route(zc0)
    e0 = []
    for c in xs:
        if c0 <= c < c1:
            m = np.zeros(zc0.shape, bool); m[ROAD_ROW - r0, c - c0] = True
            e0.append(touches_edge(drains_to(m, dc0, vc0), vc0))
    check("a window that cuts the catchments (top 20 rows removed): clip_edge = 1", all(e0), str(e0))
    wm = contributing_window(z, valid, GT, road(), max_mcells=0.012, mode="mask")
    check("mask mode: fewer valid cells than the bounding box", wm["clip_cells_valid"] < wm["clip_cells"]
          and wm["mask"].shape == (wm["window"][1] - wm["window"][0], wm["window"][3] - wm["window"][2]))


def test_parts():
    print("\n2. Pieces")
    z = np.array([[5, 4, 9, 9], [6, 1, 9, 9], [7, 7, 2, 8.0]])
    v = np.ones(z.shape, bool); v[1, 1] = False
    bm = block_min(z, v, 2)
    check("block minimum keeps valleys (NoData skipped)", bm.tolist() == [[4.0, 9.0], [7.0, 2.0]])
    m = np.zeros((7, 7), bool); m[3, 3] = True
    check("square dilation by 2 cells: 5 x 5", dilate_cells(m, 2).sum() == 25)
    d = np.full((1, 6), 0, np.int64); d[0, -1] = -1               # all flow east
    mk = np.zeros((1, 6), bool); mk[0, 3] = True
    check("drains_to: cells upstream of the mark", drains_to(mk, d, np.ones((1, 6), bool)).tolist()
          == [[True, True, True, True, False, False]])


def test_write(tmp):
    print("\n3. Window writer and tag")
    from osgeo import gdal, osr
    gdal.UseExceptions()
    src = os.path.join(tmp, "big.tif")
    ds = gdal.GetDriverByName("GTiff").Create(src, COLS, ROWS, 1, gdal.GDT_Float32)
    ds.SetGeoTransform(GT)
    s = osr.SpatialReference(); s.ImportFromEPSG(32737)
    ds.SetProjection(s.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(terrain().astype(np.float32))
    ds = None
    from ..core.raster import read_dem, raster_tags
    z, valid, info = read_dem(src)
    w = contributing_window(z, valid, info.geotransform, road(), max_mcells=0.012, mode="mask")
    out = write_clip(src, os.path.join(tmp, "clip.tif"), info, w["window"], mask=w["mask"],
                     meta={"k": w["k"]})
    zc, vc, ic = read_dem(out)
    r0, r1, c0, c1 = w["window"]
    check("clip: window size, origin and values; NoData outside the mask",
          zc.shape == (r1 - r0, c1 - c0) and math.isclose(ic.geotransform[0], GT[0] + c0 * CS)
          and np.allclose(zc[vc], terrain()[r0:r1, c0:c1][w["mask"]].astype(np.float32))
          and vc.sum() == w["mask"].sum())
    tag = json.loads(raster_tags(out)[TAG])
    check("tag QEHT_AUTOCLIP records the window", tag["window"] == [r0, r1, c0, c1] and tag["k"] == w["k"])


def main(argv=None):
    tmp = tempfile.mkdtemp(prefix="qeht_clip_")
    test_clip()
    test_parts()
    test_write(tmp)
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
