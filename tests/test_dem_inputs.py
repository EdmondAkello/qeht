# -*- coding: utf-8 -*-
"""1.0.2 input checks: truncated catchments and rasters without a CRS.
NumPy (+ GDAL for the raster check).

    python -m qeht.tests.test_dem_inputs
"""

import os
import sys
import tempfile

import numpy as np

from ..core.conditioning.fill import fill_depressions
from ..core.flow.direction import d8_direction
from ..core.flow.accumulation import flow_accumulation
from ..core.watershed.statistics import catchment_truncation, catchment_characteristics
from ..core.interop.heas_exchange import build_exchange_records

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def test_truncation_mask():
    print("\n1. Truncation flag on masks")
    valid = np.ones((20, 20), dtype=bool)
    m = np.zeros_like(valid)
    m[5:10, 5:10] = True
    check("interior catchment: not truncated", catchment_truncation(m, valid, outlet=(9, 7)) == 0)
    m2 = np.zeros_like(valid)
    m2[0:8, 4:9] = True                                   # divide runs along the top edge
    n = catchment_truncation(m2, valid, outlet=(7, 6))
    check("catchment reaching the top edge away from its outlet: truncated", n == 5, f"{n} cells")
    m3 = np.zeros_like(valid)
    m3[15:20, 8:10] = True                                # outlet at the bottom edge, nothing else on it
    m3[19, 8:10] = False
    m3[19, 9] = True
    check("outlet at the edge only: not truncated", catchment_truncation(m3, valid, outlet=(19, 9)) == 0)
    v4 = valid.copy()
    v4[3:6, 3:6] = False                                  # NoData hole beside the catchment
    m4 = np.zeros_like(valid)
    m4[6:12, 2:8] = True
    check("catchment beside a NoData area: truncated", catchment_truncation(m4, v4, outlet=(11, 5)) > 0)
    z = np.tile(np.arange(20, 0, -1, dtype=float)[:, None], (1, 20))
    ch = catchment_characteristics(m2, z, valid, 30.0, 30.0, outlet=(7, 6))
    check("catchment_characteristics carries truncated / truncated_cells",
          ch.get("truncated") == 1 and ch.get("truncated_cells") == 5, str({k: ch.get(k) for k in ("truncated", "truncated_cells")}))
    ch = catchment_characteristics(m, z, valid, 30.0, 30.0, outlet=(9, 7))
    check("interior catchment: truncated = 0", ch.get("truncated") == 0)


def test_truncation_in_exchange():
    print("\n2. Truncation in the exchange records (issue text, field)")
    rows, cols, cs = 60, 40, 30.0
    r, c = np.mgrid[0:rows, 0:cols]
    dem = 200.0 - 2.0 * r + 0.5 * np.abs(c - cols / 2)    # valley draining south, divides at the edges
    valid = np.ones_like(dem, dtype=bool)
    gt = (500000.0, cs, 0.0, 9900000.0, 0.0, -cs)
    filled, _, _ = fill_depressions(dem, valid, cell_width=cs, cell_height=cs)
    direction, _ = d8_direction(filled, valid, cs, cs)
    acc, _ = flow_accumulation(direction, valid)
    xy = lambda row, col: (gt[0] + (col + 0.5) * cs, gt[3] - (row + 0.5) * cs)
    pts = [dict(zip(("x", "y"), xy(50, cols // 2)), fid=1, source_id=None)]
    cr, catch, fp, issues, _ = build_exchange_records(direction, valid, acc, dem, gt, pts,
                                                      snap_radius_cells=2, id_prefix="T")
    flag = catch[0][1].get("truncated") if catch else None
    check("catchment whose divide is the DEM edge: truncated = 1 and an issue says area is a lower bound",
          flag == 1 and any("lower bound" in i for i in issues), f"truncated {flag}; {issues[:1]}")


def test_raster_without_crs():
    print("\n3. A factor raster without a CRS is refused")
    try:
        from osgeo import gdal, osr
    except ImportError:
        print("  [SKIP] GDAL not available")
        return
    from ..core.raster import RasterInfo, warp_to_grid
    gdal.UseExceptions()
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32737)
    gt = (500000.0, 30.0, 0.0, 9900000.0, 0.0, -30.0)
    info = RasterInfo(gt, srs.ExportToWkt(), 10, 10)
    with tempfile.TemporaryDirectory() as tmp:
        for name, wkt in (("with_crs.tif", srs.ExportToWkt()), ("no_crs.tif", None)):
            ds = gdal.GetDriverByName("GTiff").Create(os.path.join(tmp, name), 10, 10, 1, gdal.GDT_Float32)
            ds.SetGeoTransform(gt)
            if wkt:
                ds.SetProjection(wkt)
            ds.GetRasterBand(1).WriteArray(np.full((10, 10), 0.03, dtype=np.float32))
            ds = None
        a, _ = warp_to_grid(os.path.join(tmp, "with_crs.tif"), info)
        check("raster with a CRS: read onto the DEM grid", np.allclose(a, 0.03))
        try:
            warp_to_grid(os.path.join(tmp, "no_crs.tif"), info)
            ok, msg = False, "no error"
        except ValueError as e:
            ok, msg = "no coordinate reference system" in str(e), str(e)[:80]
        check("raster without a CRS: refused with a clear message", ok, msg)


def main():
    test_truncation_mask()
    test_truncation_in_exchange()
    test_raster_without_crs()
    print(f"\n{'ALL PASS' if not FAILURES else f'{len(FAILURES)} FAILURE(S)'}")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
