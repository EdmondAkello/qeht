# -*- coding: utf-8 -*-
"""1.0.2 DEM audit: vertical storage step, average drop per pixel and the
sloped-tie indicator of nearest-neighbour reprojection. NumPy + GDAL.

    python -m qeht.tests.test_dem_precision
"""

import sys

import numpy as np

from ..core.raster import (vertical_quantum, mean_steepest_drop, precision_stats,
                           audit_precision, audit_resampling, resampling_stats)

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def smooth_surface(n=300, relief=40.0, base=1500.0, seed=3):
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:n, 0:n] / n
    z = base + relief * np.sin(3 * xx) * np.cos(2 * yy) + rng.normal(0, 0.05, (n, n))
    return z, np.ones_like(z, dtype=bool)


def test_quantum():
    print("\n1. Vertical storage step")
    z, v = smooth_surface()
    check("continuous float surface: no step", vertical_quantum(z, v) is None)
    check("rounded to whole metres: 1 m", vertical_quantum(np.round(z), v) == 1.0)
    check("whole metres stored as float32: 1 m",
          vertical_quantum(np.round(z).astype(np.float32), v) == 1.0)
    check("rounded to 0.25 m: 0.25 m", vertical_quantum(np.round(z / 0.25) * 0.25, v) == 0.25)
    cm = (np.round(z * 100) / 100).astype(np.float32)          # FABDEM-like 1 cm steps near 1500 m
    check("1 cm steps stored as float32 at ~1500 m: 0.01 m", vertical_quantum(cm, v) == 0.01,
          f"got {vertical_quantum(cm, v)}")
    check("1 cm steps plus 1 mm noise: no step",
          vertical_quantum(cm + np.random.default_rng(1).normal(0, 0.001, z.shape), v) is None)
    vv = v.copy()
    vv[:, :50] = False
    zz = np.round(z)
    zz[:, :50] = -9999.5                                          # invalid cells are ignored
    check("invalid cells are ignored", vertical_quantum(zz, vv) == 1.0)


def test_drop():
    print("\n2. Mean steepest drop per pixel")
    s, cols = 0.37, 40
    z = np.tile(-s * np.arange(cols, dtype=float), (20, 1))   # plane falling to the east, 1 unit per cell
    v = np.ones_like(z, dtype=bool)
    expect = s * (cols - 1) / cols                            # last column has no lower neighbour
    got = mean_steepest_drop(z, v)
    check("plane: mean drop = slope x (cols-1)/cols", abs(got - expect) < 1e-9, f"{got:.4f} vs {expect:.4f}")
    check("flat grid: mean drop 0", mean_steepest_drop(np.full((10, 10), 5.0), np.ones((10, 10), bool)) == 0.0)
    st = precision_stats(np.round(z * 10) / 10, v)
    check("precision_stats: step, drop and ratio", st["quantum"] == 0.1 and st["ratio"] is not None
          and abs(st["ratio"] - st["mean_drop"] / 0.1) < 1e-9, str(st))


def test_audit_precision():
    print("\n3. Precision audit (Gyasi-Agyei ratio)")
    low, v = smooth_surface(relief=2.0)                       # ~0.02-0.04 m drop per cell
    mid, _ = smooth_surface(relief=300.0)                     # a few metres per cell
    high, _ = smooth_surface(relief=5000.0)
    w = audit_precision(np.round(low), v)
    check("low relief in whole metres: below ratio 1, warned",
          len(w) == 1 and "below the adequacy ratio of 1" in w[0] and "Gyasi-Agyei" in w[0])
    w = audit_precision(np.round(mid), v)
    r = precision_stats(np.round(mid), v)["ratio"]
    check("moderate relief in whole metres: ratio 1-10, caution",
          len(w) == 1 and "necessary rather than sufficient" in w[0] and 1 <= r < 10, f"ratio {r:.2f}")
    r = precision_stats(np.round(high), v)["ratio"]
    check("high relief in whole metres: no warning", not audit_precision(np.round(high), v), f"ratio {r:.1f}")
    check("continuous DEM: no warning", not audit_precision(low, v))
    w = audit_precision((np.round(low * 100) / 100).astype(np.float32), v)
    check("1 cm steps reported in centimetres", (not w) or "1 cm" in w[0])


def _warp(z, resample, rotate_deg=0.0):
    """Warp a geographic float grid (1 arcsec near 1 S, 36.5 E) to UTM 37S at 30 m."""
    from osgeo import gdal
    gdal.UseExceptions()
    n = z.shape[0]
    src = gdal.GetDriverByName("MEM").Create("", n, n, 1, gdal.GDT_Float32)
    src.SetGeoTransform((36.5, 1 / 3600, 0.0, -1.0, 0.0, -1 / 3600))
    src.SetProjection('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],'
                      'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4326"]]')
    src.GetRasterBand(1).WriteArray(z.astype(np.float32))
    src.GetRasterBand(1).SetNoDataValue(-9999)
    out = gdal.Warp("", src, format="MEM", dstSRS="EPSG:32737", xRes=30, yRes=30,
                    targetAlignedPixels=True, resampleAlg=resample, srcNodata=-9999, dstNodata=-9999)
    a = out.GetRasterBand(1).ReadAsArray().astype(np.float64)
    return a, a != -9999


def test_sloped_ties():
    print("\n4. Nearest-neighbour reprojection onto a rotated grid (sloped ties)")
    try:
        from osgeo import gdal  # noqa: F401
    except ImportError:
        print("  [SKIP] GDAL not available")
        return
    z, _ = smooth_surface(n=400, relief=60.0)
    nn, vn = _warp(z, "near")
    bl, vb = _warp(z, "bilinear")
    stn, stb = resampling_stats(nn, vn), resampling_stats(bl, vb)
    check("nearest to UTM: no whole repeated rows (the 1.0.1 test misses it)",
          max(stn["duplicated_row_share"], stn["duplicated_col_share"]) < 0.002,
          f"rows {stn['duplicated_row_share']:.4f}, cols {stn['duplicated_col_share']:.4f}")
    w = audit_resampling(nn, vn)
    check("nearest to UTM: sloped ties flagged", len(w) == 1 and "sloping ground" in w[0],
          f"sloped ties {stn['sloped_tied_share']:.3f}")
    check("bilinear to UTM: no warning", not audit_resampling(bl, vb),
          f"sloped ties {stb['sloped_tied_share']:.4f}")
    rn, rv = _warp(np.round(z), "near")
    check("whole-metre DEM reprojected with nearest: left to the precision audit, no tie warning",
          not any("sloping ground" in x for x in audit_resampling(rn, rv))
          and len(audit_precision(rn, rv)) <= 1)
    flat = z.copy()
    flat[100:300, 100:300] = 1500.0                           # a lake: large exactly flat patch
    check("continuous DEM with a flat lake: not flagged", not audit_resampling(flat, np.ones_like(flat, bool)))


def main():
    test_quantum()
    test_drop()
    test_audit_precision()
    test_sloped_ties()
    print(f"\n{'ALL PASS' if not FAILURES else f'{len(FAILURES)} FAILURE(S)'}")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
