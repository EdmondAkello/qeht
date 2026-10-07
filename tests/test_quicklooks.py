# -*- coding: utf-8 -*-
"""v0.18 validation: raster quicklooks (A6). NumPy; GDAL for the end-to-end part.

    python -m qeht.tests.test_quicklooks
"""

import json
import math
import os
import struct
import sys
import tempfile
import zlib

import numpy as np

from ..core.report.quicklooks import (write_png, read_png_size, write_world_file, downsample,
                                      factor_for, render_classified, render_continuous,
                                      make_quicklook, standard_items, INDEX_FIELDS, hillshade)

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def decode_png(path):
    """Minimal decoder for the files write_png makes (RGBA, filter 0)."""
    with open(path, "rb") as f:
        b = f.read()
    assert b[:8] == b"\x89PNG\r\n\x1a\n"
    o, idat = 8, b""
    w = h = None
    while o < len(b):
        n = struct.unpack(">I", b[o:o + 4])[0]
        tag, data = b[o + 4:o + 8], b[o + 8:o + 8 + n]
        crc = struct.unpack(">I", b[o + 8 + n:o + 12 + n])[0]
        assert crc == zlib.crc32(tag + data) & 0xFFFFFFFF
        if tag == b"IHDR":
            w, h = struct.unpack(">II", data[:8])
        elif tag == b"IDAT":
            idat += data
        o += 12 + n
    raw = zlib.decompress(idat)
    rows = [raw[r * (w * 4 + 1) + 1:(r + 1) * (w * 4 + 1)] for r in range(h)]
    return np.frombuffer(b"".join(rows), np.uint8).reshape(h, w, 4)


def test_basics():
    print("\n1. PNG, world file, downsampling, colours")
    tmp = tempfile.mkdtemp()
    a = np.random.default_rng(1).integers(0, 256, (7, 5, 4)).astype(np.uint8)
    p = write_png(os.path.join(tmp, "a.png"), a)
    check("PNG round trip (CRCs valid, pixels identical)", np.array_equal(decode_png(p), a)
          and read_png_size(p) == (5, 7))
    w = write_world_file(os.path.join(tmp, "a.pgw"), (1000.0, 30.0, 0.0, 2000.0, 0.0, -30.0))
    v = [float(x) for x in open(w).read().split()]
    check("world file: pixel size, then the CENTRE of the upper-left pixel",
          v == [30.0, 0.0, 0.0, -30.0, 1015.0, 1985.0])
    cls = np.array([[1, 1, 2, 2, 3], [1, 1, 2, 2, 3], [4, 4, 5, 5, 3]])
    d = downsample(cls, 2, "classified")
    check("classified downsampling keeps class values (no blending)",
          set(np.unique(d)) <= set(np.unique(cls)) and d.shape == (2, 3))
    c = np.arange(16, dtype=float).reshape(4, 4)
    c[0, 0] = np.nan
    m = downsample(c, 2, "continuous")
    check("continuous downsampling: block mean of the finite cells",
          math.isclose(m[0, 0], (1 + 4 + 5) / 3) and math.isclose(m[1, 1], (10 + 11 + 14 + 15) / 4))
    check("factor keeps the long side within the limit",
          factor_for((10000, 3000), 4096) == 3 and factor_for((100, 100), 4096) == 1)
    rgba = render_classified(np.array([[0, 1], [2, 9]]), {1: (10, 20, 30), 2: (40, 50, 60)})
    check("class colours exact; NoData and unknown classes transparent",
          tuple(rgba[0, 1]) == (10, 20, 30, 255) and rgba[0, 0, 3] == 0 and rgba[1, 1, 3] == 0)
    rc = render_continuous(np.array([[0.0, 99.0, np.nan]]), "blues", 0.0, 99.0, log=True)
    check("continuous ramp: ends of the stretch get the ramp ends; NaN transparent",
          tuple(rc[0, 0, :3]) == (240, 246, 251) and tuple(rc[0, 1, :3]) == (8, 48, 107)
          and rc[0, 2, 3] == 0)
    z = np.add.outer(np.zeros(5), np.arange(5.0))      # plane rising to the east
    hs = hillshade(z, 1.0, 1.0)
    exp = 0.5 + 0.5 * math.cos(math.radians(45))       # 45 deg west-facing plane, sun NW 45 deg
    check("hillshade of a 45-degree west-facing plane: uniform 0.854 inside (ESRI convention)",
          np.allclose(hs[1:-1, 1:-1], exp), f"{hs[2, 2]:.4f}")


def _write_tif(path, arr, gt, dtype, nodata, ct=None, names=None):
    from osgeo import gdal, osr
    gdal.UseExceptions()
    ds = gdal.GetDriverByName("GTiff").Create(path, arr.shape[1], arr.shape[0], 1, dtype)
    ds.SetGeoTransform(gt)
    srs = osr.SpatialReference(); srs.ImportFromEPSG(32737)
    ds.SetProjection(srs.ExportToWkt())
    b = ds.GetRasterBand(1)
    b.SetNoDataValue(nodata)
    if ct:
        t = gdal.ColorTable()
        for k, rgb in ct.items():
            t.SetColorEntry(k, tuple(rgb) + (255,))
        b.SetColorTable(t)
    if names:
        b.SetCategoryNames(names)
    b.WriteArray(arr)
    ds = None


def test_end_to_end():
    print("\n2. Quicklooks from GeoTIFFs")
    try:
        from osgeo import gdal  # noqa: F401
    except ImportError:
        print("  [SKIP] GDAL not available")
        return
    tmp = tempfile.mkdtemp()
    gt = (500000.0, 10.0, 0.0, 9900000.0, 0.0, -10.0)
    cls = np.zeros((300, 500), np.uint8)
    cls[:, :250] = 1
    cls[:, 250:] = 4
    cls[100:200, 100:400] = 5
    ct = {1: (26, 150, 65), 4: (253, 174, 97), 5: (215, 25, 28)}
    names = ["", "Very low", "", "", "High", "Very high"]
    _write_tif(os.path.join(tmp, "spi_class.tif"), cls, gt, 1, 0, ct, names)
    row = make_quicklook(dict(name="spi_class", title="SPI classes",
                              path=os.path.join(tmp, "spi_class.tif"), kind="classified"),
                         os.path.join(tmp, "ql"), max_px=128)
    leg = json.load(open(row["legend_json"]))
    check("legend classes = the colour table entries present, with their names",
          [(c["value"], c["label"], c["colour"]) for c in leg["classes"]]
          == [(1, "Very low", "#1a9641"), (4, "High", "#fdae61"), (5, "Very high", "#d7191c")])
    img = decode_png(row["png"])
    k = leg["downsample_factor"]
    check("PNG within the size limit; only the class colours appear",
          max(img.shape[:2]) <= 128 and {tuple(px[:3]) for px in img.reshape(-1, 4)} <= set(ct.values()),
          f"{img.shape[1]} x {img.shape[0]}, factor {k}")
    v = [float(x) for x in open(row["world_file"]).read().split()]
    x0 = v[4] - v[0] / 2
    y0 = v[5] - v[3] / 2
    check("world file + PNG size reproduce the source extent (to one PNG pixel)",
          math.isclose(x0, gt[0]) and math.isclose(y0, gt[3])
          and abs((x0 + v[0] * img.shape[1]) - (gt[0] + 10 * 500)) <= v[0]
          and abs((y0 + v[3] * img.shape[0]) - (gt[3] - 10 * 300)) <= abs(v[3])
          and row["xmin"] == gt[0] and row["ymax"] == gt[3] and row["crs_epsg"] == "32737")
    acc = np.exp(np.random.default_rng(2).uniform(0, 10, (200, 200)))
    _write_tif(os.path.join(tmp, "fac.tif"), acc.astype(np.float32), gt, 6, -1.0)
    dem = 1500 + np.add.outer(np.arange(200.0), np.arange(200.0))
    _write_tif(os.path.join(tmp, "dem.tif"), dem.astype(np.float32), gt, 6, -9999.0)
    items = standard_items(os.path.join(tmp, "dem.tif"), os.path.join(tmp, "fac.tif"))
    rows = [make_quicklook(it, os.path.join(tmp, "ql"), 4096) for it in items]
    legs = [json.load(open(r["legend_json"])) for r in rows]
    check("relief and log-stretched accumulation written; stretch recorded in the legend",
          [r["name"] for r in rows] == ["relief", "flow_accumulation"]
          and legs[1]["stretch_rule"].startswith("log10") and legs[0]["ramp"] == "terrain"
          and all(r["kind"] == "continuous" for r in rows))
    check("index row has every INDEX_FIELDS column", set(r for r, _ in INDEX_FIELDS) == set(rows[0]))
    from ..core.interop.field_dictionary import OPTIONAL_LAYERS, METADATA_KEYS
    check("'rasters' table and quicklook_params_json in the contract",
          [f[0] for f in OPTIONAL_LAYERS["rasters"][1]] == [n for n, _ in INDEX_FIELDS]
          and "quicklook_params_json" in [k for k, _ in METADATA_KEYS])


def main(argv=None):
    test_basics()
    test_end_to_end()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
