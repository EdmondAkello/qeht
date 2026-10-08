# -*- coding: utf-8 -*-
"""v0.21 validation: Prepare DEM for hydrology (F9).
NumPy + the GDAL bindings (the warp is GDAL's), no QGIS. Synthetic 1-arc-second
tiles of an analytic surface near 37.5 E, 1 S.

    python -m qeht.tests.test_prepare_dem
"""

import json
import math
import os
import sys
import tempfile

import numpy as np

from ..core.conditioning.prepare import (prepare_dem, utm_epsg, metres_per_degree,
                                         transform_bounds, OUT_NODATA)
from ..core.raster import audit_resampling, raster_tags, read_dem

FAILURES = []
ARC1 = 1.0 / 3600.0
N = 240                                   # cells per tile side


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def surface(lon, lat):
    """Smooth analytic terrain (m) in lon / lat degrees."""
    x = (lon - 37.4) * 111300.0
    y = (lat + 1.0) * 110600.0
    return 1500.0 + 25.0 * np.sin(x / 900.0) + 18.0 * np.cos(y / 700.0) + 0.004 * x - 0.003 * y


def write_tile(path, lon0, lat0, n=N, data=None, nodata=None):
    from osgeo import gdal, osr
    gdal.UseExceptions()
    ds = gdal.GetDriverByName("GTiff").Create(path, n, n, 1, gdal.GDT_Float32)
    ds.SetGeoTransform((lon0, ARC1, 0.0, lat0, 0.0, -ARC1))
    s = osr.SpatialReference(); s.ImportFromEPSG(4326)
    ds.SetProjection(s.ExportToWkt())
    if data is None:
        c, r = np.meshgrid(np.arange(n) + 0.5, np.arange(n) + 0.5)
        data = surface(lon0 + c * ARC1, lat0 - r * ARC1)
    ds.GetRasterBand(1).WriteArray(data.astype(np.float32))
    if nodata is not None:
        ds.GetRasterBand(1).SetNoDataValue(nodata)
    ds = None
    return path


def to_lonlat(path):
    """Cell-centre lon / lat of a projected raster and its values (NaN = NoData)."""
    from osgeo import osr
    a, v, info = read_dem(path)
    gt = info.geotransform
    s = osr.SpatialReference(); s.SetFromUserInput(info.projection_wkt)
    g = osr.SpatialReference(); g.ImportFromEPSG(4326)
    for x in (s, g):
        x.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    t = osr.CoordinateTransformation(s, g)
    c, r = np.meshgrid(np.arange(info.cols) + 0.5, np.arange(info.rows) + 0.5)
    X, Y = gt[0] + c * gt[1], gt[3] + r * gt[5]
    ll = np.array(t.TransformPoints(np.c_[X.ravel(), Y.ravel()].tolist()))
    return (ll[:, 0].reshape(X.shape), ll[:, 1].reshape(X.shape),
            np.where(v, a, np.nan), info)


def tiles(tmp):
    lon0, lat0 = 37.40, -0.96
    a = write_tile(os.path.join(tmp, "tile_a.tif"), lon0, lat0)
    b = write_tile(os.path.join(tmp, "tile_b.tif"), lon0 + N * ARC1, lat0)
    return a, b, (lon0, lat0 - N * ARC1, lon0 + 2 * N * ARC1, lat0)


def test_merge(tmp):
    print("\n1. Two adjacent tiles merge seamlessly; auto UTM")
    a, b, ext = tiles(tmp)
    out = os.path.join(tmp, "merged.tif")
    log = prepare_dem([a, b], out, log_path=os.path.join(tmp, "merged.json"))
    lon, lat, z, info = to_lonlat(out)
    inner = ((lon > ext[0] + 3 * ARC1) & (lon < ext[2] - 3 * ARC1)
             & (lat > ext[1] + 3 * ARC1) & (lat < ext[3] - 3 * ARC1))
    seam = inner & (np.abs(lon - (ext[0] + N * ARC1)) < 3 * ARC1)
    err = np.abs(z - surface(lon, lat))
    check("no NoData inside the merged extent (no seam rows or columns)",
          np.isfinite(z[inner]).all() and seam.sum() > 0, f"{int(np.isnan(z[inner]).sum())} NaN")
    check("bilinear values follow the analytic surface (max error < 0.5 m, also at the seam)",
          np.nanmax(err[inner]) < 0.5 and np.nanmax(err[seam]) < 0.5,
          f"max {np.nanmax(err[inner]):.3f} m, seam {np.nanmax(err[seam]):.3f} m")
    check("auto CRS: UTM 37S (EPSG:32737) for a centre near 37.5 E, 1 S",
          "EPSG:32737" in log["target_crs"] and "32737" in info.projection_wkt.replace('"', ""))
    m_lon, m_lat = metres_per_degree(-1.0)
    check("native cell in metres: 1 arc-second -> max(30.9, 30.7) rounded to 0.1 m",
          log["cell_m"] == round(max(ARC1 * m_lon, ARC1 * m_lat), 1)
          and math.isclose(info.cell_width, log["cell_m"]), f"{log['cell_m']}")
    gt = info.geotransform
    check("target-aligned pixels (origin a multiple of the cell size)",
          math.isclose(gt[0] / gt[1], round(gt[0] / gt[1]), abs_tol=1e-6)
          and math.isclose(gt[3] / gt[1], round(gt[3] / gt[1]), abs_tol=1e-6))
    with open(os.path.join(tmp, "merged.json"), encoding="utf-8") as f:
        j = json.load(f)
    check("one zone: no zone note", j["zone_note"] == "")
    check("log: inputs, CRS, cell, resampling, both audits, native estimate",
          j["inputs"] == [a, b] and j["resampling"] == "bilinear" and not j["audit_input"]["warnings"]
          and not j["audit_output"]["warnings"] and len(j["native_estimate_m"]) == 2)
    check("utm_epsg anywhere: 37.5 E 1 S -> 32737; 15 E 50 N -> 32633; 74 W 41 N -> 32618; "
          "151 E 34 S -> 32756; 179.9 W -> zone 1",
          utm_epsg(37.5, -1.0) == 32737 and utm_epsg(15.0, 50.0) == 32633
          and utm_epsg(-74.0, 41.0) == 32618 and utm_epsg(151.0, -34.0) == 32756
          and utm_epsg(-179.9, 10.0) == 32601)
    check("UTM exceptions: Norway 32V (5 E 60 N), Svalbard 33X (15 E 78 N); UPS beyond 84 N / 80 S",
          utm_epsg(5.0, 60.0) == 32632 and utm_epsg(15.0, 78.0) == 32633
          and utm_epsg(0.0, 85.0) == 32661 and utm_epsg(0.0, -82.0) == 32761)


def test_resampling(tmp):
    print("\n2. Resampling audit: bilinear passes, nearest fails")
    from osgeo import gdal
    a, b, _ = tiles(tmp)
    out = os.path.join(tmp, "fine.tif")
    log = prepare_dem([a, b], out, cell_m=20.0)
    check("bilinear to a 20 m grid (finer than the 30.9 m native): no repeat warning",
          not log["audit_output"]["warnings"], "; ".join(log["audit_output"]["warnings"]))
    nn = os.path.join(tmp, "nearest.tif")
    gdal.Warp(nn, [a, b], options=gdal.WarpOptions(dstSRS="EPSG:32737", xRes=20.0, yRes=20.0,
                                                   resampleAlg="near", dstNodata=OUT_NODATA))
    z, v, _ = read_dem(nn)
    check("the same warp with nearest neighbour (control) fails audit_resampling",
          len(audit_resampling(z, v)) == 1)
    try:
        prepare_dem([a], os.path.join(tmp, "x.tif"), resampling="near")
        refused = False
    except ValueError as e:
        refused = "nearest" in str(e)
    check("nearest is refused with the reason", refused)
    log2 = prepare_dem([nn], os.path.join(tmp, "from_nn.tif"))
    check("an input that already repeats rows is reported as damaged upstream (get the native "
          "product); the native estimate is coarser than the input cell",
          log2["input_damaged_upstream"] and "native product" in log2["advice"]
          and log2["native_estimate_m"][0] > 20.0 * 1.2, str(log2["native_estimate_m"]))


def test_nodata_clip_tags(tmp):
    print("\n3. NoData override, clip, tags")
    lon0, lat0 = 37.40, -0.96
    c, r = np.meshgrid(np.arange(N) + 0.5, np.arange(N) + 0.5)
    z = surface(lon0 + c * ARC1, lat0 - r * ARC1)
    z[:, :60] = 0.0                                      # undeclared NoData zeros on the west
    t = write_tile(os.path.join(tmp, "zeros.tif"), lon0, lat0, data=z)
    log0 = prepare_dem([t], os.path.join(tmp, "zeros_kept.tif"))
    check("without the override the input audit warns about the zeros",
          any("undeclared NoData" in w for w in log0["audit_input"]["warnings"]))
    out = os.path.join(tmp, "zeros_nd.tif")
    log = prepare_dem([t], out, nodata_override=0.0, source_text="synthetic 1\" tiles")
    lon, lat, zz, _ = to_lonlat(out)
    west = lon < lon0 + 55 * ARC1
    east = (lon > lon0 + 70 * ARC1) & (lon < lon0 + (N - 3) * ARC1) & (lat < lat0 - 3 * ARC1) \
        & (lat > lat0 - (N - 3) * ARC1)
    check("override 0: the zeros become NoData, the terrain is kept, no zero in the output",
          np.isnan(zz[west]).all() and np.isfinite(zz[east]).all()
          and np.nanmin(zz) > 1400.0, f"min {np.nanmin(zz):.1f}")
    tags = raster_tags(out)
    check("tags QEHT_DEM_SOURCE and QEHT_DEM_PREP written",
          tags.get("QEHT_DEM_SOURCE") == "synthetic 1\" tiles"
          and tags.get("QEHT_DEM_PREP", "").startswith("QEHT prepare DEM")
          and "override 0" in tags["QEHT_DEM_PREP"])
    a, b, ext = tiles(tmp)
    box = (37.45, -1.0, 37.48, -0.98)
    cl = os.path.join(tmp, "clip.tif")
    log = prepare_dem([a, b], cl, clip_bounds=box, clip_crs="EPSG:4326", buffer_m=500.0)
    _, _, info = read_dem(cl)
    tb = transform_bounds(box, "EPSG:4326", "EPSG:32737")
    gt = info.geotransform
    ob = (gt[0], gt[3] + info.rows * gt[5], gt[0] + info.cols * gt[1], gt[3])
    cell = info.cell_width
    check("clip: the box in UTM plus the 500 m buffer, to within one cell (-tap)",
          all(abs(o - (t_ + s * 500.0)) <= cell + 1e-6
              for o, t_, s in zip(ob, tb, (-1, -1, 1, 1))),
          f"{[round(o - (t_ + s * 500.0), 1) for o, t_, s in zip(ob, tb, (-1, -1, 1, 1))]}")
    w0 = write_tile(os.path.join(tmp, "zone_w.tif"), 35.97, -0.96)
    lz = prepare_dem([w0], os.path.join(tmp, "zones.tif"))
    check("an extent across a zone boundary (36 E) gets a zone note", "2 UTM zones" in lz["zone_note"],
          lz["zone_note"])
    try:
        prepare_dem([a], os.path.join(tmp, "geo.tif"), target="EPSG:4326")
        refused = False
    except ValueError as e:
        refused = "geographic" in str(e)
    check("a geographic target CRS is refused", refused)
    log = prepare_dem([a], os.path.join(tmp, "arc.tif"), target=21037)
    _, _, ai = read_dem(os.path.join(tmp, "arc.tif"))
    check("Arc 1960 / UTM 37S (EPSG:21037) as the target", "Arc 1960" in ai.projection_wkt)


def main(argv=None):
    tmp = tempfile.mkdtemp(prefix="qeht_prep_")
    test_merge(tmp)
    test_resampling(tmp)
    test_nodata_clip_tags(tmp)
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
