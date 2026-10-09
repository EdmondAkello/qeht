"""Demonstration site for the QEHT user guide (invented terrain, names and coordinates).

    python3.12 qeht/tools/guide/make_demo_site.py <output folder>

Built on examples/example_dem.tif (synthetic, EPSG:32737, 30 m). Writes:
  demo_dem.tif            the example DEM with an old lake bed on the south-east plain (flat ground);
                          the undeclared NoData strip is kept: the guide shows the check
  demo_road.gpkg          "Demo Road": across the slope, 2.5 km beside a main valley, then over the plain
  demo_landcover.tif      ESA WorldCover codes from elevation, slope and noise
  demo_landcover_2040.tif scenario: cropland within 1 km of the road becomes built-up
  demo_soils.gpkg         three soil units with sand, silt, clay, oc, drain and hsg fields
  demo_rain_zones.gpkg    two rainfall zones (field zone)
  demo_map.tif            mean annual rainfall, 760-1,140 mm/yr, wetter to the north and uphill
  demo_r.tif              RUSLE R, 3,200-3,800 MJ mm/(ha h yr)
  demo_rivers.gpkg        "mapped" rivers (field name), >= 8 km2, traced on the DEM plus a smooth error
                          field so that they depart from the DEM streams where a real map would
Nothing describes a real place.
"""

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.dirname(REPO))
from qeht.core.conditioning.fill import fill_depressions  # noqa: E402
from qeht.core.flow.direction import d8_direction  # noqa: E402
from qeht.core.flow.accumulation import flow_accumulation  # noqa: E402

RIVER_NAMES = ["Mto Kijani", "Mto Mawe", "Mto Jua", "Mto Nyota", "Mto Upepo", "Mto Mvua"]


def read(path):
    from osgeo import gdal
    gdal.UseExceptions()
    ds = gdal.Open(path)
    return ds.GetRasterBand(1).ReadAsArray().astype(np.float64), ds.GetGeoTransform(), ds.GetProjection()


def write_raster(path, a, gt, wkt, dtype="float32", nodata=None):
    from osgeo import gdal
    t = {"float32": gdal.GDT_Float32, "uint8": gdal.GDT_Byte}[dtype]
    ds = gdal.GetDriverByName("GTiff").Create(path, a.shape[1], a.shape[0], 1, t,
                                              options=["COMPRESS=DEFLATE", "TILED=YES"])
    ds.SetGeoTransform(gt); ds.SetProjection(wkt)
    b = ds.GetRasterBand(1)
    if nodata is not None:
        b.SetNoDataValue(nodata)
    b.WriteArray(a.astype(np.uint8 if dtype == "uint8" else np.float32))
    ds = None


def write_vector(path, layer, kind, fields, rows, wkt):
    from osgeo import ogr, osr
    ogr.UseExceptions()
    if os.path.exists(path):
        os.remove(path)
    ds = ogr.GetDriverByName("GPKG").CreateDataSource(path)
    srs = osr.SpatialReference(); srs.ImportFromWkt(wkt)
    lyr = ds.CreateLayer(layer, srs, {"line": ogr.wkbLineString, "polygon": ogr.wkbPolygon}[kind])
    for name, t in fields:
        lyr.CreateField(ogr.FieldDefn(name, {"text": ogr.OFTString, "real": ogr.OFTReal}[t]))
    for geom, attrs in rows:
        f = ogr.Feature(lyr.GetLayerDefn())
        for k, v in attrs.items():
            f.SetField(k, v)
        if kind == "line":
            g = ogr.Geometry(ogr.wkbLineString)
            for x, y in geom:
                g.AddPoint_2D(float(x), float(y))
        else:
            g = ogr.Geometry(ogr.wkbPolygon); ring = ogr.Geometry(ogr.wkbLinearRing)
            for x, y in geom + [geom[0]]:
                ring.AddPoint_2D(float(x), float(y))
            g.AddGeometry(ring)
        f.SetGeometry(g); lyr.CreateFeature(f)
    ds = None


def main(out):
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(2026)
    z, gt, wkt = read(os.path.join(REPO, "examples", "example_dem.tif"))
    rows, cols = z.shape
    cs = gt[1]
    valid = z > 0
    x0, y0 = gt[0], gt[3]
    W, H = cols * cs, rows * cs
    z = flatten_plain(z, valid, gt)
    write_raster(os.path.join(out, "demo_dem.tif"), z, gt, wkt)

    # drainage of the DEM: the road follows one main valley, so the guide can show a parallel reach
    filled, _, _ = fill_depressions(np.where(valid, z, 0.0), valid, cell_width=cs, cell_height=cs)
    dirs, _ = d8_direction(filled, valid, cs, cs, flat_method="barnes")
    acc, _ = flow_accumulation(dirs, valid)
    trunk = max(trace_rivers(acc, dirs, valid, gt, 8.0), key=lambda ln: sum(
        1 for x, y in ln if x0 + 6500 < x < x0 + 9000 and y0 - 14000 < y < y0 - 9500))
    road = make_road(trunk, x0, y0)
    write_vector(os.path.join(out, "demo_road.gpkg"), "road", "line", [("name", "text")],
                 [(road, {"name": "Demo Road"})], wkt)

    # land cover (WorldCover codes): trees high and on the cone, shrubland on mid slopes,
    # grassland, cropland on the gentle lower ground, bare ground on the steepest slopes and
    # the crater, a built-up town where the road crosses the plain
    zz = np.where(valid, z, np.nan)
    fillv = float(np.nanmean(zz))
    gy, gx = np.gradient(np.where(valid, z, fillv), cs)
    slope = np.hypot(gx, gy)
    noise = blur(rng.standard_normal(z.shape), 6.0)
    noise = noise / noise.std()
    zr = (np.where(valid, z, fillv) - np.nanpercentile(zz, 5)) / (np.nanpercentile(zz, 98) - np.nanpercentile(zz, 5))
    lc = np.full(z.shape, 30, np.uint8)                               # grassland
    lc[(zr + 0.12 * noise) < 0.42] = 40                               # cropland, lower ground
    lc[((zr + 0.12 * noise) > 0.55) & (slope < 0.25)] = 20            # shrubland
    lc[(zr + 0.15 * noise) > 0.72] = 10                               # tree cover
    lc[slope > 0.35] = 60                                             # bare / sparse
    rr, cc = np.mgrid[0:rows, 0:cols]
    X = x0 + (cc + 0.5) * cs
    Y = y0 - (rr + 0.5) * cs
    tx, ty = road[45]
    town = np.hypot(X - tx, Y - ty) < 900 + 250 * noise
    lc[town] = 50
    lc[~valid] = 0
    write_raster(os.path.join(out, "demo_landcover.tif"), lc, gt, wkt, "uint8", 0)
    from osgeo import ogr
    rl = ogr.Geometry(ogr.wkbLineString)
    for x, y in road:
        rl.AddPoint_2D(float(x), float(y))
    # distance to the road (cells): sample the road densely
    pts = np.array([rl.GetPoint_2D(i) for i in range(rl.GetPointCount())])
    dense = np.vstack([np.linspace(a, b, 40) for a, b in zip(pts[:-1], pts[1:])])
    d = np.full(z.shape, np.inf)
    for x, y in dense[::3]:
        d = np.minimum(d, np.hypot(X - x, Y - y))
    lc2 = lc.copy()
    lc2[(d < 1000) & (lc == 40)] = 50
    write_raster(os.path.join(out, "demo_landcover_2040.tif"), lc2, gt, wkt, "uint8", 0)

    # soils: three units across the slope
    yA, yB = y0 - 0.35 * H, y0 - 0.70 * H
    soils = [([(x0, y0), (x0 + W, y0), (x0 + W, yA + 900), (x0, yA - 600)],
              dict(unit="Volcanic loam", sand=42.0, silt=38.0, clay=20.0, oc=2.4, drain="W", hsg="B")),
             ([(x0, yA - 600), (x0 + W, yA + 900), (x0 + W, yB - 400), (x0, yB + 700)],
              dict(unit="Red clay loam", sand=30.0, silt=30.0, clay=40.0, oc=1.5, drain="MW", hsg="C")),
             ([(x0, yB + 700), (x0 + W, yB - 400), (x0 + W, y0 - H), (x0, y0 - H)],
              dict(unit="Plain vertisol", sand=18.0, silt=27.0, clay=55.0, oc=1.1, drain="P", hsg="D"))]
    write_vector(os.path.join(out, "demo_soils.gpkg"), "soils", "polygon",
                 [("unit", "text"), ("sand", "real"), ("silt", "real"), ("clay", "real"),
                  ("oc", "real"), ("drain", "text"), ("hsg", "text")], soils, wkt)

    # rainfall zones and mean annual rainfall
    yz = y0 - 0.45 * H
    zones = [([(x0, y0), (x0 + W, y0), (x0 + W, yz - 1200), (x0, yz + 800)], dict(zone="Zone I (highland)")),
             ([(x0, yz + 800), (x0 + W, yz - 1200), (x0 + W, y0 - H), (x0, y0 - H)], dict(zone="Zone II (plain)"))]
    write_vector(os.path.join(out, "demo_rain_zones.gpkg"), "zones", "polygon", [("zone", "text")], zones, wkt)
    mapv = 760.0 + 260.0 * (1.0 - (rr / rows)) + 120.0 * np.clip(zr, 0, 1)
    write_raster(os.path.join(out, "demo_map.tif"), np.where(valid, mapv, -9999.0), gt, wkt, nodata=-9999.0)
    rv = 3200.0 + 600.0 * (1.0 - rr / rows)
    write_raster(os.path.join(out, "demo_r.tif"), np.where(valid, rv, -9999.0), gt, wkt, nodata=-9999.0)

    # "mapped" rivers: traced on the DEM plus a smooth error field, the way a map digitised from
    # other sources departs from the DEM's own streams (mostly on the plain)
    e = blur(rng.standard_normal(z.shape), 8.0)
    zp = np.where(valid, z + 6.0 * e / e.std(), 0.0)
    fp, _, _ = fill_depressions(zp, valid, cell_width=cs, cell_height=cs)
    dp, _ = d8_direction(fp, valid, cs, cs, flat_method="barnes")
    ap, _ = flow_accumulation(dp, valid)
    lines = trace_rivers(ap, dp, valid, gt, 8.0)
    rivers = [(ln, dict(name=RIVER_NAMES[i % len(RIVER_NAMES)])) for i, ln in enumerate(lines)]
    write_vector(os.path.join(out, "demo_rivers.gpkg"), "rivers", "line", [("name", "text")], rivers, wkt)
    print(f"demo site written to {out}: road {len(road)} vertices, {road_len(road) / 1000:.1f} km, {len(rivers)} mapped rivers, "
          f"land cover classes {sorted(set(np.unique(lc).tolist()) - {0})}")


def flatten_plain(z, valid, gt):
    """An old lake bed on the south-east plain: inside an ellipse the ground blends into a plane
    falling 0.15 % to the south-east, stored in whole metres (so it has real flats). The road
    crosses it, which gives the guide a flat stretch and flat-method sensitivity."""
    rows, cols = z.shape
    cs, x0, y0 = gt[1], gt[0], gt[3]
    rr, cc = np.mgrid[0:rows, 0:cols]
    X = x0 + (cc + 0.5) * cs
    Y = y0 - (rr + 0.5) * cs
    cx, cy, ax, ay = x0 + 17600, y0 - 17500, 3400.0, 1700.0
    q = np.hypot((X - cx) / ax, (Y - cy) / ay)
    w = np.clip((1.25 - q) / 0.45, 0, 1)
    w = w * w * (3 - 2 * w)                                       # smoothstep
    inside = valid & (q < 1.0)
    zref = float(np.median(z[inside]))
    plane = zref - 0.0015 * ((X - cx) * 0.7071 - (Y - cy) * 0.7071)
    out = np.where(valid, (1 - w) * z + w * plane, z)
    return np.where(valid, np.round(out), z)


def trace_rivers(acc, dirs, valid, gt, km2):
    """Lines along the cells draining at least km2, from each channel head down, largest first."""
    from qeht.core.grid import DROW, DCOL
    rows, cols = acc.shape
    cs, x0, y0 = gt[1], gt[0], gt[3]
    big = valid & ((acc + 1) * cs * cs / 1e6 >= km2)
    heads = []
    for r_, c_ in zip(*np.where(big)):
        up = False
        for k in range(8):
            rr2, cc2 = r_ - DROW[k], c_ - DCOL[k]
            if 0 <= rr2 < rows and 0 <= cc2 < cols and big[rr2, cc2] and dirs[rr2, cc2] == k:
                up = True
                break
        if not up:
            heads.append((r_, c_))
    lines, seen = [], np.zeros(acc.shape, bool)
    for r_, c_ in sorted(heads, key=lambda t: -acc[t]):
        line = []
        while 0 <= r_ < rows and 0 <= c_ < cols and big[r_, c_]:
            line.append((x0 + (c_ + 0.5) * cs, y0 - (r_ + 0.5) * cs))
            if seen[r_, c_]:
                break
            seen[r_, c_] = True
            k = dirs[r_, c_]
            if k < 0:
                break
            r_, c_ = r_ + DROW[k], c_ + DCOL[k]
        if len(line) > 8:
            lines.append(line[::2] + [line[-1]])
    return lines


def make_road(trunk, x0, y0):
    """West to east across the slope, then 2.5 km beside the trunk valley (a parallel reach),
    then down onto the southern plain and east across it (flat ground)."""
    t = np.array(trunk)
    seg = t[(t[:, 1] < y0 - 11000) & (t[:, 1] > y0 - 13500)]           # trunk, north to south
    d = np.diff(seg, axis=0)
    n = np.column_stack([d[:, 1], -d[:, 0]]) / np.hypot(d[:, 0], d[:, 1])[:, None]
    n = np.vstack([n, n[-1:]])
    side = seg + 25.0 * n * np.sign(n[:, :1].mean() or 1.0)             # 25 m east of the stream
    west = [(x0 + 1200, y0 - 11000), (x0 + 3000, y0 - 11150), (x0 + 5000, y0 - 10750),
            (seg[0, 0] - 600, seg[0, 1] + 250)]
    east = [(seg[-1, 0] + 900, seg[-1, 1] - 900), (x0 + 12000, y0 - 16600),
            (x0 + 15000, y0 - 17200), (x0 + 18000, y0 - 17300), (x0 + 20200, y0 - 17100)]
    pts = np.vstack([chaikin(np.array(west + [tuple(side[0])]), 3)[:-1], side,
                     chaikin(np.array([tuple(side[-1])] + east), 3)[1:]])
    keep = [0]
    for i in range(1, len(pts)):                                          # vertices about 60 m apart
        if np.hypot(*(pts[i] - pts[keep[-1]])) >= 60.0 or i == len(pts) - 1:
            keep.append(i)
    return [tuple(p) for p in pts[keep]]


def road_len(r):
    a = np.array(r)
    return float(np.hypot(*np.diff(a, axis=0).T).sum())


def chaikin(p, n):
    for _ in range(n):
        q = np.empty((2 * len(p) - 2, 2))
        q[0::2] = 0.75 * p[:-1] + 0.25 * p[1:]
        q[1::2] = 0.25 * p[:-1] + 0.75 * p[1:]
        p = np.vstack([p[:1], q, p[-1:]])
    return p


def blur(a, sigma):
    ky = np.fft.fftfreq(a.shape[0])[:, None]
    kx = np.fft.rfftfreq(a.shape[1])[None, :]
    g = np.exp(-2 * (np.pi * sigma) ** 2 * (kx ** 2 + ky ** 2))
    return np.fft.irfft2(np.fft.rfft2(a) * g, s=a.shape)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "demo_site")
