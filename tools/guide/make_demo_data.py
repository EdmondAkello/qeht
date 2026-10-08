"""Synthetic demonstration dataset for the QEHT user guide (invented names and coordinates).

    python3.12 qeht/tools/guide/make_demo_data.py <output folder>

Writes, in EPSG:32737 (WGS 84 / UTM 37S), on a 10 m grid of 1,500 x 1,200 cells:
  demo_dem.tif           terrain: a main valley with tributaries (Gaussian valleys on a tilted
                         surface with Gaussian hills), a flat floodplain reach stored in whole
                         metres (so the flat-method check has something to show), three closed
                         depressions, and a road embankment raised 3 m along a curved centreline
  demo_road.gpkg         road centreline (layer road)
  demo_pour_points.gpkg  hand-placed pour points near the main crossings (field culvert)
  demo_landcover.tif     ESA WorldCover codes (tree on the hills, grass, crops in the valleys,
                         a built-up patch by the road, bare ground on the steepest slopes)
  demo_landcover_2040.tif  scenario: crops near the road become built-up
  demo_soils.gpkg        soil polygons with sand, silt, clay, oc, drain and hsg fields
  demo_rain_zones.gpkg   two rainfall zones (field zone)
  demo_map.tif           mean annual rainfall, 700 to 1,100 mm/yr from south to north
  demo_r.tif             RUSLE R, constant 3,500 MJ mm/(ha h yr)
  demo_waterways.gpkg    "mapped" waterways traced from the generator's valley lines, with names

Every figure in the guide is made from this dataset; nothing in it describes a real place.
"""

import math
import os
import sys

import numpy as np

ROWS, COLS, CS = 1500, 1200, 10.0
X0, Y0 = 512000.0, 9915000.0               # invented false origin (upper-left corner)
SEED = 2026


def valley_lines():
    """Valley centre lines in local metres (x east from the west edge, y north from the south)."""
    ys = np.linspace(0, ROWS * CS, 150)
    main = [(6000 + 1500 * math.sin(y / 4000.0) + 400 * math.sin(y / 1300.0), y) for y in ys]
    tribs = []
    for (y0, side, length, ang) in ((5200, -1, 4200, 0.55), (8800, 1, 4800, 0.5), (11800, -1, 3800, 0.6),
                                    (2600, 1, 3600, 0.45)):
        xm = 6000 + 1500 * math.sin(y0 / 4000.0) + 400 * math.sin(y0 / 1300.0)
        pts = [(xm + side * t * math.cos(ang), y0 + t * math.sin(ang) + 120 * math.sin(t / 700.0))
               for t in np.linspace(0, length, 50)]
        tribs.append(pts)
    return main, tribs


def dist_to_line(X, Y, line):
    a = np.asarray(line)
    d = np.full(X.shape, np.inf)
    for (x0, y0), (x1, y1) in zip(a[:-1], a[1:]):
        dx, dy = x1 - x0, y1 - y0
        L2 = dx * dx + dy * dy
        t = np.clip(((X - x0) * dx + (Y - y0) * dy) / L2, 0, 1)
        d = np.minimum(d, np.hypot(X - (x0 + t * dx), Y - (y0 + t * dy)))
    return d


def road_line():
    xs = np.linspace(150, COLS * CS - 150, 120)
    return [(x, 3200 + 700 * math.sin(x / 2600.0) + 150 * math.sin(x / 900.0)) for x in xs]


def terrain():
    rng = np.random.default_rng(SEED)
    c, r = np.meshgrid(np.arange(COLS), np.arange(ROWS))
    X = (c + 0.5) * CS
    Y = (ROWS - r - 0.5) * CS                       # north up: row 0 is the north edge
    z = 1400.0 + 0.030 * Y + 0.004 * np.abs(X - 6000)
    for _ in range(26):                              # hills
        hx, hy = rng.uniform(0, COLS * CS), rng.uniform(2000, ROWS * CS)
        z += rng.uniform(15, 55) * np.exp(-((X - hx) ** 2 + (Y - hy) ** 2) / (2 * rng.uniform(400, 1000) ** 2))
    main, tribs = valley_lines()
    dm = dist_to_line(X, Y, main)
    z -= 70.0 * np.exp(-(dm / 900.0) ** 2) + 0.004 * np.maximum(1500 - dm, 0)
    for t in tribs:
        dt = dist_to_line(X, Y, t)
        z -= 28.0 * np.exp(-(dt / 380.0) ** 2) + 0.004 * np.maximum(500 - dt, 0)
    # flat floodplain reach on the main valley, stored in whole metres
    rows_ = np.arange(ROWS)
    centre = np.argmin(dm, axis=1)                         # valley centre cell of every row
    bottom = z[rows_, centre]
    bottom = np.convolve(np.pad(bottom, 100, mode="edge"), np.ones(201) / 201, mode="valid")[:, None]
    shape = (bottom - 1.0 + 0.002 * dm + 0.12 * np.maximum(dm - 350.0, 0)    # 700 m floor, 1:8 sides
             + 0.08 * np.maximum(2000.0 - Y, 0) + 0.08 * np.maximum(Y - 4100.0, 0))   # tapered ends
    fp = (dm < 1500) & (Y < 5500) & (shape < z)
    floor_ = (dm <= 350) & (Y > 2000) & (Y < 4100)
    z = np.where(fp, np.where(floor_, np.round(shape), shape), z)
    for (px, py, depth) in ((3000, 9000, 5.0), (9300, 12500, 4.0), (8200, 6200, 3.5)):   # depressions
        z -= depth * np.exp(-((X - px) ** 2 + (Y - py) ** 2) / (2 * 120.0 ** 2))
    # road embankment, 3 m high, 16 m wide
    dr = dist_to_line(X, Y, road_line())
    z = np.where(dr <= 8.0, z + 3.0, z)
    return np.round(z, 2), X, Y, dm


def write_tif(path, a, gdal_type=None, nodata=None):
    from osgeo import gdal, osr
    gdal.UseExceptions()
    srs = osr.SpatialReference(); srs.ImportFromEPSG(32737)
    ds = gdal.GetDriverByName("GTiff").Create(path, COLS, ROWS, 1, gdal_type or gdal.GDT_Float32,
                                              options=["COMPRESS=DEFLATE", "TILED=YES"])
    ds.SetGeoTransform((X0, CS, 0.0, Y0, 0.0, -CS))
    ds.SetProjection(srs.ExportToWkt())
    b = ds.GetRasterBand(1)
    b.WriteArray(a)
    if nodata is not None:
        b.SetNoDataValue(nodata)
    ds = None


def write_vector(path, layer, geom_type, fields, features):
    from osgeo import ogr, osr
    ogr.UseExceptions()
    srs = osr.SpatialReference(); srs.ImportFromEPSG(32737)
    if os.path.exists(path):
        os.remove(path)
    ds = ogr.GetDriverByName("GPKG").CreateDataSource(path)
    lyr = ds.CreateLayer(layer, srs, geom_type)
    for name, kind in fields:
        lyr.CreateField(ogr.FieldDefn(name, {"text": ogr.OFTString, "real": ogr.OFTReal,
                                             "int": ogr.OFTInteger}[kind]))
    for coords, attrs in features:
        f = ogr.Feature(lyr.GetLayerDefn())
        for k, v in attrs.items():
            f.SetField(k, v)
        if geom_type == ogr.wkbPoint:
            g = ogr.Geometry(ogr.wkbPoint); g.AddPoint_2D(*coords)
        elif geom_type == ogr.wkbLineString:
            g = ogr.Geometry(ogr.wkbLineString)
            for x, y in coords:
                g.AddPoint_2D(x, y)
        else:
            g = ogr.Geometry(ogr.wkbPolygon)
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for x, y in coords + [coords[0]]:
                ring.AddPoint_2D(x, y)
            g.AddGeometry(ring)
        f.SetGeometry(g)
        lyr.CreateFeature(f)
    ds = None


def world(pts):
    """Local metres (y north from the south edge) -> EPSG:32737 coordinates."""
    return [(X0 + x, Y0 - ROWS * CS + y) for x, y in pts]


def main(out):
    from osgeo import ogr
    os.makedirs(out, exist_ok=True)
    z, X, Y, dm = terrain()
    write_tif(os.path.join(out, "demo_dem.tif"), z.astype(np.float32))
    # land cover from elevation, slope and distance to the valley
    gy, gx = np.gradient(z, CS)
    slope = np.hypot(gx, gy)
    rel = z - np.percentile(z, 50)
    lc = np.full(z.shape, 30, np.uint8)                              # grassland
    lc[(dm < 1200) & (slope < 0.06)] = 40                            # cropland in the valleys
    lc[rel > 110] = 10                                                # tree cover on the hills
    lc[slope > 0.18] = 60                                            # bare on the steepest slopes
    dr = dist_to_line(X, Y, road_line())
    town = (dr < 250) & (np.abs(X - 3000) < 600)
    lc[town] = 50                                                    # built-up by the road
    write_tif(os.path.join(out, "demo_landcover.tif"), lc, __import__("osgeo.gdal").gdal.GDT_Byte, 0)
    scn = lc.copy()
    scn[(dr < 600) & (lc == 40)] = 50                                # crops near the road built over
    write_tif(os.path.join(out, "demo_landcover_2040.tif"), scn, __import__("osgeo.gdal").gdal.GDT_Byte, 0)
    write_tif(os.path.join(out, "demo_map.tif"), (700.0 + 400.0 * Y / (ROWS * CS)).astype(np.float32))
    write_tif(os.path.join(out, "demo_r.tif"), np.full(z.shape, 3500.0, np.float32))
    W, H = COLS * CS, ROWS * CS
    write_vector(os.path.join(out, "demo_road.gpkg"), "road", ogr.wkbLineString, [("name", "text")],
                 [(world(road_line()), {"name": "Demo Road"})])
    soils = [([(0, 0), (W * 0.45, 0), (W * 0.40, H), (0, H)],
              {"unit": "DS1", "sand": 55.0, "silt": 25.0, "clay": 20.0, "oc": 1.4, "drain": "W", "hsg": "B"}),
             ([(W * 0.45, 0), (W, 0), (W, H * 0.55), (W * 0.42, H * 0.6)],
              {"unit": "DS2", "sand": 30.0, "silt": 30.0, "clay": 40.0, "oc": 1.1, "drain": "M", "hsg": "C"}),
             ([(W * 0.42, H * 0.6), (W, H * 0.55), (W, H), (W * 0.40, H)],
              {"unit": "DS3", "sand": 20.0, "silt": 35.0, "clay": 45.0, "oc": 2.0, "drain": "I", "hsg": "D"})]
    write_vector(os.path.join(out, "demo_soils.gpkg"), "soils", ogr.wkbPolygon,
                 [("unit", "text"), ("sand", "real"), ("silt", "real"), ("clay", "real"), ("oc", "real"),
                  ("drain", "text"), ("hsg", "text")],
                 [(world(p), a) for p, a in soils])
    write_vector(os.path.join(out, "demo_rain_zones.gpkg"), "zones", ogr.wkbPolygon, [("zone", "text")],
                 [(world([(0, 0), (W, 0), (W, H * 0.5), (0, H * 0.5)]), {"zone": "Zone I (lowland)"}),
                  (world([(0, H * 0.5), (W, H * 0.5), (W, H), (0, H)]), {"zone": "Zone II (upland)"})])
    main_v, tribs = valley_lines()
    names = ["Demo River", "Kijito East", "Mto Mdogo", "Kijito North", "Lower Stream"]
    lines = [main_v] + tribs
    feats = []
    for nm, ln in zip(names, lines):
        pts = [(x + 15 * math.sin(i / 7.0), y) for i, (x, y) in enumerate(ln[::3])]   # map-like wobble
        feats.append((world(pts), {"name": nm}))
    write_vector(os.path.join(out, "demo_waterways.gpkg"), "waterways", ogr.wkbLineString,
                 [("name", "text")], feats)
    # pour points: where the road meets the main valley and the southern tributary, a few metres off
    rd = np.asarray(road_line())
    pts = []
    for k, ln in enumerate((main_v, tribs[3])):
        dd = dist_to_line(rd[:, 0], rd[:, 1], ln)
        i = int(np.argmin(dd))
        pts.append((world([(rd[i, 0] + 12.0, rd[i, 1] + 25.0)])[0], {"culvert": f"C-{k + 1:02d}"}))
    write_vector(os.path.join(out, "demo_pour_points.gpkg"), "pour_points", ogr.wkbPoint,
                 [("culvert", "text")], pts)
    print("demo data in", out, "| DEM", z.shape, f"{np.nanmin(z):.0f}-{np.nanmax(z):.0f} m")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "demo")
