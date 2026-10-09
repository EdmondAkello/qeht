from osgeo import ogr, osr
import math, sys
ogr.UseExceptions()
uid = sys.argv[1] if len(sys.argv) > 1 else "X014"
ds = ogr.Open("run/layers/flowpaths.gpkg"); l = ds.GetLayer(); l.SetAttributeFilter(f"outlet_uid = '{uid}'")
f = next(iter(l)); g = f.GetGeometryRef(); pts = g.GetPoints()
cds = ogr.Open("run/layers/crossings.gpkg"); cr = cds.GetLayer(); cr.SetAttributeFilter(f"outlet_uid = '{uid}'")
c = next(iter(cr)); ox, oy = c["outlet_x"], c["outlet_y"]
if math.hypot(pts[0][0] - ox, pts[0][1] - oy) > math.hypot(pts[-1][0] - ox, pts[-1][1] - oy):
    pts = pts[::-1]                                   # outlet first
cum = [0.0]
for a, b in zip(pts[:-1], pts[1:]):
    cum.append(cum[-1] + math.hypot(b[0] - a[0], b[1] - a[1]))
L = cum[-1]
def at(d):
    for i in range(1, len(cum)):
        if cum[i] >= d:
            t = (d - cum[i - 1]) / (cum[i] - cum[i - 1] or 1)
            return (pts[i - 1][0] + t * (pts[i][0] - pts[i - 1][0]), pts[i - 1][1] + t * (pts[i][1] - pts[i - 1][1]))
    return pts[-1][:2]
L10, L85, ov = f["lfp_L10_m"], f["lfp_L85_m"], f["lfp_overland_m"]
print("geom length", round(L), "lfp_length_m", f["lfp_length_m"], "L10", L10, "L85", L85, "overland", ov)
rows = [("10", "10 %", at(L10)), ("85", "85 %", at(L85)), ("head", "channel head", at(f["lfp_length_m"] - ov))]
import os
os.makedirs("meta", exist_ok=True)
if os.path.exists("meta/lfp_points.gpkg"):
    os.remove("meta/lfp_points.gpkg")
drv = ogr.GetDriverByName("GPKG"); out = drv.CreateDataSource("meta/lfp_points.gpkg")
srs = osr.SpatialReference(); srs.ImportFromEPSG(32737)
ol = out.CreateLayer("pts", srs, ogr.wkbPoint)
ol.CreateField(ogr.FieldDefn("kind", ogr.OFTString)); ol.CreateField(ogr.FieldDefn("lab", ogr.OFTString))
for k, lab, (x, y) in rows:
    ft = ogr.Feature(ol.GetLayerDefn()); ft["kind"] = k; ft["lab"] = lab
    p = ogr.Geometry(ogr.wkbPoint); p.AddPoint_2D(x, y); ft.SetGeometry(p); ol.CreateFeature(ft)
out = None
