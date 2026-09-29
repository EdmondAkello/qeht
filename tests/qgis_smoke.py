# -*- coding: utf-8 -*-
"""QGIS-level smoke test: every QEHT Processing tool, run inside QGIS.

Unlike test_core / test_interop this needs a QGIS installation. It runs
the real Processing chain on examples/example_dem.tif:

    fill -> flow direction -> accumulation -> streams -> delineate
         -> longest flow path -> characteristics -> Build HEAS exchange

and checks that every tool completes, outputs load, outlet_uid is
consistent across layers, and the exchange package validates. It also
confirms that the exchange tool refuses a geographic DEM.

Run from the folder that CONTAINS the `qeht` plugin folder:

    QT_QPA_PLATFORM=offscreen python -m qeht.tests.qgis_smoke

or paste into the QGIS Python console:

    from qeht.tests import qgis_smoke; qgis_smoke.main(in_qgis=True)
"""

import os
import sys
import tempfile

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def _start_qgis():
    from qgis.core import QgsApplication
    app = QgsApplication([], False)
    app.initQgis()
    for p in ("/usr/share/qgis/python/plugins",
              os.path.join(QgsApplication.prefixPath(), "python", "plugins")):
        if os.path.isdir(p) and p not in sys.path:
            sys.path.append(p)
    from processing.core.Processing import Processing
    Processing.initialize()
    return app


def main(in_qgis=False):
    app = None if in_qgis else _start_qgis()
    from qgis.core import (QgsApplication, QgsVectorLayer, QgsFeature,
                           QgsGeometry, QgsPointXY, QgsRasterLayer, QgsVectorFileWriter,
                           QgsCoordinateTransformContext, QgsProcessingException)
    import processing
    from ..processing_provider.provider import QehtProvider
    from ..core.interop.heas_exchange import validate_exchange
    from ..core.interop import gpkg

    reg = QgsApplication.processingRegistry()
    if reg.providerById("qeht") is None:
        reg.addProvider(QehtProvider())
    algs = sorted(a.id() for a in reg.providerById("qeht").algorithms())
    print("QEHT algorithms:", ", ".join(algs))
    check("provider loads with 13 algorithms incl. exchange, crossings, burn, relink, soils",
          len(algs) == 13 and all(a in algs for a in (
              "qeht:buildheasexchange", "qeht:crossingcandidates", "qeht:burncrossings",
              "qeht:renumberrelink", "qeht:soilparameters")))

    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dem = os.path.join(here, "examples", "example_dem.tif")
    geo_dem = os.path.join(here, "examples", "synthetic_dem.tif")
    tmp = tempfile.mkdtemp(prefix="qeht_smoke_")
    out = lambda n: os.path.join(tmp, n)

    r = processing.run("qeht:filldepressions", {"DEM": dem, "MIN_SLOPE": 0.0,
                       "OUTPUT": out("fill.tif"), "DEPTH": out("depth.tif")})
    check("fill depressions", os.path.exists(r["OUTPUT"]))
    r = processing.run("qeht:flowdirection", {"DEM": out("fill.tif"), "RESOLVE_FLATS": True,
                       "FLAT_METHOD": 0, "OUTPUT": out("fdr.tif")})
    check("flow direction", os.path.exists(r["OUTPUT"]))
    r = processing.run("qeht:flowaccumulation", {"FDR": out("fdr.tif"), "QUANTITY": 0,
                       "OUTPUT": out("fac.tif")})
    check("flow accumulation", os.path.exists(r["OUTPUT"]))
    r = processing.run("qeht:streamnetwork", {"FDR": out("fdr.tif"), "FAC": out("fac.tif"),
                       "MODE": 0, "THRESHOLD": 500, "STREAMS": out("str.tif"),
                       "ORDER": out("ord.tif")})
    check("stream network + Strahler", os.path.exists(r["STREAMS"]))

    # three pour points near large-accumulation stream cells, with an ID field
    import numpy as np
    from ..core.raster import read_dem
    acc, valid, info = read_dem(out("fac.tif"))
    strm, _, _ = read_dem(out("str.tif"))
    cand = np.argwhere((strm > 0) & (acc > 3000))
    cand = cand[(cand[:, 0] > 20) & (cand[:, 1] > 20) & (cand[:, 0] < info.rows - 20)
                & (cand[:, 1] < info.cols - 20)]
    pick = cand[np.random.default_rng(5).choice(len(cand), 3, replace=False)]
    # field declared in the URI: no QVariant/QMetaType enum (Qt5/Qt6 neutral)
    pts = QgsVectorLayer(f"Point?crs={QgsRasterLayer(dem).crs().authid()}&field=culvert:string",
                         "pts", "memory")
    feats = []
    for k, (rr, cc) in enumerate(pick):
        x, y = info.rowcol_to_xy(int(rr), int(cc))
        f = QgsFeature(pts.fields())
        f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x + 25.0, y - 20.0)))  # off-channel
        f["culvert"] = f"CV{10 + k}"
        feats.append(f)
    pts.dataProvider().addFeatures(feats)
    opts = QgsVectorFileWriter.SaveVectorOptions(); opts.driverName = "GPKG"
    QgsVectorFileWriter.writeAsVectorFormatV3(pts, out("pts.gpkg"),
                                              QgsCoordinateTransformContext(), opts)
    ptsfile = out("pts.gpkg")

    r = processing.run("qeht:delineatecatchment", {"FDR": out("fdr.tif"), "FAC": out("fac.tif"),
                       "POINTS": ptsfile, "SNAP": 5, "SNAP_THRESHOLD": 200, "ID_FIELD": "culvert",
                       "RASTER_OUT": out("cat.tif"), "POLY_OUT": out("cat.gpkg")})
    lyr = QgsVectorLayer(r["POLY_OUT"], "c", "ogr")
    uids = sorted(f["outlet_uid"] for f in lyr.getFeatures())
    check("delineate: polygons carry outlet_uid from the ID field",
          uids == ["CV10", "CV11", "CV12"], str(uids))

    r = processing.run("qeht:longestflowpath", {"FDR": out("fdr.tif"), "DEM": out("fill.tif"),
                       "FAC": out("fac.tif"), "POINTS": ptsfile, "SNAP": 5,
                       "SNAP_THRESHOLD": 200, "NESTED": False, "ID_PREFIX": "KE-",
                       "OUTPUT": out("lfp.gpkg")})
    lyr = QgsVectorLayer(r["OUTPUT"], "l", "ogr")
    uids = sorted(f["outlet_uid"] for f in lyr.getFeatures())
    check("longest flow path: sequential outlet_uid with prefix", uids == ["KE-001", "KE-002", "KE-003"],
          str(uids))

    r = processing.run("qeht:catchmentcharacteristics", {
        "FDR": out("fdr.tif"), "DEM": out("fill.tif"), "RAW_DEM": dem, "FAC": out("fac.tif"),
        "POINTS": ptsfile, "SNAP": 5, "SNAP_THRESHOLD": 200, "NESTED": False,
        "CATCH_OUT": out("ch_c.gpkg"), "PATH_OUT": out("ch_p.gpkg")})
    cl = QgsVectorLayer(r["CATCH_OUT"], "c", "ogr"); pl = QgsVectorLayer(r["PATH_OUT"], "p", "ogr")
    names = [f.name() for f in cl.fields()]
    check("characteristics: new + legacy slope fields present",
          all(n in names for n in ("outlet_uid", "catch_slope_horn", "catch_relief_ratio",
                                   "slope_mean", "slope_relief_ratio")), str(names))
    pn = [f.name() for f in pl.fields()]
    check("characteristics: lfp_z10_m / lfp_z85_m on flow paths",
          "lfp_z10_m" in pn and "lfp_L85_m" in pn)
    ok = all(abs(f.geometry().area() - f["area_km2"] * 1e6) < 1e-3 * f["area_km2"] * 1e6
             for f in cl.getFeatures())
    check("characteristics: polygon area == area_km2 for every catchment", ok)
    check("characteristics: catchment and path uids match",
          sorted(f["outlet_uid"] for f in cl.getFeatures())
          == sorted(f["outlet_uid"] for f in pl.getFeatures()))

    r = processing.run("qeht:buildheasexchange", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "RAW_DEM": dem, "ORDER": out("ord.tif"),
        "POINTS": ptsfile, "ID_FIELD": "culvert", "ID_PREFIX": "", "ID_ORDER": 0, "SNAP": 5,
        "SNAP_THRESHOLD": 200, "LOCAL": False, "DEM_SOURCE": "example_dem.tif",
        "CONDITIONING": "QEHT fill, min_slope 0", "FLAT_METHOD": "toward", "CSV": True,
        "OUTPUT": out("smoke_qeht_exchange.gpkg")})
    xp = r["OUTPUT"]
    errors, warnings = validate_exchange(xp)
    check("exchange: package validates", not errors and not warnings, "; ".join(errors + warnings))
    for layer in ("crossings", "catchments", "flowpaths"):
        ql = QgsVectorLayer(f"{xp}|layername={layer}", layer, "ogr")
        check(f"exchange: '{layer}' opens in QGIS with 3 features and EPSG:21037",
              ql.isValid() and ql.featureCount() == 3 and ql.crs().authid() == "EPSG:21037",
              f"{ql.featureCount()} features, {ql.crs().authid()}")
    md = {row["key"]: row["value"] for row in gpkg.read_table(xp, "qeht_run_metadata")}
    check("exchange: metadata records DEM hash, flat method and id attribute",
          md["dem_sha256"].startswith("sha256:") and md["flat_method"] == "toward"
          and md["id_attribute"] == "culvert")
    check("exchange: CSV sidecars written",
          os.path.exists(os.path.splitext(xp)[0] + "_csv"))
    try:
        processing.run("qeht:buildheasexchange", {
            "FDR": geo_dem, "FAC": geo_dem, "RAW_DEM": geo_dem, "POINTS": ptsfile,
            "OUTPUT": out("geo.gpkg")})
        refused = False
    except QgsProcessingException as e:
        refused = "geographic" in str(e).lower()
    check("exchange: geographic DEM refused with a clear message", refused)

    # ---- v0.10: road crossings workflow --------------------------------
    import shutil
    from osgeo import ogr
    # a road across the middle of the DEM, from west to east
    y_mid = info.geotransform[3] + info.geotransform[5] * info.rows * 0.5 + 7.0
    x0 = info.geotransform[0] + 5 * info.geotransform[1]
    x1 = info.geotransform[0] + (info.cols - 5) * info.geotransform[1]
    road = QgsVectorLayer(f"LineString?crs={QgsRasterLayer(dem).crs().authid()}", "road", "memory")
    rf = QgsFeature(); rf.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(x0, y_mid),
                                                                  QgsPointXY(x1, y_mid)]))
    road.dataProvider().addFeatures([rf])
    QgsVectorFileWriter.writeAsVectorFormatV3(road, out("road.gpkg"),
                                              QgsCoordinateTransformContext(), opts)
    r = processing.run("qeht:crossingcandidates", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "STREAMS": out("str.tif"),
        "ORDER": out("ord.tif"), "ROAD": out("road.gpkg"), "START": 1000.0, "REVERSE": False,
        "MIN_AREA": 0.05, "HALFWIDTH": 30.0, "MIN_PARALLEL": 100.0, "MERGE": 0.0,
        "PREFIX": "C", "CANDIDATES": out("cand.gpkg"), "PARALLEL": out("par.gpkg")})
    cl = QgsVectorLayer(r["CANDIDATES"], "cand", "ogr")
    feats = list(cl.getFeatures())
    rec = [f for f in feats if f["recommended"] == 1]
    chs = [f["chainage_m"] for f in feats]
    check("crossing candidates: found, chainage from 1000 m, sorted, >= 1 recommended",
          len(feats) >= 1 and len(rec) >= 1 and min(chs) >= 1000.0 and chs == sorted(chs),
          f"{len(feats)} candidates, {len(rec)} recommended")
    check("crossing candidates: status 'candidate', angle and area filled",
          all(f["status"] == "candidate" and 0 <= f["crossing_angle_deg"] <= 90
              and f["acc_km2"] >= 0.05 for f in feats))

    r = processing.run("qeht:buildheasexchange", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "RAW_DEM": dem, "ORDER": out("ord.tif"),
        "POINTS": out("cand.gpkg"), "ROAD": out("road.gpkg"), "START": 1000.0,
        "ID_PREFIX": "X", "ID_ORDER": 0, "SNAP": 5, "SNAP_THRESHOLD": 200, "LOCAL": False,
        "OUTPUT": out("cand_exchange.gpkg")})
    xp2 = r["OUTPUT"]
    errors, warnings = validate_exchange(xp2)
    tabs = dict(gpkg.list_tables(xp2))
    cr = gpkg.read_table(xp2, "crossings")
    check("exchange from candidates: validates, carries candidates + road layers",
          not errors and "crossing_candidates" in tabs and "road_alignment" in tabs,
          "; ".join(errors))
    check("exchange from candidates: recommended crossings, numbered along chainage",
          len(cr) == len(rec) and [c["outlet_uid"] for c in sorted(cr, key=lambda c: c["chainage_m"])]
          == [f"X{k + 1:03d}" for k in range(len(cr))] and all(c["snap_dist_m"] < 45 for c in cr))

    r = processing.run("qeht:burncrossings", {
        "DEM": dem, "CROSSINGS": out("cand.gpkg"), "HALF": 30.0, "SEARCH": 2,
        "OUTPUT": out("burned.tif"), "LOG": out("breach.gpkg")})
    bl = QgsVectorLayer(r["LOG"], "b", "ogr")
    check("burn crossings: burned DEM + one breach per recommended crossing",
          os.path.exists(r["OUTPUT"]) and bl.featureCount() == len(rec),
          f"{bl.featureCount()} breaches")

    # edit the package: delete the first crossing, add one on a stream elsewhere
    edited = out("edited.gpkg"); shutil.copy(xp2, edited)
    ds = ogr.Open(edited, 1); L = ds.GetLayerByName("crossings")
    first = L.GetNextFeature(); del_uid = first.GetField("outlet_uid"); L.DeleteFeature(first.GetFID())
    rr, cc = [int(v) for v in pick[0]]
    nx, ny = info.rowcol_to_xy(rr, cc)
    nf = ogr.Feature(L.GetLayerDefn()); g = ogr.Geometry(ogr.wkbPoint); g.AddPoint_2D(nx + 10, ny)
    nf.SetGeometry(g); L.CreateFeature(nf); ds = None
    r = processing.run("qeht:renumberrelink", {
        "PACKAGE": edited, "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "RAW_DEM": dem,
        "ROAD": out("road.gpkg"), "START": 1000.0, "PREFIX": "", "SNAP": 5,
        "SNAP_THRESHOLD": 200, "OUTPUT": out("relinked.gpkg")})
    rel = r["OUTPUT"]
    errors, _ = validate_exchange(rel)
    log = gpkg.read_table(rel, "renumber_log")
    changes = sorted(x["change"] for x in log)
    new_ids = sorted(x["new_uid"] for x in log if x["new_uid"])
    check("renumber and relink: validates, logs deleted + new, gapless ids",
          not errors and "deleted" in changes and "new" in changes
          and new_ids == [f"X{k + 1:03d}" for k in range(len(new_ids))]
          and any(x["old_uid"] == del_uid and x["change"] == "deleted" for x in log),
          f"{changes}")

    # ---- v0.11: soils --------------------------------------------------
    gt = info.geotransform
    xmid = gt[0] + gt[1] * info.cols * 0.5
    xa, xb = gt[0] - 100, gt[0] + gt[1] * info.cols + 100
    ya, yb = gt[3] + gt[5] * info.rows - 100, gt[3] + 100
    soil = QgsVectorLayer(f"Polygon?crs={QgsRasterLayer(dem).crs().authid()}"
                          "&field=sand:double&field=silt:double&field=clay:double"
                          "&field=oc:double&field=drain:string", "soil", "memory")
    sf = []
    for (x0, x1), vals in (((xa, xmid), [35.0, 30.0, 35.0, 1.8, "W"]),
                           ((xmid, xb), [60.0, 20.0, 20.0, 0.9, "M"])):
        f = QgsFeature(soil.fields())
        f.setGeometry(QgsGeometry.fromPolygonXY([[QgsPointXY(x0, ya), QgsPointXY(x1, ya),
                                                  QgsPointXY(x1, yb), QgsPointXY(x0, yb)]]))
        f.setAttributes(vals); sf.append(f)
    soil.dataProvider().addFeatures(sf)
    QgsVectorFileWriter.writeAsVectorFormatV3(soil, out("soil.gpkg"),
                                              QgsCoordinateTransformContext(), opts)
    r = processing.run("qeht:soilparameters", {
        "CATCHMENTS": out("ch_c.gpkg"), "REF": out("fdr.tif"), "SOIL_POLYGONS": out("soil.gpkg"),
        "SOIL_TOP": 0.0, "SOIL_BOTTOM": 20.0, "OUTPUT": out("ch_soil.gpkg"),
        "K_RASTER": out("k.tif")})
    sl = QgsVectorLayer(r["OUTPUT"], "s", "ogr")
    fl = list(sl.getFeatures())
    check("soil parameters: every catchment gets K, texture, full coverage; K raster written",
          len(fl) == 3 and all(f["usle_k"] and 0.005 < f["usle_k"] < 0.07 and f["soil_texture"]
                               and abs(f["soil_coverage_pct"] - 100) < 1e-6 for f in fl)
          and os.path.exists(r["K_RASTER"]) and "outlet_uid" in sl.fields().names(),
          ", ".join(f"{f['outlet_uid']} K={f['usle_k']:.4f} {f['soil_texture']}" for f in fl))
    r = processing.run("qeht:buildheasexchange", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "RAW_DEM": dem, "POINTS": ptsfile,
        "ID_FIELD": "culvert", "SNAP": 5, "SNAP_THRESHOLD": 200,
        "SOIL_POLYGONS": out("soil.gpkg"), "OUTPUT": out("soil_exchange.gpkg")})
    ca = gpkg.read_table(r["OUTPUT"], "catchments", with_geometry=False)
    md = {row["key"]: row["value"] for row in gpkg.read_table(r["OUTPUT"], "qeht_run_metadata")}
    errors, _ = validate_exchange(r["OUTPUT"])
    check("exchange with soils: soil block on catchments, dataset in metadata, validates",
          not errors and all(c["usle_k"] and c["soil_hsg_proxy"] for c in ca)
          and md["soil_depth_cm"] and md["soil_dataset"])

    print("\n" + ("ALL QGIS SMOKE CHECKS PASSED" if not FAILURES
                  else f"{len(FAILURES)} FAILURE(S): " + "; ".join(FAILURES)))
    print("outputs in", tmp)
    if app is not None:
        app.exitQgis()
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
