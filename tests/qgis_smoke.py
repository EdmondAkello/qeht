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
    check("provider loads with 16 algorithms incl. exchange, crossings, burn, relink, soils, "
          "erosion, alignment profile",
          len(algs) == 16 and all(a in algs for a in ("qeht:alignmentprofile",
              "qeht:buildheasexchange", "qeht:crossingcandidates", "qeht:burncrossings",
              "qeht:renumberrelink", "qeht:soilparameters", "qeht:erosionindices",
              "qeht:erosioncorridor")))

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
    for m, extra in ((1, {}),):
        r = processing.run("qeht:flowdirection", dict({"DEM": out("fill.tif"), "RESOLVE_FLATS": True,
                           "FLAT_METHOD": m, "OUTPUT": out(f"fdr_m{m}.tif")}, **extra))
        from ..core.grid import decode_d8
        from ..core.flow.accumulation import flow_accumulation as _fa
        from ..core.raster import read_dem as _rd
        _d8, _v, _ = _rd(r["OUTPUT"])
        _, _st = _fa(decode_d8(_d8.astype(int)), _v)
        check(f"flow direction ({['toward', 'Barnes'][m]}): no flow cycles",
              _st["cells_in_cycles"] == 0, f"{_st['cells_in_cycles']} cells in cycles")
    import numpy as _np
    fdr_default = processing.run("qeht:flowdirection", {"DEM": out("fill.tif"),
                                 "OUTPUT": out("fdr_default.tif")})
    check("flow direction: Barnes is the default (0.13)",
          _np.array_equal(_rd(fdr_default["OUTPUT"])[0], _rd(out("fdr_m1.tif"))[0]))
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
    check("characteristics: flat-method check fields (0.13.1), areas consistent",
          all(n in names for n in ("area_barnes_km2", "area_toward_km2", "flat_sensitivity_pct",
                                   "flat_sensitive"))
          and all(f["flat_sensitive"] in (0, 1) and f["area_barnes_km2"] > 0 for f in cl.getFeatures()),
          ", ".join(f"{f['outlet_uid']} {f['area_barnes_km2']:.3f}/{f['area_toward_km2']:.3f}"
                    for f in cl.getFeatures()))
    r2 = processing.run("qeht:catchmentcharacteristics", {
        "FDR": out("fdr.tif"), "DEM": out("fill.tif"), "FAC": out("fac.tif"), "POINTS": ptsfile,
        "ID_FIELD": "culvert", "ID_PREFIX": "X", "FLAT_CHECK": False,
        "CATCH_OUT": out("ch_c2.gpkg"), "PATH_OUT": out("ch_p2.gpkg")})
    c2 = QgsVectorLayer(r2["CATCH_OUT"], "c", "ogr")
    check("characteristics: ID prefix not added to IDs from a field; check can be switched off",
          sorted(f["outlet_uid"] for f in c2.getFeatures()) == ["CV10", "CV11", "CV12"]
          and "flat_sensitive" not in [f.name() for f in c2.fields()],
          str(sorted(f["outlet_uid"] for f in c2.getFeatures())))
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

    r = processing.run("qeht:alignmentprofile", {
        "RAW_DEM": dem, "ROAD": out("road.gpkg"), "FILLED": out("fill.tif"), "FDR": out("fdr.tif"),
        "FAC": out("fac.tif"), "ORDER": out("ord.tif"), "THRESHOLD": 200, "START": 1000.0,
        "STEP": 10.0, "PROFILE": out("profile.gpkg"), "CSV": out("profile.csv"),
        "CHART": out("profile.png")})
    pl = QgsVectorLayer(r["PROFILE"], "profile", "ogr")
    pf = list(pl.getFeatures())
    nst = sum(1 for f in pf if f["stream"] == 1)
    check("alignment profile: stations from 1000 m every 10 m, ground and fill filled, "
          "stream crossings flagged, CSV written",
          len(pf) > 2 and pf[0]["chainage_m"] == 1000.0 and abs(pf[1]["chainage_m"] - 1010.0) < 1e-9
          and all(f["z_dem_m"] is not None and f["z_fill_m"] >= f["z_dem_m"] - 1e-6 for f in pf)
          and nst >= 1 and os.path.exists(out("profile.csv")),
          f"{len(pf)} stations, {nst} stream crossing(s), chart {os.path.exists(out('profile.png'))}")

    r = processing.run("qeht:buildheasexchange", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "RAW_DEM": dem, "ORDER": out("ord.tif"),
        "FILLED": out("fill.tif"),
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
    ap = gpkg.read_table(xp2, "alignment_profile", with_geometry=False) \
        if "alignment_profile" in tabs else []
    mdp = {row["key"]: row["value"] for row in gpkg.read_table(xp2, "qeht_run_metadata")}
    check("package carries alignment_profile (same stations as the tool) + metadata",
          len(ap) == len(pf) and sum(a["stream"] for a in ap) == nst
          and mdp.get("alignment_step_m") == "10" and mdp.get("alignment_source"),
          f"{len(ap)} stations")
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
    r = processing.run("qeht:buildheasexchange", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "RAW_DEM": dem, "POINTS": out("cand.gpkg"),
        "BURN_LOG": out("breach.gpkg"), "CONDITIONING": "fill, min_slope 0",
        "OUTPUT": out("burn_exchange.gpkg")})
    md = {row["key"]: row["value"] for row in gpkg.read_table(r["OUTPUT"], "qeht_run_metadata")}
    errors, _ = validate_exchange(r["OUTPUT"])
    bt = gpkg.read_table(r["OUTPUT"], "burn_log", with_geometry=False)
    check("exchange records the breach log: burn_log layer + conditioning metadata",
          not errors and len(bt) == bl.featureCount() and "burn crossings:" in md["conditioning"]
          and md["conditioning"].startswith("fill, min_slope 0") and md["conditioning_burn"],
          md["conditioning"])

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
    # same soil as polygons with only a unit code + a CSV table (v0.14)
    soil2 = QgsVectorLayer(f"Polygon?crs={QgsRasterLayer(dem).crs().authid()}"
                           "&field=SU:string", "soil2", "memory")
    sf2 = []
    for (x0, x1), code in (((xa, xmid), "A1"), ((xmid, xb), "B2")):
        f = QgsFeature(soil2.fields())
        f.setGeometry(QgsGeometry.fromPolygonXY([[QgsPointXY(x0, ya), QgsPointXY(x1, ya),
                                                  QgsPointXY(x1, yb), QgsPointXY(x0, yb)]]))
        f.setAttributes([code]); sf2.append(f)
    soil2.dataProvider().addFeatures(sf2)
    QgsVectorFileWriter.writeAsVectorFormatV3(soil2, out("soil_units.gpkg"),
                                              QgsCoordinateTransformContext(), opts)
    with open(out("soil_table.csv"), "w", encoding="utf-8") as fh:
        fh.write("SU,sand,silt,clay,oc,drain\nA1,35,30,35,1.8,W\nB2,60,20,20,0.9,M\n")
    r2 = processing.run("qeht:soilparameters", {
        "CATCHMENTS": out("ch_c.gpkg"), "REF": out("fdr.tif"), "SOIL_POLYGONS": out("soil_units.gpkg"),
        "SOIL_CSV": out("soil_table.csv"), "SOIL_UNIT_FIELD": "SU", "OUTPUT": out("ch_soil_csv.gpkg")})
    k1 = sorted(round(f["usle_k"], 6) for f in fl)
    k2 = sorted(round(f["usle_k"], 6) for f in QgsVectorLayer(r2["OUTPUT"], "s2", "ogr").getFeatures())
    check("soil parameters from unit polygons + CSV table = same K as attribute polygons",
          k1 == k2, f"{k1} vs {k2}")
    # v0.15 (G2): any polygon map - fields chosen by name, OC in g/kg, silt derived, HSG field
    soil3 = QgsVectorLayer(f"Polygon?crs={QgsRasterLayer(dem).crs().authid()}"
                           "&field=SAND_P:double&field=CLAY_P:double&field=SOC_GKG:double"
                           "&field=HYD_GRP:string", "soil3", "memory")
    sf3 = []
    for (x0, x1), vals in (((xa, xmid), [35.0, 35.0, 18.0, "C"]),
                           ((xmid, xb), [60.0, 20.0, 9.0, "B/D"])):
        f = QgsFeature(soil3.fields())
        f.setGeometry(QgsGeometry.fromPolygonXY([[QgsPointXY(x0, ya), QgsPointXY(x1, ya),
                                                  QgsPointXY(x1, yb), QgsPointXY(x0, yb)]]))
        f.setAttributes(vals); sf3.append(f)
    soil3.dataProvider().addFeatures(sf3)
    QgsVectorFileWriter.writeAsVectorFormatV3(soil3, out("soil_named.gpkg"),
                                              QgsCoordinateTransformContext(), opts)
    r3 = processing.run("qeht:soilparameters", {
        "CATCHMENTS": out("ch_c.gpkg"), "REF": out("fdr.tif"), "SOIL_POLYGONS": out("soil_named.gpkg"),
        "SOIL_F_SAND": "SAND_P", "SOIL_F_CLAY": "CLAY_P", "SOIL_F_OC": "SOC_GKG",
        "SOIL_OC_UNITS": 1, "SOIL_F_HSG_FIELD": "HYD_GRP", "OUTPUT": out("ch_soil_named.gpkg")})
    f3 = list(QgsVectorLayer(r3["OUTPUT"], "s3", "ogr").getFeatures())
    k3 = sorted(round(f["usle_k"], 6) for f in f3)
    check("soils from any polygon map (fields by name, OC g/kg, silt derived) = same K; "
          "HSG field used",
          k3 == k1 and all(f["soil_hsg"] in ("C", "D") and abs(f["hsg_pct_c"] + f["hsg_pct_d"] - 100) < 1e-6
                           and "HYD_GRP" in f["soil_hsg_source"] for f in f3),
          ", ".join(f"{f['outlet_uid']} HSG {f['soil_hsg']} (C {f['hsg_pct_c']:.0f} / D {f['hsg_pct_d']:.0f}, "
                    f"dual {f['hsg_pct_dual']:.0f})" for f in f3))
    # v0.15 (G3/G4): SoilGrids folder (depth layers) + HYSOGs-coded HSG raster on the DEM grid
    from osgeo import gdal as _g
    ds0 = _g.Open(dem)
    gt0, prj0 = ds0.GetGeoTransform(), ds0.GetProjection()
    ny, nx = ds0.RasterYSize, ds0.RasterXSize
    ds0 = None
    xs = gt0[0] + (np.arange(nx) + 0.5) * gt0[1]
    left = np.tile(xs < xmid, (ny, 1))
    sgdir = out("soilgrids"); os.makedirs(sgdir, exist_ok=True)

    def _wr(path, arr, dtype=_g.GDT_Float32, nodata=-32768):
        d = _g.GetDriverByName("GTiff").Create(path, nx, ny, 1, dtype)
        d.SetGeoTransform(gt0); d.SetProjection(prj0)
        b = d.GetRasterBand(1); b.SetNoDataValue(nodata); b.WriteArray(arr); d = None
    for depth in ("0-5", "5-15", "15-30"):
        for prop, lv, rv in (("sand", 350, 600), ("clay", 350, 200), ("silt", 300, 200),
                             ("soc", 180, 90)):
            _wr(os.path.join(sgdir, f"{prop}_{depth}cm_mean.tif"), np.where(left, lv, rv).astype(np.float32))
    _wr(out("hysogs.tif"), np.where(left, 2, 14).astype(np.uint8), _g.GDT_Byte, 255)
    r4 = processing.run("qeht:soilparameters", {
        "CATCHMENTS": out("ch_c.gpkg"), "REF": out("fdr.tif"), "SOIL_SOILGRIDS": sgdir,
        "SOIL_HSG_R": out("hysogs.tif"), "OUTPUT": out("ch_soil_sg.gpkg"), "K_RASTER": out("k_sg.tif")})
    f4 = list(QgsVectorLayer(r4["OUTPUT"], "s4", "ogr").getFeatures())
    k4 = sorted(round(f["usle_k"], 6) for f in f4)
    check("soils from a SoilGrids folder (3 depth layers, g/kg -> %) = same K; HYSOGs raster "
          "gives B and dual D shares; K raster written",
          max(abs(a - b) for a, b in zip(k4, k1)) < 1e-6 and os.path.exists(out("k_sg.tif"))
          and all(abs(f["hsg_pct_b"] + f["hsg_pct_d"] - 100) < 1e-6
                  and abs(f["hsg_pct_dual"] - f["hsg_pct_d"]) < 1e-6
                  and "SoilGrids" in f["soil_dataset"] for f in f4),
          f"{k4} vs {k1}")
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
    # v0.15 (F1): curve number from WorldCover classes x HSG
    _wr(out("worldcover.tif"), np.where(left, 30, 40).astype(np.uint8), _g.GDT_Byte, 0)
    r5 = processing.run("qeht:soilparameters", {
        "CATCHMENTS": out("ch_c.gpkg"), "REF": out("fdr.tif"), "SOIL_HSG_R": out("hysogs.tif"),
        "LANDCOVER": out("worldcover.tif"), "CN_CONDITION": 0, "CN_AMC": 1,
        "OUTPUT": out("ch_cn.gpkg")})
    f5 = list(QgsVectorLayer(r5["OUTPUT"], "s5", "ogr").getFeatures())
    # left: grassland x B = 69, right: cropland x B/D (as D) = 89 (TR-55 fair)
    check("curve numbers: grass-B 69 / crops-D 89 mix, AMC III export, shares add up",
          len(f5) == 3 and all(69 - 1e-6 <= f["cn_ii"] <= 89 + 1e-6 and f["cn_amc"] == "III"
                               and f["cn_export"] > f["cn_ii"]
                               and abs(f["lc_pct_grass"] + f["lc_pct_crop"] - 100) < 1e-6
                               and abs(f["cn_ii"] - (69 * f["lc_pct_grass"] + 89 * f["lc_pct_crop"]) / 100) < 1e-6
                               for f in f5),
          ", ".join(f"{f['outlet_uid']} CN {f['cn_ii']:.1f} -> {f['cn_export']:.1f}" for f in f5))
    r = processing.run("qeht:buildheasexchange", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "RAW_DEM": dem, "POINTS": ptsfile,
        "ID_FIELD": "culvert", "SNAP": 5, "SNAP_THRESHOLD": 200, "SOIL_SOILGRIDS": sgdir,
        "SOIL_HSG_R": out("hysogs.tif"), "LANDCOVER": out("worldcover.tif"),
        "OUTPUT": out("soil_sg_exchange.gpkg")})
    ca = gpkg.read_table(r["OUTPUT"], "catchments", with_geometry=False)
    md = {row["key"]: row["value"] for row in gpkg.read_table(r["OUTPUT"], "qeht_run_metadata")}
    errors, _ = validate_exchange(r["OUTPUT"])
    check("package with SoilGrids + HYSOGs + WorldCover: soil block, HSG, CN, sources in metadata",
          not errors and all(c["usle_k"] and c["soil_hsg"] for c in ca)
          and "SoilGrids" in md["soil_dataset"] and "hysogs" in md["soil_hsg_source"].lower()
          and all(c["cn_export"] for c in ca) and '"cn_proxy": true' in md["runoff_json"])

    # ---- v0.14: erosion (WP-F) -------------------------------------------
    import json as _json
    from osgeo import gdal as _gdal
    r = processing.run("qeht:erosionindices", {
        "RAW_DEM": dem, "FAC": out("fac.tif"), "LS_METHOD": 0, "CHANNEL_CELLS": 500,
        "R_VALUE": 1800.0, "C_VALUE": 0.05, "SOIL_POLYGONS": out("soil.gpkg"),
        "SPI_SCHEME": 1, "RUSLE_SCHEME": 1, "OUTPUT": out("erosion")})
    ef = r["OUTPUT"]
    with open(os.path.join(ef, "erosion_run.json")) as fh:
        run = _json.load(fh)
    need = ["ln_spi.tif", "ls_factor.tif", "twi.tif", "rusle_soil_loss_t_ha_yr.tif",
            "spi_class.tif", "rusle_class.tif", "erosion_combined_class.tif", "class_extents.csv",
            "erosion_combined_class.qml", "k_factor.tif"]
    check("erosion indices + RUSLE: rasters, class rasters, styles, run record",
          all(os.path.exists(os.path.join(ef, n)) for n in need) and run["mode"] == "RUSLE"
          and any(f["factor"] == "K" and "Williams" in f["source"] for f in run["factors"]),
          str([n for n in need if not os.path.exists(os.path.join(ef, n))]))
    ds = _gdal.Open(os.path.join(ef, "spi_class.tif"))
    rat = ds.GetRasterBand(1).GetDefaultRAT()
    check("SPI class raster carries the A14 class names (attribute table)",
          rat is not None and [rat.GetValueAsString(i, 1) for i in range(rat.GetRowCount())]
          == ["Low", "Moderate", "High", "Severe"])
    ds = None
    r2 = processing.run("qeht:erosionindices", {"RAW_DEM": dem, "FAC": out("fac.tif"),
                                                "OUTPUT": out("erosion_ls")})
    with open(os.path.join(r2["OUTPUT"], "erosion_run.json")) as fh:
        run2 = _json.load(fh)
    check("no R/K/C: LS-only terrain potential, no soil-loss raster",
          run2["mode"].startswith("LS-only")
          and os.path.exists(os.path.join(r2["OUTPUT"], "ls_class.tif"))
          and not os.path.exists(os.path.join(r2["OUTPUT"], "rusle_soil_loss_t_ha_yr.tif")))
    r = processing.run("qeht:erosioncorridor", {
        "EROSION": ef, "ROAD": out("road.gpkg"), "START": 1000.0, "STEP": 10, "HALF_WIDTH": 50,
        "STATIONS": out("ero_st.gpkg"), "REACHES": out("ero_re.gpkg"), "CHART": out("ero.png")})
    st = QgsVectorLayer(r["STATIONS"], "s", "ogr"); rl = QgsVectorLayer(r["REACHES"], "r", "ogr")
    chs = [f["chainage"] for f in st.getFeatures()]
    check("corridor: stations from ch 1000 every 10 m, reaches cover the road",
          chs[0] == 1000.0 and all(b > a for a, b in zip(chs, chs[1:])) and rl.featureCount() >= 1
          and abs(sum(f["length_m"] for f in rl.getFeatures()) - (chs[-1] - chs[0])) < 1e-6,
          f"{st.featureCount()} stations, {rl.featureCount()} reaches, chart {os.path.exists(out('ero.png'))}")
    r = processing.run("qeht:buildheasexchange", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "RAW_DEM": dem, "POINTS": ptsfile,
        "ID_FIELD": "culvert", "SNAP": 5, "SNAP_THRESHOLD": 200, "SOIL_POLYGONS": out("soil.gpkg"),
        "EROSION": ef, "OUTPUT": out("erosion_exchange.gpkg")})
    ca = gpkg.read_table(r["OUTPUT"], "catchments", with_geometry=False)
    cr = gpkg.read_table(r["OUTPUT"], "crossings", with_geometry=False)
    md = {row["key"]: row["value"] for row in gpkg.read_table(r["OUTPUT"], "qeht_run_metadata")}
    errors, _ = validate_exchange(r["OUTPUT"])
    check("exchange with erosion: ero_* on catchments and crossings, composite, metadata, validates",
          not errors and all(c["ero_a_mean_tha"] is not None and c["ero_sy_m3_yr"] is not None
                             for c in ca)
          and all(x["ero_composite"] is not None and x["ero_impact"] for x in cr)
          and "mcdma_weights" in md["erosion_json"],
          ", ".join(f"{x['outlet_uid']} {x['ero_impact']} ({x['ero_composite']:.1f})" for x in cr))
    check("v0.15: STI raster written; every crossing has channel slopes and a deposition "
          "indicator (or a note); STI/deposition settings in metadata",
          os.path.exists(os.path.join(ef, "sti_overland.tif"))
          and all(x["ch_slope_us"] is not None and x["ch_slope_ds"] is not None
                  and (x["ero_dep_flag"] in ("deposition-prone", "neutral", "scour-prone")
                       or x["ero_dep_note"]) for x in cr)
          and '"deposition"' in md["erosion_json"] and md["channel_slope_m"] == "200",
          ", ".join(f"{x['outlet_uid']} us {x['ch_slope_us']:.3f} ds {x['ch_slope_ds']:.3f} "
                    f"{x['ero_dep_flag'] or x['ero_dep_note']}" for x in cr))

    print("\n" + ("ALL QGIS SMOKE CHECKS PASSED" if not FAILURES
                  else f"{len(FAILURES)} FAILURE(S): " + "; ".join(FAILURES)))
    print("outputs in", tmp)
    if app is not None:
        app.exitQgis()
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
