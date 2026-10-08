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

import numpy as np

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
    check("provider loads with 23 algorithms incl. exchange, crossings, burn, relink, soils, "
          "erosion, alignment profile, pipeline, run report, prepare DEM",
          len(algs) == 23 and all(a in algs for a in ("qeht:alignmentprofile", "qeht:drainagecoverage",
              "qeht:preparedem", "qeht:exportpackage", "qeht:autoclip", "qeht:demuncertainty",
              "qeht:hydrologypipeline", "qeht:runreport",
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
        check(f"exchange: '{layer}' opens in QGIS with 3 features in the DEM's CRS",
              ql.isValid() and ql.featureCount() == 3
              and ql.crs().authid() == QgsRasterLayer(dem).crs().authid(),
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
    crx = [c for c in cr if c["status"] == "existing"]
    crp = [c for c in cr if c["status"] == "proposed"]
    check("exchange from candidates: recommended crossings, numbered along chainage",
          len(crx) == len(rec) and [c["outlet_uid"] for c in sorted(crx, key=lambda c: c["chainage_m"])]
          == [f"X{k + 1:03d}" for k in range(len(crx))] and all(c["snap_dist_m"] < 45 for c in crx))
    cov = gpkg.read_table(xp2, "coverage_check", with_geometry=False) if "coverage_check" in tabs else []
    cat_p = [c for c in gpkg.read_table(xp2, "catchments", with_geometry=False) if c["status"] == "proposed"]
    check("v0.16 coverage check in the package: coverage_check / sag_points / flat_stretches "
          "layers; proposed crossings P001.. delineated with catchments",
          all(t in tabs for t in ("coverage_check", "sag_points", "flat_stretches"))
          and len(cat_p) == len(crp)
          and sorted(c["outlet_uid"] for c in crp) == [f"P{k + 1:03d}" for k in range(len(crp))]
          and all(c["proposed_reason"] in ("uncovered_stream", "sag_point") for c in crp)
          and mdp.get("n_proposed") == str(len(crp)),
          f"{len(crp)} proposed ({sorted(set(c['proposed_reason'] for c in crp))}), "
          f"{len(cov)} findings: {sorted(set(c['issue'] for c in cov))}")
    r = processing.run("qeht:drainagecoverage", {
        "RAW_DEM": dem, "FILLED": out("fill.tif"), "FDR": out("fdr.tif"), "FAC": out("fac.tif"),
        "ROAD": out("road.gpkg"), "CROSSINGS": out("cand.gpkg"), "START": 1000.0,
        "COVERAGE_OUT": out("cov.gpkg"), "SAGS_OUT": out("sags.gpkg"), "FLATS_OUT": out("flats.gpkg")})
    cv = list(QgsVectorLayer(r["COVERAGE_OUT"], "cv", "ogr").getFeatures())
    check("standalone coverage tool: same missing crossings as the package",
          sum(1 for f in cv if f["issue"] == "missing_crossing")
          == sum(1 for c in cov if c["issue"] == "missing_crossing"),
          f"{len(cv)} findings")

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
    rel_cr = gpkg.read_table(rel, "crossings", with_geometry=False)
    xs = sorted(c["outlet_uid"] for c in rel_cr if c["status"] != "proposed")
    ps = sorted(c["outlet_uid"] for c in rel_cr if c["status"] == "proposed")
    check("relink keeps status: proposed stay P (gapless), the rest X (gapless)",
          xs == [f"X{k + 1:03d}" for k in range(len(xs))]
          and ps == [f"P{k + 1:03d}" for k in range(len(ps))] and len(ps) == len(crp),
          f"{len(xs)} X, {len(ps)} P")
    new_ids = [u for u in new_ids if u.startswith("X")]
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
    # ---- v0.26 (F14): land-cover scenario (grassland on the left becomes built-up) -----
    _wr(out("worldcover_scn.tif"), np.where(left, 50, 40).astype(np.uint8), _g.GDT_Byte, 0)
    r6 = processing.run("qeht:soilparameters", {
        "CATCHMENTS": out("ch_c.gpkg"), "REF": out("fdr.tif"), "SOIL_HSG_R": out("hysogs.tif"),
        "LANDCOVER": out("worldcover.tif"), "LANDCOVER_SCN": out("worldcover_scn.tif"),
        "SCENARIO_NAME": "built-up west", "OUTPUT": out("ch_scn.gpkg")})
    f6 = list(QgsVectorLayer(r6["OUTPUT"], "s6", "ogr").getFeatures())
    check("scenario in the soil tool: grass -> built-up raises CN where there was grass (built share "
          "= former grass share), no change elsewhere",
          len(f6) == 3 and all(f["scenario_name"] == "built-up west"
                               and abs(f["lc_pct_built_scn"] - f["lc_pct_grass"]) < 1e-6
                               and (f["d_cn"] > 0) == (f["lc_pct_grass"] > 0) for f in f6),
          ", ".join(f"{f['outlet_uid']} dCN {f['d_cn']:.2f}" for f in f6))
    r = processing.run("qeht:buildheasexchange", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "RAW_DEM": dem, "POINTS": ptsfile,
        "ID_FIELD": "culvert", "SNAP": 5, "SNAP_THRESHOLD": 200, "SOIL_HSG_R": out("hysogs.tif"),
        "LANDCOVER": out("worldcover.tif"), "LANDCOVER_SCN": out("worldcover_scn.tif"),
        "SCENARIO_NAME": "built-up west", "OUTPUT": out("scn_exchange.gpkg")})
    ca6 = gpkg.read_table(r["OUTPUT"], "catchments", with_geometry=False)
    md6 = {row["key"]: row["value"] for row in gpkg.read_table(r["OUTPUT"], "qeht_run_metadata")}
    errors, _ = validate_exchange(r["OUTPUT"])
    check("scenario in the package: cn_ii_scn and d_cn on every catchment, scenario_json, validates",
          not errors and all(c["cn_ii_scn"] is not None and c["d_cn"] is not None for c in ca6)
          and "built-up west" in md6["scenario_json"])

    # ---- v0.17 (F3): rainfall zones, mean annual rainfall, R estimate ------
    zl = QgsVectorLayer(f"Polygon?crs={QgsRasterLayer(dem).crs().authid()}&field=zone:string",
                        "zones", "memory")
    zf = []
    for (x0, x1), name in (((xa, xmid), "Zone II"), ((xmid, xb), "Zone I")):
        f = QgsFeature(zl.fields())
        f.setGeometry(QgsGeometry.fromPolygonXY([[QgsPointXY(x0, ya), QgsPointXY(x1, ya),
                                                  QgsPointXY(x1, yb), QgsPointXY(x0, yb)]]))
        f.setAttributes([name]); zf.append(f)
    zl.dataProvider().addFeatures(zf)
    QgsVectorFileWriter.writeAsVectorFormatV3(zl, out("rain_zones.gpkg"),
                                              QgsCoordinateTransformContext(), opts)
    _wr(out("map.tif"), np.tile(np.linspace(600, 1400, nx), (ny, 1)).astype(np.float32))
    r6 = processing.run("qeht:soilparameters", {
        "CATCHMENTS": out("ch_c.gpkg"), "REF": out("fdr.tif"), "RAIN_ZONES": out("rain_zones.gpkg"),
        "RAIN_ZONE_FIELD": "zone", "RAIN_MAP": out("map.tif"), "RAIN_MAP_DATASET": "synthetic",
        "RAIN_R_RELATION": 1, "OUTPUT": out("ch_rain.gpkg")})
    f6 = list(QgsVectorLayer(r6["OUTPUT"], "s6", "ogr").getFeatures())
    check("rainfall only (no soil): zone, zone shares, MAP within the raster range, R estimate",
          len(f6) == 3 and all(f["rain_zone"] in ("Zone I", "Zone II") and 0 < f["rain_zone_pct"] <= 100
                               and 600 <= f["map_mm"] <= 1400 and f["rusle_r"] > 0
                               and "ESTIMATE" in f["rusle_r_method"] for f in f6),
          ", ".join(f"{f['outlet_uid']} {f['rain_zone']} {f['rain_zone_pct']:.0f}% "
                    f"{f['map_mm']:.0f} mm R {f['rusle_r']:.0f}" for f in f6))
    r = processing.run("qeht:buildheasexchange", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "RAW_DEM": dem, "POINTS": ptsfile,
        "ID_FIELD": "culvert", "SNAP": 5, "SNAP_THRESHOLD": 200,
        "RAIN_ZONES": out("rain_zones.gpkg"), "RAIN_ZONE_FIELD": "zone", "RAIN_MAP": out("map.tif"),
        "OUTPUT": out("rain_exchange.gpkg")})
    ca = gpkg.read_table(r["OUTPUT"], "catchments", with_geometry=False)
    md = {row["key"]: row["value"] for row in gpkg.read_table(r["OUTPUT"], "qeht_run_metadata")}
    errors, _ = validate_exchange(r["OUTPUT"])
    check("package with rainfall: zone + MAP on catchments, no R without a relation, "
          "rainfall_json in metadata, validates",
          not errors and all(c["rain_zone"] and c["map_mm"] and c["rusle_r"] is None for c in ca)
          and '"r_relation": "none"' in md["rainfall_json"])

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
    r3 = processing.run("qeht:erosionindices", {
        "RAW_DEM": dem, "FAC": out("fac.tif"), "R_MAP": out("map.tif"), "R_RELATION": 0,
        "K_VALUE": 0.03, "C_VALUE": 0.05, "OUTPUT": out("erosion_rmap")})
    with open(os.path.join(r3["OUTPUT"], "erosion_run.json")) as fh:
        run3 = _json.load(fh)
    rf = [f for f in run3["factors"] if f["factor"] == "R"]
    check("R from a rainfall raster: RUSLE runs, R labelled an estimate (Renard & Freimund)",
          run3["mode"] == "RUSLE" and rf and rf[0]["basis"] == "estimate from rainfall"
          and "Renard & Freimund" in rf[0]["source"], str(rf))
    r = processing.run("qeht:erosioncorridor", {
        "EROSION": ef, "ROAD": out("road.gpkg"), "START": 1000.0, "STEP": 10, "HALF_WIDTH": 50,
        "STATIONS": out("ero_st.gpkg"), "REACHES": out("ero_re.gpkg"), "CHART": out("ero.png")})
    st = QgsVectorLayer(r["STATIONS"], "s", "ogr"); rl = QgsVectorLayer(r["REACHES"], "r", "ogr")
    chs = [f["chainage"] for f in st.getFeatures()]
    check("corridor: stations from ch 1000 every 10 m, reaches cover the road",
          chs[0] == 1000.0 and all(b > a for a, b in zip(chs, chs[1:])) and rl.featureCount() >= 1
          and abs(sum(f["length_m"] for f in rl.getFeatures()) - (chs[-1] - chs[0])) < 1e-6,
          f"{st.featureCount()} stations, {rl.featureCount()} reaches, chart {os.path.exists(out('ero.png'))}")
    # ---- v0.20 (STI R3): side-drain siltation indicator --------------------
    r = processing.run("qeht:erosioncorridor", {
        "EROSION": ef, "ROAD": out("road.gpkg"), "START": 1000.0, "STEP": 10, "HALF_WIDTH": 50,
        "FDR": out("fdr.tif"), "RAW_DEM": dem,
        "STATIONS": out("sd_st.gpkg"), "REACHES": out("sd_re.gpkg"), "CHART": out("sd.png")})
    st = QgsVectorLayer(r["STATIONS"], "s", "ogr"); rl = QgsVectorLayer(r["REACHES"], "r", "ogr")
    sf, rf_ = st.fields().names(), rl.fields().names()
    need_s = ["sti_p50_lhs", "sti_p90_lhs", "n_lhs", "sti_p50_rhs", "sti_p90_rhs", "n_rhs",
              "slope_long_pct", "siltation_lhs", "siltation_rhs"]
    need_r = ["sti_p90_lhs", "sti_p90_rhs", "sti_max_lhs", "sti_max_rhs", "sti_class_lhs",
              "sti_class_rhs", "siltation_len_lhs_m", "siltation_len_rhs_m"]
    feats = list(st.getFeatures())
    differ = sum(1 for f in feats if f["n_lhs"] != f["n_rhs"])
    check("side drains: station and reach fields, LHS and RHS differ, slope joined, chart written",
          all(k in sf for k in need_s) and all(k in rf_ for k in need_r) and differ > 0
          and any(f["slope_long_pct"] is not None and str(f["slope_long_pct"]) != "NULL" for f in feats)
          and os.path.exists(out("sd.png")),
          f"{differ} of {len(feats)} stations with different LHS/RHS counts; flagged "
          f"{sum(1 for f in feats if f['siltation_lhs'] == 1)} L / "
          f"{sum(1 for f in feats if f['siltation_rhs'] == 1)} R")
    r = processing.run("qeht:buildheasexchange", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "RAW_DEM": dem, "POINTS": out("cand.gpkg"),
        "ROAD": out("road.gpkg"), "START": 1000.0, "COVERAGE": False, "EROSION": ef,
        "OUTPUT": out("sd_pkg.gpkg")})
    errors, _ = validate_exchange(r["OUTPUT"])
    cst = gpkg.read_table(r["OUTPUT"], "corridor_sti", with_geometry=False) \
        if "corridor_sti" in dict(gpkg.list_tables(r["OUTPUT"])) else []
    mds = {row["key"]: row["value"] for row in gpkg.read_table(r["OUTPUT"], "qeht_run_metadata")}
    check("package with erosion folder and road: corridor_sti reaches, settings in metadata, "
          "validates",
          not errors and cst and abs(sum(c["length_m"] for c in cst)
                                     - (cst[-1]["ch_end"] - cst[0]["ch_start"])) < 1e-6
          and any(c["sti_class_lhs"] or c["sti_class_rhs"] for c in cst)
          and "class_breaks_p50_p75_p90" in mds["corridor_sti_params_json"],
          f"{len(cst)} reaches; " + "; ".join(errors))
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

    # ---- v0.15.1: real-data findings (Site C corridor) ----------------------
    import json as _json2
    r = processing.run("qeht:crossingcandidates", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "STREAMS": out("str.tif"),
        "ROAD": out("road.gpkg"), "START": 0.0, "MIN_AREA": 0.05,
        "CANDIDATES": out("cand_shp.shp"), "PARALLEL": out("par_shp.shp")})
    with open(out("cand_shp.shp"), "rb") as fh:
        magic = fh.read(4)
    check("a .shp output is a real shapefile (not a GeoPackage named .shp)",
          magic == b"\x00\x00\x27\x0a" and os.path.exists(out("cand_shp.dbf")))
    r = processing.run("qeht:buildheasexchange", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "RAW_DEM": dem, "POINTS": out("cand_shp.shp"),
        "ROAD": out("road.gpkg"), "FILLED": out("fill.tif"), "SNAP": 5, "SNAP_THRESHOLD": 200,
        "OUTPUT": out("cand_shp_pkg.gpkg")})
    md = {row["key"]: row["value"] for row in gpkg.read_table(r["OUTPUT"], "qeht_run_metadata")}
    check("candidates read back from a shapefile (10-character field names); layer sources "
          "and raster tags recorded",
          "crossing_candidates" in dict(gpkg.list_tables(r["OUTPUT"]))
          and _json2.loads(md["parameters_json"])["FDR"].endswith("fdr.tif")
          and md["flat_method"] in ("barnes", "toward") and "fill" in md["conditioning"],
          f"flat {md['flat_method']}, conditioning {md['conditioning']}")
    r = processing.run("qeht:catchmentcharacteristics", {
        "FDR": out("fdr.tif"), "DEM": out("fill.tif"), "RAW_DEM": dem, "FAC": out("fac.tif"),
        "POINTS": out("cand_shp.shp"), "CATCH_OUT": out("cc.shp"), "PATH_OUT": out("cc_lfp.shp")})
    check("characteristics takes a candidate layer (recommended outlets, no snapping) and "
          "writes shapefiles",
          QgsVectorLayer(out("cc.shp"), "cc", "ogr").featureCount() ==
          sum(1 for f in QgsVectorLayer(out("cand_shp.shp"), "c", "ogr").getFeatures()
              if f["recommende"] == 1))
    from osgeo import gdal as _g2
    _ds = _g2.Open(out("fill.tif"))
    _wc = _g2.GetDriverByName("GTiff").Create(out("wc_classes.tif"), _ds.RasterXSize, _ds.RasterYSize, 1, _g2.GDT_Byte)
    _wc.SetGeoTransform(_ds.GetGeoTransform()); _wc.SetProjection(_ds.GetProjection())
    _wc.GetRasterBand(1).WriteArray(np.where(np.arange(_ds.RasterXSize)[None, :] % 2 == 0, 30, 40)
                                    * np.ones((_ds.RasterYSize, 1), np.uint8))
    _wc = None; _ds = None
    try:
        processing.run("qeht:erosionindices", {"RAW_DEM": dem, "FAC": out("fac.tif"),
                                               "R_VALUE": 3000.0, "K_VALUE": 0.03,
                                               "C_RASTER": out("wc_classes.tif"), "OUTPUT": out("ero_bad")})
        refused = False
    except Exception as e:
        refused = "class codes" in str(e)
    check("a land-cover class raster given as the C raster is refused with a hint", refused)



    # ---- v0.18 (A4): floodplain width indicator -----------------------------
    r = processing.run("qeht:buildheasexchange", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "RAW_DEM": dem, "POINTS": out("cand.gpkg"),
        "ROAD": out("road.gpkg"), "START": 1000.0, "FILLED": out("fill.tif"), "COVERAGE": False,
        "FP_MIN_AREA": 1.0, "OUTPUT": out("fp_pkg.gpkg")})
    crf = gpkg.read_table(r["OUTPUT"], "crossings", with_geometry=False)
    mdf = {row["key"]: row["value"] for row in gpkg.read_table(r["OUTPUT"], "qeht_run_metadata")}
    big = [c for c in crf if c["acc_at_outlet_km2"] >= 1.0]
    check("floodplain width at crossings >= 1 km2: profile + HAND widths grow with dz, "
          "smaller crossings empty, settings in metadata",
          big and all(c["fp_method"] == "profile+HAND"
                      and c["fp_w_0p5_m"] <= c["fp_w_1p0_m"] <= c["fp_w_2p0_m"] for c in big)
          and all(c["fp_method"] is None for c in crf if c["acc_at_outlet_km2"] < 1.0)
          and '"min_area_km2": 1.0' in mdf["fp_params_json"],
          ", ".join(f"{c['outlet_uid']} {c['acc_at_outlet_km2']:.1f} km2: {c['fp_w_1p0_m']:.0f} m "
                    f"(HAND {c['fp_hand_w_1p0_m'] if c['fp_hand_w_1p0_m'] is None else round(c['fp_hand_w_1p0_m'])})"
                    for c in big))
    # ---- v0.19 (F5): channel cross-sections at crossings ----------------------
    xst = gpkg.read_table(r["OUTPUT"], "xs_transects")
    grow = all((c["xs_w_0p5_m"] is None or c["xs_w_1p0_m"] is None or c["xs_w_0p5_m"] <= c["xs_w_1p0_m"])
               and (c["xs_w_1p0_m"] is None or c["xs_w_2p0_m"] is None or c["xs_w_1p0_m"] <= c["xs_w_2p0_m"])
               for c in crf)
    errors, _ = validate_exchange(r["OUTPUT"])
    check("channel sections: every crossing has xs_quality, widths grow with dz, one transect "
          "per crossing, settings in metadata, package validates",
          not errors and crf and all(c["xs_quality"] in ("high", "medium", "low") for c in crf)
          and grow and len(xst) == len(crf) and {t["outlet_uid"] for t in xst} == {c["outlet_uid"] for c in crf}
          and '"dist_m": 30.0' in mdf["xs_params_json"],
          ", ".join(f"{c['outlet_uid']} {c['xs_quality']} w1 {c['xs_w_1p0_m']}" for c in crf))
    r = processing.run("qeht:buildheasexchange", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "RAW_DEM": dem, "POINTS": out("cand.gpkg"),
        "ROAD": out("road.gpkg"), "START": 1000.0, "COVERAGE": False, "XS": False,
        "OUTPUT": out("xs_off_pkg.gpkg")})
    check("channel sections can be switched off (no layer, fields empty)",
          "xs_transects" not in dict(gpkg.list_tables(r["OUTPUT"]))
          and all(c["xs_quality"] is None for c in gpkg.read_table(r["OUTPUT"], "crossings",
                                                                   with_geometry=False)))
    # ---- v0.21 (F9): Prepare DEM for hydrology --------------------------------
    import json as _jp
    r = processing.run("qeht:preparedem", {"DEMS": [geo_dem], "SOURCE": "synthetic test DEM",
                                           "OUTPUT": out("prepared.tif"), "LOG": out("prepared.json")})
    from ..core.raster import raster_tags as _rt, read_dem as _rdp
    _, _, pinf = _rdp(r["OUTPUT"])
    with open(r["LOG"], encoding="utf-8") as fh:
        plog = _jp.load(fh)
    check("prepare DEM: geographic synthetic DEM -> auto UTM, native cell in metres, tags and "
          "log with both audits",
          "UTM zone 26N" in pinf.projection_wkt and abs(pinf.cell_width - plog["cell_m"]) < 1e-9
          and _rt(r["OUTPUT"]).get("QEHT_DEM_SOURCE") == "synthetic test DEM"
          and "audit_input" in plog and not plog["audit_output"]["warnings"],
          f"{plog['target_crs']}, cell {plog['cell_m']} m, {r['SUMMARY']}")
    from osgeo import gdal as _gp
    _gp.UseExceptions()
    _gp.Translate(out("dem_tagged.tif"), dem)
    _ds = _gp.Open(out("dem_tagged.tif"), _gp.GA_Update)
    _ds.SetMetadata({"QEHT_DEM_SOURCE": "FABDEM v1.2 (test tag)", "QEHT_DEM_PREP": "QEHT prepare DEM test"})
    _ds = None
    r = processing.run("qeht:buildheasexchange", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "RAW_DEM": out("dem_tagged.tif"),
        "POINTS": out("cand.gpkg"), "XS": False, "OUTPUT": out("tag_pkg.gpkg")})
    mdt = {row["key"]: row["value"] for row in gpkg.read_table(r["OUTPUT"], "qeht_run_metadata")}
    check("package reads dem_source and dem_prep from the DEM tags when left blank",
          mdt["dem_source"].startswith("FABDEM v1.2 (test tag)") and mdt["dem_prep"] ==
          "QEHT prepare DEM test", mdt["dem_source"])
    # ---- v0.27 (F15): DEM uncertainty at crossings ------------------------------------
    try:
        processing.run("qeht:demuncertainty", {"DEM": dem, "CROSSINGS": out("cand.gpkg"),
                                               "OUTPUT": out("unc_x.gpkg")})
        refused = False
    except Exception as e:
        refused = "no default" in str(e).lower() or "There is no default" in str(e)
    check("DEM uncertainty refuses to run without a chosen vertical error", refused)
    r = processing.run("qeht:demuncertainty", {
        "DEM": dem, "CROSSINGS": out("cand.gpkg"), "ID_FIELD": "status", "UNC_PRESET": 2,
        "UNC_N": 3, "UNC_SENS": True, "OUTPUT": out("unc.gpkg"), "LOG": out("unc.json")})
    ul = list(QgsVectorLayer(r["OUTPUT"], "u", "ogr").getFeatures())
    import json as _ju
    with open(r["LOG"], encoding="utf-8") as fh:
        uj = _ju.loads(fh.read())
    check("DEM uncertainty tool (FABDEM preset, N 3): P10 <= P50 <= P90 per crossing, settings, "
          "source and correlation-length sensitivity recorded",
          ul and all(f["unc_area_p10"] <= f["unc_area_p50"] <= f["unc_area_p90"] for f in ul
                     if f["unc_area_p50"] is not None and str(f["unc_area_p50"]) != "NULL")
          and uj["sigma_m"] == 2.5 and "Hawker" in uj["sigma_source"] and uj["n"] == 3
          and len(uj["corr_len_sensitivity"]["outlets"]) == min(3, len(ul)),
          r["SUMMARY"])
    r = processing.run("qeht:buildheasexchange", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "RAW_DEM": dem, "POINTS": out("cand.gpkg"),
        "XS": False, "UNC": True, "UNC_PRESET": 3, "UNC_SIGMA": 1.0, "UNC_N": 2, "UNC_SENS": False,
        "OUTPUT": out("unc_pkg.gpkg")})
    crq = gpkg.read_table(r["OUTPUT"], "crossings", with_geometry=False)
    mdq = {row["key"]: row["value"] for row in gpkg.read_table(r["OUTPUT"], "qeht_run_metadata")}
    errors, _ = validate_exchange(r["OUTPUT"])
    check("DEM uncertainty in the package (custom 1 m, N 2): unc_* on crossings, uncertainty_json, "
          "validates",
          not errors and all(c["unc_n"] == 2 for c in crq) and '"preset": "custom"' in mdq["uncertainty_json"])
    # ---- v0.26 (F16): Graphical Modeler example ------------------------------------
    from qgis.core import QgsProcessingModelAlgorithm
    model = QgsProcessingModelAlgorithm()
    loaded = model.fromFile(os.path.join(here, "examples", "qeht_corridor.model3"))
    r = processing.run(model, {"dem": dem, "road": out("road.gpkg"), "threshold": 200,
                               "design_hydrology_package": out("model_pkg.gpkg")})
    pkm = r.get("design_hydrology_package") or out("model_pkg.gpkg")
    errors, _ = validate_exchange(pkm)
    check("Graphical Modeler example: loads, runs on the example DEM and road, package validates",
          loaded and not errors and len(gpkg.read_table(pkm, "crossings", with_geometry=False)) > 0,
          "; ".join(errors))
    # ---- v0.17 (F7/F8): one-click pipeline and run report --------------------
    import json as _json3
    from ..core.raster import read_dem as _rd3
    _, _, _inf = _rd3(dem)
    km2_200 = 200 * _inf.cell_width * _inf.cell_height / 1e6
    pf = out("pipe_review")
    r = processing.run("qeht:hydrologypipeline", {
        "DEM": dem, "ROAD": out("road.gpkg"), "START": 1000.0, "STAGE": 1, "STREAM_KM2": km2_200,
        "RUN_NAME": "review", "LOAD": False, "OUTPUT": pf})
    check("pipeline, stop for review: candidates, streams and settings written, no report yet",
          os.path.exists(os.path.join(pf, "layers", "crossing_candidates.gpkg"))
          and os.path.exists(os.path.join(pf, "layers", "streams.gpkg"))
          and os.path.exists(os.path.join(pf, "settings.json"))
          and not os.path.exists(os.path.join(pf, "report", "run_report.html"))
          and "Stopped for review" in r["SUMMARY"])
    pf2 = out("pipe_full")
    r = processing.run("qeht:hydrologypipeline", {
        "SETTINGS": os.path.join(pf, "settings.json"),
        "CROSSINGS": os.path.join(pf, "layers", "crossing_candidates.gpkg"), "STAGE": 0,
        "OUTPUT": pf2})
    # the settings file sets everything else; add inputs by editing a copy of it
    with open(os.path.join(pf, "settings.json"), encoding="utf-8") as fh:
        st = _json3.load(fh)
    st["parameters"].update({"SOIL_HSG_R": out("hysogs.tif"), "LANDCOVER": out("worldcover.tif"),
                             "RAIN_ZONES": out("rain_zones.gpkg"), "RAIN_ZONE_FIELD": "zone",
                             "RAIN_MAP": out("map.tif"), "RAIN_R_RELATION": 1, "PACKAGE": True,
                             "BURN": True, "RUN_NAME": "full", "NODATA_OVERRIDE": 0,
                             "TC_P2": 60.0})
    with open(out("settings_edit.json"), "w", encoding="utf-8") as fh:
        _json3.dump(st, fh)
    pf3 = out("pipe_rich")
    r3 = processing.run("qeht:hydrologypipeline", {
        "SETTINGS": out("settings_edit.json"),
        "CROSSINGS": os.path.join(pf, "layers", "crossing_candidates.gpkg"), "STAGE": 0,
        "OUTPUT": pf3})
    need = ["layers/crossings.gpkg", "layers/catchments.gpkg", "layers/flowpaths.gpkg",
            "layers/streams.gpkg", "layers/burn_log.gpkg", "rasters/filled.tif",
            "rasters/flow_direction.tif", "rasters/flow_accumulation.tif", "rasters/dem_burned.tif",
            "rasters/erosion/ls_factor.tif", "tables/catchment_characteristics.csv",
            "tables/crossings.csv", "report/run_report.html", "settings.json",
            "package/design_hydrology.gpkg"]
    missing = [n for n in need if not os.path.exists(os.path.join(pf3, n))]
    check("full run from edited settings: fixed folder layout, burn, erosion, package kept",
          not missing and not os.path.exists(os.path.join(pf3, "_work")), str(missing))
    pkg3 = os.path.join(pf3, "package", "design_hydrology.gpkg")
    errors, _ = validate_exchange(pkg3)
    ca3 = gpkg.read_table(pkg3, "catchments", with_geometry=False)
    md3 = {row["key"]: row["value"] for row in gpkg.read_table(pkg3, "qeht_run_metadata")}
    with open(os.path.join(pf3, "tables", "catchment_characteristics.csv"), encoding="utf-8") as fh:
        lines = fh.read().strip().splitlines()
    head = lines[0].split(",")
    check("package validates; CN, rainfall and R estimate on catchments; one table row per "
          "crossing with the flat-method check",
          not errors and all(c["cn_export"] and c["rain_zone"] and c["rusle_r"] for c in ca3)
          and len(lines) - 1 == len(ca3) and "flat_sensitive" in head and "cn_export" in head
          and "breach" in (md3.get("conditioning_burn") or "") + (md3.get("conditioning") or "")
          and all(c["elev_min_m"] > 1000 for c in ca3)
          and os.path.exists(os.path.join(pf3, "rasters", "dem_input.tif")),
          f"{len(ca3)} catchments; conditioning: {md3.get('conditioning')}")
    with open(os.path.join(pf3, "report", "run_report.html"), encoding="utf-8") as fh:
        rep = fh.read()
    cr3 = gpkg.read_table(pkg3, "crossings", with_geometry=False)
    with open(os.path.join(pf3, "tables", "time_of_concentration.csv"), encoding="utf-8") as fh:
        tcl = fh.read().strip().splitlines()
    md_tc = _json3.loads(md3["tc_params_json"])
    n_fp = sum(1 for c in cr3 if c["tc_note"] != "no flow path")
    n5 = sum(1 for c in cr3 if all(c[f"tc_{m}_min"] is not None for m in
                                   ("kirpich", "kerby_kirpich", "scs_lag", "bransby_williams")))
    n_tr = sum(1 for c in cr3 if c["tc_tr55_min"] is not None)
    import zipfile as _zf
    kz = os.path.join(pf3, "exports", "design_hydrology.kmz")
    xz = os.path.join(pf3, "exports", "design_hydrology.xlsx")
    kml = _zf.ZipFile(kz).read("doc.kml").decode("utf-8") if os.path.exists(kz) else ""
    from ..core.report.xlsx import read_xlsx_cells as _rx
    xc = _rx(xz, 1) if os.path.exists(xz) else {}
    check("pipeline exports: KMZ with every crossing, relief overlay and legend; XLSX with one "
          "row per crossing",
          kml.count("<styleUrl>#x_existing</styleUrl>") + kml.count("<styleUrl>#x_proposed</styleUrl>") == len(cr3) and "gx:LatLonQuad" in kml
          and "files/legend.png" in _zf.ZipFile(kz).namelist()
          and sum(1 for k in xc if k.startswith("A") and k[1:].isdigit() and int(k[1:]) >= 3) == len(cr3),
          f"{kml.count('<Placemark>')} placemarks, {len(xc)} cells")
    # ---- v0.25 (F13): corridor auto-clip ------------------------------------------
    r = processing.run("qeht:autoclip", {"DEM": dem, "ROAD": out("road.gpkg"), "MARGIN": 500.0,
                                         "OUTPUT": out("dem_clip.tif")})
    from ..core.raster import read_dem as _rdc, raster_tags as _rtc
    _, _, cinf = _rdc(r["OUTPUT"])
    _, _, finf = _rdc(dem)
    ctag = _json3.loads(_rtc(r["OUTPUT"]).get("QEHT_AUTOCLIP", "{}"))
    check("auto-clip: smaller grid on the same cell size, window and memory recorded in the tag",
          cinf.rows * cinf.cols < finf.rows * finf.cols and abs(cinf.cell_width - finf.cell_width) < 1e-9
          and ctag.get("clip_peak_gb", 9) < ctag.get("full_peak_gb", 0) and len(ctag["window"]) == 4,
          r["SUMMARY"])
    pfc = out("pipe_clip")
    r = processing.run("qeht:hydrologypipeline", {
        "DEM": dem, "ROAD": out("road.gpkg"), "START": 1000.0, "STREAM_KM2": km2_200, "AUTO_CLIP": True,
        "EROSION": False, "QUICKLOOKS": False, "EXPORTS": False, "PACKAGE": True, "LOAD": False,
        "OUTPUT": pfc})
    pkc = os.path.join(pfc, "package", "design_hydrology.gpkg")
    cac = gpkg.read_table(pkc, "catchments", with_geometry=False)
    mdc = {row["key"]: row["value"] for row in gpkg.read_table(pkc, "qeht_run_metadata")}
    errors, _ = validate_exchange(pkc)
    check("pipeline with the auto-clip: runs on rasters/dem_clip.tif, every catchment has "
          "clip_edge, autoclip_json recorded, validates",
          not errors and os.path.exists(os.path.join(pfc, "rasters", "dem_clip.tif"))
          and cac and all(c["clip_edge"] in (0, 1) for c in cac) and mdc.get("autoclip_json"),
          f"{len(cac)} catchments, {sum(c['clip_edge'] for c in cac)} at the clip edge")
    # ---- v0.24 (F12): check against mapped drainage (the DEM's own streams as the map) --
    r = processing.run("qeht:buildheasexchange", {
        "FDR": os.path.join(pf, "rasters", "flow_direction.tif"),
        "FAC": os.path.join(pf, "rasters", "flow_accumulation.tif"), "RAW_DEM": dem,
        "POINTS": os.path.join(pf, "layers", "crossing_candidates.gpkg"), "ROAD": out("road.gpkg"),
        "START": 1000.0, "COVERAGE": False, "SNAP_THRESHOLD": 200,
        "MAPPED": os.path.join(pf, "layers", "streams.gpkg"), "MAPPED_KM2": km2_200,
        "MAPPED_SOURCE": "DEM streams (self-check)", "OUTPUT": out("mapped_pkg.gpkg")})
    errors, _ = validate_exchange(r["OUTPUT"])
    crm = gpkg.read_table(r["OUTPUT"], "crossings", with_geometry=False)
    cam = gpkg.read_table(r["OUTPUT"], "catchments", with_geometry=False)
    mdm = _json3.loads({row["key"]: row["value"] for row in
                        gpkg.read_table(r["OUTPUT"], "qeht_run_metadata")}["mapped_drainage_json"])
    tabsm = dict(gpkg.list_tables(r["OUTPUT"]))
    agr = [c["map_agrees"] for c in crm if c["map_agrees"] is not None]
    check("mapped drainage against the DEM's own streams: precision and recall near 1, every "
          "crossing on the map, per-catchment scores, layers and metadata, validates",
          not errors and mdm["precision"] > 0.95 and mdm["recall"] > 0.95 and agr and all(agr)
          and all(c["map_precision"] is not None for c in cam)
          and "drainage_divergence" in tabsm and "mapped_rivers_used" in tabsm,
          f"P {mdm['precision']:.3f} R {mdm['recall']:.3f} F1 {mdm['f1']:.3f}; {len(agr)} of "
          f"{len(crm)} crossings compared")
    r = processing.run("qeht:exportpackage", {"PACKAGE": pkg3, "TITLE": "Site C <test>",
                                              "KMZ": out("exp.kmz"), "XLSX": out("exp.xlsx")})
    check("export tool on a package: KMZ and XLSX written",
          os.path.exists(r["KMZ"]) and os.path.exists(r["XLSX"])
          and "Site C &lt;test&gt;" in _zf.ZipFile(r["KMZ"]).read("doc.kml").decode("utf-8"))
    check("time of concentration: Kirpich, Kerby, SCS lag, Bransby-Williams on every crossing "
          "with a flow path (land cover + CN), TR-55 where a bank-full section exists, flags set, table and report "
          "section, P2 recorded",
          n5 == n_fp and n_fp >= len(cr3) - 2 and n_tr >= 1 and all(c["tc_kirpich_flag"] for c in cr3 if c["tc_note"] != "no flow path")
          and len(tcl) - 1 == len(cr3) and "Time of concentration" in rep
          and md_tc["p2_mm"] == 60.0 and md_tc["sheet_n_lookup"].startswith("PROXY")
          and all(c["tc_tr55_min"] is not None or "TR-55" in (c["tc_note"] or "") for c in cr3
                  if c["tc_note"] != "no flow path"),
          f"{n5}/{n_fp} with four methods, {n_tr} with TR-55; e.g. "
          + ", ".join(f"{c['outlet_uid']} K {c['tc_kirpich_min'] or 0:.0f} SCS "
                      f"{c['tc_scs_lag_min'] or 0:.0f} min" for c in cr3[:3]))
    check("run report: title, every crossing, plan, PROXY / ESTIMATE labels, settings path",
          "<h1>full</h1>" in rep and all(c["outlet_uid"] in rep for c in ca3) and "<svg" in rep
          and "PROXY" in rep and "R ESTIMATE" in rep and "settings.json" in rep)
    ql = gpkg.read_table(pkg3, "rasters", with_geometry=False)
    check("quicklooks: relief, accumulation and erosion classes as PNG + world file + legend, "
          "indexed in the package; terrain under the report plan; floodplain widths in the table",
          {"relief", "flow_accumulation", "spi_class"} <= {q["name"] for q in ql}
          and all(os.path.exists(os.path.normpath(os.path.join(os.path.dirname(pkg3), q[k])))
                  for q in ql for k in ("png", "world_file", "legend_json"))
          and "data:image/png;base64," in rep and "fp_w_1p0_m" in head
          and not os.path.exists(os.path.join(pf3, "_work_bg")),
          ", ".join(q["name"] for q in ql))
    pf4 = out("pipe_points")
    r4 = processing.run("qeht:hydrologypipeline", {
        "DEM": dem, "CROSSINGS": ptsfile, "FLAT_METHOD": 0, "STREAM_KM2": km2_200,
        "EROSION": False, "LOAD": False, "OUTPUT": pf4})
    direct = processing.run("qeht:buildheasexchange", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "RAW_DEM": dem, "POINTS": ptsfile,
        "SNAP": 5, "SNAP_THRESHOLD": 200, "OUTPUT": out("direct_for_pipe.gpkg")})
    a_pipe = sorted(round(f["area_km2"], 6) for f in
                    QgsVectorLayer(os.path.join(pf4, "layers", "catchments.gpkg"), "c", "ogr").getFeatures())
    a_dir = sorted(round(c["area_km2"], 6) for c in
                   gpkg.read_table(direct["OUTPUT"], "catchments", with_geometry=False))
    check("pour points, package off: no package kept, areas identical to the tools run one by one",
          r4["PACKAGE_FILE"] == "" and not os.path.exists(os.path.join(pf4, "package"))
          and not os.path.exists(os.path.join(pf4, "_work")) and a_pipe == a_dir and a_pipe,
          f"{a_pipe} vs {a_dir}")
    r5 = processing.run("qeht:runreport", {"PACKAGE": direct["OUTPUT"], "TITLE": "direct",
                                           "OUTPUT": out("direct_report.html")})
    check("run report tool on any package (+ characteristics CSV)",
          os.path.exists(r5["OUTPUT"]) and os.path.exists(r5["CSV"]))

    print("\n" + ("ALL QGIS SMOKE CHECKS PASSED" if not FAILURES
                  else f"{len(FAILURES)} FAILURE(S): " + "; ".join(FAILURES)))
    print("outputs in", tmp)
    if app is not None:
        app.exitQgis()
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
