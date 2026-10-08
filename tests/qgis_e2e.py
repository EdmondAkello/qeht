# -*- coding: utf-8 -*-
"""End-to-end check of 0.19-0.27 on a fully synthetic site, through QGIS Processing.

A 2 m grid with a trapezoidal channel running south (bed, banks and floodplain
from the generator parameters below), a road across it, a mapped waterway on
the channel, land cover (grass west, crops east), soil groups and a land-cover
scenario. Every expected value is derived from the generator parameters or
recomputed with the core functions from the package's own fields - nothing is
typed in.

    QT_QPA_PLATFORM=offscreen python3.12 -m qeht.tests.qgis_e2e
"""

import json
import math
import os
import sys
import tempfile
import zipfile

import numpy as np

FAILURES = []

# generator parameters (the expected values are derived from these)
CS = 2.0                       # cell size (m)
ROWS, COLS = 300, 201
C0 = 100                       # channel centre column
BED_W, BANK_HV, BANK_H = 12.0, 2.0, 4.0
FP_SLOPE = 0.05                # floodplain cross slope (1:20)
LONG_S = 0.005                 # along-channel slope
THALWEG = 0.01                 # bed cross-fall to the centre line (a thalweg, so D8 keeps one stream)
ROAD_ROW = 150
X0, Y0 = 400000.0, 9800000.0


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def profile(u):
    """Height above the thalweg at offset u (m)."""
    u = np.abs(u)
    a = u - BED_W / 2
    edge = THALWEG * BED_W / 2
    return np.where(a <= 0, THALWEG * u, edge + np.where(a <= BANK_HV * BANK_H, a / BANK_HV,
                                                         BANK_H + (a - BANK_HV * BANK_H) * FP_SLOPE))


def width_at(level):
    """Width below `level` above the thalweg, solved on the generator profile (bisection)."""
    lo, hi = 0.0, 1000.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if profile(mid) <= level else (lo, mid)
    return 2 * lo


def build(tmp):
    from osgeo import gdal, osr, ogr
    gdal.UseExceptions(); ogr.UseExceptions()
    srs = osr.SpatialReference(); srs.ImportFromEPSG(32737)
    gt = (X0, CS, 0.0, Y0, 0.0, -CS)
    r, c = np.mgrid[0:ROWS, 0:COLS].astype(float)
    dem = 1000.0 + LONG_S * (ROWS - 1 - r) * CS + profile((C0 - c) * CS)

    def tif(name, a, dtype=gdal.GDT_Float32, nd=None):
        p = os.path.join(tmp, name)
        ds = gdal.GetDriverByName("GTiff").Create(p, COLS, ROWS, 1, dtype)
        ds.SetGeoTransform(gt); ds.SetProjection(srs.ExportToWkt())
        b = ds.GetRasterBand(1); b.WriteArray(a)
        if nd is not None:
            b.SetNoDataValue(nd)
        ds = None
        return p
    paths = {"dem": tif("dem.tif", dem.astype(np.float32)),
             "lc": tif("lc.tif", np.where(c < C0, 30, 40).astype(np.uint8), gdal.GDT_Byte, 0),
             "lc_scn": tif("lc_scn.tif", np.where(c < C0, np.where(r < 60, 50, 30), 40).astype(np.uint8),
                           gdal.GDT_Byte, 0),
             "hsg": tif("hsg.tif", np.full((ROWS, COLS), 2, np.uint8), gdal.GDT_Byte, 255)}

    def lines(name, coords_list):
        p = os.path.join(tmp, name)
        ds = ogr.GetDriverByName("GPKG").CreateDataSource(p)
        lyr = ds.CreateLayer(name[:-5], srs, ogr.wkbLineString)
        lyr.CreateField(ogr.FieldDefn("name", ogr.OFTString))
        for co in coords_list:
            f = ogr.Feature(lyr.GetLayerDefn()); f.SetField("name", "Synthetic river")
            g = ogr.Geometry(ogr.wkbLineString)
            for x, y in co:
                g.AddPoint_2D(x, y)
            f.SetGeometry(g); lyr.CreateFeature(f)
        ds = None
        return p
    xc = X0 + (C0 + 0.5) * CS
    yr = Y0 - (ROAD_ROW + 0.5) * CS
    paths["road"] = lines("road.gpkg", [[(X0 + 21.0, yr), (X0 + (COLS - 10) * CS, yr)]])
    paths["mapped"] = lines("mapped.gpkg", [[(xc, Y0 - 0.5 * CS), (xc, Y0 - (ROWS - 0.5) * CS)]])
    return paths, dem, gt


def main():
    from qgis.core import QgsApplication
    app = QgsApplication([], False)
    app.initQgis()
    sys.path.append("/usr/share/qgis/python/plugins")
    import processing
    from processing.core.Processing import Processing
    Processing.initialize()
    from ..processing_provider.provider import QehtProvider
    provider = QehtProvider()
    QgsApplication.processingRegistry().addProvider(provider)
    from ..core.interop import gpkg
    from ..core.interop.heas_exchange import validate_exchange
    from ..core.runoff.tc import kirpich_min, bransby_williams_min, scs_lag_min, tc_block
    from ..core.runoff.curve_number import tr55_lookup, amc_convert
    from ..core.report.xlsx import read_xlsx_cells, col_letter
    from ..core.report.xlsx_layout import CHARACTERISTICS

    tmp = tempfile.mkdtemp(prefix="qeht_e2e_")
    out = lambda n: os.path.join(tmp, n)   # noqa: E731
    paths, dem, gt = build(tmp)
    print("Synthetic site in", tmp)
    thr = 400.0
    processing.run("qeht:filldepressions", {"DEM": paths["dem"], "OUTPUT": out("fill.tif")})
    processing.run("qeht:flowdirection", {"DEM": out("fill.tif"), "OUTPUT": out("fdr.tif")})
    processing.run("qeht:flowaccumulation", {"FDR": out("fdr.tif"), "OUTPUT": out("fac.tif")})
    processing.run("qeht:streamnetwork", {"FDR": out("fdr.tif"), "FAC": out("fac.tif"), "MODE": 0,
                                          "THRESHOLD": thr, "STREAMS": out("str.tif"),
                                          "ORDER": out("ord.tif")})
    processing.run("qeht:crossingcandidates", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "STREAMS": out("str.tif"), "THRESHOLD": thr,
        "ROAD": paths["road"], "CANDIDATES": out("cand.gpkg"), "PARALLEL": out("par.gpkg")})
    processing.run("qeht:erosionindices", {
        "RAW_DEM": paths["dem"], "FAC": out("fac.tif"), "CHANNEL_CELLS": thr, "R_VALUE": 3000.0,
        "K_VALUE": 0.03, "WORLDCOVER": paths["lc"], "OUTPUT": out("erosion")})
    p2 = 55.0
    r = processing.run("qeht:buildheasexchange", {
        "FDR": out("fdr.tif"), "FAC": out("fac.tif"), "RAW_DEM": paths["dem"], "ORDER": out("ord.tif"),
        "POINTS": out("cand.gpkg"), "ROAD": paths["road"], "SNAP_THRESHOLD": thr, "COVERAGE": False,
        "FILLED": out("fill.tif"), "SOIL_HSG_R": paths["hsg"], "LANDCOVER": paths["lc"],
        "LANDCOVER_SCN": paths["lc_scn"], "SCENARIO_NAME": "north-west built", "EROSION": out("erosion"),
        "MAPPED": paths["mapped"], "MAPPED_NAME_FIELD": "name", "MAPPED_KM2": thr * CS * CS / 1e6,
        "MAPPED_TOL": 2 * CS, "TC_P2": p2, "UNC": True, "UNC_PRESET": 3, "UNC_SIGMA": 0.0,
        "UNC_N": 2, "UNC_SENS": False, "OUTPUT": out("pkg.gpkg")})
    pkg = r["OUTPUT"]
    errors, _ = validate_exchange(pkg)
    check("package validates", not errors, "; ".join(errors))
    cr = gpkg.read_table(pkg, "crossings", with_geometry=False)
    ca = {a["outlet_uid"]: a for a in gpkg.read_table(pkg, "catchments", with_geometry=False)}
    fp = {a["outlet_uid"]: a for a in gpkg.read_table(pkg, "flowpaths", with_geometry=False)}
    md = {row["key"]: row["value"] for row in gpkg.read_table(pkg, "qeht_run_metadata")}
    main_x = [x for x in cr if abs(x["outlet_x"] - (X0 + (C0 + 0.5) * CS)) < 2 * CS]
    check("one crossing on the channel", len(main_x) == 1, f"{len(cr)} crossing(s)")
    if not main_x:
        return 1
    x = main_x[0]
    uid = x["outlet_uid"]
    a, f = ca[uid], fp[uid]

    print("\nF5 channel section (expected from the generator)")
    d_bf = float(profile(BED_W / 2 + BANK_HV * BANK_H))           # bank top above the thalweg
    w_bf = width_at(d_bf)
    step = CS / 2
    check(f"bank-full width {w_bf:.2f} m within one station step", abs(x["xs_bankfull_w_m"] - w_bf) <= step,
          f"{x['xs_bankfull_w_m']:.2f}")
    check(f"bank-full depth {d_bf:.2f} m", abs(x["xs_bankfull_d_m"] - d_bf) <= 0.05,
          f"{x['xs_bankfull_d_m']:.3f}")
    check(f"side slopes {BANK_HV:g} H:V", abs(x["xs_side_slope_l"] - BANK_HV) <= 0.15
          and abs(x["xs_side_slope_r"] - BANK_HV) <= 0.15,
          f"{x['xs_side_slope_l']:.2f} / {x['xs_side_slope_r']:.2f}")
    exp_w1 = width_at(1.0)
    check(f"width at bed + 1 m = {exp_w1:.2f} m", abs(x["xs_w_1p0_m"] - exp_w1) <= step, f"{x['xs_w_1p0_m']:.2f}")
    row_xs = ROAD_ROW + round(x["xs_dist_m"] / CS)
    check("bed level = generator bed at the section row",
          abs(x["xs_bed_m"] - (1000.0 + LONG_S * (ROWS - 1 - row_xs) * CS)) < 1e-3,
          f"{x['xs_bed_m']:.3f}")
    check("quality high", x["xs_quality"] == "high", str(x["xs_note"]))

    print("\nF10 time of concentration (recomputed from the package's own fields)")
    check("Kirpich = kirpich_min(lfp_length_m, lfp_slope_1085)",
          math.isclose(x["tc_kirpich_min"], kirpich_min(f["lfp_length_m"], f["lfp_slope_1085"])))
    check("Bransby-Williams from length, area and 10-85 slope",
          math.isclose(x["tc_bransby_williams_min"],
                       bransby_williams_min(f["lfp_length_m"] / 1000, a["area_km2"], 1000 * f["lfp_slope_1085"])))
    check("SCS lag from length, cn_ii and the Horn slope",
          math.isclose(x["tc_scs_lag_min"], scs_lag_min(f["lfp_length_m"], a["cn_ii"], 100 * a["catch_slope_horn"])))
    basis = json.loads(x["tc_basis_json"])["tr55"]
    again = tc_block(f, a, x, sheet_n=basis["sheet_n"], kerby_n=json.loads(x["tc_basis_json"])
                     ["kerby_kirpich"]["N"], p2_mm=p2)
    check("TR-55 and Kerby reproduce from the package fields and the section", x["tc_tr55_min"] is not None
          and math.isclose(again["tc_tr55_min"], x["tc_tr55_min"])
          and math.isclose(again["tc_kerby_kirpich_min"], x["tc_kerby_kirpich_min"]),
          f"TR-55 {x['tc_tr55_min']:.1f} min")
    check("P2 recorded as given", json.loads(md["tc_params_json"])["p2_mm"] == p2)

    print("\nF12 mapped drainage (the mapped line is the channel)")
    mj = json.loads(md["mapped_drainage_json"])
    check("crossing on the mapped river, named", x["map_agrees"] == 1 and x["map_river_dist_m"] <= CS
          and x["map_river_name"] == "Synthetic river")
    check("catchment precision = 1 (every DEM stream cell is on the map)",
          math.isclose(a["map_precision"], 1.0), f"{a['map_precision']}")
    div = gpkg.read_table(pkg, "drainage_divergence", with_geometry=False)
    check("overall precision = 1, recall above 0.95 (the map runs above the channel head), "
          "no divergence reach",
          math.isclose(mj["precision"], 1.0) and mj["recall"] > 0.95 and not div,
          f"P {mj['precision']:.3f} R {mj['recall']:.3f}, {len(div)} reach(es)")

    print("\nF14 land-cover scenario (shares and lookup)")
    lut = tr55_lookup("fair")
    gB = lambda k: lut[k][1]   # noqa: E731  HSG B
    exp_cn = (a["lc_pct_grass"] * gB(30) + a["lc_pct_crop"] * gB(40)) / 100
    check("baseline cn_ii = shares x TR-55 B values", math.isclose(a["cn_ii"], exp_cn, rel_tol=1e-9),
          f"{a['cn_ii']:.3f} vs {exp_cn:.3f}")
    exp_d = a["lc_pct_built_scn"] * (gB(50) - gB(30)) / 100
    check("d_cn = built share x (CN built - CN grass)", math.isclose(a["d_cn"], exp_d, abs_tol=1e-9),
          f"{a['d_cn']:.3f} vs {exp_d:.3f}")
    check("cn_export_scn at the baseline AMC", math.isclose(a["cn_export_scn"],
                                                         amc_convert(a["cn_ii_scn"], a["cn_amc"])))

    print("\nR3 side-drain siltation (north drains to the road, south away)")
    cs = gpkg.read_table(pkg, "corridor_sti", with_geometry=False)
    check("left (north, upstream) side has values; right (south) side none",
          cs and any(c["sti_p90_lhs"] is not None for c in cs)
          and all(c["sti_p90_rhs"] is None for c in cs))

    print("\nF15 DEM uncertainty with sigma 0")
    check("P10 = P50 = P90 = the package area; nothing lost or switched",
          all(math.isclose(x[k], a["area_km2"]) for k in ("unc_area_p10", "unc_area_p50", "unc_area_p90"))
          and x["unc_lost_pct"] == 0 and (x["unc_switch_pct"] or 0) == 0)
    check("P50 LFP length = the package LFP length", math.isclose(x["unc_lfp_p50"], f["lfp_length_m"]))

    print("\nF11 exports")
    r = processing.run("qeht:exportpackage", {"PACKAGE": pkg, "KMZ": out("e.kmz"), "XLSX": out("e.xlsx")})
    kml = zipfile.ZipFile(r["KMZ"]).read("doc.kml").decode()
    check("KMZ: one crossing placemark per crossing",
          kml.count("<styleUrl>#x_existing</styleUrl>") + kml.count("<styleUrl>#x_proposed</styleUrl>") == len(cr))
    cells = read_xlsx_cells(r["XLSX"], 1)
    heads = [cells.get(f"{col_letter(i)}1") for i in range(len(CHARACTERISTICS))]
    ai = heads.index("Area")
    row = next(rr for rr in range(3, 3 + len(cr)) if cells.get(f"A{rr}") == uid)
    check("XLSX area = package area", math.isclose(cells[f"{col_letter(ai)}{row}"], a["area_km2"]))
    ti = heads.index("Tc Kirpich")
    check("XLSX Tc Kirpich = package value", math.isclose(cells[f"{col_letter(ti)}{row}"], x["tc_kirpich_min"]))

    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): " + "; ".join(FAILURES))
        return 1
    print("ALL END-TO-END CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
