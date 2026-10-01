# -*- coding: utf-8 -*-
"""v0.11 validation: soil parameters (WP-C). Bare Python + NumPy.

    python -m qeht.tests.test_soils
"""

import math
import os
import sqlite3
import sys
import tempfile

import numpy as np

from ..core.soils.usle_k import (williams_k, dg_k, cfrg, usda_texture, hsg_proxy,
                                 K_US_TO_SI)
from ..core.soils.sotwis import load_sotwis, load_attributes
from ..core.soils.catchment import soil_block, SOIL_FIELDS
from ..core.geometry.rasterize import rasterize_polygons

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def test_williams():
    print("\n1. Williams/EPIC K")
    # hand computation, sand 40 silt 40 clay 20 OC 1.5 %
    f_csand = 0.2 + 0.3 * math.exp(-0.0256 * 40 * (1 - 40 / 100))
    f_clsi = (40 / 60) ** 0.3
    f_orgc = 1 - 0.25 * 1.5 / (1.5 + math.exp(3.72 - 2.95 * 1.5))
    sn1 = 0.6
    f_hisand = 1 - 0.7 * sn1 / (sn1 + math.exp(-5.51 + 22.9 * sn1))
    k_us, k_si, f = williams_k(40, 40, 20, 1.5)
    check("hand-computed factors reproduced",
          all(math.isclose(a, b, rel_tol=1e-12) for a, b in (
              (f["f_csand"], f_csand), (f["f_clsi"], f_clsi), (f["f_orgc"], f_orgc),
              (f["f_hisand"], f_hisand))),
          f"f_csand {f_csand:.4f} f_clsi {f_clsi:.4f} f_orgc {f_orgc:.4f} f_hisand {f_hisand:.4f}")
    check("K_US = product, K_SI = 0.1317 K_US",
          math.isclose(k_us, f_csand * f_clsi * f_orgc * f_hisand)
          and math.isclose(k_si, k_us * K_US_TO_SI), f"K_US {k_us:.4f}, K_SI {k_si:.5f}")
    _, _, f = williams_k(60, 20, 20, 1.0)
    check("f_csand uses 0.0256, not the misprinted 0.256",
          math.isclose(f["f_csand"], 0.2 + 0.3 * math.exp(-1.2288)) and f["f_csand"] > 0.28,
          f"{f['f_csand']:.4f} (0.256 would give 0.2000)")
    check("more organic carbon -> lower K",
          williams_k(40, 40, 20, 3.0)[1] < williams_k(40, 40, 20, 0.5)[1])
    check("silt loam more erodible than sand",
          williams_k(20, 65, 15, 1.0)[1] > williams_k(90, 5, 5, 1.0)[1])
    check("missing input -> NaN", math.isnan(williams_k(-1, 40, 20, 1.0)[1]))
    kd = [dg_k(s, si, 100 - s - si) for s, si in ((90, 5), (40, 40), (20, 65), (10, 30))]
    check("Dg-based K within Renard's range (<= 0.0439) and silt highest",
          max(kd) <= 0.04391 and kd[2] == max(kd), ", ".join(f"{k:.4f}" for k in kd))
    check("CFRG = exp(-0.053 rock)", math.isclose(cfrg(10), math.exp(-0.53)) and cfrg(0) == 1.0)


def test_texture_hsg():
    print("\n2. USDA texture and hydrologic-group proxy")
    cases = [((90, 5, 5), "sand"), ((80, 10, 10), "loamy sand"), ((60, 30, 10), "sandy loam"),
             ((40, 40, 20), "loam"), ((20, 65, 15), "silt loam"), ((5, 85, 10), "silt"),
             ((60, 15, 25), "sandy clay loam"), ((30, 35, 35), "clay loam"),
             ((10, 55, 35), "silty clay loam"), ((50, 5, 45), "sandy clay"),
             ((5, 50, 45), "silty clay"), ((20, 20, 60), "clay")]
    wrong = [(c, e, usda_texture(*c)) for c, e in cases if usda_texture(*c) != e]
    check("12 USDA classes", not wrong, str(wrong))
    check("HSG proxy from texture", hsg_proxy("sand") == "A" and hsg_proxy("loam") == "B"
          and hsg_proxy("sandy clay loam") == "C" and hsg_proxy("clay") == "D")
    check("poor drainage -> D, imperfect moves A/B to C",
          hsg_proxy("sand", "P") == "D" and hsg_proxy("loam", "I") == "C"
          and hsg_proxy("clay", "W") == "D")


def _make_db(path):
    con = sqlite3.connect(path)
    con.executescript("""
CREATE TABLE SOTERparameterEstimates (CLAF, PRID, Drain, Layer, TopDep, BotDep, CFRAG,
  SDTO, STPC, CLPC, PSCL, BULK, TAWK, TAWC, CECS, BSAT, CECc, PHAQ, TCEQ, GYPS, ELCO,
  TOTC, TOTN, ECEC, ALSA, ESP);
CREATE TABLE SOTERunitComposition (ISOC, SUID, NEWSUID, SoilMapUnit,
  SOIL1, PROP1, PRID1, SOIL2, PROP2, PRID2, SOIL3, PROP3, PRID3);
CREATE TABLE SOTERflagTTRrules (CLAF, PRID, Layer, Newtopdep, Newbotdep, TTRsub, TTRmain, TTRfinal);
""")
    rows = [("A", "P1", "W", "D1", 0, 20, 10, 40, 40, 20, "M", 1.3, 30, 20),   # TOTC g/kg=20 -> 2 %
            ("A", "P1", "W", "D2", 20, 40, 20, 20, 50, 30, "M", 1.4, 30, 10),
            ("B", "P2", "P", "D1", 0, 10, 0, 70, 20, 10, "C", 1.5, 30, 5),
            ("C", "P3", "W", "D1", 0, 20, 0, -1, 40, -1, "M", 1.3, 30, 10)]      # missing
    for r in rows:
        con.execute("INSERT INTO SOTERparameterEstimates (CLAF, PRID, Drain, Layer, TopDep, "
                    "BotDep, CFRAG, SDTO, STPC, CLPC, PSCL, BULK, TAWC, TOTC) VALUES "
                    "(?,?,?,?,?,?,?,?,?,?,?,?,?,?)", r)
    con.execute("INSERT INTO SOTERunitComposition (NEWSUID, PRID1, PROP1, PRID2, PROP2) "
                "VALUES ('U1', 'P1', 60, 'P2', 40)")
    con.execute("INSERT INTO SOTERunitComposition (NEWSUID, PRID1, PROP1) VALUES ('U2', 'P3', 100)")
    con.execute("INSERT INTO SOTERunitComposition (NEWSUID, PRID1, PROP1, PRID2, PROP2) "
                "VALUES ('U3', 'P1', 50, 'P3', 50)")
    con.execute("INSERT INTO SOTERflagTTRrules VALUES ('A','P1','D1',0,20,'x','A3','1a')")
    con.execute("INSERT INTO SOTERflagTTRrules VALUES ('A','P1','D2',20,40,'x','B2','1a')")
    con.commit(); con.close()


def test_sotwis_loader():
    print("\n3. SOTWIS loader (synthetic SQLite with the SOTWIS schema)")
    db = os.path.join(tempfile.mkdtemp(), "t.db")
    _make_db(db)
    units, info = load_sotwis(db)
    u1 = units["U1"]
    k1 = williams_k(40, 40, 20, 2.0)[1]
    k2 = williams_k(70, 20, 10, 0.5)[1]
    check("0-20 cm: component-weighted K (60/40), K per component",
          math.isclose(u1.k_si, 0.6 * k1 + 0.4 * k2, rel_tol=1e-12), f"{u1.k_si:.5f}")
    check("TOTC g/kg -> OC %", math.isclose(u1.oc_pct, 0.6 * 2.0 + 0.4 * 0.5))
    check("shallow D1 (0-10 cm) still used for 0-20", math.isclose(u1.sand, 0.6 * 40 + 0.4 * 70))
    check("drainage = class with the largest share", u1.drain == "W")
    check("dominant profile and its top TTR rule", u1.dominant_prid == "P1" and u1.ttr_main == "A3")
    u30, _ = load_sotwis(db, 0, 30)
    s = u30["U1"].sand
    exp = 0.6 * (40 * 20 + 20 * 10) / 30 + 0.4 * 70
    check("0-30 cm: depth-weighted across D1/D2 (20/10 cm)", math.isclose(s, exp), f"{s:.3f} vs {exp:.3f}")
    check("-1 = missing: unit with no data has no K and 0 % data",
          (units["U2"].k_si is None or math.isnan(units["U2"].k_si)) and units["U2"].share_with_data == 0)
    check("partly missing unit: renormalised, share_with_data 50 %",
          math.isclose(units["U3"].k_si, k1) and units["U3"].share_with_data == 50.0)
    check("provenance", info["soil_depth_cm"] == "0-20" and "SOTWIS" in info["soil_dataset"])
    g, gi = load_attributes([{"id": 7, "sand": 40, "silt": 40, "clay": 20, "oc": 20}],
                            oc_is_gkg=True)
    check("generic attribute loader", math.isclose(g[7].k_si, k1))


def test_csv_loader():
    print("\n3b. Soil table CSV joined on a unit code (v0.14)")
    import tempfile as _tf
    from ..core.soils.sotwis import load_csv_units, load_attributes
    d = _tf.mkdtemp()
    p = os.path.join(d, "soils.csv")
    with open(p, "w", encoding="utf-8") as f:
        f.write("Unit;Sand;Silt;Clay;OC;bulk\n12;40;30;30;1.5;1.3\nKE7;70;15;15;0.6;1.5\n")
    units, info = load_csv_units(p, "unit")
    ref, _ = load_attributes([{"id": "12", "sand": 40, "silt": 30, "clay": 30, "oc": 1.5,
                               "bulk": 1.3}], bulk="bulk")
    check("CSV (semicolon, any header case) keyed by unit code as text, same K as attributes",
          set(units) == {"12", "KE7"} and abs(units["12"].k_si - ref["12"].k_si) < 1e-12
          and "soils.csv" in info["soil_dataset"])
    with open(p, "w", encoding="utf-8") as f:
        f.write("unit,sand,silt,clay,oc\n1,40,30,30,1\n1,50,25,25,1\n")
    try:
        load_csv_units(p, "unit"); dup = False
    except ValueError as e:
        dup = "more than once" in str(e)
    with open(p, "w", encoding="utf-8") as f:
        f.write("unit,sand,silt,clay\n1,40,30,30\n")
    try:
        load_csv_units(p, "unit"); miss = False
    except ValueError as e:
        miss = "oc" in str(e)
    check("duplicate units and missing columns are refused with a message", dup and miss)


def test_catchment_block():
    print("\n4. Rasterised soil map -> catchment block (exact weights, coverage)")
    gt = (0.0, 10.0, 0.0, 100.0, 0.0, -10.0)          # 10 x 30 cells of 10 m
    shape = (10, 30)
    strips = [([[(0, 0), (100, 0), (100, 100), (0, 100)]], 1),
              ([[(100, 0), (200, 0), (200, 100), (100, 100)]], 2),
              ([[(200, 0), (280, 0), (280, 100), (200, 100)]], 3)]   # last 2 columns: no soil
    grid = rasterize_polygons(strips, gt, shape)
    check("strips rasterised by cell centre", (grid == 1).sum() == 100 and (grid == 2).sum() == 100
          and (grid == 3).sum() == 80 and (grid == 0).sum() == 20)
    units, info = load_attributes([
        {"id": "A", "sand": 80, "silt": 10, "clay": 10, "oc": 1.0, "drain": "W"},
        {"id": "B", "sand": 20, "silt": 60, "clay": 20, "oc": 2.0, "drain": "M"},
        {"id": "C", "sand": 30, "silt": 30, "clay": 40, "oc": 1.5, "drain": "P"}],
        drain="drain")
    idx = {1: "A", 2: "B", 3: "C"}
    mask = np.zeros(shape, bool); mask[:, 5:] = True       # 50 A, 100 B, 80 C, 20 none
    b = soil_block(mask, grid, units, idx, info)
    w = {"A": 50, "B": 100, "C": 80}
    exp_k = sum(w[u] * units[u].k_si for u in w) / 230
    check("K = cell-weighted unit K (50/100/80)", math.isclose(b["usle_k"], exp_k), f"{b['usle_k']:.5f}")
    check("sand exact weights", math.isclose(b["soil_sand_pct"], (50 * 80 + 100 * 20 + 80 * 30) / 230))
    check("coverage = 230/250 = 92 %", math.isclose(b["soil_coverage_pct"], 92.0))
    check("dominant unit B, drainage M, 3 units", b["soil_dominant_unit"] == "B"
          and b["soil_drain_class"] == "M" and b["soil_units"] == 3)
    check("texture and HSG proxy filled", b["soil_texture"] is not None and b["soil_hsg_proxy"] in "ABCD")
    far = rasterize_polygons([([[(5000, 5000), (5100, 5000), (5100, 5100)]], 1)], gt, shape)
    b2 = soil_block(mask, far, units, idx, info)
    check("CRS mismatch / no overlap -> coverage 0 %, values empty",
          b2["soil_coverage_pct"] == 0 and b2["usle_k"] is None and b2["soil_dataset"])
    check("block has exactly the SOIL_FIELDS keys", set(b) == {f for f, _ in SOIL_FIELDS})


def test_pipeline_soil():
    print("\n5. Soil block on exchange catchments")
    from .test_interop import synthetic_run, synthetic_case
    from ..core.interop.heas_exchange import build_exchange_records
    (cr, ca, fp, _, _), g = synthetic_run()
    dem, valid, gt, cs = synthetic_case()
    rows, cols = dem.shape
    half = cols * cs / 2
    polys = [([[(gt[0], gt[3] - rows * cs), (gt[0] + half, gt[3] - rows * cs),
                (gt[0] + half, gt[3]), (gt[0], gt[3])]], 1),
             ([[(gt[0] + half, gt[3] - rows * cs), (gt[0] + cols * cs, gt[3] - rows * cs),
                (gt[0] + cols * cs, gt[3]), (gt[0] + half, gt[3])]], 2)]
    grid = rasterize_polygons(polys, gt, dem.shape)
    units, info = load_attributes([
        {"id": "W", "sand": 30, "silt": 40, "clay": 30, "oc": 1.2, "drain": "W"},
        {"id": "E", "sand": 60, "silt": 20, "clay": 20, "oc": 0.8, "drain": "M"}],
        drain="drain")
    pts = [dict(x=a["input_x"], y=a["input_y"], fid=a["outlet_id"], source_id=None) for _, a in cr]
    from ..core.watershed.delineate import extract_streams
    cr2, ca2, _, _, _ = build_exchange_records(
        g["direction"], g["valid"], g["accum"], g["dem"], g["gt"], pts, snap_radius_cells=3,
        stream_mask=extract_streams(g["accum"], g["valid"], threshold_cells=100),
        soil=(grid, {1: "W", 2: "E"}, units, info))
    ok = all(a["soil_coverage_pct"] == 100 and a["usle_k"] and a["soil_dataset"] for _, a in ca2)
    check("every catchment gets the soil block, 100 % coverage", ok,
          ", ".join(f"{a['outlet_uid']} K={a['usle_k']:.4f}" for _, a in ca2))
    check("hydrology unchanged by the soil option",
          [a["area_km2"] for _, a in ca2] == [a["area_km2"] for _, a in ca])


def main(argv=None):
    print("=" * 62)
    print("QEHT SOILS - v0.11 (WP-C: SOTWIS, USLE K, texture, HSG proxy)")
    print("=" * 62)
    test_williams()
    test_texture_hsg()
    test_sotwis_loader()
    test_csv_loader()
    test_catchment_block()
    test_pipeline_soil()
    print("\n" + "=" * 62)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S):")
        for f in FAILURES:
            print("   -", f)
        return 1
    print("ALL CHECKS PASSED")
    print("=" * 62)
    return 0


if __name__ == "__main__":
    sys.exit(main())
