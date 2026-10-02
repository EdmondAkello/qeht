# -*- coding: utf-8 -*-
"""v0.15 validation: curve number and Rational C per catchment (F1).
Bare Python + NumPy, no QGIS, no GDAL.

    python -m qeht.tests.test_runoff
"""

import math
import os
import sys
import tempfile

import numpy as np

from ..core.runoff.curve_number import (tr55_lookup, value_grid, amc_convert, runoff_block,
                                        RunoffInputs, read_lookup_csv, RUNOFF_FIELDS)

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def test_mosaic():
    print("\n1. 2x2 cover / HSG mosaic -> exact area-weighted CN")
    # quadrants of a 4x4 grid: tree/B, grass/C, crop/D, built/A
    lc = np.array([[10, 10, 30, 30], [10, 10, 30, 30], [40, 40, 50, 50], [40, 40, 50, 50]], float)
    hsg = np.array([[2, 2, 3, 3], [2, 2, 3, 3], [4, 4, 1, 1], [4, 4, 1, 1]], np.int8)
    t = tr55_lookup("fair")
    cn = value_grid(lc, hsg, t)
    b = runoff_block(np.ones(lc.shape, bool), lc, cn, amc="II", lookup_name="tr55")
    exp = (60 + 79 + 89 + 89) / 4.0
    check("CN = (woods-B 60 + pasture-C 79 + row crops-D 89 + commercial-A 89) / 4",
          math.isclose(b["cn_ii"], exp) and b["cn_coverage_pct"] == 100.0, f"{b['cn_ii']:.2f}")
    check("land-cover shares 25 % each", all(math.isclose(b[f"lc_pct_{g}"], 25.0)
                                             for g in ("tree", "grass", "crop", "built")))
    m = np.zeros(lc.shape, bool); m[:2, :] = True
    check("half mask: (60 + 79) / 2", math.isclose(runoff_block(m, lc, cn)["cn_ii"], 69.5))
    good = runoff_block(np.ones(lc.shape, bool), lc, value_grid(lc, hsg, tr55_lookup("good")))
    poor = runoff_block(np.ones(lc.shape, bool), lc, value_grid(lc, hsg, tr55_lookup("poor")))
    check("condition: good < fair < poor",
          good["cn_ii"] < b["cn_ii"] < poor["cn_ii"],
          f"{good['cn_ii']:.2f} < {b['cn_ii']:.2f} < {poor['cn_ii']:.2f}")


def test_unknowns_and_dual():
    print("\n2. Dual groups, unknown HSG and classes without a CN")
    lc = np.array([[30, 30, 90, 0]], float)
    hsg = np.array([[14, 0, 2, 2]], np.int8)
    cn = value_grid(lc, hsg, tr55_lookup("fair"))
    check("dual group (B/D) used as D: grassland-D 84",
          cn[0, 0] == 84 and np.isnan(cn[0, 1]) and np.isnan(cn[0, 2]) and np.isnan(cn[0, 3]))
    b = runoff_block(np.ones(lc.shape, bool), lc, cn)
    check("cells without a CN are reported, not invented (coverage 25 %)",
          b["cn_ii"] == 84 and b["cn_coverage_pct"] == 25.0)
    check("wetland share counts only classified cells (1 of 3 = 33.3 %)",
          math.isclose(b["lc_pct_wetland"], 100 / 3) and math.isclose(b["lc_pct_grass"], 200 / 3))


def test_amc():
    print("\n3. AMC conversions (Chow, Maidment & Mays 1988)")
    check("CN_II 70 -> CN_III 1610/19.1 = 84.29 and CN_I 294/5.94 = 49.49",
          abs(amc_convert(70, "III") - 1610 / 19.1) < 1e-9 and abs(amc_convert(70, "I") - 294 / 5.94) < 1e-9,
          f"III {amc_convert(70, 'III'):.2f}, I {amc_convert(70, 'I'):.2f}")
    check("AMC II unchanged; 100 stays 100 in III",
          amc_convert(73.5, "II") == 73.5 and abs(amc_convert(100, "III") - 100) < 1e-9)
    ro = RunoffInputs(np.full((2, 2), 30.0), hsg=np.full((2, 2), 3, np.int8), amc="III")
    b = ro.block(np.ones((2, 2), bool))
    check("exported at AMC III from the composite AMC II value",
          b["cn_ii"] == 79 and b["cn_amc"] == "III" and math.isclose(b["cn_export"], amc_convert(79, "III")))


def test_lookups():
    print("\n4. User lookups (CN replacement, Rational C) and provenance")
    tmp = tempfile.mkdtemp()
    p = os.path.join(tmp, "cn.csv")
    with open(p, "w", encoding="utf-8") as f:
        f.write("class;name;A;B;C;D\n30;grass local;50;70;80;85\n")
    lut = read_lookup_csv(p)
    lc = np.full((2, 2), 30.0)
    hsg = np.array([[1, 2], [3, 4]], np.int8)
    ro = RunoffInputs(lc, hsg=hsg, cn_lookup=lut, cn_lookup_path=p)
    b = ro.block(np.ones((2, 2), bool))
    check("user CN lookup replaces TR-55 and is recorded",
          math.isclose(b["cn_ii"], (50 + 70 + 80 + 85) / 4) and p in b["cn_lookup_id"]
          and ro.meta["cn_proxy"] is False)
    q = os.path.join(tmp, "rc.csv")
    with open(q, "w", encoding="utf-8") as f:
        f.write("class,C\n30,0.35\n40,0.5\n")
    rc = read_lookup_csv(q, "C")
    ro2 = RunoffInputs(np.array([[30, 40]], float), hsg=None, rc_lookup=rc, rc_lookup_path=q)
    b2 = ro2.block(np.ones((1, 2), bool))
    check("one-value Rational C per class works without soils",
          math.isclose(b2["rational_c"], 0.425) and b2["cn_ii"] is None)
    ro3 = RunoffInputs(lc, hsg=None, rc_lookup=lut, rc_lookup_path=p)
    check("per-HSG Rational C without soils: left empty with a note",
          ro3.rc is None and "soil" in ro3.note)
    d = RunoffInputs(lc, hsg=hsg)
    check("default lookup flagged as a TR-55 proxy in the metadata",
          d.meta["cn_proxy"] is True and "TR-55" in d.meta["cn_lookup_id"]
          and '"30": [49, 69, 79, 84]' in d.meta_json().replace("(", "[").replace(")", "]"))
    check("no Rational C default ships", d.rc is None and d.block(np.ones((2, 2), bool))["rational_c"] is None)
    check("block keys = RUNOFF_FIELDS", set(d.block(np.ones((2, 2), bool))) == set(RUNOFF_FIELDS))


def test_pipeline():
    print("\n5. CN block on package catchments")
    from .test_interop import synthetic_run
    from ..core.interop.heas_exchange import build_exchange_records
    (cr, ca, fp, _, _), g = synthetic_run()
    shape = g["dem"].shape
    lc = np.where(np.arange(shape[1])[None, :] < shape[1] // 2, 10.0, 40.0) * np.ones(shape)
    ro = RunoffInputs(lc, hsg=np.full(shape, 2, np.int8))
    pts = [dict(x=a["outlet_x"], y=a["outlet_y"], fid=i, source_id=None,
                outlet_x=a["outlet_x"], outlet_y=a["outlet_y"]) for i, (_, a) in enumerate(cr)]
    _, ca2, _, _, _ = build_exchange_records(g["direction"], g["valid"], g["accum"], g["dem"],
                                             g["gt"], pts, snap_radius_cells=0, runoff=ro)
    ok = all(c["cn_ii"] is not None and 60 - 1e-9 <= c["cn_ii"] <= 78 + 1e-9
             and math.isclose(c["lc_pct_tree"] + c["lc_pct_crop"], 100.0) for _, c in ca2)
    check("every catchment between woods-B 60 and row crops-B 78, shares add to 100", ok,
          ", ".join(f"{c['outlet_uid']} CN {c['cn_ii']:.1f}" for _, c in ca2))


def main(argv=None):
    test_mosaic()
    test_unknowns_and_dual()
    test_amc()
    test_lookups()
    test_pipeline()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
