# -*- coding: utf-8 -*-
"""v0.26 validation: land-cover scenarios (F14).
Bare Python + NumPy, no QGIS, no GDAL.

    python -m qeht.tests.test_scenario
"""

import json
import math
import sys

import numpy as np

from ..core.runoff.curve_number import RunoffInputs, tr55_lookup
from ..core.runoff.scenario import Scenario, SCN_FIELDS
from ..core.erosion.rusle import WORLDCOVER_C as _WC, c_from_worldcover

WORLDCOVER_C = {k: (v[0] if isinstance(v, (tuple, list)) else float(v)) for k, v in _WC.items()}

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def base(classes, hsg_code=2, rc=None):
    hsg = np.full(classes.shape, hsg_code, np.int8)                 # 1..4 = A..D
    return RunoffInputs(classes, hsg=hsg, condition="fair", amc="III", rc_lookup=rc,
                        rc_lookup_path="user.csv" if rc else None)


def test_identical():
    print("\n1. Identical scenario")
    cl = np.where(np.arange(100).reshape(10, 10) % 3 == 0, 30.0, 40.0)
    b = base(cl, rc={30: (0.3,) * 4, 40: (0.5,) * 4})
    c, _ = c_from_worldcover(cl)
    a = 2.0 * c
    s = Scenario("same", cl.copy(), b, soil_loss=a, c_grid=c)
    out = s.block(np.ones(cl.shape, bool))
    check("identical land cover: d_cn = d_rational_c = d_ero_a = 0 exactly",
          out["d_cn"] == 0.0 and out["d_rational_c"] == 0.0 and out["d_ero_a_mean_tha"] == 0.0,
          str({k: out[k] for k in ("d_cn", "d_rational_c", "d_ero_a_mean_tha")}))
    check("keys = SCN_FIELDS, name carried", set(out) == set(SCN_FIELDS) and out["scenario_name"] == "same")


def test_grass_to_built():
    print("\n2. Grass to built-up on 30 % of the catchment")
    cl = np.full((10, 10), 30.0)
    sc = cl.copy()
    sc[:3, :] = 50.0
    lut = tr55_lookup("fair")
    b = base(cl)
    s = Scenario("development", sc, b)
    m = np.ones(cl.shape, bool)
    out = s.block(m)
    exp = 0.3 * (lut[50][1] - lut[30][1])
    check(f"d_cn = 0.3 x (CN built B {lut[50][1]:g} - CN grass B {lut[30][1]:g}) = {exp:.2f}",
          math.isclose(out["d_cn"], exp, abs_tol=1e-9), f"{out['d_cn']:.4f}")
    check("cn_ii_scn and shares: 70 % grass, 30 % built under the scenario",
          math.isclose(out["cn_ii_scn"], 0.7 * lut[30][1] + 0.3 * lut[50][1])
          and math.isclose(out["lc_pct_grass_scn"], 70.0) and math.isclose(out["lc_pct_built_scn"], 30.0))
    check("exported at the baseline AMC (III)", out["cn_export_scn"] > out["cn_ii_scn"])
    half = np.zeros(cl.shape, bool); half[:5] = True
    oh = s.block(half)
    check("per catchment: the upper half has 60 % built", math.isclose(oh["lc_pct_built_scn"], 60.0)
          and math.isclose(oh["d_cn"], 0.6 * (lut[50][1] - lut[30][1])))
    check("no user C lookup: Rational C stays empty", out["rational_c_scn"] is None and out["d_rational_c"] is None)


def test_erosion():
    print("\n3. Soil loss with the scenario C")
    cl = np.full((10, 10), 30.0)
    sc = cl.copy(); sc[:5] = 40.0
    c, _ = c_from_worldcover(cl)
    a = 3.0 * c                                                   # R K LS P = 3
    s = Scenario("crop", sc, base(cl), soil_loss=a, c_grid=c)
    out = s.block(np.ones(cl.shape, bool))
    exp_c = 0.5 * (WORLDCOVER_C[30] + WORLDCOVER_C[40])
    check("ero_c_mean_scn from the WorldCover lookup", math.isclose(out["ero_c_mean_scn"], exp_c))
    check("A_scn = A x C_scn / C: mean 3 x C_scn, change 3 x (C_scn - C)",
          math.isclose(out["ero_a_mean_tha_scn"], 3.0 * exp_c)
          and math.isclose(out["d_ero_a_mean_tha"], 3.0 * (exp_c - WORLDCOVER_C[30])))
    mj = json.loads(s.meta_json("C from WorldCover"))
    check("metadata: name, lookups, relation, land cover only", mj["name"] == "crop"
          and "C_scn" in mj["erosion"] and "land cover only" in mj["note"])
    from ..core.interop.field_dictionary import field_names, METADATA_KEYS
    check("contract: scenario fields on catchments, scenario_json in the metadata",
          all(f in field_names("catchments") for f in SCN_FIELDS)
          and "scenario_json" in [k for k, _ in METADATA_KEYS])


def main(argv=None):
    test_identical()
    test_grass_to_built()
    test_erosion()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
