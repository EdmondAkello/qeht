# -*- coding: utf-8 -*-
"""v0.15 validation: soils from any source (G1-G4). Bare Python + NumPy.

    python -m qeht.tests.test_soils_any
"""

import math
import os
import random
import sys
import tempfile

import numpy as np

from ..core.geometry.rasterize import rasterize_polygons
from ..core.soils.sotwis import load_attributes, load_sotwis, load_csv_units
from ..core.soils.catchment import soil_block, SOIL_FIELDS
from ..core.soils.usle_k import (williams_k, dg_k, usda_texture, hsg_proxy, williams_k_arr,
                                 dg_k_arr, usda_texture_arr, hsg_proxy_arr, TEXTURES, HSG_CODES)
from ..core.soils.sources import (SoilGrid, hsg_from_codes, soilgrids_layers, depth_weighted,
                                  SOILGRIDS_FACTORS)
from .test_soils import _make_db

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def _strips():
    gt = (0.0, 10.0, 0.0, 100.0, 0.0, -10.0)
    shape = (10, 30)
    strips = [([[(0, 0), (100, 0), (100, 100), (0, 100)]], 1),
              ([[(100, 0), (200, 0), (200, 100), (100, 100)]], 2),
              ([[(200, 0), (280, 0), (280, 100), (200, 100)]], 3)]
    return gt, shape, rasterize_polygons(strips, gt, shape)


def test_vectorised_equal_scalar():
    print("\n1. Vectorised K, texture and HSG proxy = the scalar functions")
    rnd = random.Random(7)
    bad = 0
    for _ in range(5000):
        sa = rnd.uniform(0, 100); cl = rnd.uniform(0, 100 - sa); si = 100 - sa - cl
        oc = rnd.uniform(0, 6)
        t = usda_texture(sa, si, cl)
        tc = int(usda_texture_arr(sa, si, cl))
        k = williams_k(sa, si, cl, oc)[1]
        ka = float(williams_k_arr(sa, si, cl, oc))
        bad += TEXTURES[tc - 1] != t
        bad += not ((math.isnan(k) and math.isnan(ka)) or abs(k - ka) < 1e-12)
        bad += abs(dg_k(sa, si, cl) - float(dg_k_arr(sa, si, cl))) > 1e-12
        for dr in ("", "P", "I", "W", "V"):
            bad += hsg_proxy(t, dr) != HSG_CODES[int(hsg_proxy_arr(tc, np.array(dr))) - 1]
    check("5000 random textures x 5 drainage classes: no mismatch", bad == 0, f"{bad} mismatches")
    check("missing input -> NaN K, texture 0",
          math.isnan(float(williams_k_arr(np.nan, 30, 30, 1))) and int(usda_texture_arr(-1, 30, 30)) == 0)


NEW_FIELDS = {"soil_hsg", "hsg_pct_a", "hsg_pct_b", "hsg_pct_c", "hsg_pct_d", "hsg_pct_dual",
              "soil_hsg_source", "usle_k_source"}


def test_units_regression():
    print("\n2. Unit sources on the per-cell model reproduce the v0.11-0.14 block")
    gt, shape, grid = _strips()
    units, info = load_attributes([
        {"id": "A", "sand": 80, "silt": 10, "clay": 10, "oc": 1.0, "drain": "W"},
        {"id": "B", "sand": 20, "silt": 60, "clay": 20, "oc": 2.0, "drain": "M"},
        {"id": "C", "sand": 30, "silt": 30, "clay": 40, "oc": 1.5, "drain": "P"}], drain="drain")
    idx = {1: "A", 2: "B", 3: "C"}
    sg = SoilGrid.from_units(grid, idx, units, info)
    worst = 0.0
    same = True
    for m in range(6):
        mask = np.zeros(shape, bool)
        mask[:, m * 5:] = True
        mask[: m + 1, :] = False if m % 2 else mask[: m + 1, :]
        old = soil_block(mask, grid, units, idx, info)
        new = sg.block(mask)
        for k, v in old.items():
            if k in NEW_FIELDS:
                continue
            if isinstance(v, float) and isinstance(new[k], float):
                worst = max(worst, abs(v - new[k]) / max(1.0, abs(v)))
            elif v != new[k]:
                same = False
                print("     differs:", k, v, new[k])
    check("every v0.14 field identical (floats to 1e-12)", same and worst < 1e-12,
          f"max rel diff {worst:.1e}")
    tmp = tempfile.mkdtemp()
    db = os.path.join(tmp, "s.db")
    _make_db(db)
    su, si = load_sotwis(db, 0, 20)
    keys = sorted(su)
    g2 = np.zeros(shape, np.int32)
    g2[:, :15] = 1; g2[:, 15:] = 2
    idx2 = {1: keys[0], 2: keys[-1]}
    old = soil_block(np.ones(shape, bool), g2, su, idx2, si)
    new = SoilGrid.from_units(g2, idx2, su, si).block(np.ones(shape, bool))
    check("SOTWIS database units: identical block",
          all((abs(old[k] - new[k]) < 1e-12) if isinstance(old[k], float) else old[k] == new[k]
              for k in old if k not in NEW_FIELDS),
          f"K {new['usle_k']}, HSG {new['soil_hsg']} ({new['soil_hsg_source']})")
    check("block has exactly the SOIL_FIELDS keys", set(new) == {f for f, _ in SOIL_FIELDS})


def test_texture_rasters():
    print("\n3. Texture rasters: SoilGrids units, depth weighting, per-cell K")
    shape = (4, 5)
    # SoilGrids-stored values for a uniform loam: sand 400 g/kg, clay 200, silt 400, soc 150 dg/kg
    sand = np.full(shape, 400.0) / SOILGRIDS_FACTORS["sand"]
    clay = np.full(shape, 200.0) / SOILGRIDS_FACTORS["clay"]
    oc = np.full(shape, 150.0) / SOILGRIDS_FACTORS["soc"]
    check("SoilGrids conversion: 400 g/kg -> 40 %, 150 dg/kg -> 1.5 %",
          sand[0, 0] == 40.0 and abs(oc[0, 0] - 1.5) < 1e-12)
    sg = SoilGrid.from_texture(sand, clay, oc, info={"soil_dataset": "test", "soil_depth_cm": "0-20"})
    b = sg.block(np.ones(shape, bool))
    k_ref = williams_k(40, 40, 20, 1.5)[1]
    check("silt = 100 - sand - clay, K = scalar Williams/EPIC",
          abs(b["soil_silt_pct"] - 40) < 1e-12 and abs(b["usle_k"] - k_ref) < 1e-12,
          f"K {b['usle_k']:.5f}")
    check("texture loam, HSG proxy B on every cell, 100 % B",
          b["soil_texture"] == "loam" and b["soil_hsg"] == "B" and b["hsg_pct_b"] == 100.0
          and b["soil_units"] is None and b["soil_coverage_pct"] == 100.0)
    # depth weighting 0-5, 5-15, 15-30 -> 0-20: weights 5, 10, 5
    tmp = tempfile.mkdtemp()
    vals = {"sand_0-5cm_mean.tif": 300.0, "sand_5-15cm_mean.tif": 400.0,
            "sand_15-30cm_mean.tif": 600.0, "sand_15-30cm_Q0.5.tif": 0.0,
            "clay_0-5cm_mean.tif": 200.0, "readme.txt": 0}
    for name in vals:
        open(os.path.join(tmp, name), "w").close()
    lay = soilgrids_layers(tmp)
    check("folder scan: properties, depths, 'mean' preferred",
          sorted(lay) == ["clay", "sand"] and [(t, b_) for t, b_, _ in lay["sand"]] ==
          [(0, 5), (5, 15), (15, 30)] and lay["sand"][2][2].endswith("mean.tif"))
    read = lambda p: np.full(shape, vals[os.path.basename(p)])
    a, used = depth_weighted(lay["sand"], 0, 20, read)
    check("0-20 cm = (5*300 + 10*400 + 5*600) / 20 = 425",
          np.allclose(a, 425.0) and used == ["0-5", "5-15", "15-30"], f"{a[0, 0]}")
    hole = {"sand_0-5cm_mean.tif": np.where(np.eye(4, 5) > 0, np.nan, 300.0)}
    read2 = lambda p: hole.get(os.path.basename(p), np.full(shape, vals[os.path.basename(p)]))
    a2, _ = depth_weighted(lay["sand"], 0, 20, read2)
    check("NoData in one layer: averaged over the layers present",
          abs(a2[0, 0] - (10 * 400 + 5 * 600) / 15) < 1e-9 and abs(a2[0, 1] - 425) < 1e-9)


def test_hsg_and_k_overrides():
    print("\n4. HSG and K overrides (HYSOGs250m codes, letters, precedence)")
    codes = np.array([[1, 2, 3, 4, 11], [12, 13, 14, 255, np.nan]], dtype=float)
    h = hsg_from_codes(codes)
    check("HYSOGs codes: 1-4 -> A-D, 11-14 -> dual (14), 255/NaN -> unknown",
          h.tolist() == [[1, 2, 3, 4, 14], [14, 14, 14, 0, 0]])
    lt = hsg_from_codes(np.array(["a", "B ", "C/D", None, "x"], dtype=object), "letters")
    check("letters: A, B, C/D, missing, unknown", lt.tolist() == [1, 2, 14, 0, 0])
    shape = (2, 5)
    sg = SoilGrid.from_texture(np.full(shape, 40.0), np.full(shape, 20.0), np.full(shape, 1.5),
                               info={"soil_dataset": "t"})
    sg.override_hsg(h, "HYSOGs250m test")
    b = sg.block(np.ones(shape, bool))
    # known: 8 direct + 2 proxy (B) cells: A1 B1+2 C1 D1+4dual
    check("shares over known cells; dual counted in D and reported",
          math.isclose(b["hsg_pct_a"], 10) and math.isclose(b["hsg_pct_b"], 30)
          and math.isclose(b["hsg_pct_d"], 50) and math.isclose(b["hsg_pct_dual"], 40)
          and b["soil_hsg"] == "D" and b["soil_hsg_source"] == "HYSOGs250m test",
          f"A {b['hsg_pct_a']:.0f} B {b['hsg_pct_b']:.0f} C {b['hsg_pct_c']:.0f} "
          f"D {b['hsg_pct_d']:.0f} (dual {b['hsg_pct_dual']:.0f})")
    k_before = sg.block(np.ones(shape, bool))["usle_k"]
    kk = np.full(shape, np.nan); kk[0, :] = 0.05
    sg.override_k(kk, "K raster test")
    b2 = sg.block(np.ones(shape, bool))
    check("direct K where given, computed K elsewhere; source tagged",
          math.isclose(b2["usle_k"], (5 * 0.05 + 5 * k_before) / 10) and
          b2["usle_k_source"] == "K raster test")
    row1 = np.zeros(shape, bool); row1[1, :] = True
    check("catchment without direct K says so",
          "no direct K" in sg.block(row1)["usle_k_source"])
    only = SoilGrid(shape, {"soil_dataset": "none (HSG / K only)"})
    only.override_hsg(h, "HYSOGs250m test")
    b3 = only.block(np.ones(shape, bool))
    check("HSG raster alone: shares filled, texture and K empty, coverage 0",
          b3["soil_hsg"] == "D" and b3["usle_k"] is None and b3["soil_coverage_pct"] == 0)


def test_unit_raster_csv():
    print("\n5. Unit-code raster + CSV table (e.g. HWSD v2 mapping units)")
    tmp = tempfile.mkdtemp()
    p = os.path.join(tmp, "hwsd.csv")
    with open(p, "w", encoding="utf-8") as f:
        f.write("HWSD2_SMU_ID;sand;silt;clay;oc\n11000;45;35;20;1.2\n11001;20;30;50;0.8\n")
    units, info = load_csv_units(p, "HWSD2_SMU_ID")
    codes = np.array([[11000, 11000, 11001], [11001, np.nan, 99999]], dtype=float)
    keys = {"11000": 1, "11001": 2, "99999": 3}
    grid = np.zeros(codes.shape, np.int32)
    for k, i in keys.items():
        grid[codes == float(k)] = i
    sg = SoilGrid.from_units(grid, {i: k for k, i in keys.items()}, units, info)
    b = sg.block(np.ones(codes.shape, bool))
    exp = (2 * units["11000"].k_si + 2 * units["11001"].k_si) / 4
    check("codes joined as text, unknown code and NoData = no data",
          math.isclose(b["usle_k"], exp) and math.isclose(b["soil_coverage_pct"], 4 / 6 * 100)
          and b["soil_units"] == 2, f"K {b['usle_k']:.5f}, coverage {b['soil_coverage_pct']:.1f} %")


def main(argv=None):
    test_vectorised_equal_scalar()
    test_units_regression()
    test_texture_rasters()
    test_hsg_and_k_overrides()
    test_unit_raster_csv()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
