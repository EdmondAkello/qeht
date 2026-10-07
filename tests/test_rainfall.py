# -*- coding: utf-8 -*-
"""v0.17 validation: rainfall zone, mean annual rainfall and R estimate per
catchment (F3). Bare Python + NumPy, no QGIS, no GDAL.

    python -m qeht.tests.test_rainfall
"""

import json
import math
import sys

import numpy as np

from ..core.runoff.rainfall import (rain_block, dominant_zone, r_from_map, check_map,
                                    RainfallInputs, RAIN_FIELDS, R_RELATIONS)

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def test_zones():
    print("\n1. Zone shares, dominant zone and the tie rule")
    z = np.zeros((4, 4), np.int32)
    z[:, :2] = 1
    z[:, 2:] = 2
    names = {1: "Inland", 2: "Coastal"}
    b = rain_block(np.ones(z.shape, bool), z, names)
    shares = json.loads(b["rain_zones_json"])
    check("half-and-half zones give 50 / 50", shares == {"Coastal": 50.0, "Inland": 50.0},
          b["rain_zones_json"])
    check("tie -> the zone name sorting first (Coastal), share 50 %",
          b["rain_zone"] == "Coastal" and b["rain_zone_pct"] == 50.0)
    b2 = rain_block(np.ones(z.shape, bool), np.where(z == 1, 2, 1), {1: "Coastal", 2: "Inland"})
    check("tie result does not depend on zone order in the layer", b2["rain_zone"] == "Coastal")
    check("tie rule is case-insensitive", dominant_zone({"b": 2, "A": 2})[0] == "A")
    m = np.zeros(z.shape, bool)
    m[:, 1:] = True                          # 4 Inland + 8 Coastal cells
    b3 = rain_block(m, z, names)
    check("one third / two thirds", b3["rain_zone"] == "Coastal"
          and math.isclose(b3["rain_zone_pct"], 200 / 3), f"{b3['rain_zone_pct']:.2f}")
    z2 = z.copy()
    z2[0, :] = 0                              # a row outside every zone
    b4 = rain_block(np.ones(z.shape, bool), z2, names)
    check("cells outside every zone lower the coverage, not the shares",
          b4["rain_zone_coverage_pct"] == 75.0 and b4["rain_zone_pct"] == 50.0)
    z3 = np.ones((2, 2), np.int32)
    z3[0, 0] = 3                              # two polygons with the same name, two indices
    b5 = rain_block(np.ones((2, 2), bool), z3, {1: "Inland", 3: "Inland"})
    check("polygons sharing a name are one zone", b5["rain_zone"] == "Inland"
          and b5["rain_zone_pct"] == 100.0 and json.loads(b5["rain_zones_json"]) == {"Inland": 100.0})
    empty = rain_block(np.zeros(z.shape, bool), z, names)
    check("empty mask gives NULLs", all(v is None for v in empty.values()))


def test_qa():
    print("\n2. QA warning on mixed catchments")
    z = np.ones((10, 10), np.int32)
    z[:, :3] = 2                              # 70 % zone 1
    ri = RainfallInputs(z, {1: "Zone I", 2: "Zone II"}, "zones.gpkg field ZONE")
    b = ri.block(np.ones(z.shape, bool))
    check("70 % share flagged at the default 80 % threshold",
          b["rain_zone"] == "Zone I" and ri.qa_issue("X001", b) is not None)
    ri2 = RainfallInputs(z, {1: "Zone I", 2: "Zone II"}, zone_threshold_pct=60)
    check("not flagged at a 60 % threshold", ri2.qa_issue("X001", ri2.block(np.ones(z.shape, bool))) is None)
    check("no zones -> no QA warning", RainfallInputs(map_grid=np.full((2, 2), 900.0))
          .qa_issue("X", {"rain_zone_pct": None}) is None)


def test_map():
    print("\n3. Mean annual rainfall")
    g = np.full((5, 5), 812.5)
    b = rain_block(np.ones(g.shape, bool), map_grid=g, map_dataset="uniform")
    check("uniform raster gives the exact mean", b["map_mm"] == 812.5 and b["map_coverage_pct"] == 100.0)
    g2 = np.arange(25, dtype=float).reshape(5, 5) * 40 + 400
    g2[0, 0] = np.nan
    b2 = rain_block(np.ones(g2.shape, bool), map_grid=g2)
    exp = float(np.nanmean(g2))
    check("NoData cells are left out of the mean and lower the coverage",
          math.isclose(b2["map_mm"], exp) and b2["map_coverage_pct"] == 96.0, f"{b2['map_mm']:.2f}")
    for bad, label in ((np.full((2, 2), 3.2), "mm/day"), (np.full((2, 2), -5.0), "negative"),
                       (np.full((2, 2), 25000.0), "tenths of mm")):
        try:
            check_map(bad)
            ok = False
        except ValueError:
            ok = True
        check(f"implausible annual rainfall refused ({label})", ok)
    check("a 60-90 mm/yr raster passes with a warning", len(check_map(np.full((2, 2), 75.0))) == 1)
    check("a plausible raster passes silently", check_map(np.full((2, 2), 950.0)) == [])


def test_r():
    print("\n4. R from rainfall (estimate)")
    rf = lambda p: float(r_from_map(np.array([p]), "renard_freimund_1994")[0])
    check("Renard & Freimund, P = 500: 0.0483 x 500^1.61",
          math.isclose(rf(500.0), 0.0483 * 500.0 ** 1.61), f"{rf(500.0):.1f}")
    check("Renard & Freimund, P = 1000: 587.8 - 1219 + 4105 = 3473.8",
          math.isclose(rf(1000.0), 3473.8, rel_tol=1e-12), f"{rf(1000.0):.1f}")
    check("the two branches meet at 850 mm (within 0.3 %)",
          abs(rf(850.0) - rf(850.0001)) / rf(850.0) < 0.003, f"{rf(850.0):.0f} / {rf(850.0001):.0f}")
    lo = float(r_from_map(np.array([1000.0]), "lo_1985")[0])
    check("Lo et al., P = 1000: 38.46 + 3480 = 3518.46", math.isclose(lo, 3518.46), f"{lo:.2f}")
    p = np.array([[600.0, 1400.0]])
    ri = RainfallInputs(map_grid=p, map_dataset="test", r_relation="renard_freimund_1994")
    b = ri.block(np.ones(p.shape, bool))
    exp = (rf(600.0) + rf(1400.0)) / 2
    check("catchment R = mean of cell R, not R of the mean rainfall",
          math.isclose(b["rusle_r"], exp) and not math.isclose(b["rusle_r"], rf(1000.0)),
          f"{b['rusle_r']:.0f} vs R(mean P) {rf(1000.0):.0f}")
    check("rusle_r_method says ESTIMATE and names the relation",
          b["rusle_r_method"].startswith("ESTIMATE") and "Renard & Freimund" in b["rusle_r_method"])
    check("NaN rainfall gives NaN R", np.isnan(r_from_map(np.array([np.nan]), "lo_1985")[0]))
    ri2 = RainfallInputs(zone_grid=np.ones((2, 2), np.int32), zone_names={1: "A"},
                         r_relation="lo_1985")
    check("a relation without a rainfall raster warns and leaves R empty",
          ri2.r_grid is None and len(ri2.warnings) == 1
          and ri2.block(np.ones((2, 2), bool))["rusle_r"] is None)
    try:
        r_from_map(np.array([900.0]), "hurni")
        ok = False
    except ValueError:
        ok = True
    check("unknown relation refused", ok)
    meta = json.loads(ri.meta_json())
    check("metadata records relation, estimate flag and tie rule",
          meta["r_relation"] == "renard_freimund_1994" and meta["r_is_estimate"] is True
          and "ties" in meta["zone_tie_rule"] and set(R_RELATIONS) == {"renard_freimund_1994", "lo_1985"})
    check("block keys = RAIN_FIELDS", set(b) == set(RAIN_FIELDS))


def test_dictionary():
    print("\n5. Field dictionary and metadata (additive under qeht-heas-1)")
    from ..core.interop.field_dictionary import field_names, METADATA_KEYS, SCHEMA_VERSION
    names = field_names("catchments")
    check("every rainfall field is in the catchments contract", all(f in names for f in RAIN_FIELDS))
    check("rainfall_json is a metadata key", "rainfall_json" in [k for k, _ in METADATA_KEYS])
    check("schema stays qeht-heas-1", SCHEMA_VERSION == "qeht-heas-1")


def test_pipeline():
    print("\n6. Rainfall block on package catchments")
    from .test_interop import synthetic_run
    from ..core.interop.heas_exchange import build_exchange_records
    (cr, ca, fp, _, _), g = synthetic_run()
    shape = g["dem"].shape
    cols = np.arange(shape[1])[None, :] * np.ones(shape)
    zones = np.where(cols < shape[1] // 2, 1, 2).astype(np.int32)
    rain = 700.0 + 2.0 * cols                 # wetter to the east
    ri = RainfallInputs(zones, {1: "West", 2: "East"}, "synthetic", rain, "synthetic gradient",
                        "renard_freimund_1994")
    pts = [dict(x=a["outlet_x"], y=a["outlet_y"], fid=i, source_id=None,
                outlet_x=a["outlet_x"], outlet_y=a["outlet_y"]) for i, (_, a) in enumerate(cr)]
    _, ca2, _, issues, _ = build_exchange_records(g["direction"], g["valid"], g["accum"], g["dem"],
                                                  g["gt"], pts, snap_radius_cells=0, rainfall=ri)
    ok = all(c["rain_zone"] in ("West", "East") and 0 < c["rain_zone_pct"] <= 100
             and 700 <= c["map_mm"] <= 700 + 2 * shape[1] and c["rusle_r"] > 0 for _, c in ca2)
    check("every catchment has a zone, rainfall within the raster range and an R estimate", ok,
          ", ".join(f"{c['outlet_uid']} {c['rain_zone']} {c['rain_zone_pct']:.0f}% "
                    f"{c['map_mm']:.0f} mm" for _, c in ca2))
    mixed = [c for _, c in ca2 if c["rain_zone_pct"] < 80]
    check("mixed catchments (< 80 %) are reported as issues",
          len([i for i in issues if "rainfall zone" in i]) == len(mixed), f"{len(mixed)} mixed")
    _, ca3, _, _, _ = build_exchange_records(g["direction"], g["valid"], g["accum"], g["dem"],
                                             g["gt"], pts, snap_radius_cells=0)
    check("without rainfall inputs the catchments carry no rainfall values",
          all(c.get("rain_zone") is None and c.get("map_mm") is None for _, c in ca3))


def main(argv=None):
    test_zones()
    test_qa()
    test_map()
    test_r()
    test_dictionary()
    test_pipeline()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
