# -*- coding: utf-8 -*-
"""v0.22 validation: time of concentration table (F10).
Bare Python + NumPy, no QGIS, no GDAL. Hand-computed values, unit
conversions against the US-customary forms, missing inputs and the
validity-range edges.

    python -m qeht.tests.test_tc
"""

import json
import math
import sys

from ..core.runoff.tc import (kirpich_min, kerby_min, scs_lag_min, tr55_sheet_h, tr55_shallow_h,
                              manning_v, bransby_williams_min, section_hydraulics, tc_block,
                              along_path_value, params_json, TC_FIELDS, SHEET_N, KERBY_N)

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def close(a, b, rel=1e-9):
    return a is not None and math.isclose(a, b, rel_tol=rel)


def test_formulas():
    print("\n1. Formulas: hand values and unit forms")
    v = kirpich_min(1000.0, 0.05)
    check("Kirpich L 1000 m, S 0.05: 0.0195 x 1000^0.77 x 0.05^-0.385 = 12.62 min",
          close(v, 0.0195 * 1000 ** 0.77 * 0.05 ** -0.385) and abs(v - 12.62) < 0.01, f"{v:.3f}")
    ft = 0.0078 * (1000.0 / 0.3048) ** 0.77 * 0.05 ** -0.385
    check("Kirpich metric = feet form 0.0078 L_ft^0.77 S^-0.385 (within 0.5 %)",
          abs(v / ft - 1) < 0.005, f"{v:.3f} vs {ft:.3f}")
    k = kerby_min(100.0, 0.4, 0.02)
    k_ft = 0.83 * (100.0 / 0.3048 * 0.4 / math.sqrt(0.02)) ** 0.467
    check("Kerby metric 1.44 (L N)^0.467 S^-0.235 = feet form 0.83 (L_ft N / S^0.5)^0.467 "
          "(within 1 %)", abs(k / k_ft - 1) < 0.01, f"{k:.3f} vs {k_ft:.3f}")
    s = scs_lag_min(2000.0, 75.0, 5.0)
    l_ft, S = 2000.0 / 0.3048, 1000.0 / 75.0 - 10.0
    hand = 60 * (l_ft ** 0.8 * (S + 1) ** 0.7 / (1900 * 5.0 ** 0.5)) / 0.6
    check("SCS lag L 2000 m (6,562 ft), CN 75, Y 5 %: Tc = lag / 0.6 = 74.33 min",
          close(s, hand) and abs(s - 74.33) < 0.01, f"{s:.2f}")
    t = 60 * tr55_sheet_h(0.15, 30.0, 88.9, 0.02)
    hand = 60 * 0.007 * (0.15 * 30.0 / 0.3048) ** 0.8 / (3.5 ** 0.5 * 0.02 ** 0.4)
    check("TR-55 sheet n 0.15, 30 m (98.4 ft), P2 88.9 mm (3.5 in), s 0.02: 9.25 min",
          close(t, hand) and abs(t - 9.25) < 0.01, f"{t:.3f}")
    sh = 60 * tr55_shallow_h(300.0, 0.02)
    check("TR-55 shallow unpaved 300 m at s 0.02: V 2.28 ft/s, 7.19 min",
          abs(sh - (300 / 0.3048) / (16.1345 * 0.02 ** 0.5) / 60) < 1e-9 and abs(sh - 7.19) < 0.01,
          f"{sh:.3f}")
    b = bransby_williams_min(5.0, 10.0, 10.0)
    check("Bransby-Williams L 5 km, A 10 km2, S 10 m/km: 73 / 10^0.3 = 36.59 min",
          abs(b - 73.0 / 10 ** 0.3) < 1e-9 and abs(b - 36.59) < 0.01, f"{b:.3f}")
    st = [[x, (0.0 if abs(x) <= 3 else (abs(x) - 3) / 2 if abs(x) <= 7 else 2 + (abs(x) - 7) * 0.05)]
          for x in [i * 0.5 for i in range(-60, 61)]]
    a, p = section_hydraulics(st, 2.0)
    check("bank-full section of the F5 trapezoid (6 m bed, 1:2, 2 m): A 20 m2, P 14.94 m",
          abs(a - 20.0) < 1e-9 and abs(p - (6 + 2 * math.hypot(4, 2))) < 1e-9, f"{a}, {p:.3f}")
    vm = manning_v(a, p, 0.01, 0.035)
    check("Manning V = R^2/3 s^1/2 / n", close(vm, (a / p) ** (2 / 3) * 0.1 / 0.035), f"{vm:.3f} m/s")


def inputs():
    fp = {"lfp_length_m": 1000.0, "lfp_slope_1085": 0.05, "lfp_overland_m": 200.0,
          "lfp_overland_slope": 0.04, "lfp_channel_m": 800.0, "lfp_channel_slope": 0.03}
    ca = {"area_km2": 0.3, "catch_slope_horn": 0.08, "cn_ii": 75.0}
    st = [[x, (0.0 if abs(x) <= 3 else (abs(x) - 3) / 2 if abs(x) <= 7 else 2 + (abs(x) - 7) * 0.05)]
          for x in [i * 0.5 for i in range(-60, 61)]]
    xs = {"xs_station_elev_json": json.dumps(st), "xs_bankfull_d_m": 2.0, "xs_quality": "high"}
    return fp, ca, xs


def test_block():
    print("\n2. Per-crossing block")
    fp, ca, xs = inputs()
    b = tc_block(fp, ca, xs, sheet_n=0.15, kerby_n=0.4, p2_mm=60.0)
    check("keys = TC_FIELDS", set(b) == set(TC_FIELDS))
    check("all five methods computed, none chosen",
          all(b[f"tc_{m}_min"] is not None for m in ("kirpich", "kerby_kirpich", "scs_lag", "tr55",
                                                      "bransby_williams")) and b["tc_note"] is None,
          str(b["tc_note"]))
    exp = kerby_min(200.0, 0.4, 0.04) + kirpich_min(800.0, 0.03)
    check("Kerby + Kirpich = overland Kerby + channel Kirpich", close(b["tc_kerby_kirpich_min"], exp))
    a, p = section_hydraulics(json.loads(xs["xs_station_elev_json"]), 2.0)
    exp = 60 * (tr55_sheet_h(0.15, 30.0, 60.0, 0.04) + tr55_shallow_h(170.0, 0.04)) \
        + 800.0 / manning_v(a, p, 0.03, 0.035) / 60.0
    check("TR-55 = sheet 30 m + shallow 170 m + channel 800 m (Manning on the section)",
          close(b["tc_tr55_min"], exp), f"{b['tc_tr55_min']:.2f}")
    bas = json.loads(b["tc_basis_json"])
    check("basis records the inputs per method (S basis, CN field, P2, n)",
          bas["kirpich"]["S_basis"] == "lfp_slope_1085" and bas["scs_lag"]["CN_field"] == "cn_ii"
          and bas["tr55"]["P2_mm"] == 60.0 and bas["tr55"]["channel_n"] == 0.035)
    m = tc_block(fp, ca, xs, sheet_n=0.15, kerby_n=0.4, p2_mm=None)
    check("missing P2: TR-55 NULL with the note 'give P2', others unchanged",
          m["tc_tr55_min"] is None and "give P2" in m["tc_note"]
          and close(m["tc_kirpich_min"], b["tc_kirpich_min"]))
    n = tc_block(fp, ca, None, sheet_n=None, kerby_n=None, p2_mm=60.0)
    check("no land cover and no section: Kerby and TR-55 NULL with reasons",
          n["tc_kerby_kirpich_min"] is None and n["tc_tr55_min"] is None
          and "retardance" in n["tc_note"] and "channel section" in n["tc_note"])
    c = tc_block(fp, {"area_km2": 0.3, "catch_slope_horn": 0.08}, xs, 0.15, 0.4, 60.0)
    check("no CN: SCS lag NULL with a note", c["tc_scs_lag_min"] is None and "curve number" in c["tc_note"])
    nc = dict(fp, lfp_overland_m=1000.0, lfp_channel_m=0.0, lfp_channel_slope=None)
    w = tc_block(nc, ca, None, 0.15, 0.4, 60.0)
    check("whole path overland (no channel): Kerby only overland; TR-55 needs no section",
          close(w["tc_kerby_kirpich_min"], kerby_min(1000.0, 0.4, 0.04))
          and w["tc_tr55_min"] is not None)


def test_flags():
    print("\n3. Validity flags switch at the range edges")
    fp, ca, xs = inputs()

    def flag(m, fp_=None, ca_=None):
        return tc_block(fp_ or fp, ca_ or ca, xs, 0.15, 0.4, 60.0)[f"tc_{m}_flag"]
    check("Kirpich area 0.45 within, 0.46 outside",
          flag("kirpich", ca_=dict(ca, area_km2=0.45)) == "within"
          and flag("kirpich", ca_=dict(ca, area_km2=0.46)).startswith("outside: area"))
    check("Kirpich slope 3 % within, 2.9 % outside; 10 % within, 10.1 % outside",
          flag("kirpich", fp_=dict(fp, lfp_slope_1085=0.03)) == "within"
          and "slope" in flag("kirpich", fp_=dict(fp, lfp_slope_1085=0.029))
          and flag("kirpich", fp_=dict(fp, lfp_slope_1085=0.10)) == "within"
          and "slope" in flag("kirpich", fp_=dict(fp, lfp_slope_1085=0.101)))
    check("SCS lag area 8 within, 8.1 outside; CN 95 within, 96 outside",
          flag("scs_lag", ca_=dict(ca, area_km2=8.0)) == "within"
          and "area" in flag("scs_lag", ca_=dict(ca, area_km2=8.1))
          and flag("scs_lag", ca_=dict(ca, cn_ii=95.0)) == "within"
          and "CN" in flag("scs_lag", ca_=dict(ca, cn_ii=96.0)))
    check("Kerby overland 365 m within, 366 m outside",
          flag("kerby_kirpich", fp_=dict(fp, lfp_overland_m=365.0)) == "within"
          and "overland" in flag("kerby_kirpich", fp_=dict(fp, lfp_overland_m=366.0)))
    check("Bransby-Williams always 'rural catchments'; TR-55 low section flagged",
          flag("bransby_williams") == "rural catchments"
          and tc_block(fp, ca, dict(xs, xs_quality="low"), 0.15, 0.4, 60.0)["tc_tr55_flag"]
          .startswith("outside"))


def test_lookups():
    print("\n4. Land-cover lookups along the path")
    coords = [(float(x), 0.0) for x in range(0, 201, 10)]          # 200 m, divide at x = 0

    def classes_at(xs, ys):
        return [30 if x < 50 else 10 for x in xs]
    check("sheet n over the first 30 m (grass): 0.15",
          close(along_path_value(coords, 30.0, classes_at, SHEET_N), 0.15))
    check("Kerby N over 100 m: 50 m grass 0.40 + 50 m tree 0.60 -> 0.50",
          close(along_path_value(coords, 100.0, classes_at, KERBY_N), 0.5))
    check("water only: None", along_path_value(coords, 30.0, lambda x, y: [80] * len(x), SHEET_N) is None)
    pj = json.loads(params_json(None, "", 0.035, "PROXY: x", "PROXY: y", "cn_ii"))
    check("params JSON: ranges with sources, P2 not given, PROXY lookups, no method chosen",
          "Kirpich (1940)" in pj["ranges"]["kirpich"]["source"] and "not given" in pj["p2_source"]
          and pj["sheet_n_lookup"].startswith("PROXY") and "does not choose" in pj["note"])
    from ..core.interop.field_dictionary import field_names, METADATA_KEYS
    check("contract: tc fields on crossings, tc_params_json in the metadata",
          all(f in field_names("crossings") for f in TC_FIELDS)
          and "tc_params_json" in [k for k, _ in METADATA_KEYS])


def main(argv=None):
    test_formulas()
    test_block()
    test_flags()
    test_lookups()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
