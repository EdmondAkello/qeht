# -*- coding: utf-8 -*-
"""v0.18 validation: floodplain width indicator at major crossings (A4).
Bare Python + NumPy, no QGIS, no GDAL.

    python -m qeht.tests.test_floodplain
"""

import math
import sys

import numpy as np

from ..core.network.floodplain import (profile_widths, hand_grid, hand_widths, floodplain_block,
                                       FP_FIELDS)
from ..core.flow.direction import d8_direction

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def valley(ch, bed_from=480.0, bed_to=520.0, z_bed=100.0, s_left=0.05, s_right=0.05):
    ch = np.asarray(ch, float)
    return np.where(ch < bed_from, z_bed + (bed_from - ch) * s_left,
                    np.where(ch > bed_to, z_bed + (ch - bed_to) * s_right, z_bed))


def test_profile():
    print("\n1. Profile method on analytic valleys")
    ch = np.arange(0.0, 1000.1, 10.0)
    out, note = profile_widths(ch, valley(ch), 500.0)
    exp = {0.5: 60.0, 1.0: 80.0, 2.0: 120.0}               # 40 + 2 dz / 0.05
    check("trapezoid 40 m bed, 1:20 banks: 60 / 80 / 120 m",
          all(math.isclose(out[f"fp_w_{k:.1f}_m".replace(".", "p")], v) for k, v in exp.items())
          and out["fp_z_bed_m"] == 100.0 and not note,
          f"{out['fp_w_0p5_m']}, {out['fp_w_1p0_m']}, {out['fp_w_2p0_m']}")
    check("extent at bed + 1 m: ch 460 to 540", out["fp_ch_from_m"] == 460.0 and out["fp_ch_to_m"] == 540.0)
    ch7 = np.arange(3.0, 1000.0, 7.0)                       # stations off the breakpoints
    o7, _ = profile_widths(ch7, valley(ch7), 500.0)
    check("stations off the bank breakpoints: within one station of the analytic 80 m",
          abs(o7["fp_w_1p0_m"] - 80.0) <= 7.0, f"{o7['fp_w_1p0_m']:.2f}")
    a, _ = profile_widths(ch, valley(ch, s_left=0.05, s_right=0.1), 500.0)
    check("asymmetric banks 1:20 / 1:10: 40 + 20 + 10 = 70 m at bed + 1 m",
          math.isclose(a["fp_w_1p0_m"], 70.0) and math.isclose(a["fp_w_2p0_m"], 100.0))
    v, _ = profile_widths(ch, valley(ch, 500.0, 500.0, s_left=0.1, s_right=0.1), 500.0)
    check("no floodplain (V, banks 1:10): width = 2 dz / 0.1 (10 / 20 / 40 m)",
          [round(v[k], 9) for k in ("fp_w_0p5_m", "fp_w_1p0_m", "fp_w_2p0_m")] == [10.0, 20.0, 40.0])
    b, _ = profile_widths(ch, valley(ch), 470.0)
    check("crossing 30 m off the low point: the bed search (50 m) finds the same floor",
          math.isclose(b["fp_w_1p0_m"], 80.0))
    t, note = profile_widths(ch[ch <= 530], valley(ch[ch <= 530]), 500.0)
    check("profile ending inside the floodplain: lower bound and a note",
          math.isclose(t["fp_w_1p0_m"], 530.0 - 460.0) and "end of the profile" in note)
    z = valley(ch)
    z[ch == 560] = np.nan
    n, note = profile_widths(ch, z, 500.0)
    check("NoData station at ch 560 stops the bed + 2 m run at ch 550 (lower bound, noted); "
          "bed + 1 m ends before it",
          math.isclose(n["fp_w_2p0_m"], 550.0 - 440.0) and math.isclose(n["fp_w_1p0_m"], 80.0)
          and "dz 2" in note and "dz 1 " not in note)
    e, note = profile_widths(ch, valley(ch), 5000.0)
    check("no profile near the crossing: empty with a reason", not e and "no ground profile" in note)


def test_hand():
    print("\n2. Height above nearest drainage")
    rows, cols, cs = 40, 41, 30.0
    c0 = 20
    cc = np.arange(cols)[None, :] * np.ones((rows, 1))
    rr = np.arange(rows)[:, None] * np.ones((1, cols))
    dem = 200.0 + np.abs(cc - c0) * cs * 0.1 - rr * cs * 0.001     # 1:10 sides, 1:1000 along
    valid = np.ones(dem.shape, bool)
    direction, _ = d8_direction(dem, valid, cs, cs)
    stream = np.zeros(dem.shape, bool)
    stream[:, c0] = True
    h = hand_grid(direction, valid, dem, stream)
    exp = np.abs(cc - c0) * cs * 0.1
    check("V valley: HAND = lateral height above the stream (cells drain sideways)",
          np.allclose(h[1:-1, :], exp[1:-1, :]), f"max error {np.nanmax(np.abs(h - exp)[1:-1]):.2e}")
    check("stream cells have HAND 0", np.all(h[:, c0] == 0.0))
    ramp = (np.arange(500, dtype=float)[::-1] * 0.5)[None, :]
    d2, _ = d8_direction(ramp, np.ones(ramp.shape, bool), 10.0, 10.0)
    s2 = np.zeros(ramp.shape, bool); s2[0, -1] = True
    h2 = hand_grid(d2, np.ones(ramp.shape, bool), ramp, s2)
    check("500-cell path (pointer jumping): HAND at the head = total drop 249.5 m",
          math.isclose(h2[0, 0], 249.5))
    pit = np.ones((3, 3)) * 10.0; pit[1, 1] = 5.0
    d3, _ = d8_direction(pit, np.ones(pit.shape, bool), 1.0, 1.0)
    h3 = hand_grid(d3, np.ones(pit.shape, bool), pit, np.zeros(pit.shape, bool))
    check("no stream reached: NaN, not zero", np.all(np.isnan(h3)))
    ch = np.arange(0.0, 1000.1, 10.0)
    hw = hand_widths(ch, np.abs(ch - 500.0) * 0.1, 500.0)
    check("HAND widths along a road: 2 dz / 0.1 (10 / 20 / 40 m)",
          [round(hw[k], 9) for k in ("fp_hand_w_0p5_m", "fp_hand_w_1p0_m", "fp_hand_w_2p0_m")]
          == [10.0, 20.0, 40.0])


def test_block():
    print("\n3. Per-crossing block")
    ch = np.arange(0.0, 1000.1, 10.0)
    prof = [{"chainage_m": c, "z_dem_m": z} for c, z in zip(ch, valley(ch))]
    crs = [{"chainage_m": 500.0, "acc_at_outlet_km2": 25.0},
           {"chainage_m": 200.0, "acc_at_outlet_km2": 2.0},
           {"chainage_m": None, "acc_at_outlet_km2": 50.0}]
    b = floodplain_block(crs, prof, hand_at_stations=np.abs(ch - 500.0) * 0.05)
    check("large crossing: profile and HAND widths, method profile+HAND",
          math.isclose(b[0]["fp_w_1p0_m"], 80.0) and math.isclose(b[0]["fp_hand_w_1p0_m"], 40.0)
          and b[0]["fp_method"] == "profile+HAND" and b[0]["fp_note"] is None)
    check("below the area limit or without chainage: all empty",
          all(v is None for v in b[1].values()) and all(v is None for v in b[2].values()))
    b2 = floodplain_block(crs[:1], list(reversed(prof)))
    check("profile order does not matter; without HAND the method is 'profile'",
          math.isclose(b2[0]["fp_w_1p0_m"], 80.0) and b2[0]["fp_method"] == "profile"
          and b2[0]["fp_hand_w_1p0_m"] is None)
    check("keys = FP_FIELDS", set(b[0]) == set(FP_FIELDS))
    from ..core.interop.field_dictionary import field_names, METADATA_KEYS
    check("fields in the crossings contract, fp_params_json in the metadata",
          all(f in field_names("crossings") for f in FP_FIELDS)
          and "fp_params_json" in [k for k, _ in METADATA_KEYS])


def main(argv=None):
    test_profile()
    test_hand()
    test_block()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
