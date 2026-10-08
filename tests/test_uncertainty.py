# -*- coding: utf-8 -*-
"""v0.27 validation: DEM-error sensitivity per crossing (F15).
Bare Python + NumPy, no QGIS, no GDAL.

    python -m qeht.tests.test_uncertainty
"""

import json
import math
import sys

import numpy as np

from ..core.watershed.uncertainty import (correlated_field, run, measure, params_json, PRESETS,
                                          UNC_FIELDS, upstream_lengths)

FAILURES = []
CS = 10.0
GT = (500000.0, CS, 0.0, 9900000.0, 0.0, -CS)


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def v_valley(rows=50, cols=41, side=0.5, along=0.02):
    r, c = np.mgrid[0:rows, 0:cols].astype(float)
    return 300.0 + side * np.abs(c - cols // 2) * CS - along * r * CS


def flat_tie(rows=50, cols=41):
    """Two valleys (cols 10 and 30) below a dead-flat plateau (rows 0-24) that drains to
    either one: which way it goes is decided by the error field."""
    r, c = np.mgrid[0:rows, 0:cols].astype(float)
    valleys = 200.0 + 0.3 * np.minimum(np.abs(c - 10), np.abs(c - 30)) * CS - 0.05 * (r - 25) * CS
    return np.where(r < 25, 210.0, valleys)


def test_field():
    print("\n1. Correlated error field")
    rng = np.random.default_rng(1)
    f = np.stack([correlated_field((200, 200), 2.5, 50.0, CS, rng) for _ in range(6)])
    sd = float(f.std())
    check("sample standard deviation within 5 % of sigma (2.5 m)", abs(sd / 2.5 - 1) < 0.05, f"{sd:.3f}")

    def corr(lag):
        a, b = f[:, :, :-lag].ravel(), f[:, :, lag:].ravel()
        return float(np.corrcoef(a, b)[0, 1])
    c5, c10, c15 = corr(5), corr(10), corr(15)
    check("correlation exp(-r^2/L^2): 1/e at L = 50 m (5 cells), e^-4 at 2 L, ~0 at 3 L",
          abs(c5 - math.exp(-1)) < 0.07 and abs(c10 - math.exp(-4)) < 0.05 and abs(c15) < 0.05,
          f"{c5:.3f}, {c10:.3f}, {c15:.3f}")
    check("sigma 0: zero field", not correlated_field((10, 10), 0.0, 90, CS, rng).any())


def test_runs():
    print("\n2. Monte Carlo at the crossings")
    z = v_valley()
    valid = np.ones(z.shape, bool)
    out_rc = [(45, 20)]
    det = measure(z, valid, CS, CS, out_rc, threshold_cells=20)[0]
    s0, _ = run(z, valid, GT, out_rc, 0.0, n=4, threshold_cells=20)
    b = s0[0]
    check("sigma 0: p10 = p50 = p90 = the deterministic value (area, LFP, slope, Tc)",
          all(math.isclose(b[f"{p}_{q}"], det[k]) for k, p in (("area", "unc_area"),
                                                              ("lfp_length_m", "unc_lfp"),
                                                              ("lfp_slope_1085", "unc_s1085"),
                                                              ("tc_kirpich_min", "unc_tc"))
              for q in ("p10", "p50", "p90")) and b["unc_area_cv"] == 0.0 and b["unc_lost_pct"] == 0.0)
    a1, _ = run(z, valid, GT, out_rc, 1.0, 30.0, n=6, seed=7, threshold_cells=20)
    a2, _ = run(z, valid, GT, out_rc, 1.0, 30.0, n=6, seed=7, threshold_cells=20)
    check("fixed seed: identical results", a1 == a2)
    check("keys = UNC_FIELDS", set(a1[0]) == set(UNC_FIELDS))
    sv, _ = run(z, valid, GT, out_rc, 0.5, 30.0, n=15, seed=3, threshold_cells=20)
    zf = flat_tie()
    sf, info = run(zf, np.ones(zf.shape, bool), GT, [(48, 10)], 0.5, 30.0, n=15, seed=3,
                   threshold_cells=20)
    # residual variation: the snapped outlet can move one cell along the channel (one row of
    # this small catchment is about 2 % of its area)
    check("steep V valley: area CV small (< 5 %)", sv[0]["unc_area_cv"] < 0.05,
          f"CV {sv[0]['unc_area_cv']:.4f}")
    check("flat plateau with a tie: area CV large (> 20 %) and switching reported",
          sf[0]["unc_area_cv"] > 0.2 and sf[0]["unc_switch_pct"] > 0 and "switching" in sf[0]["unc_note"],
          f"CV {sf[0]['unc_area_cv']:.3f}, switch {sf[0]['unc_switch_pct']:.0f} %")
    zl = v_valley()
    lost, _ = run(zl, np.ones(zl.shape, bool), GT, [(45, 2)], 0.5, 30.0, n=3, threshold_cells=20,
                  snap_cells=1)
    check("pour point far from any stream: lost in every realisation, noted",
          lost[0]["unc_lost_pct"] == 100.0 and lost[0]["unc_area_p50"] is None and lost[0]["unc_note"])
    pj = json.loads(params_json(info, "FABDEM", PRESETS["FABDEM"][1]))
    check("presets carry their sources; parameters JSON records N, sigma, L, seed, timing",
          PRESETS["AW3D30"][0] == 4.4 and "Tadono" in PRESETS["AW3D30"][1]
          and PRESETS["FABDEM"][0] == 2.5 and "Hawker" in PRESETS["FABDEM"][1]
          and "ESTIMATE" in PRESETS["FABDEM"][1]
          and pj["n"] == 15 and pj["seed"] == 3 and pj["sigma_m"] == 0.5 and "seconds_per_realisation" in pj)
    from ..core.interop.field_dictionary import field_names, METADATA_KEYS
    check("contract: unc fields on crossings, uncertainty_json in the metadata",
          all(f in field_names("crossings") for f in UNC_FIELDS)
          and "uncertainty_json" in [k for k, _ in METADATA_KEYS])


def test_lengths():
    print("\n3. One-pass longest flow lengths")
    from .test_interop import synthetic_case
    from ..core.conditioning.fill import fill_depressions
    from ..core.flow.direction import d8_direction
    from ..core.watershed.delineate import longest_flow_path
    z, v, gt, cs = synthetic_case()
    f, _, _ = fill_depressions(z, v, cell_width=cs, cell_height=cs)
    d, _ = d8_direction(f, v, cs, cs)
    L, fc = upstream_lengths(d, v, cs, cs)
    pts = [(55, 20), (20, 20), (25, 38), (59, 20), (40, 10)]
    ref = [longest_flow_path(d, v, rc, cell_width=cs, cell_height=cs) for rc in pts]
    check("longest flow length at any cell = longest_flow_path to that outlet (5 outlets)",
          all(math.isclose(a["length"], L[r * z.shape[1] + c]) for a, (r, c) in zip(ref, pts)))
    r, c = pts[0]
    k, n = r * z.shape[1] + c, 1
    while fc[k] >= 0:
        k, n = fc[k], n + 1
    check("traced path starts at a divide and has as many cells as the reference path",
          n == len(ref[0]["cells"]))


def main(argv=None):
    test_lengths()
    test_field()
    test_runs()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
