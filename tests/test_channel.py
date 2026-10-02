# -*- coding: utf-8 -*-
"""v0.15 validation: approach / exit channel slopes (A5), overland STI and the
crossing deposition indicator (STI advisory R1, R2). Bare Python + NumPy.

    python -m qeht.tests.test_channel
"""

import math
import sys

import numpy as np

from ..core.erosion.terrain import erosion_indices, sti_overland, channel_mask
from ..core.watershed.channel import (main_stem_upstream, downstream_path, channel_slopes,
                                      deposition_indicator)
from ..core.flow.direction import d8_direction
from ..core.flow.accumulation import flow_accumulation
from .test_crossings import valley_dem, route, gt_for

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def test_sti():
    print("\n1. STI = Moore & Burch LS with m = 0.6 on overland cells (R1)")
    rows, cols, cs = 30, 30, 10.0
    r, c = np.mgrid[0:rows, 0:cols].astype(float)
    dem = 200.0 - 0.1 * r * cs                 # plane, 10 % fall to the south
    v = np.ones((rows, cols), bool)
    d, _ = d8_direction(dem, v, cs, cs)
    a, _ = flow_accumulation(d, v)
    ind06 = erosion_indices(dem, v, a, cs, cs, m_exponent=0.6)
    ind = erosion_indices(dem, v, a, cs, cs)
    chan = channel_mask(a, v, 1e9)             # no channels
    sti = sti_overland(ind["spec_catch_area"], ind["tan_beta"], chan, cap_m=0)
    inner = np.s_[2:-2, 2:-2]
    check("identity: STI == LS(m = 0.6) without a cap",
          np.nanmax(np.abs(sti[inner] - ind06["ls"][inner])) < 1e-9)
    a_s = 50.0
    beta = math.atan(0.1)
    exact = (a_s / 22.13) ** 0.6 * (math.sin(beta) / 0.0896) ** 1.3
    one = sti_overland(np.array([[a_s]]), np.array([[0.1]]), np.zeros((1, 1), bool), cap_m=0)
    check("exact value: A_s 50 m, 10 % slope", abs(one[0, 0] - exact) < 1e-12, f"{one[0, 0]:.4f}")
    capped = sti_overland(np.array([[400.0]]), np.array([[0.1]]), np.zeros((1, 1), bool), cap_m=100)
    ref = sti_overland(np.array([[100.0]]), np.array([[0.1]]), np.zeros((1, 1), bool), cap_m=0)
    check("A_s capped at 100 m", abs(capped[0, 0] - ref[0, 0]) < 1e-12)
    chan2 = channel_mask(a, v, 10)
    sti2 = sti_overland(ind["spec_catch_area"], ind["tan_beta"], chan2)
    check("channel cells NoData, overland cells defined",
          np.isnan(sti2[chan2]).all() and np.isfinite(sti2[~chan2 & np.isfinite(ind['ls'])]).all())


def _straight_channel(slope_up, slope_dn, n=120, cs=10.0, brk=60):
    """A single column channel draining south (row increases), with a
    side wall so every cell drains into it; bed slope slope_up above the
    break row and slope_dn below."""
    rows, cols = n, 21
    z_bed = np.zeros(rows)
    for i in range(1, rows):
        z_bed[i] = z_bed[i - 1] - (slope_up if i <= brk else slope_dn) * cs
    z = np.tile(z_bed[:, None], (1, cols)) + 0.5 * np.abs(np.arange(cols) - 10)[None, :] * cs
    v = np.ones(z.shape, bool)
    d, _ = d8_direction(z, v, cs, cs)
    a, _ = flow_accumulation(d, v)
    return z, v, d, a


def test_channel_slopes():
    print("\n2. Approach / exit channel slopes (A5)")
    cs = 10.0
    z, v, d, a = _straight_channel(0.02, 0.02)
    s = channel_slopes(d, v, a, z, 60, 10, cs, cs, distance_m=200.0)
    check("constant 2 % channel: up and down slope exact, 200 m each",
          abs(s["ch_slope_us"] - 0.02) < 1e-12 and abs(s["ch_slope_ds"] - 0.02) < 1e-12
          and s["ch_len_us_m"] == 200.0 and s["ch_len_ds_m"] == 200.0, str(s))
    z2, v2, d2, a2 = _straight_channel(0.04, 0.01)
    s2 = channel_slopes(d2, v2, a2, z2, 60, 10, cs, cs, distance_m=200.0)
    check("slope break at the crossing: 4 % upstream, 1 % downstream",
          abs(s2["ch_slope_us"] - 0.04) < 1e-12 and abs(s2["ch_slope_ds"] - 0.01) < 1e-12,
          f"{s2['ch_slope_us']:.3f} / {s2['ch_slope_ds']:.3f}")
    up = main_stem_upstream(d, v, a, 60, 10, 1e9, cs, cs)
    check("main stem stays on the channel column up to the head, never a side tributary",
          all(cc == 10 for rr, cc, _ in up if rr > 0) and [p for p in up if p[0] == 0],
          f"{len(up)} cells")
    s3 = channel_slopes(d, v, a, z, 1, 10, cs, cs, distance_m=200.0)
    check("at the divide: shorter length reported", s3["ch_len_us_m"] < 200.0,
          f"{s3['ch_len_us_m']:.0f} m")
    dn = downstream_path(d, v, 110, 10, 500.0, cs, cs)
    check("downstream path stops at the grid edge", dn[-1][0] == 119 and dn[-1][2] == 90.0)


def test_deposition():
    print("\n3. Deposition indicator: change in stream power into the crossing (R2)")
    cs = 10.0
    path = [(60 - k, 10, k * cs) for k in range(0, 61)]      # 600 m straight up
    chan = np.ones((121, 21), bool)
    tb = np.full((121, 21), 0.02)
    spi = np.full((121, 21), 5.0)
    out = deposition_indicator(path, spi, tb, chan, 100, 500)
    check("uniform channel: ratio 1, neutral", out["ero_dep_ratio"] == 1.0
          and out["ero_dep_flag"] == "neutral")
    spi2 = spi.copy()
    spi2[50:60, 10] = 2.0                                     # 10-100 m upstream: low power
    out2 = deposition_indicator(path, spi2, tb, chan, 100, 500)
    check("capacity falls into the crossing: ratio 0.4, deposition-prone",
          abs(out2["ero_dep_ratio"] - 0.4) < 1e-12 and out2["ero_dep_flag"] == "deposition-prone")
    spi3 = spi.copy(); spi3[50:60, 10] = 8.0
    check("capacity rises: scour-prone",
          deposition_indicator(path, spi3, tb, chan, 100, 500)["ero_dep_flag"] == "scour-prone")
    z, v, d, a = _straight_channel(0.04, 0.01, brk=49)
    up = main_stem_upstream(d, v, a, 60, 10, 500.0, cs, cs)
    ind = erosion_indices(z, v, a, cs, cs)
    ch = channel_mask(a, v, 5)
    # constant A (as on a short reach of a large stream): SPI ratio = slope ratio
    spi_c = 1000.0 * np.maximum(ind["tan_beta"], 1e-3)
    o = deposition_indicator(up, spi_c, ind["tan_beta"], ch, 100, 500, elevation=z)
    near = [spi_c[a_, b_] for a_, b_, dd in up[1:] if dd <= 100]
    far = [spi_c[a_, b_] for a_, b_, dd in up[1:] if 100 < dd <= 500]
    check("slope break, constant A: ratio = median slope near / far < 1, deposition-prone",
          abs(o["ero_dep_ratio"] - np.median(near) / np.median(far)) < 1e-12
          and o["ero_dep_ratio"] < 0.7 and o["ero_dep_flag"] == "deposition-prone"
          and abs(o["ero_slope_app_far_pct"] - (39 * 4.0 + 1.0) / 40) < 1e-9
          and abs(o["ero_slope_app_near_pct"] - 1.0) < 1e-9,
          f"ratio {o['ero_dep_ratio']:.2f}, near {o['ero_slope_app_near_pct']:.2f} %, "
          f"far {o['ero_slope_app_far_pct']:.2f} %")
    o_real = deposition_indicator(up, np.exp(ind["ln_spi"]), ind["tan_beta"], ch, 100, 500)
    check("with the real A growth the ratio reflects both (area doubles here: > slope ratio)",
          o_real["ero_dep_ratio"] > o["ero_dep_ratio"], f"{o_real['ero_dep_ratio']:.2f}")
    short = deposition_indicator(path[:20], spi, tb, chan, 100, 500)
    check("fewer than far-reach channel cells upstream: NULL + note",
          short["ero_dep_ratio"] is None and "needs 500" in short["ero_dep_note"])
    chan2 = chan.copy(); chan2[30, 10] = False                 # channel ends 300 m up
    o2 = deposition_indicator(path, spi, tb, chan2, 100, 500)
    check("main stem leaving the channel network stops the reach",
          o2["ero_dep_ratio"] is None and "290 m" in o2["ero_dep_note"])


def main(argv=None):
    test_sti()
    test_channel_slopes()
    test_deposition()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
