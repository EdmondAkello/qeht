# -*- coding: utf-8 -*-
"""v0.14 WP-F erosion: analytic tests and A14 (Project B) formula anchors.

    python -m qeht.tests.test_erosion
"""
import math
import os
import sys
import tempfile

import numpy as np

from ..core.flow.direction import d8_direction
from ..core.flow.accumulation import flow_accumulation
from ..core.erosion.terrain import (erosion_indices, mccool_m, mccool_s, m_by_slope_class,
                                    channel_mask, BETA_MIN_TAN)
from ..core.erosion.rusle import (Factor, rusle, c_from_worldcover, WORLDCOVER_C, WORLDCOVER_P)
from ..core.erosion.classes import (get_scheme, combine_scores, class_extents, write_class_raster,
                                    ClassScheme, SCHEMES)
from ..core.erosion.summary import (ErosionInputs, catchment_erosion, crossing_erosion,
                                    composite_score, impact_level, sdr_area)
from ..core.erosion.corridor import sample_corridor
from ..core.network.alignment import Alignment

FAILURES = []


def check(name, condition, detail=""):
    status = "PASS" if condition else "FAIL"
    if not condition:
        FAILURES.append(name)
    print(f"  [{status}] {name}" + (f"   {detail}" if detail else ""))


def plane(n=60, d=10.0, s=0.12):
    """Plane dipping south (row index increases downhill) with slope s."""
    z = np.tile(((n - np.arange(n)) * d * s)[:, None], (1, n)).astype(float)
    v = np.ones_like(z, dtype=bool)
    dirn, _ = d8_direction(z, v, d, d)
    acc, _ = flow_accumulation(dirn, v)
    return z, v, dirn, acc


def test_plane():
    print("\n1. Inclined plane (analytic)")
    n, d, s = 60, 10.0, 0.12
    z, v, _, acc = plane(n, d, s)
    r = erosion_indices(z, v, acc, d, d)
    i, j = 30, 30
    b = math.atan(s)
    a_s = (i + 1) * d
    check("accumulation = cells upslope in the column", acc[i, j] == i)
    check("tan(beta) = plane slope (Horn)", abs(r["tan_beta"][i, j] - s) < 1e-12)
    check("A_s = (acc + 1) x cell size", abs(r["spec_catch_area"][i, j] - a_s) < 1e-9)
    check("SPI = A_s tan(beta)", abs(r["spi"][i, j] - a_s * s) < 1e-9)
    check("ln SPI", abs(r["ln_spi"][i, j] - math.log(a_s * s)) < 1e-12)
    check("TWI = ln(A_s / tan(beta))", abs(r["twi"][i, j] - math.log(a_s / s)) < 1e-12)
    ls = (a_s / 22.13) ** 0.4 * (math.sin(b) / 0.0896) ** 1.3
    check("LS Moore & Burch", abs(r["ls"][i, j] - ls) < 1e-9, f"{r['ls'][i, j]:.6f}")
    rd = erosion_indices(z, v, acc, d, d, ls_method="desmet_govers")
    m = float(mccool_m(math.sin(b)))
    S = float(mccool_s(math.sin(b), s))
    ok = True
    for k in range(1, 50):                            # interior rows: A_in = k cells
        Lk = (((k + 1) * d * d) ** (m + 1) - (k * d * d) ** (m + 1)) / (d ** (m + 2) * 22.13 ** m)
        ok &= abs(rd["ls"][k, j] - Lk * S) < 1e-9 * max(1.0, Lk * S)
    check("LS Desmet & Govers per cell (McCool m and S)", ok, f"m={m:.4f} S={S:.4f}")
    telesc = sum((((k + 1) * d * d) ** (m + 1) - (k * d * d) ** (m + 1)) / (d ** (m + 2) * 22.13 ** m)
                 for k in range(0, 40))
    check("D&G L summed down the slope = (lambda/22.13)^m x n cells",
          abs(telesc / 40 - (40 * d / 22.13) ** m) < 1e-9)
    # m by slope class with the same formula
    rm = erosion_indices(z, v, acc, d, d, m_exponent="slope_classes")
    mc = 0.5                                           # 12 % is in the 10-20 % class
    ls2 = (a_s / 22.13) ** mc * (math.sin(b) / 0.0896) ** 1.3
    check("LS with m by slope class (12 % -> m 0.5)", abs(rm["ls"][i, j] - ls2) < 1e-9)
    # flat plane: beta floored inside the indices only
    zf = np.full((20, 20), 100.0); vf = np.ones_like(zf, dtype=bool)
    rf = erosion_indices(zf, vf, np.zeros_like(zf), 10, 10)
    check("flat: slope 0 reported, indices use tan(beta) floor 0.001",
          rf["tan_beta"][5, 5] == 0 and abs(rf["spi"][5, 5] - 10 * BETA_MIN_TAN) < 1e-12
          and np.isfinite(rf["twi"][5, 5]))


def test_valley():
    print("\n2. V-shaped valley: indices concentrate in the channel")
    n, d = 61, 10.0
    y, x = np.mgrid[0:n, 0:n]
    z = 200 - 0.05 * y * d + 0.2 * np.abs(x - 30) * d
    v = np.ones_like(z, dtype=bool)
    dirn, _ = d8_direction(z, v, d, d)
    acc, _ = flow_accumulation(dirn, v)
    r = erosion_indices(z, v, acc, d, d)
    ch = channel_mask(acc, v, 50)
    check("channel mask = the thalweg column (lower half)", ch[40:, 30].all() and ch.sum() < 3 * n)
    check("ln SPI and TWI peak on the thalweg",
          r["ln_spi"][50, 30] > np.nanmax(r["ln_spi"][50, :29])
          and r["twi"][50, 30] == np.nanmax(r["twi"][50, :]))


def test_a14_anchors():
    print("\n3. A14 corridor formula anchors (Akello & Omosa 2025, AJERI)")
    # Table 14: RUSLE, SY, SPI scores -> composite and impact level
    rows = [(1, 5, 3, 3, 3.6, "High"), (2, 5, 5, 3, 4.2, "Severe"), (4, 3, 5, 5, 4.4, "Severe"),
            (10, 4, 3, 3, 3.3, "High"), (11, 5, 3, 2, 3.2, "High"), (12, 4, 3, 3, 3.3, "High"),
            (18, 4, 3, 3, 3.3, "High"), (21, 4, 5, 3, 3.9, "High"), (23, 5, 3, 3, 3.6, "High"),
            (31, 4, 5, 5, 4.7, "Severe"), (37, 2, 5, 5, 4.1, "Severe"), (61, 1, 5, 5, 3.8, "High"),
            (62, 1, 5, 5, 3.8, "High"), (74, 1, 5, 5, 3.8, "High")]
    bad = [cid for cid, ru, sy, sp, comp, lvl in rows
           if abs(composite_score(sp, ru, sy) - comp) > 1e-9 or impact_level(comp)[0] != lvl]
    check("Table 14: 14 crossings, composite 0.4 SPI + 0.3 RUSLE + 0.3 SY and impact level",
          not bad, f"mismatches {bad}")
    check("impact thresholds 2.0 / 3.0 / 4.0",
          [impact_level(x)[0] for x in (1.0, 2.0, 2.1, 3.0, 3.1, 4.0, 4.1, 5.0)] ==
          ["Low", "Low", "Moderate", "Moderate", "High", "High", "Severe", "Severe"])
    spi = get_scheme("spi_a14")
    codes = spi.classify(np.array([-6.0, -0.1, 0.0, 4.9, 5.0, 9.9, 10.0, 15.0]), np.ones(8, bool))
    check("Table 4 SPI classes (ln SPI 0 / 5 / 10) and scores 1, 2, 3, 5",
          [spi.names[c - 1] for c in codes] == ["Low", "Low", "Moderate", "Moderate", "High",
                                                  "High", "Severe", "Severe"]
          and list(spi.score_of(codes)) == [1, 1, 2, 2, 3, 3, 5, 5])
    ru = get_scheme("rusle_a14")
    codes = ru.classify(np.array([0.0, 4.9, 5.0, 9.9, 10, 19.9, 20, 39.9, 40, 150]), np.ones(10, bool))
    check("Table 8 RUSLE classes 5 / 10 / 20 / 40 t/ha/yr",
          [ru.names[c - 1] for c in codes] == ["Slight", "Slight", "Moderate", "Moderate", "High",
                                                 "High", "Very high", "Very high", "Severe", "Severe"])
    vol = get_scheme("volume_a14")
    codes = vol.classify(np.array([1.0, 999, 1000, 4999, 5000, 14999, 15000, 815092]),
                         np.ones(8, bool))
    check("Table 9 sediment volume classes and scores 1, 2, 3, 5",
          list(vol.score_of(codes)) == [1, 1, 2, 2, 3, 3, 5, 5])
    tb = np.array([0.005, 0.02, 0.04, 0.07, 0.15, 0.25, 0.45])
    check("Table 5 m by slope class", list(m_by_slope_class(tb)) == [0.2, 0.3, 0.4, 0.45, 0.5,
                                                                        0.55, 0.6])
    table6 = {10: (0.001, 0.05), 20: (0.01, 0.05), 30: (0.01, 0.15), 40: (0.1, 0.4),
              50: (0.05, 0.2), 60: (0.4, 0.6)}
    check("Table 6 C lookup = mean of each ESA class range",
          all(abs(WORLDCOVER_C[k][0] - (lo + hi) / 2) < 1e-12 for k, (lo, hi) in table6.items()))
    check("Table 7 P lookup", [WORLDCOVER_P[k] for k in (10, 20, 30, 40, 50, 60)] ==
          [0.6, 0.8, 0.9, 0.7, 1.0, 1.0])
    check("SDR = 0.565 A^-0.125: within the reported 0.2-0.5 for 3-4,000 km2",
          abs(sdr_area(1.0) - 0.565) < 1e-12 and 0.2 < sdr_area(4000) < sdr_area(3) < 0.5)


def test_rusle_and_classes():
    print("\n4. RUSLE factors, fallback, classes and class rasters")
    rng = np.random.default_rng(4)
    ls = rng.gamma(1.5, 1.0, (40, 50)); v = np.ones_like(ls, dtype=bool)
    wc = rng.choice([10, 20, 30, 40, 60, 80], size=ls.shape)
    c, info = c_from_worldcover(wc, v)
    check("C from WorldCover: every class mapped, water = 0", np.isfinite(c).all()
          and (c[wc == 80] == 0).all() and not info["unknown_classes"])
    R = Factor("R", value=1800.0, source="GloREDa 2023")
    K = Factor("K", grid=np.full(ls.shape, 0.03), source="SOTWIS Williams")
    C = Factor("C", grid=c, source="ESA WorldCover 2021", proxy=True)
    a, inf = rusle(ls, v, R, K, C)
    check("A = R K LS C P (P defaults to 1.0)", np.allclose(a, 1800 * 0.03 * ls * c)
          and inf["mode"] == "RUSLE" and inf["proxy_factors"] == ["C"])
    a2, inf2 = rusle(ls, v, R, Factor("K"), C)
    check("missing K -> no soil loss, LS-only fallback stated", a2 is None
          and inf2["mode"].startswith("LS-only") and "K" in inf2["note"])
    sch = get_scheme("spi_percentile", ls, v)
    codes = sch.classify(ls, v)
    shares = [100.0 * (codes == k).mean() for k in range(1, 6)]
    check("percentile classes ~ 50 / 25 / 15 / 7 / 3 %",
          np.allclose(shares, [50, 25, 15, 7, 3], atol=0.6), str(np.round(shares, 1)))
    check("percentile breaks reproducible", get_scheme("spi_percentile", ls, v).breaks == sch.breaks)
    flat = np.zeros((10, 10)); flat[0, 0] = 1
    fb = get_scheme("ls_percentile", flat, np.ones_like(flat, bool)).breaks
    check("plateau: breaks kept strictly increasing", all(y > x for x, y in zip(fb, fb[1:])))
    s1 = np.array([[1, 5, 3]], dtype=np.uint8); s2 = np.array([[4, 2, 0]], dtype=np.uint8)
    check("combined = higher score (default matrix); NoData stays 0",
          combine_scores(s1, s2).tolist() == [[4, 5, 0]])
    m = [[1] * 5 for _ in range(5)]
    check("custom matrix honoured", combine_scores(s1, s2, m).tolist() == [[1, 1, 0]])
    ext = class_extents(codes, sch, 900.0, channel=ls > 3)
    check("class extents sum to 100 % and channel + hillslope = area",
          abs(sum(r["pct"] for r in ext) - 100) < 1e-9
          and all(abs(r["channel_km2"] + r["hillslope_km2"] - r["area_km2"]) < 1e-12 for r in ext))
    try:
        ClassScheme((1, 1, 2), ["a", "b", "c", "d"], [1, 2, 3, 4], "")
        refused = False
    except ValueError:
        refused = True
    check("non-increasing breaks refused", refused)
    from osgeo import gdal
    tmp = tempfile.mkdtemp()
    p = os.path.join(tmp, "spi_class.tif")
    write_class_raster(p, codes, sch, (0, 30, 0, 0, 0, -30), "", "SPI class")
    ds = gdal.Open(p); b = ds.GetRasterBand(1)
    rat = b.GetDefaultRAT()
    check("class raster: uint8, NoData 0, colour table, attribute table with names, .qml",
          b.DataType == gdal.GDT_Byte and b.GetNoDataValue() == 0 and b.GetColorTable() is not None
          and rat is not None and [rat.GetValueAsString(i, 1) for i in range(5)] == sch.names
          and os.path.exists(p[:-4] + ".qml") and "paletteEntry" in open(p[:-4] + ".qml").read())
    ds = None


def _synthetic_catchments():
    n, d = 61, 10.0
    y, x = np.mgrid[0:n, 0:n]
    z = 200 - 0.05 * y * d + 0.2 * np.abs(x - 30) * d
    v = np.ones_like(z, dtype=bool)
    dirn, _ = d8_direction(z, v, d, d)
    acc, _ = flow_accumulation(dirn, v)
    return z, v, dirn, acc, d


def test_blocks():
    print("\n5. Catchment and crossing erosion blocks")
    from ..core.watershed.delineate import delineate_catchment
    z, v, dirn, acc, d = _synthetic_catchments()
    r = erosion_indices(z, v, acc, d, d)
    ch = channel_mask(acc, v, 50)
    wc = np.where(np.abs(np.mgrid[0:61, 0:61][1] - 30) < 5, 30, 20)
    c, _ = c_from_worldcover(wc, v)
    a, _ = rusle(r["ls"], v, Factor("R", value=2000), Factor("K", value=0.03),
                 Factor("C", grid=c, proxy=True))
    inp = ErosionInputs(r["ln_spi"], r["ls"], r["twi"], r["tan_beta"], ch,
                        get_scheme("spi_a14"), get_scheme("rusle_a14"), soil_loss=a,
                        k=np.full(z.shape, 0.03), c=c, p=None, bulk_kgm3=1400)
    rr, cc = 59, 30
    mask = delineate_catchment(dirn, v, [(rr, cc)]) > 0
    cb = catchment_erosion(mask, inp, d * d)
    area_km2 = mask.sum() * d * d / 1e6
    gross = cb["ero_a_mean_tha"] * area_km2 * 100
    check("gross soil loss = mean A x area (ha)", abs(cb["ero_gross_t_yr"] - gross) < 1e-9)
    check("SDR and sediment volume follow the formulae",
          abs(cb["ero_sdr"] - min(1, 0.565 * area_km2 ** -0.125)) < 1e-12
          and abs(cb["ero_sy_m3_yr"] - gross * cb["ero_sdr"] * 1000 / 1400) < 1e-9)
    check("combined-score shares sum to 100 %",
          abs(sum(cb[f"ero_pct_s{s}"] for s in range(1, 6)) - 100) < 1e-9)
    check("K and C means are area-weighted", abs(cb["ero_k_mean"] - 0.03) < 1e-12
          and abs(cb["ero_c_mean"] - c[mask].mean()) < 1e-12)
    xb = crossing_erosion(rr, cc, mask, dirn, inp, cb)
    comp = 0.4 * xb["ero_spi_score"] + 0.3 * cb["ero_rusle_score"] + 0.3 * cb["ero_sy_score"]
    check("crossing composite from the three scores", abs(xb["ero_composite"] - comp) < 1e-12
          and xb["ero_impact"] == impact_level(comp)[0], f"{comp:.2f} {xb['ero_impact']}")
    check("approach channel ln SPI: p50 <= p90 <= outlet 3x3 max",
          xb["ero_lnspi_app_p50"] <= xb["ero_lnspi_app_p90"] <= xb["ero_lnspi_max3x3"] + 1e-12)
    inp2 = ErosionInputs(r["ln_spi"], r["ls"], r["twi"], r["tan_beta"], ch,
                         get_scheme("spi_a14"), get_scheme("ls_percentile", r["ls"], v))
    cb2 = catchment_erosion(mask, inp2, d * d)
    xb2 = crossing_erosion(rr, cc, mask, dirn, inp2, cb2)
    check("LS-only mode: no soil loss, no composite, stated", cb2["ero_mode"] == "LS-only"
          and cb2["ero_a_mean_tha"] is None and xb2["ero_composite"] is None
          and "LS-only" in xb2["ero_mitigation"])


def test_corridor():
    print("\n6. Corridor sampler")
    z, v, dirn, acc, d = _synthetic_catchments()
    r = erosion_indices(z, v, acc, d, d)
    score = np.where(np.abs(np.mgrid[0:61, 0:61][1] - 30) <= 2, 5, 1).astype(np.uint8)
    road = Alignment([[(5.0, -300.0), (605.0, -300.0)]], start_chainage=1000.0)
    st, re = sample_corridor(road, {"ln_spi": r["ln_spi"], "ls": r["ls"]}, score,
                             (0, d, 0, 0, 0, -d), step=10, half_width=30, offset_step=10)
    chs = [s["chainage"] for s in st]
    check("stations every 10 m from the start chainage, monotonic",
          chs[0] == 1000.0 and all(b > a for a, b in zip(chs, chs[1:])) and chs[-1] == 1600.0)
    check("reach lengths add up to the alignment length",
          abs(sum(x["length_m"] for x in re) - 600.0) < 1e-9)
    worst = [x for x in re if x["worst_score"] == 5]
    check("the valley crossing is one worst-score reach at ~ch 1300",
          len(worst) == 1 and worst[0]["ch_start"] < 1300 < worst[0]["ch_end"],
          f"{[(round(x['ch_start']), round(x['ch_end']), x['worst_score']) for x in re]}")
    check("left/right statistics present", "ln_spi_L_max" in st[0] and "ls_R_mean" in st[0])


def test_factor_checks():
    print("\n7. RUSLE factor range checks (v0.15.1, Site C run)")
    from ..core.erosion.rusle import check_factor, FactorError
    def refused(name, g=None, v=None, word=""):
        try:
            check_factor(name, g, v)
            return False
        except FactorError as e:
            return word in str(e)
    wc = np.where(np.arange(400)[None, :] % 3 == 0, 40, 30) * np.ones((50, 1))
    wc_bil = wc.astype(float); wc_bil[::7, 1::5] = 35.0      # bilinear mixing at edges
    check("WorldCover classes as C refused with a hint (also after bilinear resampling)",
          refused("C", wc, word="class codes") and refused("C", wc_bil, word="class codes"))
    check("C or P above 1 refused", refused("C", np.array([[0.2, 3.0]])) and refused("P", v=1.5))
    check("K in US units refused with the 0.1317 factor", refused("K", np.array([[0.32]]), word="0.1317"))
    check("negative values refused", refused("R", v=-5.0))
    w = check_factor("R", value=100.0)
    check("R = 100 accepted with a warning (too low / US units?)", len(w) == 1 and "17.02" in w[0])
    check("plausible factors pass silently",
          check_factor("R", value=3500.0) == [] and check_factor("K", np.array([[0.02, 0.04]])) == []
          and check_factor("C", np.array([[0.01, 0.3]])) == [] and check_factor("P", value=1.0) == [])


def main():
    test_plane()
    test_valley()
    test_a14_anchors()
    test_rusle_and_classes()
    test_blocks()
    test_corridor()
    test_factor_checks()
    print("\n" + "=" * 62)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S):")
        for f in FAILURES:
            print("   -", f)
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
