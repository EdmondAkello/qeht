# -*- coding: utf-8 -*-
"""v0.16 validation: drainage coverage check, sag points, flat stretches and
proposed crossings (A2, A3). Bare Python + NumPy, no QGIS, no GDAL.

    python -m qeht.tests.test_coverage
"""

import math
import sys

import numpy as np

from ..core.network.alignment import Alignment
from ..core.network.profile import profile_with_crossings
from ..core.network.coverage import (cluster_cover, missing_crossings, sag_points, flat_stretches,
                                     walled_accumulation, sag_areas, run_coverage, smooth,
                                     alignment_wall)
from ..core.interop.heas_exchange import build_exchange_records
from .test_crossings import valley_dem, route, gt_for

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def _profile(z, step=10.0):
    return [{"chainage_m": k * step, "x": k * step, "y": 0.0, "z_dem_m": float(v),
             "stream": 0, "acc_km2": None, "pond_depth_m": 0.0} for k, v in enumerate(z)]


def test_missing():
    print("\n1. Missing crossings: two streams, one existing crossing")
    rows, cols, cs = 80, 90, 10.0
    gt = gt_for(rows, cs)
    dem = valley_dem([[(255, 800), (255, 0)], [(655, 800), (655, 0)]], rows, cols, cs)
    v, f, d, a, st = route(dem, cs, threshold=40)
    al = Alignment([[(5, 395), (895, 395)]])
    prof, cr = profile_with_crossings(al, dem, gt, direction=d, valid=v, accumulation=a,
                                      filled=f, stream_threshold_cells=40, step=10.0)
    pts = [{"x": 255.0, "y": 395.0, "fid": 1, "source_id": None, "chainage": 250.0}]
    ex, ca, fp, _, _ = build_exchange_records(d, v, a, dem, gt, pts, snap_radius_cells=3,
                                              stream_mask=st, id_prefix="X",
                                              channel_threshold_cells=40)
    miss = missing_crossings(prof, [(c[1]["outlet_uid"], c[1]["chainage_m"]) for c in ex],
                             min_area_km2=0.0)
    unc = [m for m in miss if not m["covered"]]
    check("two stream crossings found, one covered by X001, one uncovered at ~650 m",
          len(miss) == 2 and len(unc) == 1 and abs(unc[0]["chainage_m"] - 650) <= 10,
          f"{[(round(m['chainage_m']), m['covered']) for m in miss]}")
    far = missing_crossings(prof, [("X001", 250.0)], min_area_km2=1e6)
    check("area threshold honoured (huge minimum -> none)", far == [])
    tight = missing_crossings(prof, [("X001", 290.0)], search_m=30.0)
    check("search distance honoured (crossing 40 m away no longer covers)",
          sum(1 for m in tight if not m["covered"]) == 2)
    dup = [dict(r) for r in prof]
    for r in dup:
        if abs(r["chainage_m"] - 670) < 1e-9:
            r.update(stream=1, acc_km2=0.01)
    merged = missing_crossings(dup, [], merge_m=30.0)
    split = missing_crossings(dup, [], merge_m=5.0)
    check("merge distance respected (20 m apart: merged at 30 m, split at 5 m)",
          len(merged) == 2 and len(split) == 3)
    res = run_coverage(d, v, a, dem, f, gt, al, prof, ex,
                       build_kwargs=dict(channel_threshold_cells=40), stream_threshold_cells=40,
                       sag_min_depth_m=1e9)
    pc = res["crossings"]
    check("one proposed crossing P001, delineated and characterised",
          len(pc) == 1 and pc[0][1]["outlet_uid"] == "P001" and pc[0][1]["status"] == "proposed"
          and pc[0][1]["proposed_reason"] == "uncovered_stream" and pc[0][1]["nearest_uid"] == "X001"
          and len(res["catchments"]) == 1 and res["catchments"][0][1]["area_km2"] > 0.05
          and res["flowpaths"][0][1]["lfp_length_m"] > 100,
          f"{pc[0][1]['outlet_uid'] if pc else '-'} A={res['catchments'][0][1]['area_km2']:.3f} km2"
          if pc else "")
    check("existing crossings default to status 'existing'",
          ex[0][1]["status"] == "existing" and ca[0][1]["status"] == "existing")
    check("coverage layer lists the missing crossing with its proposed ID",
          any(r["issue"] == "missing_crossing" and r["uid"] == "P001" for _, r in res["coverage"]))
    try:
        run_coverage(d, v, a, dem, f, gt, al, prof,
                     [((0, 0), dict(ex[0][1], outlet_uid="P001"))], stream_threshold_cells=40,
                     sag_min_depth_m=1e9)
        clash = False
    except Exception as e:
        clash = "already used" in str(e)
    check("a proposed ID never collides with an existing one", clash)


def test_sags():
    print("\n2. Sag points on a synthetic profile")
    ch = np.arange(0, 1001, 10.0)
    z = 100.0 + 0.0 * ch
    z -= 1.0 * np.exp(-((ch - 300) / 30.0) ** 2)      # 1.0 m sag at 300
    z -= 0.5 * np.exp(-((ch - 700) / 30.0) ** 2)      # 0.5 m sag at 700
    z -= 0.1 * np.exp(-((ch - 900) / 20.0) ** 2)      # 0.1 m dip: below the threshold
    prof = _profile(z)
    s = sag_points(prof, smooth_m=0.0, min_depth_m=0.3)
    check("two sags at 300 and 700 m with exact depths 1.0 and 0.5 m",
          len(s) == 2 and [x["chainage_m"] for x in s] == [300.0, 700.0]
          and abs(s[0]["sag_depth_m"] - 1.0) < 1e-6 and abs(s[1]["sag_depth_m"] - 0.5) < 1e-6,
          f"{[(x['chainage_m'], round(x['sag_depth_m'], 4)) for x in s]}")
    check("depth threshold honoured", len(sag_points(prof, 0.0, 0.6)) == 1
          and len(sag_points(prof, 0.0, 0.05)) == 3)
    slope = _profile(100.0 - 0.01 * ch)
    check("a steady grade has no sag", sag_points(slope, 30.0, 0.0) == [])
    asym = _profile(np.r_[np.linspace(105, 100, 51), np.linspace(100, 101.2, 50)])
    sa = sag_points(asym, 0.0, 0.3)
    check("depth = the LOWER side (prominence): 1.2 m, not 5 m",
          len(sa) == 1 and abs(sa[0]["sag_depth_m"] - 1.2) < 1e-6, f"{sa}")
    from ..core.network.coverage import ponding_sags
    p2 = [dict(r) for r in prof]
    p2[30]["stream"] = 1                                   # a stream crosses at the 300 m sag
    ps = ponding_sags(p2, [("X001", 690.0)], 50.0, 0.0, 0.3)
    check("sags at a stream crossing or an existing crossing are not ponding sags",
          ps == [], f"{[x['chainage_m'] for x in ps]}")
    sm = smooth(np.r_[0, 0, 3, 0, 0.0], np.arange(5) * 10.0, 30.0)
    check("smoothing = centred moving average", abs(sm[2] - 1.0) < 1e-12)


def test_sag_area():
    print("\n3. Local area draining to a sag with the alignment as a wall")
    rows, cols, cs = 60, 80, 10.0
    gt = gt_for(rows, cs)
    r, c = np.mgrid[0:rows, 0:cols].astype(float)
    x = (c + 0.5) * cs
    y = rows * cs - (r + 0.5) * cs
    # ground falls to the south (towards the road at y=300) and towards x=400 along it;
    # south of the road it falls further south (so flow would cross the road)
    z = 100.0 + 0.02 * np.abs(y - 300.0) * np.where(y > 300, 1, -1) + 0.01 * np.abs(x - 405.0)
    v = np.ones(z.shape, bool)
    al = Alignment([[(5, 305), (795, 305)]])
    acc_w, wall, wd = walled_accumulation(z, v, al, gt, cs, cs)
    prof = [{"chainage_m": k * 10.0, "x": 5 + k * 10.0, "y": 305.0,
             "z_dem_m": 100.0 + 0.01 * abs(5 + k * 10.0 - 405)} for k in range(80)]
    s = sag_points(prof, 0.0, 0.3)
    sag_areas(s, acc_w, wall, gt, al, direction=wd)
    north = int((y > 305.0 + 5).sum())         # all cells north of the wall drain to the sag
    check("one sag at x = 405, contributing side = left (north)",
          len(s) == 1 and abs(s[0]["x"] - 405) <= 10 and s[0]["side"] == "left", f"{s}")
    check("sag area = the whole north side (wall stops it crossing), within 5 %",
          s and abs(s[0]["sag_area_km2"] / (north * cs * cs / 1e6) - 1) < 0.05,
          f"{s[0]['sag_area_km2']:.4f} vs {north * cs * cs / 1e6:.4f} km2" if s else "")
    acc_g, wall_g, wdg = walled_accumulation(z, v, al, gt, cs, cs, gaps=[(655.0, 305.0)])
    s2 = sag_points(prof, 0.0, 0.3)
    sag_areas(s2, acc_g, wall_g, gt, al, direction=wdg)
    check("a culvert (gap) at 655 m lets part of the water through: the sag gets less",
          s2 and 0.5 * s[0]["sag_area_km2"] < s2[0]["sag_area_km2"] < 0.97 * s[0]["sag_area_km2"],
          f"{s2[0]['sag_area_km2']:.4f} vs {s[0]['sag_area_km2']:.4f} km2" if s2 else "")
    w = alignment_wall(Alignment([[(5, 5), (595, 595)]]), (60, 60), gt_for(60, 10.0))
    # a diagonal wall must have no diagonal gap: every 2x2 block on it has >= 3 wall cells or none crossing
    gaps = int(np.sum(w[:-1, :-1] & w[1:, 1:] & ~w[:-1, 1:] & ~w[1:, :-1]))
    check("diagonal wall is closed (no diagonal leak)", gaps == 0)


def test_flats():
    print("\n4. Flat stretches")
    rows, cols, cs = 40, 120, 10.0
    gt = gt_for(rows, cs)
    r, c = np.mgrid[0:rows, 0:cols].astype(float)
    x = (c + 0.5) * cs
    # flat (0.1 %) between x = 300 and 800, 3 % elsewhere along x; no cross-fall
    z = np.where(x < 300, 100 - 0.03 * x,
                 np.where(x <= 800, 91 - 0.001 * (x - 300), 90.5 - 0.03 * (x - 800)))
    al = Alignment([[(5, 195), (1195, 195)]])
    prof, _ = profile_with_crossings(al, z, gt, step=10.0)
    fl = flat_stretches(prof, z, gt, al, flat_slope_pct=0.5, crossfall_m=100.0, min_len_m=300.0,
                        existing=[("X001", 500.0)])
    check("one flat stretch ~300-800 m, length within one step of 500 m",
          len(fl) == 1 and abs(fl[0]["length_m"] - 500.0) <= 20.0
          and fl[0]["n_crossings_within"] == 1, f"{[(f['chainage_m'], f['chainage_to_m']) for f in fl]}")
    check("minimum length honoured", flat_stretches(prof, z, gt, al, 0.5, 100.0, 600.0) == [])
    tilted = z + 0.01 * ((rows * cs - (r + 0.5) * cs) - 200)      # 1 % cross-fall
    pt, _ = profile_with_crossings(al, tilted, gt, step=10.0)
    check("cross-fall above the limit -> not flat", flat_stretches(pt, tilted, gt, al, 0.5, 100.0, 300.0) == [])


def test_cluster_cover():
    print("\n[Clustered candidates cover their stream]")
    cands = [("c1", 100.0), ("c1", 160.0), ("c1", 240.0), ("c2", 900.0), (None, 1500.0)]
    used = [("c1", "X001")]
    cov = cluster_cover(cands, used)
    check("every member of a used cluster is covered by that crossing; other clusters are not",
          cov == [("X001", 100.0), ("X001", 160.0), ("X001", 240.0)], str(cov))
    # a stream beside the road crosses the centreline at 100, 160 and 240 m; the crossing used
    # is at 100 m, so with a 50 m search the crossings at 160 and 240 m count as missing
    # unless the cluster covers them
    prof = [{"chainage_m": c, "x": c, "y": 0.0, "acc_km2": 30.0, "stream": 1} for c in (100.0, 160.0, 240.0)]
    plain = missing_crossings(prof, [("X001", 100.0)], 1.0, 30.0, 50.0)
    clus = missing_crossings(prof, [("X001", 100.0)] + cov, 1.0, 30.0, 50.0)
    check("parallel stream: missing without the cluster, covered with it",
          sum(not m["covered"] for m in plain) == 2 and all(m["covered"] for m in clus))


def main(argv=None):
    test_missing()
    test_sags()
    test_sag_area()
    test_flats()
    test_cluster_cover()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
