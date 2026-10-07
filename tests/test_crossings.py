# -*- coding: utf-8 -*-
"""v0.10 validation: crossing candidates (WP-A), burn crossings (WP-D),
renumber and relink (D3). Bare Python + NumPy, no QGIS, no GDAL.

    python -m qeht.tests.test_crossings
"""

import math
import sys

import numpy as np

from ..core.conditioning.fill import fill_depressions, depression_depth
from ..core.conditioning.burn import burn_crossings
from ..core.flow.direction import d8_direction
from ..core.flow.accumulation import flow_accumulation
from ..core.grid import DROW, DCOL
from ..core.watershed.delineate import extract_streams
from ..core.network.alignment import Alignment
from ..core.network.crossings import find_crossing_candidates, select_crossings
from ..core.interop.heas_exchange import build_exchange_records
from ..core.linking.relink import renumber_log

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def valley_dem(paths, rows, cols, cs, slope=0.02, side=0.5, base=100.0):
    """V-valleys along polylines (lists of (x, y)): the floor falls at `slope`
    toward each path's end and the sides rise at `side`; lowest surface wins."""
    r, c = np.mgrid[0:rows, 0:cols].astype(float)
    x = (c + 0.5) * cs
    y = rows * cs - (r + 0.5) * cs
    z = np.full((rows, cols), np.inf)
    for p in paths:
        al = Alignment([p])
        ch, off, _ = al.locate(x.ravel(), y.ravel())
        z = np.minimum(z, (base + slope * (al.length - ch) + side * np.abs(off)).reshape(rows, cols))
    return z + 1e-4 * c


def route(dem, cs, threshold=40):
    v = np.ones(dem.shape, bool)
    f, _, _ = fill_depressions(dem, v, cell_width=cs, cell_height=cs)
    d, _ = d8_direction(f, v, cs, cs)
    a, _ = flow_accumulation(d, v)
    return v, f, d, a, extract_streams(a, v, threshold_cells=threshold)


def gt_for(rows, cs):
    return (0.0, cs, 0.0, rows * cs, 0.0, -cs)


# ------------------------------------------------------------------
def test_alignment():
    print("\n1. Linear referencing")
    al = Alignment([[(0, 0), (100, 0)]], start_chainage=1000)
    ch, off, _ = al.locate([30, 30, 150], [5, -5, 0])
    check("chainage = start + distance", ch[0] == 1030 and ch[1] == 1030)
    check("offset signed: + left, - right", off[0] == 5 and off[1] == -5)
    check("beyond the end: clamped chainage, distance to end", ch[2] == 1100 and abs(off[2]) == 50)
    al2 = Alignment([[(0, 0), (100, 0)], [(200, 0), (300, 0)]])
    check("multi-part: gap between parts not counted", al2.locate([250], [1])[0][0] == 150
          and al2.length == 200)
    al3 = Alignment([[(0, 0), (100, 0)]], reverse=True)
    c3, o3, _ = al3.locate([30], [5])
    check("reverse direction", c3[0] == 70 and o3[0] == -5)
    check("point_at", al.point_at(1025) == (25.0, 0.0))
    i = al.intersect([50, 10], [-10, -10], [50, 10], [10, -5])
    check("segment intersection found once", i[0].size == 1 and i[4][0] == 50 and i[5][0] == 0)


def test_perpendicular_and_angle():
    print("\n2. Single crossings: chainage, side, angle")
    rows, cols, cs = 80, 60, 10.0
    dem = valley_dem([[(305, 800), (305, 0)]], rows, cols, cs)
    v, f, d, a, st = route(dem, cs)
    gt = gt_for(rows, cs)
    y = 400.0 + 3.0
    cand, par, s = find_crossing_candidates(d, v, a, st, gt, Alignment([[(0, y), (600, y)]],
                                            start_chainage=1000), corridor_halfwidth_m=30)
    check("one candidate", len(cand) == 1, str(s))
    at = cand[0][1]
    check("chainage exact (start 1000 + 305 m)", math.isclose(at["chainage_m"], 1305.0), f"{at['chainage_m']}")
    check("square crossing = 90 deg", math.isclose(at["crossing_angle_deg"], 90.0))
    check("flow comes from the left (north of an eastbound road)", at["side_in"] == "L" and at["side_out"] == "R")
    check("recommended, cluster 1, status candidate",
          at["recommended"] == 1 and at["cluster_id"] == 1 and at["status"] == "candidate")
    check("contributing area = (acc+1) x cell at the upstream cell",
          math.isclose(at["acc_km2"], (a[at["outlet_row"], at["outlet_col"]] + 1) * 100 / 1e6))
    check("no parallel reach on a square crossing", not par and at["parallel_reach"] == 0)
    cand, _, _ = find_crossing_candidates(d, v, a, st, gt, Alignment([[(0, 100), (600, 700)]]),
                                          corridor_halfwidth_m=30)
    angles = [c[1]["crossing_angle_deg"] for c in cand if c[1]["acc_km2"] > 0.05]
    check("45 deg road over a N-S channel -> 45 deg", len(angles) == 1 and math.isclose(angles[0], 45.0),
          str(angles))


def _weaving_case():
    rows, cols, cs, yr = 100, 90, 10.0, 500.0
    wob = [(x, yr + 12 * math.sin(2 * math.pi * (x - 100) / 120)) for x in np.arange(100, 701, 10)]
    dem = valley_dem([[(100, 1000)] + wob + [(700, 0)]], rows, cols, cs)
    v, f, d, a, st = route(dem, cs)
    return d, v, a, st, gt_for(rows, cs), Alignment([[(0, yr + 3), (900, yr + 3)]])


def test_parallel_flow():
    print("\n3. Stream running alongside the road, then crossing (N raw -> 1 cluster)")
    d, v, a, st, gt, al = _weaving_case()
    cand, par, s = find_crossing_candidates(d, v, a, st, gt, al, min_area_km2=0.02,
                                            corridor_halfwidth_m=30, min_parallel_m=200)
    n = len(cand)
    check("several raw intersections", n >= 5, f"{n} candidates")
    check("all in one cluster", len({c[1]["cluster_id"] for c in cand}) == 1)
    rec = [c[1] for c in cand if c[1]["recommended"]]
    exitc = max(cand, key=lambda c: c[1]["acc_km2"])[1]
    check("exactly one recommended = the exit (largest area)", len(rec) == 1 and rec[0] is exitc,
          f"{rec[0]['cand_id']} at ch {rec[0]['chainage_m']:.0f} m, {rec[0]['acc_km2']:.3f} km2")
    check("all flagged parallel_reach", all(c[1]["parallel_reach"] == 1 for c in cand))
    check("a parallel reach >= 500 m exported with its cluster",
          any(p[1]["length_m"] >= 500 and p[1]["cluster_id"] == 1 for p in par),
          ", ".join(f"{p[1]['length_m']:.0f} m" for p in par))
    chs = [c[1]["chainage_m"] for c in cand]
    check("chainage monotonic", chs == sorted(chs))
    check("dist_prev sums to first->last",
          math.isclose(sum(c[1]["dist_prev_m"] for c in cand[1:]), chs[-1] - chs[0]))
    check("ids along chainage", [c[1]["cand_id"] for c in cand] == [f"C{k + 1:03d}" for k in range(n)])


def test_close_tributaries():
    print("\n4. Two tributaries 100 m apart: merge distance decides")
    rows, cols, cs = 80, 70, 10.0
    dem = valley_dem([[(305, 800), (305, 0)], [(405, 800), (405, 0)]], rows, cols, cs)
    v, f, d, a, st = route(dem, cs)
    gt = gt_for(rows, cs)
    al = Alignment([[(0, 403), (700, 403)]])
    c1, _, _ = find_crossing_candidates(d, v, a, st, gt, al, min_area_km2=0.05, merge_distance_m=50)
    c2, _, _ = find_crossing_candidates(d, v, a, st, gt, al, min_area_km2=0.05, merge_distance_m=150)
    check("merge 50 m -> 2 clusters, both recommended",
          len(c1) == 2 and len({c[1]["cluster_id"] for c in c1}) == 2
          and sum(c[1]["recommended"] for c in c1) == 2)
    check("merge 150 m -> 1 cluster, one recommended",
          len({c[1]["cluster_id"] for c in c2}) == 1 and sum(c[1]["recommended"] for c in c2) == 1)


def test_winding_double_crossing():
    print("\n4b. A winding centreline crossing one stream link twice within metres (v0.15.1)")
    rows, cols, cs = 80, 70, 10.0
    dem = valley_dem([[(305, 800), (305, 0)]], rows, cols, cs)
    v, f, d, a, st = route(dem, cs)
    gt = gt_for(rows, cs)
    al = Alignment([[(0, 401), (312, 401), (312, 403), (298, 403), (298, 404.5), (700, 404.5)]])
    c, _, _ = find_crossing_candidates(d, v, a, st, gt, al, min_area_km2=0.05)
    check("three intersections, one cluster, one recommended (auto merge < 2 cell diagonals)",
          len(c) == 3 and len({x[1]["cluster_id"] for x in c}) == 1
          and sum(x[1]["recommended"] for x in c) == 1,
          f"{len(c)} candidates, clusters {[x[1]['cluster_id'] for x in c]}")
    two, _, _ = find_crossing_candidates(d, v, a, st, gt,
                                         Alignment([[(0, 401), (700, 401)], [(0, 395), (700, 395)]]),
                                         min_area_km2=0.05)
    check("separate parts far apart in chainage stay separate", len({x[1]["cluster_id"] for x in two}) == 2)


def test_selection_and_exchange():
    print("\n5. Candidates -> crossings -> exchange records")
    feats = [{"status": "candidate", "recommended": 1}, {"status": "accepted", "recommended": 0},
             {"status": "rejected", "recommended": 1}]
    idx, rule = select_crossings(feats)
    check("any 'accepted' -> exactly the accepted", idx == [1], rule)
    idx, _ = select_crossings([{"status": "candidate", "recommended": 1},
                               {"status": "rejected", "recommended": 1},
                               {"status": None, "recommended": 0}])
    check("else recommended, never rejected", idx == [0])

    rows, cols, cs = 80, 70, 10.0
    dem = valley_dem([[(305, 800), (305, 0)], [(405, 800), (405, 0)]], rows, cols, cs)
    v, f, d, a, st = route(dem, cs)
    gt = gt_for(rows, cs)
    cand, _, _ = find_crossing_candidates(d, v, a, st, gt, Alignment([[(700, 403), (0, 403)]],
                                          start_chainage=2000), min_area_km2=0.05)
    outlets = [dict(x=g[0], y=g[1], fid=k, source_id=None, outlet_x=at["outlet_x"],
                    outlet_y=at["outlet_y"], chainage=at["chainage_m"])
               for k, (g, at) in enumerate(cand)]
    cr, ca, fp, iss, info = build_exchange_records(d, v, a, dem, gt, outlets, id_order="chainage",
                                                   snap_radius_cells=5, stream_mask=st)
    by = {at["outlet_uid"]: at for _, at in cr}
    check("chainage carried to crossings", all(at["chainage_m"] is not None for at in by.values()))
    check("X001 = smallest chainage (numbered along the alignment)",
          by["X001"]["chainage_m"] < by["X002"]["chainage_m"], info.get("id_order"))
    check("fixed outlets are not snapped: outlet = candidate cell, <= 1 cell from the crossing",
          all(at["outlet_x"] == cand[at["outlet_id"]][1]["outlet_x"]
              and at["snap_dist_m"] <= cs * math.sqrt(2) for at in by.values()),
          ", ".join(f"{at['snap_dist_m']:.1f} m" for at in by.values()))
    ok = all(math.isclose(c[1]["area_km2"], cand[k][1]["acc_km2"])
             for k, c in enumerate(sorted(ca, key=lambda z: z[1]["outlet_id"])))
    check("catchment area == candidate acc_km2", ok)


def test_burn():
    print("\n6. Burn crossings through an embankment (WP-D)")
    rows, cols, cs = 90, 60, 10.0
    dem = valley_dem([[(305, 900), (305, 0)]], rows, cols, cs)
    road_y = 450.0
    r_emb = [int((rows * cs - road_y) / cs) + k for k in (-2, -1, 0, 1, 2)]
    col_c = 30
    crest = dem[r_emb[2], col_c] + 5.0
    emb = dem.copy()
    for r in r_emb:                     # embankment: flat crest 5 m above the valley floor
        emb[r, :] = np.maximum(emb[r, :], crest)
    v = np.ones(emb.shape, bool)
    f0, _, _ = fill_depressions(emb, v, cell_width=cs, cell_height=cs)
    pond0 = float(depression_depth(emb, f0, v).sum())
    burned, log = burn_crossings(emb, v, gt_for(rows, cs),
                                 [{"x": 305.0, "y": road_y, "tx": 1.0, "ty": 0.0, "id": "X001"}],
                                 half_length_m=40.0, search_radius_cells=2)
    f1, _, _ = fill_depressions(burned, v, cell_width=cs, cell_height=cs)
    pond1 = float(depression_depth(burned, f1, v).sum())
    check("embankment ponds the valley before burning", pond0 > 10, f"sum of fill depth {pond0:.1f} m")
    check("pond removed after burning", pond1 < 1e-6 * max(pond0, 1), f"{pond1:.4f} m")
    check("breach only lowers ground", np.all(burned <= emb + 1e-12))
    changed = np.argwhere(burned < emb)
    check("breach confined to the crossing (<= 1 cell off the culvert line)",
          changed.size > 0 and np.all(np.abs(changed[:, 1] - col_c) <= 1), f"{len(changed)} cells")
    e = log[0]
    check("log: cells, max cut ~ 5 m, grade falls downstream",
          e["cells"] == len(changed) and 4.0 < e["max_cut_m"] <= 5.5 and e["z_down"] < e["z_up"],
          f"cells {e['cells']}, max cut {e['max_cut_m']:.2f} m, {e['z_up']:.2f} -> {e['z_down']:.2f}")
    _, log2 = burn_crossings(emb, v, gt_for(rows, cs),
                             [{"x": 305.0, "y": road_y, "azimuth_deg": 90.0, "id": "X002"}],
                             half_length_m=10.0, search_radius_cells=6)
    check("breach ends stay on opposite sides (short breach, wide search)",
          log2[0]["cells"] > 0 and not log2[0]["note"], f"{log2[0]['cells']} cells {log2[0]['note']}")
    d1, _ = d8_direction(f1, v, cs, cs)
    r, c = 5, col_c                      # follow the flow from far upstream
    crossed_at = None
    for _ in range(400):
        k = d1[r, c]
        if k < 0:
            break
        r, c = r + int(DROW[k]), c + int(DCOL[k])
        if r == r_emb[2]:
            crossed_at = c
            break
    check("flow passes through the breach at the culvert", crossed_at is not None
          and abs(crossed_at - col_c) <= 1, f"crossed at col {crossed_at}")


def test_relink():
    print("\n7. Renumber and relink (D3)")
    crossings_new = [((10.0, 0.0), {"outlet_id": 1, "outlet_uid": "X001", "chainage_m": 10.0}),
                     ((50.0, 0.0), {"outlet_id": 7, "outlet_uid": "X002", "chainage_m": 50.0}),
                     ((90.0, 0.0), {"outlet_id": 3, "outlet_uid": "X003", "chainage_m": 90.0})]
    old_by_fid = {1: {"uid": "X001", "x": 10.0, "y": 0.0},
                  7: {"uid": None, "x": None, "y": None},          # added by the user
                  3: {"uid": "X004", "x": 90.0, "y": 30.0}}       # moved and now third
    rows = renumber_log(crossings_new, old_by_fid, ["X001", "X002", "X003", "X004"])
    ch = {r["old_uid"]: r for r in rows if r["old_uid"]}
    check("unchanged id kept", ch["X001"]["change"] == "unchanged")
    check("added point -> 'new'", any(r["change"] == "new" and r["new_uid"] == "X002" for r in rows))
    check("renumbered with distance moved", ch["X004"]["change"] == "renumbered"
          and ch["X004"]["new_uid"] == "X003" and math.isclose(ch["X004"]["moved_m"], 30.0))
    check("deleted ids listed", sorted(r["old_uid"] for r in rows if r["change"] == "deleted")
          == ["X002", "X003"])
    check("new ids gapless", sorted(r["new_uid"] for r in rows if r["new_uid"]) == ["X001", "X002", "X003"])


def main(argv=None):
    print("=" * 62)
    print("QEHT CROSSINGS - v0.10 (WP-A candidates, WP-D burn, D3 relink)")
    print("=" * 62)
    test_alignment()
    test_perpendicular_and_angle()
    test_parallel_flow()
    test_close_tributaries()
    test_winding_double_crossing()
    test_selection_and_exchange()
    test_burn()
    test_relink()
    print("\n" + "=" * 62)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S):")
        for f in FAILURES:
            print("   -", f)
        return 1
    print("ALL CHECKS PASSED")
    print("=" * 62)
    return 0


if __name__ == "__main__":
    sys.exit(main())
