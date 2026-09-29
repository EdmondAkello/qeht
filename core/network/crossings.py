# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Road x drainage crossing candidates (WP-A).

Method
------
1. Every D8 flow link i -> j of the stream network near the road (a straight
   segment between cell centres - exactly what QEHT's stream polylines are
   made of) is intersected with the alignment. Because the D8 network is
   continuous, a crossing cannot be missed through a digitising gap; this
   replaces both a vector stream x road intersection and a raster
   near-miss check with one exact operation.
2. Each intersection is a candidate with chainage, crossing point, crossing
   angle (90 = square), contributing area at the upstream cell i
   ((accumulation + 1) x cell area), Strahler order, stream reach id and
   the side the flow comes from (L/R of the alignment direction).
3. Parallel flow: stream cells within the corridor half-width of the road
   are traced downstream into runs. A run at least `min_parallel_m` long
   is a parallel reach (exported as a line - a side-drain hint). All
   candidates touching the same connected parallel system share a
   cluster; so do candidates closer than `merge_distance_m` in chainage.
4. Recommendation, not deletion (decision D3): in each cluster the most
   downstream candidate (largest contributing area) gets recommended = 1.
   Every candidate is written with status = 'candidate'; the engineer sets
   accepted / rejected.

No QGIS imports, no GDAL.
"""

import math

import numpy as np

from ..grid import DROW, DCOL, NO_RECEIVER, neighbour_distances

STATUS_VALUES = ("candidate", "accepted", "rejected")

CANDIDATE_FIELDS = [
    ("cand_id", "text"), ("chainage_m", "real"), ("status", "text"),
    ("recommended", "int"), ("cluster_id", "int"), ("parallel_reach", "int"),
    ("acc_km2", "real"), ("stream_order", "int"), ("reach_id", "int"),
    ("crossing_angle_deg", "real"), ("road_azimuth_deg", "real"),
    ("side_in", "text"), ("side_out", "text"),
    ("dist_prev_m", "real"), ("alignment_part", "int"),
    ("outlet_x", "real"), ("outlet_y", "real"), ("outlet_row", "int"), ("outlet_col", "int"),
]

PARALLEL_FIELDS = [
    ("pr_id", "int"), ("length_m", "real"), ("chainage_from_m", "real"),
    ("chainage_to_m", "real"), ("side", "text"), ("acc_max_km2", "real"),
    ("cluster_id", "int"),
]


class _UnionFind(object):
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[max(ra, rb)] = min(ra, rb)


def _dilate(mask, k):
    """Chebyshev (square) dilation by k cells, separable, no wrap-around."""
    out = mask.copy()
    for axis in (0, 1):
        acc = out.copy()
        for s in range(1, k + 1):
            if axis == 0:
                acc[s:, :] |= out[:-s, :]; acc[:-s, :] |= out[s:, :]
            else:
                acc[:, s:] |= out[:, :-s]; acc[:, :-s] |= out[:, s:]
        out = acc
    return out


def find_crossing_candidates(direction, valid, accumulation, stream_mask, geotransform,
                             alignment, min_area_km2=0.0, corridor_halfwidth_m=30.0,
                             min_parallel_m=100.0, merge_distance_m=0.0,
                             stream_order=None, reach_ids=None, prefix="C"):
    """Crossing candidates and parallel reaches along `alignment`.

    Parameters
    ----------
    direction, valid : internal D8 direction grid and validity mask
    accumulation : accumulation in cells (QEHT convention, excludes the cell)
    stream_mask : bool grid of the stream network
    geotransform : GDAL-style, north-up
    alignment : core.network.alignment.Alignment in the DEM CRS
    min_area_km2 : candidates draining less than this are dropped
    corridor_halfwidth_m : distance either side of the centreline treated
        as "alongside the road" for parallel-flow detection
    min_parallel_m : along-flow length inside the corridor that makes a
        stream run a parallel reach
    merge_distance_m : candidates closer than this in chainage share a
        cluster (0 = merge only through parallel reaches)
    stream_order, reach_ids : optional grids for the attributes

    Returns (candidates, parallel_reaches, summary) where candidates is a
    list of ((x, y), attrs) sorted by chainage and parallel_reaches a list
    of ([(x, y), ...], attrs).
    """
    gt = tuple(geotransform)
    cw, chh = abs(gt[1]), abs(gt[5])
    cell_area = cw * chh
    rows, cols = direction.shape
    hw = float(corridor_halfwidth_m)

    # -- window around the road ---------------------------------------
    xs = np.concatenate([alignment.x0, alignment.x1])
    ys = np.concatenate([alignment.y0, alignment.y1])
    pad = int(math.ceil(max(hw, 0.0) / min(cw, chh))) + 3
    c_lo = max(0, int(math.floor((xs.min() - gt[0]) / gt[1])) - pad)
    c_hi = min(cols, int(math.floor((xs.max() - gt[0]) / gt[1])) + pad + 1)
    r_lo = max(0, int(math.floor((ys.max() - gt[3]) / gt[5])) - pad)
    r_hi = min(rows, int(math.floor((ys.min() - gt[3]) / gt[5])) + pad + 1)
    if c_lo >= c_hi or r_lo >= r_hi:
        return [], [], {"n_candidates": 0, "n_clusters": 0, "n_parallel": 0,
                        "note": "alignment lies outside the grid"}

    # rasterise the centreline (dense sampling) and dilate to a pre-corridor
    road = np.zeros((r_hi - r_lo, c_hi - c_lo), dtype=bool)
    step = min(cw, chh) / 4.0
    for k in range(alignment.x0.size):
        n = max(2, int(math.ceil(alignment.seg_len[k] / step)) + 1)
        t = np.linspace(0.0, 1.0, n)
        px = alignment.x0[k] + t * alignment.dx[k]
        py = alignment.y0[k] + t * alignment.dy[k]
        cc = np.floor((px - gt[0]) / gt[1]).astype(np.int64) - c_lo
        rr = np.floor((py - gt[3]) / gt[5]).astype(np.int64) - r_lo
        ok = (rr >= 0) & (rr < road.shape[0]) & (cc >= 0) & (cc < road.shape[1])
        road[rr[ok], cc[ok]] = True
    k_cells = int(math.ceil((hw + 1.5 * max(cw, chh)) / min(cw, chh)))
    near = _dilate(road, k_cells)

    win_stream = (np.asarray(stream_mask, dtype=bool) & valid)[r_lo:r_hi, c_lo:c_hi] & near
    wr, wc = np.nonzero(win_stream)
    gr, gc = wr + r_lo, wc + c_lo                       # global row/col
    if gr.size == 0:
        return [], [], {"n_candidates": 0, "n_clusters": 0, "n_parallel": 0,
                        "note": "no stream cells near the alignment"}
    cx = gt[0] + (gc + 0.5) * gt[1]
    cy = gt[3] + (gr + 0.5) * gt[5]
    chain, off, _ = alignment.locate(cx, cy)

    d = direction[gr, gc]
    has = d >= 0
    jr = np.where(has, gr + DROW[np.where(has, d, 0)], -1)
    jc = np.where(has, gc + DCOL[np.where(has, d, 0)], -1)
    inside = has & (jr >= 0) & (jr < rows) & (jc >= 0) & (jc < cols)
    jx = gt[0] + (jc + 0.5) * gt[1]
    jy = gt[3] + (jr + 0.5) * gt[5]

    # -- 1. intersections of flow links with the alignment -------------
    li = np.flatnonzero(inside)
    link, seg, tpar, upar, hx, hy = alignment.intersect(cx[li], cy[li], jx[li], jy[li])
    raw = []
    for q in range(link.size):
        i = int(li[int(link[q])]); k = int(seg[q])
        r, c = int(gr[i]), int(gc[i])
        acc_km2 = (float(accumulation[r, c]) + 1.0) * cell_area / 1.0e6
        if acc_km2 < float(min_area_km2):
            continue
        fx, fy = jx[i] - cx[i], jy[i] - cy[i]
        fl = math.hypot(fx, fy)
        tx, ty = alignment.dx[k] / alignment.seg_len[k], alignment.dy[k] / alignment.seg_len[k]
        cosang = abs(fx * tx + fy * ty) / fl
        angle = math.degrees(math.acos(min(1.0, cosang)))
        side_cross = tx * (cy[i] - alignment.y0[k]) - ty * (cx[i] - alignment.x0[k])
        side_in = "L" if side_cross > 0 else "R"
        so = None
        if stream_order is not None:
            v = stream_order[r, c]
            so = int(v) if np.isfinite(v) and v > 0 else None
        rid = None
        if reach_ids is not None:
            v = reach_ids[r, c]
            rid = int(v) if np.isfinite(v) and v > 0 else None
        raw.append(dict(
            x=float(hx[q]), y=float(hy[q]),
            chainage_m=float(alignment.ch0[k] + upar[q] * alignment.seg_len[k]),
            acc_km2=acc_km2, stream_order=so, reach_id=rid,
            crossing_angle_deg=angle,
            road_azimuth_deg=math.degrees(math.atan2(tx, ty)) % 360.0,
            side_in=side_in,
            side_out="R" if side_in == "L" else "L",
            alignment_part=int(alignment.part[k]) + 1,
            outlet_row=r, outlet_col=c,
            outlet_x=float(cx[i]), outlet_y=float(cy[i]),
            _i=i, _j_rc=(int(jr[i]), int(jc[i]))))
    raw.sort(key=lambda a: (a["chainage_m"], -a["acc_km2"]))

    # -- 2. parallel runs inside the corridor --------------------------
    dist = neighbour_distances(cw, chh)
    in_corr = np.abs(off) <= hw
    idx_of = {(int(gr[i]), int(gc[i])): i for i in range(gr.size)}
    corr_idx = np.flatnonzero(in_corr)
    down = {}
    has_up = set()
    for i in corr_idx:
        if inside[i]:
            j = idx_of.get((int(jr[i]), int(jc[i])))
            if j is not None and in_corr[j]:
                down[int(i)] = j
                has_up.add(j)
    run_of = {}
    runs = []            # list of dict(cells=[i...], length, joins=set())
    for i in corr_idx:
        i = int(i)
        if i in has_up:
            continue
        rid = len(runs)
        cells, length, joins = [], 0.0, set()
        cur = i
        while True:
            if cur in run_of:
                joins.add(run_of[cur])
                break
            run_of[cur] = rid
            cells.append(cur)
            nxt = down.get(cur)
            if nxt is None:
                break
            length += float(dist[int(direction[gr[cur], gc[cur]])])
            cur = nxt
        runs.append({"cells": cells, "length": length, "joins": joins})
    uf_run = _UnionFind(len(runs))
    for rid, run in enumerate(runs):
        for j in run["joins"]:
            uf_run.union(rid, j)
    comp_parallel = {}
    for rid, run in enumerate(runs):
        root = uf_run.find(rid)
        comp_parallel[root] = comp_parallel.get(root, False) or run["length"] >= float(min_parallel_m)

    # -- 3. clusters ---------------------------------------------------
    n = len(raw)
    uf = _UnionFind(n)
    comp_first = {}
    for q, a in enumerate(raw):
        comps = set()
        for cell in (a["_i"], idx_of.get(a["_j_rc"])):
            if cell is not None and cell in run_of:
                comps.add(uf_run.find(run_of[cell]))
        a["_comps"] = comps
        a["parallel_reach"] = int(any(comp_parallel.get(cp, False) for cp in comps))
        for cp in comps:
            if comp_parallel.get(cp, False):
                if cp in comp_first:
                    uf.union(q, comp_first[cp])
                else:
                    comp_first[cp] = q
    if merge_distance_m and merge_distance_m > 0:
        for q in range(1, n):
            if raw[q]["chainage_m"] - raw[q - 1]["chainage_m"] < float(merge_distance_m):
                uf.union(q, q - 1)
    roots = {}
    for q in range(n):                      # cluster ids in chainage order
        r = uf.find(q)
        if r not in roots:
            roots[r] = len(roots) + 1
        raw[q]["cluster_id"] = roots[r]
    best = {}
    for q, a in enumerate(raw):
        key = (a["acc_km2"], a["chainage_m"])
        cid = a["cluster_id"]
        if cid not in best or key > best[cid][0]:
            best[cid] = (key, q)

    width = max(3, len(str(n)))
    candidates = []
    prev = None
    for q, a in enumerate(raw):
        attrs = {k: v for k, v in a.items() if not k.startswith("_") and k not in ("x", "y")}
        attrs.update(cand_id=f"{prefix}{q + 1:0{width}d}", status="candidate",
                     recommended=int(best[a["cluster_id"]][1] == q),
                     dist_prev_m=None if prev is None else a["chainage_m"] - prev)
        prev = a["chainage_m"]
        candidates.append(((a["x"], a["y"]), attrs))

    # -- 4. parallel reach lines ----------------------------------------
    comp_cluster = {}
    for a in raw:
        for cp in a["_comps"]:
            comp_cluster.setdefault(cp, a["cluster_id"])
    parallel = []
    for rid, run in enumerate(runs):
        if run["length"] < float(min_parallel_m) or len(run["cells"]) < 2:
            continue
        cells = run["cells"]
        coords = [(float(cx[i]), float(cy[i])) for i in cells]
        sides = {("L" if off[i] > 0 else "R") for i in cells if off[i] != 0}
        accs = [(float(accumulation[gr[i], gc[i]]) + 1.0) * cell_area / 1e6 for i in cells]
        parallel.append((coords, dict(
            pr_id=len(parallel) + 1, length_m=run["length"],
            chainage_from_m=float(min(chain[i] for i in cells)),
            chainage_to_m=float(max(chain[i] for i in cells)),
            side="".join(sorted(sides)) if sides else "C",
            acc_max_km2=max(accs),
            cluster_id=comp_cluster.get(uf_run.find(rid)))))

    summary = {"n_candidates": len(candidates), "n_clusters": len(best),
               "n_recommended": sum(a["recommended"] for _, a in candidates),
               "n_parallel": len(parallel),
               "alignment_length_m": alignment.length}
    return candidates, parallel, summary


def select_crossings(features, status_field="status", recommended_field="recommended"):
    """Which candidate features become crossings.

    If any feature is 'accepted', exactly the accepted ones; otherwise the
    recommended ones. Rejected features are never used. Returns
    (selected_indices, rule_text).
    """
    status = [str(f.get(status_field) or "").strip().lower() for f in features]
    acc = [k for k, s in enumerate(status) if s == "accepted"]
    if acc:
        return acc, "status = accepted"
    rec = [k for k, f in enumerate(features)
           if int(f.get(recommended_field) or 0) == 1 and status[k] != "rejected"]
    return rec, "recommended = 1 (no candidate marked accepted)"
