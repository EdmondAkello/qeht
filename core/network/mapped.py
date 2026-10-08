# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Check the DEM drainage against mapped drainage (F12, v0.24).

The independent check the flat-sensitivity finding points to: did the DEM
put the river where the map does? Mapped waterways (OSM, HydroRIVERS, a
national layer) are compared with the DEM streams at a comparison
threshold matched to the map's scale (default 1 km2; about 10 km2 for
HydroRIVERS).

1. Network agreement. A DEM stream cell agrees when a mapped line passes
   within the tolerance (default 60 m) of its centre; a point on a mapped
   line agrees when a DEM stream cell centre lies within the tolerance.
   precision = agreeing DEM stream length / DEM stream length (D8 link
   lengths); recall = agreeing mapped length / mapped length inside the
   grid (samples every half cell); F1 = 2PR / (P + R). Overall and per
   catchment (the catchment's cells).
2. Per crossing: map_river_dist_m (outlet to the nearest mapped line),
   map_river_name, map_agrees = 1 within the tolerance, 0 within 1 km,
   NULL beyond 1 km or for a crossing smaller than the comparison
   threshold (the map is not expected to show it; map_note says so).
3. Road: a mapped river crossing the road with no crossing within 50 m of
   chainage is a coverage finding mapped_river_uncovered (the designer
   decides; no proposed crossing is created).
4. Divergence reaches: a run of DEM stream cells that leaves agreement
   (its upstream stream cell agrees, it does not) and stays out for at
   least `min_cells` cells along D8. Length, cells and the area that
   flows down it (accumulation at its first cell) - how a catchment that
   switches between crossings shows up.

No QGIS imports, no GDAL.
"""

import json
import math

import numpy as np

from ..grid import DROW, DCOL, receivers_from_direction, neighbour_distances

MAP_FIELDS_X = ["map_river_dist_m", "map_river_name", "map_agrees", "map_note"]
MAP_FIELDS_C = ["map_precision", "map_recall", "map_f1"]
DIVERGENCE_FIELDS = [("length_m", "real"), ("cells", "int"), ("area_km2", "real"),
                     ("rejoins", "int"), ("from_x", "real"), ("from_y", "real")]
MAPPED_FIELDS = [("name", "text"), ("length_m", "real")]


def densify(line, step):
    """Points every <= step along a polyline (vertices kept)."""
    a = np.asarray(line, float).reshape(-1, 2)
    out = [a[:1]]
    for p, q in zip(a[:-1], a[1:]):
        n = max(int(math.ceil(math.hypot(*(q - p)) / step)), 1)
        t = (np.arange(1, n + 1) / n)[:, None]
        out.append(p + t * (q - p))
    return np.vstack(out)


def _cells(gt, shape, xy):
    rows, cols = shape
    c = np.floor((xy[:, 0] - gt[0]) / gt[1]).astype(np.int64)
    r = np.floor((xy[:, 1] - gt[3]) / gt[5]).astype(np.int64)
    ok = (r >= 0) & (r < rows) & (c >= 0) & (c < cols)
    return r, c, ok


def dilate(mask, radius_cells, cw=1.0, ch=1.0):
    """Cells whose centre lies within radius (in units of cw) of a True cell centre."""
    m = np.asarray(mask, bool)
    out = m.copy()
    rr = int(math.floor(radius_cells * cw / ch)) + 1
    rc = int(math.floor(radius_cells)) + 1
    rows, cols = m.shape
    for dr in range(-rr, rr + 1):
        for dc in range(-rc, rc + 1):
            if (dr == 0 and dc == 0) or math.hypot(dr * ch, dc * cw) > radius_cells * cw + 1e-9:
                continue
            src = m[max(0, -dr):rows - max(0, dr), max(0, -dc):cols - max(0, dc)]
            out[max(0, dr):rows - max(0, -dr), max(0, dc):cols - max(0, -dc)] |= src
    return out


def mapped_mask(lines, gt, shape):
    """Cells touched by the mapped lines (sampled every quarter cell)."""
    step = 0.25 * min(abs(gt[1]), abs(gt[5]))
    m = np.zeros(shape, bool)
    for ln in lines:
        r, c, ok = _cells(gt, shape, densify(ln, step))
        m[r[ok], c[ok]] = True
    return m


def link_lengths(direction, valid, cw, ch):
    """Length of each cell's D8 link to its receiver (half a cell where none)."""
    dist = neighbour_distances(cw, ch)
    d = np.asarray(direction)
    return np.where((d >= 0) & valid, dist[np.clip(d, 0, 7)], 0.5 * (cw + ch) / 2.0)


def prepare(stream, direction, valid, lines, gt, tol_m):
    """Grids and line samples shared by every region's score."""
    shape = stream.shape
    cw, ch = abs(gt[1]), abs(gt[5])
    s = stream & valid
    ctx = {"stream": s, "near_map": dilate(mapped_mask(lines, gt, shape), tol_m / cw, cw, ch),
           "ll": link_lengths(direction, valid, cw, ch)}
    near_stream = dilate(s, tol_m / cw, cw, ch)
    step = 0.5 * min(cw, ch)
    rr, cc, ww, hit = [], [], [], []
    for ln in lines:
        pts = densify(ln, step)
        seg = np.r_[np.hypot(*np.diff(pts, axis=0).T), 0.0]
        w = 0.5 * (seg + np.r_[0.0, seg[:-1]])                 # length per sample
        r, c, ok = _cells(gt, shape, pts)
        ok[ok] &= valid[r[ok], c[ok]]
        rr.append(r[ok]); cc.append(c[ok]); ww.append(w[ok])
        hit.append(near_stream[r[ok], c[ok]])
    cat = (lambda a, t: np.concatenate(a) if a else np.zeros(0, t))
    ctx.update(r=cat(rr, np.int64), c=cat(cc, np.int64), w=cat(ww, float), hit=cat(hit, bool))
    return ctx


def scores(ctx, region=None):
    """(precision, recall, f1) inside region (bool grid) or overall."""
    s = ctx["stream"] if region is None else ctx["stream"] & region
    tot = float(ctx["ll"][s].sum())
    prec = float(ctx["ll"][s & ctx["near_map"]].sum()) / tot if tot > 0 else None
    sel = np.ones(ctx["w"].size, bool) if region is None else region[ctx["r"], ctx["c"]]
    got = float(ctx["w"][sel].sum())
    rec = float(ctx["w"][sel & ctx["hit"]].sum()) / got if got > 0 else None
    f1 = 2 * prec * rec / (prec + rec) if prec is not None and rec is not None and prec + rec > 0 else None
    return prec, rec, f1


def agreement(stream, direction, valid, lines, gt, tol_m, region=None):
    """(precision, recall, f1, agree_stream_mask) over `region` (bool grid, or all)."""
    ctx = prepare(stream, direction, valid, lines, gt, tol_m)
    reg = None if region is None else np.asarray(region, bool)
    p, r, f = scores(ctx, reg)
    return p, r, f, ctx["stream"] & ctx["near_map"]


def nearest_line(x, y, lines, names=None):
    """(distance, index) of the nearest mapped line to (x, y)."""
    best, bi = math.inf, None
    for i, ln in enumerate(lines):
        a = np.asarray(ln, float).reshape(-1, 2)
        if len(a) < 2:
            continue
        p, q = a[:-1], a[1:]
        d = q - p
        L2 = (d * d).sum(axis=1)
        t = np.clip(((x - p[:, 0]) * d[:, 0] + (y - p[:, 1]) * d[:, 1]) / np.where(L2 > 0, L2, 1), 0, 1)
        dist = np.hypot(p[:, 0] + t * d[:, 0] - x, p[:, 1] + t * d[:, 1] - y).min()
        if dist < best:
            best, bi = float(dist), i
    return best, bi


def crossing_fields(crossings, lines, names, tol_m, threshold_km2, far_m=1000.0):
    """MAP_FIELDS_X for each crossing attrs dict (outlet_x / outlet_y / acc_at_outlet_km2)."""
    out = []
    for a in crossings:
        b = {k: None for k in MAP_FIELDS_X}
        if a.get("outlet_x") is None or not lines:
            out.append(b)
            continue
        d, i = nearest_line(a["outlet_x"], a["outlet_y"], lines)
        if i is None:
            out.append(b)
            continue
        b["map_river_dist_m"] = d
        b["map_river_name"] = (names[i] if names and names[i] not in (None, "") else None)
        area = a.get("acc_at_outlet_km2")
        if area is not None and threshold_km2 and area < threshold_km2:
            b["map_note"] = (f"{area:.3g} km2 is below the comparison threshold "
                             f"({threshold_km2:g} km2): the map is not expected to show it")
        elif d > far_m:
            b["map_note"] = f"no mapped line within {far_m:g} m"
        else:
            b["map_agrees"] = int(d <= tol_m)
            if not b["map_agrees"]:
                b["map_note"] = f"nearest mapped line {d:.0f} m away (tolerance {tol_m:g} m)"
        out.append(b)
    return out


def uncovered_rivers(lines, names, alignment, crossing_chainages, search_m=50.0):
    """Mapped rivers crossing the road with no crossing within search_m of chainage.

    -> list of ((x, y), coverage attrs)."""
    out = []
    chs = np.array([c for c in crossing_chainages if c is not None], float)
    for i, ln in enumerate(lines):
        a = np.asarray(ln, float).reshape(-1, 2)
        if len(a) < 2:
            continue
        _, seg, _, u, xs, ys = alignment.intersect(a[:-1, 0], a[:-1, 1], a[1:, 0], a[1:, 1])
        for s_, u_, x_, y_ in zip(seg, u, xs, ys):
            ch_ = float(alignment.ch0[int(s_)] + u_ * alignment.seg_len[int(s_)])
            near = float(np.abs(chs - ch_).min()) if chs.size else math.inf
            if near <= search_m:
                continue
            nm = names[i] if names and names[i] not in (None, "") else "unnamed"
            out.append(((float(x_), float(y_)), {
                "issue": "mapped_river_uncovered", "chainage_m": ch_, "chainage_to_m": None,
                "area_km2": None, "uid": None, "nearest_uid": None,
                "nearest_m": near if math.isfinite(near) else None,
                "note": f"mapped river '{nm}' crosses the road with no crossing within "
                        f"{search_m:g} m (designer decides)"}))
    return out


def divergence_reaches(stream, agree, direction, valid, accumulation, gt, min_cells=10):
    """Runs of disagreeing stream cells that start below an agreeing stream cell.

    -> list of (line [(x, y)], attrs DIVERGENCE_FIELDS)."""
    shape = stream.shape
    rows, cols = shape
    cw, ch = abs(gt[1]), abs(gt[5])
    dist = neighbour_distances(cw, ch)
    d = np.where(valid, direction, -1)
    rec = receivers_from_direction(d, shape)
    s = (stream & valid).ravel()
    ag = (agree & s.reshape(shape)).ravel()
    dis = s & ~ag
    # heads: disagreeing stream cells with an agreeing stream donor
    donors_agree = np.zeros(s.size, bool)
    idx = np.flatnonzero(ag & (rec >= 0))
    donors_agree[rec[idx]] = True
    heads = np.flatnonzero(dis & donors_agree)
    acc = np.asarray(accumulation, float).ravel()
    cell_km2 = cw * ch / 1e6
    seen = np.zeros(s.size, bool)
    out = []
    for h in heads[np.argsort(-acc[heads])]:            # largest first: a main stem keeps its run
        path, cur, length = [], h, 0.0
        while cur >= 0 and dis[cur] and not seen[cur]:
            seen[cur] = True
            path.append(cur)
            nxt = rec[cur]
            if nxt >= 0:
                length += float(dist[int(d.ravel()[cur])])
            cur = nxt
        if len(path) < min_cells:
            continue
        rejoins = int(cur >= 0 and ag[cur])
        r_, c_ = np.divmod(np.array(path + ([cur] if cur >= 0 else [])), cols)
        line = [(gt[0] + (c + 0.5) * gt[1], gt[3] + (r + 0.5) * gt[5]) for r, c in zip(r_, c_)]
        out.append((line, {"length_m": length, "cells": len(path),
                           "area_km2": (acc[h] + 1.0) * cell_km2, "rejoins": rejoins,
                           "from_x": line[0][0], "from_y": line[0][1]}))
    return out


def clip_lines(lines, names, gt, shape):
    """Parts of the mapped lines inside the grid extent (vertex-level) -> layer rows."""
    rows, cols = shape
    x0, x1 = gt[0], gt[0] + cols * gt[1]
    y1, y0 = gt[3], gt[3] + rows * gt[5]
    step = 0.5 * min(abs(gt[1]), abs(gt[5]))
    out = []
    for i, ln in enumerate(lines):
        pts = densify(ln, step)
        inside = (pts[:, 0] >= x0) & (pts[:, 0] <= x1) & (pts[:, 1] >= y0) & (pts[:, 1] <= y1)
        run = []
        for p, ok in zip(pts, inside):
            if ok:
                run.append((float(p[0]), float(p[1])))
            elif run:
                if len(run) >= 2:
                    out.append((run, i))
                run = []
        if len(run) >= 2:
            out.append((run, i))
    rows_ = []
    for run, i in out:
        a = np.asarray(run)
        rows_.append((run[::4] + [run[-1]] if len(run) > 8 else run,
                      {"name": names[i] if names else None,
                       "length_m": float(np.hypot(*np.diff(a, axis=0).T).sum())}))
    return rows_


def params_json(source, threshold_km2, tol_m, min_cells, overall, n_lines, name_field):
    p, r, f = overall
    return json.dumps({"source": source, "name_field": name_field or "",
                       "comparison_threshold_km2": threshold_km2, "tolerance_m": tol_m,
                       "divergence_min_cells": min_cells, "lines": n_lines,
                       "precision": p, "recall": r, "f1": f, "far_m": 1000.0,
                       "road_search_m": 50.0}, sort_keys=True)
