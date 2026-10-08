# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Approach / exit channel cross-section at a crossing (F5, v0.19).

A terrain-only first look at the natural channel downstream of a crossing,
for the tailwater rating, the bridge waterway and the sediment-continuity
check. Indicative: on a 30 m DEM a small channel is below the grid
resolution; use survey where available.

1. Location. Walk down the D8 receivers from the outlet for `dist_m`
   (cell or cell x sqrt 2 per step) and stop at the grid edge or NoData.
   The section centre is the centre of the cell reached.
2. Local flow direction. The vector from the cell k = 2 steps upstream to
   k = 2 steps downstream of the section cell (whatever exists of them),
   which smooths the 45 degree steps of D8. The section is perpendicular.
3. Transect. The raw DEM sampled bilinearly every `step_m` from -half_m to
   +half_m; offset positive on the RIGHT looking downstream.
4. Bed. The lowest sample within one cell of the centre.
5. Bank tops by break of slope, walking outward from the bed on each side:
   at each station the slope over the two cells behind it (towards the bed)
   is compared with the slope over the two cells ahead. A station is a
   break when the slope behind is at least `bank_slope` and the slope ahead
   is below `bank_slope` or less than half the slope behind. The bank top
   is the sharpest break (largest drop in slope) of the first run of
   breaks. Bank-full depth comes from the lower bank; the width is the
   contiguous run below that level (ends interpolated, as in A4).
6. Widths at bed + 0.5, 1 and 2 m by the same run.
7. Side slopes H:V by least squares of offset on height over the samples
   between the bed and the bank top (or bed + 1 m without a bank).
8. Quality: low for a bank-full width under 3 cells or no bank on either
   side; medium for one bank only or a transect cut by NoData or the grid
   edge (or a walk cut short); high otherwise. xs_note says why.

No QGIS imports, no GDAL.
"""

import json
import math

import numpy as np

from ..grid import DROW, DCOL, NO_RECEIVER, neighbour_distances
from ..network.profile import bilinear
from ..network.floodplain import run_width

DZ = (0.5, 1.0, 2.0)
XS_FIELDS = ["xs_dist_m", "xs_bed_m", "xs_bed_offset_m", "xs_bankfull_w_m", "xs_bankfull_d_m",
             "xs_w_0p5_m", "xs_w_1p0_m", "xs_w_2p0_m", "xs_side_slope_l", "xs_side_slope_r",
             "xs_station_elev_json", "xs_quality", "xs_note"]
HELP = ("Indicative only; on a 30 m DEM a small channel is below the grid resolution. "
        "Use survey where available.")


def _key(dz):
    return f"{dz:.1f}".replace(".", "p")


def _walk(direction, valid, r, c, n_steps):
    """Up to n_steps cells down the D8 receivers from (r, c), excluding (r, c)."""
    rows, cols = direction.shape
    out, a, b, seen = [], r, c, {(r, c)}
    for _ in range(n_steps):
        k = int(direction[a, b])
        if k == NO_RECEIVER or k < 0:
            break
        na, nb = a + int(DROW[k]), b + int(DCOL[k])
        if not (0 <= na < rows and 0 <= nb < cols) or not valid[na, nb] or (na, nb) in seen:
            break
        a, b = na, nb
        seen.add((a, b))
        out.append((a, b))
    return out


def section_location(direction, valid, r, c, dist_m, cell_width, cell_height):
    """Walk down from (r, c) for dist_m -> (path [(row, col, dist)], reached_m, cut_short)."""
    rows, cols = direction.shape
    dist = neighbour_distances(cell_width, cell_height)
    path = [(r, c, 0.0)]
    d, a, b, seen = 0.0, r, c, {(r, c)}
    cut = False
    while d < dist_m - 1e-9:
        k = int(direction[a, b])
        na, nb = (a + int(DROW[k]), b + int(DCOL[k])) if k >= 0 else (-1, -1)
        if k == NO_RECEIVER or k < 0 or not (0 <= na < rows and 0 <= nb < cols) \
                or not valid[na, nb] or (na, nb) in seen:
            cut = True
            break
        d += float(dist[k])
        a, b = na, nb
        seen.add((a, b))
        path.append((a, b, d))
    return path, d, cut


def _cell_xy(gt, r, c):
    return gt[0] + (c + 0.5) * gt[1] + (r + 0.5) * gt[2], gt[3] + (c + 0.5) * gt[4] + (r + 0.5) * gt[5]


def _bank(s, z, i0, side, base_n, bank_slope):
    """Bank top index on one side of the bed index i0 (side -1 = left, +1 = right).

    -> (index or None, hit_nodata)."""
    n = len(s)
    run, best, best_drop = False, None, -math.inf
    j = i0 + side
    while 0 <= j < n:
        if not np.isfinite(z[j]):
            return best, best is None
        jb = j - side * base_n                       # behind: towards the bed, not past it
        jb = max(jb, i0) if side > 0 else min(jb, i0)
        ja = j + side * base_n                       # ahead: outward, clamped to the data
        ja = min(ja, n - 1) if side > 0 else max(ja, 0)
        k = j
        while k != ja and np.isfinite(z[k + side]):
            k += side
        if k == j:                                    # nothing ahead: end of transect or NoData
            return best, best is None and 0 <= j + side < n
        behind = (z[j] - z[jb]) / abs(s[j] - s[jb])
        ahead = (z[k] - z[j]) / abs(s[k] - s[j])
        brk = behind >= bank_slope and (ahead < bank_slope or ahead < 0.5 * behind)
        if brk:
            run = True
            if behind - ahead > best_drop:
                best, best_drop = j, behind - ahead
        elif run:
            return best, False
        j += side
    return best, False


def _side_slope(s, z, i0, side, top, z_bed):
    """H:V by least squares of |offset| on height, bed -> top (exclusive of the flat bed)."""
    lim = 0.01 * max(top - z_bed, 1e-9)
    xs, zs = [], []
    j = i0 + side
    while 0 <= j < len(s) and np.isfinite(z[j]) and z[j] <= top + 1e-9:
        if z[j] > z_bed + lim:
            xs.append(abs(s[j] - s[i0]))
            zs.append(z[j])
        j += side
    if len(zs) < 2:
        return None
    zs, xs = np.array(zs), np.array(xs)
    dz = zs - zs.mean()
    den = float((dz * dz).sum())
    if den <= 0:
        return None
    hv = float((dz * (xs - xs.mean())).sum() / den)
    return abs(hv) if math.isfinite(hv) else None


def section_at_outlet(direction, valid, raw, geotransform, r, c, dist_m=30.0, half_m=150.0,
                      step_m=None, bank_slope=0.05, dzs=DZ):
    """F5 fields for the crossing whose outlet cell is (r, c).

    raw: raw DEM as float with NaN for NoData. Returns (fields, transect)
    where transect is {"line": [(x0, y0), (x1, y1)], "stations": [[s, z], ...],
    "centre": (x, y), "flow": (fx, fy)}, or None when no section can be cut.
    """
    gt = geotransform
    cw, ch = abs(gt[1]), abs(gt[5])
    cell = min(cw, ch)
    step = float(step_m) if step_m else 0.5 * cell
    out = {k: None for k in XS_FIELDS}
    notes = []
    path, reached, cut = section_location(direction, valid, r, c, dist_m, cw, ch)
    out["xs_dist_m"] = float(reached)
    if cut and reached < dist_m - 1e-9:
        notes.append(f"walk stopped after {reached:.0f} m of {dist_m:g} m (grid edge, NoData "
                     "or no receiver)")
    sr, sc = path[-1][0], path[-1][1]
    up = [(a, b) for a, b, _ in path[:-1]][-2:]
    down = _walk(direction, valid, sr, sc, 2)
    p0 = up[0] if up else (sr, sc)
    p1 = down[-1] if down else (sr, sc)
    x0, y0 = _cell_xy(gt, *p0)
    x1, y1 = _cell_xy(gt, *p1)
    fx, fy = x1 - x0, y1 - y0
    norm = math.hypot(fx, fy)
    if norm == 0:
        out["xs_quality"] = "low"
        out["xs_note"] = "; ".join(notes + ["no flow direction at the section cell"])
        return out, None
    fx, fy = fx / norm, fy / norm
    rx, ry = fy, -fx                                   # right of the flow (map y north)
    cx, cy = _cell_xy(gt, sr, sc)
    n_half = int(round(half_m / step))
    s = np.arange(-n_half, n_half + 1, dtype=float) * step
    z = bilinear(raw, gt, cx + s * rx, cy + s * ry)
    ic = n_half
    near = np.where((np.abs(s) <= cell + 1e-9) & np.isfinite(z))[0]
    transect = {"line": [(cx + s[0] * rx, cy + s[0] * ry), (cx + s[-1] * rx, cy + s[-1] * ry)],
                "stations": [[round(float(a), 2), round(float(b), 2)]
                             for a, b in zip(s, z) if np.isfinite(b)],
                "centre": (cx, cy), "flow": (fx, fy)}
    out["xs_station_elev_json"] = json.dumps(transect["stations"], separators=(",", ":"))
    if near.size == 0:
        out["xs_quality"] = "low"
        out["xs_note"] = "; ".join(notes + ["no DEM value at the section centre"])
        return out, transect
    i0 = int(near[np.lexsort((np.abs(s[near]), z[near]))[0]])
    z_bed = float(z[i0])
    out["xs_bed_m"] = z_bed
    out["xs_bed_offset_m"] = float(s[i0])
    base_n = max(1, int(round(2.0 * cell / step)))
    il, nd_l = _bank(s, z, i0, -1, base_n, bank_slope)
    ir, nd_r = _bank(s, z, i0, +1, base_n, bank_slope)
    truncated = nd_l or nd_r

    def width(level, what):
        nonlocal truncated
        a, b, tl, tr = run_width(s, z, i0, level)
        if a is None:
            return None
        for t, end in ((tl, a), (tr, b)):
            if t:
                if abs(end) >= s[-1] - 1e-9:
                    notes.append(f"{what} reaches the transect end (+/-{half_m:g} m; lower bound)")
                else:
                    truncated = True
                    notes.append(f"{what} stops at NoData or the grid edge (lower bound)")
        return b - a

    tops = [float(z[i]) for i in (il, ir) if i is not None]
    if tops:
        d = min(tops) - z_bed
        out["xs_bankfull_d_m"] = d
        out["xs_bankfull_w_m"] = width(z_bed + d, "bank-full width")
    for dz in dzs:
        out[f"xs_w_{_key(dz)}_m"] = width(z_bed + dz, f"width at bed + {dz:g} m")
    out["xs_side_slope_l"] = _side_slope(s, z, i0, -1, float(z[il]) if il is not None
                                         else z_bed + 1.0, z_bed)
    out["xs_side_slope_r"] = _side_slope(s, z, i0, +1, float(z[ir]) if ir is not None
                                         else z_bed + 1.0, z_bed)
    # quality
    q = "high"
    if il is None and ir is None:
        q = "low"
        notes.append("no bank found on either side (no break of slope)")
    elif out["xs_bankfull_w_m"] is not None and out["xs_bankfull_w_m"] < 3.0 * cell:
        q = "low"
        notes.append(f"bank-full width {out['xs_bankfull_w_m']:.1f} m is under 3 cells "
                     "(sub-grid channel)")
    else:
        if il is None or ir is None:
            q = "medium"
            notes.append(f"bank found on the {'right' if il is None else 'left'} side only")
        if truncated:
            q = "medium"
            if not any("NoData" in x for x in notes):
                notes.append("transect cut by NoData or the grid edge")
        if cut and reached < dist_m - 1e-9:
            q = "medium"
    out["xs_quality"] = q
    out["xs_note"] = "; ".join(dict.fromkeys(notes)) or None
    return out, transect


def outlet_cell(geotransform, x, y, shape):
    """Row, column of the cell containing (x, y) (north-up grid), or None outside."""
    gt = geotransform
    c = int(math.floor((x - gt[0]) / gt[1]))
    r = int(math.floor((y - gt[3]) / gt[5]))
    if 0 <= r < shape[0] and 0 <= c < shape[1]:
        return r, c
    return None


def section_block(crossings, direction, valid, raw, geotransform, dist_m=30.0, half_m=150.0,
                  step_m=None, bank_slope=0.05, dzs=DZ):
    """XS fields for each crossing attrs dict (needs outlet_x / outlet_y).

    -> (list of field dicts in the order of `crossings`, transect layer rows
    [(line, {outlet_uid, xs_dist_m, xs_station_elev_json, xs_quality})])."""
    out, lines = [], []
    for a in crossings:
        b = {k: None for k in XS_FIELDS}
        rc = None
        if a.get("outlet_x") is not None and a.get("outlet_y") is not None:
            rc = outlet_cell(geotransform, a["outlet_x"], a["outlet_y"], direction.shape)
        if rc is None or not valid[rc]:
            b["xs_quality"] = "low"
            b["xs_note"] = "outlet outside the grid"
            out.append(b)
            continue
        f, t = section_at_outlet(direction, valid, raw, geotransform, rc[0], rc[1], dist_m,
                                 half_m, step_m, bank_slope, dzs)
        b.update(f)
        out.append(b)
        if t is not None:
            lines.append((t["line"], {"outlet_uid": a.get("outlet_uid"),
                                      "xs_dist_m": b["xs_dist_m"],
                                      "xs_station_elev_json": b["xs_station_elev_json"],
                                      "xs_quality": b["xs_quality"]}))
    return out, lines


def params_json(dist_m, half_m, step_m, bank_slope, dzs=DZ):
    return json.dumps({"dist_m": dist_m, "half_m": half_m,
                       "step_m": step_m if step_m else "half a cell",
                       "bank_slope": bank_slope, "bank_rule": "break of slope over two cells: "
                       "slope ahead < bank_slope or < 0.5 x slope behind; lower bank governs",
                       "dz_m": list(dzs), "note": HELP}, sort_keys=True)
