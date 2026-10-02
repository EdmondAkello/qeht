# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Flow-path segments (F2) and basin shape / network indices (F6), v0.15.

F2 - overland / channel split of the longest flow path
------------------------------------------------------
Walking the LFP from the divide, the CHANNEL HEAD is the first cell whose
accumulation reaches the channel threshold - the same threshold as the
stream network (stream_threshold_cells in the run metadata). The path is
split there:

    lfp_overland_m        divide -> channel head (planimetric)
    lfp_overland_slope    (z_head - z_channel_head) / lfp_overland_m
    lfp_channel_m         channel head -> outlet
    lfp_channel_slope     (z_channel_head - z_outlet) / lfp_channel_m
    lfp_channel_slope_1085  10-85 slope of the channel part, outlet-referenced
    lfp_sheet_m           min(lfp_overland_m, sheet cap); TR-55 practice
                          caps sheet flow at about 100 m
    lfp_shallow_m         lfp_overland_m - lfp_sheet_m (shallow concentrated)
    lfp_no_channel        1 when no path cell reaches the threshold (the
                          whole path is overland; channel fields 0 / empty)

lfp_overland_m + lfp_channel_m = lfp_length_m exactly. Slopes use the raw
DEM at the segment end cells (headwater, channel head, outlet).

F6 - basin shape and network indices (information only)
-------------------------------------------------------
    perimeter_km       outline length: 3x3-smoothed mask contoured at 0.5
                       (marching squares, interpolated) - within about 1-2 %
                       of the true outline, not the cell-edge staircase
                       (+27 % on round shapes)
    form_factor        A / L^2           (Horton 1932), L = LFP length
    elongation_ratio   (2 / L) sqrt(A/pi) (Schumm 1956)
    circularity_ratio  4 pi A / P^2       (Miller 1953)
    drainage_density   channel length / A (km/km2), channel cells at the
                       channel threshold, length along D8 links inside the
                       catchment
    stream_frequency   channel links / A (1/km2); a link runs from a source
                       or a confluence to the next confluence or the outlet
    max_strahler       highest Strahler order of the channel cells

No QGIS imports, no GDAL.
"""

import math

import numpy as np

from ..grid import DROW, DCOL, NO_RECEIVER, neighbour_distances
from .statistics import path_10_85, _cumulative_from_divide

NAN = float("nan")

LFP_SPLIT_FIELDS = ["lfp_overland_m", "lfp_overland_slope", "lfp_channel_m",
                    "lfp_channel_slope", "lfp_channel_slope_1085", "lfp_sheet_m",
                    "lfp_shallow_m", "lfp_threshold_km2", "lfp_no_channel"]
SHAPE_FIELDS = ["perimeter_km", "form_factor", "elongation_ratio", "circularity_ratio",
                "drainage_density", "stream_frequency", "max_strahler"]


def lfp_split(cells, elevation, accumulation, threshold_cells, cell_width, cell_height,
              sheet_cap_m=100.0):
    """F2 fields for one longest flow path (`cells` divide -> outlet)."""
    out = {k: None for k in LFP_SPLIT_FIELDS}
    out["lfp_threshold_km2"] = float(threshold_cells) * abs(cell_width * cell_height) / 1e6
    if not cells or len(cells) < 2:
        return out
    cum = _cumulative_from_divide(cells, cell_width, cell_height)
    total = float(cum[-1])
    acc = np.array([accumulation[rc] for rc in cells], dtype=float)
    z = np.array([elevation[rc] for rc in cells], dtype=float)
    hits = np.flatnonzero(np.nan_to_num(acc, nan=-1.0) >= float(threshold_cells))
    if hits.size == 0:
        h = len(cells) - 1
        out["lfp_no_channel"] = 1
    else:
        h = int(hits[0])
        out["lfp_no_channel"] = 0
    over = float(cum[h])
    chan = total - over
    out["lfp_overland_m"] = over
    out["lfp_channel_m"] = chan
    if over > 0 and np.isfinite(z[0]) and np.isfinite(z[h]):
        out["lfp_overland_slope"] = (z[0] - z[h]) / over
    if chan > 0 and np.isfinite(z[h]) and np.isfinite(z[-1]):
        out["lfp_channel_slope"] = (z[h] - z[-1]) / chan
        if len(cells) - h >= 2:
            s = path_10_85(cells[h:], elevation, cell_width, cell_height, reference="outlet")
            out["lfp_channel_slope_1085"] = s["slope"] if np.isfinite(s["slope"]) else None
    cap = max(float(sheet_cap_m), 0.0)
    out["lfp_sheet_m"] = min(over, cap)
    out["lfp_shallow_m"] = over - out["lfp_sheet_m"]
    return out


def contour_perimeter(mask, cell_width, cell_height):
    """Outline length of a bool mask, close to the true boundary length.

    The mask is smoothed with a 3x3 mean and contoured at 0.5 by marching
    squares with linear interpolation along the cell edges. A plain cell-edge
    staircase overstates a round outline by about 27 % (4/pi) and an
    unsmoothed mid-point contour by about 6 %; this estimate is within about
    1 % on round outlines and 1.5 % on squares (corners are rounded) a few
    dozen cells across. Holes add their
    own outline, as a polygon perimeter would.
    """
    m = np.pad(np.asarray(mask, bool).astype(float), 2)
    k = np.zeros_like(m)
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            k += np.roll(np.roll(m, dr, 0), dc, 1)
    f = k / 9.0
    a = f[:-1, :-1]; b = f[:-1, 1:]; c = f[1:, 1:]; d = f[1:, :-1]
    lv = 0.5
    code = ((a >= lv) * 8 + (b >= lv) * 4 + (c >= lv) * 2 + (d >= lv)).astype(np.int8)

    def _t(p, q):
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (lv - p) / (q - p)
        return np.clip(np.nan_to_num(t, nan=0.5), 0.0, 1.0)
    # crossing points on the four edges, unit square (x right, y down)
    top = (_t(a, b), np.zeros_like(a))
    right = (np.ones_like(a), _t(b, c))
    bottom = (_t(d, c), np.ones_like(a))
    left = (np.zeros_like(a), _t(a, d))
    edges = {"T": top, "R": right, "B": bottom, "L": left}
    pairs = {1: ["LB"], 2: ["BR"], 3: ["LR"], 4: ["TR"], 5: ["LT", "BR"], 6: ["TB"],
             7: ["LT"], 8: ["LT"], 9: ["TB"], 10: ["TR", "LB"], 11: ["TR"], 12: ["LR"],
             13: ["BR"], 14: ["LB"]}
    cw, ch = abs(cell_width), abs(cell_height)
    total = 0.0
    for cd, segs in pairs.items():
        sel = code == cd
        if not sel.any():
            continue
        for e1, e2 in segs:
            x1, y1 = edges[e1][0][sel], edges[e1][1][sel]
            x2, y2 = edges[e2][0][sel], edges[e2][1][sel]
            total += float(np.hypot((x2 - x1) * cw, (y2 - y1) * ch).sum())
    return total


def channel_network_stats(mask, direction, valid, channel_mask, cell_width, cell_height,
                          strahler=None):
    """(channel length m, number of links, max Strahler) inside `mask`."""
    rows, cols = direction.shape
    chm = np.asarray(channel_mask, bool) & mask & valid
    rr, cc = np.nonzero(chm)
    if rr.size == 0:
        return 0.0, 0, None
    dist = neighbour_distances(cell_width, cell_height)
    d = direction[rr, cc].astype(np.int64)
    ok = d != NO_RECEIVER
    tr = np.where(ok, rr + DROW[np.clip(d, 0, 7)], -1)
    tc = np.where(ok, cc + DCOL[np.clip(d, 0, 7)], -1)
    inside = ok & (tr >= 0) & (tr < rows) & (tc >= 0) & (tc < cols)
    inside[inside] = chm[tr[inside], tc[inside]]
    length = float(np.sum(np.where(inside, np.asarray(dist)[np.clip(d, 0, 7)], 0.0)))
    donors = np.zeros((rows, cols), dtype=np.int32)
    np.add.at(donors, (tr[inside], tc[inside]), 1)
    nd = donors[rr, cc]
    links = int(np.sum(nd == 0) + np.sum(nd >= 2))
    mx = None
    if strahler is not None:
        so = np.asarray(strahler, float)[rr, cc]
        so = so[np.isfinite(so) & (so > 0)]
        mx = int(so.max()) if so.size else None
    return length, links, mx


def shape_indices(mask, area_km2, lfp_length_m, cell_width, cell_height, direction=None,
                  valid=None, channel_mask=None, strahler=None):
    """F6 fields for one catchment (see the module docstring)."""
    out = {k: None for k in SHAPE_FIELDS}
    if not area_km2 or area_km2 <= 0:
        return out
    a_m2 = area_km2 * 1e6
    p = contour_perimeter(mask, cell_width, cell_height)
    out["perimeter_km"] = p / 1000.0
    if p > 0:
        out["circularity_ratio"] = 4.0 * math.pi * a_m2 / p ** 2
    if lfp_length_m and lfp_length_m > 0:
        out["form_factor"] = a_m2 / lfp_length_m ** 2
        out["elongation_ratio"] = (2.0 / lfp_length_m) * math.sqrt(a_m2 / math.pi)
    if direction is not None and channel_mask is not None:
        v = valid if valid is not None else np.ones(direction.shape, bool)
        length, links, mx = channel_network_stats(mask, direction, v, channel_mask,
                                                  cell_width, cell_height, strahler)
        out["drainage_density"] = (length / 1000.0) / area_km2
        out["stream_frequency"] = links / area_km2
        out["max_strahler"] = mx
    return out
