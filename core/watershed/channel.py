# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Approach and exit channel along the D8 network (A5, v0.15).

main_stem_upstream walks UP from a crossing: at each cell the upstream
neighbour (a cell draining into it) with the largest accumulation is taken,
i.e. the main stem, never a side tributary. downstream_path follows the D8
receivers DOWN from the crossing. Both stop at `max_len_m` of planimetric
length (diagonal steps count as cell diagonals), at the grid edge, at NoData
or - upstream - where no cell drains in.

channel_slopes gives the representative approach (upstream) and exit
(downstream) slopes a downstream floodplain or tailwater calculation needs:

    ch_slope_us = (z at the far upstream end - z at the crossing) / length
    ch_slope_ds = (z at the crossing - z at the far downstream end) / length

with z from the raw DEM. ch_len_us_m / ch_len_ds_m give the length actually
used (shorter than requested near the divide or the grid edge).

No QGIS imports, no GDAL.
"""

import math

import numpy as np

from ..grid import DROW, DCOL, NO_RECEIVER, neighbour_distances

CHANNEL_FIELDS = ["ch_slope_us", "ch_slope_ds", "ch_len_us_m", "ch_len_ds_m", "ch_slope_dist_m"]


def main_stem_upstream(direction, valid, accumulation, r, c, max_len_m, cell_width,
                       cell_height, mask=None):
    """[(row, col, distance_from_start_m), ...] starting with (r, c) at 0."""
    rows, cols = direction.shape
    dist = neighbour_distances(cell_width, cell_height)
    out = [(r, c, 0.0)]
    d = 0.0
    a, b = r, c
    while d < max_len_m:
        best, best_acc, best_k = None, -1.0, None
        for k in range(8):
            ua, ub = a - int(DROW[k]), b - int(DCOL[k])
            if not (0 <= ua < rows and 0 <= ub < cols) or not valid[ua, ub]:
                continue
            if mask is not None and not mask[ua, ub]:
                continue
            if int(direction[ua, ub]) != k:
                continue
            acc = float(accumulation[ua, ub])
            if acc > best_acc:
                best, best_acc, best_k = (ua, ub), acc, k
        if best is None:
            break
        d += float(dist[best_k])
        a, b = best
        out.append((a, b, d))
    return out


def downstream_path(direction, valid, r, c, max_len_m, cell_width, cell_height):
    rows, cols = direction.shape
    dist = neighbour_distances(cell_width, cell_height)
    out = [(r, c, 0.0)]
    d = 0.0
    a, b = r, c
    seen = {(a, b)}
    while d < max_len_m:
        k = int(direction[a, b])
        if k == NO_RECEIVER or k < 0:
            break
        na, nb = a + int(DROW[k]), b + int(DCOL[k])
        if not (0 <= na < rows and 0 <= nb < cols) or not valid[na, nb] or (na, nb) in seen:
            break
        d += float(dist[k])
        a, b = na, nb
        seen.add((a, b))
        out.append((a, b, d))
    return out


def _slope(path, elevation, up):
    if len(path) < 2:
        return None, 0.0
    z0 = float(elevation[path[0][0], path[0][1]])
    z1 = float(elevation[path[-1][0], path[-1][1]])
    L = path[-1][2]
    if L <= 0 or not (math.isfinite(z0) and math.isfinite(z1)):
        return None, L
    return ((z1 - z0) if up else (z0 - z1)) / L, L


def channel_slopes(direction, valid, accumulation, elevation, r, c, cell_width, cell_height,
                   distance_m=200.0, mask=None):
    """A5 fields for one crossing at (r, c)."""
    up = main_stem_upstream(direction, valid, accumulation, r, c, distance_m, cell_width,
                            cell_height, mask)
    dn = downstream_path(direction, valid, r, c, distance_m, cell_width, cell_height)
    s_up, l_up = _slope(up, elevation, True)
    s_dn, l_dn = _slope(dn, elevation, False)
    return {"ch_slope_us": s_up, "ch_slope_ds": s_dn, "ch_len_us_m": l_up,
            "ch_len_ds_m": l_dn, "ch_slope_dist_m": float(distance_m)}


def deposition_indicator(up_path, spi, tan_beta, channel, near_m=100.0, far_m=500.0,
                         breaks=(0.7, 1.3), elevation=None):
    """STI advisory R2: change in stream power into the crossing.

    up_path from main_stem_upstream (starting at the crossing). Near reach:
    channel cells with 0 < distance <= near_m; far reach: near_m < distance
    <= far_m. ratio = median SPI near / median SPI far; below breaks[0]
    capacity falls into the crossing (deposition-prone), above breaks[1] it
    rises (scour-prone). Slopes of both reaches are exported too (drop /
    length on the raw DEM when given, else mean Horn slope), because over a
    short reach A changes little and the ratio is mostly the slope break.
    NULL with a note when the main stem leaves the channel network or the
    grid before far_m.
    """
    out = {"ero_spi_app_near": None, "ero_spi_app_far": None, "ero_dep_ratio": None,
           "ero_dep_flag": None, "ero_slope_app_near_pct": None, "ero_slope_app_far_pct": None,
           "ero_dep_note": None}
    # stop at the first non-channel cell on the main stem
    stem = []
    for a, b, d in up_path[1:]:
        if not channel[a, b]:
            break
        stem.append((a, b, d))
    near = [(a, b, d) for a, b, d in stem if d <= near_m]
    far = [(a, b, d) for a, b, d in stem if near_m < d <= far_m]
    reach_end = stem[-1][2] if stem else 0.0
    if not near or not far or reach_end < 0.9 * far_m:
        out["ero_dep_note"] = (f"main channel upstream only {reach_end:.0f} m "
                               f"(needs {far_m:.0f} m) - no ratio")
        return out
    sn = np.array([spi[a, b] for a, b, _ in near], float)
    sf = np.array([spi[a, b] for a, b, _ in far], float)
    sn, sf = sn[np.isfinite(sn)], sf[np.isfinite(sf)]
    if not sn.size or not sf.size:
        out["ero_dep_note"] = "no SPI on the approach channel"
        return out
    mn, mf = float(np.median(sn)), float(np.median(sf))
    out["ero_spi_app_near"], out["ero_spi_app_far"] = mn, mf

    def slope(cells, start):
        if elevation is not None:
            z0 = float(elevation[start[0], start[1]])
            z1 = float(elevation[cells[-1][0], cells[-1][1]])
            L = cells[-1][2] - start[2]
            if L > 0 and math.isfinite(z0) and math.isfinite(z1):
                return 100.0 * (z1 - z0) / L
        t = np.array([tan_beta[a, b] for a, b, _ in cells], float)
        t = t[np.isfinite(t)]
        return float(100.0 * t.mean()) if t.size else None
    out["ero_slope_app_near_pct"] = slope(near, up_path[0])
    out["ero_slope_app_far_pct"] = slope(far, near[-1])
    if mf > 0:
        ratio = mn / mf
        out["ero_dep_ratio"] = ratio
        lo, hi = breaks
        out["ero_dep_flag"] = ("deposition-prone" if ratio < lo else
                               "scour-prone" if ratio > hi else "neutral")
    return out
