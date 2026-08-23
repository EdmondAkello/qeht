# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Convergent flat resolution (Barnes, Lehman & Mulla 2014) - v0.8.

Source-faithful reconstruction of RichDEM's flat_resolution.hpp, with the
one structural change QEHT needs: iteration to convergence.

Why iteration
-------------
Barnes assumes flow directions on non-flat cells (including flat outlets)
are already assigned before flats are resolved. RichDEM guarantees this by
running its flow-direction pass first. In QEHT the outlet of a flat may
itself be an unrouted flat cell of a LOWER flat, so a single pass starves
any flat whose spill point is not yet resolved. The fix: repeat the Barnes
pass; each iteration resolves flats whose outlets became available in the
previous one. Converges in a handful of passes (depth of the flat
hierarchy), and every pass is O(N).

Definitions (from RichDEM find_flat_edges, docstring verbatim):
  low edge  : a flat cell adjacent to LOWER terrain, or to an equal-
              elevation cell that already drains (its spill point).
  high edge : a flat cell adjacent to HIGHER terrain.

Gradients (from BuildAwayGradient / BuildTowardsCombinedGradient):
  away   : BFS inward from high edges; per-flat max recorded as flat_height.
  toward : BFS inward from low edges, combined as
             flat_mask = 2*toward + (flat_height - away)
  Each flat cell then flows to the same-label neighbour of lowest flat_mask,
  cardinal preferred over diagonal on ties (d8_masked_FlowDir).

Reference: Barnes, R., Lehman, C., Mulla, D. (2014). Computers &
Geosciences 62, 128-135. doi:10.1016/j.cageo.2013.01.009
"""

import numpy as np
from collections import deque

from ..grid import DROW, DCOL

_IS_CARDINAL = np.array([True, False, True, False, True, False, True, False])


def _find_flat_edges(elev, valid, direction):
    """low_edges, high_edges: flat (no-flow) cells adjacent to lower /
    higher valid terrain respectively (RichDEM find_flat_edges)."""
    rows, cols = elev.shape
    no_flow = valid & (direction < 0)
    low, high = [], []
    for fid in np.flatnonzero(no_flow.reshape(-1)):
        r, c = divmod(int(fid), cols)
        e = elev[r, c]
        is_low = is_high = False
        for k in range(8):
            nr, nc = r + int(DROW[k]), c + int(DCOL[k])
            if nr < 0 or nr >= rows or nc < 0 or nc >= cols:
                continue
            if not valid[nr, nc]:
                continue
            if (not is_low) and (elev[nr, nc] < e - 1e-9 or
                                 (direction[nr, nc] >= 0 and abs(elev[nr, nc] - e) <= 1e-9)):
                low.append((r, c)); is_low = True
            if (not is_high) and elev[nr, nc] > e + 1e-9:
                high.append((r, c)); is_high = True
            if is_low and is_high:
                break
    return low, high


def _label_flats(elev, valid, direction, low_edges):
    """Flood-fill flats FROM low edges over equal-elevation no-flow cells.
    Only flats with a low edge (an outlet) get labelled; the rest stay 0
    and remain unrouted this pass. (RichDEM label_this, seeded from
    low_edges.)"""
    rows, cols = elev.shape
    no_flow = valid & (direction < 0)
    labels = np.zeros((rows, cols), dtype=np.int64)
    n = 0
    for (r0, c0) in low_edges:
        if labels[r0, c0] != 0:
            continue
        n += 1
        e0 = elev[r0, c0]
        q = deque([(r0, c0)])
        labels[r0, c0] = n
        while q:
            r, c = q.popleft()
            for k in range(8):
                nr, nc = r + int(DROW[k]), c + int(DCOL[k])
                if nr < 0 or nr >= rows or nc < 0 or nc >= cols:
                    continue
                if labels[nr, nc] != 0 or not no_flow[nr, nc]:
                    continue
                if abs(elev[nr, nc] - e0) <= 1e-9:
                    labels[nr, nc] = n
                    q.append((nr, nc))
    return labels, n


def _away_gradient(labels, direction, high_edges, n_labels):
    """BFS inward from high edges. RichDEM BuildAwayGradient."""
    rows, cols = labels.shape
    away = np.zeros((rows, cols), dtype=np.int64)
    flat_height = np.zeros(n_labels + 1, dtype=np.int64)
    edges = deque((r, c) for (r, c) in high_edges if labels[r, c] != 0)
    MARK = (-1, -1); edges.append(MARK); loops = 1
    while len(edges) > 1:
        r, c = edges.popleft()
        if r == -1:
            loops += 1; edges.append(MARK); continue
        if away[r, c] > 0:
            continue
        away[r, c] = loops
        flat_height[labels[r, c]] = loops
        for k in range(8):
            nr, nc = r + int(DROW[k]), c + int(DCOL[k])
            if nr < 0 or nr >= rows or nc < 0 or nc >= cols:
                continue
            if labels[nr, nc] == labels[r, c] and away[nr, nc] == 0:
                edges.append((nr, nc))
    return away, flat_height


def _combined_gradient(labels, direction, low_edges, away, flat_height):
    """BFS inward from low edges, superimposed on the away-gradient.
    RichDEM BuildTowardsCombinedGradient:
        flat_mask starts as -away; on visit,
        flat_mask = (flat_height[label] + flat_mask) + 2*loops   (if !=0)
                  = 2*loops                                       (if ==0)
    """
    rows, cols = labels.shape
    flat_mask = -away.astype(np.int64)
    edges = deque((r, c) for (r, c) in low_edges if labels[r, c] != 0)
    MARK = (-1, -1); edges.append(MARK); loops = 1
    while len(edges) > 1:
        r, c = edges.popleft()
        if r == -1:
            loops += 1; edges.append(MARK); continue
        if flat_mask[r, c] > 0:
            continue
        if flat_mask[r, c] != 0:
            flat_mask[r, c] = (flat_height[labels[r, c]] + flat_mask[r, c]) + 2 * loops
        else:
            flat_mask[r, c] = 2 * loops
        for k in range(8):
            nr, nc = r + int(DROW[k]), c + int(DCOL[k])
            if nr < 0 or nr >= rows or nc < 0 or nc >= cols:
                continue
            if labels[nr, nc] == labels[r, c] and flat_mask[nr, nc] <= 0:
                edges.append((nr, nc))
    return flat_mask


def _assign_flat_dirs(elev, valid, labels, direction, flat_mask):
    """Assign flow directions to flat cells.

    Primary rule (RichDEM d8_masked_FlowDir): flow to the same-label
    neighbour of lowest flat_mask, cardinal over diagonal on ties.

    Boundary rule: a flat cell that is the local flat_mask minimum has no
    lower same-label neighbour - it is a low-edge cell whose true receiver
    lies OUTSIDE the flat (the spill point). For such cells, route to the
    lowest neighbour that is either strictly lower or already draining,
    exactly as an ordinary D8 outlet would. Without this, the flat's
    minimum becomes a false sink and the whole flat piles into it.
    """
    rows, cols = labels.shape
    resolved = 0
    for fid in np.flatnonzero(labels.reshape(-1) > 0):
        r, c = divmod(int(fid), cols)
        e = elev[r, c]
        best_val = flat_mask[r, c]
        best_k = -1
        for k in range(8):
            nr, nc = r + int(DROW[k]), c + int(DCOL[k])
            if nr < 0 or nr >= rows or nc < 0 or nc >= cols:
                continue
            if labels[nr, nc] != labels[r, c]:
                continue
            v = flat_mask[nr, nc]
            if v < best_val or (v == best_val and best_k >= 0
                                and (not _IS_CARDINAL[best_k]) and _IS_CARDINAL[k]):
                best_val = v
                best_k = k
        if best_k >= 0:
            direction[r, c] = best_k
            resolved += 1
            continue
        # Local minimum: seek an exit out of the flat (lower or draining).
        exit_k = -1
        for k in range(8):
            nr, nc = r + int(DROW[k]), c + int(DCOL[k])
            if nr < 0 or nr >= rows or nc < 0 or nc >= cols:
                continue
            if not valid[nr, nc]:
                continue
            lower = elev[nr, nc] < e - 1e-9
            drains_equal = (direction[nr, nc] >= 0 and abs(elev[nr, nc] - e) <= 1e-9)
            if lower or drains_equal:
                if exit_k < 0 or (not _IS_CARDINAL[exit_k]) and _IS_CARDINAL[k]:
                    exit_k = k
        if exit_k >= 0:
            direction[r, c] = exit_k
            resolved += 1
    return resolved


def _one_pass(elev, valid, direction):
    low, high = _find_flat_edges(elev, valid, direction)
    if not low:
        return 0, 0
    labels, n = _label_flats(elev, valid, direction, low)
    away, flat_height = _away_gradient(labels, direction, high, n)
    flat_mask = _combined_gradient(labels, direction, low, away, flat_height)
    resolved = _assign_flat_dirs(elev, valid, labels, direction, flat_mask)
    return resolved, n


def resolve_flats_barnes(elev, valid, direction, max_iterations=50):
    """Iterated Barnes flat resolution. Modifies `direction` in place.

    Each iteration is a full Barnes pass. Because a pass only resolves
    flats that have an outlet by the current direction state, repeating it
    lets flats drain in hierarchy order (lowest first). Iteration stops
    when a pass resolves nothing more.
    """
    elev = np.asarray(elev, dtype=np.float64)
    total_resolved = 0
    passes = 0
    regions_first = 0
    for it in range(max_iterations):
        resolved, n = _one_pass(elev, valid, direction)
        if it == 0:
            regions_first = n
        passes += 1
        total_resolved += resolved
        if resolved == 0:
            break

    outletless = int((valid & (direction < 0)).sum())
    return direction, {
        "flat_cells_resolved": total_resolved,
        "flat_regions_first_pass": regions_first,
        "iterations": passes,
        "outletless_flat_cells": outletless,
        "method": "Barnes 2014 (iterated)",
    }
