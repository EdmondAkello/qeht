# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""REFERENCE implementations (v0.8.3), kept verbatim as test oracles.

The production code (flats.py, direction.py) is vectorised from v0.12. The
per-cell Python loops below are the validated v0.8.3 algorithms; the test
suite requires the vectorised versions to reproduce them cell for cell.
Do not optimise or "fix" this file - its value is that it does not change.
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


def reference_resolve_flats_barnes(elev, valid, direction, max_iterations=50):
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


def reference_resolve_flats_toward(work, valid, direction, progress=None):
    """Route flat cells by BFS distance toward the nearest flat outlet."""
    rows, cols = work.shape
    flat = valid & (direction < 0)
    if not flat.any():
        return direction, 0

    INF = np.iinfo(np.int32).max
    bfs = np.full((rows, cols), INF, dtype=np.int32)
    queue = deque()

    # Seed: flat cells adjacent to an equal-elevation cell that already
    # has a direction (a genuine flat outlet).
    flat_ids = np.flatnonzero(flat.reshape(-1))
    for fid in flat_ids:
        r, c = divmod(int(fid), cols)
        for k in range(8):
            nr, nc = r + int(DROW[k]), c + int(DCOL[k])
            if nr < 0 or nr >= rows or nc < 0 or nc >= cols:
                continue
            if not valid[nr, nc]:
                continue
            if direction[nr, nc] >= 0 and work[nr, nc] <= work[r, c]:
                bfs[r, c] = 1
                queue.append((r, c))
                break

    while queue:
        r, c = queue.popleft()
        d = bfs[r, c]
        for k in range(8):
            nr, nc = r + int(DROW[k]), c + int(DCOL[k])
            if nr < 0 or nr >= rows or nc < 0 or nc >= cols:
                continue
            if not flat[nr, nc] or bfs[nr, nc] <= d + 1:
                continue
            if work[nr, nc] != work[r, c]:
                continue
            bfs[nr, nc] = d + 1
            queue.append((nr, nc))

    # Each flat cell drains to the neighbour with the smallest BFS
    # distance, preferring an already-routed cell at equal or lower
    # elevation (the flat outlet itself).
    resolved = 0
    for fid in flat_ids:
        r, c = divmod(int(fid), cols)
        best_k, best_d = -1, INF
        for k in range(8):
            nr, nc = r + int(DROW[k]), c + int(DCOL[k])
            if nr < 0 or nr >= rows or nc < 0 or nc >= cols:
                continue
            if not valid[nr, nc]:
                continue
            if work[nr, nc] > work[r, c]:
                continue
            nd = 0 if (direction[nr, nc] >= 0 and not flat[nr, nc]) else int(bfs[nr, nc])
            if nd < best_d:
                best_d, best_k = nd, k
        if best_k >= 0 and best_d < INF:
            direction[r, c] = best_k
            resolved += 1

    if progress is not None:
        progress(1.0, "Resolving flats")
    return direction, resolved


def reference_snap_to_stream(row, col, accumulation, valid, search_radius_cells, stream_mask):
    """v0.8.3 nearest-stream snapping (ring by ring), kept as an oracle."""
    rows, cols = accumulation.shape
    for radius in range(0, int(search_radius_cells) + 1):
        best = None
        for dr in range(-radius, radius + 1):
            for dc in range(-radius, radius + 1):
                if max(abs(dr), abs(dc)) != radius:
                    continue
                r, c = row + dr, col + dc
                if not (0 <= r < rows and 0 <= c < cols):
                    continue
                if not (stream_mask[r, c] and valid[r, c]):
                    continue
                acc = float(accumulation[r, c])
                if best is None or acc > best[0]:
                    best = (acc, r, c)
        if best is not None:
            _, br, bc = best
            moved = int(round(np.hypot(br - row, bc - col)))
            return br, bc, moved, float(accumulation[br, bc])
    return row, col, 0, float(accumulation[row, col])


_TIE_PRIORITY = [2, 4, 6, 0, 1, 3, 5, 7]
_TIE_RANK = np.array([_TIE_PRIORITY.index(k) for k in range(8)], dtype=np.int64)


def reference_d8_core(elev, valid, cell_width, cell_height, tie_tol=1e-6):
    """v0.8.3 D8 pass (rows x cols x 8 arrays), before flat resolution."""
    from ..grid import neighbour_distances, NO_RECEIVER
    elev = np.asarray(elev, dtype=np.float64)
    rows, cols = elev.shape
    dist = neighbour_distances(cell_width, cell_height)
    work = np.where(valid, elev, np.inf)
    all_slopes = np.full((rows, cols, 8), -np.inf, dtype=np.float64)
    for k in range(8):
        dr, dc = int(DROW[k]), int(DCOL[k])
        r0_s, r1_s = max(0, -dr), rows - max(0, dr)
        c0_s, c1_s = max(0, -dc), cols - max(0, dc)
        r0_n, r1_n = max(0, dr), rows - max(0, -dr)
        c0_n, c1_n = max(0, dc), cols - max(0, -dc)
        with np.errstate(invalid="ignore"):
            slope = (work[r0_s:r1_s, c0_s:c1_s] - work[r0_n:r1_n, c0_n:c1_n]) / dist[k]
        all_slopes[r0_s:r1_s, c0_s:c1_s, k] = np.where(np.isfinite(slope), slope, -np.inf)
    best_slope = all_slopes.max(axis=2)
    with np.errstate(invalid="ignore"):
        tied = all_slopes >= (best_slope[:, :, None] * (1.0 - tie_tol))
    ranked = np.where(tied, _TIE_RANK[None, None, :], 99)
    direction = np.argmin(ranked, axis=2).astype(np.int64)
    no_flow = ~(valid & np.isfinite(best_slope) & (best_slope > 0))
    direction[no_flow] = NO_RECEIVER
    direction[~valid] = NO_RECEIVER
    return direction


import heapq
from ..conditioning.fill import valid_boundary_mask
from ..grid import receivers_from_direction, NO_RECEIVER


def reference_fill_depressions(dem, valid, min_slope=0.0, cell_width=1.0, cell_height=1.0,
                     progress=None):
    """Fill depressions, returning (filled_dem, n_filled_cells, seed_count).

    Parameters
    ----------
    dem : 2-D float array, elevations. Values under ~valid are ignored.
    valid : 2-D bool array, True where the cell carries real elevation.
    min_slope : optional gradient (m per m) imposed across filled flats.
        0.0 gives a classic flat fill (matches common `Fill` tool behaviour). A small
        positive value (1e-4 is typical) tilts filled surfaces so that
        downstream routing has a defined direction instead of relying on
        flat-resolution. Applied as an epsilon increment per traversal
        step scaled by cell size.
    progress : optional callable(fraction_0_to_1, message).

    Returns
    -------
    filled : 2-D float array. Cells outside `valid` are set to NaN.
    n_filled : int, count of cells whose elevation was raised.
    seed_count : int, number of boundary outlet cells seeded.
    """
    dem = np.asarray(dem, dtype=np.float64)
    rows, cols = dem.shape
    filled = np.full((rows, cols), np.nan, dtype=np.float64)
    closed = np.zeros((rows, cols), dtype=bool)

    mean_cell = 0.5 * (abs(cell_width) + abs(cell_height))
    epsilon = float(min_slope) * mean_cell

    heap = []
    counter = 0
    seed_count = 0

    boundary = valid_boundary_mask(valid)
    seeds = np.flatnonzero(boundary.reshape(-1))
    for flat in seeds:
        r, c = divmod(int(flat), cols)
        closed[r, c] = True
        filled[r, c] = dem[r, c]
        heapq.heappush(heap, (float(dem[r, c]), counter, r, c))
        counter += 1
        seed_count += 1

    if not heap:
        # Degenerate case: fully interior valid region with no edge and no
        # NoData contact. Seed the single lowest valid cell.
        flat_ids = np.flatnonzero(valid.reshape(-1))
        if flat_ids.size == 0:
            return filled, 0, 0
        lowest = flat_ids[np.argmin(dem.reshape(-1)[flat_ids])]
        r, c = divmod(int(lowest), cols)
        closed[r, c] = True
        filled[r, c] = dem[r, c]
        heapq.heappush(heap, (float(dem[r, c]), counter, r, c))
        counter += 1
        seed_count = 1

    total = max(1, int(valid.sum()))
    processed = 0
    n_filled = 0

    while heap:
        elev, _seq, r, c = heapq.heappop(heap)
        processed += 1
        if progress is not None and processed % 20000 == 0:
            progress(processed / total, "Filling depressions")

        for k in range(8):
            nr = r + int(DROW[k])
            nc = c + int(DCOL[k])
            if nr < 0 or nr >= rows or nc < 0 or nc >= cols:
                continue
            if closed[nr, nc] or not valid[nr, nc]:
                continue

            closed[nr, nc] = True
            raw = float(dem[nr, nc])
            spill = elev + epsilon
            if raw < spill:
                new_elev = spill
                n_filled += 1
            else:
                new_elev = raw
            filled[nr, nc] = new_elev
            heapq.heappush(heap, (new_elev, counter, nr, nc))
            counter += 1

    if progress is not None:
        progress(1.0, "Filling depressions")
    return filled, n_filled, seed_count



def reference_flow_accumulation(direction, valid, weights=None, progress=None):
    """Accumulated upstream contribution per cell.

    Parameters
    ----------
    direction : 2-D internal direction index array (from d8_direction).
    valid : 2-D bool array.
    weights : optional 2-D float array of per-cell contribution. Default
        is 1.0 per cell, giving the standard D8 convention where accumulation is
        a count of upstream cells and EXCLUDES the cell itself. Pass a
        cell-area array for contributing area, or a rainfall/runoff grid
        for weighted accumulation.

    Returns
    -------
    accum : 2-D float array. NaN outside `valid`.
    stats : dict.
    """
    rows, cols = direction.shape
    n = rows * cols
    valid_flat = valid.reshape(-1)

    receiver = receivers_from_direction(direction, (rows, cols))
    # A receiver that is itself invalid is treated as off-grid.
    has_recv = receiver != NO_RECEIVER
    bad = has_recv.copy()
    bad[has_recv] = ~valid_flat[receiver[has_recv]]
    receiver[bad] = NO_RECEIVER

    if weights is None:
        weight_flat = np.ones(n, dtype=np.float64)
    else:
        weight_flat = np.asarray(weights, dtype=np.float64).reshape(-1).copy()
    weight_flat[~valid_flat] = 0.0

    # Standard D8 convention: a cell's accumulation excludes its own weight.
    total = np.zeros(n, dtype=np.float64)

    indegree = np.zeros(n, dtype=np.int32)
    live = valid_flat & (receiver != NO_RECEIVER)
    np.add.at(indegree, receiver[live], 1)

    queue = deque(int(i) for i in np.flatnonzero(valid_flat & (indegree == 0)))

    processed = 0
    n_valid = max(1, int(valid_flat.sum()))
    while queue:
        i = queue.popleft()
        processed += 1
        if progress is not None and processed % 50000 == 0:
            progress(processed / n_valid, "Accumulating flow")

        j = receiver[i]
        if j != NO_RECEIVER:
            total[j] += total[i] + weight_flat[i]
            indegree[j] -= 1
            if indegree[j] == 0:
                queue.append(int(j))

    unresolved = int((indegree > 0).sum())

    accum = np.full(n, np.nan, dtype=np.float64)
    accum[valid_flat] = total[valid_flat]
    accum = accum.reshape(rows, cols)

    stats = {
        "cells_processed": processed,
        "cells_in_cycles": unresolved,   # >0 means the graph has a loop
        "max_accumulation": float(np.nanmax(accum)) if processed else 0.0,
    }
    if progress is not None:
        progress(1.0, "Accumulating flow")
    return accum, stats




from ..grid import neighbour_distances


def reference_delineate_catchment(direction, valid, outlet_rc):
    """Upstream contributing area for one or more outlets.

    Returns an int32 label grid: 0 = outside, 1..n = catchment index in
    the order outlets were supplied. Where catchments nest, the FIRST
    outlet in the list wins, so supply downstream outlets last if you
    want nested subcatchments carved out.

    Implemented as a reverse traversal over the receiver graph, seeded at
    the outlets. O(N), no recursion.
    """
    rows, cols = direction.shape
    n = rows * cols
    valid_flat = valid.reshape(-1)
    receiver = receivers_from_direction(direction, (rows, cols))

    # Build upstream adjacency as a CSR-style structure - far cheaper in
    # memory than a dict of lists on a large DEM.
    live = valid_flat & (receiver != NO_RECEIVER)
    src = np.flatnonzero(live)
    dst = receiver[live]

    counts = np.bincount(dst, minlength=n)
    starts = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(counts, out=starts[1:])
    order = np.argsort(dst, kind="stable")
    children = src[order]

    labels = np.zeros(n, dtype=np.int32)
    if isinstance(outlet_rc, tuple) and len(outlet_rc) == 2 and np.isscalar(outlet_rc[0]):
        outlet_rc = [outlet_rc]

    for idx, (orow, ocol) in enumerate(outlet_rc, start=1):
        seed = int(orow) * cols + int(ocol)
        if not valid_flat[seed] or labels[seed] != 0:
            continue
        stack = [seed]
        labels[seed] = idx
        while stack:
            i = stack.pop()
            for p in range(starts[i], starts[i + 1]):
                child = int(children[p])
                if labels[child] == 0:
                    labels[child] = idx
                    stack.append(child)

    return labels.reshape(rows, cols)



def reference_longest_flow_path(direction, valid, outlet_rc, elevation=None,
                      cell_width=1.0, cell_height=1.0, catchment_mask=None):
    """Longest flow path from the catchment divide to the outlet.

    Returns a dict with:
        cells        list of (row, col) from the most remote cell to the outlet
        length       planimetric length in map units
        drop         elevation drop along the path (if `elevation` given)
        slope        drop / length (if `elevation` given)

    Method: propagate cumulative downstream distance upward through the
    receiver tree, take the maximum, then walk back down. Distances are
    accumulated per link using the true diagonal weighting, so the result
    is a real length rather than a cell count.

    The length and slope returned here are exactly the inputs required by
    Kirpich, Bransby-Williams and the TRRL time-of-concentration methods.
    """
    rows, cols = direction.shape
    n = rows * cols
    valid_flat = valid.reshape(-1)
    dist = neighbour_distances(cell_width, cell_height)

    if catchment_mask is None:
        catchment_mask = reference_delineate_catchment(direction, valid, [outlet_rc]) > 0
    in_catch = np.asarray(catchment_mask).reshape(-1)

    receiver = receivers_from_direction(direction, (rows, cols))
    dir_flat = direction.reshape(-1)

    # Topological order: sources first (indegree 0 within the catchment).
    active = in_catch & valid_flat
    live = active & (receiver != NO_RECEIVER)
    live_recv = receiver.copy()
    ok = live.copy()
    ok[live] = in_catch[receiver[live]]
    live_recv[~ok] = NO_RECEIVER

    indegree = np.zeros(n, dtype=np.int32)
    contributing = active & (live_recv != NO_RECEIVER)
    np.add.at(indegree, live_recv[contributing], 1)

    upstream_len = np.zeros(n, dtype=np.float64)
    from_cell = np.full(n, -1, dtype=np.int64)

    queue = deque(int(i) for i in np.flatnonzero(active & (indegree == 0)))
    while queue:
        i = queue.popleft()
        j = live_recv[i]
        if j == NO_RECEIVER:
            continue
        step = dist[int(dir_flat[i])]
        candidate = upstream_len[i] + step
        if candidate > upstream_len[j]:
            upstream_len[j] = candidate
            from_cell[j] = i
        indegree[j] -= 1
        if indegree[j] == 0:
            queue.append(int(j))

    outlet_flat = int(outlet_rc[0]) * cols + int(outlet_rc[1])
    path = [outlet_flat]
    cur = outlet_flat
    while from_cell[cur] >= 0:
        cur = int(from_cell[cur])
        path.append(cur)
    path.reverse()   # divide -> outlet

    cells = [(int(p) // cols, int(p) % cols) for p in path]
    result = {
        "cells": cells,
        "length": float(upstream_len[outlet_flat]),
        "n_cells": len(cells),
    }

    if elevation is not None and len(cells) >= 2:
        top = float(elevation[cells[0]])
        bottom = float(elevation[cells[-1]])
        drop = top - bottom
        result["drop"] = drop
        result["slope"] = drop / result["length"] if result["length"] > 0 else 0.0
    return result


def reference_strahler_order(direction, valid, stream_mask, progress=None):
    """Strahler order computed over an EXTRACTED STREAM NETWORK.

    This distinction is deliberate and important. Strahler order computed
    over every cell of a DEM is not the same quantity as Strahler order of
    the channel network, and the two are routinely confused. Orders here
    are defined only on cells where `stream_mask` is True; everything else
    is 0.

    Rule: a cell with no upstream stream cells is order 1. Otherwise, if
    the maximum upstream order occurs more than once, order = max + 1;
    if it occurs exactly once, order = max.
    """
    rows, cols = direction.shape
    n = rows * cols
    valid_flat = valid.reshape(-1)
    stream_flat = np.asarray(stream_mask).reshape(-1) & valid_flat

    receiver = receivers_from_direction(direction, (rows, cols))
    has_recv = receiver != NO_RECEIVER
    bad = has_recv.copy()
    bad[has_recv] = ~stream_flat[receiver[has_recv]]
    receiver = receiver.copy()
    receiver[bad] = NO_RECEIVER

    order = np.zeros(n, dtype=np.int32)
    max_up = np.zeros(n, dtype=np.int32)
    count_at_max = np.zeros(n, dtype=np.int32)

    indegree = np.zeros(n, dtype=np.int32)
    live = stream_flat & (receiver != NO_RECEIVER)
    np.add.at(indegree, receiver[live], 1)

    queue = deque(int(i) for i in np.flatnonzero(stream_flat & (indegree == 0)))

    while queue:
        i = queue.popleft()
        if max_up[i] == 0:
            order[i] = 1                       # headwater
        elif count_at_max[i] >= 2:
            order[i] = max_up[i] + 1           # confluence of equals
        else:
            order[i] = max_up[i]               # single dominant tributary

        j = receiver[i]
        if j != NO_RECEIVER:
            if order[i] > max_up[j]:
                max_up[j] = order[i]
                count_at_max[j] = 1
            elif order[i] == max_up[j]:
                count_at_max[j] += 1
            indegree[j] -= 1
            if indegree[j] == 0:
                queue.append(int(j))

    if progress is not None:
        progress(1.0, "Strahler ordering")
    return order.reshape(rows, cols)


# ---------------------------------------------------------------------------
# v0.13: per-cell port of Barnes' own implementation (RichDEM,
# include/richdem/flats/flat_resolution.hpp: find_flat_edges, label_this,
# BuildAwayGradient, BuildTowardsCombinedGradient, d8_masked_FlowDir), kept
# as the oracle for core/flow/flats.py at w = 2. RichDEM neighbour n = 1..8
# is W, NW, N, NE, E, SE, S, SW; QEHT index for each:
_RD_N = (4, 5, 6, 7, 0, 1, 2, 3)


def reference_barnes_richdem(elev, valid, direction, drain=None):
    """Barnes 2014 exactly as RichDEM runs it, on QEHT arrays.

    direction : QEHT indices, -1 = no flow. drain : cells without descent
    that drain off the grid (RichDEM gives edge cells an off-grid flow
    direction, so they count as routed). Modifies `direction` in place.
    """
    from collections import deque
    from ..grid import DROW, DCOL
    rows, cols = elev.shape
    drain = np.zeros(elev.shape, bool) if drain is None else drain
    noflow = valid & (direction < 0) & ~drain
    routed = valid & ~noflow

    def inside(r, c):
        return 0 <= r < rows and 0 <= c < cols

    low, high = deque(), deque()
    for c in range(cols):                    # RichDEM loops x (col) then y (row)
        for r in range(rows):
            if not valid[r, c]:
                continue
            for k in _RD_N:
                nr, nc = r + int(DROW[k]), c + int(DCOL[k])
                if not inside(nr, nc) or not valid[nr, nc]:
                    continue
                if routed[r, c] and noflow[nr, nc] and elev[nr, nc] == elev[r, c]:
                    low.append((r, c)); break
                elif noflow[r, c] and elev[r, c] < elev[nr, nc]:
                    high.append((r, c)); break
    if not low:
        return direction
    labels = np.zeros(elev.shape, np.int64)
    group = 1
    for (r0, c0) in low:
        if labels[r0, c0]:
            continue
        q = deque([(r0, c0)]); target = elev[r0, c0]
        while q:
            r, c = q.popleft()
            if not valid[r, c] or elev[r, c] != target or labels[r, c] > 0:
                continue
            labels[r, c] = group
            for k in _RD_N:
                nr, nc = r + int(DROW[k]), c + int(DCOL[k])
                if inside(nr, nc):
                    q.append((nr, nc))
        group += 1
    high = deque(h for h in high if labels[h])
    mask = np.zeros(elev.shape, np.int64)
    fh = np.zeros(group, np.int64)

    def bfs(edges, combine):
        loops = 1
        q = deque(edges); q.append(None)
        while len(q) != 1:
            cell = q.popleft()
            if cell is None:
                loops += 1; q.append(None); continue
            r, c = cell
            if mask[r, c] > 0:
                continue
            combine(r, c, loops)
            for k in _RD_N:
                nr, nc = r + int(DROW[k]), c + int(DCOL[k])
                if inside(nr, nc) and labels[nr, nc] == labels[r, c] and noflow[nr, nc]:
                    q.append((nr, nc))

    def away(r, c, loops):
        mask[r, c] = loops; fh[labels[r, c]] = loops

    def toward(r, c, loops):
        mask[r, c] = (fh[labels[r, c]] + mask[r, c]) + 2 * loops if mask[r, c] != 0 else 2 * loops

    bfs(high, away)
    mask *= -1
    bfs(low, toward)
    for r in range(rows):
        for c in range(cols):
            if not noflow[r, c] or labels[r, c] == 0:
                continue
            best, fd, fd_rd = mask[r, c], -1, 0
            for n, k in enumerate(_RD_N, start=1):
                nr, nc = r + int(DROW[k]), c + int(DCOL[k])
                if not inside(nr, nc) or labels[nr, nc] != labels[r, c]:
                    continue
                v = mask[nr, nc]
                if v < best or (v == best and fd_rd > 0 and fd_rd % 2 == 0 and n % 2 == 1):
                    best, fd, fd_rd = v, k, n
            if fd >= 0:
                direction[r, c] = fd
    return direction
