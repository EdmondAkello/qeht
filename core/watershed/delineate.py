# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Stream extraction, pour-point snapping, catchment delineation and
longest flow path.

Longest flow path is one of the harder tools to find a good open equivalent for,
and the one that matters most for road drainage: it is the length term in
almost every time-of-concentration formula you will use.
"""

import numpy as np
from collections import deque

from ..grid import DROW, DCOL, NO_RECEIVER, neighbour_distances, receivers_from_direction


def extract_streams(accumulation, valid, threshold_cells=None,
                    threshold_area=None, cell_area=None):
    """Boolean stream mask from an accumulation threshold.

    Give either a cell count or an area (in map units squared); area is
    the more defensible parameter in an engineering report because it
    survives a change of DEM resolution.
    """
    if threshold_area is not None:
        if cell_area is None:
            raise ValueError("cell_area is required when using threshold_area")
        threshold_cells = float(threshold_area) / float(cell_area)
    if threshold_cells is None:
        raise ValueError("Provide threshold_cells or threshold_area")

    mask = np.zeros(accumulation.shape, dtype=bool)
    with np.errstate(invalid="ignore"):
        mask[valid] = accumulation[valid] >= float(threshold_cells)
    return mask


def snap_pour_point(row, col, accumulation, valid, search_radius_cells=5,
                    stream_mask=None):
    """Move a pour point onto the drainage network.

    Two strategies. Use the stream one.

    `stream_mask` given  -> NEAREST cell on the extracted stream network,
        ties at equal distance broken by higher accumulation. This
        preserves which tributary the point sits on.

    `stream_mask` None   -> highest accumulation within the radius. Simple,
        and WRONG near confluences: validated against a reference hydrology
        toolset on the
        Site A DEM, this pulled 5 of 16 road-crossing points off their
        small tributary and onto the adjacent trunk stream, inflating one
        catchment by 40x. Kept only for backward compatibility. Do not use
        it for crossing points spaced closer than the snap radius.

    Returns (row, col, moved_cells, accumulation_at_result).
    """
    rows, cols = accumulation.shape

    if stream_mask is not None:
        # Expand ring by ring so the nearest stream cell wins on distance
        # first, accumulation only as a tie-break within the same ring.
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
        # No stream cell within the radius: leave the point untouched
        # rather than silently grabbing whatever is nearby.
        return row, col, 0, float(accumulation[row, col])

    r0 = max(0, row - search_radius_cells)
    r1 = min(rows, row + search_radius_cells + 1)
    c0 = max(0, col - search_radius_cells)
    c1 = min(cols, col + search_radius_cells + 1)

    window = accumulation[r0:r1, c0:c1].copy()
    window_valid = valid[r0:r1, c0:c1]
    window[~window_valid] = -np.inf
    window = np.nan_to_num(window, nan=-np.inf)

    if not np.isfinite(window).any():
        return row, col, 0, float("nan")

    local = np.unravel_index(np.argmax(window), window.shape)
    best_r, best_c = r0 + int(local[0]), c0 + int(local[1])
    moved = int(round(np.hypot(best_r - row, best_c - col)))
    return best_r, best_c, moved, float(accumulation[best_r, best_c])


def delineate_catchment(direction, valid, outlet_rc):
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


def longest_flow_paths_by_catchment(direction, valid, outlets, elevation=None,
                                    cell_width=1.0, cell_height=1.0,
                                    nested=False, progress=None):
    """Longest flow path for EVERY outlet, each bounded by its own catchment.

    v0.1/v0.2 delineated one catchment at a time and wrote only the first
    path, which is why only a single line ever appeared in the output.
    This computes one path per outlet, with each path confined to that
    outlet's own contributing area.

    `nested=False` (default): each outlet gets its FULL upstream area, so
        catchments of outlets on the same channel overlap. This is what
        you want for a road crossing - the design flow at a culvert comes
        from everything upstream of it.

    `nested=True`: outlets are processed in the order given and each one
        only claims cells not already taken, producing non-overlapping
        local catchments. Supply outlets upstream-first. This reproduces
        the "Catchment" (local drainage area) concept used by commercial
        hydrology toolsets.

    Returns a list of dicts, one per outlet, each with the keys returned
    by longest_flow_path plus 'outlet_index', 'catchment_cells' and
    'catchment_area'.
    """
    rows, cols = direction.shape
    cell_area = abs(cell_width * cell_height)
    results = []

    if nested:
        labels = delineate_catchment(direction, valid, list(outlets))

    for i, (orow, ocol) in enumerate(outlets):
        if progress is not None:
            progress(i / max(1, len(outlets)), "Tracing longest flow paths")

        if nested:
            mask = labels == (i + 1)
        else:
            mask = delineate_catchment(direction, valid, [(orow, ocol)]) > 0

        n_cells = int(mask.sum())
        if n_cells == 0:
            results.append({"outlet_index": i, "cells": [], "length": 0.0,
                            "n_cells": 0, "catchment_cells": 0,
                            "catchment_area": 0.0})
            continue

        lfp = longest_flow_path(direction, valid, (orow, ocol),
                                elevation=elevation, cell_width=cell_width,
                                cell_height=cell_height, catchment_mask=mask)
        lfp["outlet_index"] = i
        lfp["catchment_cells"] = n_cells
        lfp["catchment_area"] = n_cells * cell_area
        results.append(lfp)

    if progress is not None:
        progress(1.0, "Tracing longest flow paths")
    return results


def longest_flow_path(direction, valid, outlet_rc, elevation=None,
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
        catchment_mask = delineate_catchment(direction, valid, [outlet_rc]) > 0
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
