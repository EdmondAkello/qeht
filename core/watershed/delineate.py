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
        # Nearest stream cell by ring (Chebyshev distance); within the
        # nearest ring the highest accumulation wins, and on equal
        # accumulation the first cell in row-major order (the v0.8.3
        # ring-by-ring scan, vectorised over the window in v0.12).
        rad = int(search_radius_cells)
        r0, r1 = max(0, row - rad), min(rows, row + rad + 1)
        c0, c1 = max(0, col - rad), min(cols, col + rad + 1)
        cand = stream_mask[r0:r1, c0:c1] & valid[r0:r1, c0:c1]
        if cand.any():
            rr, cc = np.nonzero(cand)
            gr, gc = rr + r0, cc + c0
            ring = np.maximum(np.abs(gr - row), np.abs(gc - col))
            acc = accumulation[gr, gc].astype(np.float64)
            order = np.lexsort((np.arange(rr.size), -acc, ring))   # ring, then -acc, then scan
            k = order[0]
            br, bc = int(gr[k]), int(gc[k])
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

    Reverse traversal over the receiver graph from each outlet, one whole
    frontier at a time (vectorised in v0.12; same result as the v0.8.3
    per-cell stack, kept in core/flow/_reference.py). O(N).
    """
    rows, cols = direction.shape
    n = rows * cols
    valid_flat = valid.reshape(-1)
    receiver = receivers_from_direction(direction, (rows, cols))

    live = valid_flat & (receiver != NO_RECEIVER)
    src = np.flatnonzero(live)
    dst = receiver[live]
    counts = np.bincount(dst, minlength=n)
    starts = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(counts, out=starts[1:])
    children = src[np.argsort(dst, kind="stable")]

    labels = np.zeros(n, dtype=np.int32)
    if isinstance(outlet_rc, tuple) and len(outlet_rc) == 2 and np.isscalar(outlet_rc[0]):
        outlet_rc = [outlet_rc]

    for idx, (orow, ocol) in enumerate(outlet_rc, start=1):
        seed = int(orow) * cols + int(ocol)
        if not valid_flat[seed] or labels[seed] != 0:
            continue
        labels[seed] = idx
        front = np.array([seed], dtype=np.int64)
        while front.size:
            cnt = counts[front]
            tot = int(cnt.sum())
            if tot == 0:
                break
            offs = np.repeat(starts[front] - (np.cumsum(cnt) - cnt), cnt) + np.arange(tot)
            ch = children[offs]
            ch = ch[labels[ch] == 0]
            labels[ch] = idx
            front = ch

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

    v0.12 processes the topological order in waves (vectorised). Where two
    branches are equally long, the branch chosen is the one the v0.8.3
    first-in-first-out queue met first - reproduced exactly by tracking
    each cell's wave and position - so paths, lengths and every slope
    derived from them are unchanged.

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

    front = np.flatnonzero(active & (indegree == 0))        # FIFO order = this order
    while front.size:
        j = live_recv[front]
        okj = j != NO_RECEIVER
        src, dst = front[okj], j[okj]
        pos = np.flatnonzero(okj)                           # FIFO position in this wave
        if src.size:
            cand = upstream_len[src] + dist[dir_flat[src]]
            # best child per receiver this wave: max length, then earliest position
            o = np.lexsort((pos, -cand, dst))
            d_s, c_s, s_s = dst[o], cand[o], src[o]
            first = np.ones(d_s.size, dtype=bool)
            first[1:] = d_s[1:] != d_s[:-1]
            d1, c1, s1 = d_s[first], c_s[first], s_s[first]
            better = c1 > upstream_len[d1]                  # strict: earlier waves keep ties
            upstream_len[d1[better]] = c1[better]
            from_cell[d1[better]] = s1[better]
            np.subtract.at(indegree, dst, 1)
            # next wave in FIFO order: by the position of each cell's LAST child
            ready = indegree[dst] == 0
            rd, rp = dst[ready], pos[ready]
            o = np.lexsort((-rp, rd))
            last = np.ones(rd.size, dtype=bool)
            last[1:] = rd[o][1:] != rd[o][:-1]
            rd, rp = rd[o][last], rp[o][last]
            front = rd[np.argsort(rp, kind="stable")]
        else:
            front = front[:0]

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
