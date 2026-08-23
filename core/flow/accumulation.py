# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Flow accumulation and Strahler ordering by topological traversal.

O(N) after flow direction is established. No recursion - a recursive
downstream walk overflows the stack on any real DEM, and Python's
recursion limit would be hit long before the memory limit.

Algorithm: count how many upstream cells drain into each cell (indegree),
queue every cell with indegree 0 (a ridge/source cell), then repeatedly
pop a cell, push its load downstream, and decrement the receiver's
indegree - enqueueing it once it reaches 0. Each cell is popped exactly
once and its total is complete when popped, so there is no double
counting.
"""

import numpy as np
from collections import deque

from ..grid import NO_RECEIVER, receivers_from_direction


def flow_accumulation(direction, valid, weights=None, progress=None):
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


def strahler_order(direction, valid, stream_mask, progress=None):
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
