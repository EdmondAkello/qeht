# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Vectorise a raster stream network into polyline reaches.

A reach runs from a source or a junction downstream to the next junction
or to the network outlet. That is the "Drainage Line" concept used by
commercial hydrology toolsets
and the unit an engineer actually works with - one feature per channel
segment, carrying its own order, length and slope.

This is deliberately NOT a cell-by-cell polygonize of the stream raster,
which yields one tangled multipart feature, nor one feature per cell,
which yields hundreds of thousands of two-point lines.
"""

import numpy as np
from collections import deque

from ..grid import DROW, DCOL, NO_RECEIVER, neighbour_distances, receivers_from_direction


def vectorize_streams(direction, valid, stream_mask, cell_width=1.0,
                      cell_height=1.0, elevation=None, order=None,
                      progress=None):
    """Trace the stream raster into reaches.

    Returns a list of dicts, one per reach:
        cells      ordered [(row, col), ...] upstream -> downstream
        length     planimetric length in map units (true diagonal weighting)
        order      Strahler order of the reach (if `order` supplied)
        drop       elevation drop (if `elevation` supplied)
        slope      drop / length
        head_type  'source' or 'junction'
        reach_id   sequential integer
    """
    rows, cols = direction.shape
    n = rows * cols
    dist = neighbour_distances(cell_width, cell_height)

    stream = np.asarray(stream_mask, dtype=bool) & valid
    stream_flat = stream.reshape(-1)
    dir_flat = direction.reshape(-1)

    receiver = receivers_from_direction(direction, (rows, cols))
    # Restrict the graph to the stream network.
    has = receiver != NO_RECEIVER
    keep = has.copy()
    keep[has] = stream_flat[receiver[has]]
    receiver = receiver.copy()
    receiver[~keep] = NO_RECEIVER
    receiver[~stream_flat] = NO_RECEIVER

    # Count upstream stream contributors per cell.
    upstream_count = np.zeros(n, dtype=np.int32)
    live = stream_flat & (receiver != NO_RECEIVER)
    np.add.at(upstream_count, receiver[live], 1)

    # A reach starts at a source (no upstream stream cell) or immediately
    # below a junction (a cell whose receiver has >= 2 contributors).
    is_source = stream_flat & (upstream_count == 0)
    is_junction = stream_flat & (upstream_count >= 2)

    starts = []
    for flat_id in np.flatnonzero(is_source):
        starts.append((int(flat_id), "source"))
    # Cells draining into a junction start a new reach at the junction.
    for flat_id in np.flatnonzero(is_junction):
        starts.append((int(flat_id), "junction"))

    reaches = []
    total = max(1, len(starts))

    for idx, (start, head_type) in enumerate(starts):
        if progress is not None and idx % 500 == 0:
            progress(idx / total, "Tracing stream reaches")

        cells = [start]
        cur = start
        length = 0.0

        while True:
            nxt = receiver[cur]
            if nxt == NO_RECEIVER:
                break
            length += dist[int(dir_flat[cur])]
            cells.append(int(nxt))
            cur = int(nxt)
            # Stop when we arrive at a junction: the next reach starts there.
            if upstream_count[cur] >= 2:
                break

        if len(cells) < 2:
            continue

        reach = {
            "reach_id": len(reaches) + 1,
            "cells": [(int(c) // cols, int(c) % cols) for c in cells],
            "length": float(length),
            "head_type": head_type,
        }
        if order is not None:
            reach["order"] = int(order.reshape(-1)[cells[0]])
        if elevation is not None:
            elev_flat = np.asarray(elevation).reshape(-1)
            top = float(elev_flat[cells[0]])
            bottom = float(elev_flat[cells[-1]])
            reach["drop"] = top - bottom
            reach["slope"] = (top - bottom) / length if length > 0 else 0.0
        reaches.append(reach)

    if progress is not None:
        progress(1.0, "Tracing stream reaches")
    return reaches


def network_summary(reaches, catchment_area=None):
    """Morphometric summary of a vectorised network."""
    if not reaches:
        return {}
    lengths = np.array([r["length"] for r in reaches])
    total = float(lengths.sum())
    summary = {
        "reaches": len(reaches),
        "total channel length": total,
        "mean reach length": float(lengths.mean()),
        "longest reach": float(lengths.max()),
    }
    orders = [r.get("order") for r in reaches if r.get("order")]
    if orders:
        summary["maximum Strahler order"] = int(max(orders))
    if catchment_area:
        summary["drainage density (1/map unit)"] = total / catchment_area
    return summary
