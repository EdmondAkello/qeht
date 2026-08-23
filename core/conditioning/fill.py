# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Depression filling by priority-flood.

Deliberately separated from D8 routing. DDM HydroLogic fuses the two -
its priority-flood pass emits the receiver graph directly, so the flow
direction it produces is the traversal parent rather than the steepest
descent neighbour. That is fine for connectivity tracing but it is not
compatible with standard D8 routing and it ignores diagonal distance entirely.

Here, conditioning produces ONLY a conditioned elevation surface. Routing
is a separate decision made in `qeht.core.flow.direction`. That separation
is what lets a user document which conditioning method was applied - which
is the difference between a defensible drainage design and a black box.

Reference: Barnes, Lehman & Mulla (2014), "Priority-flood: an optimal
depression-filling and watershed-labeling algorithm for digital elevation
models", Computers & Geosciences 62:117-127. O(n log n).
"""

import heapq

import numpy as np

from ..grid import DROW, DCOL


def valid_boundary_mask(valid):
    """Valid cells that touch the grid edge or touch a NoData cell.

    A clipped DEM usually arrives with a NoData collar, so the outer row
    and column are not the hydrological boundary. Seeding only the raster
    edge would force the whole grid to drain to one global low point.
    """
    rows, cols = valid.shape
    boundary = np.zeros_like(valid, dtype=bool)

    boundary[0, :] |= valid[0, :]
    boundary[-1, :] |= valid[-1, :]
    boundary[:, 0] |= valid[:, 0]
    boundary[:, -1] |= valid[:, -1]

    # Valid cells orthogonally or diagonally adjacent to NoData.
    invalid = ~valid
    padded = np.zeros((rows + 2, cols + 2), dtype=bool)
    padded[1:-1, 1:-1] = invalid
    touches_nodata = np.zeros_like(valid, dtype=bool)
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue
            touches_nodata |= padded[1 + dr:rows + 1 + dr, 1 + dc:cols + 1 + dc]
    boundary |= valid & touches_nodata
    return boundary


def fill_depressions(dem, valid, min_slope=0.0, cell_width=1.0, cell_height=1.0,
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


def depression_depth(dem, filled, valid):
    """Per-cell fill depth - a QA/QC product, not an intermediate.

    Large contiguous depths in a satellite DEM usually mark a road
    embankment, a dam wall, or a culvert the DEM does not know about.
    Reviewing this raster before accepting the conditioning is the step
    that catches a catchment draining the wrong side of a road.
    """
    depth = np.full(dem.shape, np.nan, dtype=np.float64)
    depth[valid] = filled[valid] - dem[valid]
    depth[valid & (depth < 0)] = 0.0
    return depth
