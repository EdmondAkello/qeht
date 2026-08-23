# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Catchment and longest-flow-path characteristics.

These are the numbers that feed a design flood calculation, so each one
is defined explicitly here rather than left to the reader's assumption.

On "catchment slope"
--------------------
Two definitions are in common use and they are NOT interchangeable:

  * **Mean terrain slope** - the average gradient of the ground surface
    over every cell in the catchment, computed by Horn's 3x3 method (a
    standard slope algorithm used by most GIS raster toolsets). Runoff coefficient and
    curve number tables generally assume this one.

  * **Relief ratio** - (highest elevation - lowest elevation) divided by
    the longest flow path length. This is the "catchment slope" of most
    road drainage manuals and is the term that goes into Kirpich.

Both are returned, under distinct names. Do not substitute one for the
other in a design calculation without saying which you used.

On flow path slope
------------------
The plain slope is drop over planimetric length, divided over the whole
path. The 10-85 slope excludes the top 10% and bottom 15% of the path
length, which removes the steep headwater and the flat outlet reach that
otherwise distort the average. It is required by several UK and TRRL
methods and is usually the more defensible figure for a design flood.
"""

import numpy as np

from ..grid import DROW, DCOL, neighbour_distances


def horn_slope(elevation, valid, cell_width, cell_height):
    """Terrain slope per cell, in m/m, by Horn's 3x3 method.

    This is a standard slope algorithm used across GIS raster toolsets, so
    values are directly comparable to a slope raster produced by one of them.
    """
    z = np.asarray(elevation, dtype=np.float64)
    rows, cols = z.shape

    # Pad by edge replication so border cells still get a value.
    p = np.pad(np.where(valid, z, np.nan), 1, mode="edge")

    a = p[0:-2, 0:-2]; b = p[0:-2, 1:-1]; c = p[0:-2, 2:]
    d = p[1:-1, 0:-2];                    f = p[1:-1, 2:]
    g = p[2:,   0:-2]; h = p[2:,   1:-1]; i = p[2:,   2:]

    with np.errstate(invalid="ignore"):
        dzdx = ((c + 2.0 * f + i) - (a + 2.0 * d + g)) / (8.0 * cell_width)
        dzdy = ((g + 2.0 * h + i) - (a + 2.0 * b + c)) / (8.0 * cell_height)
        slope = np.hypot(dzdx, dzdy)

    slope[~valid] = np.nan
    return slope


def path_slope_10_85(cells, elevation, cell_width, cell_height, direction=None):
    """Slope between the 10% and 85% points along a flow path.

    `cells` runs from the catchment divide to the outlet.
    """
    if len(cells) < 3:
        return float("nan")

    dist = neighbour_distances(cell_width, cell_height)
    cumulative = [0.0]
    for (r0, c0), (r1, c1) in zip(cells[:-1], cells[1:]):
        step = np.hypot((c1 - c0) * cell_width, (r1 - r0) * cell_height)
        cumulative.append(cumulative[-1] + step)
    total = cumulative[-1]
    if total <= 0:
        return float("nan")

    cumulative = np.array(cumulative)
    lo_target, hi_target = 0.10 * total, 0.85 * total
    lo_idx = int(np.searchsorted(cumulative, lo_target))
    hi_idx = int(np.searchsorted(cumulative, hi_target))
    lo_idx = min(lo_idx, len(cells) - 1)
    hi_idx = min(hi_idx, len(cells) - 1)
    if hi_idx <= lo_idx:
        return float("nan")

    z_lo = float(elevation[cells[lo_idx]])
    z_hi = float(elevation[cells[hi_idx]])
    span = cumulative[hi_idx] - cumulative[lo_idx]
    if span <= 0:
        return float("nan")
    # cells run divide -> outlet, so the upstream point is higher
    return (z_lo - z_hi) / span


def catchment_characteristics(mask, elevation, valid, cell_width, cell_height,
                              flow_path=None, slope_raster=None):
    """Morphometry for one catchment.

    Parameters
    ----------
    mask : bool array marking the catchment.
    elevation : the DEM to measure. Use the RAW DEM, not the filled one,
        so reported elevations are real ground rather than fill surface.
    flow_path : optional dict from longest_flow_path(), used for the
        relief ratio and the path statistics.
    slope_raster : optional precomputed Horn slope, to avoid recomputing
        it once per catchment on a multi-outlet run.

    Returns a dict of named characteristics in metres, km and m/m.
    """
    cell_area = abs(cell_width * cell_height)
    live = mask & valid
    n = int(live.sum())
    if n == 0:
        return {}

    z = elevation[live]
    z = z[np.isfinite(z)]
    if z.size == 0:
        return {}

    z_max = float(z.max())
    z_min = float(z.min())
    relief = z_max - z_min

    if slope_raster is None:
        slope_raster = horn_slope(elevation, valid, cell_width, cell_height)
    cell_slopes = slope_raster[live]
    cell_slopes = cell_slopes[np.isfinite(cell_slopes)]

    out = {
        "area_km2": n * cell_area / 1.0e6,
        "area_m2": n * cell_area,
        "cells": n,
        "elev_max_m": z_max,
        "elev_min_m": z_min,
        "elev_mean_m": float(z.mean()),
        "relief_m": relief,
        "slope_mean": float(cell_slopes.mean()) if cell_slopes.size else float("nan"),
        "slope_median": float(np.median(cell_slopes)) if cell_slopes.size else float("nan"),
    }

    if flow_path and flow_path.get("cells"):
        cells = flow_path["cells"]
        length = float(flow_path.get("length", 0.0))
        z_path = np.array([elevation[rc] for rc in cells], dtype=np.float64)
        z_path = z_path[np.isfinite(z_path)]

        lfp_max = float(z_path.max()) if z_path.size else float("nan")
        lfp_min = float(z_path.min()) if z_path.size else float("nan")
        lfp_drop = lfp_max - lfp_min

        out.update({
            "lfp_length_km": length / 1000.0,
            "lfp_length_m": length,
            "lfp_elev_max_m": lfp_max,
            "lfp_elev_min_m": lfp_min,
            "lfp_drop_m": lfp_drop,
            "lfp_slope": lfp_drop / length if length > 0 else float("nan"),
            "lfp_slope_1085": path_slope_10_85(cells, elevation,
                                               cell_width, cell_height),
            # Relief ratio: the "catchment slope" of most road drainage
            # manuals. Distinct from slope_mean above.
            "slope_relief_ratio": relief / length if length > 0 else float("nan"),
        })

    return out


# Field order and types for vector output. Kept here so the raster core
# and the Processing layer cannot drift apart.
CATCHMENT_FIELDS = [
    ("outlet_id", "int"),
    ("area_km2", "float"),
    ("elev_max_m", "float"),
    ("elev_min_m", "float"),
    ("elev_mean_m", "float"),
    ("relief_m", "float"),
    ("slope_mean", "float"),
    ("slope_relief_ratio", "float"),
    ("lfp_length_km", "float"),
]

FLOWPATH_FIELDS = [
    ("outlet_id", "int"),
    ("lfp_length_km", "float"),
    ("lfp_elev_max_m", "float"),
    ("lfp_elev_min_m", "float"),
    ("lfp_drop_m", "float"),
    ("lfp_slope", "float"),
    ("lfp_slope_1085", "float"),
    ("area_km2", "float"),
]
