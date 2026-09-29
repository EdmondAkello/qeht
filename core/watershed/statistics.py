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
path. The 10-85 slope is measured between the points at 10% and 85% of
the path length MEASURED FROM THE OUTLET upstream (the conventional
definition), so it excludes the bottom 10% (flat outlet reach) and the
top 15% (steep headwater):

    S_10-85 = (z85 - z10) / (L85 - L10),  L10 = 0.10 L,  L85 = 0.85 L

Elevations are interpolated linearly along the path at exactly those
distances. QEHT <= 0.8.3 measured the percentages from the divide instead
(points at 90% and 15% from the outlet); see CHANGELOG 0.9.0.

Field names (v0.9, HEAS exchange schema qeht-heas-1)
---------------------------------------------------
Four slope domains, never merged:
    catch_slope_horn    mean Horn 3x3 terrain slope        (was slope_mean)
    catch_relief_ratio  relief / LFP length                 (was slope_relief_ratio)
    lfp_slope           drop / length over the whole LFP
    lfp_slope_1085      10-85 slope along the LFP (outlet-referenced)
The old names are still written as aliases for one release.
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


def _cumulative_from_divide(cells, cell_width, cell_height):
    """Cumulative planimetric distance along `cells` (divide -> outlet)."""
    cum = np.zeros(len(cells), dtype=np.float64)
    for k in range(1, len(cells)):
        (r0, c0), (r1, c1) = cells[k - 1], cells[k]
        cum[k] = cum[k - 1] + np.hypot((c1 - c0) * cell_width,
                                       (r1 - r0) * cell_height)
    return cum


def path_10_85(cells, elevation, cell_width, cell_height, reference="outlet"):
    """10-85 slope along a flow path, with the points used to compute it.

    `cells` runs from the catchment divide to the outlet (as returned by
    longest_flow_path). Distances are the true planimetric path distance
    (diagonal steps = cell diagonal). Elevations at the 10% and 85% points
    are interpolated linearly between cell centres, so the result does not
    depend on which cell happens to be "next".

    reference="outlet" (default, v0.9+): L10/L85 measured from the outlet
        upstream - the conventional definition.
    reference="divide": the mirrored convention used by QEHT <= 0.8.3
        (kept only so the change can be demonstrated and tested).

    Returns a dict with slope, L10_m, L85_m (distance FROM THE OUTLET of
    the two points, whatever the reference), z10_m, z85_m (elevation at the
    points labelled 10 and 85 under the chosen reference) and length_m.
    Values are NaN when the path is shorter than two cells.
    """
    nan = float("nan")
    out = {"slope": nan, "L10_m": nan, "L85_m": nan, "z10_m": nan,
           "z85_m": nan, "length_m": nan, "reference": reference}
    if reference not in ("outlet", "divide"):
        raise ValueError("reference must be 'outlet' or 'divide'")
    if len(cells) < 2:
        return out

    cum = _cumulative_from_divide(cells, cell_width, cell_height)
    total = float(cum[-1])
    out["length_m"] = total
    if total <= 0:
        return out

    z = np.array([elevation[rc] for rc in cells], dtype=np.float64)
    if reference == "outlet":
        pos10, pos85 = 0.90 * total, 0.15 * total   # measured from divide
        L10, L85 = 0.10 * total, 0.85 * total       # measured from outlet
    else:
        pos10, pos85 = 0.10 * total, 0.85 * total
        L10, L85 = 0.90 * total, 0.15 * total
    z10 = float(np.interp(pos10, cum, z))
    z85 = float(np.interp(pos85, cum, z))
    span = 0.75 * total
    if reference == "outlet":
        slope = (z85 - z10) / span          # 85% point is upstream
    else:
        slope = (z10 - z85) / span          # 10% point is upstream
    out.update({"slope": slope, "L10_m": L10, "L85_m": L85,
                "z10_m": z10, "z85_m": z85})
    return out


def path_slope_10_85(cells, elevation, cell_width, cell_height,
                     direction=None, reference="outlet"):
    """10-85 slope only (see path_10_85). Outlet-referenced from v0.9."""
    return path_10_85(cells, elevation, cell_width, cell_height,
                      reference=reference)["slope"]


def path_slope_10_85_v083(cells, elevation, cell_width, cell_height):
    """The exact QEHT <= 0.8.3 computation, kept for regression comparison.

    Divide-referenced points, elevation taken at the first cell at or past
    each target distance (no interpolation). Do not use for design.
    """
    if len(cells) < 3:
        return float("nan")
    cumulative = _cumulative_from_divide(cells, cell_width, cell_height)
    total = cumulative[-1]
    if total <= 0:
        return float("nan")
    lo_idx = min(int(np.searchsorted(cumulative, 0.10 * total)), len(cells) - 1)
    hi_idx = min(int(np.searchsorted(cumulative, 0.85 * total)), len(cells) - 1)
    if hi_idx <= lo_idx:
        return float("nan")
    span = cumulative[hi_idx] - cumulative[lo_idx]
    if span <= 0:
        return float("nan")
    return (float(elevation[cells[lo_idx]]) - float(elevation[cells[hi_idx]])) / span


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
        "catch_slope_horn": float(cell_slopes.mean()) if cell_slopes.size else float("nan"),
        "slope_mean": float(cell_slopes.mean()) if cell_slopes.size else float("nan"),  # legacy alias
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

        s1085 = path_10_85(cells, elevation, cell_width, cell_height,
                           reference="outlet")
        relief_ratio = relief / length if length > 0 else float("nan")
        out.update({
            "lfp_length_km": length / 1000.0,
            "lfp_length_m": length,
            "lfp_elev_max_m": lfp_max,
            "lfp_elev_min_m": lfp_min,
            "lfp_drop_m": lfp_drop,
            "lfp_slope": lfp_drop / length if length > 0 else float("nan"),
            "lfp_slope_1085": s1085["slope"],
            "lfp_L10_m": s1085["L10_m"],
            "lfp_L85_m": s1085["L85_m"],
            "lfp_z10_m": s1085["z10_m"],
            "lfp_z85_m": s1085["z85_m"],
            # Relief ratio (relief / LFP length). Distinct from the Horn
            # mean slope; HEAS "catchment relief ratio" domain.
            "catch_relief_ratio": relief_ratio,
            "slope_relief_ratio": relief_ratio,     # legacy alias (<= 0.8.3)
        })

    return out


# Field order and types for the per-tool vector outputs. Kept here so the
# raster core and the Processing layer cannot drift apart. The HEAS
# exchange package has its own, fuller field list (core/interop).
# outlet_uid is the stable identifier (text); outlet_id is the pour-point
# feature id, kept for backward compatibility and NOT stable.
CATCHMENT_FIELDS = [
    ("outlet_uid", "text"),
    ("outlet_id", "int"),
    ("area_km2", "float"),
    ("elev_max_m", "float"),
    ("elev_min_m", "float"),
    ("elev_mean_m", "float"),
    ("relief_m", "float"),
    ("catch_slope_horn", "float"),
    ("catch_relief_ratio", "float"),
    ("lfp_length_km", "float"),
    ("slope_mean", "float"),            # legacy alias of catch_slope_horn
    ("slope_relief_ratio", "float"),    # legacy alias of catch_relief_ratio
]

FLOWPATH_FIELDS = [
    ("outlet_uid", "text"),
    ("outlet_id", "int"),
    ("lfp_length_km", "float"),
    ("lfp_elev_max_m", "float"),
    ("lfp_elev_min_m", "float"),
    ("lfp_drop_m", "float"),
    ("lfp_slope", "float"),
    ("lfp_slope_1085", "float"),
    ("lfp_L10_m", "float"),
    ("lfp_L85_m", "float"),
    ("lfp_z10_m", "float"),
    ("lfp_z85_m", "float"),
    ("area_km2", "float"),
]
