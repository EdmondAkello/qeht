# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Flat-method sensitivity of catchment areas (v0.13.1).

Where a DEM has flats or exactly tied cells on a drainage line, the flat
resolution method - not the terrain - decides which way the flow goes. On
a real road project one tied cell on a river moved ~176 km2 between two
culverts depending on the method. This check routes the conditioned DEM
with both methods (Barnes 2014 and toward lower terrain) and reports the
contributing area at each outlet cell under each, so crossings whose area
depends on an arbitrary choice are flagged for verification against
mapped drainage or a site visit.
"""

import numpy as np

from .direction import d8_direction
from .accumulation import flow_accumulation


def flat_method_sensitivity(conditioned, valid, cell_width, cell_height, outlets_rc,
                            tolerance=0.10, search_cells=1):
    """Contributing area at each outlet under both flat methods.

    conditioned : filled (conditioned) DEM used for routing
    outlets_rc : list of (row, col) outlet cells (e.g. the snapped outlets)
    tolerance : relative difference above which an outlet is flagged
    search_cells : the area under each method is the largest within this
        many cells of the outlet (default 1, a 3x3 window). The two methods
        can route the same drainage line one cell apart on a flat; reading
        the outlet cell alone then reports a few cells instead of the
        catchment (seen at every crossing of a 30 m FABDEM corridor, 0.15).
        A real diversion to another crossing still shows, because the flow
        leaves the window altogether.

    Returns a list of dicts, one per outlet: area_barnes_km2,
    area_toward_km2, flat_sensitivity_pct (|difference| / larger area,
    in %), flat_sensitive (1 when above the tolerance, else 0).
    """
    cell_km2 = abs(cell_width * cell_height) / 1e6
    areas = {}
    for method in ("barnes", "toward"):
        direction, _ = d8_direction(conditioned, valid, cell_width, cell_height,
                                    resolve_flats=True, flat_method=method)
        acc, _ = flow_accumulation(direction, valid)
        rows, cols = acc.shape
        k = int(max(search_cells, 0))
        vals = []
        for r, c in outlets_rc:
            if not valid[r, c]:
                vals.append(float("nan"))
                continue
            win = acc[max(r - k, 0):r + k + 1, max(c - k, 0):c + k + 1]
            wv = valid[max(r - k, 0):r + k + 1, max(c - k, 0):c + k + 1]
            vals.append((float(win[wv].max()) + 1.0) * cell_km2)
        areas[method] = vals
    out = []
    for ab, at in zip(areas["barnes"], areas["toward"]):
        big = max(ab, at)
        rel = abs(ab - at) / big if big > 0 else 0.0
        out.append({"area_barnes_km2": ab, "area_toward_km2": at,
                    "flat_sensitivity_pct": 100.0 * rel,
                    "flat_sensitive": int(rel > tolerance)})
    return out
