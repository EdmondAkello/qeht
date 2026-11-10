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
                            tolerance=0.10, search_cells=1, stream_threshold_cells=None,
                            flat_shares=True):
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
    stream_threshold_cells : stream threshold of the run (cells). With it,
        an outlet is "moved" when one method has a stream within the window
        and the other has none: the channel itself runs elsewhere, so the
        crossing point, not only its area, depends on the flat method
        (1.1). Without it, "moved" means the smaller area is under 10 % of
        the larger.
    flat_shares : also report how much of each catchment and of its longest
        flow path (Barnes routing) lies on resolved flats (1.1).

    Returns a list of dicts, one per outlet: area_barnes_km2,
    area_toward_km2, flat_sensitivity_pct (|difference| / larger area,
    in %), flat_sensitive (1 when above the tolerance, else 0),
    flat_area_ratio (toward / Barnes), flat_moved (0/1) and, with
    flat_shares, flat_catch_share_pct and flat_lfp_share_pct (cells with
    no downhill neighbour after filling, whose direction the flat method
    chose).
    """
    cell_km2 = abs(cell_width * cell_height) / 1e6
    k = int(max(search_cells, 0))
    areas, cells_max, where, dirs = {}, {}, {}, {}
    for method in ("barnes", "toward"):
        direction, _ = d8_direction(conditioned, valid, cell_width, cell_height,
                                    resolve_flats=True, flat_method=method)
        acc, _ = flow_accumulation(direction, valid)
        dirs[method] = direction
        vals, mx, pos = [], [], []
        for r, c in outlets_rc:
            if not valid[r, c]:
                vals.append(float("nan")); mx.append(float("nan")); pos.append(None)
                continue
            r0, c0 = max(r - k, 0), max(c - k, 0)
            win = np.where(valid[r0:r + k + 1, c0:c + k + 1], acc[r0:r + k + 1, c0:c + k + 1], -1.0)
            i = np.unravel_index(int(np.argmax(win)), win.shape)
            m = float(win[i])
            vals.append((m + 1.0) * cell_km2)
            mx.append(m + 1.0)
            pos.append((r0 + int(i[0]), c0 + int(i[1])))
        areas[method], cells_max[method], where[method] = vals, mx, pos
    flats = None
    if flat_shares:
        from ..conditioning.fill import valid_boundary_mask
        base, _ = d8_direction(conditioned, valid, cell_width, cell_height, resolve_flats=False)
        flats = valid & (base < 0) & ~valid_boundary_mask(valid)
    out = []
    for j, (ab, at) in enumerate(zip(areas["barnes"], areas["toward"])):
        big = max(ab, at)
        rel = abs(ab - at) / big if big > 0 else 0.0
        nb, nt = cells_max["barnes"][j], cells_max["toward"][j]
        if not (np.isfinite(nb) and np.isfinite(nt)):
            moved = None
        elif stream_threshold_cells:
            moved = int((nb >= stream_threshold_cells) != (nt >= stream_threshold_cells))
        else:
            moved = int(min(nb, nt) < 0.1 * max(nb, nt))
        rec = {"area_barnes_km2": ab, "area_toward_km2": at,
               "flat_sensitivity_pct": 100.0 * rel,
               "flat_sensitive": int(rel > tolerance),
               "flat_area_ratio": (at / ab) if ab and np.isfinite(ab) and ab > 0 else None,
               "flat_moved": moved}
        if flats is not None:
            rec["flat_catch_share_pct"] = rec["flat_lfp_share_pct"] = None
            o = where["barnes"][j]
            if o is not None:
                from ..watershed.delineate import delineate_catchment, longest_flow_path
                mask = delineate_catchment(dirs["barnes"], valid, [o]) == 1
                n = int(mask.sum())
                if n:
                    rec["flat_catch_share_pct"] = 100.0 * float((flats & mask).sum()) / n
                    lfp = longest_flow_path(dirs["barnes"], valid, o, cell_width=cell_width,
                                            cell_height=cell_height, catchment_mask=mask)
                    cells = lfp.get("cells") or []
                    if cells:
                        rr, cc = zip(*cells)
                        rec["flat_lfp_share_pct"] = 100.0 * float(flats[list(rr), list(cc)].mean())
        out.append(rec)
    return out
