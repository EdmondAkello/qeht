# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Fill-depression review (1.2).

Filling raises every closed depression to its spill level. Most of these
are DEM noise, but a deep fill against a road embankment usually marks a
place where the DEM treats the road as a dam: there may be a culvert the
DEM cannot see, or a missing one. This module turns the fill-depth grid
into one row per depression so the deep ones near the road can be listed,
mapped and checked.

A depression is an 8-connected group of cells raised by more than
`extent_tol` (default 1 cm, so the tiny gradient a minimum-slope fill adds
across flats does not join separate depressions). For each one:

    cells, area_m2, max_depth_m, mean_depth_m, volume_m3 (sum of depth x
    cell area), x / y of the deepest cell, and, with a road alignment,
    dist_road_m (nearest depression cell to the centreline), chainage_m
    (of the deepest cell), touches_road (some cell within the road
    half-width), nearest crossing and its chainage distance, and a review
    flag.

Flags (defaults are review settings, not design criteria - change them to
suit the DEM and the road):

    deep_at_road   touches the road and max depth >= deep_m
    deep_near_road within near_m of the centreline and max depth >= deep_m
    ""             anything else

`no_crossing_within_m` on a flagged row is the chainage distance to the
nearest crossing; a deep fill against the road with no crossing close by
is the case most worth a look.

No QGIS imports, no GDAL.
"""

import numpy as np

from ..flow.flats import _components
from ..geometry.polygonize import mask_to_polygons

FLAG_AT_ROAD = "deep_at_road"
FLAG_NEAR_ROAD = "deep_near_road"

DEPRESSION_FIELDS = [
    ("dep_id", "int"), ("cells", "int"), ("area_m2", "float"), ("max_depth_m", "float"),
    ("mean_depth_m", "float"), ("volume_m3", "float"), ("deep_x", "float"), ("deep_y", "float"),
    ("dist_road_m", "float"), ("chainage_m", "float"), ("touches_road", "int"),
    ("nearest_uid", "str"), ("nearest_ch_dist_m", "float"), ("flag", "str"), ("note", "str"),
]


def fill_depth(dem, filled, valid):
    depth = np.where(valid, filled.astype(np.float64) - dem.astype(np.float64), 0.0)
    depth[~np.isfinite(depth)] = 0.0
    return np.maximum(depth, 0.0)


def label_depressions(depth, extent_tol=0.01):
    """-> (labels int32 grid, 0 = none, 1..n), n."""
    mask = depth > extent_tol
    rows, cols = mask.shape
    idx = np.full(mask.shape, -1, dtype=np.int64)
    n = int(mask.sum())
    labels = np.zeros(mask.shape, dtype=np.int32)
    if n == 0:
        return labels, 0
    idx[mask] = np.arange(n)
    a, b = [], []
    for dr, dc in ((0, 1), (1, -1), (1, 0), (1, 1)):          # each undirected link once
        r0, r1 = max(0, -dr), rows - max(0, dr)
        c0, c1 = max(0, -dc), cols - max(0, dc)
        src = idx[r0:r1, c0:c1]
        dst = idx[r0 + dr:r1 + dr, c0 + dc:c1 + dc]
        ok = (src >= 0) & (dst >= 0)
        a.append(src[ok]); b.append(dst[ok])
    parent = _components(n, np.concatenate(a), np.concatenate(b))
    _, lab = np.unique(parent, return_inverse=True)
    labels[mask] = (lab + 1).astype(np.int32)
    return labels, int(lab.max()) + 1


def _boundary(mask):
    """Cells of mask with a 4-neighbour outside it (the mask's edge cells)."""
    p = np.pad(mask, 1, constant_values=False)
    inner = p[:-2, 1:-1] & p[2:, 1:-1] & p[1:-1, :-2] & p[1:-1, 2:]
    return mask & ~inner


def depression_table(dem, filled, valid, geotransform, alignment=None, crossings=None,
                     min_depth=0.5, deep_m=1.0, road_half_width_m=15.0, near_m=100.0,
                     extent_tol=0.01, max_rows=5000, with_polygons=True):
    """-> (rows, depth grid, summary).

    rows: one dict per depression with max depth >= min_depth (deepest
    first, at most max_rows), fields as DEPRESSION_FIELDS, plus '_geom'
    (list of polygons) when with_polygons.
    alignment: core.network.alignment.Alignment or None.
    crossings: list of (outlet_uid, chainage_m) or None.
    """
    gt = geotransform
    cw, ch = abs(gt[1]), abs(gt[5])
    cell_area = cw * ch
    depth = fill_depth(dem, filled, valid)
    labels, n = label_depressions(depth, extent_tol)
    summary = {"depressions": n, "reported": 0, "raised_cells": int((labels > 0).sum()),
               "volume_m3": float(depth[labels > 0].sum() * cell_area),
               "flag_at_road": 0, "flag_near_road": 0}
    if n == 0:
        return [], depth, summary
    flat = labels.ravel()
    d = depth.ravel()
    sel = flat > 0
    L = flat[sel]
    cells = np.bincount(L, minlength=n + 1)
    vol = np.bincount(L, weights=d[sel], minlength=n + 1)
    mx = np.zeros(n + 1)
    np.maximum.at(mx, L, d[sel])
    keep = np.flatnonzero(mx >= min_depth)
    keep = keep[keep > 0]
    keep = keep[np.argsort(-mx[keep], kind="stable")][:max_rows]
    summary["reported"] = int(keep.size)
    if keep.size == 0:
        return [], depth, summary
    # deepest cell of each kept depression
    lin = np.flatnonzero(sel)
    order = np.lexsort((-d[lin], flat[lin]))
    first = lin[order][np.r_[0, np.flatnonzero(np.diff(flat[lin][order])) + 1]]
    deep_of = {int(flat[i]): int(i) for i in first}
    # bounding boxes
    rr, cc = np.divmod(lin, labels.shape[1])
    rmin = np.full(n + 1, labels.shape[0]); rmax = np.full(n + 1, -1)
    cmin = np.full(n + 1, labels.shape[1]); cmax = np.full(n + 1, -1)
    np.minimum.at(rmin, L, rr); np.maximum.at(rmax, L, rr)
    np.minimum.at(cmin, L, cc); np.maximum.at(cmax, L, cc)

    xing = sorted(crossings or [], key=lambda t: t[1])
    xch = np.array([c for _, c in xing], dtype=np.float64)
    rows = []
    for k, lab in enumerate(keep, start=1):
        lab = int(lab)
        r0, r1, c0, c1 = int(rmin[lab]), int(rmax[lab]) + 1, int(cmin[lab]), int(cmax[lab]) + 1
        m = labels[r0:r1, c0:c1] == lab
        i = deep_of[lab]
        dr, dc = divmod(i, labels.shape[1])
        x = gt[0] + (dc + 0.5) * gt[1] + (dr + 0.5) * gt[2]
        y = gt[3] + (dc + 0.5) * gt[4] + (dr + 0.5) * gt[5]
        row = {"dep_id": k, "cells": int(cells[lab]), "area_m2": float(cells[lab] * cell_area),
               "max_depth_m": float(mx[lab]), "mean_depth_m": float(vol[lab] / cells[lab]),
               "volume_m3": float(vol[lab] * cell_area), "deep_x": float(x), "deep_y": float(y),
               "dist_road_m": None, "chainage_m": None, "touches_road": None,
               "nearest_uid": None, "nearest_ch_dist_m": None, "flag": "", "note": ""}
        if alignment is not None:
            br, bc = np.nonzero(_boundary(m))
            if br.size > 4000:                          # enough to place a large depression
                step = int(np.ceil(br.size / 4000))
                br, bc = br[::step], bc[::step]
            bx = gt[0] + (bc + c0 + 0.5) * gt[1]
            by = gt[3] + (br + r0 + 0.5) * gt[5]
            _, off, _ = alignment.locate(np.r_[bx, x], np.r_[by, y])
            # a cell centre is up to half a diagonal from its edge
            dist = max(0.0, float(np.abs(off[:-1]).min()) - 0.5 * float(np.hypot(cw, ch)))
            chg, _, _ = alignment.locate([x], [y])
            row["dist_road_m"] = dist
            row["chainage_m"] = float(chg[0])
            row["touches_road"] = int(dist <= road_half_width_m)
            if xch.size:
                j = int(np.argmin(np.abs(xch - chg[0])))
                row["nearest_uid"] = str(xing[j][0])
                row["nearest_ch_dist_m"] = float(abs(xch[j] - chg[0]))
            if row["max_depth_m"] >= deep_m and row["touches_road"]:
                row["flag"] = FLAG_AT_ROAD
                summary["flag_at_road"] += 1
            elif row["max_depth_m"] >= deep_m and dist <= near_m:
                row["flag"] = FLAG_NEAR_ROAD
                summary["flag_near_road"] += 1
            if row["flag"]:
                row["note"] = _note(row, near_m)
        if with_polygons:
            sub_gt = (gt[0] + c0 * gt[1], gt[1], gt[2], gt[3] + r0 * gt[5], gt[4], gt[5])
            row["_geom"] = mask_to_polygons(m, sub_gt)
        rows.append(row)
    return rows, depth, summary


def _note(row, near_m):
    where = "against the road" if row["flag"] == FLAG_AT_ROAD else f"{row['dist_road_m']:.0f} m from the road"
    text = (f"Filled {row['max_depth_m']:.2f} m deep, {row['area_m2']:,.0f} m2, "
            f"{row['volume_m3']:,.0f} m3, {where}")
    if row["nearest_uid"] is None:
        return text + "; no crossing on the road. Check whether the road ponds water here."
    if row["nearest_ch_dist_m"] > near_m:
        return (text + f"; nearest crossing {row['nearest_uid']} is {row['nearest_ch_dist_m']:,.0f} m "
                "along the road. Check for a culvert the DEM does not show, or a missing one.")
    return (text + f"; crossing {row['nearest_uid']} is {row['nearest_ch_dist_m']:,.0f} m along "
            "the road. Check that its catchment drains this depression.")
