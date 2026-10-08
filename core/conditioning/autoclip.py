# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Clip a large DEM to the road's contributing area (F13, v0.25).

QEHT holds the whole grid in memory; a national DEM does not fit on a work
laptop, and clipping by hand risks cutting catchments. The clip is found
on a coarse copy of the DEM:

1. Downsample by k (default: the smallest k that gives at most 4 Mcells),
   taking the MINIMUM of each block so valleys are kept.
2. Fill, D8 (Barnes) and receivers on the coarse grid.
3. Mark the coarse cells within the road buffer (default 200 m) and
   propagate "drains to the road" upstream: every cell whose receiver chain
   reaches a marked cell (pointer jumping, as hand_grid).
4. Dilate by the margin (default 1 km plus 2 coarse cells); take the
   bounding box (or the mask).
5. Read that window of the full-resolution DEM (gdal.Translate srcWin, in
   the processing layer) into rasters/dem_clip.tif, tagged QEHT_AUTOCLIP.
6. After the full-resolution run, a catchment that touches the clip edge
   (the grid edge, or NoData outside the mask) gets clip_edge = 1.

No QGIS imports; GDAL imports are function-local.
"""

import json
import math

import numpy as np

from ..grid import receivers_from_direction

TAG = "QEHT_AUTOCLIP"


def coarse_factor(rows, cols, max_mcells=4.0):
    return max(1, int(math.ceil(math.sqrt(rows * cols / (max_mcells * 1e6)))))


def block_min(z, valid, k):
    """Minimum of each k x k block (NaN where a block has no valid cell)."""
    rows, cols = z.shape
    R, Cc = -(-rows // k), -(-cols // k)
    pad = np.full((R * k, Cc * k), np.inf)
    pad[:rows, :cols] = np.where(valid, z, np.inf)
    m = pad.reshape(R, k, Cc, k).min(axis=(1, 3))
    return np.where(np.isfinite(m), m, np.nan)


def drains_to(marked, direction, valid):
    """Bool grid: cells whose D8 receiver chain reaches a marked cell (pointer jumping)."""
    shape = direction.shape
    rec = receivers_from_direction(np.where(valid, direction, -1), shape)
    idx = np.arange(rec.size, dtype=np.int64)
    mk = np.asarray(marked, bool).ravel() & np.asarray(valid, bool).ravel()
    ptr = np.where(mk | (rec < 0), idx, rec)
    for _ in range(64):
        nxt = ptr[ptr]
        if np.array_equal(nxt, ptr):
            break
        ptr = nxt
    return (mk[ptr] & np.asarray(valid, bool).ravel()).reshape(shape)


def mark_road(alignment, gt, shape, buffer_m):
    """Cells of the grid (geotransform gt) within buffer_m of the alignment (dense sampling)."""
    rows, cols = shape
    cell = min(abs(gt[1]), abs(gt[5]))
    step = cell / 3.0
    m = np.zeros(shape, bool)
    offs = np.arange(-buffer_m, buffer_m + 1e-9, step) if buffer_m > 0 else np.array([0.0])
    offs = np.union1d(offs, [0.0])
    for k in range(len(alignment.x0)):
        n = max(int(math.ceil(alignment.seg_len[k] / step)), 1)
        u = np.linspace(0.0, 1.0, n + 1)
        x = alignment.x0[k] + u * alignment.dx[k]
        y = alignment.y0[k] + u * alignment.dy[k]
        tx, ty = alignment.tangent(k)
        X = (x[:, None] - offs[None, :] * ty).ravel()
        Y = (y[:, None] + offs[None, :] * tx).ravel()
        c = np.floor((X - gt[0]) / gt[1]).astype(np.int64)
        r = np.floor((Y - gt[3]) / gt[5]).astype(np.int64)
        ok = (r >= 0) & (r < rows) & (c >= 0) & (c < cols)
        m[r[ok], c[ok]] = True
    return m


def dilate_cells(mask, n):
    """Square dilation by n cells (separable: rows, then columns)."""
    out = np.asarray(mask, bool).copy()
    if n <= 0:
        return out
    acc = out.copy()
    for s in range(1, n + 1):
        acc[s:, :] |= out[:-s, :]
        acc[:-s, :] |= out[s:, :]
    out = acc.copy()
    for s in range(1, n + 1):
        out[:, s:] |= acc[:, :-s]
        out[:, :-s] |= acc[:, s:]
    return out


def contributing_window(z, valid, gt, alignment, road_buffer_m=200.0, margin_m=1000.0,
                        max_mcells=4.0, k=None, mode="bbox", pad_cells=2):
    """The clip on a coarse copy of the full-resolution grid.

    z: full-resolution DEM (float, NaN or invalid where NoData), valid: bool.
    Returns dict: window (row0, row1, col0, col1) in full-resolution cells
    (end exclusive), mask (full-resolution bool inside the window, for
    mode 'mask'), k, coarse shape, cells and memory estimate.
    """
    from .fill import fill_depressions
    from ..flow.direction import d8_direction
    rows, cols = z.shape
    k = k or coarse_factor(rows, cols, max_mcells)
    zc = block_min(np.asarray(z, float), valid, k)
    vc = np.isfinite(zc)
    cgt = (gt[0], gt[1] * k, 0.0, gt[3], 0.0, gt[5] * k)
    cw, ch = abs(cgt[1]), abs(cgt[5])
    filled, _, _ = fill_depressions(np.nan_to_num(zc, nan=0.0), vc, cell_width=cw, cell_height=ch)
    d, _ = d8_direction(filled, vc, cw, ch)
    road = mark_road(alignment, cgt, zc.shape, road_buffer_m) & vc
    up = drains_to(road, d, vc) | road
    n = int(math.ceil(margin_m / min(cw, ch))) + int(pad_cells)
    keep = dilate_cells(up, n) & vc
    if not keep.any():
        raise ValueError("The road does not cross the DEM; nothing to clip.")
    rr, cc = np.nonzero(keep)
    r0, r1 = int(rr.min()) * k, min((int(rr.max()) + 1) * k, rows)
    c0, c1 = int(cc.min()) * k, min((int(cc.max()) + 1) * k, cols)
    out = {"window": (r0, r1, c0, c1), "k": k, "coarse_shape": list(zc.shape),
           "coarse_cells_kept": int(keep.sum()), "mode": mode,
           "full_cells": rows * cols, "clip_cells": (r1 - r0) * (c1 - c0)}
    if mode == "mask":
        fm = np.repeat(np.repeat(keep, k, axis=0), k, axis=1)[:rows, :cols]
        out["mask"] = fm[r0:r1, c0:c1]
        out["clip_cells_valid"] = int(out["mask"].sum())
    return out


def memory_note(full_cells, clip_cells):
    from ..memory import estimate_peak_memory
    a = estimate_peak_memory(1, full_cells)["peak_gb"]
    b = estimate_peak_memory(1, clip_cells)["peak_gb"]
    return {"full_peak_gb": a, "clip_peak_gb": b}


def write_clip(src_path, out_path, info, window, mask=None, meta=None):
    """Window of the full-resolution DEM -> GeoTIFF (gdal.Translate srcWin), NoData outside
    the mask in 'mask' mode; tagged QEHT_AUTOCLIP."""
    from osgeo import gdal
    gdal.UseExceptions()
    r0, r1, c0, c1 = window
    ds = gdal.Translate(out_path, src_path, options=gdal.TranslateOptions(
        srcWin=[c0, r0, c1 - c0, r1 - r0], outputType=gdal.GDT_Float32, noData=-9999.0,
        creationOptions=["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"]))
    b = ds.GetRasterBand(1)
    a = b.ReadAsArray().astype(np.float32)
    nd = info.nodata
    if nd is not None:
        a[a == nd] = -9999.0
    if mask is not None:
        a[~mask] = -9999.0
    b.WriteArray(a)
    b.SetNoDataValue(-9999.0)
    md = dict(ds.GetMetadata() or {})
    md[TAG] = json.dumps(dict(meta or {}, window=list(window)))
    ds.SetMetadata(md)
    ds.FlushCache()
    ds = None
    return out_path


def touches_edge(mask, valid):
    """1 when a catchment mask touches the grid edge or a NoData cell (8-neighbourhood)."""
    m = np.asarray(mask, bool)
    if not m.any():
        return 0
    if m[0].any() or m[-1].any() or m[:, 0].any() or m[:, -1].any():
        return 1
    inv = ~np.asarray(valid, bool)
    near = np.zeros_like(m)
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            near |= np.roll(np.roll(inv, dr, 0), dc, 1)
    return int((m & near).any())
