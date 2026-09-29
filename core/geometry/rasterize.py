# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Polygons -> integer label grid, pure NumPy (cell-centre rule).

A cell takes a polygon's value when its CENTRE lies inside the polygon
(even-odd rule, so holes work). This is the default rule of
gdal.RasterizeLayer (without ALL_TOUCHED), and the test suite checks the
two agree. Later polygons overwrite earlier ones where they overlap.

No QGIS imports, no GDAL.
"""

import numpy as np


def rasterize_polygons(polygons, geotransform, shape, fill=0, dtype=np.int32):
    """Burn polygons into a grid.

    polygons : list of (rings, value); rings = [outer, hole, ...], each a
        sequence of (x, y) (closed or not). Multipolygons: pass each part
        as its own entry with the same value.
    geotransform : GDAL-style, north-up. shape : (rows, cols).
    Returns the label grid.
    """
    gt = tuple(geotransform)
    rows, cols = shape
    out = np.full((rows, cols), fill, dtype=dtype)
    yc = gt[3] + (np.arange(rows) + 0.5) * gt[5]          # row centre y
    xc0 = gt[0] + 0.5 * gt[1]
    for rings, value in polygons:
        segs = []
        for ring in rings:
            a = np.asarray(ring, dtype=np.float64).reshape(-1, 2)
            if len(a) < 3:
                continue
            if not np.array_equal(a[0], a[-1]):
                a = np.vstack([a, a[:1]])
            segs.append(np.hstack([a[:-1], a[1:]]))
        if not segs:
            continue
        s = np.vstack(segs)                                # x0 y0 x1 y1
        ymin, ymax = s[:, [1, 3]].min(), s[:, [1, 3]].max()
        r_idx = np.flatnonzero((yc >= ymin) & (yc <= ymax))
        if r_idx.size == 0:
            continue
        # edges straddling each scanline (half-open in y: no double count at vertices)
        y0, y1 = s[:, 1][None, :], s[:, 3][None, :]
        yy = yc[r_idx][:, None]
        hit = (y0 <= yy) != (y1 <= yy)
        with np.errstate(divide="ignore", invalid="ignore"):
            xint = s[:, 0][None, :] + (yy - y0) * (s[:, 2] - s[:, 0])[None, :] / (y1 - y0)
        for k, r in enumerate(r_idx):
            xs = np.sort(xint[k][hit[k]])
            for a, b in zip(xs[0::2], xs[1::2]):
                # cells whose centre x is in [a, b)
                c0 = int(np.ceil((a - xc0) / gt[1] - 1e-12))
                c1 = int(np.ceil((b - xc0) / gt[1] - 1e-12)) - 1
                c0, c1 = max(c0, 0), min(c1, cols - 1)
                if c0 <= c1:
                    out[r, c0:c1 + 1] = value
    return out
