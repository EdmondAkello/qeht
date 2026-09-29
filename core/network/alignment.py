# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Linear referencing of a road alignment (pure NumPy).

Chainage = start chainage + planimetric distance along the alignment.
Multi-part alignments are chained in the order given; the gap between the
end of one part and the start of the next is NOT counted, so chainage runs
continuously along the carriageway parts. Offsets are signed: positive to
the LEFT of the direction of increasing chainage.

No QGIS imports, no GDAL.
"""

import numpy as np


class Alignment(object):
    """A road centreline as ordered straight segments with chainage.

        al = Alignment([[(x0, y0), (x1, y1), ...], ...], start_chainage=0.0)
        ch, off, seg = al.locate(xs, ys)
    """

    def __init__(self, parts, start_chainage=0.0, reverse=False):
        clean = []
        for p in parts:
            a = np.asarray(p, dtype=np.float64).reshape(-1, 2)
            # drop repeated vertices (zero-length segments)
            if len(a) > 1:
                keep = np.ones(len(a), dtype=bool)
                keep[1:] = np.any(np.diff(a, axis=0) != 0, axis=1)
                a = a[keep]
            if len(a) >= 2:
                clean.append(a)
        if not clean:
            raise ValueError("The alignment has no segment of non-zero length.")
        if reverse:
            clean = [a[::-1] for a in clean[::-1]]
        self.parts = clean
        x0, y0, x1, y1, part = [], [], [], [], []
        for k, a in enumerate(clean):
            x0.append(a[:-1, 0]); y0.append(a[:-1, 1])
            x1.append(a[1:, 0]); y1.append(a[1:, 1])
            part.append(np.full(len(a) - 1, k))
        self.x0 = np.concatenate(x0); self.y0 = np.concatenate(y0)
        self.x1 = np.concatenate(x1); self.y1 = np.concatenate(y1)
        self.part = np.concatenate(part)
        self.dx = self.x1 - self.x0
        self.dy = self.y1 - self.y0
        self.seg_len = np.hypot(self.dx, self.dy)
        self.ch0 = float(start_chainage) + np.concatenate(
            [[0.0], np.cumsum(self.seg_len)[:-1]])
        self.start_chainage = float(start_chainage)
        self.length = float(self.seg_len.sum())
        # last segment of each part: its end point belongs to it (u == 1)
        self.part_end = np.zeros(len(self.x0), dtype=bool)
        self.part_end[np.r_[np.flatnonzero(np.diff(self.part)), len(self.x0) - 1]] = True

    @property
    def end_chainage(self):
        return self.start_chainage + self.length

    def tangent(self, seg):
        """Unit direction of segment(s) `seg`."""
        seg = np.asarray(seg)
        return self.dx[seg] / self.seg_len[seg], self.dy[seg] / self.seg_len[seg]

    def locate(self, x, y, chunk=4096):
        """Nearest-point projection of points onto the alignment.

        Returns (chainage, signed_offset, segment_index) arrays; offset is
        the perpendicular distance, positive left of the alignment direction.
        """
        x = np.atleast_1d(np.asarray(x, dtype=np.float64))
        y = np.atleast_1d(np.asarray(y, dtype=np.float64))
        n = x.size
        ch = np.empty(n); off = np.empty(n); seg = np.empty(n, dtype=np.int64)
        L2 = self.seg_len ** 2
        for s in range(0, n, chunk):
            px = x[s:s + chunk, None]; py = y[s:s + chunk, None]
            t = ((px - self.x0) * self.dx + (py - self.y0) * self.dy) / L2
            t = np.clip(t, 0.0, 1.0)
            qx = self.x0 + t * self.dx; qy = self.y0 + t * self.dy
            d2 = (px - qx) ** 2 + (py - qy) ** 2
            k = np.argmin(d2, axis=1)
            rows = np.arange(k.size)
            tk = t[rows, k]
            ch[s:s + chunk] = self.ch0[k] + tk * self.seg_len[k]
            cross = (self.dx[k] * (py[:, 0] - self.y0[k])
                     - self.dy[k] * (px[:, 0] - self.x0[k]))
            off[s:s + chunk] = np.where(cross >= 0, 1.0, -1.0) * np.sqrt(d2[rows, k])
            seg[s:s + chunk] = k
        return ch, off, seg

    def point_at(self, chainage):
        """(x, y) at a chainage (clamped to the alignment)."""
        c = float(chainage)
        k = int(np.clip(np.searchsorted(self.ch0, c, side="right") - 1, 0, len(self.ch0) - 1))
        u = np.clip((c - self.ch0[k]) / self.seg_len[k], 0.0, 1.0)
        return float(self.x0[k] + u * self.dx[k]), float(self.y0[k] + u * self.dy[k])

    def intersect(self, ax, ay, bx, by, chunk=2048):
        """Intersections of line segments a->b with the alignment.

        A crossing is counted when the segment parameter t is in [0, 1) (so a
        path of consecutive links crossing exactly at a shared vertex counts
        once) and the alignment parameter u is in [0, 1) (u == 1 allowed only
        at the end of a part). Collinear overlaps are ignored.

        Returns arrays (link_index, seg_index, t, u, x, y).
        """
        ax = np.asarray(ax, dtype=np.float64); ay = np.asarray(ay, dtype=np.float64)
        bx = np.asarray(bx, dtype=np.float64); by = np.asarray(by, dtype=np.float64)
        out = [[] for _ in range(6)]
        for s in range(0, ax.size, chunk):
            rx = (bx[s:s + chunk] - ax[s:s + chunk])[:, None]
            ry = (by[s:s + chunk] - ay[s:s + chunk])[:, None]
            qpx = self.x0[None, :] - ax[s:s + chunk, None]
            qpy = self.y0[None, :] - ay[s:s + chunk, None]
            denom = rx * self.dy[None, :] - ry * self.dx[None, :]
            with np.errstate(divide="ignore", invalid="ignore"):
                t = (qpx * self.dy[None, :] - qpy * self.dx[None, :]) / denom
                u = (qpx * ry - qpy * rx) / denom
            eps = 1e-12
            u_hi = np.where(self.part_end[None, :], u <= 1.0 + eps, u < 1.0 - eps)
            hit = (denom != 0) & (t >= -eps) & (t < 1.0 - eps) & (u >= -eps) & u_hi
            li, sj = np.nonzero(hit)
            tt = t[li, sj]; uu = np.clip(u[li, sj], 0.0, 1.0)
            out[0].append(li + s); out[1].append(sj); out[2].append(tt); out[3].append(uu)
            out[4].append(ax[li + s] + tt * rx[li, 0]); out[5].append(ay[li + s] + tt * ry[li, 0])
        return tuple(np.concatenate(o) if o else np.array([]) for o in out)
