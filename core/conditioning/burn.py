# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Burn crossings through road embankments (WP-D).

A DSM or a LiDAR/photogrammetric DEM shows a road embankment as a dam:
the culvert under it is invisible. Filling then ponds the valley behind the
embankment and routes the flow over the crest wherever the crest happens to
be lowest - often a long way from the real culvert - which corrupts the
catchment boundary, the longest flow path and every slope along it.

At each accepted crossing this module cuts a short, straight breach across
the road: from the lowest cell near a point `half_length_m` upstream of the
centreline to the lowest cell near a point the same distance downstream,
perpendicular to the alignment. Along the breach the ground is lowered to a
straight grade between the two end elevations (never raised: z = min(z,
grade)). Nothing else in the DEM changes. Every breach is reported (cells,
maximum and total cut) so it can be recorded and checked.

Run it on the RAW DEM before filling. No QGIS imports, no GDAL.
"""

import math

import numpy as np


def _line_cells(r0, c0, r1, c1):
    """8-connected cells from (r0, c0) to (r1, c1) inclusive (Bresenham)."""
    cells = []
    dr, dc = abs(r1 - r0), abs(c1 - c0)
    sr = 1 if r1 >= r0 else -1
    sc = 1 if c1 >= c0 else -1
    err = dc - dr
    r, c = r0, c0
    while True:
        cells.append((r, c))
        if r == r1 and c == c1:
            break
        e2 = 2 * err
        if e2 > -dr:
            err -= dr; c += sc
        if e2 < dc:
            err += dc; r += sr
    return cells


def _lowest_near(dem, valid, r, c, radius, side_ok=None):
    """Lowest valid cell within `radius` of (r, c); `side_ok(rows, cols)`
    optionally restricts the search to one side of the road."""
    rows, cols = dem.shape
    r0, r1 = max(0, r - radius), min(rows, r + radius + 1)
    c0, c1 = max(0, c - radius), min(cols, c + radius + 1)
    ok = valid[r0:r1, c0:c1].copy()
    if side_ok is not None:
        rr, cc = np.mgrid[r0:r1, c0:c1]
        ok &= side_ok(rr, cc)
    win = np.where(ok, dem[r0:r1, c0:c1], np.inf)
    if not np.isfinite(win).any():
        return None
    k = np.unravel_index(int(np.argmin(win)), win.shape)
    return r0 + int(k[0]), c0 + int(k[1])


def burn_crossings(dem, valid, geotransform, crossings, half_length_m=30.0,
                   search_radius_cells=2, min_grade=1e-4):
    """Carve a breach through the road at every crossing.

    Parameters
    ----------
    dem, valid : raw elevation grid and validity mask (not modified)
    geotransform : GDAL-style, north-up
    crossings : list of dicts {x, y, tx, ty, id} - crossing point on the
        centreline and the unit road tangent (tx, ty) at that point
    half_length_m : breach extends this far either side of the centreline
        (embankment half-width + a margin; 20-40 m typical)
    search_radius_cells : each breach end snaps to the lowest cell within
        this radius, so it starts and ends in the channel
    min_grade : minimum fall (m/m) imposed along the breach when both ends
        are at the same height, so it drains without flat resolution

    Returns (burned_dem, log) where log has one dict per crossing:
    id, cells, max_cut_m, total_cut_m3 (volume), z_up, z_down, length_m,
    path (breach cell centres, upstream -> downstream) and note.
    """
    gt = tuple(geotransform)
    cw, ch = abs(gt[1]), abs(gt[5])
    rows, cols = dem.shape
    out = np.array(dem, dtype=np.float64, copy=True)
    log = []

    def rc(x, y):
        return int(math.floor((y - gt[3]) / gt[5])), int(math.floor((x - gt[0]) / gt[1]))

    for cr in crossings:
        if cr.get("tx") is None and cr.get("azimuth_deg") is not None:
            az = math.radians(float(cr["azimuth_deg"]))
            cr = dict(cr, tx=math.sin(az), ty=math.cos(az))
        tx, ty = float(cr["tx"]), float(cr["ty"])
        tl = math.hypot(tx, ty) or 1.0
        nx, ny = -ty / tl, tx / tl                    # left normal
        L = float(half_length_m)
        ends = []
        px, py = float(cr["x"]), float(cr["y"])
        margin = 0.5 * max(cw, ch)
        for sgn in (1.0, -1.0):
            r, c = rc(px + sgn * L * nx, py + sgn * L * ny)
            r = min(max(r, 0), rows - 1); c = min(max(c, 0), cols - 1)

            # each end stays on its own side of the centreline, so the two
            # search windows can never pick the same cell
            def side_ok(rr, cc, sgn=sgn):
                x = gt[0] + (cc + 0.5) * gt[1]
                y = gt[3] + (rr + 0.5) * gt[5]
                return sgn * ((x - px) * nx + (y - py) * ny) > margin
            ends.append(_lowest_near(out, valid, r, c, int(search_radius_cells), side_ok))
        entry = {"id": cr.get("id"), "cells": 0, "max_cut_m": 0.0, "total_cut_m3": 0.0,
                 "z_up": None, "z_down": None, "length_m": 0.0, "note": ""}
        if ends[0] is None or ends[1] is None:
            entry["note"] = "no valid cells at one end; not burned"
            log.append(entry); continue
        (ra, ca), (rb, cb) = ends
        za, zb = out[ra, ca], out[rb, cb]
        (ru, cu, zu), (rd, cd, zd) = ((ra, ca, za), (rb, cb, zb)) if za >= zb else \
            ((rb, cb, zb), (ra, ca, za))
        cells = _line_cells(ru, cu, rd, cd)
        steps = [0.0]
        for (r0, c0), (r1, c1) in zip(cells[:-1], cells[1:]):
            steps.append(steps[-1] + math.hypot((c1 - c0) * cw, (r1 - r0) * ch))
        length = steps[-1]
        if length <= 0:
            entry["note"] = "breach ends coincide; not burned"
            log.append(entry); continue
        z_end = min(zd, zu - min_grade * length)
        cut_max, cut_sum, n = 0.0, 0.0, 0
        for (r, c), s in zip(cells, steps):
            if not valid[r, c]:
                continue
            grade = zu + (z_end - zu) * (s / length)
            if out[r, c] > grade:
                cut = out[r, c] - grade
                out[r, c] = grade
                cut_max = max(cut_max, cut); cut_sum += cut; n += 1
        centre = lambda r, c: (gt[0] + (c + 0.5) * gt[1], gt[3] + (r + 0.5) * gt[5])
        entry.update(cells=n, max_cut_m=cut_max, total_cut_m3=cut_sum * cw * ch,
                     z_up=float(zu), z_down=float(z_end), length_m=length,
                     path=[centre(r, c) for r, c in cells])
        log.append(entry)
    return out, log
