# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Boolean cell mask -> (multi)polygon, in pure NumPy / Python.

Why not gdal.Polygonize()
-------------------------
gdal.Polygonize works on a single label raster, so overlapping catchments
(the default for road crossings: each culvert gets its FULL upstream
area) cannot be represented - in v0.8.3 the characteristics tool burned
catchments into one grid first-come-first-served, and a downstream
catchment listed after an upstream one lost the upstream cells from its
polygon (its area attribute was still right). Tracing each catchment's
own mask removes that failure, needs no GDAL (so the exchange package is
testable on bare Python), and gives identical geometry on every platform.

Method
------
Every cell side between a mask cell and a non-mask cell becomes a directed
edge with the mask on its LEFT (y-up map sense). Edges are chained into
closed rings; where two cells meet only at a corner the tracer turns left,
so diagonally-touching cells become separate parts (4-connectivity, same
as gdal.Polygonize's default). Collinear vertices are dropped. Outer rings
come out counter-clockwise, holes clockwise; each hole is assigned to the
smallest outer ring containing a point just inside it.

Vertices sit on cell corners, so area = cell count x cell area exactly.
No QGIS imports, no GDAL.
"""

import numpy as np

# Direction codes for edges, as (drow, dcol) in grid space (row down).
_DIRS = {(0, -1): 0, (0, 1): 1, (1, 0): 2, (-1, 0): 3}   # W, E, S, N
_DVEC = [(0, -1), (0, 1), (1, 0), (-1, 0)]


def _left_turn(d):
    dr, dc = _DVEC[d]
    return _DIRS[(-dc, dr)]


def _edges(mask):
    """Directed boundary edges: arrays (r0, c0, dir, out_r, out_c)."""
    m = np.pad(np.asarray(mask, dtype=bool), 1)
    core = m[1:-1, 1:-1]
    out = []
    # top side: neighbour above absent; edge (r, c+1) -> (r, c), heading W
    rr, cc = np.nonzero(core & ~m[:-2, 1:-1])
    out.append((rr, cc + 1, np.full(rr.size, 0), rr - 1, cc))
    # bottom side: (r+1, c) -> (r+1, c+1), heading E
    rr, cc = np.nonzero(core & ~m[2:, 1:-1])
    out.append((rr + 1, cc, np.full(rr.size, 1), rr + 1, cc))
    # left side: (r, c) -> (r+1, c), heading S
    rr, cc = np.nonzero(core & ~m[1:-1, :-2])
    out.append((rr, cc, np.full(rr.size, 2), rr, cc - 1))
    # right side: (r+1, c+1) -> (r, c+1), heading N
    rr, cc = np.nonzero(core & ~m[1:-1, 2:])
    out.append((rr + 1, cc + 1, np.full(rr.size, 3), rr, cc + 1))
    cat = [np.concatenate([o[k] for o in out]).astype(np.int64) for k in range(5)]
    return cat


def _signed_area(ring):
    """Shoelace area of a closed ring given as (row, col) corners, y-up."""
    a = np.asarray(ring, dtype=np.float64)
    x, y = a[:, 1], -a[:, 0]
    return 0.5 * float(np.sum(x[:-1] * y[1:] - x[1:] * y[:-1]))


def _point_in_ring(pr, pc, ring):
    """Even-odd test for point (row, col) against a ring of (row, col)."""
    a = np.asarray(ring, dtype=np.float64)
    r0, c0 = a[:-1, 0], a[:-1, 1]
    r1, c1 = a[1:, 0], a[1:, 1]
    cross = (r0 > pr) != (r1 > pr)
    with np.errstate(divide="ignore", invalid="ignore"):
        cx = c0 + (pr - r0) * (c1 - c0) / (r1 - r0)
    return int(np.count_nonzero(cross & (pc < cx))) % 2 == 1


def trace_rings(mask):
    """Simple rings in grid-corner coordinates. Returns list of (ring, pt).

    ring : list of (row, col) corners, closed (first == last), never
           touching itself (a boundary that pinches at a corner is split
           there into separate rings, so every ring is OGC-simple)
    pt   : (row, col) of a cell centre on the ring's non-mask side, next
           to one of its edges - inside the hole for a hole ring
    """
    r0, c0, d, orr, occ = _edges(mask)
    n = r0.size
    if n == 0:
        return []
    by_start = {}
    for i in range(n):
        by_start.setdefault((int(r0[i]), int(c0[i])), []).append(i)
    used = np.zeros(n, dtype=bool)
    rings = []
    for start in range(n):
        if used[start]:
            continue
        # 1) follow edges until we are back at the start edge
        seq = []
        e = start
        while True:
            used[e] = True
            seq.append(e)
            dr, dc = _DVEC[int(d[e])]
            end = (int(r0[e]) + dr, int(c0[e]) + dc)
            outs = by_start[end]
            if len(outs) == 1:
                nxt = outs[0]
            else:
                # corner where cells touch diagonally: turn left, i.e.
                # keep hugging the current cell (4-connected regions)
                want = _left_turn(int(d[e]))
                nxt = [k for k in outs if int(d[k]) == want][0]
            if nxt == start:
                break
            if used[nxt]:
                raise RuntimeError("polygon tracer revisited an edge")
            e = nxt
        # 2) split at repeated start corners -> simple loops
        stack, pos = [], {}
        for e in seq:
            v = (int(r0[e]), int(c0[e]))
            if v in pos:
                k = pos[v]
                loop = stack[k:]
                del stack[k:]
                for ee in loop:
                    pos.pop((int(r0[ee]), int(c0[ee])), None)
                rings.append(_loop_to_ring(loop, r0, c0, d, orr, occ))
            pos[v] = len(stack)
            stack.append(e)
        if stack:
            rings.append(_loop_to_ring(stack, r0, c0, d, orr, occ))
    return rings


def _loop_to_ring(loop, r0, c0, d, orr, occ):
    """Edge-id loop -> (closed simplified ring, non-mask-side point)."""
    verts = [(int(r0[e]), int(c0[e])) for e in loop]
    dirs = [int(d[e]) for e in loop]
    m = len(loop)
    # keep a vertex only where the direction changes
    ring = [verts[k] for k in range(m) if dirs[k] != dirs[k - 1]]
    if not ring:                      # cannot happen for a closed loop
        ring = verts[:1]
    ring.append(ring[0])
    e = loop[0]
    return ring, (int(orr[e]) + 0.5, int(occ[e]) + 0.5)


def mask_to_polygons(mask, geotransform):
    """Trace `mask` and return a list of polygons in map coordinates.

    Each polygon is [outer_ring, hole_ring, ...]; each ring is a list of
    (x, y) tuples, closed. Outer rings are counter-clockwise, holes
    clockwise (in a north-up map). An empty mask returns [].
    """
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return []
    # Crop to the bounding box for speed on large grids.
    rows = np.flatnonzero(mask.any(axis=1))
    cols = np.flatnonzero(mask.any(axis=0))
    ra, rb, ca, cb = rows[0], rows[-1] + 1, cols[0], cols[-1] + 1
    rings = trace_rings(mask[ra:rb, ca:cb])

    outers, holes = [], []
    for ring, pt in rings:
        area = _signed_area(ring)
        (outers if area > 0 else holes).append((ring, pt, abs(area)))
    polys = [[o[0]] for o in outers]
    for ring, pt, _ in holes:
        best, best_area = None, None
        for k, (oring, _, oarea) in enumerate(outers):
            if _point_in_ring(pt[0], pt[1], oring):
                if best is None or oarea < best_area:
                    best, best_area = k, oarea
        if best is not None:
            polys[best].append(ring)

    gt = geotransform
    def to_xy(ring):
        return [(gt[0] + (c + ca) * gt[1] + (r + ra) * gt[2],
                 gt[3] + (c + ca) * gt[4] + (r + ra) * gt[5]) for r, c in ring]
    # deterministic order: by outer ring's first vertex (row, col)
    polys.sort(key=lambda p: (p[0][0][0], p[0][0][1]))
    return [[to_xy(r) for r in poly] for poly in polys]


def polygons_area(polys):
    """Planar area of a list of polygons from mask_to_polygons."""
    total = 0.0
    for poly in polys:
        for k, ring in enumerate(poly):
            a = np.asarray(ring, dtype=np.float64)
            s = 0.5 * float(np.sum(a[:-1, 0] * a[1:, 1] - a[1:, 0] * a[:-1, 1]))
            total += s   # outer +, holes - (orientation carries the sign)
    return total
