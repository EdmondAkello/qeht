# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Alignment ground profile (A1, v0.15).

Stations every `step` metres of chainage along a road centreline. At each
station:

    z_dem_m        raw DEM ground, bilinear between cell centres
    z_fill_m       filled (conditioned) DEM, bilinear; z_fill - z_dem is
                   the ponding depth the fill removed at the station
    acc_km2        contributing area: the largest (accumulation + 1) x cell
                   area within `search_cells` of the station's cell, so a
                   stream line one cell off the centreline is not missed
    stream         1 at the station nearest to each place where a D8
                   stream link crosses the alignment (the same exact link x
                   alignment intersection as 'Road crossing candidates'),
                   else 0
    strahler       Strahler order of that stream (when a stream-order grid
                   is given)
    slope_long_pct longitudinal ground slope along the alignment, centred
                   difference of z_dem_m over the neighbouring stations
                   (one-sided at the ends), positive = rising chainage

It is a terrain profile for provisional culvert and relief levels before
the geometric design exists, and a check on DEM height bias once the
design existing-ground profile is available. Not a design long section.

No QGIS imports, no GDAL.
"""

import math

import numpy as np

PROFILE_FIELDS = [
    ("align_name", "text"), ("chainage_m", "real"), ("x", "real"), ("y", "real"),
    ("z_dem_m", "real"), ("z_fill_m", "real"), ("pond_depth_m", "real"),
    ("acc_km2", "real"), ("stream", "int"), ("strahler", "int"),
    ("slope_long_pct", "real"), ("alignment_part", "int"),
]


def stations(alignment, step):
    """Chainages every `step` m from the start, plus the end chainage."""
    if step <= 0:
        raise ValueError("Station spacing must be positive.")
    total = alignment.length
    n = int(np.floor(total / step + 1e-9)) + 1
    chs = alignment.start_chainage + np.arange(n) * float(step)
    if total - (n - 1) * step > 1e-6:
        chs = np.append(chs, alignment.end_chainage)
    return chs


def bilinear(grid, gt, xs, ys):
    """Bilinear interpolation between cell centres; NaN outside or next to NoData.

    `grid` is float with NaN for NoData. Points within half a cell of the
    edge use the nearest edge cells (clamped), so the value is defined over
    the whole raster extent.
    """
    rows, cols = grid.shape
    fc = (np.asarray(xs, float) - gt[0]) / gt[1] - 0.5
    fr = (np.asarray(ys, float) - gt[3]) / gt[5] - 0.5
    inside = (fc >= -0.5) & (fc <= cols - 0.5) & (fr >= -0.5) & (fr <= rows - 0.5)
    fc = np.clip(fc, 0.0, cols - 1.0)
    fr = np.clip(fr, 0.0, rows - 1.0)
    c0 = np.minimum(np.floor(fc).astype(np.int64), max(cols - 2, 0))
    r0 = np.minimum(np.floor(fr).astype(np.int64), max(rows - 2, 0))
    c1 = np.minimum(c0 + 1, cols - 1)
    r1 = np.minimum(r0 + 1, rows - 1)
    u = fc - c0
    v = fr - r0
    z = ((1 - u) * (1 - v) * grid[r0, c0] + u * (1 - v) * grid[r0, c1]
         + (1 - u) * v * grid[r1, c0] + u * v * grid[r1, c1])
    return np.where(inside, z, np.nan)


def _window_max(grid, gt, xs, ys, k):
    rows, cols = grid.shape
    c = np.floor((np.asarray(xs) - gt[0]) / gt[1]).astype(np.int64)
    r = np.floor((np.asarray(ys) - gt[3]) / gt[5]).astype(np.int64)
    out = np.full(len(c), np.nan)
    for i, (rr, cc) in enumerate(zip(r, c)):
        r0, r1 = max(rr - k, 0), min(rr + k + 1, rows)
        c0, c1 = max(cc - k, 0), min(cc + k + 1, cols)
        if r0 >= r1 or c0 >= c1:
            continue
        w = grid[r0:r1, c0:c1]
        if np.isfinite(w).any():
            out[i] = float(np.nanmax(w))
    return out


def longitudinal_slope_pct(chainage, z):
    """Centred-difference slope in % (one-sided at the ends, NaN-safe)."""
    ch = np.asarray(chainage, float)
    z = np.asarray(z, float)
    n = ch.size
    s = np.full(n, np.nan)
    if n < 2:
        return s
    for i in range(n):
        a, b = max(i - 1, 0), min(i + 1, n - 1)
        if b == a or ch[b] == ch[a] or not (np.isfinite(z[a]) and np.isfinite(z[b])):
            continue
        s[i] = 100.0 * (z[b] - z[a]) / (ch[b] - ch[a])
    return s


def alignment_profile(alignment, raw, geotransform, accumulation=None, valid=None,
                      filled=None, step=10.0, search_cells=1, crossings=None,
                      align_name=""):
    """Ground profile along `alignment` (see the module docstring).

    raw, filled : float grids with NaN for NoData (filled optional)
    accumulation : accumulation in CELLS (QEHT convention, excludes the cell);
        optional - acc_km2 is empty without it
    valid : bool grid for the accumulation (cells outside are ignored)
    crossings : optional list of dicts {chainage_m, acc_km2, stream_order}
        from core.network.crossings.find_crossing_candidates - the exact D8
        stream x alignment intersections. Each flags the station nearest
        its chainage (stream = 1); that station's acc_km2 is raised to the
        crossing's area if larger and takes its Strahler order.

    Returns a list of dicts (PROFILE_FIELDS keys).
    """
    gt = tuple(geotransform)
    cell_area = abs(gt[1] * gt[5])
    chs = stations(alignment, float(step))
    pts = [alignment.point_at(c) for c in chs]
    xs = np.array([p[0] for p in pts])
    ys = np.array([p[1] for p in pts])
    _, _, seg = alignment.locate(xs, ys)
    part = alignment.part[seg]

    z = bilinear(np.asarray(raw, float), gt, xs, ys)
    zf = bilinear(np.asarray(filled, float), gt, xs, ys) if filled is not None \
        else np.full(xs.size, np.nan)
    if accumulation is not None:
        acc = np.asarray(accumulation, float)
        if valid is not None:
            acc = np.where(valid, acc, np.nan)
        a = _window_max(acc, gt, xs, ys, int(search_cells))
        acc_km2 = (a + 1.0) * cell_area / 1.0e6
    else:
        acc_km2 = np.full(xs.size, np.nan)

    stream = np.zeros(xs.size, dtype=int)
    order = [None] * xs.size
    for c in (crossings or []):
        ch = c.get("chainage_m")
        if ch is None or not np.isfinite(ch):
            continue
        i = int(np.argmin(np.abs(chs - ch)))
        stream[i] = 1
        ca = c.get("acc_km2")
        if ca is not None and np.isfinite(ca) and not (np.isfinite(acc_km2[i]) and acc_km2[i] >= ca):
            acc_km2[i] = ca
        so = c.get("stream_order")
        if so is not None and (order[i] is None or so > order[i]):
            order[i] = int(so)

    # slope within each alignment part only (no difference across a gap)
    slope = np.full(xs.size, np.nan)
    for k in np.unique(part):
        sel = np.flatnonzero(part == k)
        slope[sel] = longitudinal_slope_pct(chs[sel], z[sel])

    def _f(v):
        return float(v) if np.isfinite(v) else None

    rows = []
    for i in range(xs.size):
        pond = zf[i] - z[i] if np.isfinite(zf[i]) and np.isfinite(z[i]) else float("nan")
        rows.append({
            "align_name": align_name, "chainage_m": float(chs[i]),
            "x": float(xs[i]), "y": float(ys[i]),
            "z_dem_m": _f(z[i]), "z_fill_m": _f(zf[i]),
            "pond_depth_m": _f(max(pond, 0.0)) if np.isfinite(pond) else None,
            "acc_km2": _f(acc_km2[i]), "stream": int(stream[i]), "strahler": order[i],
            "slope_long_pct": _f(slope[i]), "alignment_part": int(part[i]) + 1,
        })
    return rows


def profile_summary(rows):
    """Short text summary for the run log."""
    if not rows:
        return "no stations"
    z = [r["z_dem_m"] for r in rows if r["z_dem_m"] is not None]
    n_st = sum(r["stream"] for r in rows)
    pond = [r["pond_depth_m"] for r in rows if r["pond_depth_m"]]
    txt = (f"{len(rows)} stations, ch {rows[0]['chainage_m']:,.0f} - {rows[-1]['chainage_m']:,.0f} m; "
           f"{n_st} stream crossing(s)")
    if z:
        txt += f"; ground {min(z):,.1f} - {max(z):,.1f} m"
    if pond:
        txt += f"; {len(pond)} station(s) in filled depressions (max {max(pond):.2f} m)"
    return txt


def profile_chart(path, rows, title="Alignment ground profile", max_labels=12):
    """PNG of the ground profile with stream crossings marked (matplotlib optional).

    Only the `max_labels` largest crossings get an area label.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return False
    ch = np.array([r["chainage_m"] for r in rows]) / 1000.0
    z = np.array([np.nan if r["z_dem_m"] is None else r["z_dem_m"] for r in rows])
    zf = np.array([np.nan if r["z_fill_m"] is None else r["z_fill_m"] for r in rows])
    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(ch, z, color="#5b4636", lw=1.2, label="DEM ground")
    if np.isfinite(zf).any():
        ax.plot(ch, zf, color="#2b7bb9", lw=0.8, ls="--", label="filled DEM")
    st = np.array([r["stream"] for r in rows]) == 1
    if st.any():
        ax.scatter(ch[st], z[st], marker="v", color="#1f78b4", zorder=3, label="stream crossing")
        idx = [i for i in np.flatnonzero(st)
               if rows[i]["acc_km2"] is not None and math.isfinite(rows[i]["acc_km2"])]
        # label only the largest crossings, so the chart stays readable
        idx = sorted(sorted(idx, key=lambda i: -rows[i]["acc_km2"])[:max_labels])
        for n, i in enumerate(idx):
            a = rows[i]["acc_km2"]
            if True:
                ax.annotate(f"{a:.2g} km²", (ch[i], z[i]), textcoords="offset points",
                            xytext=(0, -14 - 10 * (n % 2)), ha="center", fontsize=7, color="#1f78b4")
    ax.set_xlabel("Chainage (km)")
    ax.set_ylabel("Elevation (m)")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return True


def profile_with_crossings(alignment, raw, geotransform, direction=None, valid=None,
                           accumulation=None, filled=None, stream_threshold_cells=None,
                           stream_order=None, step=10.0, search_cells=1, align_name=""):
    """Profile plus the exact D8 stream x alignment crossings.

    With direction, valid, accumulation and a stream threshold the stream
    network is extracted at that threshold and intersected with the
    alignment (as 'Road crossing candidates' does, without a minimum area);
    otherwise stream stays 0 everywhere. Returns (rows, crossings).
    """
    crossings = []
    if direction is not None and accumulation is not None and stream_threshold_cells:
        from ..watershed.delineate import extract_streams
        from .crossings import find_crossing_candidates
        mask = extract_streams(accumulation, valid, threshold_cells=float(stream_threshold_cells))
        cand, _, _ = find_crossing_candidates(
            direction, valid, accumulation, mask, geotransform, alignment,
            min_area_km2=0.0, corridor_halfwidth_m=0.0, min_parallel_m=1e12,
            stream_order=stream_order)
        crossings = [a for _, a in cand]
    rows = alignment_profile(alignment, raw, geotransform, accumulation=accumulation,
                             valid=valid, filled=filled, step=step,
                             search_cells=search_cells, crossings=crossings,
                             align_name=align_name)
    return rows, crossings
