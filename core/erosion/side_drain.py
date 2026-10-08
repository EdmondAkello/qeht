# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Side-drain siltation indicator along a road corridor (STI advisory R3, v0.20).

Side drains silt up where a low-gradient drain is fed by a slope with a
high sediment-transport capacity. The overland STI (Moore & Wilson 1992;
channel cells NoData, A_s capped) is a relative index of that capacity.

1. Stations every `step` m of chainage (the stations of sample_corridor).
   On each side cells are sampled at offsets offset_step ... half_width
   (nearest cell). LHS / RHS are left / right looking up-chainage, the
   sign convention of Alignment.locate.
2. A sampled cell is kept when its D8 path reaches the road strip (the
   cells within road_half_m of the centreline, or the cell the centreline
   passes through) within ceil(half_width / cell) + 1 steps without
   entering a channel cell. Cells in the strip itself are the road and
   are not sampled. Repeated cells at one station count once.
3. Per station and side: sti_p50 and sti_p90 of the kept cells' STI and
   n (cells kept); NULL when none is kept.
4. Relative classes from the 50 / 75 / 90th percentiles of this
   corridor's station p90 values (both sides pooled): low (<= P50),
   moderate, high, very high (> P90). "Relative transport capacity",
   never severity; scaling every STI value leaves the classes unchanged.
5. siltation = 1 where the class is high or very high and the
   |longitudinal ground slope| is below drain_slope_pct (default 1 %).
   A screening flag only.

Not used in the composite score, RUSLE or any class weight (rule R4).

No QGIS imports, no GDAL.
"""

import json
import math

import numpy as np

from ..grid import receivers_from_direction
from .corridor import _stations

CLASSES = ("low", "moderate", "high", "very high")
STATION_FIELDS = [("sti_p50_lhs", "real"), ("sti_p90_lhs", "real"), ("n_lhs", "int"),
                  ("sti_p50_rhs", "real"), ("sti_p90_rhs", "real"), ("n_rhs", "int"),
                  ("sti_class_lhs", "text"), ("sti_class_rhs", "text"),
                  ("slope_long_pct", "real"), ("siltation_lhs", "int"), ("siltation_rhs", "int")]
REACH_FIELDS = [("sti_p90_lhs", "real"), ("sti_p90_rhs", "real"), ("sti_max_lhs", "real"),
                ("sti_max_rhs", "real"), ("sti_class_lhs", "text"), ("sti_class_rhs", "text"),
                ("siltation_len_lhs_m", "real"), ("siltation_len_rhs_m", "real")]
NOTE = ("relative transport capacity of the slope draining to each side drain; screening only, "
        "not a severity and not used in the composite or RUSLE")


def _flat(gt, shape, xs, ys):
    rows, cols = shape
    c = np.floor((np.asarray(xs) - gt[0]) / gt[1]).astype(np.int64)
    r = np.floor((np.asarray(ys) - gt[3]) / gt[5]).astype(np.int64)
    ok = (r >= 0) & (r < rows) & (c >= 0) & (c < cols)
    return np.where(ok, r * cols + c, -1)


def road_strip(alignment, geotransform, shape, road_half_m=5.0):
    """Bool grid: cells within road_half_m of the centreline (sampled densely)."""
    gt = tuple(geotransform)
    cell = min(abs(gt[1]), abs(gt[5]))
    d = cell / 4.0
    strip = np.zeros(shape, bool).ravel()
    offs = np.union1d(np.arange(-road_half_m, road_half_m + 1e-9, d), [0.0]) \
        if road_half_m > 0 else np.array([0.0])
    for k in range(len(alignment.x0)):
        n = max(int(math.ceil(alignment.seg_len[k] / d)), 1)
        u = np.linspace(0.0, 1.0, n + 1)
        x = alignment.x0[k] + u * alignment.dx[k]
        y = alignment.y0[k] + u * alignment.dy[k]
        tx, ty = alignment.tangent(k)
        X = (x[:, None] + offs[None, :] * -ty).ravel()
        Y = (y[:, None] + offs[None, :] * tx).ravel()
        f = _flat(gt, shape, X, Y)
        strip[f[f >= 0]] = True
    return strip.reshape(shape)


def station_sti(alignment, sti, direction, valid, geotransform, channel=None, step=10.0,
                half_width=50.0, offset_step=10.0, road_half_m=5.0):
    """Per-station STI statistics on the cells that drain to each side of the road.

    sti: float grid, NaN on channel cells and NoData. direction: internal D8
    indices (NO_RECEIVER = -1). channel: bool grid (optional; NaN STI is not
    treated as channel, a channel mask stops the paths). Returns a list of
    dicts with chainage, x, y and the station fields of steps 2-3.
    """
    gt = tuple(geotransform)
    shape = direction.shape
    cell = min(abs(gt[1]), abs(gt[5]))
    rec = receivers_from_direction(np.where(valid, direction, -1), shape)
    strip = road_strip(alignment, gt, shape, road_half_m).ravel()
    chan = np.zeros(rec.size, bool) if channel is None else np.asarray(channel, bool).ravel()
    sti_f = np.asarray(sti, float).ravel()
    chs = _stations(alignment, float(step))
    k = np.clip(np.searchsorted(alignment.ch0, chs, side="right") - 1, 0, len(alignment.ch0) - 1)
    u = np.clip((chs - alignment.ch0[k]) / alignment.seg_len[k], 0.0, 1.0)
    x = alignment.x0[k] + u * alignment.dx[k]
    y = alignment.y0[k] + u * alignment.dy[k]
    tx, ty = alignment.tangent(k)
    nx, ny = -ty, tx                                     # left normal
    offs = np.arange(offset_step, half_width + 1e-9, offset_step) if half_width > 0 else np.array([])
    max_steps = int(math.ceil(half_width / cell)) + 1
    out = [{"chainage": float(chs[i]), "x": float(x[i]), "y": float(y[i])} for i in range(chs.size)]
    for side, sgn in (("lhs", 1.0), ("rhs", -1.0)):
        if offs.size == 0:
            for st in out:
                st.update({f"sti_p50_{side}": None, f"sti_p90_{side}": None, f"n_{side}": 0})
            continue
        X = x[:, None] + sgn * offs[None, :] * nx[:, None]
        Y = y[:, None] + sgn * offs[None, :] * ny[:, None]
        f = _flat(gt, shape, X, Y)
        start = f.ravel()
        ok = start >= 0
        ok[ok] &= ~strip[start[ok]] & ~chan[start[ok]]
        cur = np.where(ok, start, -1)
        kept = np.zeros(start.size, bool)
        active = ok.copy()
        for _ in range(max_steps):
            if not active.any():
                break
            nxt = np.where(active, rec[np.maximum(cur, 0)], -1)
            dead = active & (nxt < 0)
            active &= ~dead
            ii = np.flatnonzero(active)
            nn = nxt[ii]
            hit = strip[nn]
            kept[ii[hit]] = True
            active[ii[hit]] = False
            blocked = chan[nn] & ~hit
            active[ii[blocked]] = False
            cur[ii] = nn
        kept = kept.reshape(f.shape)
        for i, st in enumerate(out):
            cells = np.unique(f[i][kept[i]])
            v = sti_f[cells]
            v = v[np.isfinite(v)]
            st[f"n_{side}"] = int(v.size)
            if v.size:
                st[f"sti_p50_{side}"] = float(np.percentile(v, 50))
                st[f"sti_p90_{side}"] = float(np.percentile(v, 90))
            else:
                st[f"sti_p50_{side}"] = None
                st[f"sti_p90_{side}"] = None
    return out


def class_breaks(stations):
    """P50 / P75 / P90 of the station p90 values over both sides (None if none)."""
    v = [s[k] for s in stations for k in ("sti_p90_lhs", "sti_p90_rhs") if s.get(k) is not None]
    if not v:
        return None
    return [float(np.percentile(v, q)) for q in (50, 75, 90)]


def sti_class(v, breaks):
    if v is None or breaks is None:
        return None
    return CLASSES[int(np.searchsorted(np.asarray(breaks), v, side="left"))]


def classify(stations, breaks, slope_long_pct=None, drain_slope_pct=1.0):
    """Adds sti_class_*, slope_long_pct and siltation_* to each station in place.

    slope_long_pct: one value per station (None / NaN when unknown); without
    it the siltation flags are NULL."""
    for i, st in enumerate(stations):
        sl = None
        if slope_long_pct is not None:
            v = slope_long_pct[i]
            sl = float(v) if v is not None and np.isfinite(v) else None
        st["slope_long_pct"] = sl
        for side in ("lhs", "rhs"):
            c = sti_class(st.get(f"sti_p90_{side}"), breaks)
            st[f"sti_class_{side}"] = c
            if c is None or sl is None:
                st[f"siltation_{side}"] = None
            else:
                st[f"siltation_{side}"] = int(c in ("high", "very high") and abs(sl) < drain_slope_pct)
    return stations


def reach_sti(reaches, stations, breaks):
    """Reach fields from the stations inside each reach (sample_corridor reaches:
    ch_start / ch_end; station spacing lengths from the half-way bounds)."""
    chs = np.array([s["chainage"] for s in stations], float)
    if chs.size:
        bounds = np.r_[chs[0], 0.5 * (chs[1:] + chs[:-1]), chs[-1]]
        seg = np.diff(bounds)
    out = []
    for r in reaches:
        sel = [i for i in range(chs.size) if r["ch_start"] - 1e-9 <= chs[i] <= r["ch_end"] + 1e-9
               and (chs[i] < r["ch_end"] - 1e-9 or i == chs.size - 1)]
        b = {}
        for side in ("lhs", "rhs"):
            v = [stations[i][f"sti_p90_{side}"] for i in sel
                 if stations[i].get(f"sti_p90_{side}") is not None]
            b[f"sti_p90_{side}"] = float(np.percentile(v, 90)) if v else None
            b[f"sti_max_{side}"] = float(max(v)) if v else None
            b[f"sti_class_{side}"] = sti_class(b[f"sti_p90_{side}"], breaks)
            flags = [stations[i].get(f"siltation_{side}") for i in sel]
            b[f"siltation_len_{side}_m"] = (None if all(f is None for f in flags) else
                                           float(sum(seg[i] for i in sel
                                                     if stations[i].get(f"siltation_{side}") == 1)))
        out.append(b)
    return out


def corridor_sti(alignment, sti, direction, valid, geotransform, reaches, channel=None,
                 slope_at=None, step=10.0, half_width=50.0, offset_step=10.0, road_half_m=5.0,
                 drain_slope_pct=1.0):
    """Stations, reach fields and class breaks in one call.

    slope_at: callable chainage array -> longitudinal slope % array (or None)."""
    st = station_sti(alignment, sti, direction, valid, geotransform, channel, step, half_width,
                     offset_step, road_half_m)
    br = class_breaks(st)
    sl = slope_at(np.array([s["chainage"] for s in st])) if slope_at is not None and st else None
    classify(st, br, sl, drain_slope_pct)
    return st, reach_sti(reaches, st, br), br


def params_json(step, half_width, offset_step, road_half_m, drain_slope_pct, breaks,
                slope_source):
    return json.dumps({"step_m": step, "half_width_m": half_width, "offset_step_m": offset_step,
                       "road_half_m": road_half_m, "drain_slope_pct": drain_slope_pct,
                       "class_breaks_p50_p75_p90": breaks, "classes": list(CLASSES),
                       "slope_source": slope_source, "note": NOTE}, sort_keys=True)
