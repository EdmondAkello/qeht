# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Floodplain width indicator at major crossings (A4, v0.18).

A terrain-only first look at how wide the valley floor is along the road
at a large crossing, to inform the split of the check flood between the
main structure and relief culverts. NOT a flood level: use it with a
hydraulic calculation.

Profile method (always). From the alignment ground profile (A1, raw DEM,
bilinear): z_bed is the lowest ground within `bed_window_m` of the
crossing's chainage. For each dz the width is the contiguous run of the
profile around that low point with z <= z_bed + dz; the two ends are
interpolated linearly between the stations either side, so a linear bank is
measured exactly.

HAND method (optional, second estimate). Height Above Nearest Drainage:
every cell's elevation minus the elevation of the first stream cell reached
along its D8 path (Renno et al. 2008; Nobre et al. 2011). Sampled at the
stations, the width is the contiguous run around the crossing with
HAND <= dz. HAND follows the drainage line, so on a valley that slopes
along the road it differs from the profile method; the gap is information.

Widths that reach the end of the alignment, or a NoData station, are
truncated and say so in fp_note.

No QGIS imports, no GDAL.
"""

import math

import numpy as np

from ..grid import receivers_from_direction

DZ = (0.5, 1.0, 2.0)
FP_FIELDS = ["fp_w_0p5_m", "fp_w_1p0_m", "fp_w_2p0_m", "fp_ch_from_m", "fp_ch_to_m",
             "fp_z_bed_m", "fp_hand_w_0p5_m", "fp_hand_w_1p0_m", "fp_hand_w_2p0_m",
             "fp_method", "fp_note"]


def _key(dz):
    return f"{dz:.1f}".replace(".", "p")


def run_width(ch, v, i0, level):
    """Contiguous run around station i0 with v <= level.

    Returns (from_ch, to_ch, truncated_left, truncated_right). Ends are
    interpolated between the last station inside and the first outside.
    A NaN station or the end of the profile stops the run (truncated).
    """
    n = len(ch)
    if not (0 <= i0 < n) or not np.isfinite(v[i0]) or v[i0] > level:
        return None, None, False, False
    lo = i0
    while lo - 1 >= 0 and np.isfinite(v[lo - 1]) and v[lo - 1] <= level:
        lo -= 1
    hi = i0
    while hi + 1 < n and np.isfinite(v[hi + 1]) and v[hi + 1] <= level:
        hi += 1
    tl = lo == 0 or not np.isfinite(v[lo - 1])
    tr = hi == n - 1 or not np.isfinite(v[hi + 1])
    a = ch[lo] if tl else ch[lo - 1] + (level - v[lo - 1]) / (v[lo] - v[lo - 1]) * (ch[lo] - ch[lo - 1])
    b = ch[hi] if tr else ch[hi] + (level - v[hi]) / (v[hi + 1] - v[hi]) * (ch[hi + 1] - ch[hi])
    return float(a), float(b), tl, tr


def profile_widths(chainage, z, crossing_ch, dzs=DZ, bed_window_m=50.0):
    """Floodplain widths by the profile method -> dict of FP fields (profile part)."""
    ch = np.asarray(chainage, float)
    zz = np.asarray(z, float)
    out = {}
    win = np.where(np.abs(ch - crossing_ch) <= bed_window_m)[0]
    win = win[np.isfinite(zz[win])]
    if win.size == 0:
        return out, "no ground profile near the crossing"
    i0 = int(win[np.argmin(zz[win])])
    z_bed = float(zz[i0])
    out["fp_z_bed_m"] = z_bed
    notes = []
    for dz in dzs:
        a, b, tl, tr = run_width(ch, zz, i0, z_bed + dz)
        out[f"fp_w_{_key(dz)}_m"] = None if a is None else b - a
        if math.isclose(dz, 1.0):
            out["fp_ch_from_m"], out["fp_ch_to_m"] = a, b
        if tl or tr:
            notes.append(f"dz {dz:g} m reaches the {'start' if tl else 'end'} of the profile "
                         "(width is a lower bound)")
    return out, "; ".join(notes)


def hand_grid(direction, valid, elevation, stream_mask):
    """Height above nearest drainage (m) on the D8 grid; NaN where no stream
    is reached (path ends in a pit or off the grid without a stream cell).

    Pointer jumping: every cell points to its receiver until a stream cell,
    which points to itself; doubling the pointers converges in log2(path
    length) passes.
    """
    shape = direction.shape
    rec = receivers_from_direction(np.where(valid, direction, -1), shape)
    n = rec.size
    idx = np.arange(n, dtype=np.int64)
    st = np.asarray(stream_mask, bool).ravel() & np.asarray(valid, bool).ravel()
    ptr = np.where(st | (rec < 0), idx, rec)
    for _ in range(64):
        nxt = ptr[ptr]
        if np.array_equal(nxt, ptr):
            break
        ptr = nxt
    z = np.asarray(elevation, float).ravel()
    reached = st[ptr]
    hand = np.where(reached & np.asarray(valid, bool).ravel(), z - z[ptr], np.nan)
    return hand.reshape(shape)


def hand_widths(chainage, hand_at_stations, crossing_ch, dzs=DZ, bed_window_m=50.0):
    """Widths with HAND <= dz around the station of least HAND near the crossing."""
    ch = np.asarray(chainage, float)
    h = np.asarray(hand_at_stations, float)
    out = {}
    win = np.where(np.abs(ch - crossing_ch) <= bed_window_m)[0]
    win = win[np.isfinite(h[win])]
    if win.size == 0:
        return out
    i0 = int(win[np.argmin(h[win])])
    for dz in dzs:
        a, b, _, _ = run_width(ch, h, i0, dz)
        out[f"fp_hand_w_{_key(dz)}_m"] = None if a is None else b - a
    return out


def floodplain_block(crossings, profile_rows, min_area_km2=10.0, hand_at_stations=None,
                     bed_window_m=50.0, dzs=DZ):
    """FP fields for each crossing attrs dict (needs chainage_m and
    acc_at_outlet_km2). profile_rows: the A1 rows (chainage_m, z_dem_m).
    Returns a list of dicts in the order of `crossings`."""
    ch = np.array([r["chainage_m"] for r in profile_rows], float)
    z = np.array([np.nan if r.get("z_dem_m") is None else r["z_dem_m"] for r in profile_rows], float)
    order = np.argsort(ch, kind="stable")
    ch, z = ch[order], z[order]
    hs = None if hand_at_stations is None else np.asarray(hand_at_stations, float)[order]
    out = []
    for c in crossings:
        b = {k: None for k in FP_FIELDS}
        area = c.get("acc_at_outlet_km2")
        cch = c.get("chainage_m")
        if area is None or cch is None or area < min_area_km2 or ch.size < 2:
            out.append(b)
            continue
        pw, note = profile_widths(ch, z, cch, dzs, bed_window_m)
        b.update(pw)
        b["fp_method"] = "profile"
        if hs is not None:
            b.update(hand_widths(ch, hs, cch, dzs, bed_window_m))
            b["fp_method"] = "profile+HAND"
        b["fp_note"] = note or None
        out.append(b)
    return out
