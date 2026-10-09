# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Drainage coverage along the road: missing crossings, sag points, flat
stretches (A2, A3; v0.16).

All three read the alignment ground profile (core.network.profile).

Missing crossings (A2)
    Stations where a D8 stream crosses (stream = 1) with contributing area
    >= min_area_km2 are merged within `merge_m` of chainage into one point
    at the largest area. A point with no existing crossing within
    `search_m` is UNCOVERED; the design hydrology package delineates and
    characterises it in the same run as a PROPOSED crossing (outlet_uid
    prefix P), so the designer can adopt it without another QEHT run.

Sag points (A3)
    Local minima of the ground profile after a moving-average smoothing over
    `smooth_m`, kept when their depth - the lower of the highest smoothed
    ground to the left and to the right before the profile drops below the
    minimum again (topographic prominence) - is at least `min_depth_m`.
    Water ponds there against the embankment even without a flow path:
    the place for an equaliser / relief culvert. The local area draining to
    each sag is measured with the alignment as a wall on the filled DEM
    (`walled_accumulation`), on each side; the larger side is reported.

Flat stretches (A3)
    Runs of stations where |longitudinal slope| and the cross-fall to
    both sides over `crossfall_m` are below `flat_slope_pct`, longer than
    `min_len_m`. They mark where nominal relief culverts for unaccounted
    sheet flow may be needed (common practice: single pipes at a maximum
    spacing, e.g. 200 m). QEHT does not place or size them.

No QGIS imports, no GDAL.
"""

import math

import numpy as np

from .profile import bilinear

COVERAGE_FIELDS = [("issue", "text"), ("chainage_m", "real"), ("chainage_to_m", "real"),
                   ("area_km2", "real"), ("uid", "text"), ("nearest_uid", "text"),
                   ("nearest_m", "real"), ("note", "text")]
SAG_FIELDS = [("chainage_m", "real"), ("z_dem_m", "real"), ("sag_depth_m", "real"),
              ("sag_area_km2", "real"), ("sag_area_left_km2", "real"),
              ("sag_area_right_km2", "real"), ("side", "text"), ("pond_depth_m", "real"),
              ("uid", "text")]
FLAT_FIELDS = [("chainage_m", "real"), ("chainage_to_m", "real"), ("length_m", "real"),
               ("slope_long_pct", "real"), ("crossfall_pct", "real"),
               ("n_crossings_within", "int"), ("uids_within", "text")]


def _nearest(ch, existing):
    """(uid, distance) of the nearest existing crossing by chainage."""
    best = (None, None)
    for uid, c in existing:
        if c is None:
            continue
        d = abs(float(c) - ch)
        if best[1] is None or d < best[1]:
            best = (uid, d)
    return best


def cluster_cover(candidates, used):
    """Chainages covered through a candidate cluster.

    A stream that runs beside the road crosses the centreline many times on a DEM; Road crossing
    candidates puts those intersections in one cluster and recommends one of them. The others
    are the same drainage line, so they are covered by the crossing used for the cluster.

    candidates : [(cluster_id, chainage_m)] for every candidate
    used : [(cluster_id, uid)] for the candidates used as crossings
    Returns [(uid, chainage_m)] for every candidate whose cluster holds a used crossing.
    """
    by_cluster = {}
    for cl, uid in used:
        if cl is not None:
            by_cluster.setdefault(cl, uid)
    return [(by_cluster[cl], float(ch)) for cl, ch in candidates
            if cl is not None and ch is not None and cl in by_cluster]


def missing_crossings(profile, existing, min_area_km2=0.0, merge_m=30.0, search_m=50.0):
    """Uncovered stream crossings.

    profile : rows from alignment_profile (chainage_m, x, y, acc_km2, stream)
    existing : [(uid, chainage_m)] of the crossings already in the design
    Returns a list of dicts: chainage_m, x, y, area_km2, nearest_uid,
    nearest_m, covered (bool). Only stations with stream = 1 count.
    """
    pts = [r for r in profile if r.get("stream") == 1 and r.get("acc_km2") is not None
           and r["acc_km2"] >= float(min_area_km2)]
    groups = []
    for r in sorted(pts, key=lambda r: r["chainage_m"]):
        if groups and r["chainage_m"] - groups[-1][-1]["chainage_m"] <= float(merge_m):
            groups[-1].append(r)
        else:
            groups.append([r])
    out = []
    for g in groups:
        best = max(g, key=lambda r: r["acc_km2"])
        uid, d = _nearest(best["chainage_m"], existing)
        out.append({"chainage_m": best["chainage_m"], "x": best["x"], "y": best["y"],
                    "area_km2": best["acc_km2"], "nearest_uid": uid, "nearest_m": d,
                    "covered": d is not None and d <= float(search_m),
                    "strahler": best.get("strahler")})
    return out


def smooth(values, chainage, window_m):
    """Centred moving average over `window_m` of chainage (NaN-aware)."""
    v = np.asarray(values, float)
    ch = np.asarray(chainage, float)
    if window_m <= 0 or v.size < 3:
        return v.copy()
    out = np.full(v.size, np.nan)
    half = window_m / 2.0
    lo = np.searchsorted(ch, ch - half, side="left")
    hi = np.searchsorted(ch, ch + half, side="right")
    for i in range(v.size):
        w = v[lo[i]:hi[i]]
        w = w[np.isfinite(w)]
        if w.size:
            out[i] = w.mean()
    return out


def sag_points(profile, smooth_m=30.0, min_depth_m=0.3):
    """Sags of the ground profile: [{index, chainage_m, x, y, z_dem_m, sag_depth_m}].

    Depth = prominence of the smoothed minimum: min(highest smoothed ground to
    the left, highest to the right) minus the minimum, each searched until the
    profile drops lower than the minimum again (or the alignment ends).
    """
    ch = np.array([r["chainage_m"] for r in profile], float)
    z = np.array([np.nan if r.get("z_dem_m") is None else r["z_dem_m"] for r in profile], float)
    zs = smooth(z, ch, smooth_m)
    n = zs.size
    out = []
    for i in range(1, n - 1):
        if not np.isfinite(zs[i]):
            continue
        # local minimum (plateaus: take the first cell of an equal run)
        if not (zs[i] < zs[i - 1] and zs[i] <= zs[i + 1]):
            continue
        j = i + 1
        while j < n and np.isfinite(zs[j]) and zs[j] == zs[i]:
            j += 1
        if j < n and np.isfinite(zs[j]) and zs[j] < zs[i]:
            continue
        k = i - 1
        while k >= 0 and not (np.isfinite(zs[k]) and zs[k] < zs[i]):
            k -= 1
        left = np.nanmax(zs[k + 1:i]) if i - (k + 1) > 0 else zs[i]
        k = i + 1
        while k < n and not (np.isfinite(zs[k]) and zs[k] < zs[i]):
            k += 1
        right = np.nanmax(zs[i + 1:k]) if k - (i + 1) > 0 else zs[i]
        depth = min(left, right) - zs[i]
        if depth >= float(min_depth_m):
            # the lowest raw station inside the plateau / near the smoothed minimum
            lo_i = max(0, i - 3)
            hi_i = min(n, j + 3)
            seg = z[lo_i:hi_i]
            m = lo_i + int(np.nanargmin(seg)) if np.isfinite(seg).any() else i
            r = profile[m]
            out.append({"index": m, "chainage_m": r["chainage_m"], "x": r["x"], "y": r["y"],
                        "z_dem_m": r.get("z_dem_m"), "sag_depth_m": float(depth),
                        "pond_depth_m": r.get("pond_depth_m")})
    return out


def ponding_sags(profile, existing, search_m=50.0, smooth_m=30.0, min_depth_m=0.3):
    """Sag points with no flow path: a low point within `search_m` of a stream
    crossing or an existing crossing is a drainage crossing, not a sag."""
    stream_ch = [r["chainage_m"] for r in profile if r.get("stream") == 1]
    ex_ch = [c for _, c in existing if c is not None]
    near = lambda c, lst: any(abs(c - x) <= float(search_m) for x in lst)
    return [s for s in sag_points(profile, smooth_m, min_depth_m)
            if not near(s["chainage_m"], stream_ch) and not near(s["chainage_m"], ex_ch)]


def flat_stretches(profile, raw, geotransform, alignment, flat_slope_pct=0.5,
                   crossfall_m=100.0, min_len_m=300.0, existing=()):
    """Flat / floodplain stretches (see the module docstring)."""
    ch = np.array([r["chainage_m"] for r in profile], float)
    x = np.array([r["x"] for r in profile], float)
    y = np.array([r["y"] for r in profile], float)
    z = np.array([np.nan if r.get("z_dem_m") is None else r["z_dem_m"] for r in profile], float)
    sl = np.array([np.nan if r.get("slope_long_pct") is None else r["slope_long_pct"]
                   for r in profile], float)
    _, _, seg = alignment.locate(x, y)
    tx, ty = alignment.tangent(seg)
    nx, ny = -ty, tx                                   # left normal
    zl = bilinear(np.asarray(raw, float), geotransform, x + nx * crossfall_m, y + ny * crossfall_m)
    zr = bilinear(np.asarray(raw, float), geotransform, x - nx * crossfall_m, y - ny * crossfall_m)
    cf = np.fmax(np.abs(zl - z), np.abs(zr - z)) / crossfall_m * 100.0
    flat = (np.abs(sl) < flat_slope_pct) & (cf < flat_slope_pct)
    out = []
    i = 0
    n = flat.size
    while i < n:
        if not flat[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and flat[j + 1]:
            j += 1
        c0, c1 = ch[i], ch[j]
        # boundaries half-way to the neighbouring stations
        a = c0 - (0.5 * (ch[i] - ch[i - 1]) if i > 0 else 0.0)
        b = c1 + (0.5 * (ch[j + 1] - ch[j]) if j + 1 < n else 0.0)
        if b - a >= float(min_len_m):
            within = [u for u, c in existing if c is not None and a <= float(c) <= b]
            out.append({"i0": i, "i1": j, "chainage_m": float(a), "chainage_to_m": float(b),
                        "length_m": float(b - a),
                        "slope_long_pct": float(np.nanmean(np.abs(sl[i:j + 1]))),
                        "crossfall_pct": float(np.nanmean(cf[i:j + 1])),
                        "n_crossings_within": len(within), "uids_within": ",".join(within)})
        i = j + 1
    return out


def alignment_wall(alignment, shape, geotransform, step=None, gaps=(), gap_cells=1):
    """Bool grid of the cells the centreline passes through (8-connected
    sampling every half cell), used as a wall for the sag areas.

    gaps : (x, y) points where the embankment is open - existing crossings and
        stream crossings. Wall cells within `gap_cells` of each are removed,
        so drainage that has (or will have) a culvert passes as before and only
        water with no way through ponds at the sags."""
    gt = geotransform
    rows, cols = shape
    cs = min(abs(gt[1]), abs(gt[5]))
    step = step or cs / 2.0
    wall = np.zeros(shape, bool)
    n = int(math.ceil(alignment.length / step)) + 1
    for k in range(n):
        c = alignment.start_chainage + min(k * step, alignment.length)
        px, py = alignment.point_at(c)
        cc = int(math.floor((px - gt[0]) / gt[1]))
        rr = int(math.floor((py - gt[3]) / gt[5]))
        if 0 <= rr < rows and 0 <= cc < cols:
            wall[rr, cc] = True
    # close diagonal gaps so flow cannot slip between two wall cells
    d = wall.copy()
    diag = (wall[:-1, :-1] & wall[1:, 1:] & ~wall[:-1, 1:] & ~wall[1:, :-1])
    d[:-1, 1:] |= diag
    anti = (wall[:-1, 1:] & wall[1:, :-1] & ~wall[:-1, :-1] & ~wall[1:, 1:])
    d[:-1, :-1] |= anti
    k = int(gap_cells)
    for gx, gy in gaps:
        cc = int(math.floor((gx - gt[0]) / gt[1]))
        rr = int(math.floor((gy - gt[3]) / gt[5]))
        d[max(rr - k, 0):rr + k + 1, max(cc - k, 0):cc + k + 1] = False
    return d


def walled_accumulation(filled, valid, alignment, geotransform, cell_width, cell_height,
                        wall_height_m=1000.0, flat_method="barnes", gaps=()):
    """Accumulation with the alignment as a wall.

    `filled` is the depression-filled DEM (spurious pits removed). The
    centreline cells are raised by `wall_height_m` (open at `gaps`: existing
    and stream crossings) and D8 is computed WITHOUT
    filling again, so water reaching the embankment runs along it to the
    lowest point and stops there: the accumulation at a sag is the local
    area draining to it on that side. Returns (acc_cells, wall, direction).
    """
    from ..flow.direction import d8_direction
    from ..flow.accumulation import flow_accumulation
    wall = alignment_wall(alignment, filled.shape, geotransform, gaps=gaps)
    z = np.asarray(filled, float)
    zw = np.where(wall, z + wall_height_m, z)
    v = valid & np.isfinite(zw)
    d, _ = d8_direction(np.nan_to_num(zw, nan=0.0), v, cell_width, cell_height,
                        resolve_flats=True, flat_method=flat_method)
    a, _ = flow_accumulation(d, v)
    return a, wall, d


def sag_areas(sags, acc_walled, wall, geotransform, alignment, offset_cells=1.5,
              direction=None):
    """Local area ponding at each sag on the left and right of the wall: the
    largest (accumulation + 1) x cell area within a 3x3 window centred
    `offset_cells` cells off the centreline on each side.

    With `direction` (the walled D8 grid) only cells where the water STOPS
    (no receiver: a pit against the embankment) count, so a river running
    alongside the road past the low point is not taken for ponding water."""
    from ..grid import NO_RECEIVER
    gt = geotransform
    rows, cols = acc_walled.shape
    cell_km2 = abs(gt[1] * gt[5]) / 1e6
    cs = math.sqrt(abs(gt[1] * gt[5]))
    for s in sags:
        _, _, seg = alignment.locate([s["x"]], [s["y"]])
        tx, ty = alignment.tangent(seg)
        nx, ny = -float(ty[0]), float(tx[0])
        res = {}
        for side, sgn in (("left", 1.0), ("right", -1.0)):
            px = s["x"] + sgn * nx * offset_cells * cs
            py = s["y"] + sgn * ny * offset_cells * cs
            c = int(math.floor((px - gt[0]) / gt[1]))
            r = int(math.floor((py - gt[3]) / gt[5]))
            best = None
            for rr in range(r - 2, r + 3):
                for cc in range(c - 2, c + 3):
                    if 0 <= rr < rows and 0 <= cc < cols and not wall[rr, cc] and \
                            (direction is None or direction[rr, cc] == NO_RECEIVER):
                        v = acc_walled[rr, cc]
                        if np.isfinite(v) and (best is None or v > best):
                            best = float(v)
            res[side] = (best + 1.0) * cell_km2 if best is not None else None
        s["sag_area_left_km2"], s["sag_area_right_km2"] = res["left"], res["right"]
        cand = [(v, k) for k, v in res.items() if v is not None]
        if cand:
            v, k = max(cand)
            s["sag_area_km2"], s["side"] = v, k
        else:
            s["sag_area_km2"], s["side"] = None, None
    return sags


# -- the whole check, as the design hydrology package runs it ----------------

def run_coverage(direction, valid, accumulation, elevation, filled, geotransform, alignment,
                 profile, crossings, build_kwargs=None, min_area_km2=None, merge_m=30.0,
                 search_m=50.0, small_area_km2=0.01, sag_smooth_m=30.0, sag_min_depth_m=0.3,
                 sag_min_area_km2=0.05, flat_slope_pct=0.5, crossfall_m=100.0,
                 flat_min_len_m=300.0, proposed_prefix="P", stream_threshold_cells=200.0,
                 progress=None, cluster_chainages=None):
    """Coverage check + sags + flat stretches, and the proposed crossings
    delineated and characterised exactly like the existing ones.

    crossings : the existing crossings already built ((xy, attrs) with
        outlet_uid, chainage_m, acc_at_outlet_km2)
    build_kwargs : keyword arguments passed on to build_exchange_records
        (soil, erosion, runoff, channel threshold ...)
    Returns dict(crossings, catchments, flowpaths, coverage, sags, flats,
    issues) - the first three hold only the PROPOSED features, tagged
    status = 'proposed'.
    cluster_chainages : [(uid, chainage_m)] from cluster_cover(): the other
        intersections of a clustered stream, covered by the crossing used for it.
    """
    from ..interop.heas_exchange import build_exchange_records, ExchangeError
    gt = tuple(geotransform)
    cw, chh = abs(gt[1]), abs(gt[5])
    if min_area_km2 is None:
        min_area_km2 = (float(stream_threshold_cells) + 1.0) * cw * chh / 1e6
    existing = [(a["outlet_uid"], a.get("chainage_m")) for _, a in crossings]
    issues, coverage = [], []

    # A2 missing crossings
    covered = existing + list(cluster_chainages or [])
    miss = missing_crossings(profile, covered, min_area_km2, merge_m, search_m)
    uncovered = [m for m in miss if not m["covered"]]
    for _, a in crossings:
        ar = a.get("acc_at_outlet_km2")
        if ar is not None and ar < float(small_area_km2):
            coverage.append(((a["outlet_x"], a["outlet_y"]), {
                "issue": "small_area", "chainage_m": a.get("chainage_m"), "chainage_to_m": None,
                "area_km2": ar, "uid": a["outlet_uid"], "nearest_uid": None, "nearest_m": None,
                "note": f"contributing area {ar:.4f} km2 below {small_area_km2:g} km2 - "
                        "redundant or misplaced? (review only)"}))

    # A3 sags (walled routing) and flat stretches
    sags = ponding_sags(profile, covered, search_m, sag_smooth_m, sag_min_depth_m)
    walled_dir = walled_acc = wall = None
    if sags:
        if progress is not None:
            progress(0.1, "Routing with the alignment as a wall (sag areas)")
        gaps = [(a["outlet_x"], a["outlet_y"]) for _, a in crossings] + \
            [(r["x"], r["y"]) for r in profile if r.get("stream") == 1]
        acc_w, wall, walled_dir = walled_accumulation(filled, valid, alignment, gt, cw, chh,
                                                      gaps=gaps)
        sag_areas(sags, acc_w, wall, gt, alignment, direction=walled_dir)
        walled_acc = acc_w
    flats = flat_stretches(profile, elevation, gt, alignment, flat_slope_pct, crossfall_m,
                           flat_min_len_m, existing)

    # proposed crossings: uncovered streams (normal routing) ...
    bk = dict(build_kwargs or {})
    bk.pop("id_prefix", None); bk.pop("id_order", None); bk.pop("id_scheme", None)
    from ..watershed.delineate import extract_streams
    smask = extract_streams(accumulation, valid, threshold_cells=float(stream_threshold_cells))
    out_cr, out_ca, out_fp = [], [], []
    taken = {u for u, _ in existing}
    outlets = []
    for k, m in enumerate(uncovered):
        outlets.append({"x": m["x"], "y": m["y"], "fid": -(k + 1), "source_id": None,
                        "chainage": m["chainage_m"], "_reason": "uncovered_stream",
                        "_nearest": (m["nearest_uid"], m["nearest_m"])})
    sag_outlets = []
    for k, s in enumerate(sags):
        if s.get("sag_area_km2") is not None and s["sag_area_km2"] >= float(sag_min_area_km2):
            # outlet = the cell with the largest walled accumulation on the
            # contributing side, next to the embankment
            _, _, seg = alignment.locate([s["x"]], [s["y"]])
            tx, ty = alignment.tangent(seg)
            sgn = 1.0 if s["side"] == "left" else -1.0
            px = s["x"] + sgn * (-float(ty[0])) * 1.5 * cw
            py = s["y"] + sgn * float(tx[0]) * 1.5 * chh
            c = int(math.floor((px - gt[0]) / gt[1])); r = int(math.floor((py - gt[3]) / gt[5]))
            from ..grid import NO_RECEIVER
            best = None
            for rr in range(r - 2, r + 3):
                for cc in range(c - 2, c + 3):
                    if 0 <= rr < walled_acc.shape[0] and 0 <= cc < walled_acc.shape[1] \
                            and not wall[rr, cc] and valid[rr, cc] \
                            and walled_dir[rr, cc] == NO_RECEIVER:
                        if best is None or walled_acc[rr, cc] > walled_acc[best]:
                            best = (rr, cc)
            if best is None:
                continue
            ox = gt[0] + (best[1] + 0.5) * gt[1]; oy = gt[3] + (best[0] + 0.5) * gt[5]
            uid, d = _nearest(s["chainage_m"], existing)
            sag_outlets.append({"x": ox, "y": oy, "outlet_x": ox, "outlet_y": oy,
                                "fid": -(1000 + k), "source_id": None,
                                "chainage": s["chainage_m"], "_reason": "sag_point",
                                "_nearest": (uid, d), "_sag": s})
    # dedupe outlet cells (two sags or a sag and a stream on one cell)
    seen = set()
    for grp in (outlets, sag_outlets):
        keep = []
        for o in grp:
            key = (round(o.get("outlet_x", o["x"]), 3), round(o.get("outlet_y", o["y"]), 3))
            if key not in seen:
                seen.add(key); keep.append(o)
        grp[:] = keep
    for n, o in enumerate(outlets + sag_outlets, start=1):
        o["_tmp"] = f"__tmp{n:05d}"
    built = []
    for group, dgrid, agrid, snap in ((outlets, direction, accumulation, 2),
                                      (sag_outlets, walled_dir, walled_acc, 0)):
        if not group:
            continue
        try:
            cr, ca, fp, iss, _ = build_exchange_records(
                dgrid, valid, agrid, elevation, gt,
                [dict(o, source_id=o["_tmp"]) for o in group], snap_radius_cells=snap,
                stream_mask=smask if snap else None, id_scheme="attribute", id_prefix="",
                **{k: v for k, v in bk.items() if k not in ("snap_radius_cells", "stream_mask")})
        except ExchangeError as e:
            issues.append(f"Proposed crossings ({group[0]['_reason']}): {e}")
            continue
        issues += iss
        by_tmp = {o["_tmp"]: o for o in group}
        for c in cr:
            built.append(by_tmp[c[1]["outlet_uid"]])
        built_layers = (cr, ca, fp)
        for lst, dest in zip(built_layers, (out_cr, out_ca, out_fp)):
            for geom, a in lst:
                o = by_tmp[a["outlet_uid"]]
                a.update(status="proposed", proposed_reason=o["_reason"],
                         nearest_uid=o["_nearest"][0], nearest_m=o["_nearest"][1])
                dest.append((geom, a))
    # one gapless numbering along the chainage for the proposed crossings built
    built.sort(key=lambda o: (o["chainage"], o["_tmp"]))
    rename = {}
    for n, o in enumerate(built, start=1):
        uid = f"{proposed_prefix}{n:03d}"
        if uid in taken:
            raise ExchangeError(f"Proposed ID {uid} already used by an existing crossing; "
                                "choose another proposed prefix.")
        rename[o["_tmp"]] = uid
        o["_uid"] = uid
        if o["_reason"] == "sag_point":
            o["_sag"]["uid"] = uid
    for lst in (out_cr, out_ca, out_fp):
        for _, a in lst:
            u = rename[a["outlet_uid"]]
            a["outlet_uid"] = u
            for k in ("crossing_id", "catchment_id", "flowpath_id", "source_id"):
                if k in a:
                    a[k] = u
    for lst in (out_cr, out_ca, out_fp):
        lst.sort(key=lambda t: t[1]["outlet_uid"])
    for m in miss:
        uid = next((o.get("_uid") for o in outlets if o["chainage"] == m["chainage_m"]), None)
        if m["covered"]:
            continue
        coverage.append(((m["x"], m["y"]), {
            "issue": "missing_crossing", "chainage_m": m["chainage_m"], "chainage_to_m": None,
            "area_km2": m["area_km2"], "uid": uid, "nearest_uid": m["nearest_uid"],
            "nearest_m": m["nearest_m"],
            "note": f"stream of {m['area_km2']:.3g} km2 crosses the road with no crossing within "
                    f"{search_m:g} m" + (f" - proposed {uid}" if uid else "")}))
    for s in sags:
        coverage.append(((s["x"], s["y"]), {
            "issue": "sag_point", "chainage_m": s["chainage_m"], "chainage_to_m": None,
            "area_km2": s.get("sag_area_km2"), "uid": s.get("uid"),
            "nearest_uid": _nearest(s["chainage_m"], existing)[0],
            "nearest_m": _nearest(s["chainage_m"], existing)[1],
            "note": (f"low point {s['sag_depth_m']:.2f} m deep; water ponds against the embankment "
                     f"({s['side']} side, {s['sag_area_km2']:.3g} km2)" if s.get("sag_area_km2") else
                     f"low point {s['sag_depth_m']:.2f} m deep in the road profile; no water ponds "
                     "against the embankment here (the ground drains past) - check the road "
                     "long section")
                    + (f" - proposed {s['uid']}" if s.get("uid") else "")}))
    for f in flats:
        mid = alignment.point_at(0.5 * (f["chainage_m"] + f["chainage_to_m"]))
        coverage.append((mid, {
            "issue": "flat_stretch", "chainage_m": f["chainage_m"], "chainage_to_m": f["chainage_to_m"],
            "area_km2": None, "uid": None, "nearest_uid": None, "nearest_m": None,
            "note": f"{f['length_m']:.0f} m flat (long {f['slope_long_pct']:.2f} %, cross-fall "
                    f"{f['crossfall_pct']:.2f} %), {f['n_crossings_within']} crossing(s) within: "
                    "relief culverts for sheet flow may be needed (not placed or sized by QEHT)"}))
    coverage.sort(key=lambda t: (t[1]["chainage_m"] is None, t[1]["chainage_m"] or 0))
    return {"crossings": out_cr, "catchments": out_ca, "flowpaths": out_fp,
            "coverage": coverage, "sags": sags, "flats": flats, "issues": issues,
            "n_missing": len(uncovered)}
