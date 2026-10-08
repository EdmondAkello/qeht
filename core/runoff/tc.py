# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Time of concentration by five published methods, side by side (F10, v0.22).

QEHT reports every method with its inputs and a validity flag and does not
pick one: the choice is the engineer's. A method missing a user input is
NULL with a note, never filled with an invented value.

Kirpich (1940)          Tc [min] = 0.0195 L^0.77 S^-0.385
                        L flow path (m), S 10-85 slope (m/m)
Kerby (1959) + Kirpich  t_ov [min] = 1.44 (L_ov N)^0.467 S_ov^-0.235 (L_ov in m,
                        N retardance), plus Kirpich on the channel part
SCS / NRCS lag          lag [h] = l^0.8 (S + 1)^0.7 / (1900 Y^0.5), l in ft,
(NEH 630 ch. 15)        S = 1000 / CN - 10 (in), Y average watershed slope (%);
                        Tc = lag / 0.6
TR-55 segmental (1986)  sheet  Tt [h] = 0.007 (n L)^0.8 / (P2^0.5 s^0.4), L in ft
                        (<= 100 ft = 30 m), P2 the 2-yr 24-h rainfall (in);
                        shallow V = 16.1345 s^0.5 ft/s (unpaved);
                        channel V = (1 / n) R^(2/3) s^(1/2) (Manning, SI) on the
                        bank-full section of F5
Bransby-Williams        Tc [min] = 14.6 L / (A^0.1 S^0.2), L km, A km2, S m/km

Validity ranges (flags only; values are always reported) and their
sources are in RANGES and in the run metadata tc_params_json.

No QGIS imports, no GDAL.
"""

import json
import math

FT = 0.3048
IN = 25.4

METHODS = ("kirpich", "kerby_kirpich", "scs_lag", "tr55", "bransby_williams")
TC_FIELDS = [f"tc_{m}_min" for m in METHODS] + [f"tc_{m}_flag" for m in METHODS] + \
    ["tc_basis_json", "tc_note"]

# WorldCover class -> TR-55 (1986) Table 3-1 sheet-flow Manning n (PROXY)
SHEET_N = {10: 0.40, 20: 0.24, 30: 0.15, 40: 0.06, 50: 0.011, 60: 0.011,
           70: None, 80: None, 90: None, 95: None, 100: None}
SHEET_N_ID = ("PROXY: TR-55 (1986) Table 3-1 matched to WorldCover - tree 0.40 (woods, light "
              "underbrush), shrub 0.24 (dense grasses), grass 0.15 (short grass prairie), crop "
              "0.06 (cultivated, residue <= 20 %), built-up and bare 0.011 (smooth surfaces); "
              "water, wetland, snow, moss none")
# WorldCover class -> Kerby (1959) retardance N (PROXY)
KERBY_N = {10: 0.60, 20: 0.40, 30: 0.40, 40: 0.20, 50: 0.02, 60: 0.10,
           70: None, 80: None, 90: None, 95: None, 100: None}
KERBY_N_ID = ("PROXY: Kerby (1959) matched to WorldCover - built-up 0.02 (smooth impervious), "
              "bare 0.10 (smooth bare packed soil), crop 0.20 (poor grass / cultivated rows), "
              "grass and shrub 0.40 (pasture / average grass), tree 0.60 (deciduous forest); "
              "water, wetland, snow, moss none")
CHANNEL_N = 0.035                     # clean natural channel (Chow 1959); editable default
SHEET_MAX_M = 30.0                    # 100 ft (NEH 630 ch. 15, 2010)

RANGES = {
    "kirpich": {"area_km2": [0.004, 0.45], "slope_pct": [3.0, 10.0],
                "source": "Kirpich (1940), Tennessee data"},
    "scs_lag": {"area_km2_max": 8.0, "cn": [50.0, 95.0],
                "source": "NEH 630 ch. 15 (about 2,000 acres)"},
    "kerby_kirpich": {"overland_m_max": 365.0, "source": "Kerby (1959)"},
    "tr55": {"sheet_m_max": SHEET_MAX_M,
             "source": "TR-55 (1986); sheet flow <= 100 ft per NEH 630 ch. 15 (2010)"},
    "bransby_williams": {"note": "rural catchments", "source": "no hard limit"},
}


def kirpich_min(L_m, S):
    return 0.0195 * L_m ** 0.77 * S ** -0.385


def kerby_min(L_m, N, S):
    return 1.44 * (L_m * N) ** 0.467 * S ** -0.235


def scs_lag_min(L_m, cn, Y_pct):
    s = 1000.0 / cn - 10.0
    lag_h = (L_m / FT) ** 0.8 * (s + 1.0) ** 0.7 / (1900.0 * Y_pct ** 0.5)
    return 60.0 * lag_h / 0.6


def tr55_sheet_h(n, L_m, p2_mm, s):
    return 0.007 * (n * L_m / FT) ** 0.8 / ((p2_mm / IN) ** 0.5 * s ** 0.4)


def tr55_shallow_h(L_m, s, paved=False):
    v_fts = (20.3282 if paved else 16.1345) * s ** 0.5
    return (L_m / FT) / v_fts / 3600.0


def manning_v(area, perimeter, s, n):
    r = area / perimeter
    return r ** (2.0 / 3.0) * s ** 0.5 / n


def bransby_williams_min(L_km, A_km2, S_m_per_km):
    return 14.6 * L_km / (A_km2 ** 0.1 * S_m_per_km ** 0.2)


def section_hydraulics(stations, depth):
    """Flow area and wetted perimeter of the contiguous section below
    z_bed + depth around the lowest station ([[s, z], ...] from F5)."""
    if not stations or depth is None or depth <= 0:
        return None, None
    s = [float(a) for a, _ in stations]
    z = [float(b) for _, b in stations]
    i0 = min(range(len(z)), key=lambda i: (z[i], abs(s[i])))
    level = z[i0] + depth
    lo = i0
    while lo > 0 and z[lo - 1] <= level:
        lo -= 1
    hi = i0
    while hi < len(z) - 1 and z[hi + 1] <= level:
        hi += 1
    area = perim = 0.0
    for i in range(max(lo - 1, 0), min(hi + 1, len(z) - 1)):
        x0, z0, x1, z1 = s[i], z[i], s[i + 1], z[i + 1]
        d0, d1 = level - z0, level - z1
        if d0 <= 0 and d1 <= 0:
            continue
        if d0 < 0 or d1 < 0:              # partly wet: cut at the water line
            t = d0 / (d0 - d1)
            xw = x0 + t * (x1 - x0)
            if d0 < 0:
                x0, z0, d0 = xw, level, 0.0
            else:
                x1, z1, d1 = xw, level, 0.0
        area += 0.5 * (d0 + d1) * abs(x1 - x0)
        perim += math.hypot(x1 - x0, z1 - z0)
    return (area, perim) if area > 0 and perim > 0 else (None, None)


def _ok(v):
    return v is not None and isinstance(v, (int, float)) and math.isfinite(v)


def _range_flag(checks):
    out = [msg for bad, msg in checks if bad]
    return "outside: " + "; ".join(out) if out else "within"


def tc_block(fp, ca, xs=None, sheet_n=None, kerby_n=None, p2_mm=None, channel_n=CHANNEL_N,
             cn_field="cn_ii"):
    """TC fields for one crossing.

    fp: flowpath attrs (lfp_*), ca: catchment attrs (area_km2, catch_slope_horn,
    cn_ii / cn_export), xs: crossing attrs with the F5 section (optional).
    sheet_n / kerby_n: area-weighted values along the sheet / overland part
    (None when no land cover). p2_mm: 2-yr 24-h rainfall (mm) or None.
    """
    out = {k: None for k in TC_FIELDS}
    notes, basis = [], {}
    L = fp.get("lfp_length_m")
    S1085 = fp.get("lfp_slope_1085")
    A = ca.get("area_km2")
    # Kirpich
    if _ok(L) and _ok(S1085) and S1085 > 0:
        out["tc_kirpich_min"] = kirpich_min(L, S1085)
        r = RANGES["kirpich"]
        out["tc_kirpich_flag"] = _range_flag([
            (_ok(A) and not r["area_km2"][0] <= A <= r["area_km2"][1],
             f"area {A:.3g} km2 not in {r['area_km2'][0]}-{r['area_km2'][1]}"),
            (not r["slope_pct"][0] <= 100 * S1085 <= r["slope_pct"][1],
             f"slope {100 * S1085:.2g} % not in {r['slope_pct'][0]:g}-{r['slope_pct'][1]:g}")])
        basis["kirpich"] = {"L_m": L, "S_m_m": S1085, "S_basis": "lfp_slope_1085"}
    else:
        notes.append("Kirpich: no flow-path length or 10-85 slope")
    # Kerby overland + Kirpich channel
    Lov, Sov = fp.get("lfp_overland_m"), fp.get("lfp_overland_slope")
    Lch, Sch = fp.get("lfp_channel_m") or 0.0, fp.get("lfp_channel_slope")
    if not _ok(kerby_n):
        notes.append("Kerby: no retardance N (give a land-cover raster or an N lookup)")
    elif not (_ok(Lov) and _ok(Sov) and Sov > 0):
        notes.append("Kerby: no overland length or slope")
    elif Lch > 0 and not (_ok(Sch) and Sch > 0):
        notes.append("Kerby: channel part without a slope")
    else:
        t_ov = kerby_min(Lov, kerby_n, Sov) if Lov > 0 else 0.0
        t_ch = kirpich_min(Lch, Sch) if Lch > 0 else 0.0
        out["tc_kerby_kirpich_min"] = t_ov + t_ch
        out["tc_kerby_kirpich_flag"] = _range_flag([
            (Lov > RANGES["kerby_kirpich"]["overland_m_max"],
             f"overland {Lov:.0f} m > {RANGES['kerby_kirpich']['overland_m_max']:g} m")])
        basis["kerby_kirpich"] = {"L_ov_m": Lov, "S_ov": Sov, "N": kerby_n, "t_ov_min": t_ov,
                                  "L_ch_m": Lch, "S_ch": Sch, "t_ch_min": t_ch}
    # SCS lag
    cn = ca.get(cn_field)
    Y = ca.get("catch_slope_horn")
    if not _ok(cn):
        notes.append(f"SCS lag: no curve number ({cn_field})")
    elif not (_ok(L) and _ok(Y) and Y > 0):
        notes.append("SCS lag: no flow-path length or catchment slope")
    else:
        out["tc_scs_lag_min"] = scs_lag_min(L, cn, 100.0 * Y)
        r = RANGES["scs_lag"]
        out["tc_scs_lag_flag"] = _range_flag([
            (_ok(A) and A > r["area_km2_max"], f"area {A:.3g} km2 > {r['area_km2_max']:g}"),
            (not r["cn"][0] <= cn <= r["cn"][1], f"CN {cn:.0f} not in {r['cn'][0]:.0f}-{r['cn'][1]:.0f}")])
        basis["scs_lag"] = {"L_m": L, "CN": cn, "CN_field": cn_field, "Y_pct": 100.0 * Y,
                            "Y_basis": "catch_slope_horn (Horn mean)"}
    # TR-55 segmental
    tr = []
    if not _ok(p2_mm):
        tr.append("give P2")
    if not _ok(sheet_n) and _ok(Lov) and Lov > 0:
        tr.append("no sheet-flow n (land cover)")
    if not (_ok(Lov) and _ok(Sov) and Sov > 0) and _ok(Lov) and Lov > 0:
        tr.append("no overland slope")
    v_ch, area, perim = None, None, None
    if Lch > 0:
        if xs is not None and xs.get("xs_station_elev_json") and _ok(xs.get("xs_bankfull_d_m")):
            area, perim = section_hydraulics(json.loads(xs["xs_station_elev_json"]),
                                             xs["xs_bankfull_d_m"])
        if area is None:
            tr.append("no bank-full channel section (F5)")
        elif not (_ok(Sch) and Sch > 0):
            tr.append("no channel slope")
        else:
            v_ch = manning_v(area, perim, Sch, channel_n)
    if tr:
        notes.append("TR-55: " + ", ".join(tr))
    else:
        Lov_ = Lov if _ok(Lov) else 0.0
        Ls = min(Lov_, SHEET_MAX_M)
        Lsh = Lov_ - Ls
        t_s = tr55_sheet_h(sheet_n, Ls, p2_mm, Sov) if Ls > 0 else 0.0
        t_sh = tr55_shallow_h(Lsh, Sov) if Lsh > 0 else 0.0
        t_c = (Lch / v_ch / 3600.0) if Lch > 0 else 0.0
        out["tc_tr55_min"] = 60.0 * (t_s + t_sh + t_c)
        out["tc_tr55_flag"] = "within" if xs is None or xs.get("xs_quality") != "low" else \
            "outside: channel section quality low (velocity indicative)"
        basis["tr55"] = {"P2_mm": p2_mm, "sheet_m": Ls, "sheet_n": sheet_n, "s_overland": Sov,
                         "shallow_m": Lsh, "shallow_surface": "unpaved", "channel_m": Lch,
                         "channel_n": channel_n, "channel_s": Sch,
                         "channel_area_m2": area, "channel_perimeter_m": perim,
                         "channel_v_ms": v_ch, "t_sheet_min": 60 * t_s,
                         "t_shallow_min": 60 * t_sh, "t_channel_min": 60 * t_c}
    # Bransby-Williams
    if _ok(L) and _ok(A) and A > 0 and _ok(S1085) and S1085 > 0:
        out["tc_bransby_williams_min"] = bransby_williams_min(L / 1000.0, A, 1000.0 * S1085)
        out["tc_bransby_williams_flag"] = "rural catchments"
        basis["bransby_williams"] = {"L_km": L / 1000.0, "A_km2": A, "S_m_km": 1000.0 * S1085}
    else:
        notes.append("Bransby-Williams: no length, area or slope")
    out["tc_basis_json"] = json.dumps(basis, sort_keys=True, default=float)
    out["tc_note"] = "; ".join(notes) or None
    return out


def along_path_value(coords, length_m, classes_at, lookup):
    """Length-weighted lookup value over the first length_m of a path
    (coords from the divide). classes_at(xs, ys) -> class codes. None
    when no part of it has a value."""
    if not coords or length_m is None or length_m <= 0:
        return None
    tot = wsum = 0.0
    d = 0.0
    segs = []
    for (x0, y0), (x1, y1) in zip(coords[:-1], coords[1:]):
        seg = math.hypot(x1 - x0, y1 - y0)
        if d >= length_m or seg <= 0:
            d += seg
            continue
        use = min(seg, length_m - d)
        segs.append((0.5 * (x0 + x0 + (x1 - x0) * use / seg),
                     0.5 * (y0 + y0 + (y1 - y0) * use / seg), use))
        d += seg
    if not segs:
        return None
    cls = classes_at([s[0] for s in segs], [s[1] for s in segs])
    for (_, _, w), c in zip(segs, cls):
        try:
            v = lookup.get(int(c)) if c is not None and math.isfinite(float(c)) else None
        except (TypeError, ValueError):
            v = None
        if v is None:
            continue
        tot += w
        wsum += w * v
    return wsum / tot if tot > 0 else None


def read_value_lookup(path):
    """CSV class,value -> {class: value} (a user n or N lookup)."""
    from .curve_number import read_lookup_csv
    return {k: v[0] for k, v in read_lookup_csv(path, "value").items()}


def params_json(p2_mm, p2_source, channel_n, sheet_id, kerby_id, cn_field):
    return json.dumps({"methods": list(METHODS), "ranges": RANGES,
                       "p2_mm": p2_mm, "p2_source": p2_source or "not given (TR-55 NULL)",
                       "channel_n": channel_n,
                       "channel_n_source": "default 0.035, clean natural channel (Chow 1959)"
                       if channel_n == CHANNEL_N else "user value",
                       "sheet_n_lookup": sheet_id, "kerby_n_lookup": kerby_id,
                       "scs_cn": cn_field, "shallow_surface": "unpaved",
                       "sheet_max_m": SHEET_MAX_M,
                       "note": "every method reported; QEHT does not choose one"}, sort_keys=True)
