# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""DEM-error sensitivity per crossing (F15, v0.27).

Turns each deterministic catchment value into a range. For each of N
realisations a spatially correlated Gaussian error field is added to the
raw DEM, which is then filled, routed (D8, Barnes) and accumulated; every
crossing's fixed pour point is snapped to the stream within the snap
radius (or recorded as lost) and its area, longest flow path length,
10-85 slope and Kirpich Tc are measured.

Error field: white noise N(0, 1) convolved (FFT, on a grid padded by three
correlation lengths so the periodic wrap does not correlate opposite
edges) with the kernel k(r) = exp(-2 r^2 / L^2), scaled by 1 / sqrt(sum k^2)
so its variance is 1, then multiplied by sigma. Its correlation is
rho(r) = exp(-r^2 / L^2): 1/e at the correlation length L.

Per crossing: p10, p50, p90 and the coefficient of variation of each
value over the realisations where the outlet kept a stream; unc_lost_pct
(the outlet lost its stream) and unc_switch_pct (the area moved by more
than 25 % from the deterministic run: the catchment switched).

sigma presets (vertical error, one standard deviation, metres):
  AW3D30  4.4  Tadono et al. (2016), JAXA validation of AW3D30: 4.40 m RMSE
               at 5,121 check points worldwide.
  FABDEM  2.5  ESTIMATE from Hawker et al. (2022), Environ. Res. Lett. 17
               024016: mean absolute error 1.12 m (built-up) to 2.88 m
               (forest) against LiDAR / ICESat, i.e. about 1.4-3.6 m RMSE
               for Gaussian errors (RMSE = 1.25 MAE); error varies by land
               cover - choose 'custom' for a site-specific value.
There is no silent default: the user picks a preset or a value.

No QGIS imports, no GDAL.
"""

import json
import math
import time

import numpy as np

PRESETS = {
    "AW3D30": (4.4, "Tadono et al. (2016), AW3D30 validation (JAXA): 4.40 m RMSE at 5,121 "
                    "independent check points"),
    "FABDEM": (2.5, "ESTIMATE from Hawker et al. (2022), Environ. Res. Lett. 17 024016: MAE "
                    "1.12 m (built-up) to 2.88 m (forest) vs LiDAR / ICESat, about 1.4-3.6 m "
                    "RMSE; mid-range value"),
}
DEFAULT_CORR_M = 90.0
DEFAULT_N = 50
DEFAULT_SEED = 2026
QUANT = [("area", "unc_area"), ("lfp_length_m", "unc_lfp"), ("lfp_slope_1085", "unc_s1085"),
         ("tc_kirpich_min", "unc_tc")]
UNC_FIELDS = [f"{p}_{s}" for _, p in QUANT for s in ("p10", "p50", "p90", "cv")] + \
    ["unc_lost_pct", "unc_switch_pct", "unc_n", "unc_note"]


def correlated_field(shape, sigma, corr_len_m, cell, rng):
    """Gaussian field with standard deviation sigma and correlation exp(-r^2 / L^2)."""
    rows, cols = shape
    if sigma == 0:
        return np.zeros(shape)
    L = max(float(corr_len_m) / float(cell), 1e-6)        # in cells
    pad = int(math.ceil(3 * L)) + 1
    R, C = rows + 2 * pad, cols + 2 * pad
    noise = rng.standard_normal((R, C))
    dy = np.minimum(np.arange(R), R - np.arange(R))[:, None]
    dx = np.minimum(np.arange(C), C - np.arange(C))[None, :]
    k = np.exp(-2.0 * (dx * dx + dy * dy) / (L * L))
    f = np.fft.irfft2(np.fft.rfft2(noise) * np.fft.rfft2(k), s=(R, C))
    f /= math.sqrt(float((k * k).sum()))
    return sigma * f[pad:pad + rows, pad:pad + cols]


def upstream_lengths(direction, valid, cw, ch):
    """Longest flow length arriving at every cell, in one topological pass.

    For a full (overlapping) catchment the longest flow path to an outlet
    is the longest path arriving at that cell, so one pass serves every
    outlet. Returns (length grid flattened, from_cell: the upstream
    neighbour on that path, -1 at a divide)."""
    from ..grid import receivers_from_direction, neighbour_distances
    shape = direction.shape
    rec = receivers_from_direction(np.where(valid, direction, -1), shape)
    dist = neighbour_distances(cw, ch)
    dflat = np.asarray(direction).ravel()
    n = rec.size
    live = np.asarray(valid).ravel() & (rec >= 0)
    indeg = np.bincount(rec[live], minlength=n)
    length = np.zeros(n)
    from_cell = np.full(n, -1, dtype=np.int64)
    front = np.flatnonzero(np.asarray(valid).ravel() & (indeg == 0))
    while front.size:
        f = front[live[front]]
        if f.size == 0:
            break
        tgt = rec[f]
        cand = length[f] + dist[dflat[f]]
        order = np.lexsort((f, cand))            # ascending: the last write per target wins
        tgt_o, cand_o, f_o = tgt[order], cand[order], f[order]
        better = cand_o > length[tgt_o]
        length[tgt_o[better]] = cand_o[better]
        from_cell[tgt_o[better]] = f_o[better]
        np.subtract.at(indeg, tgt, 1)
        nxt = np.unique(tgt)
        front = nxt[indeg[nxt] == 0]
    return length, from_cell


def measure(dem, valid, cw, ch, outlets, threshold_cells, snap_cells=5):
    """Route dem and measure each outlet -> list of dicts (None = lost)."""
    from ..conditioning.fill import fill_depressions
    from ..flow.direction import d8_direction
    from ..flow.accumulation import flow_accumulation
    from .delineate import extract_streams, snap_pour_point
    from .statistics import path_10_85
    from ..runoff.tc import kirpich_min
    z = np.where(valid, dem, 0.0)
    filled, _, _ = fill_depressions(z, valid, cell_width=cw, cell_height=ch)
    d, _ = d8_direction(filled, valid, cw, ch)
    acc, _ = flow_accumulation(d, valid)
    streams = extract_streams(acc, valid, threshold_cells=threshold_cells)
    length, from_cell = upstream_lengths(d, valid, cw, ch)
    cols = dem.shape[1]
    out = []
    for r0, c0 in outlets:
        r, c, _, _ = snap_pour_point(r0, c0, acc, valid, search_radius_cells=snap_cells,
                                     stream_mask=streams)
        if not streams[r, c]:
            out.append(None)
            continue
        path, k = [], r * cols + c
        while k >= 0:
            path.append(divmod(int(k), cols))
            k = from_cell[k]
        cells = path[::-1]                                    # divide -> outlet
        p = path_10_85(cells, dem, cw, ch)
        L = float(length[r * cols + c])
        S = p["slope"] if math.isfinite(p["slope"]) else None
        out.append({"area": (float(acc[r, c]) + 1.0) * cw * ch / 1e6, "lfp_length_m": L,
                    "lfp_slope_1085": S,
                    "tc_kirpich_min": kirpich_min(L, S) if L and S and S > 0 else None,
                    "cell": (r, c)})
    return out


def _stats(v):
    v = np.array([x for x in v if x is not None and math.isfinite(x)], float)
    if v.size == 0:
        return None, None, None, None
    m = float(v.mean())
    return (float(np.percentile(v, 10)), float(np.percentile(v, 50)), float(np.percentile(v, 90)),
            float(v.std() / m) if m != 0 else None)


def run(dem, valid, gt, outlets, sigma, corr_len_m=DEFAULT_CORR_M, n=DEFAULT_N, seed=DEFAULT_SEED,
        threshold_cells=200.0, snap_cells=5, switch_frac=0.25, progress=None, log=None):
    """Monte Carlo over n realisations -> (list of UNC field dicts per outlet, info dict)."""
    cw, ch = abs(gt[1]), abs(gt[5])
    t0 = time.time()
    base = measure(dem, valid, cw, ch, outlets, threshold_cells, snap_cells)
    per_run = time.time() - t0
    if log is not None:
        log(f"DEM uncertainty: {n} realisations at about {per_run:.1f} s each "
            f"(about {per_run * n / 60.0:.1f} min).")
    rng = np.random.default_rng(seed)
    runs = []
    for i in range(n):
        if progress is not None:
            progress(i / max(n, 1))
        e = correlated_field(dem.shape, sigma, corr_len_m, min(cw, ch), rng)
        runs.append(measure(np.where(valid, dem + e, dem), valid, cw, ch, outlets, threshold_cells,
                            snap_cells))
    out = []
    for j in range(len(outlets)):
        b = {k: None for k in UNC_FIELDS}
        vals = [r[j] for r in runs]
        kept = [v for v in vals if v is not None]
        b["unc_n"] = n
        b["unc_lost_pct"] = 100.0 * (n - len(kept)) / n if n else None
        for key, pre in QUANT:
            p10, p50, p90, cv = _stats([v[key] for v in kept])
            b[f"{pre}_p10"], b[f"{pre}_p50"], b[f"{pre}_p90"], b[f"{pre}_cv"] = p10, p50, p90, cv
        a0 = base[j]["area"] if base[j] is not None else None
        if a0 and kept:
            b["unc_switch_pct"] = 100.0 * sum(1 for v in kept if abs(v["area"] - a0) > switch_frac * a0) / n
        notes = []
        if base[j] is None:
            notes.append("no stream within the snap radius on the deterministic DEM")
        if b["unc_lost_pct"]:
            notes.append(f"outlet lost its stream in {b['unc_lost_pct']:.0f} % of realisations")
        if b["unc_switch_pct"]:
            notes.append(f"area moved > {100 * switch_frac:.0f} % in {b['unc_switch_pct']:.0f} % "
                         "(catchment switching)")
        b["unc_note"] = "; ".join(notes) or None
        out.append(b)
    info = {"n": n, "sigma_m": sigma, "corr_len_m": corr_len_m, "seed": seed,
            "threshold_cells": threshold_cells, "snap_cells": snap_cells,
            "switch_frac": switch_frac, "seconds_per_realisation": per_run,
            "deterministic": [None if x is None else {k: x[k] for k in ("area", "lfp_length_m",
                                                                         "lfp_slope_1085",
                                                                         "tc_kirpich_min")}
                              for x in base]}
    return out, info


def corr_sensitivity(dem, valid, gt, outlets, sigma, corr_len_m, n, seed, threshold_cells,
                     snap_cells=5, factors=(0.5, 2.0)):
    """Area p50 and CV at corr_len x factor for the given outlets (the three largest)."""
    res = {}
    for f in factors:
        o, _ = run(dem, valid, gt, outlets, sigma, corr_len_m * f, n, seed, threshold_cells,
                   snap_cells)
        res[f"{f:g}x"] = [{"area_p50": b["unc_area_p50"], "area_cv": b["unc_area_cv"]} for b in o]
    return res


def params_json(info, preset, source, sensitivity=None, uids=None):
    d = {k: v for k, v in info.items() if k != "deterministic"}
    d.update({"preset": preset, "sigma_source": source,
              "field": "Gaussian, correlation exp(-r^2/L^2) (1/e at the correlation length), "
                       "FFT convolution of white noise",
              "routing": "fill, D8 (Barnes), accumulation per realisation; pour point snapped "
                         "to the stream within the snap radius",
              "switch": "area differs from the deterministic run by more than switch_frac"})
    if sensitivity is not None:
        d["corr_len_sensitivity"] = {"outlets": uids, **sensitivity}
    return json.dumps(d, sort_keys=True, default=float)
