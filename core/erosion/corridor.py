# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Erosion indices sampled along a road alignment (WP-F Level 2, v0.14).

Stations every `step` metres of chainage (default 10 m). At each station
the grids are sampled on the centreline and on perpendicular offsets out
to `half_width` on each side (every `offset_step` metres, nearest cell).
Per station and side: max and mean of ln(SPI), LS, RUSLE A and TWI; the
combined severity score on the centreline and the worst score within the
buffer (either side). Consecutive stations with the same worst score form
a reach; reach boundaries sit half-way between stations, so reach lengths
add up exactly to the sampled alignment length.
"""

import numpy as np

GRIDS = ("ln_spi", "ls", "soil_loss", "twi")


def _stations(alignment, step):
    total = alignment.length
    n = int(np.floor(total / step)) + 1
    chs = alignment.start_chainage + np.arange(n) * step
    if total - (n - 1) * step > 1e-6:
        chs = np.append(chs, alignment.end_chainage)
    return chs


def _sample(grid, gt, xs, ys):
    rows, cols = grid.shape
    c = np.floor((xs - gt[0]) / gt[1]).astype(np.int64)
    r = np.floor((ys - gt[3]) / gt[5]).astype(np.int64)
    ok = (r >= 0) & (r < rows) & (c >= 0) & (c < cols)
    out = np.full(xs.shape, np.nan)
    out[ok] = grid[r[ok], c[ok]]
    return out


def sample_corridor(alignment, grids, combined_score, geotransform, step=10.0,
                    half_width=50.0, offset_step=10.0):
    """Returns (stations, reaches).

    grids : dict name -> 2-D array for any of GRIDS (missing ones skipped)
    combined_score : uint8 grid of combined severity scores (0 = none)
    stations : list of dicts (chainage, x, y, centre_score, worst_score,
               worst_side, and <grid>_<L|R>_max / _mean)
    reaches : list of dicts (ch_start, ch_end, length_m, worst_score,
              n_stations, ln_spi_max, soil_loss_max)
    """
    gt = tuple(geotransform)
    chs = _stations(alignment, float(step))
    k = np.clip(np.searchsorted(alignment.ch0, chs, side="right") - 1, 0, len(alignment.ch0) - 1)
    u = np.clip((chs - alignment.ch0[k]) / alignment.seg_len[k], 0.0, 1.0)
    x = alignment.x0[k] + u * alignment.dx[k]
    y = alignment.y0[k] + u * alignment.dy[k]
    tx, ty = alignment.tangent(k)
    nx, ny = -ty, tx                                   # left normal
    offs = np.arange(offset_step, half_width + 1e-9, offset_step) if half_width > 0 else np.array([])
    present = {g: grids[g] for g in GRIDS if grids.get(g) is not None}
    score = np.asarray(combined_score, dtype=np.float64)
    centre = _sample(score, gt, x, y)
    side_vals = {}
    for side, sgn in (("L", 1.0), ("R", -1.0)):
        o = np.concatenate([[0.0], offs]) * sgn
        X = x[:, None] + o[None, :] * nx[:, None]
        Y = y[:, None] + o[None, :] * ny[:, None]
        sc = _sample(score, gt, X, Y)
        side_vals[side] = {"score": np.nanmax(np.where(np.isfinite(sc), sc, -1), axis=1)}
        for g, grid in present.items():
            v = _sample(grid, gt, X, Y)
            with np.errstate(all="ignore"):
                import warnings
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    side_vals[side][g + "_max"] = np.nanmax(v, axis=1)
                    side_vals[side][g + "_mean"] = np.nanmean(v, axis=1)
    stations = []
    for i in range(chs.size):
        sl, sr = side_vals["L"]["score"][i], side_vals["R"]["score"][i]
        worst = max(sl, sr)
        st = {"chainage": float(chs[i]), "x": float(x[i]), "y": float(y[i]),
              "centre_score": int(centre[i]) if np.isfinite(centre[i]) and centre[i] > 0 else None,
              "worst_score": int(worst) if worst > 0 else None,
              "worst_side": ("L" if sl > sr else "R" if sr > sl else "both") if worst > 0 else None}
        for side in ("L", "R"):
            for g in present:
                for stat in ("max", "mean"):
                    v = side_vals[side][f"{g}_{stat}"][i]
                    st[f"{g}_{side}_{stat}"] = float(v) if np.isfinite(v) else None
        stations.append(st)
    # reaches
    reaches = []
    if stations:
        bounds = [chs[0]] + [0.5 * (chs[i] + chs[i + 1]) for i in range(chs.size - 1)] + [chs[-1]]
        start = 0
        for i in range(1, chs.size + 1):
            if i == chs.size or stations[i]["worst_score"] != stations[start]["worst_score"]:
                seg = stations[start:i]
                def _mx(key):
                    vals = [s.get(key) for s in seg if s.get(key) is not None]
                    return max(vals) if vals else None
                spi = [v for v in (_mx("ln_spi_L_max"), _mx("ln_spi_R_max")) if v is not None]
                sl_ = [v for v in (_mx("soil_loss_L_max"), _mx("soil_loss_R_max")) if v is not None]
                reaches.append({"ch_start": float(bounds[start]), "ch_end": float(bounds[i]),
                                "length_m": float(bounds[i] - bounds[start]),
                                "worst_score": stations[start]["worst_score"],
                                "n_stations": i - start,
                                "ln_spi_max": max(spi) if spi else None,
                                "soil_loss_max": max(sl_) if sl_ else None})
                start = i
    return stations, reaches


def profile_chart(path, stations, title, class_names=None):
    """Chainage profile PNG (ln SPI max, RUSLE A max, worst score). Needs matplotlib;
    returns False when it is not available (the tool then skips the chart)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # nosec B110 - optional dependency, chart is skipped
        return False
    ch = np.array([s["chainage"] for s in stations]) / 1000.0
    def best(key):
        return np.array([np.nanmax([s.get(f"{key}_L_max") if s.get(f"{key}_L_max") is not None
                                    else np.nan,
                                    s.get(f"{key}_R_max") if s.get(f"{key}_R_max") is not None
                                    else np.nan]) for s in stations])
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        spi = best("ln_spi"); a = best("soil_loss")
    ws = np.array([s["worst_score"] or 0 for s in stations], dtype=float)
    sd = any(s.get("sti_p90_lhs") is not None or s.get("sti_p90_rhs") is not None
             for s in stations)
    fig, axes = plt.subplots(4 if sd else 3, 1, figsize=(12, 9 if sd else 7), sharex=True)
    axes[0].plot(ch, spi, lw=0.8, color="#2b83ba"); axes[0].set_ylabel("max ln(SPI)")
    if np.isfinite(a).any():
        axes[1].plot(ch, a, lw=0.8, color="#d7191c"); axes[1].set_ylabel("max A (t/ha/yr)")
    else:
        axes[1].text(0.5, 0.5, "no RUSLE (LS-only)", ha="center", transform=axes[1].transAxes)
    colours = {0: "#cccccc", 1: "#1a9641", 2: "#a6d96a", 3: "#ffffbf", 4: "#fdae61", 5: "#d7191c"}
    axes[2].bar(ch, ws, width=(ch[1] - ch[0]) if ch.size > 1 else 0.01,
                color=[colours[int(v)] for v in ws])
    axes[2].set_ylabel("worst score"); axes[2].set_ylim(0, 5.5)
    if sd:                                           # STI R3: side-drain feed, per side
        for side, colour in (("lhs", "#5e3c99"), ("rhs", "#e66101")):
            v = np.array([s.get(f"sti_p90_{side}") if s.get(f"sti_p90_{side}") is not None
                          else np.nan for s in stations])
            axes[3].plot(ch, v, lw=0.8, color=colour, label=side.upper())
        axes[3].set_ylabel("STI p90 to drain"); axes[3].legend(loc="upper right", fontsize=8)
    axes[-1].set_xlabel("chainage (km)")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return True
