# -*- coding: utf-8 -*-
"""v0.12 validation: vectorised core == v0.8.3 reference, and the hybrid
flat-resolution family (WP-E). Bare Python + NumPy.

    python -m qeht.tests.test_flats
"""

import sys
from collections import deque

import numpy as np

from ..core.conditioning.fill import fill_depressions
from ..core.flow.direction import d8_direction, _resolve_flats_toward
from ..core.flow.flats import resolve_flats
from ..core.flow.accumulation import flow_accumulation, strahler_order
from ..core.flow import _reference as R
from ..core.grid import DROW, DCOL, receivers_from_direction, NO_RECEIVER
from ..core.watershed.delineate import (delineate_catchment, longest_flow_path,
                                        extract_streams, snap_pour_point)
from ..core.memory import estimate_peak_memory

FAILURES = []
N_RANDOM = 80


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def rand_surface(rng, kind=None):
    rows, cols = int(rng.integers(4, 45)), int(rng.integers(4, 45))
    kind = int(rng.integers(0, 4)) if kind is None else kind
    if kind == 0:                       # quantised noise: many plateaus
        z = np.round(rng.random((rows, cols)) * rng.integers(2, 6))
    elif kind == 1:                     # terraces
        r, c = np.mgrid[0:rows, 0:cols]
        z = np.floor((r + 0.3 * c + rng.random((rows, cols)) * 2) / rng.integers(2, 6))
    elif kind == 2:                     # smooth noise (fill makes flats)
        z = rng.random((rows, cols)) * 10
    else:                               # large flat, bumps, one low outlet
        z = np.full((rows, cols), 5.0)
        z[rng.random((rows, cols)) < 0.05] = 6
        z[0, int(rng.integers(0, cols))] = 4
    valid = np.ones((rows, cols), bool)
    if rng.random() < 0.3:
        valid[rng.random((rows, cols)) < 0.05] = False
    return z, valid


def unrouted_before(z, valid):
    f, _, _ = fill_depressions(z, valid)
    d, _ = d8_direction(f, valid, 1.0, 1.0, resolve_flats=False)
    return f, d


def has_cycle(direction):
    rows, cols = direction.shape
    n = rows * cols
    rec = receivers_from_direction(direction, (rows, cols))
    nxt = np.where(rec == NO_RECEIVER, n, rec)
    nxt = np.append(nxt, n)
    for _ in range(int(np.ceil(np.log2(n + 1))) + 2):
        nxt = nxt[nxt]
    return bool((nxt[:n] != n).any())


# ------------------------------------------------------------------
def test_oracles():
    print(f"\n1. Vectorised == v0.8.3 reference ({N_RANDOM} random surfaces each)")
    rng = np.random.default_rng(2026)
    bad = {k: 0 for k in ("fill", "fill_eps", "d8", "toward", "barnes", "barnes_cycles",
                          "acc", "acc_w", "delineate", "lfp", "strahler", "snap")}
    n_ref_cycles = [0]
    for t in range(N_RANDOM):
        z, v = rand_surface(rng)
        cw, ch = float(rng.choice([1.0, 30.0])), float(rng.choice([1.0, 30.0, 29.7]))
        for key, eps in (("fill", 0.0), ("fill_eps", 1e-3)):
            a = R.reference_fill_depressions(z, v, eps, cw, ch)
            b = fill_depressions(z, v, eps, cw, ch)
            if not (np.array_equal(np.nan_to_num(a[0], nan=-9), np.nan_to_num(b[0], nan=-9))
                    and a[1:] == b[1:]):
                bad[key] += 1
        f, _, _ = fill_depressions(z, v, 0.0, cw, ch)
        d0 = R.reference_d8_core(f, v, cw, ch)
        d1, _ = d8_direction(f, v, cw, ch, resolve_flats=False)
        bad["d8"] += not np.array_equal(d0, d1)
        work = np.where(v, f, np.inf)
        a, b = d1.copy(), d1.copy()
        R.reference_resolve_flats_toward(work, v, a); _resolve_flats_toward(work, v, b)
        bad["toward"] += not np.array_equal(a, b)
        a, b = d1.copy(), d1.copy()
        R.reference_resolve_flats_barnes(f, v, a); resolve_flats(f, v, b, w=2.0)
        # v0.12 fixes the v0.8.3 exit rule (it could form 2-cell cycles):
        # outputs must match except at exit cells, which must now leave the
        # flat to a lower or pre-routed neighbour.
        diff = np.argwhere(a != b)
        for r, c in diff:
            k = b[r, c]
            nr, nc = r + DROW[k], c + DCOL[k]
            ok = (d1[r, c] < 0 and k >= 0 and (f[nr, nc] < f[r, c] - 1e-9 or
                  (d1[nr, nc] >= 0 and abs(f[nr, nc] - f[r, c]) <= 1e-9)))
            if not ok:
                bad["barnes"] += 1
                break
        bad["barnes_cycles"] += has_cycle(b)
        n_ref_cycles[0] += has_cycle(a)
        d = b
        x, sx = R.reference_flow_accumulation(d, v); y, sy = flow_accumulation(d, v)
        bad["acc"] += not (np.array_equal(np.nan_to_num(x, nan=-1), np.nan_to_num(y, nan=-1)) and sx == sy)
        wts = rng.random(z.shape) * 900
        x, _ = R.reference_flow_accumulation(d, v, weights=wts); y, _ = flow_accumulation(d, v, weights=wts)
        bad["acc_w"] += not np.allclose(np.nan_to_num(x), np.nan_to_num(y), rtol=1e-12, atol=1e-9)
        ids = np.flatnonzero(v)
        outs = [tuple(int(q) for q in np.unravel_index(i, z.shape))
                for i in rng.choice(ids, min(3, ids.size), replace=False)]
        bad["delineate"] += not np.array_equal(R.reference_delineate_catchment(d, v, outs),
                                               delineate_catchment(d, v, outs))
        for o in outs:
            m = delineate_catchment(d, v, [o]) > 0
            p = R.reference_longest_flow_path(d, v, o, elevation=z, cell_width=cw,
                                              cell_height=ch, catchment_mask=m)
            q = longest_flow_path(d, v, o, elevation=z, cell_width=cw, cell_height=ch,
                                  catchment_mask=m)
            bad["lfp"] += p != q
        acc, _ = flow_accumulation(d, v)
        st = extract_streams(acc, v, threshold_cells=int(rng.integers(1, 6)))
        bad["strahler"] += not np.array_equal(R.reference_strahler_order(d, v, st),
                                              strahler_order(d, v, st))
        r0, c0 = outs[0]
        rad = int(rng.integers(0, 6))
        bad["snap"] += snap_pour_point(r0, c0, acc, v, rad, st) != \
            R.reference_snap_to_stream(r0, c0, acc, v, rad, st)
    labels = {"barnes": "barnes: identical to v0.8.3 except corrected exit cells",
              "barnes_cycles": "barnes: no flow cycles (v0.8.3 had cycles on "
                               f"{n_ref_cycles[0]} of {N_RANDOM} surfaces)"}
    for k, n in bad.items():
        check(labels.get(k, f"{k}: identical to v0.8.3"), n == 0, f"{n} failing surfaces")


def _toward_steps(f, valid, d0):
    """Min 8-steps from each no-flow cell to a low edge of its flat (1 at the edge)."""
    rows, cols = f.shape
    flat = valid & (d0 < 0)
    dist = np.zeros(f.shape, dtype=np.int64)
    q = deque()
    for r, c in zip(*np.nonzero(flat)):
        for k in range(8):
            nr, nc = r + DROW[k], c + DCOL[k]
            if 0 <= nr < rows and 0 <= nc < cols and valid[nr, nc] and (
                    f[nr, nc] < f[r, c] - 1e-9 or (d0[nr, nc] >= 0 and abs(f[nr, nc] - f[r, c]) <= 1e-9)):
                dist[r, c] = 1; q.append((r, c)); break
    while q:
        r, c = q.popleft()
        for k in range(8):
            nr, nc = r + DROW[k], c + DCOL[k]
            if 0 <= nr < rows and 0 <= nc < cols and flat[nr, nc] and dist[nr, nc] == 0 \
                    and abs(f[nr, nc] - f[r, c]) <= 1e-9:
                dist[nr, nc] = dist[r, c] + 1; q.append((nr, nc))
    return dist


def test_hybrid_family():
    print("\n2. Hybrid family flat_mask = w*toward + (flat_height - away)")
    rng = np.random.default_rng(7)
    surfaces = []
    for i in range(60):
        z, v = rand_surface(rng, kind=i % 4)
        f, d0 = unrouted_before(z, v)
        surfaces.append((f, d0, v))
    for w in (1.01, 1.5, 2.0, 5.0, 100.0):
        cyc = unres = 0
        for f, d0, v in surfaces:
            ref = d0.copy(); resolve_flats(f, v, ref, w=2.0)
            dd = d0.copy(); resolve_flats(f, v, dd, w=w)
            cyc += has_cycle(dd)
            unres += int((v & (dd < 0)).sum()) != int((v & (ref < 0)).sum())
        check(f"w={w:g}: no cycles; resolves exactly the cells Barnes (w=2) resolves",
              cyc == 0 and unres == 0, f"{cyc} with cycles, {unres} with different unrouted counts")
    bad = 0
    for f, d0, v in surfaces:
        dd = d0.copy(); resolve_flats(f, v, dd, w=1e6)
        steps = _toward_steps(f, v, d0)
        rr, cc = np.nonzero((steps > 1) & (dd >= 0))
        for r, c in zip(rr, cc):
            k = dd[r, c]
            if steps[r + DROW[k], c + DCOL[k]] != steps[r, c] - 1:
                bad += 1
    check("w -> inf: every flat cell steps along a shortest path to its outlet", bad == 0,
          f"{bad} cells off a shortest path")
    try:
        resolve_flats(surfaces[0][0], surfaces[0][2], surfaces[0][1].copy(), w=1.0)
        refused = False
    except ValueError:
        refused = True
    check("w <= 1 refused (false sinks possible)", refused)
    same = all(np.array_equal(d8_direction(f, v, 1, 1, flat_method="hybrid", flat_weight=2.0)[0],
                              d8_direction(f, v, 1, 1, flat_method="barnes")[0])
               for f, _, v in surfaces[:15])
    check("d8_direction: hybrid w=2 == barnes", same)


def test_size_switch():
    print("\n3. Size switch (toward-only below N cells)")
    rng = np.random.default_rng(11)
    bad_all = bad_off = seams = 0
    for i in range(40):
        z, v = rand_surface(rng, kind=i % 4)
        f, d0 = unrouted_before(z, v)
        a = d0.copy(); resolve_flats(f, v, a, w=2.0, small_flat_cells=0)
        b = d0.copy(); resolve_flats(f, v, b, w=2.0)
        bad_off += not np.array_equal(a, b)
        c = d0.copy(); resolve_flats(f, v, c, w=2.0, small_flat_cells=10 ** 9)
        bad_all += has_cycle(c) or int((v & (c < 0)).sum()) != int((v & (b < 0)).sum())
        m = d0.copy(); resolve_flats(f, v, m, w=2.0, small_flat_cells=20)
        seams += has_cycle(m) or int((v & (m < 0)).sum()) != int((v & (b < 0)).sum())
    check("threshold 0 == plain Barnes", bad_off == 0)
    check("all flats toward-only: no cycles, same cells resolved", bad_all == 0)
    check("mixed threshold: no cycles, no unresolved seams", seams == 0)


def test_memory_estimate():
    print("\n4. Memory estimate")
    e = estimate_peak_memory(5000, 5000)
    check("per-step estimate reported, peak is the maximum", e["peak_gb"] == max(e["steps_gb"].values())
          and 1.0 < e["peak_gb"] < 10.0, f"{e['peak_gb']:.1f} GB for 25 Mcells")


def main(argv=None):
    print("=" * 62)
    print("QEHT FLATS & PERFORMANCE - v0.12 (WP-E)")
    print("=" * 62)
    test_oracles()
    test_hybrid_family()
    test_size_switch()
    test_memory_estimate()
    print("\n" + "=" * 62)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S):")
        for f in FAILURES:
            print("   -", f)
        return 1
    print("ALL CHECKS PASSED")
    print("=" * 62)
    return 0


if __name__ == "__main__":
    sys.exit(main())
