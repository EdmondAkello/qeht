# -*- coding: utf-8 -*-
"""Time the routing core on a synthetic DEM (repository tool; not in the plugin zip).

    python -m qeht.tools.dev.bench_core [cells_per_side ...]     (default 1000 2000)

Run from the folder that contains the `qeht` folder. The terrain is
synthetic: four octaves of smooth noise on a regional slope, with the
northern third stored in whole metres so that it has real flats. Seeded,
so runs are comparable across versions and machines.
"""

import sys
import time

import numpy as np

from qeht.core.conditioning.fill import fill_depressions
from qeht.core.flow.direction import d8_direction
from qeht.core.flow.accumulation import flow_accumulation
from qeht.core.raster import resampling_stats, precision_stats


def synthetic_dem(n, seed=0):
    rng = np.random.default_rng(seed)
    z = np.zeros((n, n))
    x = np.linspace(0.0, 1.0, n)
    for k, amp in ((8, 300.0), (32, 60.0), (128, 12.0), (512, 3.0)):
        g = rng.normal(0.0, 1.0, (k + 1, k + 1))
        knots = np.linspace(0.0, 1.0, k + 1)
        rows = np.array([np.interp(x, knots, g[i]) for i in range(k + 1)])     # (k+1, n)
        z += amp * np.array([np.interp(x, knots, rows[:, j]) for j in range(n)]).T
    z += 1500.0 + np.linspace(0.0, 200.0, n)[:, None]
    z[: n // 3] = np.round(1500.0 + 0.05 * (z[: n // 3] - 1500.0))
    return z


def timed(label, fn):
    t0 = time.perf_counter()
    out = fn()
    print(f"  {label:30s} {time.perf_counter() - t0:8.2f} s", flush=True)
    return out


def main(argv):
    sizes = [int(a) for a in argv] or [1000, 2000]
    cs = 30.0
    for n in sizes:
        z = synthetic_dem(n)
        v = np.ones_like(z, dtype=bool)
        print(f"{n} x {n} = {n * n / 1e6:.1f} Mcells")
        timed("audits (resampling, precision)", lambda: (resampling_stats(z, v), precision_stats(z, v)))
        filled, _, _ = timed("fill", lambda: fill_depressions(z, v, cell_width=cs, cell_height=cs))
        d, _ = timed("D8 + Barnes flats", lambda: d8_direction(filled, v, cs, cs, flat_method="barnes"))
        timed("D8 + toward-lower flats", lambda: d8_direction(filled, v, cs, cs, flat_method="toward"))
        timed("accumulation", lambda: flow_accumulation(d, v))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
