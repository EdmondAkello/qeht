# -*- coding: utf-8 -*-
"""v0.20 validation: side-drain siltation indicator along the corridor (STI R3).
Bare Python + NumPy, no QGIS, no GDAL. Synthetic planes and channels.

    python -m qeht.tests.test_corridor_sti
"""

import json
import math
import sys

import numpy as np

from ..core.erosion.side_drain import (station_sti, class_breaks, classify, reach_sti,
                                       corridor_sti, sti_class, params_json, STATION_FIELDS,
                                       REACH_FIELDS)
from ..core.erosion.corridor import sample_corridor
from ..core.network.alignment import Alignment
from ..core.flow.direction import d8_direction

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


ROWS, COLS, CS = 60, 80, 10.0
GT = (500000.0, CS, 0.0, 9901000.0, 0.0, -CS)
ROAD_ROW = 30


def road():
    """E-W road along the centre of row 30, chainage rising eastward (LHS = north)."""
    y = GT[3] - (ROAD_ROW + 0.5) * CS
    return Alignment([[(GT[0] + 55.0, y), (GT[0] + 745.0, y)]])


def plane(south=0.1):
    r, c = np.mgrid[0:ROWS, 0:COLS].astype(float)
    return 500.0 - south * r * CS + 0.0001 * c       # falls to the south everywhere


def d8(dem):
    valid = np.ones(dem.shape, bool)
    d, _ = d8_direction(dem, valid, CS, CS)
    return d, valid


def test_plane():
    print("\n1-2. Tilted plane draining toward the road from the left (north) only")
    dem = plane()
    d, valid = d8(dem)
    sti = np.full(dem.shape, 5.0)
    st = station_sti(road(), sti, d, valid, GT, step=10.0, half_width=50.0, offset_step=10.0)
    check("LHS (uphill, draining to the road): 5 cells kept at every station",
          all(s["n_lhs"] == 5 for s in st), str(sorted({s["n_lhs"] for s in st})))
    check("RHS (draining away): n_rhs = 0 and the values are NULL",
          all(s["n_rhs"] == 0 and s["sti_p50_rhs"] is None and s["sti_p90_rhs"] is None for s in st))
    check("uniform STI = 5 on the kept cells: p50 = p90 = 5",
          all(s["sti_p50_lhs"] == 5.0 and s["sti_p90_lhs"] == 5.0 for s in st))
    g = np.where(np.arange(ROWS)[:, None] < ROAD_ROW, np.arange(ROWS)[:, None] * 1.0, 0.0) \
        * np.ones((1, COLS))                           # STI = row number north of the road
    st2 = station_sti(road(), g, d, valid, GT, step=10.0, half_width=50.0, offset_step=10.0)
    check("STI = row number (rows 25-29 kept): p50 = 27, p90 = 28.6 (linear percentiles)",
          all(math.isclose(s["sti_p50_lhs"], 27.0) and math.isclose(s["sti_p90_lhs"], 28.6)
              for s in st2))
    st3 = station_sti(road(), sti, d, valid, GT, step=10.0, half_width=20.0, offset_step=10.0)
    check("half width 20 m: 2 cells kept", all(s["n_lhs"] == 2 for s in st3))
    st4 = station_sti(road(), sti, d, valid, GT, step=50.0, half_width=50.0, offset_step=10.0)
    check("stations every 50 m from the start chainage (0 ... 690)",
          [s["chainage"] for s in st4][:3] == [0.0, 50.0, 100.0] and st4[-1]["chainage"] == 690.0)


def test_channel():
    print("\n3. A channel crossing the corridor")
    r, c = np.mgrid[0:ROWS, 0:COLS].astype(float)
    c0 = 40
    dem = plane() + np.where(np.abs(c - c0) <= 3, 0.5 * np.abs(c - c0) * CS, 1.5 * CS)
    d, valid = d8(dem)
    chan = np.zeros(dem.shape, bool)
    chan[:, c0] = True
    sti = np.where(chan, np.nan, 5.0)
    st = station_sti(road(), sti, d, valid, GT, chan, step=10.0, half_width=50.0, offset_step=10.0)
    x_of = {s["chainage"]: (s["x"] - GT[0]) / CS for s in st}
    near = [s for s in st if abs(x_of[s["chainage"]] - (c0 + 0.5)) <= 3.0]
    far = [s for s in st if abs(x_of[s["chainage"]] - (c0 + 0.5)) >= 6.0]
    check("cells draining into the channel are dropped (path stopped): n_lhs = 0 near it",
          near and all(s["n_lhs"] == 0 for s in near), str([s["n_lhs"] for s in near]))
    check("away from the channel the slope still drains to the road (5 cells)",
          far and all(s["n_lhs"] == 5 for s in far))
    st_nc = station_sti(road(), sti, d, valid, GT, None, step=10.0, half_width=50.0, offset_step=10.0)
    on = [s for s in st_nc if abs(x_of[s["chainage"]] - (c0 + 0.5)) < 0.6]
    check("without a channel mask the channel cells (STI NoData) still never count",
          on and all(s["n_lhs"] == 0 or s["sti_p90_lhs"] == 5.0 for s in on))


def test_classes():
    print("\n4. Relative classes from the corridor's own percentiles")
    rng = np.random.default_rng(7)
    st = [{"chainage": 10.0 * i, "sti_p90_lhs": float(v), "sti_p90_rhs": float(w)}
          for i, (v, w) in enumerate(rng.gamma(2.0, 3.0, size=(200, 2)))]
    br = class_breaks(st)
    allv = np.array([s[k] for s in st for k in ("sti_p90_lhs", "sti_p90_rhs")])
    check("breaks are the pooled P50 / P75 / P90", np.allclose(br, np.percentile(allv, [50, 75, 90])))
    classify(st, br)
    cls = [s[k] for s in st for k in ("sti_class_lhs", "sti_class_rhs")]
    shares = {c: cls.count(c) / len(cls) for c in ("low", "moderate", "high", "very high")}
    check("shares follow the breaks: 50 / 25 / 15 / 10 %",
          abs(shares["low"] - 0.5) < 0.01 and abs(shares["moderate"] - 0.25) < 0.01
          and abs(shares["high"] - 0.15) < 0.01 and abs(shares["very high"] - 0.10) < 0.01,
          str({k: round(v, 3) for k, v in shares.items()}))
    st2 = [{"chainage": s["chainage"], "sti_p90_lhs": 2 * s["sti_p90_lhs"],
            "sti_p90_rhs": 2 * s["sti_p90_rhs"]} for s in st]
    classify(st2, class_breaks(st2))
    check("doubled STI everywhere: identical classes (scale-free)",
          all(a["sti_class_lhs"] == b["sti_class_lhs"] and a["sti_class_rhs"] == b["sti_class_rhs"]
              for a, b in zip(st, st2)))
    check("edges: = P50 is low, just above P90 is very high, NULL stays NULL",
          sti_class(br[0], br) == "low" and sti_class(br[2] + 1e-9, br) == "very high"
          and sti_class(None, br) is None)


def test_flag_and_reaches():
    print("\n5. Siltation flag and reaches")
    br = [1.0, 2.0, 3.0]
    st = [{"chainage": 10.0 * i, "sti_p90_lhs": v, "sti_p90_rhs": w}
          for i, (v, w) in enumerate([(3.5, 0.5), (2.5, 2.5), (3.5, 3.5), (1.5, None), (2.5, 0.5)])]
    slope = [0.5, -0.8, 2.0, 0.1, float("nan")]
    classify(st, br, slope, drain_slope_pct=1.0)
    check("flag only where class >= high and |slope| < 1 %",
          [s["siltation_lhs"] for s in st] == [1, 1, 0, 0, None]
          and [s["siltation_rhs"] for s in st] == [0, 1, 0, None, None],
          f"{[s['siltation_lhs'] for s in st]} / {[s['siltation_rhs'] for s in st]}")
    check("classes", [s["sti_class_lhs"] for s in st] == ["very high", "high", "very high",
                                                          "moderate", "high"])
    reaches = [{"ch_start": 0.0, "ch_end": 25.0}, {"ch_start": 25.0, "ch_end": 40.0}]
    rr = reach_sti(reaches, st, br)
    check("reach 1 (stations 0-20): p90 of 3.5/2.5/3.5, max 3.5, flagged length 5 + 10 m on LHS",
          math.isclose(rr[0]["sti_p90_lhs"], 3.5) and rr[0]["sti_max_lhs"] == 3.5
          and rr[0]["sti_class_lhs"] == "very high" and rr[0]["siltation_len_lhs_m"] == 15.0
          and rr[0]["siltation_len_rhs_m"] == 10.0, str(rr[0]))
    check("reach 2 (stations 30-40): no flagged length; RHS p90 from the one value",
          rr[1]["siltation_len_lhs_m"] == 0.0 and rr[1]["sti_p90_rhs"] == 0.5
          and rr[1]["sti_class_rhs"] == "low", str(rr[1]))
    dem = plane()
    d, valid = d8(dem)
    al = road()
    score = np.ones(dem.shape, np.uint8)
    _, reaches = sample_corridor(al, {}, score, GT, step=10.0)
    st, rr, br = corridor_sti(al, np.full(dem.shape, 5.0), d, valid, GT, reaches,
                              slope_at=lambda ch: np.zeros(ch.size))
    check("corridor_sti on the plane: one reach, LHS p90 5, RHS NULL, no flag (uniform = low)",
          len(rr) == 1 and rr[0]["sti_p90_lhs"] == 5.0 and rr[0]["sti_p90_rhs"] is None
          and rr[0]["siltation_len_lhs_m"] == 0.0 and all(s["slope_long_pct"] == 0.0 for s in st))


def test_contract():
    print("\n6. Contract")
    from ..core.interop.field_dictionary import OPTIONAL_LAYERS, METADATA_KEYS, DOWNSTREAM_USE
    lay = OPTIONAL_LAYERS.get("corridor_sti")
    check("optional layer corridor_sti (LINESTRING) carries every reach field",
          lay is not None and lay[0] == "LINESTRING"
          and {f for f, _ in REACH_FIELDS} <= {f[0] for f in lay[1]})
    check("metadata key corridor_sti_params_json",
          "corridor_sti_params_json" in [k for k, _ in METADATA_KEYS])
    check("downstream-use term 'side-drain siltation'", "side-drain siltation" in DOWNSTREAM_USE)
    check("station fields named as specified",
          [f for f, _ in STATION_FIELDS] == ["sti_p50_lhs", "sti_p90_lhs", "n_lhs", "sti_p50_rhs",
                                             "sti_p90_rhs", "n_rhs", "sti_class_lhs",
                                             "sti_class_rhs", "slope_long_pct", "siltation_lhs",
                                             "siltation_rhs"])
    pj = json.loads(params_json(10.0, 50.0, 10.0, 5.0, 1.0, [1, 2, 3], "alignment profile"))
    check("params JSON records the breaks and says it is relative and not a severity",
          pj["class_breaks_p50_p75_p90"] == [1, 2, 3] and "not a severity" in pj["note"])


def main(argv=None):
    test_plane()
    test_channel()
    test_classes()
    test_flag_and_reaches()
    test_contract()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
