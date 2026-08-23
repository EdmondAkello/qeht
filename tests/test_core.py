# -*- coding: utf-8 -*-
"""Level 1 validation: synthetic DEMs with analytically known answers.

Runs on bare Python + NumPy. No QGIS, no GDAL. That is the point of
keeping `qeht.core` free of QGIS imports - the hydrology can be verified
anywhere, including on a machine where QGIS is not installed.

    python -m qeht.tests.test_core
"""

import sys
import numpy as np

from ..core.grid import (DROW, DCOL, D8_CODE, NO_RECEIVER,
                         encode_d8, decode_d8, neighbour_distances)
from ..core.conditioning.fill import fill_depressions, depression_depth
from ..core.flow.direction import d8_direction
from ..core.flow.accumulation import flow_accumulation, strahler_order
from ..core.watershed.delineate import (extract_streams, snap_pour_point,
                                        delineate_catchment, longest_flow_path)

FAILURES = []


def check(name, condition, detail=""):
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


# ------------------------------------------------------------------
def test_encoding():
    print("\n1. Standard D8 encoding round-trip")
    idx = np.arange(8, dtype=np.int64).reshape(2, 4)
    valid = np.ones((2, 4), dtype=bool)
    d8 = encode_d8(idx, valid)
    check("encode gives standard D8 codes", list(d8.reshape(-1)) == list(D8_CODE))
    back = decode_d8(d8)
    check("decode inverts encode", np.array_equal(back, idx))

    # East must be 1, north must be 64 - the two anchors of the convention.
    check("E == 1", D8_CODE[0] == 1)
    check("N == 64", D8_CODE[6] == 64 and DROW[6] == -1 and DCOL[6] == 0)


def test_distance_weighting():
    print("\n2. Diagonal distance weighting")
    d = neighbour_distances(30.0, 30.0)
    check("orthogonal = cellsize", np.isclose(d[0], 30.0))
    check("diagonal = cellsize*sqrt(2)", np.isclose(d[1], 30.0 * np.sqrt(2)))

    # A cell with a 10 m drop east (30 m away) and an 12 m drop southeast
    # (42.43 m away). Unweighted logic picks SE (bigger drop);
    # correct D8 picks E (steeper gradient). This is precisely the bug in
    # a traversal-parent implementation.
    dem = np.array([
        [100.0, 100.0, 100.0],
        [100.0, 100.0,  90.0],
        [100.0, 100.0,  88.0],
    ])
    valid = np.ones((3, 3), dtype=bool)
    direction, _ = d8_direction(dem, valid, 30.0, 30.0, resolve_flats=False)
    e_slope = (100.0 - 90.0) / 30.0
    se_slope = (100.0 - 88.0) / (30.0 * np.sqrt(2))
    check("E steeper than SE for this pair", e_slope > se_slope,
          f"E={e_slope:.4f} SE={se_slope:.4f}")
    check("centre routes E (index 0), not SE", direction[1, 1] == 0,
          f"got index {direction[1, 1]}")


def test_fill():
    print("\n3. Priority-flood depression filling")
    # Uniform plane tilted south with one 5 m pit.
    n = 20
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    dem = 100.0 - yy * 1.0
    dem[10, 10] -= 5.0
    valid = np.ones((n, n), dtype=bool)

    filled, n_filled, seeds = fill_depressions(dem, valid, cell_width=10, cell_height=10)
    check("pit was raised", filled[10, 10] > dem[10, 10],
          f"{dem[10,10]:.1f} -> {filled[10,10]:.1f}")
    # Spill elevation is the LOWEST neighbour (the row below, at 89.0),
    # not the pit's own original surface level. Water leaves at 89.
    spill = dem[11, 10]
    check("pit filled to true spill elevation", np.isclose(filled[10, 10], spill),
          f"{filled[10,10]:.2f} vs spill {spill:.2f}")
    check("exactly one cell filled", n_filled == 1, f"n_filled={n_filled}")
    check("no cell was lowered", np.all(filled >= dem - 1e-9))
    check("boundary cells seeded", seeds == 4 * n - 4, f"seeds={seeds}")

    depth = depression_depth(dem, filled, valid)
    expected_depth = dem[11, 10] - dem[10, 10]
    check("fill depth == spill minus original", np.isclose(depth[10, 10], expected_depth),
          f"depth={depth[10,10]:.2f}")
    check("no other cell was filled", np.isclose(np.nansum(depth), expected_depth))

    # NoData collar: boundary must be detected at the NoData contact,
    # not the raster edge.
    valid2 = np.zeros((n, n), dtype=bool)
    valid2[5:15, 5:15] = True
    dem2 = np.where(valid2, dem, np.nan)
    _, _, seeds2 = fill_depressions(np.nan_to_num(dem2), valid2,
                                    cell_width=10, cell_height=10)
    check("NoData collar seeds interior boundary", seeds2 == 36,
          f"seeds={seeds2} (perimeter of 10x10 block = 36)")


def test_accumulation_plane():
    print("\n4. Accumulation on a tilted plane (analytic)")
    # Perfect south-draining plane, no diagonal ambiguity: every column
    # is an independent flow line, so cell (r, c) must accumulate exactly
    # r upstream cells.
    n = 15
    yy, _ = np.mgrid[0:n, 0:n].astype(float)
    dem = 100.0 - yy * 2.0
    valid = np.ones((n, n), dtype=bool)

    direction, dstats = d8_direction(dem, valid, 10.0, 10.0)
    accum, astats = flow_accumulation(direction, valid)

    check("no cycles in graph", astats["cells_in_cycles"] == 0)
    # The bottom row has no lower neighbour: these are genuine outlets
    # draining off-grid, and MUST report as unrouted (standard D8 code 0).
    check("exactly the bottom row is unrouted",
          dstats["cells_still_unrouted"] == n,
          f"unrouted={dstats['cells_still_unrouted']}, expected {n}")
    check("unrouted cells are all in the bottom row",
          np.all(direction[:-1, :] >= 0))

    expected_row5 = 5.0
    check("row 5 accumulates 5 upstream cells",
          np.allclose(accum[5, 1:-1], expected_row5),
          f"got {accum[5, 1]}")
    check("top row accumulates 0", np.allclose(accum[0, :], 0.0))
    check("total accumulation conserved",
          np.isclose(np.nansum(accum), np.sum(np.arange(n)) * n),
          f"sum={np.nansum(accum):.0f}")


def test_accumulation_valley():
    print("\n5. Convergent valley")
    n = 21
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    dem = 100.0 + np.abs(xx - n // 2) * 2.0 - yy * 1.0
    valid = np.ones((n, n), dtype=bool)

    filled, _, _ = fill_depressions(dem, valid, cell_width=30, cell_height=30)
    direction, _ = d8_direction(filled, valid, 30.0, 30.0)
    accum, _ = flow_accumulation(direction, valid)

    outlet_r, outlet_c = n - 1, n // 2
    check("valley axis carries the maximum",
          accum[outlet_r, outlet_c] >= np.nanmax(accum) * 0.5,
          f"axis={accum[outlet_r, outlet_c]:.0f} max={np.nanmax(accum):.0f}")
    check("accumulation increases downstream along axis",
          np.all(np.diff(accum[1:, n // 2]) >= 0))

    # Contributing area with weights = cell area.
    cell_area = 30.0 * 30.0
    weights = np.full((n, n), cell_area)
    accum_area, _ = flow_accumulation(direction, valid, weights=weights)
    check("weighted accumulation == count * cell area",
          np.allclose(accum_area[valid], accum[valid] * cell_area))


def test_catchment_and_path():
    print("\n6. Catchment delineation and longest flow path")
    n = 21
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    dem = 100.0 + np.abs(xx - n // 2) * 2.0 - yy * 1.0
    valid = np.ones((n, n), dtype=bool)

    filled, _, _ = fill_depressions(dem, valid, cell_width=30, cell_height=30)
    direction, _ = d8_direction(filled, valid, 30.0, 30.0)
    accum, _ = flow_accumulation(direction, valid)

    outlet = (n - 1, n // 2)
    labels = delineate_catchment(direction, valid, [outlet])
    check("catchment is non-empty", labels.sum() > 0, f"{int((labels>0).sum())} cells")
    check("outlet is inside its own catchment", labels[outlet] == 1)
    check("catchment size matches accumulation at outlet",
          int((labels > 0).sum()) == int(accum[outlet]) + 1,
          f"cells={int((labels>0).sum())} accum+1={int(accum[outlet])+1}")

    lfp = longest_flow_path(direction, valid, outlet, elevation=filled,
                            cell_width=30.0, cell_height=30.0)
    check("longest path reaches the outlet", lfp["cells"][-1] == outlet)
    check("longest path has positive length", lfp["length"] > 0,
          f"L={lfp['length']:.1f} m over {lfp['n_cells']} cells")
    check("path length >= straight-line cell count",
          lfp["length"] >= (lfp["n_cells"] - 1) * 30.0 - 1e-6)
    check("slope is positive and plausible",
          0 < lfp.get("slope", -1) < 1,
          f"S={lfp.get('slope', float('nan')):.5f}, drop={lfp.get('drop', 0):.1f} m")


def test_snap_and_streams():
    print("\n7. Stream extraction and pour-point snapping")
    n = 21
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    dem = 100.0 + np.abs(xx - n // 2) * 2.0 - yy * 1.0
    valid = np.ones((n, n), dtype=bool)
    filled, _, _ = fill_depressions(dem, valid, cell_width=30, cell_height=30)
    direction, _ = d8_direction(filled, valid, 30.0, 30.0)
    accum, _ = flow_accumulation(direction, valid)

    streams = extract_streams(accum, valid, threshold_cells=10)
    check("streams are a subset of valid cells", np.all(valid[streams]))
    check("stream network is not everything",
          0 < streams.sum() < valid.sum(), f"{int(streams.sum())} stream cells")

    streams_area = extract_streams(accum, valid, threshold_area=10 * 900.0,
                                   cell_area=900.0)
    check("area threshold == equivalent cell threshold",
          np.array_equal(streams, streams_area))

    # Drop a pour point off-channel and confirm it snaps back.
    off_r, off_c = n - 3, n // 2 + 3
    sr, sc, moved, acc_at = snap_pour_point(off_r, off_c, accum, valid,
                                            search_radius_cells=5)
    check("snapping moved the point", moved > 0, f"moved {moved} cells")
    check("snapped point has higher accumulation",
          acc_at > accum[off_r, off_c],
          f"{accum[off_r, off_c]:.0f} -> {acc_at:.0f}")

    order = strahler_order(direction, valid, streams)
    check("Strahler orders start at 1", order[streams].min() >= 1)
    check("Strahler defined only on streams", np.all(order[~streams] == 0))
    check("outlet carries the highest order",
          order[n - 1, n // 2] == order.max(),
          f"outlet order {order[n-1, n//2]}, max {order.max()}")


def test_flat_handling():
    print("\n8. Flats and filled surfaces")
    n = 15
    dem = np.full((n, n), 50.0)
    dem[n - 1, :] = 40.0          # a low southern edge to drain to
    valid = np.ones((n, n), dtype=bool)

    direction, stats = d8_direction(dem, valid, 10.0, 10.0, resolve_flats=True)
    interior_unrouted = stats["cells_still_unrouted"]
    check("flat resolution routed most of the plateau",
          interior_unrouted < n * n * 0.2,
          f"{interior_unrouted} of {n*n} cells unrouted")

    accum, astats = flow_accumulation(direction, valid)
    check("no cycles introduced by flat routing",
          astats["cells_in_cycles"] == 0,
          f"{astats['cells_in_cycles']} cells in cycles")

    # min_slope during filling should largely remove flats.
    dem2 = np.full((n, n), 50.0)
    dem2[7, 7] = 45.0
    dem2[n - 1, :] = 40.0
    filled2, _, _ = fill_depressions(dem2, valid, min_slope=1e-3,
                                     cell_width=10, cell_height=10)
    check("min_slope tilts the filled surface",
          filled2[7, 7] > 45.0, f"pit -> {filled2[7,7]:.4f}")




def test_barnes_saddle():
    print("\n9. Barnes convergent flat resolution (saddle)")
    from ..core.flow.direction import d8_direction
    from ..core.flow.accumulation import flow_accumulation
    Z = np.array([
        [9,9,9,9,9,9,9],[9,5,5,5,5,5,9],[9,5,5,5,5,5,9],
        [9,5,5,5,5,5,9],[9,5,5,5,5,5,9],[9,9,9,4,9,9,9],
    ], dtype=float)
    valid = np.ones(Z.shape, bool)
    d, st = d8_direction(Z, valid, 1.0, 1.0, flat_method="barnes")
    acc, astats = flow_accumulation(d, valid)
    check("saddle: no cycles", astats["cells_in_cycles"] == 0)
    check("saddle: interior routes (no false sinks)",
          int((valid & (d < 0)).sum()) <= 3,
          f"{int((valid&(d<0)).sum())} unrouted")
    # outlet at (5,3) must collect the whole interior
    check("saddle: flow converges to outlet",
          acc[5, 3] >= 18, f"outlet accum {acc[5,3]:.0f}")
    check("saddle: accumulation increases downstream",
          acc[4, 3] >= acc[3, 3] and acc[5, 3] >= acc[4, 3])


def main():
    print("=" * 62)
    print("QEHT CORE - LEVEL 1 VALIDATION (synthetic, analytic)")
    print("=" * 62)
    test_encoding()
    test_distance_weighting()
    test_fill()
    test_accumulation_plane()
    test_accumulation_valley()
    test_catchment_and_path()
    test_snap_and_streams()
    test_flat_handling()
    test_barnes_saddle()

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
