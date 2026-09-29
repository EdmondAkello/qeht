# WP-G benchmark: flat resolution across terrain classes (QEHT 0.13)

Run on 29 September 2026. This file records the test design, what was run, the results and the decision. Raw per-run results and the scripts sit outside the repository, in the working folder under `data/benchmark/` and `data/_scripts/`. The DEM clips are not redistributable: FABDEM is CC BY-NC-SA 4.0, and the redistribution terms for ALOS AW3D30 are unconfirmed.

## Question

QEHT 0.12 offered three ways to route flow across flats: toward-lower, Barnes 2014, and a weighted hybrid. The rule for this benchmark was agreed before the full run:

- If the three give comparable results in every terrain class, keep one of them as the only method.
- Two methods count as comparable when their stream networks differ from each other by no more than two established tools differ from each other (TauDEM vs RichDEM).
- The recommended method is the QEHT method that agrees best with TauDEM. TauDEM is the independent judge because, after the fix below, QEHT Barnes *is* RichDEM's algorithm.

## Areas (sampling design, D5)

- **Frame.** The FABDEM 30 m grid, where it overlaps ALOS AW3D30 (FABDEM stops at 4.62° N). ESA WorldCover 2021 at 30 m uses the same grid.
- **Candidates.** Non-overlapping tiles of 500 × 500 cells (about 15 km on a side). A tile is eligible when more than 99 % of its cells are valid and less than 2 % is open water. The water limit exists because lake surfaces are artificial flats. 2,360 tiles are eligible.
- **Terrain class.** Each tile is classed by the median, over its ~300 m blocks, of the block's mean Horn slope of the 30 m DEM:

  | Class | Median slope | Eligible tiles |
  |---|---|---:|
  | Flat | < 3 % | 1,553 |
  | Rolling | 3–10 % | 541 |
  | Hilly | 10–25 % | 218 |
  | Mountainous | ≥ 25 % | 48 |

  Rolling, hilly and mountainous were requested. Flat was added because flat routing is what this benchmark decides.
- **Sample.** Stratified random, 20 tiles per class, seed 2026, giving 80 areas. Each area is clipped from both FABDEM (floating point, bare earth) and ALOS (whole metres, surface model), so there are 160 runs.
- **Scaling up.** Every eligible tile is kept in `aoi_index.geojson` with a `selected` flag. The same frame can later supply randomly offset windows for the planned runs of about 100,000 areas.

The median share of cells with no downhill neighbour after filling, per class and DEM:

| DEM | Flat | Rolling | Hilly | Mountainous |
|---|---:|---:|---:|---:|
| ALOS | 53 % | 8.6 % | 4.4 % | 1.9 % |
| FABDEM | 1.2 % | 0.9 % | 0.9 % | 0.7 % |

On whole-metre DEMs the choice of flat method matters a great deal. On floating-point DEMs it matters much less.

## Methods compared

| Method | What it is |
|---|---|
| QEHT toward | Shortest path to the flat's outlet (the default up to 0.12) |
| QEHT barnes | Barnes, Lehman & Mulla 2014 |
| QEHT hybrid15 / hybrid10 | w·toward + (height − away), with w = 1.5 and w = 10 |
| TauDEM (Develop, built from source) | `pitremove` + `d8flowdir`, which uses Garbrecht & Martz 1997 flats |
| RichDEM (master, built from source) | `rd_d8_flowdirs`: PriorityFlood + Barnes 2014, written by Barnes |
| MAS Spatial Analysis Tool 1.2.1 | Epsilon fill + resolve flats toward outlets + D8 |
| DDM HydroLogic 2.3 | Priority-flood receiver graph |

The results are read in two views:

- **iso (flat-isolated).** Every method shares QEHT's distance-weighted D8 away from flats. Only cells with no descent after a common flat fill take the method's own direction. This isolates the flat algorithm and covers QEHT, TauDEM and RichDEM.
- **e2e (end to end).** Each tool runs exactly as shipped.

In both views, one common accumulation routine runs on every direction grid.

## Metrics

All figures are medians over the 10 runs in each class and DEM.

- **Stream F1.** Streams at 1 km² contributing area, compared with 1-cell tolerance. Exact stream IoU is also reported.
- **Basin agreement.** For basins of 1 km² or more: area-weighted best-match IoU of the terminal basins.
- **Flat-direction agreement.** The share of flat cells that are given the same direction.
- **Integrity checks.** Flow loops and unrouted flat cells.

## Results (iso view)

Stream F1 at 1 km²:

| Pair | Flat ALOS | Flat FABDEM | Rolling ALOS | Rolling FABDEM | Hilly ALOS | Hilly FABDEM | Mount. ALOS | Mount. FABDEM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| toward vs Barnes | 0.811 | 0.852 | 0.962 | 0.972 | 0.987 | 0.985 | 0.991 | 0.988 |
| Barnes vs hybrid w=1.5 | 0.993 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Barnes vs hybrid w=10 | 0.990 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| **TauDEM vs RichDEM (reference spread)** | 0.870 | 0.940 | 0.979 | 0.985 | 0.996 | 0.996 | 0.996 | 0.990 |
| toward vs TauDEM | 0.733 | **0.954** | 0.949 | 0.983 | 0.987 | 0.986 | 0.996 | **0.996** |
| Barnes vs TauDEM | **0.870** | 0.940 | **0.979** | **0.985** | **0.996** | **0.996** | 0.996 | 0.990 |
| Barnes vs RichDEM | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

Exact stream IoU at 1 km²:

| Pair | Flat ALOS | Flat FABDEM | Rolling ALOS | Rolling FABDEM | Hilly ALOS | Hilly FABDEM | Mount. ALOS | Mount. FABDEM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| toward vs Barnes | 0.394 | 0.679 | 0.644 | 0.911 | 0.823 | 0.930 | 0.890 | 0.945 |
| TauDEM vs RichDEM | 0.658 | 0.862 | 0.920 | 0.950 | 0.978 | 0.982 | 0.989 | 0.977 |
| toward vs TauDEM | 0.312 | 0.881 | 0.646 | 0.915 | 0.820 | 0.929 | 0.906 | 0.974 |
| Barnes vs TauDEM | 0.658 | 0.862 | 0.920 | 0.950 | 0.978 | 0.982 | 0.989 | 0.977 |

Flat-direction agreement:

| Pair | Flat ALOS | Flat FABDEM | Rolling ALOS | Rolling FABDEM | Hilly ALOS | Hilly FABDEM | Mount. ALOS | Mount. FABDEM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| toward vs Barnes | 0.576 | 0.637 | 0.621 | 0.828 | 0.717 | 0.822 | 0.768 | 0.891 |
| Barnes vs TauDEM | 0.908 | 0.916 | 0.907 | 0.941 | 0.949 | 0.945 | 0.957 | 0.943 |
| toward vs TauDEM | 0.567 | 0.648 | 0.618 | 0.830 | 0.727 | 0.814 | 0.777 | 0.900 |

Basin agreement is 0.89–1.00 for every pair and stratum. It is lowest on flat terrain, where toward vs TauDEM scores 0.889 on ALOS and Barnes vs TauDEM scores 0.949.

## Results (end to end, stream F1 at 1 km²)

| Pair | Flat ALOS | Flat FABDEM | Rolling ALOS | Rolling FABDEM | Hilly ALOS | Hilly FABDEM | Mount. ALOS | Mount. FABDEM |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Barnes vs TauDEM | 0.861 | 0.935 | 0.968 | 0.979 | 0.993 | 0.994 | 0.987 | 0.988 |
| toward vs TauDEM | 0.723 | 0.948 | 0.933 | 0.980 | 0.983 | 0.985 | 0.980 | 0.993 |
| TauDEM vs RichDEM | 0.867 | 0.787 | 0.970 | 0.939 | 0.991 | 0.975 | 0.982 | 0.964 |
| Barnes vs MAS | 0.093 | 0.809 | 0.523 | 0.947 | 0.787 | 0.978 | 0.915 | 0.975 |
| Barnes vs DDM | 0.736 | 0.556 | 0.922 | 0.883 | 0.967 | 0.946 | 0.957 | 0.938 |

## Integrity checks (all 160 runs)

| Method | Runs with loops | Cells in 2-cell loops | Unrouted flat cells |
|---|---:|---:|---:|
| QEHT toward, Barnes, hybrids | 0 | 0 | 0 |
| TauDEM | 0 | 0 | 0 |
| RichDEM | 0 | 0 | 0 |
| MAS 1.2.1 | **160** | 152,406 | 0 |
| DDM 2.3 | 0 | 0 | 0 |

## Findings

1. **The hybrid adds nothing.** It reproduces Barnes: F1 is at least 0.99 everywhere, and exactly 1.000 outside flat ALOS terrain.
2. **Toward-lower and Barnes are not comparable.** Their spread (1 − F1) is 0.009–0.189. The spread between the reference tools is 0.004–0.130. The QEHT spread is larger in every one of the 8 strata. The gap is large on whole-metre flats and small in steep terrain.
3. **Barnes agrees best with the independent TauDEM.** Averaged over the 8 strata, F1 is 0.969 for Barnes and 0.948 for toward-lower. The largest gap is on ALOS flats: 0.870 vs 0.733. Toward-lower is slightly ahead on FABDEM flat (0.954 vs 0.940) and FABDEM mountainous (0.996 vs 0.990).
4. **The benchmark found two defects in QEHT 0.12. Both are fixed in 0.13.**
   - Flats that touched the grid edge or NoData had no outlet and stayed unrouted. On the first flat ALOS test area this left 75 cells unrouted. They now drain off the grid, as in TauDEM and RichDEM.
   - QEHT Barnes departed from the published algorithm at a flat's low edge. A cell there could first run along the rim of the flat, where the reference leaves the flat immediately. QEHT also scanned neighbours in a different order. Before the fix, 0.12 Barnes routed only 78 % of flat cells the same way as RichDEM. It now matches RichDEM cell for cell in all 160 runs. A per-cell port of RichDEM's code is the new test oracle, and QEHT matches it exactly on 80 random surfaces with NoData.
5. **The other tools behave differently by design, or have defects.**
   - MAS 1.2.1 lets no cell drain off the grid edge, and every run contains flow loops, so its basins are not usable.
   - DDM HydroLogic takes each cell's receiver from the order of a priority flood, not from steepest descent.
   - RichDEM's command-line D8 picks the lowest neighbour and ignores the diagonal distance. This is why its end-to-end results differ from TauDEM and QEHT on FABDEM, while the flat-isolated view agrees.
6. **The DEM matters as much as the method.** On FABDEM, which is floating point, flats are about 1 % of cells in every class and the choice of method changes little. On ALOS, which is whole metres, flats reach 53 % in flat terrain.

## Decision (Edmond, 29 Sep 2026)

- **Barnes 2014 is the default flat method.**
- **Hybrid is removed** from the UI and from `d8_direction`. `core/flow/flats.resolve_flats` keeps the `w` parameter for research and tests.
- **Toward-lower stays as a second option.** It is labelled for matching outputs of platforms that route flats that way. The enum index is unchanged, so saved models keep their method.
- **Still open:** re-run the Site B comparison against the commercial reference with the fixed Barnes code.

## Timing

Indicative only: the steps timed are not identical across tools. The figures are medians per 0.25 Mcell area.

| Method | What is timed | Median time |
|---|---|---:|
| QEHT D8 + flats | Direction on a filled DEM | 0.04 s |
| RichDEM | Fill + flats + file I/O | 0.11 s |
| TauDEM | Fill + D8 + file I/O + MPI start-up | 0.51 s |
| DDM | Receiver graph | 0.82 s |
| MAS | Fill + flats + D8, numba | 1.46 s |

## Validation of QGIS 4.2.2 (same session)

**QGIS 4.2.2-Belém do Pará**, run headless in a Debian trixie chroot with the official qgis.org packages (Qt 6, Python 3.13, GDAL 3.10.3, NumPy 2.2.4):

- All unit suites pass: test_core 47, test_interop 74, test_crossings 43, test_soils 31, test_flats 31.
- The QGIS smoke test passes all 28 checks.
- The plugin loads through `classFactory`, registers its 13 algorithms and unloads cleanly.
- The official `pyqt5_to_pyqt6.py` migration script suggested one cosmetic change, in a fallback branch that only runs on QGIS 3.22–3.36.

**QGIS 3.34.4-Prizren:** the same suites pass.

**Bandit 1.9.4:** clean.
