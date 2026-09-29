# QEHT — QGIS Engineering Hydrology Toolkit v0.9.0

Terrain and drainage analysis for QGIS, computed entirely in-process, offering
the same class of tools as commercial GIS hydrology extensions.

Author: Edmond Akello · License: GPL-2.0-or-later

## Why it exists

Standard QGIS hydrology routes through providers that spawn child processes:
`gdal:polygonize` shells out to `gdal_polygonize.py`, SAGA algorithms invoke
`saga_cmd.exe`, GRASS algorithms write and execute a batch file. On a managed
workstation with endpoint security, those calls are frequently blocked.

QEHT spawns nothing. Every computation uses NumPy and the GDAL Python bindings,
which are C++ libraries already loaded inside the QGIS process. There is no
`subprocess` import anywhere in the codebase, and no pip installation is required.

## Design

    qeht/
      core/                  ← NO QGIS IMPORTS. Pure NumPy + lazy GDAL.
        grid.py              neighbour conventions, standard D8 encode/decode
        raster.py            GDAL I/O, in-process gdal.Polygonize()
        conditioning/fill.py priority-flood depression filling
        flow/direction.py    D8 steepest descent, distance-weighted
        flow/accumulation.py topological accumulation + Strahler
        watershed/delineate.py  streams, snapping, catchments, longest path
        watershed/statistics.py morphometry, the four slope domains, 10-85 slope
        linking/ids.py       outlet_uid generation and uniqueness
        geometry/polygonize.py  cell mask -> polygon (pure NumPy)
        interop/             HEAS exchange: field dictionary, GeoPackage writer
      processing_provider/   ← the only place QGIS and core meet
      tests/test_core.py     ← runs on bare Python + NumPy
      tests/test_interop.py  ← outlet_uid, slopes, exchange package, golden fixture
      tests/qgis_smoke.py    ← every Processing tool, run inside QGIS

The core is importable without QGIS. That is what makes the hydrology testable:

    python -m qeht.tests.test_core        # 47 analytic checks
    python -m qeht.tests.test_interop     # 74 checks incl. the golden fixture

Run these after any change to the core, and before trusting any output on a
real project. `tests/qgis_smoke.py` needs a QGIS installation (see its header).

## Tools

| Tool | Equivalent to |
|---|---|
| Fill depressions | a raster toolset's `Fill` |
| D8 flow direction | a raster toolset's `Flow Direction` (D8) |
| Flow accumulation | a raster toolset's `Flow Accumulation` |
| Stream network and Strahler order | a hydrology toolset's `Stream Definition` + `Stream Order` |
| Delineate catchments from pour points | a hydrology toolset's `Watershed` / `Batch Watershed Delineation` |
| Stream network to polylines | a hydrology toolset's `Drainage Line Processing` |
| Longest flow path | a hydrology toolset's `Longest Flow Path` |
| Catchment and flow path characteristics | a hydrology toolset's `Basin/Longest Flow Path` attributes |
| Build HEAS exchange package | — (one self-describing GeoPackage for HEAS; see below) |

For Pairwise Intersect, use the built-in `native:intersection` — it is C++ and
never spawns anything. There is no reason to wrap it.

## Catchment characteristics

The "Catchment and flow path characteristics" tool delineates each catchment and
its longest flow path, and writes the morphometry a design flood calculation
needs.

Every catchment and its flow path carry **`outlet_uid`**, a text identifier
assigned once per run: taken from a pour-point attribute you choose (e.g.
existing culvert numbers), or sequential with your prefix (`X001`, `X002` ...,
`001` = most downstream). `outlet_id` (the pour-point feature id) is still
written but is not stable — it changes when the layer is edited or re-saved —
so never join on it.

**Catchment polygon fields**

| Field | Meaning |
|---|---|
| `outlet_uid` | stable identifier shared with the flow path |
| `area_km2` | contributing area |
| `elev_max_m` / `elev_min_m` | highest and lowest ground in the catchment |
| `relief_m` | elev_max - elev_min |
| `catch_slope_horn` | mean terrain gradient, Horn 3x3 (a standard slope algorithm) |
| `catch_relief_ratio` | relief / longest flow path length |
| `lfp_length_km` | longest flow path length |
| `slope_mean`, `slope_relief_ratio` | v0.8 names of the two slopes above, kept as aliases for this release |

**Flow path line fields**

| Field | Meaning |
|---|---|
| `outlet_uid` | stable identifier shared with the catchment |
| `lfp_length_km` | planimetric length |
| `lfp_elev_max_m` / `lfp_elev_min_m` | elevation at the divide and the outlet |
| `lfp_drop_m` | drop along the path |
| `lfp_slope` | drop / length over the whole path |
| `lfp_slope_1085` | slope between the points at 10% and 85% of the length, measured from the outlet |
| `lfp_L10_m`, `lfp_L85_m`, `lfp_z10_m`, `lfp_z85_m` | where those points are and their elevations, for a hand check |

### Four slope domains, never merged

QEHT reports four slopes under four names, and none substitutes for another:
`catch_slope_horn` (mean terrain gradient), `catch_relief_ratio` (relief ÷ LFP
length), `lfp_slope` (drop ÷ length along the LFP) and `lfp_slope_1085` (10–85
along the LFP).

### Two catchment slopes, deliberately

`catch_slope_horn` and `catch_relief_ratio` are both commonly called "catchment
slope" and they are not interchangeable. On the Site A catchments they differ
by a factor of 1.2 to 4.5:

| Outlet | catch_slope_horn | catch_relief_ratio | ratio |
|---|---|---|---|
| 22528 | 0.1769 | 0.0438 | 4.0x |
| 23238 | 0.0506 | 0.0387 | 1.3x |
| 24055 | 0.1762 | 0.0689 | 2.6x |

Runoff coefficient and curve number tables generally assume the mean terrain
gradient. Kirpich and most road drainage manuals want the relief ratio. Picking
the wrong one changes a design discharge substantially. **State which you used.**

`lfp_slope_1085` = (z85 − z10) / (L85 − L10), where L10 and L85 are 10% and
85% of the path length **measured from the outlet upstream**. It excludes the
flat bottom 10% and the steep top 15% of the path. Elevations are interpolated
linearly along the path at exactly those distances. It is required by several
UK and TRRL methods and is usually the more defensible design figure.

**Changed in 0.9.** QEHT 0.8.3 and earlier measured the 10% and 85% points from
the divide (i.e. 90% and 15% from the outlet) and took the next cell's elevation.
On concave profiles the conventional definition gives a lower slope: over 40
longest flow paths on a steep 30 m test area the 0.9 value was a
median 6% lower (range −20% to +5%) than the 0.8.3 value. Re-run QEHT before
reusing 10–85 slopes from ≤ 0.8.3 outputs.

### Supply the raw DEM for elevations

Routing uses the conditioned DEM; reported elevations should come from the raw
one, or heights inside filled depressions read as fill surface rather than
ground. The tool takes both and warns if the raw DEM is omitted.

## HEAS exchange package

**Build HEAS exchange package** writes one GeoPackage per run (schema
`qeht-heas-1`) that HEAS imports with no field mapping:

| Table | Content |
|---|---|
| `crossings` | snapped outlet points with `outlet_uid`, input and snapped coordinates, snap distance, contributing area, Strahler order |
| `catchments` | one polygon per crossing, with all catchment fields |
| `flowpaths` | longest flow path per crossing, with all flow-path fields |
| `qeht_run_metadata` | QEHT version, DEM path and SHA-256, CRS, cell size, thresholds, snapping, ID scheme, full parameters |
| `qeht_field_dictionary` | meaning, unit, method and HEAS target of every field |

The three layers share `outlet_uid`, so HEAS links them by key rather than by
row order or by position. Rules: a projected, metric CRS is required (a DEM in
degrees is refused); duplicate or empty IDs stop the run with a list; missing
values are written as NULL, never 0; schema changes are additive only (a
rename or removal would be `qeht-heas-2`). An optional CSV per layer carries
`outlet_uid` for spreadsheet users. The GeoPackage is written with the Python
standard library and checked against GDAL's GeoPackage validator in the tests;
`tests/fixtures/golden_exchange.gpkg` is the shared fixture that HEAS imports
in its own tests.

## Interoperability

Flow direction is written and read in the standard D8 encoding
(E=1, SE=2, S=4, SW=8, W=16, NW=32, N=64, NE=128; 0 = sink or edge outlet).

This works in both directions. You can feed a flow direction grid produced by
another GIS/hydrology toolset straight into QEHT's accumulation tool — which
is the cleanest way to validate this toolkit against a trusted production
result.

## Known divergences from commercial reference toolsets

Stated explicitly because a drainage report needs them stated:

1. **Flat resolution is one-sided.** We implement the gradient-toward-lower-
   terrain component of Garbrecht & Martz (1997), not the gradient-away-from-
   higher-terrain component. Drainage across wide flats converges into fewer
   channels than commercial reference toolsets typically produce. Mitigate by
   filling with a small minimum slope (1e-4), which largely removes flats
   before routing.
2. **Tie-breaking.** Where two neighbours give identical drop/distance, we take
   the lowest internal index (E before SE before S...). Commercial toolsets'
   exact tie handling is generally undocumented. Ours is deterministic and
   reproducible, which matters more.
3. **No tiling.** DEMs are processed whole in memory. Beyond roughly
   5000 x 5000 cells, clip to your catchment plus a buffer first — which is
   normal practice anyway.

## Validation

**Level 1 — synthetic analytic DEMs.** 47 checks, `tests/test_core.py`, plus 74 interop checks in `tests/test_interop.py`. PASS.

**Level 3 — reference hydrology toolset production output, Site A.** 718 x 775
cells @ 30.92 m, EPSG:21037, 16 road-crossing pour points with reference
catchments and longest flow paths.

| Quantity | Result |
|---|---|
| Catchment area, 16 outlets | median ratio **1.009**, 15/16 within ±30% |
| Spatial agreement (IoU) | median **0.84** |
| Longest flow path length | median ratio **0.999**, 15/16 within ±20% |
| Best individual match | 3505.7 m (reference) vs 3503.9 m (QEHT) — 0.05% |
| Flow-direction cycles | 0 |

Two defects were found by this validation and fixed:

1. **Pour-point snapping.** Snapping to maximum accumulation within the radius
   pulled 5 of 16 road crossings off their small tributary onto the adjacent
   trunk stream — one catchment inflated 40x. Fixed by snapping to the NEAREST
   cell of the extracted stream network instead. Set the snap threshold low
   enough to include the tributary the crossing sits on; a threshold of 500
   cells re-broke the case that 200 fixed.
2. **Undeclared NoData.** The Site A DEM had 11.8% of cells at 0.0 with no
   NoData set in the header, while real terrain starts at 1574 m. `audit_nodata()`
   now detects this and warns; the Fill algorithm takes a NoData override.

The one remaining outlier (0.65 area ratio) is a genuine hydrological
difference, not a snapping artefact, and is worth inspecting on its own.

### Two validation datasets, two regimes

QEHT has been checked cell-by-cell against a reference hydrology toolset on
two real DEMs, and they bound the expected accuracy:

| Dataset | Terrain | Flat cells | Direction agreement |
|---|---|---|---|
| Site A | steep, Rift floor | ~5% | **98.3%** (aligned grid) |
| Site B | flat, coastal | **34%** | **87.9%** (same-grid, G&M flats, v0.5 tie-break) |

The difference is flats, not resolution - both are ~30 m. Where a third of the
surface fills to a flat, the flat-resolution algorithm rather than steepest
descent decides most directions, and that is where QEHT and the reference
toolset diverge most. On the flat-heavy Site B DEM:

  * Garbrecht & Martz flats: **85.5%** overall, 81.6% on flats, 87.6% off flats
  * Toward-lower-only flats: 81.7% overall, 76.2% on flats

**Barnes convergent flat resolver, built and benchmarked (v0.8).** The iterated
Barnes 2014 two-gradient resolver was completed and verified:

  * Passes a synthetic saddle test - interior converges to the single outlet,
    no false sinks, accumulation increases monotonically downstream.
  * On real terrain resolves ~98% of flat cells in 2 iterations, leaving only
    genuinely outletless cells unrouted; no cycles, no uphill flow.

But on the flat Site B coastal DEM it DISPERSES flow more than the simple method:

| Method | IoU@100 | IoU@1000 | p99 accum (reference 9,071) | cells>1000 (reference 56,119) |
|---|---|---|---|---|
| Toward lower | 0.69 | 0.63 | 6,324 | 51,974 |
| Barnes iterated | 0.50 | 0.19 | 642 | 14,428 |

The conclusion, now firmly evidenced: the reference toolset's flat behaviour -
despite its documentation citing Garbrecht & Martz - empirically resembles the
simpler toward-lower method on this terrain, not the canonical convergent
Barnes result. Both Barnes and toward-lower are legitimate, cycle-free
drainage solutions; they simply differ, and toward-lower is closer to the
reference toolset here. Toward-lower remains the default. Barnes is
selectable for terrain where a strictly convergent solution is wanted. On
steep Site A the two agree to within 0.1% (98.5 vs 98.4).

This closes the flat-resolution investigation: the remaining gap between QEHT
and the reference toolset on flat terrain is a genuine algorithmic
difference, not an implementation bug.

**Accumulation is validated to near-exact agreement (v0.7.1).** The single most
diagnostic test: feed the reference toolset's OWN flow-direction grid into
QEHT's accumulation engine (same grid, no resampling, no QEHT routing).
Against that toolset's own accumulation raster over 9,000,000 Site B cells:

| Metric | Result |
|---|---|
| Exact cell match | 99.31% |
| Within 1 | 99.33% |
| Cycles | 0 |
| Stream IoU @100 / @1000 / @5000 | 0.993 / 0.976 / 0.946 |

This proves the accumulation engine is correct and isolates ALL remaining
flat-terrain discrepancy to flow *direction* - specifically flat routing. It
means a correct convergent flat resolver would close almost the entire gap,
and that nothing needs changing in accumulation, thresholding, or Strahler code.

**Flat resolution (v0.7).** Validated against a reference hydrology toolset's
flow *accumulation* (not just direction) on the 34%-flat Site B coastal DEM,
the gradient-toward-lower method reproduces the channel network well:

| Threshold | Stream IoU | Reference cells | QEHT cells |
|---|---|---|---|
| accum > 100 | 0.67 | 670,278 | 642,721 |
| accum > 1000 | 0.60 | 220,256 | 210,325 |
| accum > 5000 | 0.56 | 109,005 | 101,866 |

Accumulation percentiles track the reference toolset closely (p99: 4,438 vs
5,810; cells above 500: 19,457 vs 21,008). This corrects an earlier
over-emphasis on raw per-cell direction agreement, which is misleading on
flats - cells can disagree on direction while producing near-identical
networks, because a flat can be crossed by several equivalent paths.

A faithful port of RichDEM's Barnes 2014 two-gradient flat router was trialled
to close the residual gap. Across three implementations it produced
hydrologically valid but DISPERSED flow (every cell reaches an outlet, no loops,
no uphill drainage - but accumulation spreads instead of concentrating, stream
IoU 0.04). Matching the reference toolset exactly would require reproducing
RichDEM's entire integrated fill-and-route pipeline including its
epsilon-elevation output stage, a large rewrite with little benefit given how
well the simpler method already tracks the reference toolset on accumulation.
The Barnes code remains in the tree (`core/flow/flats.py`) for reference but
is not wired into the tools.

**Tie-breaking (v0.5).** The residual off-flat disagreement was traced to
tie-breaking: on low-relief terrain, 21% of cells have two or more neighbours at
an equal steepest drop/distance, and 99.8% of off-flat disagreements were these
ties. The reference toolset resolves them by a fixed priority (among
cardinals, S > W > N > E, diagonals below all cardinals), recovered
empirically from the Site B grid.
Applying it lifted off-flat agreement from ~84% to ~91% across four independent
windows, with no regression on steep Site A terrain (98.3% -> 98.4%). It is now
the default and applies on all terrain. For low-lying or coastal catchments, treat
the reference toolset as the benchmark and QEHT as a cross-check, not the
other way round.

The measurement itself is grid-sensitive: comparing QEHT (run on the raw DEM) to
the reference toolset across a half-cell grid offset dropped apparent
agreement to 76%. Routing the reference toolset's OWN filled DEM - identical
grid, no resampling - recovered it to 85.5%. Always compare on an aligned
grid before drawing conclusions.

### Cell-by-cell against a reference flow direction grid

The Site A DEM sits inside a larger reference-toolset run on the same grid
(EPSG:21037, 30.9187 m, aligned to the cell). Comparing the 490,148
overlapping cells:

**98.29% identical flow direction.**

Disagreement is not uniform. It concentrates almost entirely in filled
depressions:

| Location | Disagreement rate |
|---|---|
| Inside filled depressions (2.2% of DEM) | **24.03%** |
| Everywhere else | 1.20% |

That 20x concentration is the source of the near-parallel streaks visible when
streams are drawn across flat filled areas: with a single toward-lower gradient,
every cell in a flat inherits nearly the same direction and the flat drains as a
sheet of parallel lines rather than a converging network.

### Flat resolution is a choice, not a default

Measured on Site A against the reference toolset's grid and the 16 reference
catchments:

| Method | Dir. agreement | Catchment area median | Within ±30% |
|---|---|---|---|
| Garbrecht & Martz (two-gradient) | 98.29% | 0.989 | 14/16 |
| Toward lower terrain only | — | 1.009 | 15/16 |
| No flat resolution | 99.26% | 0.005 | 0/16 |

No flat resolution is not an option: 24,938 unrouted cells sever the drainage
network and catchment delineation collapses entirely, even though the cells that
DO get a direction agree slightly better.

Between the two real methods the evidence is mixed — G&M is the theoretically
correct method and should reduce the streak artefact, while the one-sided method
scored marginally better on catchment area on this particular DEM. It is exposed
as a parameter on the Flow Direction algorithm. **Run both on your own DEM and
compare before committing to one for a deliverable.**

**Level 4 — cross-check against an independent reference implementation.** The
reference tool uses bit-identical neighbour ordering, standard D8 codes and
tie-breaking, so flow direction should match cell-for-cell on sloping terrain.
Note the reference tool initialises accumulation to the cell's own weight, so
`reference == QEHT + 1`. Harness supplied separately.

## Licence

GNU General Public License v2 or later.
