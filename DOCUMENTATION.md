# QEHT — QGIS Engineering Hydrology Toolkit

## Technical Documentation

**Version:** 0.18.0
**Type:** QGIS Processing plugin for DEM-based terrain and drainage analysis
**Licence:** GNU General Public License v2 or later
**Implementation:** Python, NumPy, GDAL Python bindings, QGIS Processing API

---

## 1. Purpose and design philosophy

QEHT is a QGIS Processing plugin that takes a digital elevation model (DEM) and produces the terrain and drainage quantities an engineering hydrology study needs: conditioned elevations, flow direction, flow accumulation, stream networks, catchments, longest flow paths, and the morphometric characteristics that feed a design-flood calculation.

Three commitments shape the whole design.

**In-process computation, no subprocesses.** Every calculation runs inside the QGIS Python process using only NumPy and the GDAL Python bindings that ship with QGIS. QEHT never spawns a child process — it does not call `gdal_polygonize.exe`, `saga_cmd.exe`, or the GRASS batch wrappers, and it requires no `pip` installation. This is a practical requirement, not an aesthetic one: on managed or endpoint-secured workstations, GUI applications that spawn hidden child processes are frequently blocked or flagged, and the standard QGIS hydrology providers (GDAL scripts, SAGA, GRASS) all route through exactly that mechanism. QEHT was built to work where those tools cannot.

**Separation of the numerical core from the GIS layer.** The `qeht.core` package contains no QGIS imports. It operates on NumPy arrays and file paths, and can be imported and tested from a bare Python interpreter with only NumPy present. The QGIS Processing layer (`qeht.processing_provider`) is the only place where QGIS and the core meet: it resolves parameters, reads and writes rasters through GDAL, and adapts progress reporting. This boundary is what makes the hydrology independently testable and what allows QEHT to be compared against a reference GIS platform's hydrology toolset without conflating interface behaviour with numerical algorithms.

**Explicit, documentable numerical decisions.** DEM conditioning is a visible step, not a side effect of routing. Flow direction is a separate algorithm from accumulation. "Catchment slope" is reported under two distinct, named definitions rather than a single ambiguous value. Every choice that affects a design discharge is exposed and recorded, because the output is intended to be defensible in an engineering review.

---

## 2. Architecture

```
qeht/
  __init__.py            plugin entry point (classFactory)
  plugin.py              registers the Processing provider; nothing else
  metadata.txt           QGIS plugin metadata
  icon.svg

  core/                  NO QGIS IMPORTS - pure NumPy + lazy GDAL
    grid.py              neighbour conventions, D8 encode/decode, receivers
    raster.py            GDAL I/O, in-process gdal.Polygonize(), NoData audit
    conditioning/
      fill.py            priority-flood depression filling
    flow/
      direction.py       D8 steepest descent, distance-weighted, tie rule
      flats.py           flat resolution (toward-lower and Barnes 2014)
      accumulation.py    topological accumulation, Strahler ordering
      streamlines.py     stream-network vectorisation into reaches
    watershed/
      delineate.py       streams, snapping, catchments, longest flow path
      statistics.py      catchment and flow-path morphometry, 10-85 slope
      morphometry.py     overland/channel split of the LFP, basin shape and network indices
      channel.py         approach / exit channel slopes, crossing deposition indicator
    linking/
      ids.py             outlet_uid generation and uniqueness
      relink.py          renumber_log for renumber-and-relink (D3)
    geometry/
      polygonize.py      cell mask -> OGC-valid (multi)polygon, pure NumPy
    soils/
      usle_k.py          Williams/EPIC and Dg-based K, CFRG, USDA texture, HSG proxy (scalar + vectorised)
      sotwis.py          SOTWIS SQLite loader, generic attribute loader, CSV table loader
      sources.py         SoilGrid: any soil source on one per-cell model; SoilGrids, HYSOGs
      catchment.py       per-catchment soil block (field list)
    runoff/
      curve_number.py    curve number (TR-55 x WorldCover proxy, user lookups), AMC, Rational C
      rainfall.py        rainfall zones, mean annual rainfall, R estimate from rainfall
    geometry/
      rasterize.py       polygons -> grid, cell-centre rule
    erosion/
      terrain.py         A_s, SPI, ln SPI, TWI, LS (Moore & Burch, Desmet & Govers)
      rusle.py           source-tagged R, K, C, P; WorldCover C/P lookups; LS-only fallback
      classes.py         class schemes, severity scores, combination matrix, class rasters
      summary.py         catchment and crossing erosion blocks, SDR, impact score
      corridor.py        stations and reaches along an alignment
      io.py              the erosion output folder and its run record
    network/
      alignment.py       linear referencing: chainage, signed offset, intersections
      crossings.py       road x drainage candidates, parallel reaches, clusters
      profile.py         alignment ground profile
      coverage.py        missing crossings, sag points (walled routing), flat stretches
      floodplain.py      floodplain width at large crossings (profile, HAND)
    interop/
      field_dictionary.py  the qeht-heas-1 contract (every exchange field)
      gpkg.py            standard-library GeoPackage 1.2 writer/reader
      heas_exchange.py   pipeline, writer, validator, CSV export
    report/
      characteristics.py one row per crossing from a package (the characteristics table)
      run_report.py      HTML run report from a package (nothing recalculated)
      quicklooks.py      PNG + world file + legend for rasters (pure-Python PNG writer)

  processing_provider/   the only QGIS-aware code
    provider.py          registers the nineteen algorithms
    base.py              shared base class and helpers (pour points, IDs)
    alg_*.py             one file per algorithm

  tests/
    test_core.py         47 analytic checks, runnable without QGIS
    test_interop.py      77 checks: outlet_uid, slopes, exchange, golden fixture
    test_crossings.py    45 checks: alignment, candidates, clusters, burn, relink
    test_soils.py        33 checks: USLE K, texture, HSG, SOTWIS and CSV loaders, soil block
    test_flats.py        36 checks: v0.8.3 oracles, Barnes == RichDEM port, edge drains, DEM QA
    test_erosion.py      50 checks: analytic plane/valley, A14 2025 anchors, RUSLE, classes, corridor
    test_alignment.py    15 checks: ground profile on a plane and a valley
    test_morphometry.py  22 checks: LFP split, perimeter/shape ratios, drainage density, links
    test_soils_any.py    18 checks: vectorised = scalar, v0.14 regression, SoilGrids, HYSOGs, overrides
    test_runoff.py       17 checks: CN mosaic, dual groups, AMC, lookups, pipeline
    test_channel.py      16 checks: STI identity, channel slopes, deposition ratio
    test_coverage.py     21 checks: missing crossings, sags, walled sag areas, flat stretches
    test_rainfall.py     35 checks: zone shares and ties, rainfall mean and units, R relations
    test_report.py       15 checks: characteristics table, run report
    test_floodplain.py   19 checks: analytic valleys, truncation, HAND on a V valley
    test_quicklooks.py   14 checks: PNG round trip, world file, legends, downsampling
    fixtures/            golden_exchange.gpkg + .json (shared with HEAS)
    qgis_smoke.py        every Processing tool run inside QGIS
```

The dependency direction is strict and one-way: `processing_provider` imports `core`; `core` never imports `processing_provider` or `qgis`. GDAL imports inside `core` are function-local and lazy, so the numeric modules remain importable where GDAL is absent.

---

## 3. Processing algorithms

QEHT registers nineteen algorithms under the "Engineering Hydrology" provider.

| Algorithm | Commercial reference analogue |
|---|---|
| Fill depressions | a reference raster toolset's `Fill` |
| D8 flow direction | a reference raster toolset's `Flow Direction` (D8) |
| Flow accumulation | a reference raster toolset's `Flow Accumulation` |
| Stream network and Strahler order | a reference hydrology toolset's `Stream Definition` + `Stream Order` |
| Stream network to polylines | a reference hydrology toolset's `Drainage Line Processing` |
| Delineate catchments from pour points | a reference hydrology toolset's `Watershed` / `Batch Watershed Delineation` |
| Longest flow path | a reference hydrology toolset's `Longest Flow Path` |
| Catchment and flow path characteristics | a reference hydrology toolset's basin/LFP attribute tools |
| Build design hydrology package (formerly "Build HEAS exchange package") | none — one self-describing GeoPackage of crossings, catchments and flow paths (Section 5.1) |
| Renumber and relink exchange package | none — gapless re-issue of IDs after edits (Section 5.2) |
| Road crossing candidates | none — road × drainage crossings (Section 4.11) |
| Burn crossings through embankments | DEM reconditioning at culverts (Section 4.12) |
| Soil, rainfall and runoff parameters for catchments | zonal soil statistics, USLE K, HSG shares, curve number, rainfall (Sections 4.13, 4.17, 4.18, 4.21) |
| Erosion indices and RUSLE soil loss | terrain indices, RUSLE and severity classes (Section 4.14) |
| Sample erosion along alignment | none — erosion stations and reaches along a road (Section 4.14) |
| Alignment ground profile | none — ground, fill, area and stream crossings along a road (Section 4.15) |
| Drainage coverage check along a road | none — missing crossings, sags, flat stretches (Section 4.20) |
| Run hydrology pipeline (one click) | none — the whole chain, one output folder, run report (Section 4.22) |
| Run report from a design hydrology package | none — HTML report and characteristics table (Section 4.22) |

Each is a `QgsProcessingAlgorithm` registered through a `QgsProcessingProvider`. Exposing the tools this way — rather than as bespoke dialogs — means they gain input validation, batch mode, the Graphical Modeler, the history log, and `processing.run()` scriptability at no additional cost. Chaining tools in the Modeler is much of a commercial hydrology extension's practical value, and this design reproduces it.

---

## 4. Algorithms and equations

### 4.1 Depression filling

Filling uses the priority-flood algorithm (Barnes, Lehman & Mulla 2014): the DEM is flooded inward from its boundary using a priority queue keyed on elevation, and each cell is raised to the maximum of its own elevation and the spill level at which water reaches it. Complexity is O(N log N), a single pass, in contrast to the O(N²) behaviour of iterative neighbourhood-scanning fills.

Conditioning is deliberately a standalone step producing only a conditioned elevation surface — it does not emit flow directions. An optional minimum slope may be imposed across filled flats (an epsilon increment per traversal step) to give routing a defined gradient. An optional fill-depth raster is produced as a QA product: large contiguous fill depths typically mark a road embankment, dam, or culvert the DEM treats as a barrier, and signal that breaching or stream burn-in may be more appropriate than filling.

**NoData auditing.** QEHT inspects the DEM for undeclared NoData — for example a large block of cells at exactly 0.0 while genuine terrain begins far above it, with no NoData value set in the header. Treated as valid terrain, such a block becomes a spurious flat sink that the whole catchment drains into. The audit warns; the Fill algorithm accepts a NoData override so the user can mask it. (This was not hypothetical: it was found on a real project DEM with 11.8% of cells at 0.0 and real terrain starting at 1574 m.)

**Resampling audit (v0.13.1).** A DEM regridded or reprojected with nearest neighbour repeats whole source rows and columns wherever the output grid is finer than the source. The repeats are exact, so they create artificial flats and exactly tied neighbours, and flow routing across them is decided by the flat method rather than the terrain. Fill and D8 flow direction count rows and columns that repeat their neighbour over ≥ 99 % of their overlapping valid cells (constant rows excluded) and warn above 0.2 %. Calibration: nearest-neighbour reprojections showed 0.45–3 %; native and bilinear grids 0 %; whole-metre DEMs (ALOS) have many tied cells but no repeated rows and are not flagged. Note that a FABDEM export at a scale slightly finer than its native 1″ grid (e.g. an Earth Engine export at 0.000269°) repeats about one row and column in 33 — export in the native projection and scale, or resample bilinearly.

**Flat-method sensitivity (v0.13.1).** Catchment characteristics routes the conditioned DEM with both flat methods and reports each outlet's contributing area under each (`area_barnes_km2`, `area_toward_km2`, `flat_sensitivity_pct`); `flat_sensitive = 1` above a tolerance (default 10 %). A flagged catchment is decided by flats or ties, not terrain, and needs checking against mapped drainage or on site. On Site C (a road project) the check flagged exactly the four crossings on the section built from a nearest-neighbour DEM (areas differing by 40–100 %) and none of the six stable ones; on a bilinear resample it flagged none.

### 4.2 D8 flow direction

Flow direction selects the neighbour of steepest descent, weighted by distance:

```
S_k = (z_i - z_k) / d_k
```

with `d_k = cellsize` for orthogonal neighbours and `d_k = cellsize·√2` for diagonals. The neighbour with the largest positive `S_k` is chosen. Distance weighting matters: choosing the largest absolute drop instead of the largest drop-per-distance biases routing toward diagonals and changes both flow paths and contributing areas.

Directions are stored internally as indices 0–7 and written in the standard D8 encoding (E=1, SE=2, S=4, SW=8, W=16, NW=32, N=64, NE=128; 0 for sinks and edge outlets). The encoder and decoder are both provided, so a flow-direction grid interchanges with a reference GIS platform's hydrology tools in **both** directions — which is also how QEHT is cross-validated against a trusted production grid.

### 4.3 Tie-breaking

Low-relief DEMs contain many cells with two or more neighbours at equal drop-per-distance. Empirically, on such terrain, this accounts for essentially all off-flat disagreement with the reference toolset. QEHT applies a fixed tie priority, recovered empirically from a reference hydrology toolset's grid:

```
S > W > N > E > SE > SW > NW > NE
```

with a relative tolerance of 1e-6 for treating gradients as tied. Applying this rule raised off-flat agreement with the reference toolset from ~84% to ~91% across four independent windows, with no regression on steep terrain. This is described as an *empirically recovered* interoperability rule, not a proof of the reference platform's undocumented internal behaviour.

### 4.4 Flat resolution

After steepest-descent routing, cells in filled flats have no downhill neighbour. Since v0.13 QEHT provides two methods; Barnes is the default (WP-G benchmark, see `BENCHMARK.md`). Flats that touch the grid edge or NoData drain to it (boundary cells without descent are outlets), as in TauDEM and RichDEM.

**Method A — toward-lower BFS (option; default up to v0.12).** Flat cells are seeded from their outlets and a breadth-first distance is propagated inward; each flat cell routes toward decreasing distance-to-outlet. This is the toward-lower component of the Garbrecht & Martz conceptual approach.

**Method B — Barnes 2014 convergent resolver (default since v0.13).** Barnes' two-gradient method, reproducing his own RichDEM implementation cell for cell (tested against a per-cell port of RichDEM's code and against RichDEM itself on 160 benchmark runs): label each flat from its low edges; build a gradient *away* from high edges and a gradient *toward* low edges; combine them as `flat_mask = 2·toward + (flat_height − away)`; route each flat cell to the same-flat neighbour of lowest combined value; cells on a flat's low edge leave the flat directly (in RichDEM the routed low-edge cells carry the lowest value of all); neighbours are scanned W, NW, N, NE, E, SE, S, SW with cardinal preferred on ties. QEHT adds iteration to convergence — repeating the pass so nested flats drain in hierarchy order, since a flat whose outlet is itself a lower flat can only resolve once that lower flat has drained. The method passes synthetic saddle tests (interior converges to a single outlet, no false sinks, monotonic downstream accumulation) and resolves ~98% of flat cells in two iterations on real terrain.

**History: why toward-lower was the default up to v0.12.** On a flat coastal DEM (34% of cells flat), the canonical Barnes result *disperses* flow and matches the reference toolset's accumulation **less** well than the simpler toward-lower method:

| Method | Stream IoU @1000 | p99 accum (reference 9,071) | cells > 1000 (reference 56,119) |
|---|---|---|---|
| Toward-lower | 0.63 | 6,324 | 51,974 |
| Barnes iterated | 0.19 | 642 | 14,428 |

**Correction (v0.12).** The Barnes figures above were produced with the ≤ 0.11 code, whose outlet rule for local-minimum cells could point two cells at each other (it accepted neighbours assigned earlier in the same pass). The resulting two-cell loops sit on main channels and cut them off from their upstream area — consistent with the collapsed p99 accumulation (642 against 9,071). On a 30 m flat coastal clip the fix raised the largest accumulation from 40,117 to 794,103 cells. The comparison must be repeated with the fixed resolver before any conclusion about Barnes against the reference platform is drawn. Toward-lower was never affected.

**Removed in v0.13 — hybrid.** The WP-G benchmark showed it reproduces Barnes (stream F1 ≥ 0.99 in every stratum); it is no longer offered. `core/flow/flats.resolve_flats(w=…)` keeps the parameter for research. Original description (v0.12): `flat_mask = w·toward + (flat_height − away)`, w > 1. w = 2 is Barnes; as w → ∞ the away term only breaks ties and every cell follows a shortest path to its outlet (tested). A cell only flows to a neighbour with strictly lower flat_mask, so no loop is possible for any w; for w > 1 every cell with toward distance t > 1 has a neighbour at t − 1 whose value is lower by at least w − 1, so only low-edge cells can be local minima, and those exit to lower or pre-routed neighbours. For w ≤ 1 interior false sinks are possible and the value is refused. On a broad coastal flat w changes little (stream IoU vs toward-lower 0.394–0.399 for w = 1.5 … 10⁶, ~44 % of flat cells routed differently from toward-lower at every w): among equally short paths, Barnes' tie rules (same flat, cardinal first) and toward-lower's (first in E, SE, S … order) differ, and that choice dominates on large flats. An optional size switch applies the toward gradient alone to flats below N cells (no seams, tested); with the vectorised code it saves no measurable time.

**WP-G benchmark (v0.13).** 80 stratified random 15 km areas across Kenya (flat, rolling, hilly, mountainous; 20 each), each on ALOS AW3D30 and FABDEM, compared with TauDEM, RichDEM, MAS 1.2.1 and DDM HydroLogic 2.3. Toward-lower and Barnes are not interchangeable (their difference exceeds the TauDEM–RichDEM difference in every stratum); Barnes agrees best with the independent TauDEM (mean stream F1 0.969 vs 0.948; 0.870 vs 0.733 on whole-metre flats). Full design, tables and decision in `BENCHMARK.md`.

**Implementation (v0.12).** The resolvers are vectorised: flat edges by neighbour shifts, flats labelled by union-find with pointer jumping over equal-elevation links, BFS by whole frontiers over a compact adjacency list, assignment in the same neighbour order and tie rules. Toward-lower reproduces v0.8.3 bit for bit; Barnes reproduces it everywhere except the corrected exit cells (tests on random surfaces with plateaus, terraces, filled noise and NoData).

### 4.5 Flow accumulation

Accumulation is a topological traversal of the receiver graph. Each cell's in-degree (number of upstream contributors) is counted; cells of in-degree zero are queued; each cell is popped once, its load passed to its receiver, and the receiver's in-degree decremented and enqueued when it reaches zero. This is O(N) after flow direction is established and uses no recursion (a recursive downstream walk overflows the stack on any real DEM).

By default accumulation follows the standard D8 convention, excluding the cell's own contribution (a ridge cell reads 0). A weight raster may be supplied for contributing area or weighted accumulation. Cells trapped in a cycle — which a valid drainage tree never contains — are reported explicitly rather than silently dropped.

**Validation of the engine in isolation.** Fed the reference platform's own flow-direction grid, QEHT's accumulation reproduces its accumulation raster to **99.31% exact cell agreement** over 9,000,000 cells, with stream IoU of 0.99 / 0.98 / 0.95 at thresholds of 100 / 1000 / 5000 and zero cycles. This isolates the accumulation, thresholding, and Strahler code as correct, and localises all flat-terrain discrepancy to flow direction alone.

### 4.6 Strahler ordering

Strahler order is computed on the *extracted stream network*, not on every DEM cell — these are different quantities and are routinely conflated. A stream cell with no upstream stream cells is order 1; where the maximum upstream order occurs once, it is retained; where two or more equal maxima meet, the order increments. Non-stream cells are 0.

### 4.7 Stream extraction and vectorisation

Streams are extracted by an accumulation threshold, expressible as a cell count or as a contributing area. Area thresholds are preferable in engineering documentation because they are physically interpretable and survive a change of DEM resolution. Vectorisation traces the raster network into polyline **reaches** — one feature per segment between sources, junctions, and outlets — each carrying its Strahler order, length, elevation drop, and slope. This is the same concept as a commercial hydrology toolset's "Drainage Line", not a cell-by-cell polygonisation.

### 4.8 Catchment delineation and pour-point snapping

Catchments are delineated by reverse traversal of the receiver graph from each pour point, in O(N). Pour points are snapped to the drainage network first, because a surveyed culvert or outfall coordinate rarely lands exactly on the DEM-derived channel. Snapping moves each point to the nearest stream cell (ties broken by accumulation), which preserves the tributary it sits on. Validation against a reference hydrology toolset showed that naive snap-to-maximum-accumulation pulls road crossings onto the adjacent trunk stream and inflates catchments by up to 40×; nearest-stream snapping fixed 15 of 16 test catchments.

### 4.9 Longest flow path

The longest flow path is found by propagating cumulative downstream distance upward through the receiver tree, taking the maximum, and walking back down. Distances use true diagonal weighting, so the result is a real planimetric length rather than a cell count. One path is produced per pour point, each bounded by its own catchment, with an option for non-overlapping (local) catchments versus full upstream areas. Length and slope along this path are the direct inputs to the Kirpich, Bransby-Williams, and TRRL time-of-concentration methods.

### 4.10 Catchment and flow-path characteristics

For each catchment the tool reports area, highest/lowest/mean elevation, relief, and two distinct slopes; for each longest flow path it reports length, endpoint elevations, drop, whole-path slope, and 10–85 slope. Both carry `outlet_uid` (Section 5.1).

**Four slope domains, never merged.** `catch_slope_horn` is the mean terrain gradient over every cell by Horn's 3×3 method (the algorithm behind most commercial GIS `Slope` tools), which runoff-coefficient and curve-number tables assume. `catch_relief_ratio` is relief divided by longest-flow-path length. On the Site A catchments these two differ by factors of 1.2 to 4.5; they are not interchangeable and the report must state which was used. Along the flow path, `lfp_slope` is the drop over the whole length and `lfp_slope_1085` the 10–85 slope. Time-of-concentration methods use the flow-path slopes; `catch_relief_ratio` divides by LFP length (a variant of Schumm's basin-length relief ratio) and is reported for morphometry, not as a Tc input. `lfp_drop_m` is the highest minus the lowest raw-DEM elevation on the path; `lfp_z_head_m` / `lfp_z_outlet_m` give the end points and `lfp_nonmonotonic` flags a spike or pit (> 0.5 m) that makes the two differ. The v0.8 names `slope_mean` and `slope_relief_ratio` are written as aliases for one release.

**10–85 slope.** S₁₀₋₈₅ = (z₈₅ − z₁₀) / (L₈₅ − L₁₀), with L₁₀ = 0.10 L and L₈₅ = 0.85 L measured from the outlet upstream along the path, and z interpolated linearly between cell centres at exactly those distances. `lfp_L10_m`, `lfp_L85_m`, `lfp_z10_m` and `lfp_z85_m` are exported so the value can be checked by hand. Up to v0.8.3 the points were measured from the divide (90 % and 15 % from the outlet) with next-cell elevations; on concave profiles that gives a steeper value (median +6 % over 40 longest flow paths on a steep 30 m test area, range −5 % to +24 % relative to the 0.9 value). The test suite includes an analytic concave profile z = a·d² on which the two conventions give exactly 0.95·aL and 1.05·aL.

**Polygons.** Each catchment is traced from its own cell mask (`core/geometry/polygonize.py`), so overlapping full-upstream catchments each keep their full area. The tracer's output is OGC-valid and identical to `gdal.Polygonize` followed by a union (checked on ~3,000 random masks and on real 30 m catchments).

### 4.11 Road crossing candidates

**Alignment.** The road centreline is linearly referenced: chainage = start chainage + planimetric distance along the parts in layer order (gaps between parts not counted; optional reverse). Offsets are signed, positive to the left of increasing chainage.

**Intersections.** For every stream cell within the corridor, the D8 link from its centre to its receiver's centre is intersected with the alignment segments (parameter conventions make a link through a shared vertex count once). The D8 network is continuous, so this finds every crossing that intersecting QEHT's stream polylines would find, with no near-miss problem. Per candidate: crossing point and chainage; `crossing_angle_deg` = the acute angle between the flow link and the road (90 = square); `acc_km2` = (accumulation + 1) × cell area at the upstream cell of the link, which is also the outlet cell used later; Strahler order; reach id; `side_in`/`side_out` (L/R of the chainage direction); road azimuth.

**Parallel flow and clusters.** Stream cells with |offset| ≤ corridor half-width are traced downstream into runs (a run stops where it leaves the corridor or joins an earlier run). A connected set of runs with at least one run ≥ the minimum parallel length is a parallel system: its runs are exported as `parallel_reaches`, and every candidate touching it joins one cluster. Candidates closer than the merge distance in chainage are also joined (union-find). In each cluster the candidate with the largest contributing area is `recommended = 1` (decision D3: most downstream).

### 4.12 Burning crossings through embankments

At each crossing the breach runs perpendicular to the road between two points `half_length` either side of the centreline. Each end moves to the lowest valid cell within the search radius **on its own side of the centreline**; the higher end is upstream. Cells on the 8-connected line between the ends are lowered to z = min(z, grade), where the grade falls linearly from the upstream to the downstream end elevation (at least `min_grade` = 1e-4 so the breach drains). Nothing is raised and nothing outside the breach changes. A crossing where no ground stands above the grade (no embankment in the DEM) is reported with 0 cells. On a synthetic valley dammed by a 5 m embankment, burning removes the pond completely and the flow passes through the culvert cell.

### 4.13 Soil parameters and USLE K

**Data.** SOTWIS: `SOTERunitComposition` (NEWSUID → up to ten profiles PRID with shares), `SOTERparameterEstimates` (per profile and depth layer: sand SDTO, silt STPC, clay CLPC %, organic carbon TOTC g/kg, bulk density, coarse fragments vol %, FAO drainage), read with the standard-library `sqlite3`. Per profile, properties are depth-weighted by the overlap of each layer with the chosen interval (default 0–20 cm, decision D6); shallow first layers (e.g. 0–10 cm) are used as they are. SOTWIS missing values (−1) are excluded and component shares renormalised; `share_with_data` records the loss.

**K.** Williams/EPIC per component: K = f_csand · f_cl-si · f_orgc · f_hisand with f_csand = 0.2 + 0.3 exp[−0.0256 SAN (1 − SIL/100)], f_cl-si = (SIL/(CLA+SIL))^0.3, f_orgc = 1 − 0.25C/(C + exp(3.72 − 2.95C)), f_hisand = 1 − 0.7SN1/(SN1 + exp(−5.51 + 22.9SN1)), SN1 = 1 − SAN/100, C = organic carbon %. SI K = 0.1317 × US K. Because the equation is non-linear, K is computed per component and then weighted — not computed from weighted texture. Alternative: Renard et al. (1997) K from the geometric-mean particle diameter. CFRG = exp(−0.053 × rock %) is reported separately.

**Proxies.** USDA texture from the weighted fractions. Hydrologic soil group proxy: A (sand, loamy sand), B (sandy loam, loam, silt loam, silt), C (sandy clay loam), D (clay loams, clays); poorly/very poorly drained → D, imperfectly drained moves A/B to C. It is labelled a proxy wherever it appears.

**Aggregation.** The soil polygons are rasterised on the DEM grid by cell centre (even-odd rule — the default of `gdal.RasterizeLayer`, reproduced exactly on 1.44 million test cells), and each catchment value is the cell-weighted mean over units with data. `soil_coverage_pct` gives the covered share; the values describe that share only.

**Kenya check.** All 397 SOTWIS units give a K (median 0.029, range 0.009–0.051). Over the country's cells the median is 16–20 % above the ESDAC global K rasters with correlation 0.24–0.40; on six steep Rift-valley catchments, 0.026–0.038 against 0.027–0.028 (ESDAC Wischmeier-based). The two products are independent estimates at very different resolutions (1:1 M soil map vs 1 km model); report which one was used.

### 4.14 Erosion indices, RUSLE and erosion along a road

**Indices.** On the raw DEM (filling would zero the slope in depressions): Horn slope β; A_s = (upslope cells + 1)·cell area / cell size (D8 flow width = cell size); SPI = A_s·tan β; TWI = ln(A_s/tan β). tan β is floored at 0.001 inside the indices only; the slope raster is not floored. LS: Moore & Burch (1986) (A_s/22.13)^m·(sin β/0.0896)^1.3 with m = 0.4, or m by percent-slope class (0–1 % 0.2 … > 30 % 0.6, after Renard 1997 and McCool 1987, as in the A14 study); or Desmet & Govers (1996) per cell with McCool m and S. The A14 study used flow accumulation × cell size for A_s; QEHT adds the cell itself so headwater cells are not zero.

**RUSLE.** A = R·K·LS·C·P (t/ha/yr; R in MJ·mm/(ha·h·yr), K in SI). Every factor is a `Factor` with a source and a proxy flag. C from WorldCover is the mean of the A14 study's class ranges (tree 0.001–0.05, shrub 0.01–0.05, grass 0.01–0.15, crop 0.1–0.4, built-up 0.05–0.2, bare 0.4–0.6; water and wetlands from Panagos 2015 / Linard 2014) and is labelled a land-cover proxy. If R, K or C is missing no soil loss is computed and the classes use LS percentiles ("terrain potential").

**Classes.** A scheme is breaks, names and a severity score on a common 1–5 scale, so schemes with four and five classes combine. Shipped: SPI percentiles 50/75/90/97 (relative to the extent, D7) or fixed ln(SPI) 0/5/10 (A14); RUSLE 5/12/25/50 (D7) or 5/10/20/40 t/ha/yr (A14, after Watene et al. 2021); sediment volume 1,000/5,000/15,000 m³/yr. Combined = D7 matrix on the two scores (default: the higher). Class rasters are uint8, NoData 0, with a colour table, an attribute table (value, class, score) and a `.qml`.

**Catchment and crossing blocks.** Gross soil loss = mean A × area; SDR = 0.565·A_km²^−0.125 (FAO, capped at 1); sediment volume = yield × 1000 / bulk density (soil block, else the tool value). At the crossing: ln(SPI) max in the 3×3 outlet window and p50/p90 on the approach channel (up to 10 D8 steps), and composite = 0.4·SPI + 0.3·RUSLE + 0.3·sediment score with levels ≤ 2.0 Low, ≤ 3.0 Moderate, ≤ 4.0 High, else Severe (Akello & Omosa 2025). The test suite reproduces the paper's Table 14 composite scores and levels for all 14 crossings.

**Corridor.** Stations every Δs along the chainage; samples on perpendicular offsets to the half-width each side (nearest cell); per side max and mean of ln(SPI), LS, A and TWI; worst combined score in the buffer; reaches are runs of equal worst score with boundaries half-way between stations, so reach lengths sum to the sampled length.

### 4.15 Alignment ground profile (v0.15)

Stations every Δs (default 10 m) from the start chainage, plus the end. z_dem and z_fill by bilinear interpolation between cell centres (NaN next to NoData); pond depth = max(z_fill − z_dem, 0). Contributing area = largest (accumulation + 1)·cell area within one cell of the station. Stream crossings are the exact intersections of D8 stream links (cell centre to receiver centre, at the stream threshold) with the centreline — the method of 4.11 without a minimum area; each flags the nearest station, which takes the crossing's area if larger and its Strahler order. Longitudinal slope = centred difference of z_dem, one-sided at the ends and within each part of a multi-part alignment. Tested exact on a plane (z, ±2 % with direction, 1.9 % on an oblique line).

### 4.16 Flow-path segments and basin shape (v0.15)

**Segments.** Walking the LFP from the divide, the channel head is the first cell with accumulation ≥ the stream threshold. lfp_overland_m = distance to that cell, lfp_channel_m = the rest (they sum to lfp_length_m); slopes are end-point drops on the raw DEM over each length; lfp_channel_slope_1085 applies 4.10 to the channel part. The overland part is reported as sheet flow up to a cap (default 100 m, TR-55 practice) plus shallow concentrated flow. No channel cell → whole path overland, lfp_no_channel = 1.

**Shape.** perimeter_km from the catchment mask smoothed with a 3×3 mean and contoured at 0.5 by marching squares with interpolation (within ~1 % of a disc's circumference and ~1.5 % of a square's; the cell-edge staircase overstates a round outline by 27 %). form_factor = A/L² (Horton), elongation_ratio = (2/L)√(A/π) (Schumm), circularity_ratio = 4πA/P² (Miller), L = LFP length. drainage_density = channel length along D8 links inside the catchment / A; stream_frequency = links / A with links = sources + confluence cells (Shreve links); max_strahler from the stream-order raster or computed at the threshold.

### 4.17 Soils from any source (v0.15)

Every source is turned into a `SoilGrid`: per-cell sand, silt, clay, organic carbon, bulk density, coarse fragments, K (SI), Dg-K, drainage class and hydrologic soil group on the DEM grid. Unit sources (SOTWIS / SOTER database, polygon attributes or user-named fields, polygons + CSV, unit raster + CSV) paint each unit's description — K still computed per component, then weighted — onto its cells; the catchment block is then identical to 4.13 (tested field by field). Texture rasters give per-cell values with K, texture and HSG proxy computed per cell (vectorised versions of the 4.13 functions, tested equal on 5,000 random textures). SoilGrids 2.0 values are divided by ISRIC's factors (sand/silt/clay g/kg ÷ 10 → %, soc dg/kg ÷ 100 → %, bdod cg/cm³ ÷ 100 → g/cm³, cfvo ‰ ÷ 10 → vol %); depth layers are weighted by their overlap with the chosen interval, per cell over the layers present. Overrides: a hydrologic soil group raster (HYSOGs250m codes 1–4 = A–D, 11–14 = dual groups counted as D, 255 NoData) or polygon field replaces the proxy where known; a K raster or field replaces the computed K where finite. The catchment block adds soil_hsg (largest share), hsg_pct_a…d, hsg_pct_dual, soil_hsg_source and usle_k_source.

### 4.18 Curve number and Rational C (v0.15)

Per cell, land cover (WorldCover, nearest-neighbour onto the DEM grid) × HSG (4.17) → CN (AMC II) from a lookup: by default USDA TR-55 (1986) Table 2-2 with each WorldCover class matched to one cover (tree → woods, shrub → brush, grass → pasture/grassland/range, crop → row crops straight row, built-up → commercial 85 % impervious, bare → fallow bare soil, water → 100; wetland, mangroves, snow and moss have no default), for the chosen condition (fair/good/poor; row crops have no 'fair' row and use 'good'). The match is a proxy and is flagged in `runoff_json`; a user CSV (class, A, B, C, D) replaces it. cn_ii is the cell mean over cells with a CN (cn_coverage_pct reports the rest); cn_export converts it to AMC I or III with CN_I = 4.2CN/(10 − 0.058CN), CN_III = 23CN/(10 + 0.13CN) (Chow, Maidment & Mays 1988). Rational C works the same way from a user lookup only. Land-cover shares are given whenever a land-cover raster is.

### 4.19 STI, deposition at crossings, channel slopes (v0.15)

**STI** (sediment transport capacity index, Moore & Wilson 1992) = (A_s/22.13)^0.6 (sin β/0.0896)^1.3 — the 4.14 LS form with m = 0.6 (tested identical) — on overland cells only, channels NoData and A_s capped at 100 m. A relative indicator; never used in soil loss, classes or the composite.

**Deposition indicator.** The main stem upstream of a crossing is followed by taking, at each cell, the donor with the largest accumulation. On channel cells, median SPI over 0–100 m (near) and 100–500 m (far) upstream; ero_dep_ratio = near / far; below 0.7 deposition-prone, above 1.3 scour-prone (editable). Both reach slopes (drop / length) are exported, because over a short reach the ratio is mostly the slope break — unless the area grows fast along it, which the ratio then includes. Null with a note when the channel upstream is shorter than 90 % of the far distance.

**Channel slopes.** ch_slope_us along the main stem upstream and ch_slope_ds along the D8 receivers downstream, each over 200 m (or less at the divide or grid edge, reported in ch_len_us_m / ch_len_ds_m), from raw-DEM end-point elevations.

### 4.20 Drainage coverage check (v0.16)

From the alignment profile (4.15). **Missing crossings:** stations with stream = 1 and area ≥ the minimum (default the stream threshold) merged within 30 m into one point at the largest area; uncovered when no existing crossing lies within 50 m of chainage. The package builds each uncovered point as a proposed crossing (snapped within 2 cells to the stream) with the full pipeline (4.8–4.19), IDs P001… gapless along the chainage, `nearest_uid` / `nearest_m` recorded. **Sags:** local minima of the profile smoothed by a 30 m moving average whose prominence (the lower of the highest smoothed ground on each side before the profile drops below the minimum again) is ≥ 0.3 m, with no stream station or crossing within 50 m. On the filled DEM the centreline cells are raised 1,000 m (diagonal gaps closed; wall removed within one cell of existing and stream crossings) and D8 is computed without refilling; the sag area on each side is the largest accumulation among pit cells (no receiver) within two cells of a point 1.5 cells off the centreline. Sags of ≥ 0.05 km² become proposed crossings delineated on that walled routing. **Flat stretches:** runs of stations with |longitudinal slope| and the steeper cross-fall over 100 m both below 0.5 %, boundaries half-way between stations, ≥ 300 m long.

### 4.21 Rainfall zone, mean annual rainfall and R estimate (v0.17)

**Zones.** The zone polygons are rasterised on the DEM grid by the cell-centre rule (as the soil map; where polygons overlap the later one wins), polygons with the same name forming one zone. Per catchment, each zone's share is its cell count over the cells covered by any zone; `rain_zone` is the largest, ties going to the name that sorts first (case-insensitive), so the result does not depend on the order of the layer. `rain_zone_coverage_pct` is the covered share of the catchment.

**Mean annual rainfall.** Resampled bilinearly to the DEM grid; `map_mm` is the mean over cells with data. A raster whose 99.5th percentile is below 20 mm/yr (mm/day or a monthly mean), above 13,000 mm/yr (tenths of mm or a NoData problem) or with negative values is refused; below 100 mm/yr gives a warning.

**R estimate.** Per cell from the rainfall P (mm/yr), then the catchment mean (R is not linear in P, so this differs from R of the mean P):

- Renard & Freimund (1994): R = 0.0483 P^1.61 for P ≤ 850 mm, R = 587.8 − 1.219 P + 0.004105 P² above (the branches meet within 0.2 % at 850 mm);
- Lo et al. (1985): R = 38.46 + 3.48 P;

in MJ mm ha⁻¹ h⁻¹ yr⁻¹. Neither was derived in East Africa. The value is labelled an estimate in `rusle_r_method`, `rainfall_json` and, when the erosion tool uses it, the factor basis "estimate from rainfall".

### 4.22 One-click pipeline and run report (v0.17)

The pipeline runs the standalone tools as child algorithms with the parameters below, so it introduces no numerical code of its own: Fill (min slope 0) → D8 (chosen flat method) → accumulation → streams (threshold = km² ÷ cell area, rounded) → crossing candidates (with a road and no crossing layer) → Burn crossings and the four routing steps again on the burned DEM (when asked) → stream polylines → Erosion indices (channel threshold = stream threshold; R from the R raster, else the value, else the rainfall estimate) → Build design hydrology package (snap 5 cells to the stream threshold for pour points; candidates keep their outlets; road, coverage check, soils, CN, rainfall, erosion folder and burn log passed through). After the package, the flat-method sensitivity (4.4, v0.13.1) is computed at every crossing outlet on the filled DEM. The package layers are copied to `layers/` with in-process GDAL `VectorTranslate`, the tables to `tables/`, and the characteristics table and report are built from the package.

The **characteristics table** joins crossings, catchments and flow paths on `outlet_uid` (core columns always; groups for flow-path segments, channel slopes, flat check, soils, CN, rainfall, shape and erosion when at least one crossing has a value) and orders rows by chainage when present. The **run report** reads `qeht_run_metadata` and the layers only. Both are pure Python (`core/report`).

### 4.23 Floodplain width indicator (v0.18, A4)

At crossings whose contributing area (acc_at_outlet_km2) is at least the limit (10 km²) and that have a chainage, from the alignment profile (4.15; raw DEM, bilinear, stations every 10 m): z_bed is the lowest z_dem_m within 50 m of the crossing's chainage, at station i₀. For each Δz in 0.5, 1, 2 m the run of stations around i₀ with z ≤ z_bed + Δz is extended both ways; each end is interpolated linearly between the last station inside and the first outside, so a linear bank is measured exactly (a trapezoid with a 40 m floor and 1:20 banks gives 60 / 80 / 120 m). A NoData station or the end of the profile stops the run, and fp_note says the width is a lower bound.

HAND (Rennó et al. 2008; Nobre et al. 2011) is z minus the elevation of the first stream cell (accumulation ≥ the stream threshold) on the cell's D8 path, computed for the whole grid by pointer jumping (each pass doubles the path length followed, so a path of n cells takes log₂ n passes); cells whose path ends without reaching a stream are NoData. HAND is read at the station's cell, and the same contiguous run is taken with HAND ≤ Δz around the station of least HAND within 50 m. HAND follows the drainage, so on a valley that slopes along the road the two methods differ; report both.

### 4.24 Raster quicklooks (v0.18, A6)

Each raster is downsampled by the integer factor k = ⌈long side / limit⌉: classified rasters take the cell nearest each block centre (class values preserved), continuous ones the mean of the block's finite cells. Classified rasters are coloured from their GDAL colour table (the same colours as their .qml) with their category names in the legend; flow accumulation uses a blue ramp on log₁₀(1 + cells) from 0 to the maximum; STI a yellow–red ramp over the 2nd–98th percentile; the relief is a hypsometric tint (2nd–98th percentile of elevation) multiplied by 0.35 + 0.65 × hillshade (Horn, azimuth 315°, altitude 45°). The PNG is RGBA with NoData transparent. The world file holds the pixel size and the centre of the upper-left pixel of the downsampled grid.

---

## 5. Data handling and interoperability

- **Rasters** are read and written through the GDAL Python bindings directly. Rotated or skewed grids are rejected, because D8 assumes a north-up grid and silently accepting a rotated one produces plausible-looking but wrong directions.
- **Vectorisation** uses `gdal.Polygonize()` in-process — the call the standard workflow reaches by shelling out to `gdal_polygonize.exe`. No external executable is launched.
- **Flow-direction interchange** works in both directions via the standard D8 encoding, so a reference hydrology toolset's grid can be ingested and a QEHT grid exported.
- **CRS** is carried with each raster and pour points are reprojected into the DEM CRS as needed.

### 5.1 Design hydrology package (schema `qeht-heas-1`)

"Build design hydrology package" (named "Build HEAS exchange package" up to 0.14; algorithm id `buildheasexchange` unchanged) writes one GeoPackage per run with `crossings` (snapped outlets), `catchments` and `flowpaths`, plus two attribute tables: `qeht_run_metadata` (key/value provenance: QEHT version, run time, CRS, cell size, DEM path and SHA-256, user-declared DEM source, conditioning and flat method, stream threshold, snap radius and strategy, catchment mode, ID scheme, full parameter set) and `qeht_field_dictionary` (field, type, unit, meaning, method, HEAS target, plain-language downstream use and HEAS field name for every field, including the optional layers — the file documents itself). HEAS, a separate design tool, imports the package directly; it is not required.

- **Linking.** Every feature carries `outlet_uid`, identical across the three layers for one crossing, plus `crossing_id`/`catchment_id`/`flowpath_id` (default = `outlet_uid`), `link_method` (`pour_point`), `link_confidence` (1.0) and `link_note`. IDs come from a chosen pour-point attribute or are sequential with a prefix, numbered downstream-first (largest contributing area = 001) or in layer order. Empty or duplicate IDs, and two pour points snapping to the same cell, stop the run with a list. `outlet_id` (feature id) is kept for backward compatibility and is not stable.
- **CRS.** A projected, metric CRS is required; a geographic DEM is refused before any computation.
- **Values.** Missing values are NULL, never 0. `acc_at_outlet_km2` is read from the accumulation raster ((accumulation + 1) × cell area) and equals `area_km2` for full catchments — a built-in QA check.
- **Renumber and relink (5.2).** After crossings are deleted, moved or added in a package, "Renumber and relink exchange package" writes a new package: unmoved crossings keep their outlet cell, moved/added ones are snapped; crossings are ordered by chainage (recomputed from the road for moved/added points) or downstream-first; IDs are re-issued gaplessly with the package's prefix; catchments and flow paths are recomputed; `renumber_log` lists old → new (`unchanged`, `renumbered`, `new`, `deleted`, and `moved_m`). IDs are only re-issued when this tool is run.
- **v0.17 additions (additive).** Catchments: `rain_zone`, `rain_zone_pct`, `rain_zones_json`, `rain_zone_coverage_pct`, `map_mm`, `map_coverage_pct`, `map_dataset`, `rusle_r`, `rusle_r_method`. Metadata: `rainfall_json`. The golden fixture was regenerated after `golden_diff` showed no existing value changed.
- **v0.18 additions (additive).** Crossings: the eleven `fp_*` fields (4.23). Table `rasters` (4.24). Metadata: `fp_params_json`, `quicklook_params_json`.
- **Evolution.** New fields may be added under `qeht-heas-1`; renaming or removing a field requires `qeht-heas-2`. The validator (also used by the tool after writing) rejects unknown major versions and checks that every catchment and flow path refers to an existing crossing.
- **Erosion block (v0.14).** With an erosion folder, catchments carry `ero_*` soil loss, LS, ln(SPI), K/C/P, class shares, SDR and sediment fields and crossings carry `ero_*` ln(SPI), local terrain, class and impact fields; `erosion_json` in the metadata records the mode, factor sources and proxy flags, schemes, weights and SDR model. Additive under `qeht-heas-1`.
- **Burn log (v0.14).** Give the breach log from "Burn crossings through embankments" and the package stores it as layer `burn_log`; `conditioning` appends a summary (breaches cut, volume) to the user-declared text and `conditioning_burn` holds it on its own.
- **v0.15 additions (all additive).** Crossings: ch_slope_us/ds and lengths; with an erosion folder ero_sti_local and the deposition indicator. Catchments: shape and network indices, HSG shares and sources, curve number / Rational C / land-cover block. Flow paths: overland / channel segments. Layer `alignment_profile` with a road alignment. Metadata: alignment_source, alignment_step_m, sheet_cap_m, channel_slope_m, soil_hsg_source, soil_k_source, runoff_json. The golden fixture was regenerated; no existing value changed.
- **Road layers.** When built from a candidate layer, the package also holds `crossing_candidates` (every candidate, the audit trail) and, when a road is given, `road_alignment`. The crossings keep each candidate's outlet cell (no snapping), carry `chainage_m`, and are numbered along the chainage. Metadata records `crossing_source` and `chainage_start_m`.
- **Implementation.** The file is written with the Python standard library (`sqlite3`, `struct`), so it is identical on every platform and testable without GDAL; the test suite validates it with GDAL's GeoPackage validator. `tests/fixtures/golden_exchange.gpkg` (synthetic DEM, three crossings: two nested on one valley, one on a tributary) is the shared contract fixture: QEHT asserts it reproduces the same attributes, HEAS asserts it imports with no field mapping.

---

## 6. Validation record

QEHT is validated at three levels.

**Level 1 — synthetic analytic DEMs.** 47 checks in `tests/test_core.py`, runnable on bare Python + NumPy, covering encoding round-trips, distance weighting, fill spill levels, accumulation on analytic surfaces, catchment areas, longest-flow-path geometry, snapping, Strahler rules, and the Barnes saddle convergence test. A further 74 checks in `tests/test_interop.py` cover the 10–85 conventions on an analytic profile, the slope domains, `outlet_uid` rules, the polygon tracer, the GeoPackage writer, the CRS rule and the exchange package against the golden fixture. All pass. (Earlier documentation quoted 48 core checks; the suite has 47.)

**Road crossings.** 43 checks in `tests/test_crossings.py`: linear referencing (chainage, signed offset, multi-part, reverse); a square crossing (exact chainage, 90°, side); a 45° road (45° exactly); a stream weaving across the centreline for ~600 m (11 raw intersections → one cluster, the exit recommended, a 748 m parallel reach); two tributaries 100 m apart (merge 50 m → two clusters, 150 m → one); candidate selection rules; candidates → exchange (numbered along chainage, outlets not snapped); an embanked valley (pond removed, breach confined to the culvert line, flow through the culvert); renumber-log rules.

**Soils.** 31 checks in `tests/test_soils.py`: hand-computed Williams factors, the 0.0256 coefficient, K monotonic in organic carbon and texture, Dg-K bounds, CFRG, 12 USDA classes, HSG rules, the SOTWIS loader on a synthetic database with the SOTWIS schema (component and depth weighting, shallow layers, −1 values, TTR), exact cell weights and partial coverage on a three-unit map, no-overlap case, and the soil block on exchange catchments.

**QGIS level.** `tests/qgis_smoke.py` runs every Processing tool through `processing.run()` inside QGIS on the example DEM and checks outputs, `outlet_uid` consistency and the exchange package. Passes on QGIS 3.34.4 (headless).

**Level 3 — reference hydrology toolset production output (Site A, steep).** 718×775 cells, EPSG:21037, 16 road-crossing pour points with reference catchments and longest flow paths. Catchment area median ratio 1.009 (15/16 within ±30%), longest-flow-path length median ratio 0.999 (15/16 within ±20%), best individual match 0.05%, zero flow-direction cycles. Cell-by-cell on an aligned grid, flow direction agrees with the reference toolset on 98.5%.

**Flat-terrain benchmark (Site B, coastal).** 34%-flat coastal DEM. Same-grid flow-direction agreement 88%; the residual is concentrated in flats and is a genuine algorithmic difference (Section 4.4). Accumulation validated in isolation to 99.31% exact against the reference platform (Section 4.5).

The two datasets bound the expected accuracy: near-exact on moderate-to-steep terrain, and a documented, quantified flat-terrain difference concentrated entirely in flow-direction assignment over flats.

---

## 7. Known limitations

- **No tiling.** The whole grid is processed in memory. Measured peaks for the v0.12 core (bytes per cell, plus ~24 resident): fill 44, flow direction 119 (toward) / 165 (Barnes), accumulation, Strahler, catchment and longest flow path ~82. Every tool reports its estimate before running and warns when it approaches the free RAM; clip to the catchment plus a buffer when it does — results inside the clip are unchanged.
- **Flat terrain.** On very flat terrain QEHT's toward-lower option reproduced the commercial reference platform's output to a stream IoU of ~0.6 (Site B; to be repeated with the fixed Barnes default). This is inherent to the algorithmic difference in Section 4.4, not a defect. For flat coastal catchments where the last increment matters, a reference-platform flow-direction grid can be ingested and QEHT's accumulation, catchment, and characteristics tools run on top of it.
- **Flat method depends on the DEM.** On whole-metre DEMs (ALOS) half the cells of flat terrain can be flat and the method changes the network; on floating-point DEMs (FABDEM) flats are ~1 % of cells and it matters little.
- **Single-threaded.** No parallelism. Since v0.12 the core is vectorised NumPy (fill, direction, flats, accumulation, Strahler, delineation, longest flow path, snapping), 3.6–5.4× faster end to end on 1-megapixel clips.

---

## 8. Release readiness (GitHub and QGIS Plugin Store)

Completed as of 0.8.2:

- README version header corrected to 0.8.2.
- `LICENSE` file added (GPL-2.0-or-later), author attributed to Edmond Akello.
- `CITATION.cff` added, author attributed to Edmond Akello.
- `metadata.txt`: author, email, and real homepage/repository/tracker URLs
  (`https://github.com/EdmondAkello/qeht`) populated; `changelog` pointer
  added.
- `CHANGELOG.md` added, derived from the version history.
- Reproducible example datasets added: a real-world DEM (`examples/example_dem.tif`,
  ~2.7 MB) and a synthetic DEM with no real-world correspondence
  (`examples/synthetic_dem.tif`, ~0.4 MB), both with a documented end-to-end
  workflow (`examples/README.md`).
- Documentation pass replacing project-identifying site names with generic
  labels throughout, ahead of public distribution.

Still recommended before a 1.0 / plugin-store submission:

- Add a QGIS-level smoke test (provider loads, one algorithm runs against
  the bundled example DEM) alongside the core tests — the core test suite
  covers the numerical engine but does not exercise the QGIS Processing
  wiring itself.
- Decide whether to clear the `experimental` flag; it is honest to retain
  it until the QGIS-level smoke test above is in place, even though the
  example dataset is now bundled.
- Establish a DOI/archive (e.g. Zenodo) for the citable release.

The Processing-provider architecture QEHT uses is explicitly consistent with QGIS's recommended approach for analytical plugins, and the in-process, no-binary design aligns with the store's cross-platform and no-external-binary expectations.

---

## 9. Publication strategy

The material supports two distinct outputs, and they should not be collapsed into one.

**Software paper — QEHT itself.** A comprehensive description of architecture, algorithms, interoperability, testing, and validation. Natural venues: the Journal of Open Source Software (JOSS), which reviews the software, its documentation and tests, and expects a public open-development history of at least six months — so the GitHub repository should be established early; or SoftwareX, which pairs a short software paper with the open-source release.

**Research paper — the flat-resolution finding.** This is the stronger scholarly contribution, and it is a genuine result rather than a tool description: *the influence of DEM flat-resolution algorithm on drainage-network structure and engineering catchment characteristics*. The controlled comparison of toward-lower versus canonical Barnes across terrain regimes — quantified on flow-direction agreement, accumulation, stream IoU, catchment area, longest-flow-path length, and slope, and traceable through to time of concentration and design discharge — is a defensible, reproducible experiment. The finding that a theoretically canonical convergent algorithm reproduces a widely-used reference implementation *less* well than a simpler method, on a specific and common terrain regime, is the kind of result that is worth reporting. Natural venue: Environmental Modelling & Software, whose scope explicitly covers hydrological modelling, GIS, and quantitative comparison of alternative methods. Computers & Geosciences is possible but would require framing the contribution as more than a standard implementation.

The recommended sequencing is to establish the public GitHub repository now — as both the product repository and the experimental record, holding the benchmark datasets, scripts, and validation results — so that the open-development history accrues while the research paper's controlled benchmarks are assembled. A v0.9 focused on reproducibility and a formal validation framework, rather than new features, is the right next step; the v0.8 code already contains the substance both papers need.

---

## 10. Citation

See `CITATION.cff`. In brief: QEHT: QGIS Engineering Hydrology Toolkit, v0.14.0, GPL-2.0-or-later.

## 11. Licence

GNU General Public License v2 or later. See `LICENSE`.
