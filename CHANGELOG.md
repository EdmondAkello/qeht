# Changelog

All notable changes to QEHT are recorded here. Versions follow the
development history of the numerical core and Processing tools.

## [0.15.1] — unreleased (fixes from the Site C test run)
### Fixed
- **Vector outputs honour the file extension.** Choosing `.shp` (or `.geojson` …) wrote a GeoPackage under that name, which other software could not open. All vector writers now pick the OGR driver from the extension (GeoPackage otherwise) and set fields by index, so 10-character shapefile names work.
- **Candidate layers saved as shapefiles are read back.** Truncated field names (`recommende`, `outlet_row` …) are matched, so "Build design hydrology package" no longer stops with "No candidate is accepted or recommended".
- **Run metadata records layer files, not project layer ids** (`parameters_json`), and empty inputs as ''.

### Added
- **RUSLE factor range checks.** A land-cover class raster given as the C (or P) raster is refused with a hint to use the WorldCover input. That was the mistake behind C ≈ 30 and soil loss two orders too high in the test run. Also refused: C or P outside 0–1, K in US units (with the 0.1317 factor), negative values, R out of range. R below 200 MJ·mm/(ha·h·yr) gives a warning. Single values are labelled "user value" in `erosion_run.json`.
- **Rasters record how they were made.** Fill writes `QEHT_CONDITIONING` and D8 flow direction writes `QEHT_FLAT_METHOD` as GeoTIFF metadata. The package fills `conditioning` and `flat_method` from them when the fields are left blank. `dem_source` defaults to the DEM file name.
- **Automatic merge of near-duplicate crossings.** Consecutive candidates closer than two cell diagonals whose areas agree within 5 % (one stream crossing a winding centreline twice, 1.8 m apart on the corridor) share a cluster, so only one is recommended.
- **Catchment characteristics takes a crossing-candidate layer** (recommended / accepted outlets, no snapping), like the package.
- The flat-method check reads the largest area within one cell of each outlet under each method, so a drainage line routed one cell apart is not reported as a change.
- Tests: crossings 45 (+2), erosion 50 (+6), smoke 49 (+5).

### Notes
- On the Site C DEM (FABDEM reprojected with nearest neighbour: 1.5 % repeated rows, 14 % tied neighbours) the flat-method check still flags every crossing; toward-lower routing moves the main river between crossings. This is the 0.13.1 Site C finding again. Re-export the DEM with bilinear resampling, or use native-resolution FABDEM, before relying on flat-sensitive areas.

## [0.15.0] — unreleased (alignment profile, flow-path segments, soils anywhere, curve numbers)
### Changed
- **"Build HEAS exchange package" is now "Build design hydrology package".** The algorithm id `buildheasexchange` and the schema `qeht-heas-1` are unchanged, so saved models and HEAS imports keep working. The field dictionary gains `downstream_use` (plain language) and `heas_field` columns; `heas_target` stays. README: HEAS moves to one optional section.
- "Soil parameters for catchments" is shown as "Soil and runoff parameters for catchments" (id unchanged); its soil polygons are now optional (other sources below).

### Added
- **Alignment ground profile** (Road drainage): stations every 10 m with raw and filled DEM (bilinear), ponding depth, contributing area, exact D8 stream crossings with Strahler order and longitudinal slope; CSV and chart. The package writes the same profile as layer `alignment_profile` when a road is given (optional `FILLED`, `PROFILE_STEP`).
- **Flow-path segments (F2):** `lfp_overland_m` / `_slope`, `lfp_channel_m` / `_slope` / `_slope_1085`, `lfp_sheet_m` (cap, default 100 m), `lfp_shallow_m`, `lfp_threshold_km2`, `lfp_no_channel`, split at the stream threshold.
- **Basin shape and network (F6):** `perimeter_km` (smoothed outline), `form_factor`, `elongation_ratio`, `circularity_ratio`, `drainage_density`, `stream_frequency`, `max_strahler`.
- **Soils from any source (G1–G4):** one per-cell soil model (`core/soils/sources.py`) for SOTWIS / SOTER, polygon fields chosen by name (OC % or g/kg), polygons + CSV, a soil unit raster + CSV (e.g. HWSD v2), texture rasters or a SoilGrids 2.0 folder (unit conversion, depth weighting), plus HSG (HYSOGs250m codes, dual groups as D) and K overrides. New catchment fields `soil_hsg`, `hsg_pct_a`…`d`, `hsg_pct_dual`, `soil_hsg_source`, `usle_k_source`. Unit sources reproduce the 0.14 block exactly.
- **Curve number and Rational C (F1):** land cover (WorldCover) × HSG per cell → CN from TR-55 Table 2-2 (condition fair/good/poor; a proxy match, flagged), `cn_ii`, `cn_export` at AMC I/II/III, coverage, land-cover shares `lc_pct_*`; user CN lookup CSV; Rational C from a user lookup only (none ships). Metadata `runoff_json`.
- **Sediment transport (STI advisory R1, R2):** `sti_overland.tif` from the erosion tool (Moore & Wilson 1992, overland cells, A_s capped at 100 m; not used in soil loss, classes or composite); per crossing `ero_sti_local` and the deposition indicator `ero_spi_app_near` / `_far`, `ero_dep_ratio`, `ero_dep_flag` (0.7 / 1.3), reach slopes, `ero_dep_note`.
- **Approach / exit channel slopes (A5):** `ch_slope_us` (main stem) and `ch_slope_ds` (D8 path) over 200 m on every crossing, with the lengths used.
- README "Data sources by region" (global DEM, soil, HSG, land cover, R and K sources); help text no longer assumes Kenyan data.
- Tests: test_alignment 15, test_morphometry 22, test_soils_any 18, test_runoff 17, test_channel 16 (new); smoke test 44 (+7). Golden fixture regenerated — new fields and metadata keys only, no existing value changed.

### Notes
- The TR-55 → WorldCover match and the WorldCover → C lookup are proxies for review. Rational C needs your table (e.g. RDM Part II rows).
- Not yet: drainage coverage check, sag points and flat stretches (A2, A3), rainfall zones (F3), floodplain width (A4), cross-sections (F5), one-click pipeline and run report (F7, F8), raster quicklooks (A6), corridor STI (R3).

## [0.14.0] — unreleased (WP-F erosion; slope review)
### Changed
- Documentation: the relief ratio is described as a basin-steepness index (QEHT divides by LFP length, a variant of Schumm's basin-length ratio) and no longer as the Kirpich slope. Kirpich, Bransby-Williams and TRRL use the flow-path slope (`lfp_slope` / `lfp_slope_1085`). Statistics docstring, README, tool help, DOCUMENTATION.md and field dictionary aligned. No values change.
- Build HEAS exchange package and Catchment characteristics print one run line: QEHT version · flat method · 10–85 reference. README gains a "Re-running a design dataset" box.

### Added
- `lfp_z_head_m`, `lfp_z_outlet_m`, `lfp_nonmonotonic` on flow paths: `lfp_drop_m` is highest − lowest raw-DEM elevation on the path; the flag marks a spike or pit (> 0.5 m) that makes it exceed headwater − outlet (additive; schema still `qeht-heas-1`).
- **Erosion indices and RUSLE soil loss** (Soils and erosion group): slope, A_s, SPI, ln(SPI), TWI and LS (Moore & Burch with m = 0.4 or m by slope class; Desmet & Govers); RUSLE A = R·K·LS·C·P with each factor from a raster, a value or a dataset (K from SOTWIS; C from ESA WorldCover 2021 through an editable lookup, flagged as a proxy), resampled in-process to the DEM grid; LS-only "terrain potential" when R, K or C is missing. Severity classes (D7 and the A14 2025 schemes, or custom breaks; 5 × 5 combination matrix) as uint8 rasters with colour table, attribute table and `.qml`; `class_extents.csv`; `erosion_run.json`.
- **Sample erosion along alignment**: stations every 10 m with left/right statistics, reaches by worst class, optional chainage chart.
- **Erosion block in the HEAS exchange package** (optional erosion folder): per-catchment soil loss, LS, ln(SPI), K/C/P, class shares, gross loss, SDR (0.565·A^−0.125) and sediment volume; per-crossing ln(SPI), local terrain and the hydrodynamic impact score 0.4·SPI + 0.3·RUSLE + 0.3·sediment (Akello & Omosa 2025). Additive (`ero_*` fields, `erosion_json` metadata); schema still `qeht-heas-1`.
- `core/erosion/` (terrain, rusle, classes, summary, corridor, io); `core.raster.warp_to_grid`.
- **Breach log in the exchange package**: Build HEAS exchange package takes the log from Burn crossings; stored as layer `burn_log`, summarised in `conditioning` and the new metadata key `conditioning_burn`.
- **Soil table loader**: soil polygons with only a unit code plus a CSV of sand, silt, clay, oc [bulk, cfrag, drain], joined on the unit field (`core.soils.sotwis.load_csv_units`), on every tool with soil inputs.
- Tests: test_interop 77 (+3), test_soils 33 (+2), test_erosion 44 (analytic plane and valley; the A14 paper's Tables 4–9 and the 14 composite scores of Table 14), smoke test 37.

### Notes
- The WorldCover → C lookup is a draft for review; edit `core/erosion/rusle.WORLDCOVER_C` or supply a C raster.

## [0.13.1] — unreleased (real-road review follow-ups)
Found on Site C, a road project (old-tool outputs vs 0.13; report kept outside the repo).
### Added
- **DEM QA: nearest-neighbour resampling.** Fill and D8 flow direction warn when rows or columns exactly repeat their neighbour (`core.raster.resampling_stats` / `audit_resampling`, threshold 0.2 %). On the project DEM (FABDEM reprojected with nearest neighbour) 1.9 % of rows and 1.2 % of columns repeated and 14 % of cells were tied; one tied cell on a river moved ~176 km² between two culverts depending on the flat method.
- **Flat-method check in Catchment characteristics** (on by default, tolerance 10 %): `area_barnes_km2`, `area_toward_km2`, `flat_sensitivity_pct`, `flat_sensitive`, plus a warning listing flagged outlets (`core/flow/sensitivity.py`). On the project it flagged exactly the four unstable crossings; none on a bilinear resample.
- Tests: test_flats 36 checks (+5); smoke test 30 checks (+2).

### Changed
- **ID prefix applies to sequential IDs only.** With an ID attribute the default prefix "X" is no longer prepended (it turned NS1 into XNS1 in Catchment characteristics and Build HEAS exchange package). The core `assign_uids(scheme="attribute", prefix=...)` keeps an explicit prefix for scripted use.

### Notes
- An Earth Engine FABDEM export at 0.000269° (finer than the native 1″) repeats about one row and column in 33; every FABDEM benchmark clip in WP-G carries this (see BENCHMARK.md). ALOS clips were native.

## [0.13.0] — unreleased (WP-G benchmark, QGIS 4)
### Changed
- **Barnes 2014 is the default flat method** (D8 flow direction and `core.flow.direction.d8_direction`). Decided on the WP-G benchmark: 80 stratified random 15 km areas across Kenya (flat, rolling, hilly, mountainous), each on ALOS AW3D30 and FABDEM, against TauDEM, RichDEM, MAS 1.2.1 and DDM HydroLogic 2.3. Toward-lower and Barnes are not interchangeable; Barnes agrees best with the independent TauDEM (mean stream F1 0.969 vs 0.948; 0.870 vs 0.733 on whole-metre flats). Design, tables and decision in `BENCHMARK.md`.
- **Toward lower terrain** remains as an option (enum index unchanged, so saved models keep their method).
- `cells_still_unrouted` no longer counts boundary outlets; new stat `boundary_outlets`.

### Removed
- **Hybrid flat method** (and the FLAT_WEIGHT / SMALL_FLATS parameters): it reproduced Barnes (stream F1 ≥ 0.99 in every stratum). `core/flow/flats.resolve_flats(w=…)` keeps w for research.

### Fixed
- **Flats touching the grid edge or NoData stayed unrouted** (all methods). Boundary cells without descent are now outlets and such flats drain to them, as in TauDEM and RichDEM — relevant for clips and coastlines.
- **Barnes departed from the published algorithm at the flat's low edge**: a cell there could run along the rim before leaving; the neighbour scan order also differed. 0.12 Barnes agreed with RichDEM on 78 % of flat cells (one flat ALOS area); 0.13 matches RichDEM cell for cell on all 160 benchmark runs.

### Added
- `core/flow/_reference.reference_barnes_richdem`: per-cell port of RichDEM's Barnes code, the new test oracle (exact match on 80 random surfaces with NoData). Tests: edge/NoData drainage, Barnes default, hybrid refused (test_flats 31 checks; smoke 28).
- `BENCHMARK.md`.

### Compatibility
- **QGIS 4.2.2** (Qt 6): all unit suites and the 28-check smoke test pass headless on the official Debian trixie packages; the plugin loads, registers 13 algorithms and unloads. Also tested on QGIS 3.34.4. `qgisMaximumVersion=4.99`.

## [0.12.0] — unreleased (WP-E flats and performance)
### Fixed
- **Barnes flat resolution could create flow loops.** Its outlet rule for local-minimum flat cells accepted an equal-elevation neighbour assigned earlier in the same pass, which could point two cells at each other. The loops are small but sit on main channels and cut them off from their upstream area (on a 30 m flat coastal clip the largest accumulation was 40,117 cells; 794,103 with the fix). The exit now goes only to a lower neighbour or one routed before flat resolution. Toward-lower (the default) was not affected. The earlier Barnes-vs-reference comparison on coastal flats (stream IoU 0.19) used the faulty code and should be repeated.

### Added
- **Hybrid flat resolution**: `flat_mask = w·toward + (flat_height − away)`, w > 1 (w = 2 is Barnes), with an optional size switch (toward gradient only below N cells). New options on D8 flow direction. Loop-free for every w > 1 (tested); on broad real flats w changes little — tie rules among equally short paths matter more.
- **Memory estimate** (`core/memory.py`) from measured per-step peaks, reported by every heavy tool before it runs, with a warning when it nears the free RAM.
- `core/flow/_reference.py`: the v0.8.3 per-cell algorithms, kept verbatim as test oracles. `tests/test_flats.py` (24 checks).

### Changed
- **Vectorised core**: fill (frontier relaxation of the priority-flood solution, with and without min_slope), D8 flow direction (streamed over the 8 directions: ~176 → ~59 bytes per cell), toward-lower and Barnes flats, accumulation, Strahler, catchment delineation, longest flow path (FIFO order reproduced so equal-length ties choose the same branch) and nearest-stream snapping. Outputs identical to v0.8.3 (except the Barnes fix above; fractional-weight accumulation to 1e-12). Full chain on ~1-megapixel 30 m clips: 10.5 s → 2.9 s (steep), 25.9 s → 4.8 s (flat coastal).
- Documentation corrected: the D8 tie rule is the fixed priority S, W, N, E, SE, SW, NW, NE (not lowest index).

## [0.11.0] — unreleased (WP-C soils)
### Added
- **Soil parameters for catchments** (new "Soils and erosion" group) and an optional soil input on **Build HEAS exchange package**: a soil block on every catchment — topsoil sand, silt, clay, organic carbon, coarse fragments, bulk density, USDA texture, FAO drainage class, USLE K (Williams/EPIC, SI) and the Renard Dg-based K, CFRG, a texture/drainage hydrologic-group proxy, coverage %, dominant unit, TTR class, dataset and depth. Optional USLE K and soil-unit rasters.
- **SOTWIS preset** (standard-library `sqlite3`, never the .mdb): full unit composition, depth-weighted over a chosen interval (default 0–20 cm), K per component then weighted. Also reads polygons with SOTWIS fields (dominant soil) or plain sand/silt/clay/oc fields.
- `core/soils/` (usle_k, sotwis, catchment), `core/geometry/rasterize.py` (pure-NumPy polygon rasteriser, identical to `gdal.RasterizeLayer` on 1.44 M test cells).
- Exchange schema: 18 soil fields on `catchments` and `soil_dataset` / `soil_depth_cm` in the metadata (additive, still `qeht-heas-1`; empty when no soil data is given).
- Tests: `tests/test_soils.py` (31 checks); QGIS smoke test extended (26 checks).

### Notes
- The Williams f_csand coefficient is 0.0256 (the SWAT 2009 theory PDF misprints 0.256).
- The SOTWIS Kenya shapefile contains invalid polygons; the soil tools read them without QGIS's validity check (which would abort) and report how many were rasterised as-is.

## [0.10.0] — unreleased (WP-A crossings, WP-D burn, D3 relink)
### Added
- **Road crossing candidates** (new "Road drainage" group): intersects every D8 flow link of the stream network near a road centreline with the alignment. Candidates carry chainage (from a start value, optionally reversed; multi-part roads chained in layer order), crossing angle, contributing area, Strahler order, reach id, flow side (L/R), road azimuth and `status = candidate`. Streams running alongside the road within the corridor for at least the minimum parallel length form one cluster and are exported as `parallel_reaches` (side-drain hints); an optional merge distance also clusters candidates by chainage. In each cluster the most downstream candidate is `recommended = 1` (D3). Nothing is deleted.
- **Burn crossings through embankments**: a short straight breach across the road at each crossing, lowered to a straight grade between the upstream and downstream low points (never raised), with a per-crossing log (cells, maximum cut, volume). Takes a candidate layer (road direction stored on each candidate) or points plus the road.
- **Renumber and relink exchange package** (D3): re-sorts edited crossings (by chainage, else downstream-first), re-issues gapless `outlet_uid`s, recomputes catchments and flow paths and writes a `renumber_log` table (old → new, unchanged / renumbered / new / deleted, distance moved) into a new package.
- **Build HEAS exchange package** accepts a candidate layer: uses the accepted candidates (else the recommended ones), keeps each outlet cell (no snapping), fills `chainage_m`, numbers IDs along the chainage and stores `crossing_candidates` and `road_alignment` in the package. Optional road alignment gives chainages to hand-placed points. New metadata keys `crossing_source`, `chainage_start_m`, `relinked_from` (additive, still `qeht-heas-1`).
- `core/network/alignment.py` (linear referencing, intersections), `core/network/crossings.py`, `core/conditioning/burn.py`, `core/linking/relink.py`.
- Tests: `tests/test_crossings.py` (43 checks); `tests/qgis_smoke.py` extended to the road workflow (23 checks, QGIS 3.34.4).

### Changed
- "Sequential numbering order" in Build HEAS exchange package now defaults to *Automatic* (along the chainage when known, else downstream first).

## [0.9.0] — unreleased (WP-B, HEAS interoperability)
### Added
- **Build HEAS exchange package** (new "Interoperability" group): one GeoPackage per run, schema `qeht-heas-1`, with `crossings`, `catchments` and `flowpaths` layers linked by `outlet_uid`, plus `qeht_run_metadata` (QEHT version, DEM path and SHA-256, CRS, cell size, thresholds, snapping, ID scheme, full parameters) and `qeht_field_dictionary` (meaning, unit, method and HEAS target of every field). Refuses a geographic CRS. Duplicate or empty IDs stop the run with a list. Optional CSV per layer.
- **`outlet_uid`**: a stable text identifier assigned once per run: from a chosen pour-point attribute, or sequential with a prefix (`X001` ...; downstream-first by default, or in pour-point order). Written by the characteristics, longest-flow-path and delineation tools as well. `outlet_id` (feature id) is still written but documented as not stable.
- Flow-path fields `lfp_L10_m`, `lfp_L85_m`, `lfp_z10_m`, `lfp_z85_m` so the 10–85 slope can be checked by hand.
- Catchment fields `catch_slope_horn` and `catch_relief_ratio` (the four slope domains are now `catch_slope_horn`, `catch_relief_ratio`, `lfp_slope`, `lfp_slope_1085`). `slope_mean` and `slope_relief_ratio` are still written as aliases for this release.
- `core/geometry/polygonize.py`: pure-NumPy mask-to-polygon tracer (OGC-valid output, identical to `gdal.Polygonize` + union in the tests). `core/interop/gpkg.py`: standard-library GeoPackage 1.2 writer/reader.
- Tests: `tests/test_interop.py` (74 checks, bare Python), golden exchange fixture in `tests/fixtures/`, and `tests/qgis_smoke.py` running every Processing tool inside QGIS (passes on QGIS 3.34.4).

### Changed
- **`lfp_slope_1085` now uses the conventional outlet-referenced definition**: points at 10% and 85% of the path length measured from the outlet, elevations interpolated linearly along the path. **Values from 0.8.3 and earlier are divide-referenced** (points at 10% and 85% from the divide, next-cell elevation). On concave profiles the new value is lower: across 40 longest flow paths on steep 30 m terrain it was a median 6% lower (range −20% to +5%). Re-run before reusing old 10–85 slopes. The old computation is kept as `path_slope_10_85_v083()` for comparison only.
- "Non-overlapping (local) catchments" in the characteristics tool are now computed upstream-first (ascending contributing area), so the result no longer depends on the order of the pour-point layer.
- Missing values (NaN) are written as NULL in the characteristics outputs instead of 0.

### Fixed
- **Catchment polygons with overlapping catchments** (characteristics tool, default full-upstream mode): catchments were burned into one label grid first-come-first-served, so a downstream catchment listed after an upstream one lost the upstream cells from its polygon, and a crossing listed after a downstream catchment that contained it got no polygon at all. The attribute values were correct. Each catchment is now traced from its own mask.
- NoData cells of the raw DEM could enter the elevation statistics as a real height (e.g. −9999 as the catchment minimum). They are now excluded.

## [0.8.3]
### Changed
- **QGIS 4 / PyQt6 compatibility: fixed all "QT6 Check" enum-scoping findings** reported by the plugin repository's `pyqgis4-checker` across `alg_watershed.py`, `alg_longest_flowpath.py`, `alg_fill.py`, `alg_characteristics.py`, `alg_streams.py`, and `alg_streamlines.py`. PyQt6 requires fully-scoped enum access (e.g. `QgsProcessing.SourceType.TypeVectorPoint`, `QgsProcessingParameterNumber.Type.Integer`, `QgsProcessingParameterNumber.Type.Double`) where PyQt5 accepted the legacy unscoped form; this check is informational rather than blocking for repository approval, but fixed anyway so the plugin loads cleanly under QGIS 4 without relying on PyQt6's temporary backward-compatibility shims.

## [0.8.2]
### Changed
- Documentation pass: replaced project-identifying site names with generic
  labels throughout `README.md`, `DOCUMENTATION.md`, `CHANGELOG.md`, and
  code docstrings ahead of public distribution. No functional changes.
- Added a second, synthetic example DEM (`examples/synthetic_dem.tif`) with
  no real-world correspondence, alongside the existing real-world sample
  (renamed to `examples/example_dem.tif`); `examples/README.md` now
  documents both.
- Terminology pass: replaced named-product comparisons (ArcGIS, ArcHydro,
  Spatial Analyst, ESRI) with generic descriptions ("a reference GIS
  platform", "a reference hydrology toolset", "standard D8") throughout
  the documentation and in code identifiers/docstrings (`encode_esri`/
  `decode_esri`/`ESRI_CODE` renamed to `encode_d8`/`decode_d8`/`D8_CODE`),
  ahead of public distribution. All validation numbers, percentages, and
  methodology are unchanged; the full 48-check test suite
  (`python -m qeht.tests.test_core`) still passes with zero regressions.
- Cleared the `experimental` flag for the initial public release, now that
  the plugin has been installed and exercised in live QGIS 3.44.6.

## [0.8.1]
### Changed
- Packaging pass for public distribution: finalized author/contact metadata,
  repository/tracker/homepage links, and license attribution ahead of
  submission to the official QGIS plugin repository and GitHub. No changes
  to the numerical core or Processing algorithms; the full 48-check test
  suite (`python -m qeht.tests.test_core`) still passes.

## [0.8.0]
### Added
- Iterated Barnes 2014 convergent flat resolver, verified against a
  synthetic saddle test (interior converges to a single outlet, no false
  sinks, monotonic downstream accumulation); resolves ~98% of flat cells
  in two iterations on real terrain.
- Flat-resolution method selector on the Flow Direction algorithm
  (toward-lower default; Barnes convergent optional).
- Saddle convergence regression test (test suite now 48 checks).
### Changed
- Documented the central empirical finding: on flat coastal terrain the
  canonical Barnes result disperses flow and matches the reference
  toolset's accumulation less well than the toward-lower method (stream
  IoU 0.19 vs 0.63); toward-lower retained as the validated default.
### Fixed
- Flat low-edge detection unified (adjacent to lower terrain OR to an
  equal-elevation draining cell) and exit-routing added so a flat's
  minimum-mask cell drains out instead of becoming a false sink.

## [0.7.1]
### Fixed
- Stale help text claiming Garbrecht & Martz as the default when the code
  ran toward-lower.
- Packaging: removed stray brace-expansion directories and a `.bak` file.
### Added
- Isolation validation: reference-platform flow direction into QEHT
  accumulation reproduces the reference platform's accumulation to
  99.31% exact over 9,000,000 cells, localising all flat-terrain
  discrepancy to flow direction.

## [0.7.0]
### Changed
- Validated flat resolution against a reference hydrology toolset's
  accumulation (stream IoU 0.60-0.67 on the flat Site B DEM); corrected
  earlier over-emphasis on raw per-cell direction agreement.
- Retired a broken experimental flat option from the UI.

## [0.5.0]
### Added
- Empirically recovered standard D8 tie-breaking rule (S > W > N > E,
  diagonals after), raising off-flat agreement with the reference
  toolset from ~84% to ~91% with no regression on steep terrain.

## [0.4.0]
### Added
- Catchment and flow-path characteristics tool (area, elevations, relief,
  mean terrain slope, relief ratio, longest-flow-path length/drop/slope,
  10-85 slope).

## [0.3.0]
### Added
- Stream-network vectorisation into polyline reaches with Strahler order,
  length, drop and slope.
- Multi-outlet longest flow path (one path per pour point, each bounded
  by its catchment; local/nested option).

## [0.2.0]
### Added
- Stream-based pour-point snapping (fixes catchment inflation from naive
  snap-to-maximum-accumulation).
- Undeclared-NoData audit and NoData override on Fill.
### Changed
- Level 3 validation against a reference hydrology toolset's Site A
  outputs (area median ratio ~0.99).

## [0.1.0]
### Added
- Initial release: Fill depressions, D8 flow direction, flow
  accumulation, stream network with Strahler ordering, catchment
  delineation, longest flow path. QGIS-independent numerical core with
  48 analytic tests; in-process GDAL I/O and polygonize.
