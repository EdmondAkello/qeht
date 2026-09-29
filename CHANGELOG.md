# Changelog

All notable changes to QEHT are recorded here. Versions follow the
development history of the numerical core and Processing tools.

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
