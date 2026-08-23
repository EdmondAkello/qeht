# Changelog

All notable changes to QEHT are recorded here. Versions follow the
development history of the numerical core and Processing tools.

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
