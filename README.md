# QEHT — QGIS Engineering Hydrology Toolkit v0.20.0

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
        interop/             design hydrology package: field dictionary, GeoPackage writer
        network/             road alignment chainage, crossing candidates, ground profile
        watershed/morphometry.py  flow-path segments, basin shape and network indices
        watershed/channel.py approach / exit channel slopes, deposition indicator
        runoff/curve_number.py   curve number and Rational C per catchment
        conditioning/burn.py breach road embankments at crossings
        soils/               soils from any source (SoilGrid), SOTWIS / table / raster loaders,
                             USLE K, texture, HSG proxy, HYSOGs
        geometry/rasterize.py polygons -> grid (cell-centre rule, pure NumPy)
      processing_provider/   ← the only place QGIS and core meet
      tests/test_core.py     ← runs on bare Python + NumPy
      tests/test_interop.py  ← outlet_uid, slopes, exchange package, golden fixture
      tests/test_crossings.py ← crossing candidates, burn, renumber-and-relink
      tests/test_soils.py    ← USLE K, texture, SOTWIS loader, soil block
      tests/test_flats.py    ← vectorised core == v0.8.3 reference; Barnes == RichDEM port
      tests/test_alignment.py, test_morphometry.py, test_soils_any.py,
      test_runoff.py, test_channel.py  ← v0.15 features
      tests/qgis_smoke.py    ← every Processing tool, run inside QGIS

The core is importable without QGIS. That is what makes the hydrology testable:

    python -m qeht.tests.test_core        # 47 analytic checks
    python -m qeht.tests.test_interop     # 77 checks incl. the golden fixture
    python -m qeht.tests.test_crossings   # 45 checks: road crossings, burn, relink
    python -m qeht.tests.test_soils       # 33 checks: soils, USLE K, CSV soil table
    python -m qeht.tests.test_flats       # 36 checks: oracles, Barnes == RichDEM, edge drains, DEM QA
    python -m qeht.tests.test_erosion     # 50 checks: erosion indices, RUSLE, classes, A14 anchors
    python -m qeht.tests.test_alignment   # 15 checks: alignment ground profile
    python -m qeht.tests.test_morphometry # 22 checks: overland/channel split, basin shape
    python -m qeht.tests.test_soils_any   # 18 checks: soils from any source, HSG / K overrides
    python -m qeht.tests.test_runoff      # 17 checks: curve number, AMC, Rational C
    python -m qeht.tests.test_channel     # 16 checks: STI, deposition indicator, channel slopes
    python -m qeht.tests.test_coverage    # 21 checks: missing crossings, sags, flat stretches, proposed crossings
    python -m qeht.tests.test_rainfall    # 35 checks: rainfall zones, mean annual rainfall, R estimate
    python -m qeht.tests.test_report      # 15 checks: characteristics table, run report
    python -m qeht.tests.test_floodplain  # 19 checks: floodplain widths (profile, HAND)
    python -m qeht.tests.test_quicklooks  # 14 checks: PNG, world file, legends, downsampling

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
| Build design hydrology package | — (one self-describing GeoPackage of crossings, catchments and flow paths; see below) |
| Renumber and relink exchange package | — (re-issue IDs after editing crossings) |
| Road crossing candidates | — (road × drainage crossings with chainage, clustering) |
| Burn crossings through embankments | a DEM-reconditioning "burn culverts" step |
| Soil, rainfall and runoff parameters for catchments | zonal soil statistics, USLE K, HSG shares, curve number, rainfall zone, mean annual rainfall |
| Erosion indices and RUSLE soil loss | SPI, TWI, LS, RUSLE and severity classes |
| Sample erosion along alignment | — (erosion stations and reaches along a road) |
| Alignment ground profile | — (ground, fill, area and stream crossings every 10 m along a road) |
| Drainage coverage check along a road | — (missing crossings, sag points, flat stretches) |
| Run hydrology pipeline (one click) | — (the whole chain into one output folder, with a run report) |
| Run report from a design hydrology package | — (HTML report and characteristics table for any package) |

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
| `catch_relief_ratio` | relief / longest flow path length (variant of Schumm's relief ratio, which uses basin length) |
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
gradient. Time-of-concentration formulas (Kirpich, Bransby-Williams, TRRL) take
the slope along the longest flow path, `lfp_slope` or `lfp_slope_1085`, not
either catchment slope. The relief ratio is a basin-steepness index. Picking the
wrong slope changes a design discharge substantially. **State which you used.**

> **Re-running a design dataset.** Two changes since 0.8.3 can move design
> inputs: the 10–85 slope is measured from the outlet since 0.9.0 (values a
> median 6 % lower on steep 30 m terrain), and Barnes is the default flat method
> since 0.13.0 (flow paths and catchments can move on flat terrain). To re-run:
> 1. re-run Longest flow path and Catchment characteristics;
> 2. compare old and new `lfp_slope_1085` (`path_slope_10_85_v083()` gives the
>    old value for the same path);
> 3. check the outlets flagged `flat_sensitive`;
> 4. record the QEHT version and flat method in the design report (both are
>    printed in the tool log and stored in `qeht_run_metadata`).

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

## Design hydrology package

**Build design hydrology package** (called "Build HEAS exchange package" up to
0.14; the algorithm id `buildheasexchange` is unchanged, so saved models keep
working) writes one self-describing GeoPackage per run (schema `qeht-heas-1`)
for spreadsheets, reports or design software:

| Table | Content |
|---|---|
| `crossings` | snapped outlet points with `outlet_uid`, input and snapped coordinates, snap distance, contributing area, Strahler order, approach / exit channel slopes |
| `catchments` | one polygon per crossing, with all catchment fields |
| `flowpaths` | longest flow path per crossing, with all flow-path fields |
| `qeht_run_metadata` | QEHT version, DEM path and SHA-256, CRS, cell size, thresholds, snapping, ID scheme, full parameters |
| `qeht_field_dictionary` | meaning, unit, method, plain-language downstream use (and the HEAS field, where there is one) of every field |
| optional layers | `road_alignment`, `alignment_profile`, `crossing_candidates`, `burn_log`, `renumber_log` when their inputs are given |

The three layers share `outlet_uid`, so any reader links them by key rather
than by row order or by position. Rules: a projected, metric CRS is required
(a DEM in degrees is refused); duplicate or empty IDs stop the run with a list;
missing values are written as NULL, never 0; schema changes are additive only
(a rename or removal would be `qeht-heas-2`). An optional CSV per layer carries
`outlet_uid` for spreadsheet users. The GeoPackage is written with the Python
standard library and checked against GDAL's GeoPackage validator in the tests.

### Using with HEAS (optional)

HEAS is a separate desktop hydrology–hydraulics design tool that can import
this package directly; it is not required. `tests/fixtures/golden_exchange.gpkg`
is the shared fixture HEAS imports in its own tests (regenerated in 0.15 with
new fields only; no existing value changed).

## Road drainage workflow (v0.10)

1. **Road crossing candidates** — every D8 flow link of the stream network
   near the road is intersected with the centreline (the same segments QEHT's
   stream polylines are made of, so no crossing is lost to a digitising gap).
   Each candidate carries chainage, crossing angle (90 = square), contributing
   area, Strahler order, reach id, the side the flow comes from, and
   `status = candidate`. Streams that run alongside the road inside the
   corridor half-width for at least the minimum parallel length cross the
   centreline many times on a DEM: those candidates share a cluster and the
   stream is written as a *parallel reach* (where a side drain must carry the
   water). In each cluster the most downstream candidate (largest area) is
   `recommended = 1`. Nothing is deleted; set `status` to accepted/rejected.
2. **Burn crossings through embankments** (optional, for DSMs and survey DEMs
   that show the embankment as a dam) — at each crossing, a short straight
   breach across the road from the low point upstream to the low point
   downstream, lowered to a straight grade (never raised). Every breach is
   logged (cells, maximum cut, volume). Then fill, flow direction and
   accumulation on the burned DEM.
3. **Build design hydrology package** on the candidate layer — uses the accepted
   candidates (or the recommended ones), keeps each outlet cell exactly (no
   snapping), numbers `outlet_uid` along the chainage (decision D2), and stores
   the full candidate layer and the road alignment in the package.
4. **Renumber and relink** — after deleting, moving or adding crossings in the
   package, re-sorts them, re-issues gapless IDs, recomputes every catchment
   and flow path, and writes `renumber_log` (old → new, including deleted and
   new, and how far a moved point moved) into a new package.

## Soils (v0.11)

**Soil, rainfall and runoff parameters for catchments** (and the optional soil input on
Build design hydrology package) adds a soil block to every catchment: area-weighted topsoil
sand, silt, clay, organic carbon, coarse fragments and bulk density; USDA
texture; dominant FAO drainage class; **USLE K** (Williams/EPIC, SI units
t·ha·h/(ha·MJ·mm)) plus the Renard Dg-based alternative; the coarse-fragment
factor CFRG; a texture-and-drainage **hydrologic soil group proxy**; and the
share of the catchment covered by soil data.

- **SOTWIS (Kenya preset).** Give the SOTWIS polygons and the SQLite database
  (never the .mdb). Each soil unit's components (up to ten profiles with their
  shares) are depth-weighted over the chosen interval (default 0–20 cm) and K
  is computed per component, then weighted by share, then by catchment area.
  SOTWIS missing values (−1) are left out and the shares renormalised.
- **Other soil maps.** Polygons carrying SOTWIS fields (dominant soil only) or
  plain `sand`, `silt`, `clay`, `oc` (%) [+ `bulk`, `cfrag`, `drain`].
- **Polygons + a soil table (v0.14).** Polygons that carry only a unit code,
  plus a CSV (comma, semicolon or tab) with the unit code and `sand`, `silt`,
  `clay`, `oc` (%) [+ `bulk`, `cfrag`, `drain`]; give the unit field name.
  Duplicate units or missing columns stop the run.
- The Williams f_csand coefficient is 0.0256; the SWAT 2009 theory PDF prints
  0.256, which pins f_csand at 0.2 for sandy soils.
- The soil map is rasterised on the DEM grid by cell centre (identical to
  `gdal.RasterizeLayer` in the tests), so catchments are weighted by the same
  cells as the hydrology. Invalid soil polygons (present in the SOTWIS Kenya
  shapefile) are rasterised as they are and counted in the log.
- Across Kenya, SOTWIS-derived K has a median of 0.029 against 0.022–0.024 for
  the ESDAC global K rasters over the same cells (correlation 0.24–0.40); on
  six steep Rift-valley test catchments it was 0.026–0.038 against 0.027–0.028
  (ESDAC Wischmeier-based K). Treat K as an estimate with that spread.

## Erosion (v0.14)

**Erosion indices and RUSLE soil loss** writes one folder with the terrain
indices, an optional RUSLE soil-loss raster, severity class rasters (with
colours, class names and a QGIS style) and `erosion_run.json`, which records
every factor's source and every setting.

- **Indices.** Slope by Horn on the raw DEM; specific catchment area
  A_s = (upslope cells + 1) × cell area / cell size; SPI = A_s·tan β;
  ln(SPI); TWI = ln(A_s / tan β); LS by Moore & Burch (m = 0.4, or m by slope
  class as in the A14 study) or Desmet & Govers with McCool m and S. tan β is
  floored at 0.001 inside the indices only.
- **RUSLE** A = R·K·LS·C·P (t/ha/yr). Each factor is a raster, a single value
  or a dataset: R from a raster (GloREDa 2023 annual erosivity recommended);
  K from SOTWIS (Williams/EPIC, 0–20 cm) by default; C from ESA WorldCover 2021
  through an editable class lookup, **flagged as a land-cover proxy**; P = 1.0
  unless a raster, value or the WorldCover P lookup is given. Rasters on other
  grids are resampled in-process onto the DEM grid. **If R, K or C is missing,
  no soil loss is computed**: QEHT never invents a factor and falls back to
  LS-only classes labelled "terrain potential".
- **Classes.** SPI by percentiles of the extent (relative) or fixed ln(SPI)
  0 / 5 / 10; RUSLE by 5 / 12 / 25 / 50 or 5 / 10 / 20 / 40 t/ha/yr, or your own
  breaks; the combined class is the higher of the two severity scores unless
  you give a 5 × 5 matrix. `class_extents.csv` gives area per class, split into
  channel and hillslope cells.

**Sample erosion along alignment** reads that folder and a road centreline
(a rough digitised line is fine; chainage is provisional until the geometric
design exists). Stations every 10 m carry max and mean ln(SPI), LS, soil loss
and TWI on each side within a buffer; consecutive stations with the same worst
class form reaches whose lengths add up to the road length. Optional chainage
profile chart.

**Build design hydrology package** takes the folder as an optional input and adds
an erosion block: per catchment the mean and p90 soil loss and LS, channel
ln(SPI) p90, area-weighted K, C and P, class shares, gross soil loss, sediment
delivery ratio (SDR = 0.565·A^−0.125) and sediment volume; per crossing the
ln(SPI) at and above the outlet, local slope, LS and TWI, and the hydrodynamic
impact score 0.4·SPI + 0.3·RUSLE + 0.3·sediment (Low / Moderate / High /
Severe) with an indicative mitigation. The scheme, weights and C lookup follow
the published A14 corridor study (Akello & Omosa 2025), whose tables are the
test anchors.

## v0.15: alignment profile, flow-path segments, soils anywhere, curve numbers

**Alignment ground profile** (Road drainage) stations a road centreline every
10 m: DEM ground (bilinear), the filled-DEM level and the ponding depth the
fill removed, contributing area (largest within one cell), the stream crossings
(the exact D8 link × centreline intersections of Road crossing candidates) with
their Strahler order, and the longitudinal ground slope. Provisional levels
before the geometric design exists, and a check on DEM height bias afterwards.
Build design hydrology package writes the same profile as `alignment_profile`
when a road alignment is given.

**Flow-path segments.** Each longest flow path is split at the channel head
(the first cell reaching the stream threshold): `lfp_overland_m` / `_slope`,
`lfp_channel_m` / `_slope` / `_slope_1085`, and the overland part as sheet flow
(up to a cap, default 100 m) plus shallow concentrated flow — inputs for Kerby
and segmental (TR-55) Tc. `lfp_overland_m + lfp_channel_m = lfp_length_m`.

**Basin shape and network** (information only): `perimeter_km` (smoothed
outline, not the cell staircase), `form_factor`, `elongation_ratio`,
`circularity_ratio`, `drainage_density`, `stream_frequency`, `max_strahler`.

**Channel slopes at every crossing:** `ch_slope_us` along the main stem
upstream and `ch_slope_ds` along the D8 path downstream, over 200 m on the raw
DEM, with the lengths actually used.

**Soils from any source.** Every source becomes the same per-cell description
on the DEM grid, so the catchment soil fields mean the same thing everywhere:

- SOTWIS / SOTER database, polygon attributes, or polygons + CSV (as before;
  results identical to 0.14);
- any soil map: choose its sand, clay and organic carbon fields (OC in % or
  g/kg; silt = 100 − sand − clay when not given);
- a soil unit raster + CSV table (e.g. HWSD v2 mapping units with an exported
  attribute table);
- texture rasters, or a folder of SoilGrids 2.0 files (`sand_0-5cm_mean.tif`
  …), converted from SoilGrids units and depth-weighted over the chosen
  interval;
- a hydrologic soil group raster (HYSOGs250m codes; dual groups A/D, B/D, C/D
  counted as D) or field, and a K raster or field, which override the
  texture-based values. `soil_hsg_source` and `usle_k_source` record which was
  used; `hsg_pct_a` … `hsg_pct_d` give the shares.

**Curve number and Rational C.** With a land-cover raster (ESA WorldCover
classes): per cell land cover × hydrologic soil group → CN from TR-55
Table 2-2 (condition fair / good / poor), averaged over the catchment (`cn_ii`)
and exported at AMC I, II or III (`cn_export`). The WorldCover → TR-55 match is
a **proxy**, flagged in `runoff_json`; replace it with your own lookup CSV.
Rational C comes only from a lookup you supply — none ships. Land-cover shares
`lc_pct_*` are always given.

**Sediment transport and deposition.** The erosion tool also writes
`sti_overland.tif`, the sediment transport capacity index of Moore & Wilson
(1992) on overland cells (channels NoData, A_s capped at 100 m) — a relative
indicator of where slopes deliver sediment, not the RUSLE LS factor and not used
in soil loss, the classes or the composite. In the package each crossing gets a
**deposition indicator**: median SPI on the main channel 0–100 m vs 100–500 m
upstream; a ratio below 0.7 means transport capacity falls into the inlet
(deposition-prone), above 1.3 scour-prone.

## v0.16: drainage coverage check

With a road alignment, **Build design hydrology package** checks the road for drainage that has no crossing, from the DEM alone (the standalone **Drainage coverage check along a road** gives the same findings without delineating):

- **Missing crossings** — a stream of at least the stream-threshold area crosses the centreline with no crossing within 50 m. In the package each becomes a **proposed crossing** (`status = proposed`, `P001…`), delineated and characterised like the existing ones.
- **Sag points** — low points of the ground profile with no stream or crossing nearby. The water ponding against the embankment is measured with the centreline as a wall (open at the crossings). A sag with ≥ 0.05 km² becomes a proposed crossing (equaliser / relief).
- **Flat stretches** — ≥ 300 m flatter than 0.5 % along and across the road: where relief culverts for sheet flow may be needed. QEHT lists them; it does not place or size reliefs.

Adopt a proposed crossing by setting its `status` to `existing`, delete the ones you don't want, then run **Renumber and relink** — proposed ones keep the P prefix, adopted ones join the X sequence.

## v0.17: one-click pipeline, run report, rainfall

**Run hydrology pipeline (one click)** (group Workflow) runs the whole chain and writes one folder:

| Folder | Contents |
|---|---|
| `layers/` | crossings, catchments, flow paths, streams, crossing candidates, coverage findings, sags, flat stretches, alignment profile (GeoPackage each) |
| `rasters/` | filled DEM, flow direction, accumulation, streams, Strahler order, `erosion/` |
| `tables/` | `catchment_characteristics.csv` (one row per crossing: chainage, area, the four slopes, flow-path segments, channel slopes, flat-method check, soils, CN, rainfall, shape, erosion) and one CSV per layer |
| `quicklooks/` | PNG + world file + legend JSON of the relief, accumulation and erosion rasters (v0.18) |
| `report/` | `run_report.html`: summary, schematic plan, crossing schedule, DEM checks, flat-method check, coverage findings, inputs and provenance, warnings, CN lookup |
| `package/` | `design_hydrology.gpkg`, only when *Keep the design hydrology package* is ticked |
| `settings.json` | every parameter of the run |

Steps: DEM checks → fill → D8 (Barnes) → accumulation → streams → candidates (with a road) → optional burn and second routing → erosion → crossings, catchments, flow paths with every optional block → coverage check → flat-method check → tables and report. Each step is the standalone tool, run as a child algorithm, so results are identical to running the tools one by one (the smoke test checks this).

- **Review stop.** *Run = Stop after crossing candidates* writes the candidates and stops. Edit `layers/crossing_candidates.gpkg` (status accepted / rejected, move or add points), then run again with *Re-run from settings* = the run's `settings.json`, *Crossings* = the edited layer and *Run = Full run*.
- **Re-run from settings.** Every parameter comes from the file; a new output folder, Crossings and Run override.
- **Undeclared NoData.** *Treat this DEM value as NoData* writes `rasters/dem_input.tif` with the value declared, and every step uses it.

**Run report from a design hydrology package** writes the same report (and the characteristics CSV) for any package, e.g. one built by hand or relinked.

**Rainfall (F3)**, in *Soil, rainfall and runoff parameters*, the package and the pipeline:

- **Zones:** a polygon layer and its name field. Per catchment `rain_zone` (largest share; ties go to the name that sorts first), `rain_zone_pct`, `rain_zones_json` (every zone's share) and `rain_zone_coverage_pct`. A warning lists catchments whose zone covers less than 80 %: pick their design zone by hand.
- **Mean annual rainfall:** a raster in mm/yr (e.g. a CHIRPS climatology), bilinear to the DEM grid: `map_mm`, `map_coverage_pct`, `map_dataset`. Values that cannot be annual mm (mm/day, tenths of mm, negative) are refused.
- **R estimate (optional):** `rusle_r` per cell from the rainfall with Renard & Freimund (1994) or Lo et al. (1985), then the catchment mean, labelled an ESTIMATE in `rusle_r_method` and the metadata. Neither relation comes from East Africa; an erosivity raster such as GloREDa is better. The erosion tool can use the same estimate when no R raster or value is given.

## v0.18: floodplain width and raster quicklooks

**Floodplain width indicator (A4).** With a road, *Build design hydrology package* (and the pipeline) measures at every crossing of at least 10 km² how wide the valley floor is along the road: the run of the alignment profile around the lowest ground within 50 m of the crossing where the ground is below bed + 0.5, 1 and 2 m (`fp_w_0p5_m`, `fp_w_1p0_m`, `fp_w_2p0_m`, extent `fp_ch_from_m`–`fp_ch_to_m` at 1 m, `fp_z_bed_m`), and the same by height above nearest drainage (`fp_hand_w_*`). A width that runs off the end of the profile or into NoData is a lower bound and says so in `fp_note`. It informs the split of the check flood between the main structure and relief culverts. It is a terrain indicator, not a flood level.

**Raster quicklooks (A6).** Give a quicklooks folder (the pipeline does by default) and the package writes, for the relief (hillshade over elevation), flow accumulation (log stretch) and the erosion class and STI rasters, a downsampled PNG (≤ 4,096 px), its world file (`.pgw`) and a legend JSON, and indexes them in the table `rasters`. Classified rasters keep their class values (nearest neighbour) and the colours of their `.qml` styles; continuous ones state their ramp and stretch. The run report draws its plan over the terrain.

## v0.19: channel sections at crossings

*Build design hydrology package* (and the pipeline) cuts a cross-section of the natural channel 30 m downstream of every crossing, existing and proposed, for the tailwater rating, the bridge waterway and the sediment-continuity check. The section sits on the D8 path and runs perpendicular to the local flow direction (taken from two cells upstream to two cells downstream, which smooths D8's 45° steps). The raw DEM is sampled bilinearly every half cell over ±150 m, offset positive on the right looking downstream.

- `xs_bed_m`: lowest ground within one cell of the section centre (`xs_bed_offset_m`).
- `xs_bankfull_w_m`, `xs_bankfull_d_m`: bank-full width and depth by break of slope. The bank top on each side is where the side slope over the next two cells falls below 1:20 or to less than half the slope behind; the lower bank governs.
- `xs_w_0p5_m`, `xs_w_1p0_m`, `xs_w_2p0_m`: width below bed + 0.5, 1 and 2 m. These hold when there is no clear bank.
- `xs_side_slope_l`, `xs_side_slope_r`: side slopes as H:V, least squares from the bed to the bank top (bed + 1 m without a bank).
- `xs_station_elev_json`: the transect as `[[station, z], …]`, and the layer `xs_transects`.
- `xs_quality`: `low` for a bank-full width under 3 cells or no bank on either side; `medium` for one bank only or a transect cut by NoData or the grid edge; `high` otherwise. `xs_note` says why.

Indicative only; on a 30 m DEM a small channel is below the grid resolution. Use survey where available. On a coarse DEM the break of slope often finds the valley shoulder rather than the channel bank; a bank-full depth of several metres says which. Distance, half width and bank slope are advanced parameters (settings in `xs_params_json`); the section can be switched off. The characteristics table has a "channel section" group, and the run report lists low-quality sections under DEM checks.

## v0.20: side-drain siltation indicator

Side drains silt up where a low-gradient drain is fed by a slope with a high sediment-transport capacity. *Sample erosion along alignment* now measures that per side of the road when you also give the **flow direction** (and, for the flag, the **raw DEM**):

- At each station, the cells sampled on each side (10 … 50 m out) are kept when their D8 path reaches the road strip (5 m either side of the centreline) without entering a channel. Those cells are the slope that feeds the side drain.
- Station fields: `sti_p50_lhs`, `sti_p90_lhs`, `n_lhs` and the same for `rhs` (left / right looking up-chainage), from the overland STI raster of *Erosion indices and RUSLE*; NULL when no cell drains to that side.
- Classes `sti_class_lhs` / `_rhs` (low, moderate, high, very high) come from the 50 / 75 / 90th percentiles of this corridor's own station values. They describe relative transport capacity, not severity.
- `slope_long_pct` is the longitudinal ground slope from the alignment profile, and `siltation_lhs` / `_rhs` = 1 where the class is high or very high and the slope is under 1 %. A screening flag.
- Reach fields: `sti_p90_lhs` / `_rhs` (p90 of the station values), `sti_max_*`, `sti_class_*` and `siltation_len_lhs_m` / `_rhs_m`. The chart gains an LHS / RHS panel.

*Build design hydrology package* writes the same reaches as layer `corridor_sti` when an erosion folder and a road are given (settings and class breaks in `corridor_sti_params_json`). The indicator is not used in the erosion classes, RUSLE or the composite score.

### Data sources by region

QEHT does not ship data; it reads what you have. Global sources that work
anywhere (national data is better where it exists):

| Input | Global source | Tool input |
|---|---|---|
| DEM | FABDEM, Copernicus GLO-30, ALOS AW3D30 (30 m) | Raw DEM |
| Soil texture | SoilGrids 2.0 (ISRIC, 250 m, six depths) | SoilGrids folder or texture rasters |
| Soil units | HWSD v2.0 (FAO/IIASA, ~1 km) | Soil unit raster + CSV |
| Hydrologic soil group | HYSOGs250m (ORNL DAAC) | HSG raster |
| Land cover | ESA WorldCover 2021 (10 m) | Land cover / WorldCover |
| Rainfall erosivity R | GloREDa (ESDAC) | R raster |
| Soil erodibility K | ESDAC global K, or computed from texture | K raster |

Kenya-specific defaults stay as they were: SOTWIS for soils, the A14 study's
WorldCover C lookup. The default RUSLE class breaks (5 / 12 / 25 / 50 t/ha/yr)
come from East African practice — check them against local guidance or give
your own.

## Flats and performance (v0.12)

**Vectorised core.** Fill, D8 flow direction, both flat resolvers,
accumulation, Strahler order, catchment delineation, longest flow path and
pour-point snapping were per-cell Python loops; they are now whole-array NumPy
operations. The v0.8.3 implementations are kept in `core/flow/_reference.py`
and the tests require identical output (bit for bit; fractional-weight
accumulation to 1e-12). Measured on 30 m clips of about one million cells, the
full chain (fill, flow direction, accumulation, Strahler, catchment, longest
flow path) went from 10.5 s to 2.9 s on steep terrain and from 25.9 s to 4.8 s
on a flat coastal clip; D8 peak memory fell from ~176 to ~59 bytes per cell.

**Barnes resolver fixed.** Up to 0.11 the Barnes option could route two flat
cells into each other; fixed in 0.12.

### 0.13: Barnes is the default flat method (WP-G benchmark)

80 random 15 km areas across Kenya (20 each flat, rolling, hilly,
mountainous), each on ALOS AW3D30 and FABDEM, were routed with QEHT and with
TauDEM, RichDEM, MAS and DDM HydroLogic. Barnes agreed best with the
independent TauDEM (stream F1 0.969 vs 0.948 for toward-lower; 0.870 vs 0.733
on whole-metre flats), and QEHT's Barnes now matches Barnes' own
implementation (RichDEM) cell for cell. The hybrid reproduced Barnes and was
removed; toward-lower stays as an option. The benchmark also fixed flats that
touch the clip edge or NoData, which now drain to it. Details: `BENCHMARK.md`.

## Interoperability

Flow direction is written and read in the standard D8 encoding
(E=1, SE=2, S=4, SW=8, W=16, NW=32, N=64, NE=128; 0 = sink or edge outlet).

This works in both directions. You can feed a flow direction grid produced by
another GIS/hydrology toolset straight into QEHT's accumulation tool — which
is the cleanest way to validate this toolkit against a trusted production
result.

## Known divergences from commercial reference toolsets

Stated explicitly because a drainage report needs them stated:

1. **Flat resolution.** The default (since 0.13) is Barnes 2014, both
   gradients, identical to Barnes' RichDEM implementation; routing toward
   lower terrain only is an option (see `BENCHMARK.md`). On wide flats the
   two give different networks, and commercial reference toolsets may differ
   from both. Filling with a small minimum slope (1e-4)
   largely removes flats before routing.
2. **Tie-breaking.** Where neighbours give identical drop/distance, the fixed
   priority S, W, N, E, SE, SW, NW, NE applies - recovered empirically from a
   reference platform (off-flat agreement 84 % → 91 %). Deterministic.
3. **No tiling.** DEMs are processed whole in memory. Every tool reports an
   estimated peak memory before it runs (measured per-step figures) and warns
   when it approaches the free RAM; clip to your catchment plus a buffer when
   it does — normal practice anyway.

## Validation

**Level 1 — synthetic analytic DEMs.** 47 checks, `tests/test_core.py`, plus 77 interop checks in `tests/test_interop.py`, 43 road-crossing checks in `tests/test_crossings.py` 33 soil checks in `tests/test_soils.py` and 36 checks in `tests/test_flats.py` (vectorised core vs the v0.8.3 reference; Barnes vs a port of RichDEM; edge drainage; resampling audit and flat-method sensitivity) and 44 erosion checks in `tests/test_erosion.py` (analytic plane and valley, the A14 2025 paper's tables). PASS. `tests/qgis_smoke.py` (37 checks) runs every tool inside QGIS.

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
3. **Nearest-neighbour resampling (0.13.1, Site C road project).** A
   FABDEM clip reprojected with nearest neighbour repeated ~1–2 % of its rows
   and columns; 14 % of cells were exactly tied with a neighbour. One tied cell
   on a river sent ~176 km² to one culvert or the next depending on the flat
   method. Fill and D8 flow direction now warn when rows or columns repeat
   (`audit_resampling()`), and Catchment characteristics reports every
   outlet's area under both flat methods and flags the ones that differ
   (`flat_sensitive`). On a bilinear resample of the same DEM both methods
   agreed at every crossing.

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
