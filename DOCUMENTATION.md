# QEHT — QGIS Engineering Hydrology Toolkit

## Technical Documentation

**Version:** 0.8.2
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
      statistics.py      catchment and flow-path morphometry

  processing_provider/   the only QGIS-aware code
    provider.py          registers the eight algorithms
    base.py              shared base class and helpers
    alg_*.py             one file per algorithm

  tests/
    test_core.py         48 analytic checks, runnable without QGIS
```

The dependency direction is strict and one-way: `processing_provider` imports `core`; `core` never imports `processing_provider` or `qgis`. GDAL imports inside `core` are function-local and lazy, so the numeric modules remain importable where GDAL is absent.

---

## 3. Processing algorithms

QEHT registers eight algorithms under the "Engineering Hydrology" provider.

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

Each is a `QgsProcessingAlgorithm` registered through a `QgsProcessingProvider`. Exposing the tools this way — rather than as bespoke dialogs — means they gain input validation, batch mode, the Graphical Modeler, the history log, and `processing.run()` scriptability at no additional cost. Chaining tools in the Modeler is much of a commercial hydrology extension's practical value, and this design reproduces it.

---

## 4. Algorithms and equations

### 4.1 Depression filling

Filling uses the priority-flood algorithm (Barnes, Lehman & Mulla 2014): the DEM is flooded inward from its boundary using a priority queue keyed on elevation, and each cell is raised to the maximum of its own elevation and the spill level at which water reaches it. Complexity is O(N log N), a single pass, in contrast to the O(N²) behaviour of iterative neighbourhood-scanning fills.

Conditioning is deliberately a standalone step producing only a conditioned elevation surface — it does not emit flow directions. An optional minimum slope may be imposed across filled flats (an epsilon increment per traversal step) to give routing a defined gradient. An optional fill-depth raster is produced as a QA product: large contiguous fill depths typically mark a road embankment, dam, or culvert the DEM treats as a barrier, and signal that breaching or stream burn-in may be more appropriate than filling.

**NoData auditing.** QEHT inspects the DEM for undeclared NoData — for example a large block of cells at exactly 0.0 while genuine terrain begins far above it, with no NoData value set in the header. Treated as valid terrain, such a block becomes a spurious flat sink that the whole catchment drains into. The audit warns; the Fill algorithm accepts a NoData override so the user can mask it. (This was not hypothetical: it was found on a real project DEM with 11.8% of cells at 0.0 and real terrain starting at 1574 m.)

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

After steepest-descent routing, cells in filled flats have no downhill neighbour. QEHT provides two methods.

**Method A — toward-lower BFS (default).** Flat cells are seeded from their outlets and a breadth-first distance is propagated inward; each flat cell routes toward decreasing distance-to-outlet. This is the toward-lower component of the Garbrecht & Martz conceptual approach.

**Method B — Barnes 2014 convergent resolver (iterated).** A faithful reconstruction of RichDEM's two-gradient method: label each flat from its low edges; build a gradient *away* from high edges and a gradient *toward* low edges; combine them as `flat_mask = 2·toward + (flat_height − away)`; route each flat cell to the same-flat neighbour of lowest combined value. QEHT adds iteration to convergence — repeating the pass so nested flats drain in hierarchy order, since a flat whose outlet is itself a lower flat can only resolve once that lower flat has drained. The method passes synthetic saddle tests (interior converges to a single outlet, no false sinks, monotonic downstream accumulation) and resolves ~98% of flat cells in two iterations on real terrain.

**Why toward-lower is the default.** This is the central empirical finding of the project and is documented honestly rather than hidden. On a flat coastal DEM (34% of cells flat), the canonical Barnes result *disperses* flow and matches the reference toolset's accumulation **less** well than the simpler toward-lower method:

| Method | Stream IoU @1000 | p99 accum (reference 9,071) | cells > 1000 (reference 56,119) |
|---|---|---|---|
| Toward-lower | 0.63 | 6,324 | 51,974 |
| Barnes iterated | 0.19 | 642 | 14,428 |

Both are valid, cycle-free drainage solutions; they simply differ. The reference platform's actual flat behaviour — despite its documentation citing Garbrecht & Martz — empirically resembles the simpler method here. The remaining QEHT-vs-reference gap on flat terrain is therefore a genuine algorithmic difference, not an implementation defect. On steep terrain the two methods agree to within 0.1%.

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

For each catchment the tool reports area, highest/lowest/mean elevation, relief, and two distinct slopes; for each longest flow path it reports length, endpoint elevations, drop, whole-path slope, and 10–85 slope.

**Two catchment slopes, deliberately.** `slope_mean` is the mean terrain gradient over every cell by Horn's 3×3 method (the algorithm behind most commercial GIS `Slope` tools), which runoff-coefficient and curve-number tables assume. `slope_relief_ratio` is relief divided by longest-flow-path length, which is the "catchment slope" of most road-drainage manuals and the term Kirpich expects. On the Site A catchments these differ by factors of 1.2 to 4.5. They are not interchangeable; the report must state which was used.

---

## 5. Data handling and interoperability

- **Rasters** are read and written through the GDAL Python bindings directly. Rotated or skewed grids are rejected, because D8 assumes a north-up grid and silently accepting a rotated one produces plausible-looking but wrong directions.
- **Vectorisation** uses `gdal.Polygonize()` in-process — the call the standard workflow reaches by shelling out to `gdal_polygonize.exe`. No external executable is launched.
- **Flow-direction interchange** works in both directions via the standard D8 encoding, so a reference hydrology toolset's grid can be ingested and a QEHT grid exported.
- **CRS** is carried with each raster and pour points are reprojected into the DEM CRS as needed.

---

## 6. Validation record

QEHT is validated at three levels.

**Level 1 — synthetic analytic DEMs.** 48 checks in `tests/test_core.py`, runnable on bare Python + NumPy, covering encoding round-trips, distance weighting, fill spill levels, accumulation on analytic surfaces, catchment areas, longest-flow-path geometry, snapping, Strahler rules, and the Barnes saddle convergence test. All pass.

**Level 3 — reference hydrology toolset production output (Site A, steep).** 718×775 cells, EPSG:21037, 16 road-crossing pour points with reference catchments and longest flow paths. Catchment area median ratio 1.009 (15/16 within ±30%), longest-flow-path length median ratio 0.999 (15/16 within ±20%), best individual match 0.05%, zero flow-direction cycles. Cell-by-cell on an aligned grid, flow direction agrees with the reference toolset on 98.5%.

**Flat-terrain benchmark (Site B, coastal).** 34%-flat coastal DEM. Same-grid flow-direction agreement 88%; the residual is concentrated in flats and is a genuine algorithmic difference (Section 4.4). Accumulation validated in isolation to 99.31% exact against the reference platform (Section 4.5).

The two datasets bound the expected accuracy: near-exact on moderate-to-steep terrain, and a documented, quantified flat-terrain difference concentrated entirely in flow-direction assignment over flats.

---

## 7. Known limitations

- **No tiling.** The whole grid is processed in memory; several working copies are held during routing. A ~41-megapixel DEM needs well over 8 GB at peak. Above ~25 megapixels the tools warn and recommend clipping to the catchment plus a buffer, which is standard practice and does not affect results within the clip.
- **Flat terrain.** On very flat terrain QEHT reproduces the reference platform's output to a stream IoU of ~0.6, not near-unity. This is inherent to the algorithmic difference in Section 4.4, not a defect. For flat coastal catchments where the last increment matters, a reference-platform flow-direction grid can be ingested and QEHT's accumulation, catchment, and characteristics tools run on top of it.
- **Flat routing is one-sided by default.** The Barnes convergent option exists but is not the default because it matches the reference toolset less well on the tested terrain.
- **Single-threaded.** No parallelism; the algorithms are O(N) or O(N log N) but run on one core.

---

## 8. Release readiness (GitHub and QGIS Plugin Store)

Completed as of 0.8.2:

- README version header corrected to 0.8.2.
- `LICENSE` file added (GPL-2.0-or-later), author attributed to Edmond Akello.
- `CITATION.cff` added, author attributed to Edmond Akello.
- `metadata.txt`: author, email, and real homepage/repository/tracker URLs
  (`https://github.com/edmondakello/qeht`) populated; `changelog` pointer
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

See `CITATION.cff`. In brief: QEHT: QGIS Engineering Hydrology Toolkit, v0.8.2, GPL-2.0-or-later.

## 11. Licence

GNU General Public License v2 or later. See `LICENSE`.
