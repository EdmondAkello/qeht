# QEHT development history

Records moved out of the README. They describe how decisions were reached and are superseded by the current documentation ([README.md](README.md), [BENCHMARK.md](BENCHMARK.md), [CHANGELOG.md](CHANGELOG.md)).

## The flat-resolution investigation (v0.5 to v0.8)

The paragraphs below record how the flat method was chosen before the WP-G
benchmark. They are kept for the record. They are superseded: since 0.13 Barnes
is the default (see *Flats and performance* above and `BENCHMARK.md`), and the
v0.8 Barnes figures were affected by a resolver defect fixed in 0.12.

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
drainage solutions; they simply differ, and toward-lower was closer to the
reference toolset here (at v0.8). On
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
(At v0.8 the Barnes code was kept out of the tools; since 0.13 it is the default.)

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
