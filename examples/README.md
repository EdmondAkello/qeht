# Example datasets

Two small DEMs are bundled so the full Processing workflow can be exercised
without needing to supply your own data first.

## `example_dem.tif` — real-world terrain

A 718 × 775 cell DEM at 30.92 m resolution, EPSG:21037 (Arc 1960 / UTM zone
37S). This is the same raster used for QEHT's Level 3 validation against a
reference hydrology toolset (see `DOCUMENTATION.md`, Section 6, and
`README.md`, section "Validation") — real terrain with real flat areas,
undeclared NoData cells, and depressions, which is why it exercises the
flat-resolution and fill-depression logic more thoroughly than a synthetic
surface can.

## `synthetic_dem.tif` — synthetic terrain

A 320 × 300 cell procedurally generated DEM (EPSG:4326) with no real-world
correspondence — useful if you want to try the toolset on a dataset that
carries no location information at all. It includes a base tilt draining
toward a corner, a mix of ridges/valleys/hills, and three engineered
depressions specifically so **Fill depressions** has something to resolve.

## Suggested end-to-end workflow

Either file works with the steps below; substitute the filename you're using.

1. Load the DEM into QGIS.
2. Open the Processing Toolbox and expand **Engineering Hydrology**.
3. Run **Fill depressions** on the DEM to produce a hydrologically conditioned
   surface.
4. Run **D8 flow direction** on the filled DEM.
5. Run **Flow accumulation** on the flow-direction raster.
6. Run **Stream network (Strahler order)** with an accumulation threshold of
   roughly 500 cells (adjust to taste — lower thresholds produce a denser
   network).
7. Digitize one or more pour points near a road crossing or catchment outlet
   of interest, then run **Delineate catchment** with pour-point snapping
   enabled.
8. Run **Catchment and flow path characteristics** on the resulting catchment
   and streams to get area, relief, mean slope, relief-ratio slope, and
   longest-flow-path length/drop/slope — the inputs a time-of-concentration
   and design-discharge calculation needs.

Running this against `example_dem.tif` mirrors the exact validation workflow
described in `README.md` under "Validation", where QEHT's output on this DEM
was compared cell-by-cell and catchment-by-catchment against a reference
hydrology toolset's production results.
