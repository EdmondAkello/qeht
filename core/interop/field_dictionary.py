# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Field dictionary for the QEHT -> HEAS exchange package (schema qeht-heas-1).

This list is the contract. The writer creates the layers from it, the
validator checks against it, and it is written verbatim into the
`qeht_field_dictionary` table of every exchange GeoPackage so the file
documents itself.

Additive evolution only: new fields may be appended under qeht-heas-1.
Renaming or removing a field requires qeht-heas-2 (see heas_exchange.py).

Each entry: (field, type, unit, meaning, method, heas_target)
type is one of "text", "int", "real".
"""

SCHEMA_VERSION = "qeht-heas-1"
SCHEMA_MAJOR = 1

_LINK = [
    ("link_method", "text", "-", "how this feature was linked to its crossing",
     "'pour_point' = generated from that crossing in this run; 'spatial' / "
     "'manual' reserved for reconciled layers", "link registry"),
    ("link_confidence", "real", "0-1", "confidence of the link",
     "1.0 for pour_point links", "link registry"),
    ("link_note", "text", "-", "free-text note on the link", "", "link registry"),
]

CROSSINGS = [
    ("outlet_uid", "text", "-", "stable crossing/outlet identifier; identical on the "
     "crossing, its catchment and its longest flow path",
     "assigned once per run (see qeht_run_metadata id_scheme)", "link key"),
    ("crossing_id", "text", "-", "crossing identifier (default = outlet_uid; editable)",
     "", "crossing ID"),
    ("outlet_id", "int", "-", "pour-point feature id in the input layer - NOT stable, "
     "kept for backward compatibility only", "QGIS feature.id()", "none"),
    ("source_id", "text", "-", "value of the user-selected ID attribute, if any",
     "", "provenance"),
    ("input_x", "real", "m", "pour point as supplied (DEM CRS)", "", "QA"),
    ("input_y", "real", "m", "pour point as supplied (DEM CRS)", "", "QA"),
    ("outlet_x", "real", "m", "snapped outlet (cell centre, DEM CRS)",
     "nearest stream cell within the snap radius", "QA"),
    ("outlet_y", "real", "m", "snapped outlet (cell centre, DEM CRS)", "", "QA"),
    ("snap_dist_m", "real", "m", "distance moved by snapping", "", "QA"),
    ("acc_at_outlet_km2", "real", "km2", "contributing area at the snapped outlet "
     "read from the accumulation raster (equals area_km2 for full catchments)",
     "(accumulation + 1 cell) x cell area; QEHT accumulation excludes the cell itself",
     "QA"),
    ("stream_order", "int", "-", "Strahler order at the outlet (if a stream-order "
     "raster was supplied)", "Strahler", "info"),
    ("chainage_m", "real", "m", "position along the road alignment (from the crossing-"
     "candidate tool; empty for hand-placed pour points)", "linear referencing",
     "crossing chainage"),
] + _LINK + [
    # erosion block at the crossing (v0.14, WP-F)
    ("ero_lnspi_max3x3", "real", "ln(m)", "highest ln(SPI) in the 3x3 cells at the outlet", "",
     "gully potential"),
    ("ero_lnspi_app_p50", "real", "ln(m)", "median ln(SPI) on the approach channel",
     "channel cells up to 10 D8 steps upstream", "gully potential"),
    ("ero_lnspi_app_p90", "real", "ln(m)", "90th percentile ln(SPI) on the approach channel", "",
     "gully potential"),
    ("ero_slope_pct", "real", "%", "local slope, 3x3 mean", "Horn, raw DEM", "info"),
    ("ero_ls_local", "real", "-", "local LS, 3x3 mean", "", "info"),
    ("ero_twi_local", "real", "ln(m)", "local TWI, 3x3 mean", "", "info"),
    ("ero_spi_class", "text", "-", "SPI class of ero_lnspi_max3x3", "erosion metadata scheme",
     "erosion"),
    ("ero_spi_score", "int", "1-5", "severity score of ero_spi_class", "", "MCDMA"),
    ("ero_worst_score_3x3", "int", "1-5", "worst combined severity score in the 3x3 window", "",
     "erosion"),
    ("ero_composite", "real", "1-5", "hydrodynamic impact score",
     "0.4 SPI + 0.3 RUSLE + 0.3 sediment score (Akello & Omosa 2025; weights in metadata)",
     "MCDMA"),
    ("ero_impact", "text", "-", "impact level", "<= 2.0 Low, <= 3.0 Moderate, <= 4.0 High, "
     "else Severe", "MCDMA"),
    ("ero_mitigation", "text", "-", "indicative mitigation for the impact level", "", "MCDMA"),
]

CATCHMENTS = [
    ("outlet_uid", "text", "-", "link key (see crossings)", "", "link key"),
    ("catchment_id", "text", "-", "catchment identifier (default = outlet_uid)", "",
     "catchment ID"),
    ("outlet_id", "int", "-", "pour-point feature id - NOT stable", "", "none"),
    ("catchment_mode", "text", "-", "'full' = entire upstream area (overlapping); "
     "'local' = non-overlapping local area", "", "provenance"),
    ("area_km2", "real", "km2", "contributing area", "cell count x cell area",
     "area_km2"),
    ("cells", "int", "-", "number of DEM cells in the catchment", "", "QA"),
    ("elev_max_m", "real", "m", "highest elevation", "raw DEM", "provenance"),
    ("elev_min_m", "real", "m", "lowest elevation", "raw DEM", "provenance"),
    ("elev_mean_m", "real", "m", "mean elevation", "raw DEM", "provenance"),
    ("relief_m", "real", "m", "elev_max_m - elev_min_m", "raw DEM", "provenance"),
    ("catch_slope_horn", "real", "m/m", "mean terrain slope over the catchment",
     "Horn 3x3 on the raw DEM, cell mean", "catchment_slope"),
    ("catch_relief_ratio", "real", "m/m", "relief ratio (basin-steepness index)",
     "relief_m / LFP length (Schumm variant: LFP length, not basin length); not a Tc input",
     "catchment relief ratio"),
    ("lfp_length_km", "real", "km", "longest flow path length (planimetric)",
     "diagonal-weighted D8 path", "flow_path_km"),
    ("slope_mean", "real", "m/m", "LEGACY alias of catch_slope_horn (removed in a "
     "later release)", "", "alias"),
    ("slope_relief_ratio", "real", "m/m", "LEGACY alias of catch_relief_ratio "
     "(removed in a later release)", "", "alias"),
    # basin shape and network indices (v0.15, F6) - information only
    ("perimeter_km", "real", "km", "catchment outline length",
     "3x3-smoothed mask contoured at 0.5 (marching squares); not the cell-edge staircase",
     "morphometry"),
    ("form_factor", "real", "-", "form factor A / L^2", "Horton (1932); L = LFP length",
     "morphometry"),
    ("elongation_ratio", "real", "-", "elongation ratio (2/L) sqrt(A/pi)",
     "Schumm (1956); L = LFP length", "morphometry"),
    ("circularity_ratio", "real", "-", "circularity ratio 4 pi A / P^2", "Miller (1953)",
     "morphometry"),
    ("drainage_density", "real", "km/km2", "channel length / area",
     "channel cells at stream_threshold_cells; length along D8 links", "morphometry"),
    ("stream_frequency", "real", "1/km2", "channel links / area",
     "a link runs from a source or confluence to the next confluence or the outlet",
     "morphometry"),
    ("max_strahler", "int", "-", "highest Strahler order in the catchment",
     "stream-order raster if given, else computed at the stream threshold", "morphometry"),
] + _LINK + [
    # soil block (v0.11, WP-C) - empty when no soil dataset was supplied
    ("soil_sand_pct", "real", "%", "topsoil sand", "area-weighted over soil units; "
     "component-weighted within a unit; depth-weighted over soil_depth_cm", "sediment/CN input"),
    ("soil_silt_pct", "real", "%", "topsoil silt", "as soil_sand_pct", "sediment/CN input"),
    ("soil_clay_pct", "real", "%", "topsoil clay", "as soil_sand_pct", "sediment/CN input"),
    ("soil_oc_pct", "real", "%", "topsoil organic carbon", "as soil_sand_pct (SOTWIS TOTC g/kg / 10)",
     "sediment input"),
    ("soil_cfrag_pct", "real", "vol %", "coarse fragments", "as soil_sand_pct", "MUSLE CFRG input"),
    ("soil_bulk_gcm3", "real", "g/cm3", "bulk density", "as soil_sand_pct", "info"),
    ("soil_texture", "text", "-", "USDA texture class of the weighted texture", "USDA triangle",
     "info"),
    ("soil_drain_class", "text", "-", "FAO drainage class with the largest share "
     "(E, S, W, M, I, P, V)", "", "info"),
    ("soil_hsg_proxy", "text", "-", "hydrologic soil group PROXY (A-D) from texture and "
     "drainage - not a measured infiltration class", "see core/soils/usle_k.py", "CN input (proxy)"),
    ("usle_k", "real", "t.ha.h/(ha.MJ.mm)", "USLE/RUSLE K (SI)", "Williams/EPIC per component "
     "(f_csand coefficient 0.0256) x 0.1317, then weighted", "MUSLE K"),
    ("usle_k_dg", "real", "t.ha.h/(ha.MJ.mm)", "alternative K from geometric mean particle "
     "diameter", "Renard et al. 1997", "alternative K"),
    ("usle_cfrg", "real", "-", "coarse-fragment factor", "exp(-0.053 x soil_cfrag_pct)", "MUSLE CFRG"),
    ("soil_coverage_pct", "real", "%", "share of the catchment covered by soil units with data",
     "", "QA"),
    ("soil_dominant_unit", "text", "-", "soil unit covering the largest share", "", "provenance"),
    ("soil_units", "int", "-", "number of soil units with data in the catchment", "", "QA"),
    ("soil_ttr_main", "text", "-", "taxotransfer rule class of the dominant unit's main "
     "profile (SOTWIS TTRmain; confidence indicator)", "", "QA"),
    ("soil_dataset", "text", "-", "soil dataset and version", "", "provenance"),
    ("soil_depth_cm", "text", "cm", "depth interval the soil values describe (default 0-20, D6)",
     "", "provenance"),
    ("soil_hsg", "text", "-", "hydrologic soil group with the largest share (A-D)",
     "direct HSG input (e.g. HYSOGs250m) where given, else the per-cell texture/drainage proxy; "
     "see soil_hsg_source", "CN input"),
    ("hsg_pct_a", "real", "%", "share of the catchment in HSG A", "of cells with a known group",
     "CN input"),
    ("hsg_pct_b", "real", "%", "share in HSG B", "", "CN input"),
    ("hsg_pct_c", "real", "%", "share in HSG C", "", "CN input"),
    ("hsg_pct_d", "real", "%", "share in HSG D, including dual groups (A/D, B/D, C/D)",
     "dual groups counted as D (undrained)", "CN input"),
    ("hsg_pct_dual", "real", "%", "share in dual groups (high runoff unless drained)",
     "part of hsg_pct_d", "CN input"),
    ("soil_hsg_source", "text", "-", "source of the HSG values", "direct dataset or "
     "'proxy: texture and drainage class'", "provenance"),
    ("usle_k_source", "text", "-", "source of usle_k", "direct K input, or Williams/EPIC "
     "from the soil dataset", "provenance"),
] + [
    # erosion block (v0.14, WP-F) - empty when no erosion inputs were supplied
    ("ero_mode", "text", "-", "RUSLE, or LS-only when R, K or C was not supplied",
     "no factor is ever invented", "erosion provenance"),
    ("ero_a_mean_tha", "real", "t/ha/yr", "mean RUSLE soil loss over the catchment",
     "A = R K LS C P", "MUSLE / sediment"),
    ("ero_a_p90_tha", "real", "t/ha/yr", "90th percentile of RUSLE soil loss", "", "erosion"),
    ("ero_ls_mean", "real", "-", "mean LS factor", "Moore & Burch or Desmet & Govers "
     "(ero metadata)", "MUSLE LS"),
    ("ero_ls_p90", "real", "-", "90th percentile LS", "", "erosion"),
    ("ero_lnspi_ch_p90", "real", "ln(m)", "90th percentile ln(SPI) on channel cells",
     "SPI = A_s tan(beta)", "gully potential"),
    ("ero_k_mean", "real", "t ha h/(ha MJ mm)", "area-weighted K", "mean over catchment cells",
     "MUSLE K"),
    ("ero_c_mean", "real", "-", "area-weighted C", "see erosion metadata for proxy flags",
     "MUSLE C"),
    ("ero_p_mean", "real", "-", "area-weighted P", "", "MUSLE P"),
    ("ero_pct_s1", "real", "%", "share of cells with combined severity score 1 (lowest)",
     "D7 matrix on SPI and RUSLE (or LS) scores", "erosion"),
    ("ero_pct_s2", "real", "%", "share with combined score 2", "", "erosion"),
    ("ero_pct_s3", "real", "%", "share with combined score 3", "", "erosion"),
    ("ero_pct_s4", "real", "%", "share with combined score 4", "", "erosion"),
    ("ero_pct_s5", "real", "%", "share with combined score 5 (most severe)", "", "erosion"),
    ("ero_rusle_class", "text", "-", "RUSLE class of the mean soil loss", "erosion metadata "
     "scheme", "erosion"),
    ("ero_rusle_score", "int", "1-5", "severity score of ero_rusle_class", "", "MCDMA"),
    ("ero_gross_t_yr", "real", "t/yr", "gross soil loss", "ero_a_mean_tha x area (ha)",
     "sediment"),
    ("ero_sdr", "real", "-", "sediment delivery ratio", "0.565 A_km2^-0.125, capped at 1 (FAO)",
     "sediment"),
    ("ero_sy_t_yr", "real", "t/yr", "sediment yield to the crossing", "gross x SDR", "sediment"),
    ("ero_sy_m3_yr", "real", "m3/yr", "sediment volume to the crossing",
     "yield x 1000 / bulk density", "sediment"),
    ("ero_bulk_kgm3", "real", "kg/m3", "bulk density used for the volume",
     "soil block when present, else the tool value", "sediment"),
    ("ero_sy_class", "text", "-", "sediment volume impact class", "<1,000 / 5,000 / 15,000 m3/yr",
     "MCDMA"),
    ("ero_sy_score", "int", "1-5", "severity score of ero_sy_class", "", "MCDMA"),
]

FLOWPATHS = [
    ("outlet_uid", "text", "-", "link key (see crossings)", "", "link key"),
    ("flowpath_id", "text", "-", "flow path identifier (default = outlet_uid)", "",
     "flow path ID"),
    ("outlet_id", "int", "-", "pour-point feature id - NOT stable", "", "none"),
    ("lfp_length_km", "real", "km", "longest flow path length", "planimetric, "
     "diagonal-weighted", "flow_path_km"),
    ("lfp_length_m", "real", "m", "longest flow path length", "", "provenance"),
    ("lfp_elev_max_m", "real", "m", "highest elevation on the path", "raw DEM",
     "provenance"),
    ("lfp_elev_min_m", "real", "m", "lowest elevation on the path", "raw DEM",
     "provenance"),
    ("lfp_drop_m", "real", "m", "lfp_elev_max_m - lfp_elev_min_m",
     "highest - lowest raw-DEM elevation on the path; equals headwater - outlet on a "
     "monotonic path (see lfp_nonmonotonic)", "provenance"),
    ("lfp_slope", "real", "m/m", "whole-path slope", "lfp_drop_m / lfp_length_m",
     "flow_path_slope"),
    ("lfp_slope_1085", "real", "m/m", "10-85 slope along the LFP",
     "(z85 - z10) / (L85 - L10), L measured from the OUTLET; elevations "
     "interpolated linearly along the path (QEHT <= 0.8.3 measured from the divide)",
     "flow_path_slope_10_85"),
    ("lfp_L10_m", "real", "m", "distance of the 10% point from the outlet",
     "0.10 x lfp_length_m", "hand check"),
    ("lfp_L85_m", "real", "m", "distance of the 85% point from the outlet",
     "0.85 x lfp_length_m", "hand check"),
    ("lfp_z10_m", "real", "m", "elevation at the 10% point", "raw DEM, interpolated",
     "hand check"),
    ("lfp_z85_m", "real", "m", "elevation at the 85% point", "raw DEM, interpolated",
     "hand check"),
    ("lfp_z_head_m", "real", "m", "raw-DEM elevation of the first (headwater) path cell", "",
     "hand check"),
    ("lfp_z_outlet_m", "real", "m", "raw-DEM elevation of the last (outlet) path cell", "",
     "hand check"),
    ("lfp_nonmonotonic", "int", "0/1", "1 when lfp_drop_m exceeds lfp_z_head_m - lfp_z_outlet_m "
     "by more than 0.5 m (spike or pit on the path; lfp_slope slightly overstated)", "", "QA"),
    ("area_km2", "real", "km2", "area of the catchment this path belongs to", "",
     "provenance"),
    # overland / channel split (v0.15, F2)
    ("lfp_overland_m", "real", "m", "divide -> channel head length",
     "channel head = first path cell with accumulation >= lfp_threshold_km2",
     "Kerby overland length"),
    ("lfp_overland_slope", "real", "m/m", "slope of the overland part",
     "(z head - z channel head) / lfp_overland_m, raw DEM", "Kerby overland slope"),
    ("lfp_channel_m", "real", "m", "channel head -> outlet length",
     "lfp_overland_m + lfp_channel_m = lfp_length_m", "Kirpich / Kerby channel length"),
    ("lfp_channel_slope", "real", "m/m", "slope of the channel part",
     "(z channel head - z outlet) / lfp_channel_m", "Kirpich / Kerby channel slope"),
    ("lfp_channel_slope_1085", "real", "m/m", "10-85 slope of the channel part",
     "outlet-referenced, as lfp_slope_1085", "channel slope (10-85)"),
    ("lfp_sheet_m", "real", "m", "sheet-flow part of the overland length",
     "min(lfp_overland_m, sheet cap; default 100 m, TR-55 practice)", "TR-55 sheet flow"),
    ("lfp_shallow_m", "real", "m", "shallow concentrated part of the overland length",
     "lfp_overland_m - lfp_sheet_m", "TR-55 shallow flow"),
    ("lfp_threshold_km2", "real", "km2", "channel threshold used for the split",
     "stream threshold x cell area", "provenance"),
    ("lfp_no_channel", "int", "0/1", "1 when no path cell reaches the threshold "
     "(whole path overland)", "", "QA"),
] + _LINK

LAYERS = {
    "crossings": ("POINT", CROSSINGS,
                  "Crossing / outlet points (snapped location)"),
    "catchments": ("MULTIPOLYGON", CATCHMENTS,
                   "One catchment per crossing, with morphometry"),
    "flowpaths": ("LINESTRING", FLOWPATHS,
                  "Longest flow path per crossing, digitised divide -> outlet"),
}

# Optional layers (additive; written only when their input was given). Listed
# in qeht_field_dictionary so the package still documents itself.
OPTIONAL_LAYERS = {
    "alignment_profile": ("POINT", [
        ("align_name", "text", "-", "alignment identifier", "", "info"),
        ("chainage_m", "real", "m", "station chainage", "start chainage + distance along the "
         "alignment, every alignment_step_m", "crossing chainage"),
        ("x", "real", "m", "station easting (package CRS)", "", "info"),
        ("y", "real", "m", "station northing (package CRS)", "", "info"),
        ("z_dem_m", "real", "m", "raw DEM ground at the station", "bilinear between cell centres",
         "ground profile"),
        ("z_fill_m", "real", "m", "filled DEM at the station (empty without a filled DEM)",
         "bilinear", "ground profile"),
        ("pond_depth_m", "real", "m", "depth the fill raised the ground (ponding)",
         "max(z_fill_m - z_dem_m, 0)", "ground profile"),
        ("acc_km2", "real", "km2", "contributing area at the station",
         "largest (accumulation + 1) x cell area within 1 cell", "ground profile"),
        ("stream", "int", "0/1", "1 at the station nearest a D8 stream crossing",
         "exact stream link x alignment intersection at the stream threshold", "ground profile"),
        ("strahler", "int", "-", "Strahler order of that stream (if an order raster was given)",
         "", "info"),
        ("slope_long_pct", "real", "%", "longitudinal ground slope, positive = rising chainage",
         "centred difference of z_dem_m (one-sided at part ends)", "ground profile"),
        ("alignment_part", "int", "-", "part of a multi-part alignment", "", "info"),
    ], "Ground profile along the road alignment (stations)"),
}

METADATA_KEYS = [
    # key, meaning
    ("schema_version", "exchange schema (qeht-heas-<major>)"),
    ("qeht_version", "QEHT plugin version that wrote the file"),
    ("run_utc", "UTC time of the run (ISO 8601)"),
    ("crs_epsg", "EPSG code of the (projected, metric) CRS of every layer"),
    ("crs_name", "CRS description"),
    ("cell_size_m", "DEM cell size (x, y)"),
    ("dem_path", "flow-direction / conditioned DEM used for routing"),
    ("raw_dem_path", "DEM used for reported elevations and slopes"),
    ("dem_sha256", "SHA-256 of the raw DEM file (or size+mtime above the hash limit)"),
    ("dem_source", "user description of the DEM (e.g. ALOS 30 m, LiDAR)"),
    ("conditioning", "how the DEM was conditioned (user-declared)"),
    ("flat_method", "flat resolution used for flow direction (user-declared)"),
    ("tie_rule", "D8 tie-break rule"),
    ("stream_threshold_cells", "snap-to-stream threshold in cells"),
    ("stream_threshold_km2", "snap-to-stream threshold as area"),
    ("snap_radius_cells", "pour-point snap radius in cells"),
    ("snap_strategy", "nearest_stream | max_accumulation | none"),
    ("catchment_mode", "full (overlapping) | local (non-overlapping)"),
    ("id_scheme", "sequential | attribute"),
    ("id_prefix", "prefix used for outlet_uid"),
    ("id_order", "numbering order for sequential ids"),
    ("id_attribute", "pour-point attribute used for ids (attribute scheme)"),
    ("lfp_1085_reference", "outlet (v0.9+) | divide (<= 0.8.3)"),
    ("crossing_source", "pour points | crossing candidates | relinked package"),
    ("chainage_start_m", "start chainage of the road alignment, if one was used"),
    ("relinked_from", "package this one was renumbered/relinked from (D3)"),
    ("soil_dataset", "soil dataset used for the soil block (empty = none)"),
    ("soil_depth_cm", "topsoil depth interval of the soil block"),
    ("soil_hsg_source", "source of the hydrologic soil groups (direct dataset or texture proxy)"),
    ("soil_k_source", "source of K in the soil block (direct K or Williams/EPIC)"),
    ("conditioning_burn", "summary of the breach log when one was supplied (layer burn_log)"),
    ("erosion_json", "erosion inputs: mode, factor provenance (source, proxy flags), LS method, "
     "class schemes, MCDMA weights, SDR model (empty = no erosion block)"),
    ("sheet_cap_m", "sheet-flow cap used for lfp_sheet_m (m)"),
    ("alignment_source", "road alignment layer used for chainage and the ground profile"),
    ("alignment_step_m", "station spacing of alignment_profile (m)"),
    ("n_crossings", "number of crossings written"),
    ("parameters_json", "full parameter dictionary of the run"),
]

FIELD_DICTIONARY_COLUMNS = ("layer", "field", "type", "unit", "meaning",
                            "method", "heas_target", "downstream_use", "heas_field")

# Plain-language "downstream use" for each heas_target term (v0.15, §0 of the
# feature recommendations): the package is for spreadsheets, reports or any
# design software, not only HEAS. heas_target stays as it was (HEAS reads it;
# the table is additive under qeht-heas-1); heas_field names the HEAS input
# only where the target is one.
DOWNSTREAM_USE = {
    "link key": "joins crossings, catchments and flow paths",
    "link registry": "how the feature was linked to its crossing",
    "crossing ID": "crossing label in schedules and reports",
    "catchment ID": "catchment label in schedules and reports",
    "flow path ID": "flow path label in schedules and reports",
    "crossing chainage": "crossing position along the road",
    "none": "not for linking (unstable feature id)",
    "provenance": "record of how the value was produced",
    "QA": "quality check",
    "info": "information for the report",
    "hand check": "lets a reviewer check the value by hand",
    "alias": "old name kept for one release",
    "area_km2": "catchment area for any design flood method",
    "catchment_slope": "mean catchment slope (runoff coefficient / CN tables)",
    "catchment relief ratio": "basin-steepness index (morphometry; not a Tc input)",
    "flow_path_km": "time of concentration: flow path length",
    "flow_path_slope": "time of concentration: flow path slope (drop / length)",
    "flow_path_slope_10_85": "time of concentration: 10-85 flow path slope",
    "sediment/CN input": "soil texture for CN and sediment calculations",
    "sediment input": "soil property for sediment calculations",
    "CN input (proxy)": "hydrologic soil group for curve numbers (texture proxy)",
    "CN input": "hydrologic soil group for curve numbers",
    "MUSLE K": "soil erodibility for RUSLE / MUSLE",
    "MUSLE LS": "slope length-steepness factor for RUSLE / MUSLE",
    "MUSLE C": "cover factor for RUSLE / MUSLE",
    "MUSLE P": "support practice factor for RUSLE / MUSLE",
    "MUSLE CFRG": "coarse-fragment factor for MUSLE",
    "MUSLE CFRG input": "coarse fragments for MUSLE",
    "MUSLE / sediment": "soil loss for sediment calculations",
    "alternative K": "alternative soil erodibility",
    "gully potential": "gully / scour potential at the crossing",
    "erosion": "erosion severity",
    "erosion provenance": "how the erosion values were produced",
    "sediment": "sediment yield to the crossing",
    "MCDMA": "multi-criteria erosion impact score",
    "ground profile": "ground profile along the road (provisional levels, DEM check)",
    "morphometry": "basin shape / network index for the report (information only)",
    "Kerby overland length": "time of concentration: overland (Kerby) length",
    "Kerby overland slope": "time of concentration: overland (Kerby) slope",
    "Kirpich / Kerby channel length": "time of concentration: channel length",
    "Kirpich / Kerby channel slope": "time of concentration: channel slope",
    "channel slope (10-85)": "time of concentration: 10-85 channel slope",
    "TR-55 sheet flow": "segmental Tc: sheet-flow length",
    "TR-55 shallow flow": "segmental Tc: shallow concentrated flow length",
}

# heas_target terms that are HEAS input names (written to heas_field)
HEAS_FIELDS = {"area_km2", "catchment_slope", "flow_path_km", "flow_path_slope",
               "flow_path_slope_10_85"}


def downstream_use(target):
    return DOWNSTREAM_USE.get(target, target)


def field_names(layer):
    return [f[0] for f in LAYERS[layer][1]]


def dictionary_rows():
    """Rows for the qeht_field_dictionary table, in contract order."""
    rows = []
    for layer, (_, fields, _) in LAYERS.items():
        for name, ftype, unit, meaning, method, target in fields:
            rows.append((layer, name, ftype, unit, meaning, method, target,
                         downstream_use(target), target if target in HEAS_FIELDS else ""))
    for layer, (_, fields, _) in OPTIONAL_LAYERS.items():
        for name, ftype, unit, meaning, method, target in fields:
            rows.append((layer, name, ftype, unit, meaning, method, target,
                         downstream_use(target), ""))
    for key, meaning in METADATA_KEYS:
        rows.append(("qeht_run_metadata", key, "text", "-", meaning, "", "provenance",
                     downstream_use("provenance"), ""))
    return rows
