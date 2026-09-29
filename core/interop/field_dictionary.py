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
    ("chainage_m", "real", "m", "position along the road alignment (empty until the "
     "crossing-candidate tool supplies it)", "linear referencing", "crossing chainage"),
] + _LINK

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
    ("catch_relief_ratio", "real", "m/m", "relief ratio", "relief_m / LFP length",
     "catchment relief ratio"),
    ("lfp_length_km", "real", "km", "longest flow path length (planimetric)",
     "diagonal-weighted D8 path", "flow_path_km"),
    ("slope_mean", "real", "m/m", "LEGACY alias of catch_slope_horn (removed in a "
     "later release)", "", "alias"),
    ("slope_relief_ratio", "real", "m/m", "LEGACY alias of catch_relief_ratio "
     "(removed in a later release)", "", "alias"),
] + _LINK

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
    ("lfp_drop_m", "real", "m", "lfp_elev_max_m - lfp_elev_min_m", "", "provenance"),
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
    ("area_km2", "real", "km2", "area of the catchment this path belongs to", "",
     "provenance"),
] + _LINK

LAYERS = {
    "crossings": ("POINT", CROSSINGS,
                  "Crossing / outlet points (snapped location)"),
    "catchments": ("MULTIPOLYGON", CATCHMENTS,
                   "One catchment per crossing, with morphometry"),
    "flowpaths": ("LINESTRING", FLOWPATHS,
                  "Longest flow path per crossing, digitised divide -> outlet"),
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
    ("n_crossings", "number of crossings written"),
    ("parameters_json", "full parameter dictionary of the run"),
]

FIELD_DICTIONARY_COLUMNS = ("layer", "field", "type", "unit", "meaning",
                            "method", "heas_target")


def field_names(layer):
    return [f[0] for f in LAYERS[layer][1]]


def dictionary_rows():
    """Rows for the qeht_field_dictionary table, in contract order."""
    rows = []
    for layer, (_, fields, _) in LAYERS.items():
        for name, ftype, unit, meaning, method, target in fields:
            rows.append((layer, name, ftype, unit, meaning, method, target))
    for key, meaning in METADATA_KEYS:
        rows.append(("qeht_run_metadata", key, "text", "-", meaning, "", "provenance"))
    return rows
