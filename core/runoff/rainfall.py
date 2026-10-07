# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Rainfall zone and mean annual rainfall per catchment (F3, v0.17).

Zones: a polygon layer with a zone name (e.g. the design manual's rainfall
zones, or hydro-climatic zones), rasterised on the DEM grid by the
cell-centre rule, exactly like the soil map. Per catchment: the share of
every zone among the cells covered by a zone, the dominant zone and its
share. Tie rule: when two zones have the same number of cells, the
dominant zone is the one whose name sorts first (case-insensitive), so the
result never depends on the order of the polygons in the layer.

Mean annual rainfall (MAP): a raster in mm/yr (e.g. a CHIRPS or WorldClim
climatology) resampled bilinearly to the DEM grid; the catchment value is
the mean over the cells with data.

RUSLE R from rainfall (optional, an ESTIMATE): per cell from the MAP with a
published annual R-P relation, then the catchment mean (R is not linear in
P, so the mean of R is not R of the mean P). Units MJ mm/(ha h yr).

    renard_freimund_1994   Renard & Freimund (1994), J. Hydrol. 157:287-306
        R = 0.0483 P^1.61                      P <= 850 mm
        R = 587.8 - 1.219 P + 0.004105 P^2     P >  850 mm
    lo_1985                Lo, El-Swaify, Dangler & Shinshiro (1985),
                           Hawaii: R = 38.46 + 3.48 P

Neither relation was derived in East Africa. A measured or modelled
erosivity raster (e.g. GloREDa) is preferable; the relation is a fallback
and every output says which was used.

No QGIS imports, no GDAL.
"""

import json

import numpy as np

RAIN_FIELDS = ["rain_zone", "rain_zone_pct", "rain_zones_json", "rain_zone_coverage_pct",
               "map_mm", "map_coverage_pct", "map_dataset", "rusle_r", "rusle_r_method"]
RAIN_TEXT = ("rain_zone", "rain_zones_json", "map_dataset", "rusle_r_method")


def _renard_freimund(p):
    p = np.asarray(p, dtype=np.float64)
    return np.where(p <= 850.0, 0.0483 * np.power(np.maximum(p, 0.0), 1.61),
                    587.8 - 1.219 * p + 0.004105 * p * p)


def _lo(p):
    return 38.46 + 3.48 * np.asarray(p, dtype=np.float64)


R_RELATIONS = {
    "renard_freimund_1994": (
        "Renard & Freimund (1994): R = 0.0483 P^1.61 (P <= 850 mm), "
        "587.8 - 1.219 P + 0.004105 P^2 (P > 850 mm)", _renard_freimund),
    "lo_1985": ("Lo et al. (1985): R = 38.46 + 3.48 P", _lo),
}
R_RELATION_KEYS = ("none", "renard_freimund_1994", "lo_1985")

MAP_MIN_MM, MAP_MAX_MM = 20.0, 13000.0


def r_from_map(map_mm, relation):
    """Per-cell R (MJ mm/(ha h yr)) from mean annual rainfall in mm; NaN stays NaN."""
    if relation not in R_RELATIONS:
        raise ValueError(f"unknown R-P relation {relation!r}; one of {sorted(R_RELATIONS)}")
    p = np.asarray(map_mm, dtype=np.float64)
    out = R_RELATIONS[relation][1](p)
    return np.where(np.isfinite(p), out, np.nan)


def r_method_text(relation, map_dataset):
    if relation not in R_RELATIONS:
        return ""
    return f"ESTIMATE from mean annual rainfall ({map_dataset}): {R_RELATIONS[relation][0]}"


def check_map(grid):
    """Unit check of a mean-annual-rainfall raster. Returns warnings; raises
    ValueError when the values cannot be annual totals in mm."""
    g = np.asarray(grid, dtype=np.float64)
    v = g[np.isfinite(g)]
    if v.size == 0:
        raise ValueError("Mean annual rainfall: the raster has no data over the DEM.")
    lo, hi = float(v.min()), float(np.percentile(v, 99.5))
    if lo < 0:
        raise ValueError(f"Mean annual rainfall: negative values ({lo:g}).")
    if hi < MAP_MIN_MM:
        raise ValueError(f"Mean annual rainfall: values up to {hi:g} - this looks like mm/day "
                         "or a monthly mean, not mm/yr. Give the annual total in mm.")
    if hi > MAP_MAX_MM:
        raise ValueError(f"Mean annual rainfall: values up to {hi:g} mm/yr exceed the world "
                         "record - check the units (0.1 mm?) or the NoData value.")
    warnings = []
    if hi < 100:
        warnings.append(f"Mean annual rainfall is at most {hi:g} mm/yr - check it is an annual "
                        "total, not a monthly mean.")
    return warnings


def dominant_zone(counts):
    """counts {name: cells} -> (name, share 0-100) with the documented tie
    rule (largest count; ties -> name that sorts first, case-insensitive)."""
    total = sum(counts.values())
    if total <= 0:
        return None, None
    best = sorted(counts.items(), key=lambda kv: (-kv[1], str(kv[0]).casefold(), str(kv[0])))[0]
    return best[0], 100.0 * best[1] / total


def rain_block(mask, zone_grid=None, zone_names=None, map_grid=None, map_dataset="",
               r_grid=None, r_method=""):
    """RAIN_FIELDS for one catchment mask.

    zone_grid: int grid, 0 = no zone, k = index into zone_names {k: name}.
    """
    out = {k: None for k in RAIN_FIELDS}
    mask = np.asarray(mask, bool)
    n = int(mask.sum())
    if n == 0:
        return out
    if zone_grid is not None:
        z = np.asarray(zone_grid)[mask]
        z = z[z > 0]
        out["rain_zone_coverage_pct"] = float(100.0 * z.size / n)
        if z.size:
            idx, cnt = np.unique(z, return_counts=True)
            counts = {}
            for i, c in zip(idx.tolist(), cnt.tolist()):
                name = (zone_names or {}).get(int(i), str(int(i)))
                counts[name] = counts.get(name, 0) + int(c)
            name, pct = dominant_zone(counts)
            out["rain_zone"] = name
            out["rain_zone_pct"] = pct
            total = float(sum(counts.values()))
            out["rain_zones_json"] = json.dumps(
                {k: round(100.0 * v / total, 3) for k, v in
                 sorted(counts.items(), key=lambda kv: (-kv[1], str(kv[0]).casefold()))})
    if map_grid is not None:
        v = np.asarray(map_grid, dtype=np.float64)[mask]
        ok = np.isfinite(v)
        out["map_coverage_pct"] = float(100.0 * ok.sum() / n)
        if ok.any():
            out["map_mm"] = float(v[ok].mean())
        out["map_dataset"] = map_dataset
    if r_grid is not None:
        v = np.asarray(r_grid, dtype=np.float64)[mask]
        ok = np.isfinite(v)
        if ok.any():
            out["rusle_r"] = float(v[ok].mean())
        out["rusle_r_method"] = r_method
    return out


class RainfallInputs(object):
    """Rainfall zones, MAP and (optionally) R from MAP on the DEM grid + provenance."""

    def __init__(self, zone_grid=None, zone_names=None, zone_source="", map_grid=None,
                 map_dataset="", r_relation="none", zone_threshold_pct=80.0):
        self.zone_grid = None if zone_grid is None else np.asarray(zone_grid)
        self.zone_names = dict(zone_names or {})
        self.zone_source = zone_source
        self.map_grid = None if map_grid is None else np.asarray(map_grid, dtype=np.float64)
        self.map_dataset = map_dataset
        self.r_relation = r_relation if r_relation in R_RELATIONS else "none"
        self.zone_threshold_pct = float(zone_threshold_pct)
        self.r_grid = None
        self.r_method = ""
        self.warnings = []
        if self.map_grid is not None:
            self.warnings += check_map(self.map_grid)
            if self.r_relation != "none":
                self.r_grid = r_from_map(self.map_grid, self.r_relation)
                self.r_method = r_method_text(self.r_relation, map_dataset)
        elif self.r_relation != "none":
            self.warnings.append("An R-P relation was chosen but no rainfall raster was given; "
                                 "rusle_r left empty.")
        self.meta = {"zone_source": zone_source,
                     "zones": sorted(set(self.zone_names.values()), key=str.casefold),
                     "zone_tie_rule": "largest share; ties -> zone name sorting first",
                     "zone_qa_threshold_pct": self.zone_threshold_pct,
                     "map_dataset": map_dataset,
                     "map_resampling": "bilinear to the DEM grid" if map_grid is not None else "",
                     "r_relation": self.r_relation,
                     "r_relation_text": R_RELATIONS[self.r_relation][0]
                     if self.r_relation in R_RELATIONS else "",
                     "r_is_estimate": self.r_relation != "none"}

    @property
    def present(self):
        return self.zone_grid is not None or self.map_grid is not None

    def block(self, mask):
        return rain_block(mask, self.zone_grid, self.zone_names, self.map_grid,
                          self.map_dataset, self.r_grid, self.r_method)

    def qa_issue(self, uid, block):
        """A warning text when the dominant zone covers less than the threshold."""
        pct = block.get("rain_zone_pct")
        if pct is not None and pct < self.zone_threshold_pct:
            return (f"{uid}: rainfall zone '{block['rain_zone']}' covers only {pct:.0f} % of the "
                    f"catchment (zones: {block['rain_zones_json']}) - choose the design zone "
                    "by hand.")
        return None

    def meta_json(self):
        return json.dumps(self.meta, sort_keys=True)

