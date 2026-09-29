# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Per-catchment soil block (area-weighted over soil units).

Every value is weighted by the number of catchment cells in each soil unit.
K is averaged as K (each unit's K was already computed per component), not
recomputed from averaged texture. `soil_coverage_pct` is the share of the
catchment covered by units that have data; the other values describe only
the covered part.

No QGIS imports, no GDAL.
"""

import math

import numpy as np

NAN = float("nan")

SOIL_FIELDS = [
    ("soil_sand_pct", "real"), ("soil_silt_pct", "real"), ("soil_clay_pct", "real"),
    ("soil_oc_pct", "real"), ("soil_cfrag_pct", "real"), ("soil_bulk_gcm3", "real"),
    ("soil_texture", "text"), ("soil_drain_class", "text"), ("soil_hsg_proxy", "text"),
    ("usle_k", "real"), ("usle_k_dg", "real"), ("usle_cfrg", "real"),
    ("soil_coverage_pct", "real"), ("soil_dominant_unit", "text"),
    ("soil_units", "int"), ("soil_ttr_main", "text"), ("soil_dataset", "text"),
    ("soil_depth_cm", "text"),
]

_MAP = {"soil_sand_pct": "sand", "soil_silt_pct": "silt", "soil_clay_pct": "clay",
        "soil_oc_pct": "oc_pct", "soil_cfrag_pct": "cfrag", "soil_bulk_gcm3": "bulk",
        "usle_k": "k_si", "usle_k_dg": "k_dg"}


def soil_block(mask, unit_grid, units, index_to_unit, info=None):
    """Soil fields for one catchment.

    mask : bool catchment mask; unit_grid : int grid of unit indices (0 =
    no soil polygon); index_to_unit : {index: unit_id}; units : {unit_id:
    UnitSoil}; info : dataset provenance (soil_dataset, soil_depth_cm).
    """
    from .usle_k import usda_texture, hsg_proxy, cfrg
    out = {f: None for f, _ in SOIL_FIELDS}
    info = info or {}
    out["soil_dataset"] = info.get("soil_dataset")
    out["soil_depth_cm"] = info.get("soil_depth_cm")
    n = int(mask.sum())
    if n == 0:
        return out
    idx, counts = np.unique(unit_grid[mask], return_counts=True)
    weights = {}
    for i, cnt in zip(idx.tolist(), counts.tolist()):
        u = units.get(index_to_unit.get(i))
        if i == 0 or u is None or u.k_si is None or not math.isfinite(u.k_si):
            continue
        weights[u.unit_id] = weights.get(u.unit_id, 0) + cnt
    covered = sum(weights.values())
    out["soil_coverage_pct"] = 100.0 * covered / n
    out["soil_units"] = len(weights)
    if covered == 0:
        return out
    for f, attr in _MAP.items():
        num = den = 0.0
        for uid, w in weights.items():
            v = getattr(units[uid], attr)
            if v is not None and math.isfinite(v):
                num += w * v; den += w
        out[f] = num / den if den > 0 else None
    dom = max(weights, key=weights.get)
    out["soil_dominant_unit"] = str(dom)
    out["soil_ttr_main"] = units[dom].ttr_main
    dr = {}
    for uid, w in weights.items():
        d = units[uid].drain
        if d:
            dr[d] = dr.get(d, 0) + w
    out["soil_drain_class"] = max(dr, key=dr.get) if dr else None
    if None not in (out["soil_sand_pct"], out["soil_silt_pct"], out["soil_clay_pct"]):
        out["soil_texture"] = usda_texture(out["soil_sand_pct"], out["soil_silt_pct"],
                                           out["soil_clay_pct"])
        out["soil_hsg_proxy"] = hsg_proxy(out["soil_texture"], out["soil_drain_class"])
    if out["soil_cfrag_pct"] is not None:
        out["usle_cfrg"] = cfrg(out["soil_cfrag_pct"])
    return out
