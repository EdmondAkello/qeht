# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Review register (1.2): every finding of a run in one table.

The run already reports its findings, but in different places: DEM
warnings in the log, flags in the characteristics table, truncated
catchments in the exchange issues, deep fills in the depression layer. The
register puts them in one list a reviewer can work through, with an ID,
severity, chainage, the evidence and a suggested check, plus empty
`status` / `response` columns for the review itself.

Built only from the package and what the pipeline passes in `extra`;
nothing is recalculated. Severity:

    high    the result at that crossing may be wrong (truncated catchment,
            crossing position depends on flat routing, deep fill against
            the road far from any crossing)
    medium  the result depends on a modelling choice or on DEM error
            (flat-method area change, uncertainty switching, mapped
            drainage disagrees, DEM warnings, deep fill near the road)
    low     worth noting (outlet snapped far, Tc outside the calibration
            range of most methods, missing USLE K, non-monotonic flow path)
    info    run notes (proposed crossings, other warnings)

No QGIS imports, no GDAL.
"""

import csv
import math

from ..interop import gpkg

SEVERITY_ORDER = {"high": 0, "medium": 1, "low": 2, "info": 3}

REGISTER_FIELDS = [
    ("finding_id", "str"), ("severity", "str"), ("category", "str"), ("outlet_uid", "str"),
    ("chainage_m", "float"), ("finding", "str"), ("evidence", "str"), ("check", "str"),
    ("status", "str"), ("response", "str"),
]

# pipeline warnings already listed per crossing above
COVERED = ("between the two flat methods", "under the other flat method",
           "proposed crossing(s) from the coverage check")

TC_METHODS = ("kirpich", "kerby_kirpich", "scs_lag", "tr55", "bransby_williams")


def _num(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _f(v, nd=1):
    f = _num(v)
    return "–" if f is None else f"{f:,.{nd}f}"


def build_register(package_path, extra=None, depressions=None, near_m=100.0):
    """-> list of register rows (dicts with REGISTER_FIELDS keys), sorted by
    severity then chainage, IDs R-001...

    extra: as for the run report - warnings [str], dem_audit [str],
    flat_check {uid: {...}}. depressions: rows from
    core.review.depressions.depression_table (flagged rows are listed).
    """
    extra = dict(extra or {})
    tables = {t for t, _ in gpkg.list_tables(package_path)}
    md = {r["key"]: r["value"] for r in gpkg.read_table(package_path, "qeht_run_metadata")}
    crossings = gpkg.read_table(package_path, "crossings", with_geometry=False)
    catch = {r["outlet_uid"]: r for r in gpkg.read_table(package_path, "catchments",
                                                          with_geometry=False)}
    paths = {r["outlet_uid"]: r for r in gpkg.read_table(package_path, "flowpaths",
                                                         with_geometry=False)} \
        if "flowpaths" in tables else {}
    flat = extra.get("flat_check") or {}
    cell = _num(md.get("cell_size_m")) or 0.0
    out = []

    def add(sev, cat, uid, ch, finding, evidence, check):
        out.append({"severity": sev, "category": cat, "outlet_uid": uid or "",
                    "chainage_m": _num(ch), "finding": finding, "evidence": evidence,
                    "check": check, "status": "open", "response": ""})

    for c in crossings:
        uid, ch = c.get("outlet_uid"), c.get("chainage_m")
        k = catch.get(uid, {})
        fp = paths.get(uid, {})
        if k.get("truncated"):
            add("high", "catchment", uid, ch, "Catchment is cut off by the DEM edge or NoData",
                f"{int(k.get('truncated_cells') or 0):,} catchment cells on the edge; area "
                f"{_f(k.get('area_km2'), 3)} km²",
                "Extend the DEM upstream (or fill the NoData) and re-run; the area is a lower bound.")
        if k.get("clip_edge"):
            add("high", "catchment", uid, ch, "Catchment touches the edge of the automatic clip",
                f"area {_f(k.get('area_km2'), 3)} km²",
                "Re-run with a larger clip buffer or without clipping.")
        v = flat.get(uid) or {}
        if v.get("flat_moved"):
            add("high", "flat routing", uid, ch,
                "The channel at this crossing runs elsewhere under the other flat method",
                f"area Barnes {_f(v.get('area_barnes_km2'), 3)} km², toward-lower "
                f"{_f(v.get('area_toward_km2'), 3)} km²; flat share of catchment "
                f"{_f(v.get('flat_catch_share_pct'))} %",
                "Confirm the channel position on mapped drainage, imagery or on site.")
        elif v.get("flat_sensitive"):
            add("medium", "flat routing", uid, ch,
                "Contributing area depends on how flats are routed",
                f"area changes {_f(v.get('flat_sensitivity_pct'))} % between methods "
                f"(Barnes {_f(v.get('area_barnes_km2'), 3)}, toward-lower "
                f"{_f(v.get('area_toward_km2'), 3)} km²)",
                "Verify the catchment boundary against mapped drainage or on site.")
        sw, lost = _num(c.get("unc_switch_pct")), _num(c.get("unc_lost_pct"))
        if (sw or 0) > 0 or (lost or 0) > 0:
            add("medium", "DEM uncertainty", uid, ch,
                "Area is unstable under DEM vertical error",
                f"area moved >25 % in {_f(sw)} % of realisations, stream lost in {_f(lost)} %"
                + (f"; P10–P90 {_f(c.get('unc_area_p10'), 3)}–{_f(c.get('unc_area_p90'), 3)} km²"
                   if _num(c.get("unc_area_p10")) is not None else ""),
                "Treat the area as a range; check the divide on a better DEM or on site.")
        if c.get("map_agrees") == 0:
            add("medium", "mapped drainage", uid, ch,
                "Outlet is not on a mapped waterway",
                f"nearest mapped waterway {_f(c.get('map_river_dist_m'), 0)} m"
                + (f" ({c.get('map_river_name')})" if c.get("map_river_name") else "")
                + (f"; {c.get('map_note')}" if c.get("map_note") else ""),
                "Check whether the DEM channel or the map is right at this crossing.")
        snap = _num(c.get("snap_dist_m"))
        if snap is not None and cell and snap > 2.5 * cell:
            add("low", "outlet", uid, ch, "Outlet snapped far from the input point",
                f"{snap:,.0f} m ({snap / cell:.1f} cells)",
                "Check that the snapped outlet is on the intended channel.")
        have = [m for m in TC_METHODS if _num(c.get(f"tc_{m}_min")) is not None]
        outside = [m for m in have if str(c.get(f"tc_{m}_flag") or "").lower().startswith("outside")]
        if have and 2 * len(outside) > len(have):     # one method out of range is routine
            add("low", "time of concentration", uid, ch,
                "Catchment is outside the calibration range of most Tc methods",
                "; ".join(f"{m.replace('_', ' ')}: {c.get(f'tc_{m}_flag')}" for m in outside),
                "Choose Tc with care; compare with a velocity-based estimate.")
        if fp.get("lfp_nonmonotonic"):
            add("low", "flow path", uid, ch, "Longest flow path rises along its length",
                f"drop along path {_f(fp.get('lfp_drop_m'))} m vs head–outlet "
                f"{_f((_num(fp.get('lfp_z_head_m')) or 0) - (_num(fp.get('lfp_z_outlet_m')) or 0))} m",
                "Use the 10–85 slope; check the path over filled ground.")
        if md.get("soil_dataset") and k and _num(k.get("usle_k")) is None:
            add("low", "soil", uid, ch, "The soil source gave no USLE K for this catchment",
                f"soil coverage {_f(k.get('soil_coverage_pct'))} %",
                "Supply K from another source before using the erosion results.")
        if c.get("status") == "proposed":
            add("info", "coverage", uid, ch, "Proposed crossing from the coverage check",
                c.get("proposed_reason") or "", "Adopt or delete it, then Renumber and relink.")

    for d in depressions or []:
        if not d.get("flag"):
            continue
        far = d.get("nearest_uid") is None or (d.get("nearest_ch_dist_m") or 0) > near_m
        sev = "high" if d["flag"] == "deep_at_road" and far else "medium"
        add(sev, "fill", d.get("nearest_uid") if not far else "", d.get("chainage_m"),
            "Deep fill against the road" if d["flag"] == "deep_at_road" else "Deep fill near the road",
            f"depression {d['dep_id']}: max {d['max_depth_m']:.2f} m, {d['area_m2']:,.0f} m², "
            f"{d['volume_m3']:,.0f} m³, {_f(d.get('dist_road_m'), 0)} m from the centreline"
            + ("" if d.get("nearest_uid") is None else
               f"; nearest crossing {d['nearest_uid']} at {_f(d.get('nearest_ch_dist_m'), 0)} m"),
            "Check for a culvert the DEM does not show or a missing one; burn it in if it exists."
            if far else "Check that the nearby crossing drains this depression.")

    for w in extra.get("dem_audit") or []:
        add("medium", "DEM", "", None, "DEM check", w, "Read the DEM notes in the documentation.")
    for w in extra.get("warnings") or []:
        if w.startswith(("DEM check:", "DEM QA:")) or any(t in w for t in COVERED):
            continue
        add("info", "run", "", None, "Run warning", w, "")

    out.sort(key=lambda r: (SEVERITY_ORDER[r["severity"]],
                            r["chainage_m"] if r["chainage_m"] is not None else math.inf,
                            r["outlet_uid"]))
    for i, r in enumerate(out, start=1):
        r["finding_id"] = f"R-{i:03d}"
    return [{k: r[k] for k, _ in REGISTER_FIELDS} for r in out]


def register_counts(rows):
    n = {s: 0 for s in SEVERITY_ORDER}
    for r in rows:
        n[r["severity"]] += 1
    return n


def write_register_csv(rows, path):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=[k for k, _ in REGISTER_FIELDS])
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else
                            (f"{r[k]:.1f}" if k == "chainage_m" else r[k]))
                        for k, _ in REGISTER_FIELDS})
    return path
