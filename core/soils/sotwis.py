# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Soil datasets -> per-soil-unit topsoil properties (WP-C).

Two loaders, one output shape: {unit_id: UnitSoil}.

SOTER / SOTWIS (preset, Kenya SOTWIS v1.0)
------------------------------------------
Stdlib `sqlite3` on the SOTWIS SQLite database (never the .mdb):
    SOTERunitComposition   NEWSUID -> up to 10 components (PRIDn, PROPn %)
    SOTERparameterEstimates PRID, Layer D1..D5, TopDep/BotDep (cm), SDTO,
                           STPC, CLPC (%), TOTC (g/kg), BULK (g/cm3),
                           CFRAG (vol %), Drain (FAO class)
    SOTERflagTTRrules      taxotransfer rules applied per PRID/layer
The polygon layer (e.g. KEN_SOTWISv1_t1s1d1.shp) is only used for geometry
and its NEWSUID. For each component the properties are depth-weighted over
[depth_top, depth_bottom] (default 0-20 cm, decision D6) using the layers
that overlap the interval; K is computed PER COMPONENT (Williams is
non-linear) and then weighted by component share. Values < 0 (SOTWIS uses
-1 for missing) are treated as missing; missing components are left out
and the remaining shares renormalised, and the unit's `share_with_data`
records how much of the unit had data.

Generic polygon attributes
--------------------------
Any polygon layer whose attributes already hold sand/silt/clay % and OC %
(+ optional bulk density, rock %, drainage): one "component" per polygon.

No QGIS imports, no GDAL.
"""

import math
import os
import sqlite3

from .usle_k import williams_k, dg_k, cfrg, usda_texture, hsg_proxy

NAN = float("nan")
SOTWIS_KENYA = "ISRIC SOTWIS Kenya v1.0 (Batjes & Gicheru 2004), SQLite"

PROPS = ("sand", "silt", "clay", "oc_pct", "bulk", "cfrag")


class UnitSoil(object):
    """Topsoil description of one soil (mapping) unit."""

    __slots__ = ("unit_id", "sand", "silt", "clay", "oc_pct", "bulk", "cfrag",
                 "k_us", "k_si", "k_dg", "cfrg", "texture", "drain", "hsg",
                 "share_with_data", "n_components", "dominant_prid", "ttr_main")

    def __init__(self, unit_id):
        self.unit_id = unit_id
        for s in self.__slots__[1:]:
            setattr(self, s, None)

    def as_dict(self):
        return {s: getattr(self, s) for s in self.__slots__}


def _num(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return NAN
    return v if math.isfinite(v) and v >= 0 else NAN


def _depth_weighted(layers, top, bottom):
    """layers: list of dict(top, bot, values...). Weighted by overlap (cm)."""
    acc = {p: [0.0, 0.0] for p in PROPS}
    drain = None
    for L in layers:
        ov = min(bottom, L["bot"]) - max(top, L["top"])
        if ov <= 0:
            continue
        drain = drain or L.get("drain")
        for p in PROPS:
            v = L.get(p, NAN)
            if math.isfinite(v):
                acc[p][0] += v * ov
                acc[p][1] += ov
    out = {p: (acc[p][0] / acc[p][1] if acc[p][1] > 0 else NAN) for p in PROPS}
    out["drain"] = drain
    return out


def _component_record(props):
    k_us, k_si, _ = williams_k(props["sand"], props["silt"], props["clay"], props["oc_pct"])
    rec = dict(props)
    rec.update(k_us=k_us, k_si=k_si, k_dg=dg_k(props["sand"], props["silt"], props["clay"]))
    return rec


def _combine(unit_id, comps):
    """comps: list of (share %, component record, prid, ttr). -> UnitSoil."""
    u = UnitSoil(unit_id)
    u.n_components = len(comps)
    total = sum(s for s, _, _, _ in comps) or 0.0
    good = [(s, c, p, t) for s, c, p, t in comps if math.isfinite(c.get("k_si", NAN))]
    wsum = sum(s for s, _, _, _ in good)
    u.share_with_data = 100.0 * wsum / total if total > 0 else 0.0
    if wsum <= 0:
        return u
    for name in PROPS + ("k_us", "k_si", "k_dg"):
        num = den = 0.0
        for s, c, _, _ in good:
            v = c.get(name, NAN)
            if math.isfinite(v):
                num += s * v; den += s
        setattr(u, name, num / den if den > 0 else NAN)
    dom = max(good, key=lambda z: z[0])
    u.dominant_prid, u.ttr_main = dom[2], dom[3]
    # drainage: class with the largest share
    by = {}
    for s, c, _, _ in good:
        if c.get("drain"):
            by[c["drain"]] = by.get(c["drain"], 0.0) + s
    u.drain = max(by, key=by.get) if by else None
    u.texture = usda_texture(u.sand, u.silt, u.clay)
    u.hsg = hsg_proxy(u.texture, u.drain)
    u.cfrg = cfrg(u.cfrag)
    return u


def load_sotwis(db_path, depth_top=0.0, depth_bottom=20.0):
    """Read a SOTWIS SQLite database. Returns {NEWSUID: UnitSoil}, info."""
    if depth_bottom <= depth_top:
        raise ValueError("depth_bottom must be greater than depth_top")
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        prof = {}
        for row in con.execute(
                "SELECT PRID, Layer, TopDep, BotDep, SDTO, STPC, CLPC, TOTC, BULK, "
                "CFRAG, Drain FROM SOTERparameterEstimates"):
            prid, layer, top, bot, sdto, stpc, clpc, totc, bulk, cfr, drain = row
            prof.setdefault(prid, []).append(dict(
                top=_num(top), bot=_num(bot), sand=_num(sdto), silt=_num(stpc),
                clay=_num(clpc), oc_pct=_num(totc) / 10.0, bulk=_num(bulk),
                cfrag=_num(cfr), drain=(drain or None)))
        ttr = {}           # TTRmain of the topmost rule row per profile
        try:
            best = {}
            for prid, top, ttrmain in con.execute(
                    "SELECT PRID, Newtopdep, TTRmain FROM SOTERflagTTRrules"):
                t = _num(top)
                if prid not in best or (math.isfinite(t) and t < best[prid]):
                    best[prid] = t if math.isfinite(t) else 1e9
                    ttr[prid] = ttrmain
        except sqlite3.OperationalError:
            pass
        comp_rows = con.execute("SELECT * FROM SOTERunitComposition").fetchall()
        cols = [d[1] for d in con.execute("PRAGMA table_info(SOTERunitComposition)")]
    finally:
        con.close()

    cache = {}
    units = {}
    for row in comp_rows:
        r = dict(zip(cols, row))
        uid = r.get("NEWSUID")
        comps = []
        for n in range(1, 11):
            prid, share = r.get(f"PRID{n}"), r.get(f"PROP{n}")
            if not prid or share is None or share <= 0:
                continue
            if prid not in cache:
                cache[prid] = _component_record(
                    _depth_weighted(prof.get(prid, []), depth_top, depth_bottom))
            comps.append((float(share), cache[prid], prid, ttr.get(prid)))
        units[uid] = _combine(uid, comps)
    info = {"soil_dataset": SOTWIS_KENYA, "soil_depth_cm": f"{depth_top:g}-{depth_bottom:g}",
            "n_units": len(units), "n_profiles": len(prof)}
    return units, info


def load_attributes(records, unit_field="id", sand="sand", silt="silt", clay="clay",
                    oc_pct="oc", bulk=None, cfrag=None, drain=None, oc_is_gkg=False,
                    label="user polygon attributes"):
    """Generic dataset: records = list of dicts (one per polygon/unit)."""
    units = {}
    for r in records:
        props = {"sand": _num(r.get(sand)), "silt": _num(r.get(silt)),
                 "clay": _num(r.get(clay)),
                 "oc_pct": _num(r.get(oc_pct)) / (10.0 if oc_is_gkg else 1.0),
                 "bulk": _num(r.get(bulk)) if bulk else NAN,
                 "cfrag": _num(r.get(cfrag)) if cfrag else NAN,
                 "drain": (str(r.get(drain)).strip() if drain and r.get(drain) is not None else None)}
        uid = r.get(unit_field)
        units[uid] = _combine(uid, [(100.0, _component_record(props), None, None)])
    return units, {"soil_dataset": label, "soil_depth_cm": "as supplied",
                   "n_units": len(units)}


def load_csv_units(path, unit_field):
    """Soil properties per unit from a CSV table (v0.14), for polygons that carry
    only a unit code. Columns (case-insensitive): the unit field, sand, silt,
    clay, oc (percent), optional bulk, cfrag, drain.

    Returns load_attributes(...) output keyed by the unit code as TEXT, so a
    polygon attribute 12 and a CSV value "12" match.
    """
    import csv
    with open(path, newline="", encoding="utf-8-sig") as f:
        sample = f.read(4096)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        except csv.Error:
            dialect = csv.excel
        rows = list(csv.DictReader(f, dialect=dialect))
    if not rows:
        raise ValueError(f"{path}: no rows")
    lower = {k.strip().lower(): k for k in rows[0].keys() if k}
    need = [unit_field.lower(), "sand", "silt", "clay", "oc"]
    missing = [n for n in need if n not in lower]
    if missing:
        raise ValueError(f"{path}: missing column(s) {', '.join(missing)} "
                         "(need the unit field, sand, silt, clay, oc in percent)")
    recs, seen = [], set()
    for r in rows:
        key = str(r[lower[unit_field.lower()]]).strip()
        if not key:
            continue
        if key in seen:
            raise ValueError(f"{path}: unit '{key}' appears more than once")
        seen.add(key)
        rec = {"id": key}
        for name in ("sand", "silt", "clay", "oc", "bulk", "cfrag", "drain"):
            if name in lower:
                rec[name] = r[lower[name]]
        recs.append(rec)
    return load_attributes(recs, sand="sand", silt="silt", clay="clay", oc_pct="oc",
                           bulk="bulk" if "bulk" in lower else None,
                           cfrag="cfrag" if "cfrag" in lower else None,
                           drain="drain" if "drain" in lower else None,
                           label=f"soil table {os.path.basename(path)} (joined on {unit_field})")
