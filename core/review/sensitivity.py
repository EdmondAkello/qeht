# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Per-crossing sensitivity table (1.2).

How much does a crossing's catchment area and longest flow path move when
one modelling choice changes? Each variant changes a single setting from
the run and re-measures every crossing:

    flat method       the other flat method (Barnes <-> toward lower)
    snap radius       half and twice the run's radius (at least 1 cell)
    stream threshold  half and twice the run's threshold (snapping target)
    DEM error         P10 and P90 from the DEM-uncertainty run, when the
                      package has one (copied, not recalculated)

The base row is the run's own outlet measured on the run's routing, so it
equals the package. The snap-radius and threshold variants snap the
crossing's input point again; the flat-method variant snaps the run's
outlet onto the other method's streams. Area is the contributing
area at the snapped outlet; the flow path is the longest flow length
arriving there (the full-catchment longest flow path). `outlet_moved` is
1 when the snapped outlet is more than one cell from the base outlet;
`outlet_lost` when no stream is within the snap radius.

Long table: one row per crossing and variant. Summary per crossing: the
largest |change| in area and in flow path, and the factor behind it.

No QGIS imports, no GDAL.
"""

import math

import numpy as np

from ..flow.direction import d8_direction
from ..flow.accumulation import flow_accumulation
from ..watershed.delineate import extract_streams, snap_pour_point
from ..watershed.uncertainty import upstream_lengths

SENSITIVITY_FIELDS = [
    ("outlet_uid", "str"), ("factor", "str"), ("setting", "str"), ("area_km2", "float"),
    ("d_area_pct", "float"), ("lfp_length_m", "float"), ("d_lfp_pct", "float"),
    ("outlet_moved", "int"), ("outlet_lost", "int"),
]
SUMMARY_FIELDS = [
    ("outlet_uid", "str"), ("base_area_km2", "float"), ("base_lfp_length_m", "float"),
    ("sens_max_area_pct", "float"), ("sens_area_factor", "str"),
    ("sens_max_lfp_pct", "float"), ("sens_lfp_factor", "str"), ("sens_moved", "str"),
]


class _Routing(object):
    def __init__(self, filled, valid, cw, ch, method, direction=None):
        if direction is None:
            direction, _ = d8_direction(filled, valid, cw, ch, resolve_flats=True,
                                        flat_method=method)
        self.direction = direction
        self.acc, _ = flow_accumulation(self.direction, valid)
        self.length, _ = upstream_lengths(self.direction, valid, cw, ch)
        self._streams = {}
        self.valid = valid

    def streams(self, thr):
        key = float(thr)
        if key not in self._streams:
            self._streams[key] = extract_streams(self.acc, self.valid, threshold_cells=thr)
        return self._streams[key]

    def at(self, rc, cell_km2):
        r, c = rc
        if not self.valid[r, c]:
            return None
        return {"cell": (r, c), "area_km2": (float(self.acc[r, c]) + 1.0) * cell_km2,
                "lfp_length_m": float(self.length[r * self.acc.shape[1] + c])}

    def measure(self, rc, snap, thr, cell_km2):
        s = self.streams(thr)
        r, c, _, _ = snap_pour_point(rc[0], rc[1], self.acc, self.valid,
                                     search_radius_cells=snap, stream_mask=s)
        if snap > 0 and not s[r, c]:
            return None
        return {"cell": (r, c), "area_km2": (float(self.acc[r, c]) + 1.0) * cell_km2,
                "lfp_length_m": float(self.length[r * self.acc.shape[1] + c])}


def _pct(v, base):
    if v is None or base is None or not base:
        return None
    return 100.0 * (v - base) / base


def crossing_sensitivity(filled, valid, cell_width, cell_height, inputs_rc, uids,
                         threshold_cells, snap_cells=5, flat_method="barnes", dem_error=None,
                         base_direction=None, outlets_rc=None, progress=None):
    """-> (rows, summary).

    inputs_rc: input (pre-snap) cell of each crossing; uids: their IDs.
    dem_error: optional {uid: {"p10": (area, lfp), "p90": (area, lfp)}}.
    base_direction: the run's own D8 directions (decoded), to skip routing
    the run's flat method again.
    outlets_rc: the run's outlet cells. With them the base is measured at
    the run's outlet (so it equals the package), the flat-method variant
    snaps that outlet onto the other method's streams, and the snap-radius
    and threshold variants snap the input point. Without them every row
    snaps the input point, the base with the run's settings.
    """
    cw, ch = abs(cell_width), abs(cell_height)
    cell_km2 = cw * ch / 1e6
    other = "toward" if flat_method == "barnes" else "barnes"
    thr = max(1.0, float(threshold_cells))
    snap = int(snap_cells)
    half_snap, dbl_snap = max(1, snap // 2), max(2, 2 * snap)
    if progress:
        progress(0.0, f"Routing with {flat_method}")
    base_r = _Routing(filled, valid, cw, ch, flat_method, base_direction)
    if progress:
        progress(0.45, f"Routing with {other}")
    other_r = _Routing(filled, valid, cw, ch, other)
    if progress:
        progress(0.9, "Measuring crossings")
    variants = [
        ("flat method", other, other_r, max(1, snap), thr),
        ("snap radius", f"{half_snap} cells", base_r, half_snap, thr),
        ("snap radius", f"{dbl_snap} cells", base_r, dbl_snap, thr),
        ("stream threshold", f"{thr / 2:g} cells", base_r, snap, max(1.0, thr / 2)),
        ("stream threshold", f"{thr * 2:g} cells", base_r, snap, thr * 2),
    ]
    rows, summary = [], []
    outs = outlets_rc if outlets_rc is not None else [None] * len(uids)
    for uid, rc, orc in zip(uids, inputs_rc, outs):
        b = base_r.at(orc, cell_km2) if orc is not None else base_r.measure(rc, snap, thr, cell_km2)
        ba = b["area_km2"] if b else None
        bl = b["lfp_length_m"] if b else None
        rows.append({"outlet_uid": uid, "factor": "base",
                     "setting": (f"run outlet, {flat_method}" if orc is not None else
                                 f"{flat_method}, snap {snap} cells, threshold {thr:g} cells"),
                     "area_km2": ba, "d_area_pct": 0.0 if b else None, "lfp_length_m": bl,
                     "d_lfp_pct": 0.0 if b else None, "outlet_moved": 0, "outlet_lost": int(b is None)})
        mine = []
        for factor, setting, rt, sn, th in variants:
            start = orc if (factor == "flat method" and orc is not None) else rc
            m = rt.measure(start, sn, th, cell_km2)
            moved = 0
            if m and b:
                moved = int(max(abs(m["cell"][0] - b["cell"][0]), abs(m["cell"][1] - b["cell"][1])) > 1)
            row = {"outlet_uid": uid, "factor": factor, "setting": setting,
                   "area_km2": m["area_km2"] if m else None,
                   "d_area_pct": _pct(m["area_km2"], ba) if m else None,
                   "lfp_length_m": m["lfp_length_m"] if m else None,
                   "d_lfp_pct": _pct(m["lfp_length_m"], bl) if m else None,
                   "outlet_moved": moved, "outlet_lost": int(m is None)}
            rows.append(row)
            mine.append(row)
        for q in ("p10", "p90"):
            v = (dem_error or {}).get(uid, {}).get(q)
            if v and v[0] is not None:
                row = {"outlet_uid": uid, "factor": "DEM error", "setting": q.upper(),
                       "area_km2": v[0], "d_area_pct": _pct(v[0], ba),
                       "lfp_length_m": v[1], "d_lfp_pct": _pct(v[1], bl),
                       "outlet_moved": None, "outlet_lost": None}
                rows.append(row)
                mine.append(row)

        def worst(key):
            vals = [(abs(r[key]), r) for r in mine if r[key] is not None and math.isfinite(r[key])]
            if not vals:
                return None, ""
            v, r = max(vals, key=lambda t: t[0])
            return v, f"{r['factor']} ({r['setting']})"
        wa, fa = worst("d_area_pct")
        wl, fl = worst("d_lfp_pct")
        moved = sorted({f"{r['factor']} ({r['setting']})" for r in mine
                        if r.get("outlet_moved") or r.get("outlet_lost")})
        summary.append({"outlet_uid": uid, "base_area_km2": ba, "base_lfp_length_m": bl,
                        "sens_max_area_pct": wa, "sens_area_factor": fa,
                        "sens_max_lfp_pct": wl, "sens_lfp_factor": fl,
                        "sens_moved": "; ".join(moved)})
    if progress:
        progress(1.0, "Sensitivity done")
    return rows, summary


def package_inputs(package_path, geotransform, shape):
    """Crossing input and outlet cells, uids, run settings and DEM-error
    ranges from a package (proposed crossings are left out: the coverage
    check delineates them against the embankment, not from the routing) -> (uids, inputs_rc, settings, dem_error,
    outlets_rc)."""
    from ..interop import gpkg
    md = {r["key"]: r["value"] for r in gpkg.read_table(package_path, "qeht_run_metadata")}
    rows = gpkg.read_table(package_path, "crossings", with_geometry=False)
    gt = geotransform
    uids, rc, orc, err = [], [], [], {}

    def cell(x, y):
        r, k = int((y - gt[3]) // gt[5]), int((x - gt[0]) // gt[1])
        return (r, k) if 0 <= r < shape[0] and 0 <= k < shape[1] else None
    for c in rows:
        if c.get("status") == "proposed":      # area from the coverage check's own delineation
            continue
        o = cell(c["outlet_x"], c["outlet_y"])
        i = cell(c["input_x"], c["input_y"]) if c.get("input_x") is not None else o
        if o is None or i is None:
            continue
        uids.append(c["outlet_uid"])
        rc.append(i)
        orc.append(o)
        if c.get("unc_area_p10") is not None:
            err[c["outlet_uid"]] = {"p10": (c.get("unc_area_p10"), c.get("unc_lfp_p10")),
                                    "p90": (c.get("unc_area_p90"), c.get("unc_lfp_p90"))}

    def num(key, default):
        try:
            return float(md.get(key))
        except (TypeError, ValueError):
            return default
    method = str(md.get("flat_method") or "").lower()
    settings = {"threshold_cells": num("stream_threshold_cells", 200.0),
                "snap_cells": int(num("snap_radius_cells", 5)),
                "flat_method": "toward" if method.startswith("toward") else "barnes"}
    return uids, rc, settings, err, orc


def by_uid(rows):
    out = {}
    for r in rows:
        if r["factor"] != "base":
            out.setdefault(r["outlet_uid"], []).append(r)
    return out


def write_csv(rows, path, fields=SENSITIVITY_FIELDS):
    import csv
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=[k for k, _ in fields])
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else
                            (f"{r[k]:.10g}" if isinstance(r[k], float) else r[k]))
                        for k, _ in fields})
    return path
