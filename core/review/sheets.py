# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Crossing review sheets (1.2): one printable page per crossing.

Each page carries what a reviewer needs to accept or question one
crossing without opening a GIS: a plan of the catchment, flow path, road
and any flagged fill depression; the main characteristics; the checks that
bear on it (flat method, DEM uncertainty, mapped drainage, sensitivity
table when one was run); its review-register findings; and a sign-off
box. One HTML file, page breaks between crossings; print to PDF from a
browser.

Read from the package and the pipeline's `extra`; nothing is recalculated.
No QGIS imports, no GDAL.
"""

import datetime
import math
import zlib

from ..interop import gpkg
from ..report.run_report import CSS, _e, _fmt, _has, _lines, _table
from .register import SEVERITY_ORDER

SHEET_CSS = CSS + """
.sheet{break-after:page;page-break-after:always;padding-top:6px}
.sheet:last-child{break-after:auto;page-break-after:auto}
.grid2{display:grid;grid-template-columns:1.15fr 1fr;gap:16px;align-items:start}
.sev{display:inline-block;padding:0 6px;border-radius:4px;font-size:11px;font-weight:600;color:#fff}
.sev.high{background:var(--warn)}.sev.medium{background:var(--amber)}
.sev.low{background:var(--accent)}.sev.info{background:#7a8494}
.sign td{height:30px}
@media print{.grid2{grid-template-columns:1.15fr 1fr}}
"""

KEY_FIELDS = [
    ("area_km2", "Catchment area (km²)", "catch"),
    ("lfp_length_m", "Longest flow path (m)", "fp"),
    ("lfp_slope_1085", "10–85 slope (m/m)", "fp"),
    ("lfp_drop_m", "Flow-path drop (m)", "fp"),
    ("catch_slope_horn", "Mean catchment slope (m/m)", "catch"),
    ("relief_m", "Relief (m)", "catch"),
    ("elev_min_m", "Lowest elevation (m)", "catch"),
    ("elev_max_m", "Highest elevation (m)", "catch"),
    ("max_strahler", "Highest Strahler order", "catch"),
    ("cn_export", "Curve number", "catch"),
    ("rational_c", "Rational C", "catch"),
    ("soil_hsg", "Hydrologic soil group", "catch"),
    ("map_mm", "Mean annual rainfall (mm)", "catch"),
    ("rain_zone", "Rainfall zone", "catch"),
    ("tc_kirpich_min", "Tc Kirpich (min)", "x"),
    ("tc_kerby_kirpich_min", "Tc Kerby–Kirpich (min)", "x"),
    ("tc_scs_lag_min", "Tc SCS lag (min)", "x"),
    ("tc_tr55_min", "Tc TR-55 (min)", "x"),
    ("tc_bransby_williams_min", "Tc Bransby Williams (min)", "x"),
    ("xs_bankfull_w_m", "Bankfull width at crossing (m)", "x"),
    ("xs_quality", "Cross-section quality", "x"),
]


def _bbox(seqs):
    pts = [p for s in seqs for p in s]
    if not pts:
        return None
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    return min(xs), min(ys), max(xs), max(ys)


def crossing_svg(uid, layers, dep_rows=(), width=520, height=420):
    """Plan around one crossing: its catchment and flow path, the road,
    other crossings and flagged depressions inside the frame."""
    catch = [r for r in layers.get("catchments", []) if r.get("outlet_uid") == uid]
    fp = [r for r in layers.get("flowpaths", []) if r.get("outlet_uid") == uid]
    xing = [r for r in layers.get("crossings", []) if r.get("outlet_uid") == uid]
    own = [ln for r in catch + fp + xing for ln in _lines(r.get("_geom"))]
    bb = _bbox(own)
    if bb is None:
        return ""
    x0, y0, x1, y1 = bb
    span = max(x1 - x0, y1 - y0, 200.0)
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    half = span * 0.58
    x0, x1, y0, y1 = cx - half, cx + half, cy - half, cy + half
    pad = 14
    s = min((width - 2 * pad) / (x1 - x0), (height - 2 * pad) / (y1 - y0))
    ox = pad + ((width - 2 * pad) - (x1 - x0) * s) / 2
    oy = pad + ((height - 2 * pad) - (y1 - y0) * s) / 2

    def xy(p):
        return f"{ox + (p[0] - x0) * s:.1f},{oy + (y1 - p[1]) * s:.1f}"

    def path(seq, close=False):
        return "M" + " L".join(xy(p) for p in seq) + (" Z" if close else "")

    def inside(seq):
        b = _bbox([seq])
        return b and not (b[2] < x0 or b[0] > x1 or b[3] < y0 or b[1] > y1)
    cid = f"clip-{zlib.crc32(str(uid).encode()):08x}"
    out = [f'<svg viewBox="0 0 {width} {height}" width="100%" role="img" '
           f'aria-label="Plan of crossing {_e(uid)}">',
           f'<defs><clipPath id="{cid}"><rect x="0" y="0" width="{width}" height="{height}"/>'
           '</clipPath></defs>', f'<rect width="{width}" height="{height}" fill="#fff" '
           'stroke="#d9dee7"/>', f'<g clip-path="url(#{cid})">']
    for r in catch:
        for ring in _lines(r.get("_geom")):
            out.append(f'<path d="{path(ring, True)}" fill="#e6f1f5" fill-opacity="0.7" '
                       'stroke="#1f6f8b" stroke-width="1.2"/>')
    for d in dep_rows:
        if not d.get("flag"):
            continue
        for poly in d.get("_geom") or []:
            for ring in poly:
                if inside(ring):
                    out.append(f'<path d="{path(ring, True)}" fill="#c0392b" fill-opacity="0.35" '
                               'stroke="#c0392b" stroke-width="0.8"/>')
    for r in fp:
        for ln in _lines(r.get("_geom")):
            out.append(f'<path d="{path(ln)}" fill="none" stroke="#3a7bd5" stroke-width="1.6"/>')
    for r in layers.get("road_alignment", []):
        for ln in _lines(r.get("_geom")):
            if inside(ln):
                out.append(f'<path d="{path(ln)}" fill="none" stroke="#1f2a44" stroke-width="3"/>')
    for r in layers.get("crossings", []):
        for ln in _lines(r.get("_geom")):
            p = ln[0]
            if not (x0 <= p[0] <= x1 and y0 <= p[1] <= y1):
                continue
            x, y = (float(v) for v in xy(p).split(","))
            me = r.get("outlet_uid") == uid
            out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{6 if me else 4}" '
                       f'fill="{"#c0392b" if me else "#7a8494"}" stroke="#fff" stroke-width="1.5"/>')
            out.append(f'<text x="{x + 8:.1f}" y="{y - 6:.1f}" font-size="{12 if me else 10}" '
                       f'fill="#1d2330" paint-order="stroke" stroke="#fff" stroke-width="3">'
                       f'{_e(r.get("outlet_uid"))}</text>')
    out.append("</g>")
    out.append(f'<g transform="translate({width - 26},{30})"><path d="M0,-16 L6,5 L0,1 L-6,5 Z" '
               'fill="#1f2a44"/><text x="0" y="18" font-size="10" text-anchor="middle" '
               'fill="#1f2a44">N</text></g>')
    target = (x1 - x0) / 4
    mag = 10 ** math.floor(math.log10(target))
    bar = max(m * mag for m in (1, 2, 5) if m * mag <= target)
    lab = f"{bar / 1000:g} km" if bar >= 1000 else f"{bar:g} m"
    out.append(f'<g transform="translate({pad},{height - 12})"><rect width="{bar * s:.1f}" '
               'height="4" fill="#1f2a44"/>'
               f'<text x="{bar * s + 6:.1f}" y="5" font-size="10" fill="#1f2a44">{lab}</text></g>')
    out.append("</svg>")
    return "\n".join(out)


def _kv(rows):
    return ("<table><tbody>" + "".join(
        f"<tr><td class='kv'>{_e(k)}</td><td class='n'>{v}</td></tr>" for k, v in rows)
        + "</tbody></table>")


def build_sheets(package_path, register=(), extra=None, depressions=(), title=None):
    """-> HTML text, one .sheet section per crossing in chainage order.

    extra: flat_check {uid: {...}} and sensitivity {uid: [rows]} (rows
    with factor, setting, area_km2, lfp_length_m, d_area_pct, d_lfp_pct,
    outlet_moved), both optional.
    """
    extra = dict(extra or {})
    tables = {t for t, _ in gpkg.list_tables(package_path)}
    md = {r["key"]: r["value"] for r in gpkg.read_table(package_path, "qeht_run_metadata")}
    layers = {n: gpkg.read_table(package_path, n)
              for n in ("crossings", "catchments", "flowpaths", "road_alignment") if n in tables}
    catch = {r["outlet_uid"]: r for r in layers.get("catchments", [])}
    fps = {r["outlet_uid"]: r for r in layers.get("flowpaths", [])}
    flat = extra.get("flat_check") or {}
    sens = extra.get("sensitivity") or {}
    by_uid = {}
    for f in register or ():
        by_uid.setdefault(f.get("outlet_uid") or "", []).append(f)
    xs = sorted(layers.get("crossings", []),
                key=lambda r: (r.get("chainage_m") if _has(r.get("chainage_m")) else math.inf,
                               str(r.get("outlet_uid"))))
    gen = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    title = title or "QEHT crossing review sheets"
    P = []
    for c in xs:
        uid = c.get("outlet_uid")
        k, fp = catch.get(uid, {}), fps.get(uid, {})
        head = (f"<h1>{_e(uid)}</h1><p class='sub'>{_e(title)} · chainage "
                f"{_fmt(c.get('chainage_m'))} m · {_e(c.get('status') or 'existing')} · QEHT "
                f"{_e(md.get('qeht_version'))} · run {_e(md.get('run_utc'))} · generated {gen}</p>")
        loc = [("Crossing ID", _e(c.get("crossing_id"))),
               ("Outlet x, y", f"{_fmt(c.get('outlet_x'))}, {_fmt(c.get('outlet_y'))}"),
               ("Snapped from input (m)", _fmt(c.get("snap_dist_m"))),
               ("Stream order at outlet", _fmt(c.get("stream_order"))),
               ("Link", f"{_e(c.get('link_method'))} ({_fmt(c.get('link_confidence'))})")]
        if c.get("status") == "proposed":
            loc.append(("Proposed because", _e(c.get("proposed_reason"))))
        src = {"catch": k, "fp": fp, "x": c}
        key = [(lbl, _fmt(src[w].get(f)) if not isinstance(src[w].get(f), str)
                else _e(src[w].get(f)))
               for f, lbl, w in KEY_FIELDS if _has(src[w].get(f))]
        if k.get("truncated"):
            key.insert(0, ("Catchment", "<span class='flag'>truncated by DEM edge / NoData</span>"))
        svg = crossing_svg(uid, layers, depressions)
        P.append("<section class='sheet'>" + head + "<div class='grid2'><div>"
                 + (f"<figure>{svg}<figcaption>Catchment (blue), longest flow path, road, this "
                    "crossing (red) and nearby crossings (grey); flagged fill depressions in red "
                    "shading.</figcaption></figure>" if svg else "")
                 + "<h3>Location</h3>" + _kv(loc) + "</div><div><h3>Characteristics</h3>"
                 + _kv(key) + "</div></div>")

        checks = []
        v = flat.get(uid)
        if v:
            checks.append(("Area, Barnes / toward-lower (km²)",
                           f"{_fmt(v.get('area_barnes_km2'))} / {_fmt(v.get('area_toward_km2'))}"))
            checks.append(("Flat-method change (%)", _fmt(v.get("flat_sensitivity_pct"))
                           + (" <span class='flag'>check</span>" if v.get("flat_sensitive") else "")))
            if v.get("flat_moved"):
                checks.append(("Channel under other flat method", "<span class='flag'>moves</span>"))
            if _has(v.get("flat_catch_share_pct")):
                checks.append(("Share of catchment on resolved flats (%)",
                               _fmt(v.get("flat_catch_share_pct"))))
        if _has(c.get("unc_area_p50")):
            checks.append(("Area under DEM error P10 / P50 / P90 (km²)",
                           f"{_fmt(c.get('unc_area_p10'))} / {_fmt(c.get('unc_area_p50'))} / "
                           f"{_fmt(c.get('unc_area_p90'))}"))
            checks.append(("Tc under DEM error P10 / P90 (min)",
                           f"{_fmt(c.get('unc_tc_p10'))} / {_fmt(c.get('unc_tc_p90'))}"))
            checks.append(("Realisations switching / losing stream (%)",
                           f"{_fmt(c.get('unc_switch_pct'))} / {_fmt(c.get('unc_lost_pct'))}"))
        if c.get("map_agrees") is not None and c.get("map_agrees") != "":
            checks.append(("On mapped drainage", "yes" if c.get("map_agrees") else
                           f"<span class='flag'>no</span> ({_fmt(c.get('map_river_dist_m'))} m)"))
        if checks:
            P.append("<h3>Checks</h3>" + _kv(checks))
        srows = sens.get(uid) or []
        if srows:
            cols = ["factor", "setting", "area_km2", "d_area_pct", "lfp_length_m", "d_lfp_pct",
                    "outlet_moved"]
            P.append("<h3>Sensitivity</h3><p class='sub'>Each row changes one setting from the "
                     "run; changes are relative to the run's own result.</p>"
                     + _table(cols, [{col: (_e(r.get(col)) if isinstance(r.get(col), str)
                                            else ("yes" if col == "outlet_moved" and r.get(col)
                                                  else ("–" if col == "outlet_moved"
                                                        else _fmt(r.get(col)))))
                                      for col in cols} for r in srows],
                              {"area_km2", "d_area_pct", "lfp_length_m", "d_lfp_pct"}))
        finds = sorted(by_uid.get(uid, []), key=lambda f: SEVERITY_ORDER[f["severity"]])
        if finds:
            P.append("<h3>Review register</h3>" + _table(
                ["finding_id", "severity", "finding", "evidence", "check"],
                [{"finding_id": _e(f["finding_id"]),
                  "severity": f"<span class='sev {f['severity']}'>{f['severity']}</span>",
                  "finding": _e(f["finding"]), "evidence": _e(f["evidence"]),
                  "check": _e(f["check"])} for f in finds]))
        else:
            P.append("<h3>Review register</h3><p>No findings for this crossing.</p>")
        P.append("<h3>Review</h3><table class='sign'><tbody><tr><td class='kv'>Decision "
                 "(accept / revise / site check)</td><td></td></tr><tr><td class='kv'>Reviewer, "
                 "date</td><td></td></tr><tr><td class='kv'>Notes</td><td></td></tr>"
                 "</tbody></table></section>")
    if not P:
        P.append("<p>No crossings in the package.</p>")
    return ("<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'><title>"
            f"{_e(title)}</title><style>{SHEET_CSS}</style></head><body><main>"
            + "".join(P) + "</main></body></html>")


def write_sheets(package_path, out_html, register=(), extra=None, depressions=(), title=None):
    text = build_sheets(package_path, register, extra, depressions, title)
    with open(out_html, "w", encoding="utf-8") as f:
        f.write(text)
    return out_html


def alignment_from_package(package_path):
    """Rebuild the chainage Alignment from the package's road_alignment
    layer (None when the run had no road)."""
    from ..network.alignment import Alignment
    if "road_alignment" not in {t for t, _ in gpkg.list_tables(package_path)}:
        return None
    rows = gpkg.read_table(package_path, "road_alignment")
    parts = [ln for r in sorted(rows, key=lambda r: r.get("part") or 0)
             for ln in _lines(r.get("_geom")) if len(ln) >= 2]
    if not parts:
        return None
    start = min((r.get("chainage_from_m") for r in rows if _has(r.get("chainage_from_m"))),
                default=0.0)
    return Alignment(parts, start_chainage=start)
