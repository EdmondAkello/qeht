# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Run report (F8, v0.17): one self-contained HTML page per design
hydrology package, for a design-report appendix or a reviewer.

Built only from the package (qeht_run_metadata and its layers) plus the
optional `extra` dict the pipeline passes (warnings, DEM audit, flat-method
check); nothing is recalculated. Print it to PDF from a browser.

No QGIS imports, no GDAL.
"""

import datetime
import html
import json
import math

from ..interop import gpkg
from .characteristics import characteristics_rows

CSS = """
:root{--ink:#1d2330;--muted:#5b6575;--line:#d9dee7;--soft:#f4f6f9;--navy:#1f2a44;
--accent:#1f6f8b;--warn:#c0392b;--amber:#d9822b;--ok:#3a9e5f}
*{box-sizing:border-box}body{margin:0;background:#fff;color:var(--ink);
font:14px/1.5 Inter,"Segoe UI",Roboto,"Liberation Sans",sans-serif}
main{max-width:1100px;margin:0 auto;padding:28px 24px 48px}
h1{font-size:26px;color:var(--navy);margin:0 0 4px}h2{font-size:18px;color:var(--navy);
margin:30px 0 10px;padding-bottom:5px;border-bottom:2px solid var(--navy);break-after:avoid}
h3{font-size:15px;color:var(--navy);margin:18px 0 6px}
.sub{color:var(--muted);margin:0 0 14px}.note{background:var(--soft);border-left:3px solid
var(--accent);padding:8px 12px;color:var(--muted);font-size:12.5px;margin:10px 0}
.cards{display:grid;grid-template-columns:repeat(auto-fill,minmax(160px,1fr));gap:10px;margin:14px 0}
.card{border:1px solid var(--line);border-radius:8px;padding:9px 12px}.card .k{font-size:11.5px;
color:var(--muted)}.card .v{font-size:17px;font-weight:600;color:var(--navy)}
table{width:100%;border-collapse:collapse;font-size:12.5px;margin:6px 0 12px}
th{background:var(--navy);color:#fff;text-align:left;font-weight:600;padding:5px 7px}
td{padding:4px 7px;border-bottom:1px solid var(--line);vertical-align:top}
tbody tr:nth-child(even) td{background:#fafbfc}td.n,th.n{text-align:right;white-space:nowrap}
td.kv{color:var(--muted);width:30%}.flag{color:var(--warn);font-weight:600}
.tag{display:inline-block;padding:0 6px;border-radius:4px;font-size:11px;font-weight:600;color:#fff}
.tag.p{background:var(--amber)}.tag.e{background:var(--accent)}.tag.w{background:var(--warn)}
ul.warn li{margin-bottom:4px}code{background:var(--soft);padding:0 4px;border-radius:3px;
font-size:12px}figure{margin:8px 0 14px;border:1px solid var(--line);border-radius:8px;
padding:10px;background:#fcfcfd}figcaption{font-size:12px;color:var(--muted);margin-top:4px}
svg text{font-family:Inter,"Segoe UI",sans-serif}
@media print{main{padding:0}h2{break-after:avoid}table,figure{break-inside:avoid}}
"""

PROVENANCE = [
    ("dem_source", "DEM source (user-declared)"), ("raw_dem_path", "DEM file"),
    ("dem_sha256", "DEM SHA-256"), ("cell_size_m", "Cell size (m)"),
    ("crs_name", "CRS"), ("conditioning", "Conditioning"), ("flat_method", "Flat method"),
    ("tie_rule", "D8 tie rule"), ("stream_threshold_km2", "Stream threshold (km²)"),
    ("snap_strategy", "Outlet snapping"), ("catchment_mode", "Catchments"),
    ("lfp_1085_reference", "10–85 slope reference"), ("sheet_cap_m", "Sheet-flow cap (m)"),
    ("channel_slope_m", "Channel slope length (m)"), ("crossing_source", "Crossings from"),
    ("alignment_source", "Road alignment"), ("chainage_start_m", "Start chainage (m)"),
    ("soil_dataset", "Soil dataset"), ("soil_depth_cm", "Soil depth (cm)"),
    ("soil_hsg_source", "Hydrologic soil groups"), ("soil_k_source", "USLE K"),
    ("conditioning_burn", "Embankment breaches"), ("relinked_from", "Relinked from"),
]


def _e(v):
    return html.escape("" if v is None else str(v))


def _has(v):
    return v is not None and v != "" and not (isinstance(v, float) and not math.isfinite(v))


def _fmt(v, nd=3):
    if not _has(v):
        return "–"
    if isinstance(v, float):
        if abs(v) >= 1000:
            return f"{v:,.0f}"
        if abs(v) >= 100:
            return f"{v:,.1f}"
        return f"{v:.{nd}g}" if abs(v) < 1 else f"{v:.{nd}f}".rstrip("0").rstrip(".")
    return _e(v)


def _json(md, key):
    try:
        return json.loads(md.get(key) or "") or {}
    except (TypeError, ValueError):
        return {}


def _lines(geom):
    """Flatten a parsed geometry into a list of coordinate sequences."""
    if geom is None:
        return []
    if isinstance(geom, tuple) and len(geom) == 2 and isinstance(geom[0], (int, float)):
        return [[geom]]
    if geom and isinstance(geom[0], tuple) and isinstance(geom[0][0], (int, float)):
        return [list(geom)]
    out = []
    for g in geom:
        out += _lines(g)
    return out


def plan_svg(package_path, tables, width=1040, height=560, background=None):
    """Schematic plan: catchment outlines, longest flow paths, road and
    crossings (existing blue, proposed orange), north arrow and scale bar."""
    layers = {}
    for name in ("catchments", "flowpaths", "crossings", "road_alignment", "sag_points"):
        if name in tables:
            layers[name] = gpkg.read_table(package_path, name)
    pts = [p for lyr in layers.values() for r in lyr for ln in _lines(r.get("_geom")) for p in ln]
    if not pts:
        return ""
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    span = max(x1 - x0, y1 - y0, 1.0)
    pad = 30
    s = min((width - 2 * pad) / max(x1 - x0, span * 1e-3), (height - 2 * pad) / max(y1 - y0, span * 1e-3))
    ox = pad + ((width - 2 * pad) - (x1 - x0) * s) / 2
    oy = pad + ((height - 2 * pad) - (y1 - y0) * s) / 2

    def xy(p):
        return f"{ox + (p[0] - x0) * s:.1f},{oy + (y1 - p[1]) * s:.1f}"

    def path(seq, close=False):
        return "M" + " L".join(xy(p) for p in seq) + (" Z" if close else "")
    out = [f'<svg viewBox="0 0 {width} {height}" width="100%" role="img" '
           f'aria-label="Schematic plan">', f'<rect width="{width}" height="{height}" fill="#fff"/>']
    if background and background.get("png") and background.get("extent"):
        import base64
        bx0, by0, bx1, by1 = background["extent"]
        try:
            with open(background["png"], "rb") as f:
                data = base64.b64encode(f.read()).decode("ascii")
            px0, py0 = (float(v) for v in xy((bx0, by1)).split(","))
            px1, py1 = (float(v) for v in xy((bx1, by0)).split(","))
            out.append(f'<image href="data:image/png;base64,{data}" x="{px0:.1f}" y="{py0:.1f}" '
                       f'width="{px1 - px0:.1f}" height="{py1 - py0:.1f}" '
                       'preserveAspectRatio="none" opacity="0.85"/>')
        except OSError:
            pass
    for r in sorted(layers.get("catchments", []), key=lambda r: -(r.get("area_km2") or 0)):
        prop = r.get("status") == "proposed"
        for ring in _lines(r.get("_geom")):
            out.append(f'<path d="{path(ring, True)}" fill="{"#fbeee0" if prop else "#e6f1f5"}" '
                       f'fill-opacity="0.55" stroke="{"#d9822b" if prop else "#1f6f8b"}" '
                       'stroke-width="0.8"/>')
    for r in layers.get("flowpaths", []):
        for ln in _lines(r.get("_geom")):
            out.append(f'<path d="{path(ln)}" fill="none" stroke="#3a7bd5" stroke-width="1.1"/>')
    for r in layers.get("road_alignment", []):
        for ln in _lines(r.get("_geom")):
            out.append(f'<path d="{path(ln)}" fill="none" stroke="#1f2a44" stroke-width="3"/>')
    for r in layers.get("sag_points", []):
        for ln in _lines(r.get("_geom")):
            x, y = xy(ln[0]).split(",")
            out.append(f'<rect x="{float(x) - 4:.1f}" y="{float(y) - 4:.1f}" width="8" height="8" '
                       'fill="#c0392b" transform="rotate(45 ' + f'{x} {y})"/>')
    for r in layers.get("crossings", []):
        for ln in _lines(r.get("_geom")):
            x, y = (float(v) for v in xy(ln[0]).split(","))
            prop = r.get("status") == "proposed"
            out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" fill="{"#d9822b" if prop else "#1f6f8b"}" '
                       'stroke="#fff" stroke-width="1.5"/>')
            out.append(f'<text x="{x + 7:.1f}" y="{y - 6:.1f}" font-size="11" fill="#1d2330" '
                       f'paint-order="stroke" stroke="#fff" stroke-width="3">{_e(r.get("outlet_uid"))}</text>')
    # north arrow and scale bar
    out.append(f'<g transform="translate({width - 34},{40})"><path d="M0,-18 L7,6 L0,1 L-7,6 Z" '
               'fill="#1f2a44"/><text x="0" y="20" font-size="11" text-anchor="middle" '
               'fill="#1f2a44">N</text></g>')
    target = span / 5
    mag = 10 ** math.floor(math.log10(target))
    bar = max(m * mag for m in (1, 2, 5) if m * mag <= target)
    lab = f"{bar / 1000:g} km" if bar >= 1000 else f"{bar:g} m"
    out.append(f'<g transform="translate({pad},{height - 16})"><rect width="{bar * s:.1f}" '
               'height="4" fill="#1f2a44"/>'
               f'<text x="{bar * s + 6:.1f}" y="5" font-size="11" fill="#1f2a44">{lab}</text></g>')
    out.append("</svg>")
    return "\n".join(out)


def _table(cols, rows, num=()):
    h = "".join(f'<th class="{"n" if c in num else ""}">{_e(c)}</th>' for c in cols)
    body = []
    for r in rows:
        body.append("<tr>" + "".join(f'<td class="{"n" if c in num else ""}">{r.get(c, "")}</td>'
                                     for c in cols) + "</tr>")
    return f"<table><thead><tr>{h}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def build_report(package_path, title=None, extra=None):
    """-> HTML text. extra (all optional): warnings [str], dem_audit [str],
    flat_check {uid: {flat_sensitivity_pct, flat_sensitive, area_barnes_km2,
    area_toward_km2}}, settings_path, outputs {label: path}, run_seconds."""
    extra = dict(extra or {})
    tables = dict(gpkg.list_tables(package_path))
    md = {r["key"]: r["value"] for r in gpkg.read_table(package_path, "qeht_run_metadata")}
    flat = extra.get("flat_check") or {}
    cols, rows = characteristics_rows(package_path, flat)
    n_prop = sum(1 for r in rows if r.get("status") == "proposed")
    areas = [r["area_km2"] for r in rows if _has(r.get("area_km2"))]
    title = title or "QEHT run report"
    gen = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    P = []
    P.append(f"<h1>{_e(title)}</h1><p class='sub'>Generated {gen} from <code>"
             f"{_e(package_path.replace(chr(92), '/').split('/')[-1])}</code> · QEHT "
             f"{_e(md.get('qeht_version'))} · run {_e(md.get('run_utc'))}</p>")
    P.append("<div class='note'>Everything on this page is read from the package "
             "(<code>qeht_run_metadata</code> and its layers); nothing is recalculated. "
             "Print from a browser (Save as PDF) for an appendix.</div>")
    cards = [("Crossings", f"{len(rows) - n_prop}" + (f" + {n_prop} proposed" if n_prop else "")),
             ("Largest catchment", f"{max(areas):,.3g} km²" if areas else "–"),
             ("Smallest catchment", f"{min(areas):,.3g} km²" if areas else "–"),
             ("Cell size", f"{_e(md.get('cell_size_m'))} m"),
             ("CRS", f"EPSG:{_e(md.get('crs_epsg'))}"),
             ("Flat method", _e(md.get("flat_method"))),
             ("Stream threshold", f"{_fmt(_f(md.get('stream_threshold_km2')))} km²")]
    if flat:
        cards.append(("Flat-sensitive outlets", str(sum(1 for v in flat.values()
                                                        if v.get("flat_sensitive")))))
    P.append("<div class='cards'>" + "".join(
        f"<div class='card'><div class='k'>{k}</div><div class='v'>{v}</div></div>"
        for k, v in cards) + "</div>")

    svg = plan_svg(package_path, tables, background=extra.get("background"))
    if svg:
        P.append("<h2>Plan</h2><figure>" + svg + "<figcaption>Schematic plan from the package "
                 "geometry: catchments (blue; proposed orange), longest flow paths, road and "
                 "crossings" + (" over the terrain (hillshade of the DEM)" if extra.get("background")
                                else "") + ". Drawn to scale.</figcaption></figure>")

    # crossing schedule
    show = [c for c in ("outlet_uid", "status", "chainage_m", "area_km2", "lfp_length_m",
                        "lfp_slope_1085", "catch_slope_horn", "cn_export", "rain_zone", "map_mm",
                        "flat_sensitive") if c in cols]
    num = {"chainage_m", "area_km2", "lfp_length_m", "lfp_slope_1085", "catch_slope_horn",
           "cn_export", "map_mm"}
    srows = []
    for r in rows:
        d = {c: _fmt(r.get(c)) for c in show}
        d["status"] = (f"<span class='tag {'p' if r.get('status') == 'proposed' else 'e'}'>"
                       f"{_e(r.get('status') or 'existing')}</span>")
        if r.get("flat_sensitive"):
            d["flat_sensitive"] = "<span class='flag'>check</span>"
        elif "flat_sensitive" in show:
            d["flat_sensitive"] = "–"
        srows.append(d)
    P.append("<h2>Crossing schedule</h2>" + _table(show, srows, num))
    P.append("<p class='sub'>The full table (every characteristic) is "
             "<code>tables/catchment_characteristics.csv</code> for a pipeline run; field "
             "meanings and units are in the package's <code>qeht_field_dictionary</code>.</p>")

    # DEM checks and flat-method check
    P.append("<h2>DEM checks</h2>")
    audit = extra.get("dem_audit")
    if audit is None:
        P.append("<p>Resampling audit: not run in this session (Fill and D8 report it in their "
                 "logs).</p>")
    elif audit:
        P.append("<ul class='warn'>" + "".join(f"<li class='flag'>{_e(a)}</li>" for a in audit) + "</ul>")
    else:
        P.append("<p>Resampling audit: no repeated rows or columns, no NoData problems.</p>")
    if flat:
        frows = []
        flagged = {u: v for u, v in flat.items() if v.get("flat_sensitive")}
        for uid, v in sorted(flagged.items(), key=lambda kv: -(kv[1].get("flat_sensitivity_pct") or 0)):
            frows.append({"outlet_uid": _e(uid), "area_barnes_km2": _fmt(v.get("area_barnes_km2")),
                          "area_toward_km2": _fmt(v.get("area_toward_km2")),
                          "flat_sensitivity_pct": _fmt(v.get("flat_sensitivity_pct")),
                          "flag": "<span class='flag'>verify</span>" if v.get("flat_sensitive") else "–"})
        P.append("<h3>Flat-method check</h3><p class='sub'>Contributing area under both flat "
                 "methods (largest within one cell of the outlet). A flagged crossing's area "
                 "depends on how flats are routed: verify it against mapped drainage or on site."
                 f"</p><p>{len(flagged)} of {len(flat)} crossing(s) flagged"
                 + (" (the others agree within the tolerance; values in the characteristics "
                    "table)." if flagged else ".") + "</p>"
                 + (_table(["outlet_uid", "area_barnes_km2", "area_toward_km2",
                            "flat_sensitivity_pct", "flag"], frows,
                           {"area_barnes_km2", "area_toward_km2", "flat_sensitivity_pct"})
                    if flagged else ""))

    # time of concentration (F10)
    from .characteristics import tc_rows
    tcols, trows = tc_rows(package_path)
    if any(_has(r.get("tc_kirpich_min")) for r in trows):
        meth = [("tc_kirpich_min", "Kirpich"), ("tc_kerby_kirpich_min", "Kerby + Kirpich"),
                ("tc_scs_lag_min", "SCS lag"), ("tc_tr55_min", "TR-55"),
                ("tc_bransby_williams_min", "Bransby-Williams")]
        hdr = ["outlet_uid", "area_km2"] + [k for k, _ in meth]
        body = []
        for r in trows:
            d = {"outlet_uid": _e(r["outlet_uid"]), "area_km2": _fmt(r.get("area_km2"))}
            for k, _ in meth:
                fl = r.get(k.replace("_min", "_flag")) or ""
                d[k] = _fmt(r.get(k)) + (" <span class='flag'>*</span>" if fl.startswith("outside") else "")
            body.append(d)
        tj = _json(md, "tc_params_json") or {}
        P.append("<h2>Time of concentration (min)</h2><p class='sub'>Every method side by side; "
                 "QEHT does not choose one. * = outside the method's published calibration range "
                 "(ranges and sources in <code>tc_params_json</code>; inputs per crossing in "
                 "<code>tc_basis_json</code>). TR-55 P2: " + _e(tj.get("p2_source", "not given"))
                 + ". Sheet n and Kerby N: "
                 + ("PROXY lookups from land cover" if str(tj.get("sheet_n_lookup", "")).startswith("PROXY")
                    else "user lookups") + ".</p>" + _table(hdr, body, {"area_km2"} | {k for k, _ in meth})
                 + "<p class='sub'>Full table: <code>tables/time_of_concentration.csv</code>.</p>")

    # channel sections of low quality (F5)
    if "crossings" in tables:
        xs = [r for r in gpkg.read_table(package_path, "crossings", with_geometry=False)
              if r.get("xs_quality") == "low"]
        if xs:
            P.append("<h3>Channel sections of low quality</h3><p class='sub'>The channel "
                     "cross-section downstream of these crossings is below the grid resolution "
                     "or has no clear bank; use survey for the tailwater and waterway.</p>"
                     + _table(["outlet_uid", "xs_bankfull_w_m", "xs_w_1p0_m", "xs_note"],
                              [{k: _fmt(r.get(k)) for k in ("outlet_uid", "xs_bankfull_w_m",
                                                           "xs_w_1p0_m", "xs_note")} for r in xs],
                              {"xs_bankfull_w_m", "xs_w_1p0_m"}))

    # coverage
    if "coverage_check" in tables:
        cov = gpkg.read_table(package_path, "coverage_check", with_geometry=False)
        P.append("<h2>Drainage coverage along the road</h2>")
        if cov:
            P.append(_table(["issue", "chainage_m", "chainage_to_m", "area_km2", "uid", "note"],
                            [{k: _fmt(r.get(k)) for k in ("issue", "chainage_m", "chainage_to_m",
                                                         "area_km2", "uid", "note")} for r in cov],
                            {"chainage_m", "chainage_to_m", "area_km2"}))
        else:
            P.append("<p>No missing crossings, sags or flat stretches found.</p>")

    # provenance
    prov = [{"k": f"<span>{_e(lbl)}</span>", "v": _e(md.get(key))}
            for key, lbl in PROVENANCE if _has(md.get(key))]
    ro, rf, er = _json(md, "runoff_json"), _json(md, "rainfall_json"), _json(md, "erosion_json")
    if ro:
        prov.append({"k": "Curve numbers", "v": _e(ro.get("cn_lookup_id")) + (
            " <span class='tag w'>PROXY</span>" if ro.get("cn_proxy") else "")
            + f" · AMC {_e(ro.get('cn_amc'))} · land cover {_e(ro.get('lc_dataset'))}"})
        if ro.get("rc_lookup_id"):
            prov.append({"k": "Rational C", "v": _e(ro.get("rc_lookup_id"))})
    if rf:
        v = []
        if rf.get("zone_source"):
            v.append(f"zones {_e(rf['zone_source'])} ({_e(', '.join(rf.get('zones') or []))}); "
                     f"QA below {_fmt(rf.get('zone_qa_threshold_pct'))} %")
        if rf.get("map_dataset"):
            v.append(f"rainfall {_e(rf['map_dataset'])}")
        if rf.get("r_is_estimate"):
            v.append(f"<span class='tag w'>R ESTIMATE</span> {_e(rf.get('r_relation_text'))}")
        prov.append({"k": "Rainfall", "v": "<br>".join(v)})
    if er:
        fs = "; ".join(f"{_e(f.get('factor'))}: {_e(f.get('basis'))}"
                       + (f" ({_e(f.get('source'))})" if f.get("source") else "")
                       for f in er.get("factors", []))
        prov.append({"k": "Erosion", "v": f"{_e(er.get('mode'))}<br>{fs}"})
    P.append("<h2>Inputs and provenance</h2><table><tbody>" + "".join(
        f"<tr><td class='kv'>{r['k']}</td><td>{r['v']}</td></tr>" for r in prov) + "</tbody></table>")

    # warnings
    warns = extra.get("warnings") or []
    P.append("<h2>Warnings</h2>" + ("<ul class='warn'>" + "".join(f"<li>{_e(w)}</li>" for w in warns)
                                     + "</ul>" if warns else "<p>None.</p>"))

    # lookups
    if ro and ro.get("cn_lookup"):
        lr = [{"class": _e(k), "A": _fmt(v[0]), "B": _fmt(v[1]), "C": _fmt(v[2]), "D": _fmt(v[3])}
              for k, v in sorted(ro["cn_lookup"].items(), key=lambda kv: int(kv[0]))
              if isinstance(v, (list, tuple)) and len(v) == 4]
        P.append("<h2>Curve number lookup used</h2>" + _table(["class", "A", "B", "C", "D"], lr,
                                                              {"A", "B", "C", "D"}))
    outs = extra.get("outputs") or {}
    if outs:
        P.append("<h2>Output files</h2><table><tbody>" + "".join(
            f"<tr><td class='kv'>{_e(k)}</td><td><code>{_e(v)}</code></td></tr>"
            for k, v in outs.items()) + "</tbody></table>")
    if extra.get("settings_path"):
        P.append(f"<p class='sub'>Settings for an exact re-run: <code>{_e(extra['settings_path'])}"
                 "</code> (Run hydrology pipeline → Re-run from settings).</p>")
    return (f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' "
            f"content='width=device-width,initial-scale=1'><title>{_e(title)}</title><style>{CSS}"
            f"</style></head><body><main>{''.join(P)}</main></body></html>")


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def write_report(package_path, out_html, title=None, extra=None):
    text = build_report(package_path, title, extra)
    with open(out_html, "w", encoding="utf-8") as f:
        f.write(text)
    return out_html
