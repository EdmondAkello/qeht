# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""KMZ export of a design hydrology package (F11, v0.23).

KML 2.2 in a zip, written in pure Python (zipfile and string templates,
every text escaped). The caller supplies `to_wgs84(points) -> points`, a
transformation from the package CRS to longitude / latitude (OSR in the
processing layer), so this module never reprojects.

Folders: run record, crossings (existing blue, proposed orange; one card
per crossing with location, catchment, flow path, time of concentration,
runoff, channel and erosion values), catchments (translucent, labelled),
longest flow paths, road with kilometre posts, coverage findings, sag
points and flat stretches, channel sections, side-drain siltation and,
optionally, the relief quicklook as a ground overlay. Icons are drawn here
(NumPy, PNG writer of quicklooks); the legend and title block are drawn
with matplotlib when it is installed and left out otherwise.

No QGIS imports, no GDAL.
"""

import datetime
import html
import io
import json
import math
import os
import zipfile

import numpy as np

from ..interop import gpkg
from .quicklooks import png_bytes

# palette (hex RGB); KML colours are aabbggrr
C = {"existing": "#2b83ba", "proposed": "#f28e2b", "navy": "#1f3b57", "text": "#1d2733",
     "label": "#51606f", "rule": "#c9d3dd", "zebra": "#f3f6f9", "lfp": "#ae017e",
     "lfp_hi": "#ff3fbf", "road": "#252525", "road_case": "#fff7bc", "finding": "#d7301f",
     "sag": "#6a51a3", "flat": "#1b9e77", "xs": "#00a6ca", "kmpost": "#3b4b5c", "ok": "#1a9641",
     "warn": "#d95f02", "bad": "#d7191c", "grey": "#7f8c99", "mapped": "#4292c6"}
STI_COLOURS = {"low": "#1a9641", "moderate": "#a6d96a", "high": "#fdae61", "very high": "#d7191c"}
FONT = "font-family:Arial,Helvetica,sans-serif;font-size:12px;color:#1d2733;width:430px"


def kml_colour(hex_rgb, alpha=255):
    h = hex_rgb.lstrip("#")
    return f"{alpha:02x}{h[4:6]}{h[2:4]}{h[0:2]}"


def x(text):
    """XML text escape."""
    return html.escape("" if text is None else str(text), quote=True)


def cdata(text):
    return "<![CDATA[" + str(text).replace("]]>", "]]]]><![CDATA[>") + "]]>"


def _has(v):
    return v is not None and v != "" and not (isinstance(v, float) and not math.isfinite(v))


def fmt(v, nd=3, unit=""):
    if not _has(v):
        return None
    if isinstance(v, (int, np.integer)) and not isinstance(v, bool):
        s = f"{v:,d}"
    elif isinstance(v, (float, np.floating)):
        a = abs(v)
        s = (f"{v:,.0f}" if a >= 1000 else f"{v:,.1f}" if a >= 100 else
             f"{v:.{nd}f}".rstrip("0").rstrip(".") if a >= 1 else f"{v:.{nd}g}")
    else:
        s = str(v)
    return s + (f" {unit}" if unit else "")


def chainage(ch):
    if not _has(ch):
        return None
    km = int(ch // 1000)
    return f"{km}+{ch - 1000 * km:06.2f}"


# -- HTML cards ---------------------------------------------------------------

def badge(text, colour):
    return (f"<span style=\"background:{colour};color:#fff;padding:1px 6px;border-radius:8px;"
            f"font-size:11px\">{x(text)}</span> ")


def rows_table(rows):
    rows = [(k, v) for k, v in rows if v is not None]
    if not rows:
        return ""
    out = ["<table style=\"border-collapse:collapse;width:100%\">"]
    for i, (k, v) in enumerate(rows):
        bg = C["zebra"] if i % 2 == 0 else "#ffffff"
        out.append(f"<tr style=\"background:{bg}\"><td style=\"padding:2px 5px;color:{C['label']};"
                   f"width:48%;font-size:12px\">{x(k)}</td><td style=\"padding:2px 5px;font-weight:bold;"
                   f"font-size:12px\">"
                   f"{x(v)}</td></tr>")
    return "".join(out) + "</table>"


def grid_table(header, rows, highlight=None):
    th = "".join(f"<th style=\"background:{C['navy']};color:#fff;padding:3px\">{x(h)}</th>"
                 for h in header)
    body = []
    for i, r in enumerate(rows):
        bg = "#dbe9f6" if highlight and highlight(r) else (C["zebra"] if i % 2 == 0 else "#fff")
        body.append(f"<tr style=\"background:{bg}\">" + "".join(
            f"<td style=\"padding:2px 4px;border-bottom:1px solid #e3e8ee\">{x(c if c is not None else '–')}</td>"
            for c in r) + "</tr>")
    return (f"<table style=\"border-collapse:collapse;width:100%;text-align:center;font-size:11px\">"
            f"<tr>{th}</tr>{''.join(body)}</table>")


def note(text):
    if not text:
        return ""
    return (f"<div style=\"font-size:10px;color:{C['label']};margin:3px 2px;line-height:1.3\">"
            f"{x(text)}</div>")


def section(title, inner):
    if not inner:
        return ""
    return (f"<div style=\"margin:8px 0 3px 0;font-weight:bold;color:{C['navy']};border-bottom:1px "
            f"solid {C['rule']};font-size:12px\">{x(title)}</div>{inner}")


def card(colour, title, subtitle="", badges="", body="", foot=""):
    return (f"<div style=\"{FONT}\"><div style=\"background:{colour};color:#fff;padding:7px 9px;"
            f"border-radius:4px 4px 0 0\"><div style=\"font-size:15px;font-weight:bold\">{x(title)}"
            f"</div>" + (f"<div style=\"font-size:11px;opacity:.92\">{x(subtitle)}</div>" if subtitle
                         else "") + "</div><div style=\"padding:4px 2px\">"
            + (f"<div style=\"margin:4px 0\">{badges}</div>" if badges else "") + body
            + (f"<div style=\"margin-top:6px;font-size:10px;color:{C['label']}\">{x(foot)}</div>"
               if foot else "") + "</div></div>")


# -- icons --------------------------------------------------------------------

def _rgb(hex_rgb):
    h = hex_rgb.lstrip("#")
    return [int(h[i:i + 2], 16) for i in (0, 2, 4)]


def icon_png(shape, fill, size=64, ss=4):
    """Anti-aliased icon (circle, square, diamond, triangle) with a dark outline and a
    white halo, supersampled ss x."""
    n = size * ss
    yy, xx = (np.mgrid[0:n, 0:n] + 0.5) / n * 2.0 - 1.0
    if shape == "circle":
        d = np.hypot(xx, yy)
    elif shape == "square":
        d = np.maximum(np.abs(xx), np.abs(yy)) * 1.08
    elif shape == "diamond":
        d = (np.abs(xx) + np.abs(yy)) * 0.80
    else:                                                   # triangle, point down (rows run down)
        yy = yy + 0.2                                       # centroid a little above the middle
        d = 2.0 * np.maximum(-yy, 0.866 * np.abs(xx) + 0.5 * yy)
    halo, edge, body = d <= 0.92, d <= 0.80, d <= 0.66
    img = np.zeros((n, n, 4), float)
    img[halo] = [255, 255, 255, 235]
    img[edge] = [29, 39, 51, 255]
    img[body] = _rgb(fill) + [255]
    img = img.reshape(size, ss, size, ss, 4).mean(axis=(1, 3))
    return png_bytes(np.clip(img, 0, 255).astype(np.uint8))


# -- geometry helpers -----------------------------------------------------------

def coords_text(pts):
    return " ".join(f"{lon:.7f},{lat:.7f},0" for lon, lat in pts)


def ring_centroid(ring):
    a = cx = cy = 0.0
    for (x0, y0), (x1, y1) in zip(ring, ring[1:] + ring[:1]):
        cr = x0 * y1 - x1 * y0
        a += cr
        cx += (x0 + x1) * cr
        cy += (y0 + y1) * cr
    if abs(a) < 1e-12:
        return tuple(np.mean(np.asarray(ring), axis=0))
    return cx / (3 * a), cy / (3 * a)


def _point_in_ring(p, ring):
    xp, yp = p
    inside = False
    for (x0, y0), (x1, y1) in zip(ring, ring[1:] + ring[:1]):
        if (y0 > yp) != (y1 > yp) and xp < (x1 - x0) * (yp - y0) / (y1 - y0) + x0:
            inside = not inside
    return inside


def label_point(polys):
    """A point inside the largest polygon (centroid, or the nearest vertex mean inside)."""
    best = max(polys, key=lambda p: abs(sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1)
                                             in zip(p[0], p[0][1:] + p[0][:1]))))
    ring = list(best[0])
    c = ring_centroid(ring)
    if _point_in_ring(c, ring):
        return c
    # fall back: midpoint of the horizontal chord through the centroid's row
    ys = sorted({y for _, y in ring})
    yc = min(ys, key=lambda y: abs(y - c[1]))
    xs = sorted(x0 + (yc - y0) * (x1 - x0) / (y1 - y0) for (x0, y0), (x1, y1)
                in zip(ring, ring[1:] + ring[:1]) if (y0 > yc) != (y1 > yc))
    if len(xs) >= 2:
        return 0.5 * (xs[0] + xs[1]), yc
    return ring[0]


def along(line, d):
    """Point at distance d along a polyline."""
    acc = 0.0
    for (x0, y0), (x1, y1) in zip(line[:-1], line[1:]):
        seg = math.hypot(x1 - x0, y1 - y0)
        if acc + seg >= d and seg > 0:
            t = (d - acc) / seg
            return x0 + t * (x1 - x0), y0 + t * (y1 - y0)
        acc += seg
    return line[-1]


# -- styles -------------------------------------------------------------------

def icon_style(sid, href, scale=1.0, label=True):
    def one(suffix, sc, lc, ls):
        return (f"<Style id=\"{sid}_{suffix}\"><IconStyle><scale>{sc}</scale><Icon><href>{href}"
                f"</href></Icon><hotSpot x=\"0.5\" y=\"0.5\" xunits=\"fraction\" yunits=\"fraction\"/>"
                f"</IconStyle><LabelStyle><color>{lc}</color><scale>{ls if label else 0}</scale>"
                f"</LabelStyle><BalloonStyle><text>$[description]</text></BalloonStyle></Style>")
    return (one("n", 0.9 * scale, "ffffffff", 0.85) + one("h", 1.25 * scale, "ff00ffff", 1.05)
            + f"<StyleMap id=\"{sid}\"><Pair><key>normal</key><styleUrl>#{sid}_n</styleUrl></Pair>"
            f"<Pair><key>highlight</key><styleUrl>#{sid}_h</styleUrl></Pair></StyleMap>")


def line_style(sid, hex_rgb, width, hi_hex=None, alpha=255, fill=None, fill_alpha=56,
               label_scale=0.0):
    def one(suffix, col, w, fa):
        poly = (f"<PolyStyle><color>{kml_colour(fill, fa)}</color></PolyStyle>" if fill else "")
        icon = ("<IconStyle><scale>0</scale></IconStyle>" if fill else "")
        return (f"<Style id=\"{sid}_{suffix}\">{icon}<LabelStyle><color>{kml_colour(col)}</color>"
                f"<scale>{label_scale}</scale></LabelStyle><LineStyle><color>{kml_colour(col, alpha)}"
                f"</color><width>{w}</width></LineStyle>{poly}<BalloonStyle><text>$[description]"
                f"</text></BalloonStyle></Style>")
    return (one("n", hex_rgb, width, fill_alpha) + one("h", hi_hex or hex_rgb, width * 1.6,
                                                         min(fill_alpha * 2, 255))
            + f"<StyleMap id=\"{sid}\"><Pair><key>normal</key><styleUrl>#{sid}_n</styleUrl></Pair>"
            f"<Pair><key>highlight</key><styleUrl>#{sid}_h</styleUrl></Pair></StyleMap>")


# -- placemarks -----------------------------------------------------------------

def placemark(name, style, geometry, description="", snippet="", visible=True):
    return (f"<Placemark><name>{x(name)}</name><visibility>{int(visible)}</visibility>"
            f"<Snippet maxLines=\"{1 if snippet else 0}\">{x(snippet)}</Snippet>"
            + (f"<description>{cdata(description)}</description>" if description else "")
            + f"<styleUrl>#{style}</styleUrl>{geometry}</Placemark>")


def folder(name, items, open_=False, visible=True):
    return (f"<Folder><name>{x(name)}</name><visibility>{int(visible)}</visibility>"
            f"<open>{int(open_)}</open>{''.join(items)}</Folder>")


def point_geom(p):
    return f"<Point><coordinates>{p[0]:.7f},{p[1]:.7f},0</coordinates></Point>"


def line_geom(pts):
    return f"<LineString><tessellate>1</tessellate><coordinates>{coords_text(pts)}</coordinates></LineString>"


def polygon_geom(rings):
    outer, holes = rings[0], rings[1:]
    s = (f"<Polygon><tessellate>1</tessellate><outerBoundaryIs><LinearRing><coordinates>"
         f"{coords_text(outer)}</coordinates></LinearRing></outerBoundaryIs>")
    for h in holes:
        s += (f"<innerBoundaryIs><LinearRing><coordinates>{coords_text(h)}</coordinates>"
              f"</LinearRing></innerBoundaryIs>")
    return s + "</Polygon>"


# -- content --------------------------------------------------------------------

TC_METHODS = [("tc_kirpich", "Kirpich"), ("tc_kerby_kirpich", "Kerby + Kirpich"),
              ("tc_scs_lag", "SCS lag"), ("tc_tr55", "TR-55"),
              ("tc_bransby_williams", "Bransby-Williams")]


def crossing_card(c, ca, fp, md):
    status = c.get("status") or "existing"
    col = C["proposed"] if status == "proposed" else C["existing"]
    b = badge("Proposed - " + (c.get("proposed_reason") or "").replace("_", " ")
              if status == "proposed" else "Existing", col)
    if c.get("xs_quality"):
        b += badge(f"Channel section: {c['xs_quality']}",
                   {"high": C["ok"], "medium": C["warn"], "low": C["bad"]}[c["xs_quality"]])
    if c.get("map_agrees") is not None:
        b += badge("On the mapped waterway" if c["map_agrees"] == 1 else
                   f"{fmt(c.get('map_river_dist_m'), 0)} m from the mapped waterway",
                   C["ok"] if c["map_agrees"] == 1 else C["warn"])
    if c.get("ero_impact"):
        b += badge(f"Erosion impact: {c['ero_impact']}", C["grey"])
    if c.get("ero_dep_flag"):
        b += badge(c["ero_dep_flag"], C["grey"])
    epsg = md.get("crs_epsg") or ""
    loc = rows_table([
        ("Chainage", chainage(c.get("chainage_m"))),
        (f"Easting / Northing (EPSG:{epsg})",
         f"{c['outlet_x']:,.2f} / {c['outlet_y']:,.2f}" if _has(c.get("outlet_x")) else None),
        ("Snapped by", fmt(c.get("snap_dist_m"), 1, "m")),
        ("Strahler order", fmt(c.get("stream_order"))),
        ("Source ID", c.get("source_id")),
        ("Nearest mapped waterway", f"{c.get('map_river_name') or 'unnamed'} "
         f"({fmt(c.get('map_river_dist_m'), 1, 'm')})" if _has(c.get("map_river_dist_m")) else None),
        ("Nearest existing crossing", f"{c['nearest_uid']} ({fmt(c.get('nearest_m'), 1, 'm')})"
         if c.get("nearest_uid") else None)])
    cat = rows_table([
        ("Area", fmt(ca.get("area_km2"), 4, "km²")),
        ("Mean slope (Horn)", fmt(ca.get("catch_slope_horn"), 4, "m/m")),
        ("Elevation max / min / mean",
         " / ".join(f"{ca[k]:,.0f}" for k in ("elev_max_m", "elev_min_m", "elev_mean_m")) + " m"
         if all(_has(ca.get(k)) for k in ("elev_max_m", "elev_min_m", "elev_mean_m")) else None),
        ("Relief", fmt(ca.get("relief_m"), 1, "m")),
        ("Soil / HSG", " / ".join(str(v) for v in (ca.get("soil_texture"), ca.get("soil_hsg"))
                                  if _has(v)) or None),
        ("CN (AMC II) / exported", f"{fmt(ca.get('cn_ii'), 1)} / {fmt(ca.get('cn_export'), 1)} "
         f"(AMC {ca.get('cn_amc')})" if _has(ca.get("cn_ii")) else None),
        ("Rational C", fmt(ca.get("rational_c"), 2)),
        ("Rain zone", f"{ca['rain_zone']} ({fmt(ca.get('rain_zone_pct'), 1, '%')})"
         if ca.get("rain_zone") else None),
        ("Mean annual rainfall", fmt(ca.get("map_mm"), 0, "mm/yr"))])
    lfp = rows_table([
        ("Longest flow path", f"{fmt(fp.get('lfp_length_m') / 1000.0, 3)} km"
         if _has(fp.get("lfp_length_m")) else None),
        ("Drop", fmt(fp.get("lfp_drop_m"), 1, "m")),
        ("Slope 10-85 / drop-length", f"{fmt(fp.get('lfp_slope_1085'), 4)} / "
         f"{fmt(fp.get('lfp_slope'), 4)} m/m" if _has(fp.get("lfp_slope_1085")) else None),
        ("Overland / channel length", f"{fmt(fp.get('lfp_overland_m'), 0)} / "
         f"{fmt(fp.get('lfp_channel_m'), 0)} m" if _has(fp.get("lfp_overland_m")) else None),
        ("Channel slope up / down", f"{fmt(c.get('ch_slope_us'), 4)} / {fmt(c.get('ch_slope_ds'), 4)} m/m"
         if _has(c.get("ch_slope_us")) else None)])
    tc_rows = []
    for k, name in TC_METHODS:
        v = c.get(k + "_min")
        if _has(v):
            fl = c.get(k + "_flag") or ""
            tc_rows.append([name, fmt(v, 1), fl[len("outside: "):] if fl.startswith("outside")
                            else ("rural catchments" if fl == "rural catchments" else "within range")])
    tc = grid_table(["Method", "Tc (min)", "Validity"], tc_rows) if tc_rows else ""
    if tc:
        tc += note("Every method reported; QEHT does not choose one."
                   + (f" Not computed: {c['tc_note']}." if c.get("tc_note") else ""))
    xs = rows_table([
        ("Bed level", fmt(c.get("xs_bed_m"), 2, "m")),
        ("Bank-full width / depth", f"{fmt(c.get('xs_bankfull_w_m'), 1)} m / "
         f"{fmt(c.get('xs_bankfull_d_m'), 2)} m" if _has(c.get("xs_bankfull_w_m")) else None),
        ("Width at bed + 0.5 / 1 / 2 m", " / ".join(fmt(c.get(k), 1) or "–" for k in
                                                    ("xs_w_0p5_m", "xs_w_1p0_m", "xs_w_2p0_m")) + " m"
         if _has(c.get("xs_w_1p0_m")) else None),
        ("Floodplain width along the road (bed + 1 m)", fmt(c.get("fp_w_1p0_m"), 0, "m"))]) \
        + note(c.get("xs_note"))
    ero = rows_table([
        ("Impact score", fmt(c.get("ero_composite"), 2)),
        ("Mean soil loss", fmt(ca.get("ero_a_mean_tha"), 1, "t/ha/yr")),
        ("Sediment volume", fmt(ca.get("ero_sy_m3_yr"), 0, "m³/yr")),
        ("Deposition ratio", fmt(c.get("ero_dep_ratio"), 2))])
    body = (section("Location", loc) + section("Contributing catchment", cat)
            + section("Flow path and channel slopes", lfp) + section("Time of concentration", tc)
            + section("Channel downstream (indicative)", xs) + section("Erosion", ero))
    sub = " · ".join(v for v in (f"Ch {chainage(c.get('chainage_m'))}" if _has(c.get("chainage_m")) else "",
                                  f"{fmt(ca.get('area_km2'), 3)} km²" if _has(ca.get("area_km2")) else "")
                     if v)
    return card(col, f"Crossing {c['outlet_uid']}", sub, b, body,
                foot=f"QEHT {md.get('qeht_version', '')} · terrain values from the DEM; check "
                     "against survey")


def run_card(md, n_ex, n_prop, title):
    rows = [("Crossings", f"{n_ex} existing" + (f" + {n_prop} proposed" if n_prop else "")),
            ("DEM", md.get("dem_source")), ("Preparation", md.get("dem_prep") or None),
            ("CRS", f"EPSG:{md.get('crs_epsg')} {md.get('crs_name') or ''}".strip()),
            ("Cell size", f"{md.get('cell_size_m')} m"), ("Flat method", md.get("flat_method")),
            ("Conditioning", md.get("conditioning")),
            ("Stream threshold", f"{md.get('stream_threshold_km2')} km²"),
            ("Catchments", md.get("catchment_mode")), ("QEHT", md.get("qeht_version")),
            ("Run (UTC)", md.get("run_utc"))]
    try:
        ro = json.loads(md.get("runoff_json") or "{}")
    except ValueError:
        ro = {}
    if ro.get("cn_proxy"):
        rows.append(("Curve numbers", "PROXY lookup (TR-55 matched to WorldCover) - check before "
                                      "design use"))
    if "ESTIMATE" in (md.get("rainfall_json") or "").upper():
        rows.append(("RUSLE R", "ESTIMATE from mean annual rainfall"))
    return card(C["navy"], title, "QEHT design hydrology package - GIS summary layers", "",
                rows_table(rows), foot="Everything here is read from the package; nothing is "
                                       "recalculated.")


def _layer(pkg, name, tabs, geom=True):
    return gpkg.read_table(pkg, name, with_geometry=geom) if name in tabs else []


def build_kml(pkg, to_wgs84, title=None, relief=None):
    """-> (kml text, {zip path: bytes})."""
    tabs = dict(gpkg.list_tables(pkg))
    md = {r["key"]: r["value"] for r in gpkg.read_table(pkg, "qeht_run_metadata")}
    cr = _layer(pkg, "crossings", tabs)
    ca = {r["outlet_uid"]: r for r in _layer(pkg, "catchments", tabs)}
    fp = {r["outlet_uid"]: r for r in _layer(pkg, "flowpaths", tabs)}
    title = title or "QEHT design hydrology"
    files = {}
    for st, col in (("existing", C["existing"]), ("proposed", C["proposed"])):
        files[f"files/circle_{st}.png"] = icon_png("circle", col)
    files["files/triangle_finding.png"] = icon_png("triangle", C["finding"])
    files["files/diamond_sag.png"] = icon_png("diamond", C["sag"])
    files["files/square_kmpost.png"] = icon_png("square", C["kmpost"], size=40)
    styles = [icon_style("x_existing", "files/circle_existing.png"),
              icon_style("x_proposed", "files/circle_proposed.png", scale=1.05),
              icon_style("finding", "files/triangle_finding.png", 0.9),
              icon_style("sag", "files/diamond_sag.png", 0.85),
              icon_style("kmpost", "files/square_kmpost.png", 0.55),
              line_style("cat_existing", C["existing"], 2.0, fill=C["existing"], label_scale=0.85),
              line_style("cat_proposed", C["proposed"], 2.0, fill=C["proposed"], label_scale=0.85),
              line_style("lfp", C["lfp"], 2.6, hi_hex=C["lfp_hi"]),
              line_style("road_case", C["road_case"], 7.0),
              line_style("road", C["road"], 3.2),
              line_style("flat", C["flat"], 5.0, alpha=200),
              line_style("xs", C["xs"], 1.6),
              line_style("mapped", C["mapped"], 2.4, alpha=210),
              line_style("divergence", C["bad"], 3.4)]
    styles += [line_style("sti_" + k.replace(" ", "_"), v, 4.0, alpha=230)
               for k, v in STI_COLOURS.items()]
    # everything to WGS 84 in one call per layer
    allpts = []
    n_ex = sum(1 for r in cr if (r.get("status") or "existing") != "proposed")
    n_prop = len(cr) - n_ex
    folders = []
    # run record
    folders.append(placemark("Run record (open me)", "kmpost", "", run_card(md, n_ex, n_prop, title),
                             snippet=f"{len(cr)} crossings · {md.get('dem_source', '')}"))
    # crossings
    items = []
    pts = to_wgs84([r["_geom"] for r in cr]) if cr else []
    allpts += pts
    for r, p in sorted(zip(cr, pts), key=lambda t: (t[0].get("chainage_m") is None,
                                                    t[0].get("chainage_m") or 0, t[0]["outlet_uid"])):
        a, f = ca.get(r["outlet_uid"], {}), fp.get(r["outlet_uid"], {})
        snip = " · ".join(v for v in (
            f"Ch {chainage(r.get('chainage_m'))}" if _has(r.get("chainage_m")) else "",
            f"{fmt(a.get('area_km2'), 3)} km²" if _has(a.get("area_km2")) else "",
            f"Tc Kirpich {fmt(r.get('tc_kirpich_min'), 0)} min" if _has(r.get("tc_kirpich_min")) else "")
            if v)
        items.append(placemark(r["outlet_uid"], "x_proposed" if r.get("status") == "proposed"
                               else "x_existing", point_geom(p), crossing_card(r, a, f, md), snip))
    folders.append(folder("1 · Crossings (existing blue, proposed orange)", items, open_=True))
    # catchments
    items = []
    for uid, a in ca.items():
        polys = a.get("_geom") or []
        if not polys:
            continue
        geoms = [polygon_geom([to_wgs84(ring) for ring in poly]) for poly in polys]
        lp = to_wgs84([label_point(polys)])[0]
        st = "cat_proposed" if a.get("status") == "proposed" else "cat_existing"
        desc = card(C["proposed"] if st == "cat_proposed" else C["existing"], f"Catchment {uid}",
                    f"Drains to crossing {uid}", "", rows_table([
                        ("Area", fmt(a.get("area_km2"), 4, "km²")),
                        ("Mode", a.get("catchment_mode")),
                        ("Mean slope (Horn)", fmt(a.get("catch_slope_horn"), 4, "m/m")),
                        ("Relief", fmt(a.get("relief_m"), 1, "m")),
                        ("Perimeter", fmt(a.get("perimeter_km"), 3, "km")),
                        ("Form factor / circularity", f"{fmt(a.get('form_factor'), 3)} / "
                         f"{fmt(a.get('circularity_ratio'), 3)}" if _has(a.get("form_factor")) else None),
                        ("Drainage density", fmt(a.get("drainage_density"), 3, "km/km²")),
                        ("CN exported", fmt(a.get("cn_export"), 1)),
                        ("Land cover (largest share)", max(
                            ((k[7:], a[k]) for k in a if k.startswith("lc_pct_") and _has(a.get(k))),
                            key=lambda t: t[1], default=(None, None))[0])]))
        items.append(placemark(uid, st, "<MultiGeometry>" + point_geom(lp) + "".join(geoms)
                               + "</MultiGeometry>", desc, fmt(a.get("area_km2"), 3, "km²")))
    folders.append(folder("2 · Catchments", items))
    # flow paths
    items = []
    for uid, f in fp.items():
        line = f.get("_geom") or []
        if len(line) < 2:
            continue
        desc = card(C["lfp"], f"Longest flow path {uid}", f"Outlet at crossing {uid}", "", rows_table([
            ("Length", fmt((f.get("lfp_length_m") or 0) / 1000.0, 3, "km")),
            ("Elevation max / min", f"{fmt(f.get('lfp_elev_max_m'), 1)} / {fmt(f.get('lfp_elev_min_m'), 1)} m"),
            ("Drop", fmt(f.get("lfp_drop_m"), 1, "m")),
            ("Slope (drop / length)", fmt(f.get("lfp_slope"), 4, "m/m")),
            ("Slope 10-85", fmt(f.get("lfp_slope_1085"), 4, "m/m")),
            ("Overland / sheet / shallow", " / ".join(fmt(f.get(k), 0) or "–" for k in
                                                       ("lfp_overland_m", "lfp_sheet_m", "lfp_shallow_m")) + " m"),
            ("Channel length / slope", f"{fmt(f.get('lfp_channel_m'), 0)} m / "
             f"{fmt(f.get('lfp_channel_slope'), 4)} m/m" if _has(f.get("lfp_channel_m")) else None)]))
        items.append(placemark(uid, "lfp", line_geom(to_wgs84(line)), desc))
    folders.append(folder("3 · Longest flow paths", items))
    # road and kilometre posts
    road = _layer(pkg, "road_alignment", tabs)
    if road:
        items, posts = [], []
        for r in road:
            line = r.get("_geom") or []
            if len(line) < 2:
                continue
            ll = to_wgs84(line)
            allpts += ll
            desc = card(C["kmpost"], "Road alignment", md.get("alignment_source") or "", "", rows_table([
                ("Chainage", f"{chainage(r.get('chainage_from_m'))} to {chainage(r.get('chainage_to_m'))}"),
                ("Part", fmt(r.get("part")))]))
            items.append(placemark("Road (casing)", "road_case", line_geom(ll)))
            items.append(placemark(f"Road part {r.get('part')}", "road", line_geom(ll), desc))
            c0, c1 = r.get("chainage_from_m") or 0.0, r.get("chainage_to_m") or 0.0
            km = math.ceil(c0 / 1000.0) * 1000.0
            while km <= c1 + 1e-6:
                p = to_wgs84([along(line, km - c0)])[0]
                posts.append(placemark(chainage(km)[:-3], "kmpost", point_geom(p),
                                       card(C["kmpost"], f"Ch {chainage(km)}",
                                            md.get("alignment_source") or "")))
                km += 1000.0
        items.append(folder("Kilometre posts", posts))
        folders.append(folder("4 · Road", items))
    # coverage
    cov = _layer(pkg, "coverage_check", tabs)
    sags = _layer(pkg, "sag_points", tabs)
    flats = _layer(pkg, "flat_stretches", tabs)
    if cov or sags or flats:
        items = []
        for r in cov:
            if r.get("issue") == "flat_stretch":
                continue
            p = to_wgs84([r["_geom"]])[0]
            items.append(placemark(f"{(r.get('issue') or '').replace('_', ' ')} "
                                   f"{chainage(r.get('chainage_m')) or ''}", "finding", point_geom(p),
                                   card(C["finding"], (r.get("issue") or "").replace("_", " ").capitalize(),
                                        f"Ch {chainage(r.get('chainage_m'))}", "", rows_table([
                                            ("Area", fmt(r.get("area_km2"), 4, "km²")),
                                            ("Proposed crossing", r.get("uid")),
                                            ("Nearest existing", r.get("nearest_uid")),
                                            ("Finding", r.get("note"))])), r.get("note") or ""))
        for r in sags:
            p = to_wgs84([r["_geom"]])[0]
            items.append(placemark(f"Sag {chainage(r.get('chainage_m'))}", "sag", point_geom(p),
                                   card(C["sag"], f"Sag point Ch {chainage(r.get('chainage_m'))}",
                                        f"{r.get('side') or ''} side", "", rows_table([
                                            ("Ground level", fmt(r.get("z_dem_m"), 2, "m")),
                                            ("Sag depth", fmt(r.get("sag_depth_m"), 2, "m")),
                                            ("Local area", fmt(r.get("sag_area_km2"), 4, "km²")),
                                            ("Proposed crossing", r.get("uid"))]))))
        for r in flats:
            line = r.get("_geom") or []
            if len(line) >= 2:
                items.append(placemark(f"Flat stretch {chainage(r.get('chainage_m'))}", "flat",
                                       line_geom(to_wgs84(line)),
                                       card(C["flat"], "Flat stretch", f"Ch {chainage(r.get('chainage_m'))} "
                                            f"to {chainage(r.get('chainage_to_m'))}", "", rows_table([
                                                ("Length", fmt(r.get("length_m"), 0, "m")),
                                                ("Longitudinal slope", fmt(r.get("slope_long_pct"), 2, "%")),
                                                ("Cross-fall", fmt(r.get("crossfall_pct"), 2, "%")),
                                                ("Crossings within", r.get("uids_within"))]))))
        folders.append(folder("5 · Coverage findings, sags and flat stretches", items))
    # channel sections and side drains (off by default)
    xs = _layer(pkg, "xs_transects", tabs)
    if xs:
        folders.append(folder("6 · Channel sections (indicative)", [
            placemark(r.get("outlet_uid"), "xs", line_geom(to_wgs84(r["_geom"])),
                      card(C["xs"], f"Channel section {r.get('outlet_uid')}",
                           f"{fmt(r.get('xs_dist_m'), 0)} m downstream · quality {r.get('xs_quality')}"))
            for r in xs if r.get("_geom")], visible=False))
    dv = _layer(pkg, "drainage_divergence", tabs)
    mr = _layer(pkg, "mapped_rivers_used", tabs)
    if dv or mr:
        items = [placemark(r.get("name") or "mapped waterway", "mapped", line_geom(to_wgs84(r["_geom"])),
                           card(C["mapped"], r.get("name") or "Mapped waterway", "used for the check",
                                "", rows_table([("Length in the DEM extent", fmt(r.get("length_m"), 0, "m"))])))
                 for r in mr if r.get("_geom")]
        items += [placemark(f"Divergence {fmt(r.get('area_km2'), 3)} km²", "divergence",
                            line_geom(to_wgs84(r["_geom"])),
                            card(C["bad"], "DEM stream leaves the mapped course", "", "", rows_table([
                                ("Area flowing down it", fmt(r.get("area_km2"), 3, "km²")),
                                ("Length", fmt(r.get("length_m"), 0, "m")),
                                ("Returns to the mapped course", "yes" if r.get("rejoins") else "no")]),
                                foot="Verify on site: the catchment may switch between crossings."))
                  for r in dv if r.get("_geom")]
        folders.append(folder("6a · Mapped drainage check", items, visible=True))
    sti = _layer(pkg, "corridor_sti", tabs)
    if sti:
        items = []
        order = list(STI_COLOURS)
        for r in sti:
            cls = [c for c in (r.get("sti_class_lhs"), r.get("sti_class_rhs")) if c in STI_COLOURS]
            if not cls or not r.get("_geom"):
                continue
            worst = max(cls, key=order.index)
            items.append(placemark(f"{chainage(r.get('ch_start'))} - {chainage(r.get('ch_end'))}",
                                   "sti_" + worst.replace(" ", "_"), line_geom(to_wgs84(r["_geom"])),
                                   card(STI_COLOURS[worst], "Side-drain siltation (screening)",
                                        f"Ch {chainage(r.get('ch_start'))} to {chainage(r.get('ch_end'))}",
                                        "", rows_table([
                                            ("Left: class / STI p90", f"{r.get('sti_class_lhs') or '–'} / "
                                             f"{fmt(r.get('sti_p90_lhs'), 2) or '–'}"),
                                            ("Right: class / STI p90", f"{r.get('sti_class_rhs') or '–'} / "
                                             f"{fmt(r.get('sti_p90_rhs'), 2) or '–'}"),
                                            ("Flagged length left / right",
                                             f"{fmt(r.get('siltation_len_lhs_m'), 0) or '–'} / "
                                             f"{fmt(r.get('siltation_len_rhs_m'), 0) or '–'} m")]),
                                        foot="Relative transport capacity within this corridor, "
                                             "not severity.")))
        folders.append(folder("7 · Side-drain siltation (relative)", items, visible=False))
    # relief overlay
    if relief is not None:
        png_path, (xmin, ymin, xmax, ymax) = relief
        with open(png_path, "rb") as fh:
            files["files/relief.png"] = fh.read()
        q = to_wgs84([(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)])
        folders.append(folder("8 · Terrain", [
            f"<GroundOverlay><name>Relief (quicklook)</name><visibility>1</visibility>"
            f"<color>c8ffffff</color><drawOrder>0</drawOrder><Icon><href>files/relief.png</href></Icon>"
            f"<gx:LatLonQuad><coordinates>{coords_text(q)}</coordinates></gx:LatLonQuad></GroundOverlay>"]))
    # legend and title block
    overlays = []
    leg = legend_png(bool(n_prop), bool(cov or sags), bool(flats), bool(sti))
    tb = title_png(title, md)
    if leg:
        files["files/legend.png"] = leg
        overlays.append("<ScreenOverlay><name>Legend</name><visibility>1</visibility><Icon><href>"
                        "files/legend.png</href></Icon><overlayXY x=\"0\" y=\"1\" xunits=\"fraction\" "
                        "yunits=\"fraction\"/><screenXY x=\"0.005\" y=\"0.985\" xunits=\"fraction\" "
                        "yunits=\"fraction\"/><size x=\"0\" y=\"0\" xunits=\"pixels\" yunits=\"pixels\"/>"
                        "</ScreenOverlay>")
    if tb:
        files["files/title_block.png"] = tb
        overlays.append("<ScreenOverlay><name>Title block</name><visibility>1</visibility><Icon><href>"
                        "files/title_block.png</href></Icon><overlayXY x=\"0\" y=\"0\" xunits=\"fraction\" "
                        "yunits=\"fraction\"/><screenXY x=\"0.005\" y=\"0.035\" xunits=\"fraction\" "
                        "yunits=\"fraction\"/><size x=\"0\" y=\"0\" xunits=\"pixels\" yunits=\"pixels\"/>"
                        "</ScreenOverlay>")
    if overlays:
        folders.insert(0, folder("Legend and title", overlays))
    # view
    if allpts:
        a = np.asarray(allpts)
        lon0, lat0 = a[:, 0].mean(), a[:, 1].mean()
        span = max((a[:, 0].max() - a[:, 0].min()) * 111320 * math.cos(math.radians(lat0)),
                   (a[:, 1].max() - a[:, 1].min()) * 110574, 500.0)
        look = (f"<LookAt><longitude>{lon0:.6f}</longitude><latitude>{lat0:.6f}</latitude>"
                f"<altitude>0</altitude><heading>0</heading><tilt>0</tilt><range>{1.6 * span:.0f}"
                f"</range><altitudeMode>relativeToGround</altitudeMode></LookAt>")
    else:
        look = ""
    kml = ("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<kml xmlns=\"http://www.opengis.net/kml/2.2\" "
           "xmlns:gx=\"http://www.google.com/kml/ext/2.2\"><Document>"
           f"<name>{x(title)}</name><open>1</open>{look}"
           f"<Snippet maxLines=\"1\">QEHT {x(md.get('qeht_version', ''))} · "
           f"{x(datetime.date.today().isoformat())}</Snippet>"
           f"<description>{cdata(run_card(md, n_ex, n_prop, title))}</description>"
           + "".join(styles) + "".join(folders) + "</Document></kml>")
    return kml, files


def write_kmz(pkg, out_path, to_wgs84, title=None, relief=None):
    kml, files = build_kml(pkg, to_wgs84, title, relief)
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with zipfile.ZipFile(out_path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("doc.kml", kml.encode("utf-8"))
        for name, data in files.items():
            z.writestr(name, data)
    return out_path


# -- legend and title block (matplotlib when available) ----------------------------

def _mpl():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except Exception:  # nosec B110 - optional dependency, overlays are skipped
        return None


def _fig_png(fig, plt):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100, transparent=True)
    plt.close(fig)
    return buf.getvalue()


def legend_png(proposed=True, coverage=True, flats=True, sti=True):
    plt = _mpl()
    if plt is None:
        return None
    from matplotlib.patches import FancyBboxPatch, Rectangle
    entries = [("Crossings", None), ("circle", C["existing"], "Existing crossing")]
    if proposed:
        entries.append(("circle", C["proposed"], "Proposed crossing"))
    entries += [("Catchments", None), ("fill", C["existing"], "Catchment (existing)")]
    if proposed:
        entries.append(("fill", C["proposed"], "Catchment (proposed)"))
    entries += [("Lines", None), ("line", C["lfp"], "Longest flow path"),
                ("road", C["road"], "Road alignment"), ("square", C["kmpost"], "Kilometre post")]
    if coverage:
        entries += [("Coverage check", None), ("triangle", C["finding"], "Finding"),
                    ("diamond", C["sag"], "Sag point")]
    if flats:
        entries.append(("line", C["flat"], "Flat stretch"))
    if sti:
        entries += [("Side drains (relative)", None)] + [("line", v, k.capitalize())
                                                        for k, v in STI_COLOURS.items()]
    h = 0.55 + sum(0.40 if e[1] is None else 0.26 for e in entries) + 0.2
    fig = plt.figure(figsize=(2.6, h))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 2.6); ax.set_ylim(0, h); ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.04, 0.04), 2.52, h - 0.08, boxstyle="round,pad=0,rounding_size=0.08",
                                fc="white", ec="#7f8c99", lw=1.0))
    ax.add_patch(Rectangle((0.04, h - 0.5), 2.52, 0.46, fc=C["navy"], ec="none"))
    ax.text(1.3, h - 0.27, "QEHT DESIGN HYDROLOGY", color="white", ha="center", va="center",
            fontsize=9.5, fontweight="bold")
    y = h - 0.55
    for e in entries:
        if e[1] is None:
            y -= 0.26
            ax.text(0.14, y, e[0], fontsize=8.5, fontweight="bold", color=C["navy"], va="center")
            ax.plot([0.14, 2.46], [y - 0.12, y - 0.12], color=C["rule"], lw=0.8)
            y -= 0.14
            continue
        kind, col, label = e
        y -= 0.26
        cx = 0.32
        if kind == "circle":
            ax.scatter([cx], [y], s=90, c=col, edgecolors="#1d2733", linewidths=1.0, zorder=3)
        elif kind == "square":
            ax.scatter([cx], [y], s=55, marker="s", c=col, edgecolors="#1d2733", linewidths=1.0)
        elif kind == "diamond":
            ax.scatter([cx], [y], s=80, marker="D", c=col, edgecolors="#1d2733", linewidths=1.0)
        elif kind == "triangle":
            ax.scatter([cx], [y], s=90, marker="v", c=col, edgecolors="#1d2733", linewidths=1.0)
        elif kind == "fill":
            ax.add_patch(Rectangle((cx - 0.13, y - 0.09), 0.26, 0.18, fc=col, alpha=0.3, ec=col, lw=1.6))
        elif kind == "road":
            ax.plot([cx - 0.15, cx + 0.15], [y, y], color=C["road_case"], lw=6, solid_capstyle="butt")
            ax.plot([cx - 0.15, cx + 0.15], [y, y], color=col, lw=2.5, solid_capstyle="butt")
        else:
            ax.plot([cx - 0.15, cx + 0.15], [y, y], color=col, lw=2.5)
        ax.text(0.6, y, label, fontsize=8, va="center", color=C["text"])
    return _fig_png(fig, plt)


def title_png(title, md):
    plt = _mpl()
    if plt is None:
        return None
    from matplotlib.patches import FancyBboxPatch
    fig = plt.figure(figsize=(5.9, 0.9))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 5.9); ax.set_ylim(0, 0.9); ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.03, 0.03), 5.84, 0.84, boxstyle="round,pad=0,rounding_size=0.08",
                                fc="#2f4a63", ec="none"))
    ax.text(0.18, 0.66, title, color="white", fontsize=11, fontweight="bold", va="center")
    ax.text(0.18, 0.42, f"DEM: {md.get('dem_source', '')} · cell {md.get('cell_size_m', '')} m · "
            f"flat method {md.get('flat_method', '')}", color="#dbe6f0", fontsize=7.5, va="center")
    ax.text(0.18, 0.2, f"Package CRS EPSG:{md.get('crs_epsg', '')}, shown in WGS 84 · QEHT "
            f"{md.get('qeht_version', '')} · {datetime.date.today().isoformat()} · terrain values, "
            "check against survey", color="#9fb3c6", fontsize=6.8, va="center")
    return _fig_png(fig, plt)
