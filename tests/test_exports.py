# -*- coding: utf-8 -*-
"""v0.23 validation: KMZ and XLSX exports (F11).
NumPy + OSR for the test's own coordinate transform; the writers are pure
Python. The package is the synthetic three-crossing case of test_interop,
with a road and an ID full of XML special characters.

    python -m qeht.tests.test_exports
"""

import math
import os
import sys
import tempfile
import xml.etree.ElementTree as ET  # nosec B405 - parses our own output
import zipfile

from ..core.interop.heas_exchange import write_exchange
from ..core.report.kmz import write_kmz, icon_png, chainage
from ..core.report.xlsx import write_xlsx, Sheet, read_xlsx_cells, col_letter
from ..core.report.xlsx_layout import write_workbook, CHARACTERISTICS, TC
from .test_interop import synthetic_run, UTM37S_WKT

FAILURES = []
KML = "{http://www.opengis.net/kml/2.2}"
ODD = 'X<2>&"q\''


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def transformer():
    from osgeo import osr
    osr.UseExceptions()
    s = osr.SpatialReference(); s.ImportFromWkt(UTM37S_WKT)
    g = osr.SpatialReference(); g.ImportFromEPSG(4326)
    for r in (s, g):
        r.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    t = osr.CoordinateTransformation(s, g)
    return lambda pts: [tuple(p[:2]) for p in t.TransformPoints([list(map(float, q)) for q in pts])]


def package(tmp):
    (cr, ca, fp, issues, id_info), g = synthetic_run()
    for layer in (cr, ca, fp):
        for _, a in layer:
            if a["outlet_uid"] == "X002":
                a["outlet_uid"] = ODD
                for k in ("crossing_id", "catchment_id", "flowpath_id"):
                    if k in a:
                        a[k] = ODD
    gt = g["gt"]
    y = gt[3] - 45.5 * 10.0
    road = [((gt[0] + 5.0, y), (gt[0] + 495.0, y))]
    extra = [("road_alignment", "LINESTRING",
              [("part", "int"), ("chainage_from_m", "real"), ("chainage_to_m", "real")],
              [([road[0][0], road[0][1]], {"part": 1, "chainage_from_m": 600.0,
                                           "chainage_to_m": 1090.0})], "Road alignment used")]
    md = dict(id_info, dem_source="synthetic <test> & co", qeht_version="test",
              cell_size_m="10 x 10", stream_threshold_km2="0.01")
    path = os.path.join(tmp, "pkg.gpkg")
    write_exchange(path, cr, ca, fp, md, UTM37S_WKT, crs_name="WGS 84 / UTM zone 37S",
                   extra_layers=extra)
    return path, cr


def test_kmz(tmp):
    print("\n1. KMZ")
    pkg, cr = package(tmp)
    tw = transformer()
    out = write_kmz(pkg, os.path.join(tmp, "out.kmz"), tw, title="Site <B> & road")
    z = zipfile.ZipFile(out)
    names = z.namelist()
    check("zip: doc.kml first, icons in files/", names[0] == "doc.kml"
          and "files/circle_existing.png" in names and "files/square_kmpost.png" in names)
    root = ET.fromstring(z.read("doc.kml"))  # nosec B314
    doc = root.find(KML + "Document")
    check("document title escaped and parsed back", doc.find(KML + "name").text == "Site <B> & road")
    folders = {f.find(KML + "name").text: f for f in doc.iter(KML + "Folder")}

    def count(prefix):
        f = [v for k, v in folders.items() if k.startswith(prefix)]
        return len(f[0].findall(KML + "Placemark")) if f else -1
    check("folder counts: 3 crossings, 3 catchments, 3 flow paths",
          (count("1 ·"), count("2 ·"), count("3 ·")) == (3, 3, 3),
          str((count("1 ·"), count("2 ·"), count("3 ·"))))
    kms = [p.find(KML + "name").text for p in folders["Kilometre posts"].findall(KML + "Placemark")]
    check("kilometre post at 1+000 only (road 0+600 to 1+090)", kms == ["1+000"], str(kms))
    pms = {p.find(KML + "name").text: p for p in folders[[k for k in folders if k.startswith("1 ·")][0]]
           .findall(KML + "Placemark")}
    check("an ID with < > & \" ' survives as a placemark name", ODD in pms, str(list(pms)))
    worst = 0.0
    for _, a in cr:
        lon, lat = tw([(a["outlet_x"], a["outlet_y"])])[0]
        txt = pms[a["outlet_uid"]].find(f"{KML}Point/{KML}coordinates").text
        lo, la, _ = (float(v) for v in txt.split(","))
        worst = max(worst, abs(lo - lon), abs(la - lat))
    check("crossing coordinates within 1e-7 degree of the OSR transform", worst <= 1e-7, f"{worst:.1e}")
    desc = pms[ODD].find(KML + "description").text
    check("card HTML escapes the ID and carries the sections",
          "X&lt;2&gt;&amp;&quot;q" in desc and "Contributing catchment" in desc and "Location" in desc)
    run = doc.find(KML + "description").text
    check("run card escapes the DEM source", "synthetic &lt;test&gt; &amp; co" in run)
    poly = folders[[k for k in folders if k.startswith("2 ·")][0]].find(KML + "Placemark")
    check("catchments: label point + polygon in a MultiGeometry, style map with highlight",
          poly.find(f"{KML}MultiGeometry/{KML}Point") is not None
          and poly.find(f"{KML}MultiGeometry/{KML}Polygon") is not None
          and doc.find(f"{KML}StyleMap[@id='cat_existing']") is not None)
    look = doc.find(KML + "LookAt")
    check("LookAt over the crossings", look is not None and 36.0 < float(look.find(KML + "longitude").text) < 40.0)
    check("chainage format", chainage(1234.5) == "1+234.50" and chainage(1000.0) == "1+000.00"
          and chainage(None) is None)
    png = icon_png("circle", "#2b83ba", 32)
    check("icon is a PNG", png[:8] == b"\x89PNG\r\n\x1a\n")
    try:
        import matplotlib  # noqa: F401
        check("legend and title block overlays (matplotlib present)",
              "files/legend.png" in names and "files/title_block.png" in names)
    except ImportError:
        check("no matplotlib: no overlays, still valid", "files/legend.png" not in names)


def test_xlsx(tmp):
    print("\n2. XLSX")
    p = os.path.join(tmp, "t.xlsx")
    write_xlsx(p, [Sheet("Numbers & text", ["a", "b", "c"], [[1.5, None, 'x<&>"y'], [2, True, "z"]],
                         units=["m", "", "-"], formats=["0.00", None, None], freeze=(2, 1)),
                   Sheet("Second/sheet:name?", ["k"], [["v"]])])
    z = zipfile.ZipFile(p)
    ct = ET.fromstring(z.read("[Content_Types].xml"))  # nosec B314
    check("[Content_Types].xml parses and declares both sheets",
          sum(1 for o in ct if o.get("PartName", "").startswith("/xl/worksheets/")) == 2)
    for n in z.namelist():
        if n.endswith(".xml") or n.endswith(".rels"):
            ET.fromstring(z.read(n))  # nosec B314
    check("every part is well-formed XML", True)
    c = read_xlsx_cells(p, 1)
    check("header and units rows", c["A1"] == "a" and c["A2"] == "m" and "B2" not in c)
    check("numbers round-trip, NULL is an empty cell, text with < & > \" round-trips",
          c["A3"] == 1.5 and "B3" not in c and c["C3"] == 'x<&>"y' and c["A4"] == 2.0 and c["B4"] is True)
    wb = ET.fromstring(z.read("xl/workbook.xml"))  # nosec B314
    sn = [s.get("name") for s in wb.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}sheet")]
    check("sheet names cleaned of illegal characters", sn == ["Numbers & text", "Second_sheet_name_"], str(sn))
    s1 = z.read("xl/worksheets/sheet1.xml").decode()
    check("frozen pane at B3, column widths, number format style",
          'topLeftCell="B3"' in s1 and 'state="frozen"' in s1 and "<cols>" in s1
          and 'formatCode="0.00"' not in z.read("xl/styles.xml").decode() and ' s="3"' in s1)
    check("column letters", [col_letter(i) for i in (0, 25, 26, 701, 702)] == ["A", "Z", "AA", "ZZ", "AAA"])
    pkg, cr = package(tmp)
    out = write_workbook(pkg, os.path.join(tmp, "wb.xlsx"), extra={"X001": {"flat_sensitive": 1}})
    wbn = [s.get("name") for s in ET.fromstring(zipfile.ZipFile(out).read("xl/workbook.xml"))  # nosec B314
           .iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}sheet")]
    check("workbook: five sheets in order", wbn == ["Catchment characteristics", "Time of concentration",
                                                     "Coverage findings", "Run metadata", "Field dictionary"])
    cc = read_xlsx_cells(out, 1)
    hdr = [cc.get(f"{col_letter(i)}1") for i in range(len(CHARACTERISTICS))]
    check("characteristics columns in the layout order, units in row 2",
          hdr == [h for _, _, h, _, _ in CHARACTERISTICS] and cc["F2"] == "km²")
    ids = [cc.get(f"A{r}") for r in (3, 4, 5)]
    area = {a["outlet_uid"] for _, a in cr}
    fi = [h for _, _, h, _, _ in CHARACTERISTICS].index("Flat-sensitive")
    check("one row per crossing; the odd ID and the extra column come through",
          sorted(ids) == sorted(area) and cc.get(f"{col_letter(fi)}{3 + ids.index('X001')}") == 1.0,
          str(ids))
    md = read_xlsx_cells(out, 4)
    check("run metadata sheet carries the DEM source", any(v == "synthetic <test> & co" for v in md.values()))
    check("TC sheet header", read_xlsx_cells(out, 2)["F1"] == TC[5][2])


def main(argv=None):
    tmp = tempfile.mkdtemp(prefix="qeht_exp_")
    test_kmz(tmp)
    test_xlsx(tmp)
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
