# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Minimal XLSX (Office Open XML) writer, pure Python (F11, v0.23).

No openpyxl: pip is blocked on managed machines. Supports inline strings,
numbers, a bold header row, an italic units row, frozen panes, column
widths and simple number formats - what a hydrology workbook needs.

    write_xlsx(path, [Sheet("Name", header, rows, units=..., formats=..., widths=...,
                            freeze=(2, 1))])

No QGIS imports, no GDAL.
"""

import datetime
import math
import re
import zipfile
from html import escape as _html_escape


def escape(text, entities=None):
    """XML text escape (&, <, > and double quotes)."""
    return _html_escape(text, quote=True).replace("&#x27;", "'")

_BAD = re.compile("[\x00-\x08\x0b\x0c\x0e-\x1f￾￿]")
BUILTIN = {"General": 0, "0": 1, "0.00": 2, "#,##0": 3, "#,##0.00": 4}


class Sheet(object):
    def __init__(self, name, header, rows, units=None, formats=None, widths=None, freeze=None):
        self.name = name
        self.header = list(header)
        self.rows = [list(r) for r in rows]
        self.units = list(units) if units else None
        self.formats = list(formats) if formats else [None] * len(self.header)
        self.widths = list(widths) if widths else None
        self.freeze = freeze          # (rows, cols) kept visible, e.g. (2, 1) -> pane at B3


def col_letter(i):
    """0 -> A, 25 -> Z, 26 -> AA."""
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def _text(v):
    return escape(_BAD.sub("", str(v)))


def sheet_name(name, used):
    n = re.sub(r"[\[\]\*\?/\\:]", "_", str(name)).strip("'")[:31] or "Sheet"
    base, k = n, 2
    while n.lower() in used:
        suffix = f" ({k})"
        n = base[:31 - len(suffix)] + suffix
        k += 1
    used.add(n.lower())
    return n


def auto_width(header, units, rows, col):
    vals = [header[col]] + ([units[col]] if units else []) + [
        r[col] for r in rows[:500] if col < len(r) and r[col] is not None]
    w = max((len(f"{v:.4g}") if isinstance(v, float) else len(str(v))) for v in vals) if vals else 8
    return float(min(max(w + 2, 8), 60))


def write_xlsx(path, sheets, title="", creator="QEHT"):
    """Write the workbook; returns path."""
    # styles: 0 default, 1 header, 2 units, then one per custom/builtin number format
    fmt_codes = []
    for sh in sheets:
        for f in sh.formats:
            if f and f not in fmt_codes:
                fmt_codes.append(f)
    num_ids, custom = {}, []
    for f in fmt_codes:
        if f in BUILTIN:
            num_ids[f] = BUILTIN[f]
        else:
            num_ids[f] = 164 + len(custom)
            custom.append(f)
    style_of = {f: 3 + i for i, f in enumerate(fmt_codes)}
    xfs = ['<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>',
           '<xf numFmtId="0" fontId="1" fillId="2" borderId="1" xfId="0" applyFont="1" applyFill="1" '
           'applyBorder="1"><alignment wrapText="1" vertical="center"/></xf>',
           '<xf numFmtId="0" fontId="2" fillId="0" borderId="1" xfId="0" applyFont="1" applyBorder="1"/>']
    xfs += [f'<xf numFmtId="{num_ids[f]}" fontId="0" fillId="0" borderId="0" xfId="0" '
            f'applyNumberFormat="1"/>' for f in fmt_codes]
    styles = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
              '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
              + (f'<numFmts count="{len(custom)}">' + "".join(
                  f'<numFmt numFmtId="{164 + i}" formatCode="{escape(c)}"/>'
                  for i, c in enumerate(custom)) + "</numFmts>" if custom else "")
              + '<fonts count="3"><font><sz val="10"/><name val="Arial"/></font>'
              '<font><b/><sz val="10"/><color rgb="FFFFFFFF"/><name val="Arial"/></font>'
              '<font><i/><sz val="9"/><color rgb="FF51606F"/><name val="Arial"/></font></fonts>'
              '<fills count="3"><fill><patternFill patternType="none"/></fill>'
              '<fill><patternFill patternType="gray125"/></fill>'
              '<fill><patternFill patternType="solid"><fgColor rgb="FF1F3B57"/><bgColor indexed="64"/>'
              '</patternFill></fill></fills>'
              '<borders count="2"><border><left/><right/><top/><bottom/><diagonal/></border>'
              '<border><left/><right/><top/><bottom style="thin"><color rgb="FFC9D3DD"/></bottom>'
              '<diagonal/></border></borders>'
              '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
              f'<cellXfs count="{len(xfs)}">' + "".join(xfs) + '</cellXfs>'
              '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
              '</styleSheet>')
    used = set()
    names = [sheet_name(s.name, used) for s in sheets]
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml",
                   '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
                   '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
                   '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
                   '<Default Extension="xml" ContentType="application/xml"/>'
                   '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
                   '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
                   + "".join(f'<Override PartName="/xl/worksheets/sheet{i + 1}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
                             for i in range(len(sheets)))
                   + '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>'
                   '</Types>')
        z.writestr("_rels/.rels",
                   '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
                   '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                   '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
                   '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>'
                   '</Relationships>')
        now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        z.writestr("docProps/core.xml",
                   '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
                   '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
                   'xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" '
                   'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
                   f'<dc:title>{_text(title)}</dc:title><dc:creator>{_text(creator)}</dc:creator>'
                   f'<dcterms:created xsi:type="dcterms:W3CDTF">{now}</dcterms:created></cp:coreProperties>')
        z.writestr("xl/workbook.xml",
                   '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
                   '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
                   'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'
                   + "".join(f'<sheet name="{_text(n)}" sheetId="{i + 1}" r:id="rId{i + 1}"/>'
                             for i, n in enumerate(names)) + '</sheets></workbook>')
        z.writestr("xl/_rels/workbook.xml.rels",
                   '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
                   '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                   + "".join(f'<Relationship Id="rId{i + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i + 1}.xml"/>'
                             for i in range(len(sheets)))
                   + f'<Relationship Id="rId{len(sheets) + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
                   '</Relationships>')
        z.writestr("xl/styles.xml", styles)
        for i, sh in enumerate(sheets):
            z.writestr(f"xl/worksheets/sheet{i + 1}.xml", _sheet_xml(sh, style_of))
    return path


def _cell(ref, v, style):
    s = f' s="{style}"' if style else ""
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return ""
    if isinstance(v, bool):
        return f'<c r="{ref}"{s} t="b"><v>{int(v)}</v></c>'
    if isinstance(v, (int, float)):
        return f'<c r="{ref}"{s}><v>{repr(float(v)) if isinstance(v, float) else v}</v></c>'
    return f'<c r="{ref}"{s} t="inlineStr"><is><t xml:space="preserve">{_text(v)}</t></is></c>'


def _sheet_xml(sh, style_of):
    ncol = len(sh.header)
    out = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
           '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
           'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">']
    if sh.freeze and (sh.freeze[0] or sh.freeze[1]):
        fr, fc = sh.freeze
        tl = f"{col_letter(fc)}{fr + 1}"
        pane = "bottomRight" if fr and fc else ("bottomLeft" if fr else "topRight")
        out.append('<sheetViews><sheetView workbookViewId="0"><pane'
                   + (f' xSplit="{fc}"' if fc else "") + (f' ySplit="{fr}"' if fr else "")
                   + f' topLeftCell="{tl}" activePane="{pane}" state="frozen"/>'
                   f'<selection pane="{pane}" activeCell="{tl}" sqref="{tl}"/></sheetView></sheetViews>')
    widths = sh.widths or [auto_width(sh.header, sh.units, sh.rows, c) for c in range(ncol)]
    out.append("<cols>" + "".join(f'<col min="{c + 1}" max="{c + 1}" width="{w:.1f}" customWidth="1"/>'
                                  for c, w in enumerate(widths)) + "</cols><sheetData>")
    r = 1
    out.append(f'<row r="{r}">' + "".join(_cell(f"{col_letter(c)}{r}", h, 1)
                                          for c, h in enumerate(sh.header)) + "</row>")
    if sh.units:
        r += 1
        out.append(f'<row r="{r}">' + "".join(_cell(f"{col_letter(c)}{r}", u, 2)
                                              for c, u in enumerate(sh.units) if u) + "</row>")
    for row in sh.rows:
        r += 1
        cells = []
        for c, v in enumerate(row[:ncol]):
            f = sh.formats[c] if c < len(sh.formats) else None
            st = style_of.get(f) if f and isinstance(v, (int, float)) and not isinstance(v, bool) else 0
            cells.append(_cell(f"{col_letter(c)}{r}", v, st))
        out.append(f'<row r="{r}">' + "".join(cells) + "</row>")
    out.append("</sheetData></worksheet>")
    return "".join(out)


def read_xlsx_cells(path, sheet=1):
    """{ref: value} of one sheet (inline strings, numbers, booleans); for tests and checks."""
    import xml.etree.ElementTree as ET  # nosec B405 - parses our own output
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(path) as z:
        root = ET.fromstring(z.read(f"xl/worksheets/sheet{sheet}.xml"))  # nosec B314
    out = {}
    for c in root.iter(f"{{{ns['m']}}}c"):
        t = c.get("t")
        if t == "inlineStr":
            out[c.get("r")] = "".join(x.text or "" for x in c.iter(f"{{{ns['m']}}}t"))
        else:
            v = c.find("m:v", ns)
            if v is not None:
                out[c.get("r")] = bool(int(v.text)) if t == "b" else float(v.text)
    return out
