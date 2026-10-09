"""Assemble the user guide HTML from user_guide_src.html.

    python3 qeht/tools/guide/build_guide.py <run folder> <figure folder> [<output folder>]

<run folder> is a pipeline run with a kept package (package/design_hydrology.gpkg): its field
dictionary becomes Appendix B. Appendix A is generated from the plugin's own parameter
definitions, so QGIS must be importable (run with the Python QGIS uses). The output folder
(default qeht/docs/user-guide) gets user-guide.html and figures/ with the figures it uses.
Render the PDF with render_pdf.js.
"""

import datetime
import html
import os
import re
import shutil
import sqlite3
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.dirname(REPO))

GROUPS = ["Workflow", "Terrain and drainage", "Road drainage", "Soils and erosion", "Interoperability"]
LAYERS = ["crossings", "catchments", "flowpaths"]


def version():
    for line in open(os.path.join(REPO, "metadata.txt"), encoding="utf-8"):
        if line.startswith("version="):
            return line.split("=", 1)[1].strip()
    return "unknown"


def algorithms():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from qgis.core import QgsApplication, QgsProcessingParameterDefinition as D
    app = QgsApplication([], False)
    app.initQgis()
    from qeht.processing_provider.provider import QehtProvider
    prov = QehtProvider()
    QgsApplication.processingRegistry().addProvider(prov)
    out = []
    for a in prov.algorithms():
        params = []
        for d in a.parameterDefinitions():
            dv = d.defaultValue()
            if d.type() == "enum" and dv is not None:
                try:
                    dv = d.options()[int(dv)]
                except (ValueError, IndexError, TypeError):
                    pass
            if d.type() in ("rasterDestination", "vectorDestination", "fileDestination",
                            "folderDestination", "sink"):
                dv = "output"
            params.append((d.name(), d.description(), "" if dv in (None, "") else str(dv),
                           bool(d.flags() & D.FlagOptional), bool(d.flags() & D.FlagAdvanced)))
        out.append((a.group(), a.displayName(), a.id(), a.shortHelpString(), params))
    return out


def first_sentence(text):
    t = re.sub(r"<[^>]+>", "", text or "").strip()
    m = re.match(r"(.+?\.)(\s|$)", t, re.S)
    return (m.group(1) if m else t).replace("\n", " ")


def appendix_params(algs):
    parts = []
    for g in GROUPS:
        parts.append(f'<h2>{html.escape(g)}</h2>')
        for group, name, aid, helptext, params in sorted(a for a in algs if a[0] == g):
            parts.append(f'<h3 class="alg">{html.escape(name)} <span class="small">· <code>{aid}</code></span></h3>')
            parts.append(f'<p class="small">{html.escape(first_sentence(helptext))}</p>')
            rows = []
            for n, desc, dv, opt, adv in params:
                flag = " <span class='small'>(advanced)</span>" if adv else ""
                rows.append(f"<tr><td>{n}</td><td>{html.escape(desc)}{flag}</td><td>{html.escape(dv) or '—'}</td></tr>")
            parts.append('<table class="params"><thead><tr><th style="width:24%">Name</th><th>Label</th>'
                         '<th style="width:16%">Default</th></tr></thead><tbody>' + "".join(rows) + "</tbody></table>")
    return "\n".join(parts)


def appendix_fields(run):
    db = sqlite3.connect(os.path.join(run, "package", "design_hydrology.gpkg"))
    rows = db.execute("select layer, field, unit, meaning from qeht_field_dictionary").fetchall()
    layers = LAYERS + sorted({r[0] for r in rows} - set(LAYERS))
    parts = []
    for lyr in layers:
        rs = [r for r in rows if r[0] == lyr]
        if not rs:
            continue
        parts.append(f"<h2><code>{html.escape(lyr)}</code> <span class='small'>({len(rs)} fields)</span></h2>")
        body = "".join(f"<tr><td>{html.escape(f)}</td><td>{html.escape(u or '')}</td><td>{html.escape(m or '')}</td></tr>"
                       for _, f, u, m in rs)
        parts.append('<table class="fields"><thead><tr><th style="width:26%">Field</th><th style="width:12%">Unit</th>'
                     '<th>Meaning</th></tr></thead><tbody>' + body + "</tbody></table>")
    return "\n".join(parts)


def main(run, figs, out):
    src = open(os.path.join(HERE, "user_guide_src.html"), encoding="utf-8").read()
    v = version()
    doc = (src.replace("{{VERSION}}", v)
              .replace("{{EDITION}}", datetime.date.today().strftime("%B %Y"))
              .replace("{{APPENDIX_PARAMS}}", appendix_params(algorithms()))
              .replace("{{APPENDIX_FIELDS}}", appendix_fields(run)))
    os.makedirs(os.path.join(out, "figures"), exist_ok=True)
    used = sorted(set(re.findall(r"figures/([\w\-]+\.png)", doc)))
    missing = [f for f in used if not os.path.exists(os.path.join(figs, f))]
    if missing:
        sys.exit("missing figures: " + ", ".join(missing))
    for f in used:
        shutil.copy2(os.path.join(figs, f), os.path.join(out, "figures", f))
    with open(os.path.join(out, "user-guide.html"), "w", encoding="utf-8") as f:
        f.write(doc)
    print(f"wrote {out}/user-guide.html (QEHT {v}, {len(used)} figures)")


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) < 2:
        sys.exit(__doc__)
    main(a[0], a[1], a[2] if len(a) > 2 else os.path.join(REPO, "docs", "user-guide"))
