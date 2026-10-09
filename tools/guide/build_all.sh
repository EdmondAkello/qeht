#!/usr/bin/env bash
# Build the user guide from scratch: demo site, pipeline run, figures, HTML and PDF.
# Run from the folder that contains the qeht repository:
#   bash qeht/tools/guide/build_all.sh [build folder]
# Needs QGIS (Python 3 with qgis and osgeo), matplotlib, Pillow, Node with Playwright and
# Chromium, and pdfunite (poppler-utils). Use python3.12 on Ubuntu 24.04 (QGIS 3.34 bindings).
set -euo pipefail
PY=${PY:-python3.12}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
export QEHT_GUIDE_BUILD=$(realpath -m "${1:-guide_build}")
export QT_QPA_PLATFORM=offscreen
G=$QEHT_GUIDE_BUILD
mkdir -p "$G/fig" "$G/meta"
$PY "$REPO/tools/guide/make_demo_site.py" "$G/demo"
rm -rf "$G/run"
$PY "$REPO/tools/guide/run_pipeline.py"
(cd "$G" && $PY "$REPO/tools/guide/make_lfp_points.py" X014)
$PY "$REPO/tools/guide/make_maps.py"
$PY "$REPO/tools/guide/make_charts.py"
$PY "$REPO/tools/guide/grab_dialogs.py"
node "$REPO/tools/guide/capture.js" "$G/run/report/run_report.html" "$G/fig"
$PY "$REPO/tools/guide/build_guide.py" "$G/run" "$G/fig"
V=$(grep '^version=' "$REPO/metadata.txt" | cut -d= -f2)
node "$REPO/tools/guide/render_pdf.js" "$REPO/docs/user-guide/user-guide.html" \
     "$REPO/docs/user-guide/QEHT_${V}_User_Guide.pdf" "$V"
