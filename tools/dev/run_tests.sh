#!/bin/bash
# QEHT full test run. Run from anywhere: bash qeht/tools/dev/run_tests.sh
# - pure suites with python3.12 (needs GDAL: apt install python3-gdal or qgis)
# - QGIS 3.x smoke on the host (apt qgis python3-qgis)
# - QGIS 4.x smoke in the chroot /opt/trixie (tools/dev/build_trixie.sh), skipped if absent
# The repo folder MUST be named "qeht" (tests import qeht.*); we cd to its parent.
set -u
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
PARENT="$(dirname "$REPO")"
[ "$(basename "$REPO")" = "qeht" ] || { echo "rename/clone the repo folder as 'qeht'"; exit 2; }
cd "$PARENT"
PY=${PY:-python3.12}
SUITES="core interop crossings soils soils_any flats erosion alignment morphometry runoff channel coverage rainfall report floodplain quicklooks section corridor_sti prepare_dem tc"
for t in $SUITES; do r=$($PY -m qeht.tests.test_$t 2>&1); echo "pure $t: pass=$(echo "$r"|grep -c '\[PASS\]') fail=$(echo "$r"|grep -c '\[FAIL\]')"; done
QT_QPA_PLATFORM=offscreen $PY -m qeht.tests.qgis_smoke > smoke_3x.log 2>&1
echo "smoke QGIS 3.x: pass=$(grep -c '\[PASS\]' smoke_3x.log) fail=$(grep -c '\[FAIL\]' smoke_3x.log) done=$(grep -c 'ALL QGIS SMOKE CHECKS PASSED' smoke_3x.log)  (a segfault AFTER the final line is a known Ubuntu-3.34 teardown issue)"
if [ -d /opt/trixie/usr/share/qgis/python ]; then
  for m in proc sys dev; do mountpoint -q /opt/trixie/$m || mount --bind /$m /opt/trixie/$m; done
  mkdir -p /opt/trixie/work; mountpoint -q /opt/trixie/work || mount --bind "$PARENT" /opt/trixie/work
  chroot /opt/trixie bash -c "cd /work && PYTHONPATH=/usr/share/qgis/python:/usr/share/qgis/python/plugins QT_QPA_PLATFORM=offscreen python3 -m qeht.tests.qgis_smoke" > smoke_4x.log 2>&1
  echo "smoke QGIS 4.x: pass=$(grep -c '\[PASS\]' smoke_4x.log) fail=$(grep -c '\[FAIL\]' smoke_4x.log) done=$(grep -c 'ALL QGIS SMOKE CHECKS PASSED' smoke_4x.log)"
else
  echo "smoke QGIS 4.x: skipped (no /opt/trixie; see tools/dev/build_trixie.sh)"
fi
