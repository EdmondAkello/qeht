#!/bin/bash
# Run the pure (QGIS-free) test suites and fail if any suite errors or reports a [FAIL].
# Used by .github/workflows/tests.yml; also runs locally: bash qeht/tools/dev/ci_pure.sh
# The repo folder must be named "qeht" (tests import qeht.*).
set -u
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
[ "$(basename "$REPO")" = "qeht" ] || { echo "clone the repo into a folder named 'qeht'"; exit 2; }
cd "$(dirname "$REPO")"
PY=${PY:-python3}
status=0
total=0
for f in qeht/tests/test_*.py; do
  m=$(basename "$f" .py)
  out=$($PY -m "qeht.tests.$m" 2>&1); rc=$?
  p=$(printf '%s\n' "$out" | grep -c '\[PASS\]')
  x=$(printf '%s\n' "$out" | grep -c '\[FAIL\]')
  total=$((total + p))
  echo "$m: pass=$p fail=$x rc=$rc"
  if [ "$rc" -ne 0 ] || [ "$x" -ne 0 ]; then
    status=1
    printf '%s\n' "$out" | grep -E '\[FAIL\]|Error|Traceback' | head -20
  fi
done
echo "total checks passed: $total"
exit $status
