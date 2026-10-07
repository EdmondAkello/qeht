"""Check that a freshly built golden package only ADDS fields/keys to the
committed fixture (no existing value changed). Run from the folder that contains the qeht repo:
    python3 golden_diff.py
"""
import json, os, sys, tempfile
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
from qeht.tests import test_interop as ti

new = ti._summary(ti._golden_package(os.path.join(tempfile.mkdtemp(), "g.gpkg")))
ref = json.load(open(ti.FIXTURE_JSON, encoding="utf-8"))
changed, added, removed = [], set(), []
for layer in ref:
    if layer not in new:
        removed.append(layer); continue
    if isinstance(ref[layer], dict):
        for k, v in ref[layer].items():
            if k not in new[layer]: removed.append(f"{layer}.{k}")
            elif new[layer][k] != v and k not in ("run_time_utc", "qeht_version"):
                changed.append(f"{layer}.{k}: {v!r} -> {new[layer][k]!r}")
        added |= {f"{layer}.{k}" for k in set(new[layer]) - set(ref[layer])}
    elif ref[layer] and isinstance(ref[layer][0], dict):
        if len(ref[layer]) != len(new[layer]):
            changed.append(f"{layer}: {len(ref[layer])} -> {len(new[layer])} rows"); continue
        for i, (a, b) in enumerate(zip(ref[layer], new[layer])):
            for k, v in a.items():
                if k not in b: removed.append(f"{layer}[{i}].{k}")
                elif b[k] != v: changed.append(f"{layer}[{i}].{k}: {v!r} -> {b[k]!r}")
            added |= {f"{layer}.{k}" for k in set(b) - set(a)}
    elif ref[layer] != new[layer]:
        changed.append(f"{layer}: {ref[layer]} -> {new[layer]}")
added |= {f"{k} (layer)" for k in set(new) - set(ref)}
print("ADDED:", ", ".join(sorted(added)) or "none")
print("REMOVED:", ", ".join(removed) or "none")
print("CHANGED:", "\n  ".join(changed) or "none")
sys.exit(1 if (changed or removed) else 0)
