# -*- coding: utf-8 -*-
"""v0.17 validation: catchment characteristics table (§0) and run report (F8).
Bare Python + NumPy, no QGIS, no GDAL.

    python -m qeht.tests.test_report
"""

import csv
import os
import sys
import tempfile

import numpy as np

from .test_interop import synthetic_run, _golden_package, UTM37S_WKT
from ..core.interop.heas_exchange import build_exchange_records, write_exchange
from ..core.report.characteristics import characteristics_rows, write_characteristics_csv, CORE
from ..core.report.run_report import build_report, write_report
from ..core.runoff.curve_number import RunoffInputs
from ..core.runoff.rainfall import RainfallInputs

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def _rich_package(path):
    """Golden surface with curve number and rainfall blocks and chainages."""
    (cr, _, _, _, _), g = synthetic_run()
    shape = g["dem"].shape
    lc = np.where(np.arange(shape[1])[None, :] < shape[1] // 2, 10.0, 40.0) * np.ones(shape)
    ro = RunoffInputs(lc, hsg=np.full(shape, 2, np.int8))
    zones = np.where(np.arange(shape[1])[None, :] < shape[1] // 2, 1, 2) * np.ones(shape, np.int32)
    ri = RainfallInputs(zones.astype(np.int32), {1: "West", 2: "East"}, "synthetic",
                        np.full(shape, 900.0), "synthetic", "lo_1985")
    pts = [dict(x=a["outlet_x"], y=a["outlet_y"], fid=i, source_id=None, outlet_x=a["outlet_x"],
                outlet_y=a["outlet_y"], chainage=[300.0, 100.0, 200.0][i])
           for i, (_, a) in enumerate(cr)]
    c2, ca2, fp2, _, id_info = build_exchange_records(
        g["direction"], g["valid"], g["accum"], g["dem"], g["gt"], pts, snap_radius_cells=0,
        runoff=ro, rainfall=ri, channel_threshold_cells=100)
    md = dict(id_info, dem_path="synthetic", raw_dem_path="synthetic", dem_sha256="",
              dem_source="synthetic <surface> & test", flat_method="barnes",
              stream_threshold_km2="0.01", cell_size_m="10 x 10", qeht_version="test",
              runoff_json=ro.meta_json(), rainfall_json=ri.meta_json())
    write_exchange(path, c2, ca2, fp2, md, UTM37S_WKT, crs_name="WGS 84 / UTM zone 37S")
    return path


def test_characteristics():
    print("\n1. Catchment characteristics table")
    tmp = tempfile.mkdtemp()
    g = _golden_package(os.path.join(tmp, "golden.gpkg"))
    cols, rows = characteristics_rows(g)
    check("one row per crossing, keyed by outlet_uid",
          len(rows) == 3 and sorted(r["outlet_uid"] for r in rows) == ["X001", "X002", "X003"])
    core = [c for c, _ in CORE if c != "chainage_m"]
    check("core hydrology columns first, chainage dropped when no crossing has one",
          cols[:len(core)] == core and "chainage_m" not in cols)
    check("empty optional groups are left out (no soils, CN, rainfall, erosion)",
          not any(c in cols for c in ("soil_hsg", "cn_export", "rain_zone", "ero_a_mean_tha")))
    check("values joined from the right layer (area from catchments, 10-85 from flow paths)",
          all(r["area_km2"] > 0 and r["lfp_slope_1085"] > 0 for r in rows))
    cols2, rows2 = characteristics_rows(g, {"X001": {"flat_sensitive": 1, "flat_sensitivity_pct": 40.0}})
    check("extra values (flat-method check) join as their own group",
          "flat_sensitive" in cols2 and [r["flat_sensitive"] for r in rows2 if r["outlet_uid"] == "X001"] == [1])
    rich = _rich_package(os.path.join(tmp, "rich.gpkg"))
    cols3, rows3 = characteristics_rows(rich)
    check("curve number and rainfall groups appear when filled",
          all(c in cols3 for c in ("cn_export", "rain_zone", "map_mm", "rusle_r")))
    check("rows ordered along the chainage", [r["chainage_m"] for r in rows3] == [100.0, 200.0, 300.0])
    out, n = write_characteristics_csv(rich, os.path.join(tmp, "t", "chars.csv"))
    with open(out, newline="", encoding="utf-8") as f:
        data = list(csv.reader(f))
    check("CSV: header + one line per crossing, NULL as empty",
          n == 3 and len(data) == 4 and data[0] == cols3 and "" in data[1])


def test_report():
    print("\n2. Run report")
    tmp = tempfile.mkdtemp()
    rich = _rich_package(os.path.join(tmp, "rich.gpkg"))
    html = build_report(rich, "Demo <run>", {"warnings": ["w1 <b>"], "dem_audit": ["dup rows"],
                                             "flat_check": {"X001": {"flat_sensitive": 1,
                                                                     "flat_sensitivity_pct": 33.0,
                                                                     "area_barnes_km2": 0.3,
                                                                     "area_toward_km2": 0.2}}})
    check("self-contained page with the crossing ids and a plan",
          html.startswith("<!doctype html>") and all(u in html for u in ("X001", "X002", "X003"))
          and "<svg" in html and "<script" not in html and "http" not in html.split("<main>")[0])
    check("text is escaped (title, DEM source, warnings)",
          "Demo &lt;run&gt;" in html and "&lt;surface&gt; &amp; test" in html and "w1 &lt;b&gt;" in html)
    check("PROXY and ESTIMATE labels carried from the metadata",
          "PROXY" in html and "R ESTIMATE" in html and "Lo et al." in html)
    check("DEM audit and flat-method check shown", "dup rows" in html and "verify" in html)
    check("curve number lookup table appended", "Curve number lookup used" in html)
    g = _golden_package(os.path.join(tmp, "golden.gpkg"))
    plain = build_report(g)
    check("plain package: no rainfall / CN sections, audit 'not run', warnings 'None'",
          "Curve number lookup used" not in plain and "not run in this session" in plain
          and "<p>None.</p>" in plain)
    out = write_report(g, os.path.join(tmp, "r.html"))
    check("written to disk", os.path.getsize(out) > 3000)


def main(argv=None):
    test_characteristics()
    test_report()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: " + "; ".join(FAILURES))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
