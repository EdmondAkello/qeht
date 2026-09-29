# -*- coding: utf-8 -*-
"""v0.9 interop validation: outlet_uid, D1 slope domains, exchange package.

Runs on bare Python + NumPy (no QGIS, no GDAL), like test_core:

    python -m qeht.tests.test_interop            run the checks
    python -m qeht.tests.test_interop --regen    rewrite the golden fixture

The golden fixture (tests/fixtures/golden_exchange.gpkg + .json) is the
shared QEHT <-> HEAS contract test: QEHT asserts it reproduces the same
attributes; HEAS asserts it imports the file with no field mapping.
"""

import json
import math
import os
import sqlite3
import sys
import tempfile

import numpy as np

from ..core.conditioning.fill import fill_depressions
from ..core.flow.direction import d8_direction
from ..core.flow.accumulation import flow_accumulation, strahler_order
from ..core.watershed.delineate import extract_streams, delineate_catchment
from ..core.watershed.statistics import (path_10_85, path_slope_10_85,
                                         path_slope_10_85_v083,
                                         catchment_characteristics)
from ..core.linking.ids import (assign_uids, sequential_uids, attribute_uids,
                                OutletIdError)
from ..core.geometry.polygonize import mask_to_polygons, polygons_area
from ..core.interop import gpkg
from ..core.interop import field_dictionary as fd
from ..core.interop.heas_exchange import (build_exchange_records, write_exchange,
                                          validate_exchange, crs_check,
                                          ExchangeError, write_csvs)

FAILURES = []
HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE_GPKG = os.path.join(HERE, "fixtures", "golden_exchange.gpkg")
FIXTURE_JSON = os.path.join(HERE, "fixtures", "golden_exchange.json")
GOLDEN_TIME = "2026-09-29T00:00:00Z"

UTM37S_WKT = (
    'PROJCS["WGS 84 / UTM zone 37S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID['
    '"WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG",'
    '"6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",'
    '0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],'
    'PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],'
    'PARAMETER["central_meridian",39],PARAMETER["scale_factor",0.9996],'
    'PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],'
    'UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],'
    'AXIS["Northing",NORTH],AUTHORITY["EPSG","32737"]]')


def check(name, condition, detail=""):
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def raises(fn, exc):
    try:
        fn()
    except exc as e:
        return str(e)
    return None


# ------------------------------------------------------------------
def test_1085_convention():
    print("\n1. 10-85 slope: outlet-referenced (v0.9) vs divide-referenced (<=0.8.3)")
    # Straight N-S path of 21 cells (20 steps x 10 m = 200 m), concave-up
    # profile z = a * d^2 with d = distance FROM THE OUTLET. Nodes fall
    # exactly on 10% and 85% (2 and 17 steps), so both conventions have
    # closed-form answers:
    #   outlet-referenced : a L (0.85^2 - 0.10^2) / 0.75 = 0.95 a L
    #   divide-referenced : a L (0.90^2 - 0.15^2) / 0.75 = 1.05 a L
    a, cs, N = 0.002, 10.0, 20
    L = N * cs
    cells = [(r, 0) for r in range(N + 1)]            # divide (row 0) -> outlet
    z = np.zeros((N + 1, 1))
    for r in range(N + 1):
        d = (N - r) * cs
        z[r, 0] = 100.0 + a * d * d
    res = path_10_85(cells, z, cs, cs)
    check("outlet-referenced slope = 0.95 a L", math.isclose(res["slope"], 0.95 * a * L,
          rel_tol=1e-12), f"{res['slope']:.6f} vs {0.95 * a * L:.6f}")
    check("L10 = 0.10 L and L85 = 0.85 L from the outlet",
          math.isclose(res["L10_m"], 20.0) and math.isclose(res["L85_m"], 170.0))
    check("z10, z85 are the profile at those distances",
          math.isclose(res["z10_m"], 100 + a * 20 ** 2)
          and math.isclose(res["z85_m"], 100 + a * 170 ** 2),
          f"z10={res['z10_m']:.3f} z85={res['z85_m']:.3f}")
    div = path_10_85(cells, z, cs, cs, reference="divide")["slope"]
    check("divide-referenced slope = 1.05 a L", math.isclose(div, 1.05 * a * L, rel_tol=1e-12))
    check("the two conventions differ on a concave profile", abs(div - res["slope"]) > 1e-4,
          f"divide {div:.5f} vs outlet {res['slope']:.5f} ({100 * (div / res['slope'] - 1):+.1f}%)")
    old = path_slope_10_85_v083(cells, z, cs, cs)
    check("v0.8.3 function reproduces the divide-referenced value", math.isclose(old, div,
          rel_tol=1e-12))
    check("path_slope_10_85() now returns the outlet-referenced value",
          math.isclose(path_slope_10_85(cells, z, cs, cs), res["slope"]))

    # Linear interpolation: on a uniform-gradient path the 10-85 slope equals
    # the gradient exactly, whatever the cell count and diagonal steps.
    cells = [(k, k) for k in range(8)]                 # diagonal, 7 steps
    zz = np.full((8, 8), np.nan)
    step = cs * math.sqrt(2)
    for k in range(8):
        zz[k, k] = 50.0 + 0.03 * (7 - k) * step
    s = path_10_85(cells, zz, cs, cs)["slope"]
    check("uniform gradient recovered exactly off-node (diagonal path)",
          math.isclose(s, 0.03, rel_tol=1e-12), f"{s:.6f}")
    check("path shorter than 2 cells -> NaN", math.isnan(path_10_85([(0, 0)], zz, cs, cs)["slope"]))


def test_slope_domains():
    print("\n2. Four slope domains, new names + legacy aliases")
    n = 21
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    dem = 100.0 + np.abs(xx - n // 2) * 2.0 - yy * 1.0
    valid = np.ones((n, n), bool)
    filled, _, _ = fill_depressions(dem, valid, cell_width=30, cell_height=30)
    direction, _ = d8_direction(filled, valid, 30.0, 30.0)
    from ..core.watershed.delineate import longest_flow_path
    outlet = (n - 1, n // 2)
    mask = delineate_catchment(direction, valid, [outlet]) > 0
    lfp = longest_flow_path(direction, valid, outlet, elevation=dem,
                            cell_width=30.0, cell_height=30.0, catchment_mask=mask)
    ch = catchment_characteristics(mask, dem, valid, 30.0, 30.0, flow_path=lfp)
    for k in ("catch_slope_horn", "catch_relief_ratio", "lfp_slope", "lfp_slope_1085",
              "lfp_L10_m", "lfp_L85_m", "lfp_z10_m", "lfp_z85_m"):
        check(f"{k} present and finite", k in ch and np.isfinite(ch[k]))
    check("slope_mean alias == catch_slope_horn", ch["slope_mean"] == ch["catch_slope_horn"])
    check("slope_relief_ratio alias == catch_relief_ratio",
          ch["slope_relief_ratio"] == ch["catch_relief_ratio"])
    check("relief ratio = relief / LFP length",
          math.isclose(ch["catch_relief_ratio"], ch["relief_m"] / ch["lfp_length_m"]))
    check("hand check: (z85 - z10)/(L85 - L10) = lfp_slope_1085",
          math.isclose((ch["lfp_z85_m"] - ch["lfp_z10_m"]) / (ch["lfp_L85_m"] - ch["lfp_L10_m"]),
                       ch["lfp_slope_1085"], rel_tol=1e-12))


def test_ids():
    print("\n3. outlet_uid generation")
    u = sequential_uids(3, prefix="X", order="downstream", accumulation=[10, 500, 40])
    check("downstream-first: largest area gets X001", u == ["X003", "X001", "X002"], str(u))
    check("input order", sequential_uids(3, prefix="KE-", order="input") ==
          ["KE-001", "KE-002", "KE-003"])
    big = sequential_uids(1200, prefix="X", order="input")
    check("width grows past 999 (no ambiguous ids)", big[0] == "X0001" and big[-1] == "X1200")
    check("ties keep input order", sequential_uids(3, order="downstream",
          accumulation=[5, 5, 5]) == ["X001", "X002", "X003"])
    check("chainage order", sequential_uids(3, order="chainage", chainage=[300, 100, 200])
          == ["X003", "X001", "X002"])
    msg = raises(lambda: attribute_uids(["C1", "C2", "C1", "C3", "C3"]), OutletIdError)
    check("duplicate attribute ids stop the run and are listed",
          msg is not None and "C1" in msg and "C3" in msg, msg or "")
    msg = raises(lambda: attribute_uids(["C1", None, "NULL"]), OutletIdError)
    check("empty / NULL attribute ids are refused", msg is not None and "empty" in msg)
    uids, info = assign_uids(2, scheme="attribute", attribute_values=[" 17 ", 18], prefix="CV")
    check("attribute ids stripped and prefixed", uids == ["CV17", "CV18"], str(uids))
    check("id metadata recorded", info["id_scheme"] == "attribute")


def test_polygonize():
    print("\n4. Mask -> polygon (pure NumPy)")
    gt = (500000.0, 10.0, 0.0, 9900000.0, 0.0, -10.0)
    m = np.ones((5, 5), bool); m[2, 2] = False
    p = mask_to_polygons(m, gt)
    check("donut: one polygon with one hole", len(p) == 1 and len(p[0]) == 2)
    check("donut area = 24 cells", math.isclose(polygons_area(p), 24 * 100.0))
    m = np.zeros((3, 3), bool); m[0, 0] = m[1, 1] = True
    check("diagonal cells -> 2 parts (4-connected)", len(mask_to_polygons(m, gt)) == 2)
    m = np.zeros((4, 4), bool); m[0, :] = m[:, 0] = True; m[1, 1] = True; m[2, 2] = True
    p = mask_to_polygons(m, gt)
    check("pinched notch -> area preserved", math.isclose(polygons_area(p), m.sum() * 100.0))
    rng = np.random.default_rng(7)
    bad = 0
    for _ in range(300):
        m = rng.random((rng.integers(1, 20), rng.integers(1, 20))) < 0.55
        if m.any() and not math.isclose(polygons_area(mask_to_polygons(m, gt)), m.sum() * 100.0):
            bad += 1
    check("300 random masks: area = cells x cell area", bad == 0, f"{bad} mismatches")
    p = mask_to_polygons(np.ones((2, 3), bool), gt)
    ring = p[0][0]
    check("outer ring closed, counter-clockwise, 4 corners",
          ring[0] == ring[-1] and len(ring) == 5 and polygons_area(p) > 0)
    check("corners on the grid", ring[0] == (500000.0, 9900000.0) or
          (500000.0, 9900000.0) in ring)


def test_gpkg_roundtrip():
    print("\n5. GeoPackage writer (stdlib)")
    tmp = tempfile.mkdtemp()
    path = os.path.join(tmp, "t.gpkg")
    w = gpkg.GeoPackageWriter(path, timestamp="2026-09-29T00:00:00.000Z")
    w.add_srs(32737, "EPSG", 32737, UTM37S_WKT, "WGS 84 / UTM zone 37S")
    w.add_feature_table("pts", "POINT", 32737, [("name", "text"), ("v", "real"), ("n", "int")],
                        [((1.5, 2.5), {"name": "a", "v": float("nan"), "n": 3}),
                         ((3.0, 4.0), {"name": None, "v": 2.25, "n": None})])
    w.add_feature_table("poly", "MULTIPOLYGON", 32737, [("k", "int")],
                        [([[[(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)]]], {"k": 1})])
    w.close()
    con = sqlite3.connect(path)
    appid = con.execute("PRAGMA application_id").fetchone()[0]
    uv = con.execute("PRAGMA user_version").fetchone()[0]
    con.close()
    check("application_id = 'GPKG'", appid == 0x47504B47)
    check("user_version = 10200", uv == 10200)
    rows = gpkg.read_table(path, "pts")
    check("NaN written as NULL (not 0)", rows[0]["v"] is None and rows[1]["v"] == 2.25)
    check("point geometry round-trips", rows[1]["_geom"] == (3.0, 4.0))
    poly = gpkg.read_table(path, "poly")[0]["_geom"]
    check("multipolygon round-trips", poly[0][0][2] == (10.0, 10.0))
    tabs = dict(gpkg.list_tables(path))
    check("gpkg_contents lists both tables", tabs.get("pts") == "features" and "poly" in tabs)


def test_crs_rule():
    print("\n6. CRS rule: projected metric only")
    ok, epsg, _ = crs_check(UTM37S_WKT)
    check("UTM 37S accepted, EPSG parsed", ok and epsg == 32737, str(epsg))
    ok, _, msg = crs_check(gpkg.WGS84_WKT)
    check("geographic CRS refused with a clear message", not ok and "geographic" in msg)
    ok, _, _ = crs_check("")
    check("missing CRS refused", not ok)


# -- synthetic catchment set used by the pipeline and golden fixture -------

def synthetic_case():
    """60 x 50 cells of 10 m in UTM 37S: a main valley draining south with
    a tributary from the east. Three crossings: X on the main valley near
    the bottom, a second on the main valley upstream (nested inside the
    first), and one on the tributary."""
    rows, cols, cs = 60, 50, 10.0
    r, c = np.mgrid[0:rows, 0:cols].astype(float)
    main = 200.0 + 0.8 * (rows - 1 - r) + 1.5 * np.abs(c - 20)          # valley on col 20
    trib = 200.0 + 0.8 * (rows - 1 - 25) + 0.6 * (c - 20) + 1.2 * np.abs(r - 25)  # row 25, from E
    dem = np.where(c > 20, np.minimum(main, trib), main)
    dem += 0.001 * c                                     # break exact ties deterministically
    valid = np.ones((rows, cols), bool)
    gt = (500000.0, cs, 0.0, 9900600.0, 0.0, -cs)
    return dem, valid, gt, cs


def synthetic_run(local=False, pour_order=None):
    dem, valid, gt, cs = synthetic_case()
    filled, _, _ = fill_depressions(dem, valid, cell_width=cs, cell_height=cs)
    direction, _ = d8_direction(filled, valid, cs, cs)
    accum, _ = flow_accumulation(direction, valid)
    streams = extract_streams(accum, valid, threshold_cells=100)
    order = strahler_order(direction, valid, streams)

    def xy(row, col, dx=0.0, dy=0.0):     # a point near a cell, off-channel
        return gt[0] + (col + 0.5) * cs + dx, gt[3] - (row + 0.5) * cs + dy
    pts = [  # (row, col, offset): deliberately off the channel by a few metres
        dict(zip(("x", "y"), xy(20, 20, 12.0, 3.0)), fid=7, source_id=None),   # upstream main
        dict(zip(("x", "y"), xy(55, 20, -14.0, -2.0)), fid=3, source_id=None),  # downstream main
        dict(zip(("x", "y"), xy(25, 38, 4.0, 11.0)), fid=11, source_id=None),  # tributary
    ]
    if pour_order is not None:
        pts = [pts[i] for i in pour_order]
    recs = build_exchange_records(direction, valid, accum, dem, gt, pts,
                                  snap_radius_cells=3, stream_mask=streams,
                                  local=local, stream_order=order, id_prefix="X")
    return recs, dict(direction=direction, valid=valid, accum=accum, dem=dem, gt=gt, cs=cs)


def test_pipeline():
    print("\n7. Pipeline: crossing -> catchment -> LFP, one outlet_uid")
    (cr, ca, fp, issues, id_info), g = synthetic_run()
    check("3 crossings, 3 catchments, 3 flow paths", (len(cr), len(ca), len(fp)) == (3, 3, 3),
          f"{len(cr)}/{len(ca)}/{len(fp)} issues={issues}")
    uid = lambda layer: sorted(a["outlet_uid"] for _, a in layer)
    check("same outlet_uid set on all three layers", uid(cr) == uid(ca) == uid(fp), str(uid(cr)))
    by = {a["outlet_uid"]: a for _, a in cr}
    big = max(by.values(), key=lambda a: a["acc_at_outlet_km2"])
    check("X001 is the most downstream (largest area)", big["outlet_uid"] == "X001")
    check("legacy outlet_id kept as the input fid", sorted(a["outlet_id"] for _, a in cr) == [3, 7, 11])
    check("every crossing snapped onto the channel", all(a["snap_dist_m"] > 0 for a in by.values()))
    # the v0.8.3 overlap defect: each polygon must hold its FULL area
    cell = g["cs"] ** 2
    ok = all(math.isclose(polygons_area(geom), a["area_km2"] * 1e6, rel_tol=1e-9) for geom, a in ca)
    check("full mode: every polygon area == its area_km2 (overlap-safe)", ok)
    area = {a["outlet_uid"]: a["area_km2"] for _, a in ca}
    acc = {a["outlet_uid"]: a["acc_at_outlet_km2"] for _, a in cr}
    check("acc_at_outlet_km2 == area_km2 in full mode",
          all(math.isclose(area[k], acc[k]) for k in area))
    nested_ok = area["X001"] > area["X002"] + area["X003"] - 1e-12
    check("downstream catchment contains the upstream ones", nested_ok,
          ", ".join(f"{k}={v:.4f}" for k, v in sorted(area.items())))
    so = {a["outlet_uid"]: a["stream_order"] for _, a in cr}
    check("stream order recorded at crossings", all(v and v >= 1 for v in so.values()), str(so))
    # order-independence of local mode
    (_, ca_l, _, _, _), _ = synthetic_run(local=True)
    (_, ca_l2, _, _, _), _ = synthetic_run(local=True, pour_order=[2, 0, 1])
    al = {a["outlet_uid"]: a["area_km2"] for _, a in ca_l}
    al2 = {a["outlet_uid"]: a["area_km2"] for _, a in ca_l2}
    check("local mode: areas sum to the downstream full area",
          math.isclose(sum(al.values()), area["X001"]), f"{sum(al.values()):.4f} vs {area['X001']:.4f}")
    check("local mode is independent of pour-point order", al == al2)
    # duplicate snap target is refused
    dem, valid, gt, cs = synthetic_case()
    msg = raises(lambda: build_exchange_records(
        g["direction"], g["valid"], g["accum"], g["dem"], g["gt"],
        [dict(x=gt[0] + 205, y=gt[3] - 555, fid=1, source_id=None),
         dict(x=gt[0] + 206, y=gt[3] - 556, fid=2, source_id=None)],
        snap_radius_cells=3, stream_mask=extract_streams(g["accum"], valid, threshold_cells=100)),
        ExchangeError)
    check("two points snapping to one cell are refused", msg is not None and "same" in msg, msg or "")


def _golden_package(path):
    (cr, ca, fp, issues, id_info), g = synthetic_run()
    md = dict(id_info, dem_path="synthetic", raw_dem_path="synthetic", dem_sha256="",
              dem_source="synthetic test surface (tests/test_interop.py)",
              conditioning="priority-flood fill, min_slope 0", flat_method="toward",
              tie_rule="lowest internal index (E, SE, S ...)",
              stream_threshold_cells="100", stream_threshold_km2=str(100 * 100 / 1e6),
              snap_radius_cells="3", snap_strategy="nearest_stream",
              catchment_mode="full", qeht_version="golden",
              cell_size_m="10 x 10", parameters_json={"fixture": "golden", "rows": 60, "cols": 50})
    write_exchange(path, cr, ca, fp, md, UTM37S_WKT, crs_name="WGS 84 / UTM zone 37S",
                   timestamp=GOLDEN_TIME)
    return path


def _summary(path):
    """Attribute summary compared against the committed JSON (rounded)."""
    def rnd(v):
        return round(v, 6) if isinstance(v, float) else v
    out = {}
    for layer in ("crossings", "catchments", "flowpaths"):
        rows = gpkg.read_table(path, layer)
        out[layer] = [{k: rnd(v) for k, v in r.items() if k not in ("fid", "_geom")}
                      for r in rows]
        out[layer + "_vertices"] = [len(r["_geom"]) if layer != "catchments"
                                    else sum(len(ring) for p in r["_geom"] for ring in p)
                                    for r in rows]
    out["metadata"] = {r["key"]: r["value"] for r in gpkg.read_table(path, "qeht_run_metadata")}
    return out


def test_exchange_package(regen=False):
    print("\n8. Exchange GeoPackage + golden fixture")
    tmp = tempfile.mkdtemp()
    path = _golden_package(os.path.join(tmp, "golden_exchange.gpkg"))
    errors, warnings = validate_exchange(path)
    check("validator: no errors", not errors, "; ".join(errors))
    check("validator: no warnings", not warnings, "; ".join(warnings))
    tabs = dict(gpkg.list_tables(path))
    check("five tables", set(tabs) == {"crossings", "catchments", "flowpaths",
                                       "qeht_run_metadata", "qeht_field_dictionary"})
    for layer in ("crossings", "catchments", "flowpaths"):
        cols = [c for c, _ in gpkg.table_columns(path, layer)][2:]
        check(f"{layer}: fields exactly as the field dictionary, in order",
              cols == fd.field_names(layer))
    md = {r["key"]: r["value"] for r in gpkg.read_table(path, "qeht_run_metadata")}
    check("schema_version qeht-heas-1", md["schema_version"] == "qeht-heas-1")
    check("crs_epsg 32737 recorded", md["crs_epsg"] == "32737")
    check("lfp_1085_reference = outlet", md["lfp_1085_reference"] == "outlet")
    nd = len(gpkg.read_table(path, "qeht_field_dictionary"))
    check("field dictionary complete", nd == len(fd.dictionary_rows()), f"{nd} rows")
    # geographic refusal at write time
    (cr, ca, fp, _, _), _ = synthetic_run()
    msg = raises(lambda: write_exchange(os.path.join(tmp, "g.gpkg"), cr, ca, fp, {},
                                        gpkg.WGS84_WKT), ExchangeError)
    check("writer refuses a geographic CRS", msg is not None)
    # schema major version gate
    bad = os.path.join(tmp, "v2.gpkg"); _golden_package(bad)
    con = sqlite3.connect(bad)
    con.execute("UPDATE qeht_run_metadata SET value='qeht-heas-2' WHERE key='schema_version'")
    con.commit(); con.close()
    errors, _ = validate_exchange(bad)
    check("unknown major schema version rejected", any("not supported" in e for e in errors))
    # CSV sidecars keep the key
    csvs = write_csvs(path, os.path.join(tmp, "csv"))
    with open([c for c in csvs if c.endswith("_flowpaths.csv")][0], encoding="utf-8") as f:
        header = f.readline().strip().split(",")
    check("CSV export carries outlet_uid first", header[0] == "outlet_uid" and len(csvs) == 5)

    summary = _summary(path)
    if regen:
        os.makedirs(os.path.dirname(FIXTURE_GPKG), exist_ok=True)
        _golden_package(FIXTURE_GPKG)
        with open(FIXTURE_JSON, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=1, sort_keys=True)
        print("  [INFO] golden fixture regenerated")
    have = os.path.exists(FIXTURE_JSON) and os.path.exists(FIXTURE_GPKG)
    check("golden fixture present", have)
    if have:
        with open(FIXTURE_JSON, encoding="utf-8") as f:
            ref = json.load(f)
        check("attributes reproduce the golden fixture", summary == ref,
              "" if summary == ref else "regenerate only if the change is intended")
        check("committed golden .gpkg matches the JSON",
              _summary(FIXTURE_GPKG) == ref)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    print("=" * 62)
    print("QEHT INTEROP - v0.9 (outlet_uid, D1 slopes, HEAS exchange)")
    print("=" * 62)
    test_1085_convention()
    test_slope_domains()
    test_ids()
    test_polygonize()
    test_gpkg_roundtrip()
    test_crs_rule()
    test_pipeline()
    test_exchange_package(regen="--regen" in argv)
    print("\n" + "=" * 62)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S):")
        for f in FAILURES:
            print("   -", f)
        return 1
    print("ALL CHECKS PASSED")
    print("=" * 62)
    return 0


if __name__ == "__main__":
    sys.exit(main())
