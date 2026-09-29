# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""QEHT -> HEAS exchange package (schema qeht-heas-1).

One self-describing GeoPackage per run:

    crossings            Point        snapped outlets, with outlet_uid
    catchments           MultiPolygon one per crossing, morphometry
    flowpaths            LineString   longest flow path per crossing
    qeht_run_metadata    table        key/value provenance of the run
    qeht_field_dictionary table       every field: meaning, unit, method

Every feature in the three layers carries the SAME `outlet_uid` for one
crossing, so HEAS links them by key (link_method 'pour_point', confidence
1.0) - no spatial guessing, no shapefile/CSV round trip to misalign.

Pipeline (build_exchange_records) is pure NumPy and runs on bare Python;
the writer uses the stdlib GeoPackage writer in gpkg.py. The Processing
tool "Build HEAS exchange package" only reads the inputs and calls these.

Rules enforced here:
  * projected, metric CRS only (a geographic CRS is refused - lengths,
    areas and slopes would be in degrees);
  * outlet_uid unique and non-empty in every layer; every catchment and
    flow path refers to an existing crossing;
  * NaN is written as NULL, never as 0.
"""

import csv
import datetime as _dt
import hashlib
import json
import math
import os
import re

import numpy as np

from ..grid import NO_RECEIVER
from ..geometry.polygonize import mask_to_polygons
from ..linking.ids import assign_uids, validate_unique, OutletIdError
from ..watershed.delineate import (delineate_catchment, snap_pour_point,
                                   longest_flow_path)
from ..watershed.statistics import catchment_characteristics, horn_slope
from . import field_dictionary as fd
from .gpkg import GeoPackageWriter, read_table, list_tables, table_columns

SCHEMA_VERSION = fd.SCHEMA_VERSION


class ExchangeError(Exception):
    """The package cannot be written or does not satisfy the contract."""


# -- CRS -------------------------------------------------------------------

_GEOG = re.compile(r"^\s*(GEOGCS|GEOGCRS|GEODCRS|GEOGRAPHICCRS)\[", re.I)
_PROJ = re.compile(r"^\s*(PROJCS|PROJCRS|PROJECTEDCRS)\[", re.I)
_EPSG_TAIL = re.compile(
    r'(?:AUTHORITY\["EPSG",\s*"?(\d+)"?\]|ID\["EPSG",\s*(\d+)\])\s*\]\s*$', re.I)


def crs_check(wkt):
    """(ok, epsg_or_None, message) for a CRS given as WKT1/WKT2."""
    if not wkt or not str(wkt).strip():
        return False, None, ("The DEM has no coordinate reference system. Assign a "
                             "projected, metric CRS (e.g. UTM) before exporting.")
    w = str(wkt).strip()
    if _GEOG.match(w):
        return False, None, (
            "The DEM is in a geographic CRS (degrees). Areas, lengths and slopes "
            "in the exchange package must be metric: reproject the DEM to a "
            "projected CRS (e.g. WGS 84 / UTM zone 37S, EPSG:32737) and re-run.")
    if not _PROJ.match(w):
        return False, None, "Unrecognised CRS definition; a projected metric CRS is required."
    tail = w.split("PROJECTION")[-1] if "PROJECTION" in w.upper() else w
    if not re.search(r'UNIT\["(metre|meter|m)"', tail, re.I) and \
            not re.search(r'LENGTHUNIT\["(metre|meter)"', w, re.I):
        return False, None, "The projected CRS does not use metres; reproject to a metric CRS."
    m = _EPSG_TAIL.search(w)
    epsg = int(m.group(1) or m.group(2)) if m else None
    return True, epsg, ""


# -- provenance helpers ----------------------------------------------------

def file_fingerprint(path, hash_limit_bytes=512 * 1024 * 1024):
    """sha256 of a file, or 'size=..;mtime=..' when larger than the limit."""
    if not path or not os.path.exists(path):
        return ""
    st = os.stat(path)
    if st.st_size > hash_limit_bytes:
        return f"size={st.st_size};mtime={int(st.st_mtime)}"
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return "sha256:" + h.hexdigest()


def plugin_version():
    """Version string from the plugin's metadata.txt (or 'unknown')."""
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    try:
        with open(os.path.join(here, "metadata.txt"), encoding="utf-8") as f:
            for line in f:
                if line.strip().startswith("version="):
                    return line.split("=", 1)[1].strip()
    except OSError:
        pass
    return "unknown"


def utc_now():
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# -- pipeline --------------------------------------------------------------

def _cell_centre(gt, r, c):
    return (gt[0] + (c + 0.5) * gt[1] + (r + 0.5) * gt[2],
            gt[3] + (c + 0.5) * gt[4] + (r + 0.5) * gt[5])


def build_exchange_records(direction, valid, accumulation, elevation, geotransform,
                           outlets, snap_radius_cells=5, stream_mask=None,
                           local=False, stream_order=None, id_scheme="sequential",
                           id_prefix="X", id_width=3, id_start=1,
                           id_order="downstream", progress=None):
    """Run snap -> id -> catchment -> LFP -> characteristics for every outlet.

    Parameters
    ----------
    direction : internal 0-7 direction grid (NO_RECEIVER for none)
    valid : bool grid
    accumulation : accumulation in CELLS (as written by QEHT)
    elevation : the DEM for reported elevations/slopes (use the RAW DEM)
    geotransform : GDAL-style 6-tuple, north-up
    outlets : list of dicts {x, y, fid, source_id} in DEM CRS; `source_id`
        is the value of the user's ID attribute (or None). Optional keys:
        `outlet_x`, `outlet_y` - a fixed outlet location (e.g. from the
        crossing-candidate tool): used as-is, never snapped; `chainage` -
        chainage of the crossing (written to chainage_m, and used when
        id_order="chainage")
    stream_mask : bool grid for nearest-stream snapping (None -> snap to
        maximum accumulation, not recommended; see snap_pour_point)
    local : False = full upstream catchments (overlapping, the default for
        crossings); True = non-overlapping local catchments. In local mode
        outlets are processed upstream-first (ascending accumulation) so the
        result does not depend on the order of the input layer.
    stream_order : optional Strahler grid for the crossing attribute

    Returns (crossings, catchments, flowpaths, issues, id_info) where the
    first three are lists of (geometry, attributes) ready for the writer
    and `issues` is a list of human-readable warnings.
    """
    gt = tuple(geotransform)
    cw, ch = abs(gt[1]), abs(gt[5])
    cell_area = cw * ch
    rows, cols = direction.shape
    issues = []

    # 1. snap
    snapped = []
    for k, o in enumerate(outlets):
        fixed = o.get("outlet_x") is not None and o.get("outlet_y") is not None
        px, py = (o["outlet_x"], o["outlet_y"]) if fixed else (o["x"], o["y"])
        col = int(math.floor((px - gt[0]) / gt[1]))
        row = int(math.floor((py - gt[3]) / gt[5]))
        if not (0 <= row < rows and 0 <= col < cols):
            issues.append(f"Pour point {k + 1} (fid {o.get('fid')}) lies outside the DEM; skipped.")
            continue
        if fixed:
            r, c = row, col
        elif snap_radius_cells and snap_radius_cells > 0:
            r, c, _, _ = snap_pour_point(row, col, accumulation, valid,
                                         search_radius_cells=int(snap_radius_cells),
                                         stream_mask=stream_mask)
        else:
            r, c = row, col
        if not valid[r, c]:
            issues.append(f"Pour point {k + 1} (fid {o.get('fid')}) is on a NoData cell; skipped.")
            continue
        ox, oy = _cell_centre(gt, r, c)
        snapped.append(dict(o, row=r, col=c, outlet_x=ox, outlet_y=oy,
                            snap_dist_m=float(math.hypot(ox - o["x"], oy - o["y"])),
                            acc_cells=float(accumulation[r, c]), snapped_fixed=fixed))
    if not snapped:
        raise ExchangeError("No usable pour points. " + " ".join(issues))

    # two pour points snapping to the same cell would give identical
    # catchments under two IDs - stop and say which
    seen = {}
    for s in snapped:
        key = (s["row"], s["col"])
        if key in seen:
            raise ExchangeError(
                f"Pour points fid {seen[key]} and fid {s.get('fid')} snap to the same "
                f"cell (row {key[0]}, col {key[1]}). Remove one or reduce the snap radius.")
        seen[key] = s.get("fid")

    # 2. ids
    try:
        uids, id_info = assign_uids(
            len(snapped), scheme=id_scheme, prefix=id_prefix, width=id_width,
            start=id_start, order=id_order,
            accumulation=[s["acc_cells"] for s in snapped],
            chainage=([s.get("chainage") for s in snapped]
                      if all(s.get("chainage") is not None for s in snapped) else None),
            attribute_values=[s.get("source_id") for s in snapped])
    except OutletIdError as e:
        raise ExchangeError(str(e))
    for s, u in zip(snapped, uids):
        s["uid"] = u

    # 3. catchments
    slope_raster = horn_slope(elevation, valid, cw, ch)
    if local:
        order = sorted(range(len(snapped)), key=lambda i: (snapped[i]["acc_cells"], i))
        labels = delineate_catchment(direction, valid,
                                     [(snapped[i]["row"], snapped[i]["col"]) for i in order])
        label_of = {i: n + 1 for n, i in enumerate(order)}

    crossings, catchments, flowpaths = [], [], []
    link = {"link_method": "pour_point", "link_confidence": 1.0, "link_note": ""}
    n = len(snapped)
    for i, s in enumerate(snapped):
        if progress is not None:
            progress(i / n, f"Outlet {s['uid']}")
        mask = (labels == label_of[i]) if local else \
            (delineate_catchment(direction, valid, [(s["row"], s["col"])]) > 0)

        so = None
        if stream_order is not None:
            v = stream_order[s["row"], s["col"]]
            so = int(v) if np.isfinite(v) and v > 0 else None
        crossings.append(((s["outlet_x"], s["outlet_y"]), dict(
            link, outlet_uid=s["uid"], crossing_id=s["uid"], outlet_id=s.get("fid"),
            source_id=s.get("source_id"), input_x=s["x"], input_y=s["y"],
            outlet_x=s["outlet_x"], outlet_y=s["outlet_y"], snap_dist_m=s["snap_dist_m"],
            acc_at_outlet_km2=(s["acc_cells"] + 1.0) * cell_area / 1.0e6,
            stream_order=so, chainage_m=s.get("chainage"))))

        lfp = longest_flow_path(direction, valid, (s["row"], s["col"]),
                                elevation=elevation, cell_width=cw, cell_height=ch,
                                catchment_mask=mask)
        chs = catchment_characteristics(mask, elevation, valid, cw, ch,
                                        flow_path=lfp, slope_raster=slope_raster)
        if not chs:
            issues.append(f"{s['uid']}: empty catchment (no valid cells); no catchment "
                          "or flow path written.")
            continue
        polys = mask_to_polygons(mask, gt)
        catchments.append((polys, dict(
            chs, **link, outlet_uid=s["uid"], catchment_id=s["uid"],
            outlet_id=s.get("fid"), catchment_mode="local" if local else "full")))

        if len(lfp["cells"]) >= 2:
            line = [_cell_centre(gt, r, c) for r, c in lfp["cells"]]
            flowpaths.append((line, dict(
                chs, **link, outlet_uid=s["uid"], flowpath_id=s["uid"],
                outlet_id=s.get("fid"))))
        else:
            issues.append(f"{s['uid']}: catchment is a single cell; no flow path written.")
    if progress is not None:
        progress(1.0, "Characteristics done")
    return crossings, catchments, flowpaths, issues, id_info


# -- writer ----------------------------------------------------------------

def _check_links(crossings, catchments, flowpaths):
    for name, layer in (("crossings", crossings), ("catchments", catchments),
                        ("flowpaths", flowpaths)):
        try:
            validate_unique([a.get("outlet_uid") for _, a in layer],
                            what=f"outlet_uid in {name}")
        except OutletIdError as e:
            raise ExchangeError(str(e))
    known = {a["outlet_uid"] for _, a in crossings}
    for name, layer in (("catchments", catchments), ("flowpaths", flowpaths)):
        orphans = sorted({a["outlet_uid"] for _, a in layer} - known)
        if orphans:
            raise ExchangeError(f"{name} refer to crossings that do not exist: "
                                f"{', '.join(orphans[:20])}")


def write_exchange(path, crossings, catchments, flowpaths, metadata, crs_wkt,
                   crs_epsg=None, crs_name="", timestamp=None, csv_dir=None,
                   extra_layers=None, extra_tables=None):
    """Write the exchange GeoPackage (and optional per-layer CSVs).

    `metadata` is a dict of run-metadata values (see field_dictionary.
    METADATA_KEYS); schema_version, qeht_version, run_utc, crs_epsg and
    n_crossings are filled in here. `timestamp` fixes run_utc and the
    GeoPackage last_change (reproducible fixtures); default = now.
    `extra_layers`: optional list of (name, geom_type, fields, rows, desc),
    e.g. crossing_candidates / road_alignment / parallel_reaches.
    `extra_tables`: optional list of (name, fields, rows, desc), e.g.
    renumber_log. Both are additive under qeht-heas-1.
    Returns the path.
    """
    ok, epsg_from_wkt, msg = crs_check(crs_wkt)
    if not ok:
        raise ExchangeError(msg)
    epsg = crs_epsg or epsg_from_wkt
    _check_links(crossings, catchments, flowpaths)

    run_utc = timestamp or utc_now()
    md = {k: "" for k, _ in fd.METADATA_KEYS}
    md.update({k: ("" if v is None else v) for k, v in (metadata or {}).items()})
    md.update({"schema_version": SCHEMA_VERSION,
               "qeht_version": md.get("qeht_version") or plugin_version(),
               "run_utc": run_utc,
               "crs_epsg": str(epsg) if epsg else "",
               "crs_name": crs_name or md.get("crs_name", ""),
               "n_crossings": str(len(crossings)),
               "lfp_1085_reference": "outlet"})
    if isinstance(md.get("parameters_json"), (dict, list)):
        md["parameters_json"] = json.dumps(md["parameters_json"], sort_keys=True,
                                           default=str)

    srs_id = int(epsg) if epsg else 99999
    w = GeoPackageWriter(path, overwrite=True,
                         timestamp=run_utc.replace("Z", ".000Z") if len(run_utc) == 20 else run_utc)
    try:
        w.add_srs(srs_id, "EPSG" if epsg else "NONE", srs_id if epsg else srs_id,
                  crs_wkt, crs_name or (f"EPSG:{epsg}" if epsg else "project CRS"))
        for table, data in (("crossings", crossings), ("catchments", catchments),
                            ("flowpaths", flowpaths)):
            gtype, fields, desc = fd.LAYERS[table]
            w.add_feature_table(table, gtype, srs_id,
                                [(f[0], f[1]) for f in fields], data, description=desc)
        for name, gtype, fields, rows, desc in (extra_layers or []):
            w.add_feature_table(name, gtype, srs_id, fields, rows, description=desc)
        for name, fields, rows, desc in (extra_tables or []):
            w.add_attribute_table(name, fields, rows, description=desc)
        w.add_attribute_table(
            "qeht_run_metadata", [("key", "text"), ("value", "text")],
            [{"key": k, "value": str(md.get(k, ""))} for k, _ in fd.METADATA_KEYS]
            + [{"key": k, "value": str(v)} for k, v in sorted(md.items())
               if k not in dict(fd.METADATA_KEYS)],
            description="Run provenance (key/value), schema " + SCHEMA_VERSION)
        w.add_attribute_table(
            "qeht_field_dictionary",
            [(c, "text") for c in fd.FIELD_DICTIONARY_COLUMNS],
            [dict(zip(fd.FIELD_DICTIONARY_COLUMNS, r)) for r in fd.dictionary_rows()],
            description="Meaning, unit and method of every field")
    finally:
        w.close()

    if csv_dir:
        write_csvs(path, csv_dir)
    return path


def write_csvs(gpkg_path, csv_dir):
    """One CSV per layer/table (attributes only) for spreadsheet users.

    The GeoPackage remains the preferred import; CSVs carry outlet_uid so
    they can always be re-joined by key, never by row order.
    """
    os.makedirs(csv_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(gpkg_path))[0]
    written = []
    for table, _ in list_tables(gpkg_path):
        cols = [c for c, _ in table_columns(gpkg_path, table) if c not in ("fid", "geom")]
        rows = read_table(gpkg_path, table, with_geometry=False)
        out = os.path.join(csv_dir, f"{stem}_{table}.csv")
        with open(out, "w", newline="", encoding="utf-8") as f:
            wr = csv.writer(f)
            wr.writerow(cols)
            for r in rows:
                wr.writerow(["" if r.get(c) is None else r.get(c) for c in cols])
        written.append(out)
    return written


# -- validator (the same checks HEAS applies on import) --------------------

def validate_exchange(path):
    """Check a package against the contract. Returns (errors, warnings)."""
    errors, warnings = [], []
    tables = dict(list_tables(path))
    for t in ("crossings", "catchments", "flowpaths", "qeht_run_metadata",
              "qeht_field_dictionary"):
        if t not in tables:
            errors.append(f"missing table {t}")
    if errors:
        return errors, warnings
    md = {r["key"]: r["value"] for r in read_table(path, "qeht_run_metadata")}
    sv = md.get("schema_version", "")
    m = re.match(r"^qeht-heas-(\d+)$", sv or "")
    if not m:
        errors.append(f"unrecognised schema_version '{sv}'")
    elif int(m.group(1)) != fd.SCHEMA_MAJOR:
        errors.append(f"schema {sv} is not supported by this reader (expects "
                      f"qeht-heas-{fd.SCHEMA_MAJOR}); update the reader")
    uids = {}
    for layer in ("crossings", "catchments", "flowpaths"):
        cols = {c for c, _ in table_columns(path, layer)}
        missing = [f for f in fd.field_names(layer) if f not in cols]
        if missing:
            errors.append(f"{layer}: missing fields {', '.join(missing)}")
        vals = [r.get("outlet_uid") for r in read_table(path, layer, with_geometry=False)]
        try:
            validate_unique(vals, what=f"outlet_uid in {layer}")
        except OutletIdError as e:
            errors.append(str(e))
        uids[layer] = set(v for v in vals if v)
    for layer in ("catchments", "flowpaths"):
        orphan = uids[layer] - uids["crossings"]
        if orphan:
            errors.append(f"{layer} without a crossing: {', '.join(sorted(orphan))}")
        lonely = uids["crossings"] - uids[layer]
        if lonely:
            warnings.append(f"crossings without a {layer[:-1]}: {', '.join(sorted(lonely))}")
    return errors, warnings
