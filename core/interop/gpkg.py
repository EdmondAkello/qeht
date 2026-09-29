# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Minimal OGC GeoPackage (1.2) writer and reader - Python stdlib only.

Why not OGR
-----------
The exchange package is a contract between two codebases. Writing it with
`sqlite3` + `struct` means:
  * it can be produced and verified on bare Python (tests need no GDAL);
  * field names, types, order and NULLs are exactly what the field
    dictionary says, with no driver-side renaming or type widening;
  * the same inputs give the same file on every platform (golden fixture).
The output is validated against GDAL/OGR in the test suite and opens in
QGIS like any other GeoPackage.

Supported: POINT, LINESTRING, MULTIPOLYGON (2D) feature tables and
attribute (non-spatial) tables. No spatial index is written; QGIS and GDAL
do not need one for layers of this size.

No QGIS imports, no GDAL.
"""

import math
import os
import re
import sqlite3
import struct

GPKG_APPLICATION_ID = 0x47504B47          # "GPKG"
GPKG_USER_VERSION = 10200                 # GeoPackage 1.2

_WKB_TYPES = {"POINT": 1, "LINESTRING": 2, "POLYGON": 3, "MULTIPOLYGON": 6}
_SQL_TYPES = {"text": "TEXT", "int": "INTEGER", "real": "REAL"}

WGS84_WKT = (
    'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,'
    'AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,'
    'AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,'
    'AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],'
    'AUTHORITY["EPSG","4326"]]')


class GpkgError(Exception):
    pass


_IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]{0,62}$")


def _ident(name):
    """Validate a table/column name before it is placed in SQL.

    SQLite cannot bind identifiers as parameters, so names are checked
    against a strict pattern (letters, digits, underscore) and quoted;
    VALUES are always bound with '?'.
    """
    if not isinstance(name, str) or not _IDENT.match(name):
        raise GpkgError(f"invalid table/field name {name!r}")
    return '"' + name + '"'


# -- WKB -------------------------------------------------------------------

def _wkb_ring(ring):
    out = [struct.pack("<I", len(ring))]
    out += [struct.pack("<dd", float(x), float(y)) for x, y in ring]
    return b"".join(out)


def wkb_point(x, y):
    return struct.pack("<BIdd", 1, 1, float(x), float(y))


def wkb_linestring(coords):
    return struct.pack("<BI", 1, 2) + _wkb_ring(coords)


def wkb_polygon(rings):
    return struct.pack("<BII", 1, 3, len(rings)) + b"".join(_wkb_ring(r) for r in rings)


def wkb_multipolygon(polys):
    return (struct.pack("<BII", 1, 6, len(polys))
            + b"".join(wkb_polygon(p) for p in polys))


def _envelope(coords_iter):
    xs, ys = [], []
    for x, y in coords_iter:
        xs.append(float(x)); ys.append(float(y))
    if not xs:
        return None
    return (min(xs), max(xs), min(ys), max(ys))


def gpkg_blob(geom_type, geom, srs_id):
    """GeoPackage geometry blob: 'GP' header (+ envelope) + WKB."""
    if geom is None:
        return None
    if geom_type == "POINT":
        body = wkb_point(*geom)
        env = None
    elif geom_type == "LINESTRING":
        body = wkb_linestring(geom)
        env = _envelope(geom)
    elif geom_type == "MULTIPOLYGON":
        body = wkb_multipolygon(geom)
        env = _envelope(pt for poly in geom for pt in poly[0])
    else:
        raise GpkgError(f"unsupported geometry type {geom_type}")
    if env is None:
        flags = 0b00000001                       # little endian, no envelope
        header = b"GP" + struct.pack("<BBi", 0, flags, int(srs_id))
    else:
        flags = 0b00000011                       # little endian, xy envelope
        header = (b"GP" + struct.pack("<BBi", 0, flags, int(srs_id))
                  + struct.pack("<dddd", *env))
    return header + body


def parse_gpkg_blob(blob):
    """Inverse of gpkg_blob for the supported types. Returns (type, geom, srs)."""
    if blob is None:
        return None, None, None
    if blob[:2] != b"GP":
        raise GpkgError("not a GeoPackage geometry blob")
    flags = blob[3]
    srs = struct.unpack("<i", blob[4:8])[0]
    env_code = (flags >> 1) & 0b111
    env_len = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}[env_code]
    wkb = blob[8 + env_len:]
    geom, _ = _parse_wkb(wkb, 0)
    return geom[0], geom[1], srs


def _parse_wkb(b, o):
    order = b[o]
    fmt = "<" if order == 1 else ">"
    t = struct.unpack(fmt + "I", b[o + 1:o + 5])[0] % 1000
    o += 5

    def ring(o):
        n = struct.unpack(fmt + "I", b[o:o + 4])[0]; o += 4
        pts = [struct.unpack(fmt + "dd", b[o + 16 * k:o + 16 * k + 16]) for k in range(n)]
        return pts, o + 16 * n

    if t == 1:
        return ("POINT", struct.unpack(fmt + "dd", b[o:o + 16])), o + 16
    if t == 2:
        pts, o = ring(o)
        return ("LINESTRING", pts), o
    if t == 3:
        n = struct.unpack(fmt + "I", b[o:o + 4])[0]; o += 4
        rings = []
        for _ in range(n):
            r, o = ring(o); rings.append(r)
        return ("POLYGON", rings), o
    if t == 6:
        n = struct.unpack(fmt + "I", b[o:o + 4])[0]; o += 4
        polys = []
        for _ in range(n):
            (_, p), o = _parse_wkb(b, o); polys.append(p)
        return ("MULTIPOLYGON", polys), o
    raise GpkgError(f"unsupported WKB type {t}")


# -- writer ----------------------------------------------------------------

def _clean(value, ftype):
    """Python value -> SQLite value; NaN/inf/None -> NULL (never 0.0)."""
    if value is None:
        return None
    if ftype == "real":
        v = float(value)
        return v if math.isfinite(v) else None
    if ftype == "int":
        try:
            v = float(value)
        except (TypeError, ValueError):
            return None
        return int(v) if math.isfinite(v) else None
    return str(value)


class GeoPackageWriter(object):
    """Create a new GeoPackage file and add tables to it.

        w = GeoPackageWriter(path, timestamp="2026-09-29T00:00:00.000Z")
        w.add_srs(32737, "EPSG", 32737, wkt, "WGS 84 / UTM zone 37S")
        w.add_feature_table("crossings", "POINT", 32737, fields, rows)
        w.add_attribute_table("qeht_run_metadata", fields, rows)
        w.close()

    `fields` is a list of (name, type) with type in text/int/real.
    Feature rows are (geometry, {field: value}); attribute rows are dicts.
    `timestamp` fixes gpkg_contents.last_change (for reproducible files).
    """

    def __init__(self, path, overwrite=True, timestamp=None):
        if os.path.exists(path):
            if not overwrite:
                raise GpkgError(f"{path} exists")
            os.remove(path)
        for ext in ("-journal", "-wal", "-shm"):
            if os.path.exists(path + ext):
                os.remove(path + ext)
        self.path = path
        self.timestamp = timestamp
        self.con = sqlite3.connect(path)
        c = self.con
        c.execute(f"PRAGMA application_id = {GPKG_APPLICATION_ID}")
        c.execute(f"PRAGMA user_version = {GPKG_USER_VERSION}")
        c.executescript("""
CREATE TABLE gpkg_spatial_ref_sys (
  srs_name TEXT NOT NULL, srs_id INTEGER NOT NULL PRIMARY KEY,
  organization TEXT NOT NULL, organization_coordsys_id INTEGER NOT NULL,
  definition TEXT NOT NULL, description TEXT);
CREATE TABLE gpkg_contents (
  table_name TEXT NOT NULL PRIMARY KEY, data_type TEXT NOT NULL,
  identifier TEXT UNIQUE, description TEXT DEFAULT '',
  last_change DATETIME NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%fZ','now')),
  min_x DOUBLE, min_y DOUBLE, max_x DOUBLE, max_y DOUBLE, srs_id INTEGER,
  CONSTRAINT fk_gc_r_srs_id FOREIGN KEY (srs_id) REFERENCES gpkg_spatial_ref_sys(srs_id));
CREATE TABLE gpkg_geometry_columns (
  table_name TEXT NOT NULL, column_name TEXT NOT NULL,
  geometry_type_name TEXT NOT NULL, srs_id INTEGER NOT NULL,
  z TINYINT NOT NULL, m TINYINT NOT NULL,
  CONSTRAINT pk_geom_cols PRIMARY KEY (table_name, column_name),
  CONSTRAINT uk_gc_table_name UNIQUE (table_name),
  CONSTRAINT fk_gc_tn FOREIGN KEY (table_name) REFERENCES gpkg_contents(table_name),
  CONSTRAINT fk_gc_srs FOREIGN KEY (srs_id) REFERENCES gpkg_spatial_ref_sys (srs_id));
""")
        c.executemany(
            "INSERT INTO gpkg_spatial_ref_sys VALUES (?,?,?,?,?,?)",
            [("Undefined cartesian SRS", -1, "NONE", -1, "undefined",
              "undefined cartesian coordinate reference system"),
             ("Undefined geographic SRS", 0, "NONE", 0, "undefined",
              "undefined geographic coordinate reference system"),
             ("WGS 84 geodetic", 4326, "EPSG", 4326, WGS84_WKT,
              "longitude/latitude coordinates in decimal degrees on the WGS 84 spheroid")])

    def add_srs(self, srs_id, organization, org_id, definition, name):
        exists = self.con.execute("SELECT 1 FROM gpkg_spatial_ref_sys WHERE srs_id=?",
                                  (int(srs_id),)).fetchone()
        if exists:
            return
        self.con.execute("INSERT INTO gpkg_spatial_ref_sys VALUES (?,?,?,?,?,?)",
                         (name or f"{organization}:{org_id}", int(srs_id),
                          organization, int(org_id), definition or "undefined", None))

    def _contents(self, table, data_type, description, bbox=None, srs_id=None):
        _ident(table)
        cols = "table_name, data_type, identifier, description, min_x, min_y, max_x, max_y, srs_id"
        vals = [table, data_type, table, description or ""]
        vals += list(bbox) if bbox else [None, None, None, None]
        vals.append(srs_id)
        if self.timestamp:
            cols += ", last_change"
            vals.append(self.timestamp)
        marks = ",".join("?" * len(vals))
        # SQL identifiers come from _ident() (validated); values are bound.
        sql = f"INSERT INTO gpkg_contents ({cols}) VALUES ({marks})"  # nosec B608
        self.con.execute(sql, vals)

    @staticmethod
    def _column_sql(fields):
        for name, ftype in fields:
            if ftype not in _SQL_TYPES:
                raise GpkgError(f"field {name}: unknown type {ftype}")
        return ", ".join(f"{_ident(n)} {_SQL_TYPES[t]}" for n, t in fields)

    def add_feature_table(self, table, geom_type, srs_id, fields, rows,
                          description="", geom_column="geom"):
        if geom_type not in ("POINT", "LINESTRING", "MULTIPOLYGON"):
            raise GpkgError(f"unsupported geometry type {geom_type}")
        cols = self._column_sql(fields)
        self.con.execute(
            f"CREATE TABLE {_ident(table)} (fid INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, "
            f"{_ident(geom_column)} {geom_type}" + (f", {cols}" if cols else "") + ")")
        xs, ys = [], []
        names = [n for n, _ in fields]
        colnames = ", ".join([_ident(geom_column)] + [_ident(n) for n in names])
        marks = ",".join("?" * (len(names) + 1))
        # SQL identifiers come from _ident() (validated); values are bound.
        sql = f"INSERT INTO {_ident(table)} ({colnames}) VALUES ({marks})"  # nosec B608
        for geom, attrs in rows:
            if geom is not None:
                if geom_type == "POINT":
                    xs.append(geom[0]); ys.append(geom[1])
                elif geom_type == "LINESTRING":
                    xs += [p[0] for p in geom]; ys += [p[1] for p in geom]
                else:
                    for poly in geom:
                        xs += [p[0] for p in poly[0]]; ys += [p[1] for p in poly[0]]
            vals = [gpkg_blob(geom_type, geom, srs_id)]
            vals += [_clean(attrs.get(n), t) for n, t in fields]
            self.con.execute(sql, vals)
        bbox = (min(xs), min(ys), max(xs), max(ys)) if xs else None
        self._contents(table, "features", description, bbox, srs_id)
        self.con.execute("INSERT INTO gpkg_geometry_columns VALUES (?,?,?,?,0,0)",
                         (table, geom_column, geom_type, int(srs_id)))

    def add_attribute_table(self, table, fields, rows, description=""):
        cols = self._column_sql(fields)
        self.con.execute(f"CREATE TABLE {_ident(table)} (fid INTEGER PRIMARY KEY "
                         f"AUTOINCREMENT NOT NULL, {cols})")
        names = [n for n, _ in fields]
        colnames = ", ".join(_ident(n) for n in names)
        marks = ",".join("?" * len(names))
        # SQL identifiers come from _ident() (validated); values are bound.
        sql = f"INSERT INTO {_ident(table)} ({colnames}) VALUES ({marks})"  # nosec B608
        for attrs in rows:
            self.con.execute(sql, [_clean(attrs.get(n), t) for n, t in fields])
        self._contents(table, "attributes", description)

    def close(self):
        if self.con is not None:
            self.con.commit()
            self.con.execute("VACUUM")
            self.con.close()
            self.con = None


# -- reader ----------------------------------------------------------------

def read_table(path, table, with_geometry=True):
    """Read one table as a list of dicts (geometry decoded under '_geom')."""
    con = sqlite3.connect(path)
    try:
        con.row_factory = sqlite3.Row
        geom_col = con.execute(
            "SELECT column_name FROM gpkg_geometry_columns WHERE table_name=?",
            (table,)).fetchone()
        rows = []
        # SQL identifiers come from _ident() (validated); values are bound.
        query = f"SELECT * FROM {_ident(table)} ORDER BY fid"  # nosec B608
        for r in con.execute(query):
            d = dict(r)
            if geom_col is not None:
                blob = d.pop(geom_col[0])
                if with_geometry:
                    d["_geom"] = parse_gpkg_blob(blob)[1]
            rows.append(d)
        return rows
    finally:
        con.close()


def list_tables(path):
    con = sqlite3.connect(path)
    try:
        return [(r[0], r[1]) for r in con.execute(
            "SELECT table_name, data_type FROM gpkg_contents ORDER BY table_name")]
    finally:
        con.close()


def table_columns(path, table):
    con = sqlite3.connect(path)
    try:
        return [(r[1], r[2]) for r in con.execute("PRAGMA table_info(" + _ident(table) + ")")]
    finally:
        con.close()
