# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Raster and vector I/O through the GDAL Python bindings only.

Every GDAL import is function-local and lazy. That is not stylistic: it
keeps the numeric modules in `qeht.core` importable and testable on a
machine with nothing but NumPy, and it guarantees that no import-time
side effect can drag in a code path we have not audited.

NOTHING in this module spawns a subprocess. `gdal_polygonize.exe`,
`gdal_translate.exe` and friends are never invoked. `gdal.Polygonize()`
is a direct call into the already-loaded GDAL C++ library.
"""

import os

import numpy as np


_OGR_BY_EXT = {".shp": "ESRI Shapefile", ".geojson": "GeoJSON", ".json": "GeoJSON",
               ".gml": "GML", ".kml": "KML", ".sqlite": "SQLite", ".fgb": "FlatGeobuf"}


def raster_tags(path):
    """GeoTIFF metadata of a raster as a dict ({} if unreadable)."""
    try:
        from osgeo import gdal
        ds = gdal.Open(str(path))
        md = dict(ds.GetMetadata() or {}) if ds is not None else {}
        ds = None
        return md
    except Exception:
        return {}


def ogr_driver_for(path):
    """OGR driver matching the output file's extension (GeoPackage otherwise).

    Processing passes the file name the user chose; writing a GeoPackage into a
    file called .shp (QEHT <= 0.15) produced files other software could not open.
    """
    from osgeo import ogr
    import os
    ext = os.path.splitext(str(path))[1].lower()
    return ogr.GetDriverByName(_OGR_BY_EXT.get(ext, "GPKG"))


class RasterError(Exception):
    """Raised when a raster cannot be used for hydrological analysis."""


class RasterInfo(object):
    """Georeferencing carried alongside the array, so core functions never
    need to know about GDAL datasets."""

    def __init__(self, geotransform, projection_wkt, rows, cols, nodata=None):
        self.geotransform = tuple(geotransform)
        self.projection_wkt = projection_wkt
        self.rows = int(rows)
        self.cols = int(cols)
        self.nodata = nodata

    @property
    def cell_width(self):
        return abs(self.geotransform[1])

    @property
    def cell_height(self):
        return abs(self.geotransform[5])

    @property
    def cell_area(self):
        return self.cell_width * self.cell_height

    def rowcol_to_xy(self, row, col, centre=True):
        gt = self.geotransform
        off = 0.5 if centre else 0.0
        x = gt[0] + (col + off) * gt[1] + (row + off) * gt[2]
        y = gt[3] + (col + off) * gt[4] + (row + off) * gt[5]
        return x, y

    def xy_to_rowcol(self, x, y):
        gt = self.geotransform
        col = int((x - gt[0]) / gt[1])
        row = int((y - gt[3]) / gt[5])
        return row, col

    def __repr__(self):
        return (f"RasterInfo({self.rows}x{self.cols}, "
                f"cell {self.cell_width:g}x{self.cell_height:g})")


def audit_nodata(array, valid, declared_nodata):
    """Detect undeclared NoData - the single most damaging DEM defect.

    Found on a real project DEM: 11.8% of cells held 0.0 with no NoData
    value set in the header, while genuine terrain started at 1574 m.
    Treated as valid, those cells form a sea at 0 m that the entire
    catchment drains into, and every delineated boundary is wrong.

    Returns a list of human-readable warnings. Nothing is auto-corrected;
    the user decides, because a legitimate DEM can contain real zeros.
    """
    warnings = []
    live = array[valid]
    if live.size == 0:
        return ["Raster contains no valid cells."]

    for suspect in (0.0, -9999.0, -32768.0, -3.4028235e38):
        n = int(np.count_nonzero(array[valid] == suspect))
        if n == 0:
            continue
        rest = live[live != suspect]
        if rest.size == 0:
            continue
        share = 100.0 * n / live.size
        gap = rest.min() - suspect
        spread = float(rest.max() - rest.min())
        if share > 0.5 and spread > 0 and gap > spread * 0.25:
            warnings.append(
                f"{n:,} cells ({share:.2f}%) hold the value {suspect:g}, but "
                f"real elevations start at {rest.min():.1f}. This is almost "
                f"certainly undeclared NoData. Set the NoData override, or "
                f"every catchment boundary will be wrong.")
    return warnings


def resampling_stats(array, valid, min_overlap=50):
    """Signature of nearest-neighbour resampling (v0.13.1).

    A DEM reprojected or regridded with nearest neighbour repeats whole
    source rows and columns wherever the output grid is finer than the
    source along that axis. The repeats are exact, so they create
    artificial flats and exactly tied neighbours: flat routing then decides
    where flow goes, and a single tied cell on a river can move a large
    catchment from one crossing to another (seen on a real project:
    ~176 km2 moved between two culverts).

    Returns a dict: duplicated_rows / duplicated_cols (counts of rows or
    columns that exactly repeat the previous one over >= 99 % of their
    overlapping valid cells and are not constant), their shares,
    tied_share (share of E and S neighbour pairs with identical values) and
    sloped_tied_share (1.0.2): the share of pairs that are tied although
    both cells differ from some other neighbour by more than
    SLOPED_TIE_RELIEF_M. Reprojection with nearest neighbour onto a rotated
    or shifted grid repeats source cells without repeating whole rows, so
    the row test misses it; the copies still show up as ties on sloping
    ground, where a continuous surface has almost none. Ties inside flat
    patches (water bodies, flattened areas, quantised plains) are excluded.
    """
    a = np.asarray(array)
    v = np.asarray(valid, dtype=bool)

    def _dups(x, m):
        ov = m[1:] & m[:-1]
        n_ov = ov.sum(axis=1)
        eq = ((x[1:] == x[:-1]) & ov).sum(axis=1)
        xmask = np.where(ov, x[1:], np.nan)
        with np.errstate(all="ignore"):
            import warnings as _w
            with _w.catch_warnings():
                _w.simplefilter("ignore", RuntimeWarning)
                span = np.nanmax(xmask, axis=1) - np.nanmin(xmask, axis=1)
        dup = (n_ov >= min_overlap) & (eq >= 0.99 * n_ov) & (np.nan_to_num(span) > 0)
        considered = int((n_ov >= min_overlap).sum())
        return int(dup.sum()), considered

    dr, nr = _dups(a, v)
    dc, nc = _dups(a.T, v.T)
    pe = v[:, 1:] & v[:, :-1]
    ps = v[1:] & v[:-1]
    ties = int(((a[:, 1:] == a[:, :-1]) & pe).sum() + ((a[1:] == a[:-1]) & ps).sum())
    pairs = int(pe.sum() + ps.sum())
    aw, vw = _central_window(a, v, AUDIT_WINDOW_CELLS)
    return {"duplicated_rows": dr, "duplicated_cols": dc,
            "duplicated_row_share": dr / max(nr, 1), "duplicated_col_share": dc / max(nc, 1),
            "tied_share": ties / max(pairs, 1),
            "sloped_tied_share": _sloped_tie_share(aw, vw, SLOPED_TIE_RELIEF_M)}


SLOPED_TIE_RELIEF_M = 0.1        # a tie counts as "on sloping ground" above this local relief
SLOPED_TIE_THRESHOLD = 0.02      # share of tied pairs on sloping ground that triggers a warning
AUDIT_WINDOW_CELLS = 4_000_000   # neighbourhood indicators use a central window of large DEMs
QUANTUM_CANDIDATES = (1.0, 0.5, 0.25, 0.2, 0.1, 0.05, 0.02, 0.01)


def _central_window(array, valid, max_cells):
    """Central window of at most max_cells cells (the whole grid when smaller)."""
    rows, cols = array.shape
    if rows * cols <= max_cells:
        return array, valid
    side = int(max_cells ** 0.5)
    r0, c0 = max((rows - side) // 2, 0), max((cols - side) // 2, 0)
    return array[r0:r0 + side, c0:c0 + side], valid[r0:r0 + side, c0:c0 + side]


def _neighbour_stack(z):
    """Yield the 8 neighbour grids of z (NaN outside the grid), D8 order not needed."""
    rows, cols = z.shape
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue
            nb = np.full(z.shape, np.nan)
            nb[max(-dr, 0):rows - max(dr, 0), max(-dc, 0):cols - max(dc, 0)] = \
                z[max(dr, 0):rows - max(-dr, 0), max(dc, 0):cols - max(-dc, 0)]
            yield nb


def _sloped_tie_share(array, valid, relief_m):
    """Share of E/S neighbour pairs that are exactly tied while both cells have local relief > relief_m."""
    z = np.where(valid, np.asarray(array, dtype=np.float64), np.nan)
    relief = np.zeros(z.shape)
    for nb in _neighbour_stack(z):
        with np.errstate(invalid="ignore"):
            relief = np.fmax(relief, np.abs(z - nb))
    relief = np.nan_to_num(relief)
    pe = valid[:, 1:] & valid[:, :-1]
    ps = valid[1:] & valid[:-1]
    te = (z[:, 1:] == z[:, :-1]) & pe & (np.minimum(relief[:, 1:], relief[:, :-1]) > relief_m)
    ts = (z[1:] == z[:-1]) & ps & (np.minimum(relief[1:], relief[:-1]) > relief_m)
    return float((te.sum() + ts.sum()) / max(int(pe.sum() + ps.sum()), 1))


def vertical_quantum(array, valid, candidates=QUANTUM_CANDIDATES, need=0.999, sample=200_000):
    """Storage step of the elevations, or None when they look continuous (1.0.2).

    Returns the largest candidate step q (metres) such that at least `need`
    of the valid values lie on a multiple of q, allowing for float32
    rounding. Whole-metre products (ALOS AW3D30, SRTM) give 1.0 even when
    stored as float; FABDEM gives 0.01; a bilinear-resampled DEM gives None.
    """
    x = np.asarray(array, dtype=np.float64)[np.asarray(valid, dtype=bool)]
    if x.size == 0:
        return None
    if x.size > sample:
        x = x[np.random.default_rng(0).choice(x.size, sample, replace=False)]
    eps = 4.0 * float(np.spacing(np.float32(np.max(np.abs(x)))))
    for q in candidates:
        off = np.abs(x - np.round(x / q) * q)
        if (off <= max(0.02 * q, eps)).mean() >= need:
            return float(q)
    return None


def mean_steepest_drop(array, valid):
    """Mean over valid cells of the steepest drop (elevation units) to any of the 8 neighbours.

    Cells with no lower neighbour count as zero. This is the "average drop
    per pixel" of Gyasi-Agyei et al. (1995), computed on the DEM as given.
    """
    v = np.asarray(valid, dtype=bool)
    z = np.where(v, np.asarray(array, dtype=np.float64), np.nan)
    best = np.zeros(z.shape)
    for nb in _neighbour_stack(z):
        with np.errstate(invalid="ignore"):
            best = np.fmax(best, z - nb)
    vals = best[v]
    return float(np.nanmean(vals)) if vals.size else 0.0


def precision_stats(array, valid):
    """Vertical quantum, mean steepest drop per pixel and their ratio (None when continuous)."""
    a, v = _central_window(np.asarray(array), np.asarray(valid, dtype=bool), AUDIT_WINDOW_CELLS)
    q = vertical_quantum(a, v)
    drop = mean_steepest_drop(a, v)
    return {"quantum": q, "mean_drop": drop, "ratio": (drop / q) if q else None}


def audit_precision(array, valid, caution_ratio=10.0):
    """Warn when the elevation storage step is coarse for the relief (1.0.2).

    Gyasi-Agyei et al. (1995, Hydrol. Process. 9, 363-382) proposed that a
    DEM's vertical resolution is adequate for extracting the drainage
    network when the average drop per pixel / vertical resolution exceeds
    one. Walker and Willgoose (1999, WRR 35, 2259-2268) found that
    condition necessary but not sufficient, and note that the threshold was
    fixed only to within a factor of ten. So: below 1 is a warning; below
    caution_ratio (10) is a caution that crossings on low-relief ground may
    depend on how flats are routed. Continuous DEMs return no warning.
    """
    st = precision_stats(array, valid)
    q, r = st["quantum"], st["ratio"]
    if q is None or r is None or r >= caution_ratio:
        return []
    step = f"{q:g} m" if q >= 0.1 else f"{100 * q:g} cm"
    head = (f"Elevations are stored in steps of {step}, and the mean steepest drop per cell is "
            f"{st['mean_drop']:.2f} m (ratio {r:.1f}). ")
    if r < 1.0:
        return [head + "That is below the adequacy ratio of 1 proposed by Gyasi-Agyei et al. (1995): "
                "much of the terrain is flat at this precision, so stream positions and catchment "
                "boundaries there are set by flat routing rather than by the terrain. Use a "
                "floating-point DEM if one is available, never round a DEM to whole metres, and "
                "check crossings on low-relief ground with the flat-method sensitivity output, "
                "mapped drainage or a site visit."]
    return [head + f"That meets the ratio of 1 proposed by Gyasi-Agyei et al. (1995) but is under "
            f"{caution_ratio:g}; Walker and Willgoose (1999) found that ratio necessary rather than "
            "sufficient. Crossings on low-relief parts of the DEM may still depend on flat routing: "
            "review the flat-method sensitivity output there, and prefer a floating-point DEM."]

def audit_resampling(array, valid, threshold=0.002, tie_threshold=SLOPED_TIE_THRESHOLD):
    """Warn when the DEM looks nearest-neighbour resampled (see resampling_stats).

    threshold : share of duplicated rows or columns that triggers the
        warning (0.2 %). Calibrated on real clips: nearest-neighbour
        reprojections showed 0.45-1.5 %, native and bilinear grids 0 %.
    tie_threshold : share of tied neighbour pairs on sloping ground that
        triggers the warning when no row repeats (1.0.2). Checked only when
        the elevations are continuous or stored in steps of 1 cm or finer;
        coarser storage steps produce ties of their own (see
        audit_precision). Calibrated on FABDEM clips reprojected to UTM:
        nearest neighbour gives a few per cent, bilinear and cubic almost
        none.
    """
    st = resampling_stats(array, valid)
    worst = max(st["duplicated_row_share"], st["duplicated_col_share"])
    advice = ("The repeats create artificial flats and ties, so flow routing on them is arbitrary "
              "and a catchment can switch between crossings. Resample the source DEM with bilinear "
              "(or cubic) into the local UTM zone and run again.")
    if worst >= threshold:
        return [
            f"{st['duplicated_rows']:,} rows ({100 * st['duplicated_row_share']:.2f}%) and "
            f"{st['duplicated_cols']:,} columns ({100 * st['duplicated_col_share']:.2f}%) exactly "
            f"repeat their neighbour, and {100 * st['tied_share']:.1f}% of neighbouring cells are "
            "exactly equal. The DEM was probably resampled or reprojected with nearest neighbour. "
            + advice]
    a, v = _central_window(np.asarray(array), np.asarray(valid, dtype=bool), AUDIT_WINDOW_CELLS)
    q = vertical_quantum(a, v)
    if (q is None or q <= 0.01) and st["sloped_tied_share"] >= tie_threshold:
        return [
            f"{100 * st['sloped_tied_share']:.1f}% of neighbouring cells on sloping ground have "
            "exactly the same elevation, although no whole row or column repeats. The DEM was "
            "probably reprojected or resampled with nearest neighbour onto a rotated or shifted "
            "grid (or the source already contains copied cells). " + advice]
    return []


def read_dem(path, band=1, nodata_override=None):
    """Read a DEM into (array, valid_mask, RasterInfo).

    Rejects rotated/skewed rasters outright. A north-up grid is assumed
    by every D8 algorithm; silently accepting a rotated one produces
    plausible-looking output with wrong flow directions, which is the
    worst possible failure mode for a drainage design.
    """
    from osgeo import gdal
    gdal.UseExceptions()

    ds = gdal.Open(path)
    if ds is None:
        raise RasterError(f"Could not open raster: {path}")

    gt = ds.GetGeoTransform()
    if abs(gt[2]) > 1e-9 or abs(gt[4]) > 1e-9:
        raise RasterError(
            "Raster is rotated or skewed. D8 analysis requires a north-up "
            "grid. Reproject the DEM to a projected CRS first.")
    if gt[1] <= 0 or gt[5] >= 0:
        raise RasterError(
            "Unexpected geotransform orientation. Expected west-to-east, "
            "north-to-south pixel ordering.")

    rb = ds.GetRasterBand(band)
    array = rb.ReadAsArray().astype(np.float64)
    nodata = rb.GetNoDataValue()

    valid = np.isfinite(array)
    if nodata is not None:
        valid &= (array != nodata)
    if nodata_override is not None:
        valid &= (array != float(nodata_override))

    info = RasterInfo(gt, ds.GetProjection(), ds.RasterYSize, ds.RasterXSize,
                      nodata)
    ds = None

    if not valid.any():
        raise RasterError("Raster contains no valid elevation cells.")
    return array, valid, info


def write_raster(path, array, info, dtype=None, nodata=-9999.0,
                 valid=None, compress=True, metadata=None):
    """Write a NumPy array as a GeoTIFF using the georeferencing in `info`.

    `metadata`: optional {key: value} written as GeoTIFF metadata (QEHT_*
    tags record how a raster was made, e.g. the flat method; the design
    hydrology package reads them back when the user leaves the fields blank).
    """
    from osgeo import gdal
    gdal.UseExceptions()

    array = np.asarray(array)
    if dtype is None:
        if np.issubdtype(array.dtype, np.integer):
            gdal_type, np_type = gdal.GDT_Int32, np.int32
        else:
            gdal_type, np_type = gdal.GDT_Float32, np.float32
    elif dtype == "int32":
        gdal_type, np_type = gdal.GDT_Int32, np.int32
    else:
        gdal_type, np_type = gdal.GDT_Float32, np.float32

    out = array.astype(np_type, copy=True)
    if valid is not None:
        out[~valid] = np_type(nodata)
    out = np.where(np.isfinite(out.astype(np.float64)), out, np_type(nodata))

    options = ["TILED=YES"]
    if compress:
        options += ["COMPRESS=DEFLATE", "PREDICTOR=2"]

    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(path, info.cols, info.rows, 1, gdal_type, options=options)
    ds.SetGeoTransform(info.geotransform)
    if info.projection_wkt:
        ds.SetProjection(info.projection_wkt)
    if metadata:
        ds.SetMetadata({str(k): str(v) for k, v in metadata.items()})
    band = ds.GetRasterBand(1)
    band.WriteArray(out)
    band.SetNoDataValue(float(nodata))
    band.FlushCache()
    ds = None
    return path


def polygonize(label_array, info, output_path, layer_name="catchments",
               field_name="DN", ignore_value=0, dissolve=True):
    """Raster -> polygon, entirely in-process via gdal.Polygonize().

    This is the call that the sample scripts were trying to reach by
    shelling out to gdal_polygonize.exe with a hidden console window.
    There is no need for any of that: the function below is a direct
    binding call into the GDAL library already loaded in the process.
    """
    from osgeo import gdal, ogr, osr
    gdal.UseExceptions()

    labels = np.asarray(label_array).astype(np.int32)

    # In-memory source raster; a mask band excludes the ignore value so
    # background does not become a giant polygon.
    mem = gdal.GetDriverByName("MEM")
    src = mem.Create("", info.cols, info.rows, 1, gdal.GDT_Int32)
    src.SetGeoTransform(info.geotransform)
    if info.projection_wkt:
        src.SetProjection(info.projection_wkt)
    src.GetRasterBand(1).WriteArray(labels)

    mask_ds = mem.Create("", info.cols, info.rows, 1, gdal.GDT_Byte)
    mask_ds.SetGeoTransform(info.geotransform)
    mask_ds.GetRasterBand(1).WriteArray(
        (labels != ignore_value).astype(np.uint8) * 255)

    srs = None
    if info.projection_wkt:
        srs = osr.SpatialReference()
        srs.ImportFromWkt(info.projection_wkt)

    drv = ogr_driver_for(output_path)
    import os
    if os.path.exists(output_path):
        drv.DeleteDataSource(output_path)
    vds = drv.CreateDataSource(output_path)
    layer = vds.CreateLayer(layer_name, srs=srs, geom_type=ogr.wkbPolygon)
    layer.CreateField(ogr.FieldDefn(field_name, ogr.OFTInteger))

    gdal.Polygonize(src.GetRasterBand(1), mask_ds.GetRasterBand(1),
                    layer, 0, [], callback=None)

    n_raw = layer.GetFeatureCount()

    if dissolve:
        _dissolve_by_field(layer, field_name, info)

    n_final = layer.GetFeatureCount()
    vds = None
    src = None
    mask_ds = None
    return {"path": output_path, "layer": layer_name,
            "features_raw": n_raw, "features": n_final}


def _dissolve_by_field(layer, field_name, info):
    """Merge polygons sharing a label into one multipart feature.

    gdal.Polygonize emits one polygon per contiguous region. A catchment
    interrupted by a NoData sliver arrives as several parts; an
    engineering deliverable wants one feature per catchment.
    """
    from osgeo import ogr

    groups = {}
    layer.ResetReading()
    for feat in layer:
        key = feat.GetField(field_name)
        geom = feat.GetGeometryRef()
        if geom is None:
            continue
        groups.setdefault(key, []).append(geom.Clone())

    if not groups or all(len(v) == 1 for v in groups.values()):
        layer.ResetReading()
        return

    ids = [f.GetFID() for f in layer]
    layer.ResetReading()
    for fid in ids:
        layer.DeleteFeature(fid)

    defn = layer.GetLayerDefn()
    for key, geoms in groups.items():
        merged = geoms[0]
        for g in geoms[1:]:
            merged = merged.Union(g)
        feat = ogr.Feature(defn)
        feat.SetField(field_name, key)
        feat.SetGeometry(merged)
        layer.CreateFeature(feat)
        feat = None
    layer.ResetReading()


def polyline_from_cells(cells, info, output_path, layer_name="flowpath",
                        attributes=None):
    """Write an ordered cell list as a single polyline (for longest flow path)."""
    from osgeo import ogr, osr
    import os

    srs = None
    if info.projection_wkt:
        srs = osr.SpatialReference()
        srs.ImportFromWkt(info.projection_wkt)

    drv = ogr_driver_for(output_path)
    if os.path.exists(output_path):
        drv.DeleteDataSource(output_path)
    vds = drv.CreateDataSource(output_path)
    layer = vds.CreateLayer(layer_name, srs=srs, geom_type=ogr.wkbLineString)

    attributes = attributes or {}
    for name, value in attributes.items():
        ftype = ogr.OFTReal if isinstance(value, float) else ogr.OFTInteger
        layer.CreateField(ogr.FieldDefn(name, ftype))

    line = ogr.Geometry(ogr.wkbLineString)
    for row, col in cells:
        x, y = info.rowcol_to_xy(row, col)
        line.AddPoint_2D(x, y)

    feat = ogr.Feature(layer.GetLayerDefn())
    for name, value in attributes.items():
        feat.SetField(name, value)
    feat.SetGeometry(line)
    layer.CreateFeature(feat)
    feat = None
    vds = None
    return output_path


def warp_to_grid(path, info, resampling="bilinear", band=1):
    """Read raster `path` resampled onto the DEM grid described by `info`.

    In-process GDAL Warp to memory (no subprocess): same CRS, extent, cell
    size and shape as the DEM. resampling: "bilinear" for continuous
    factors (R, K, C, P), "near" for class rasters (WorldCover).
    Returns (array float64 with NaN for NoData, source description).
    """
    from osgeo import gdal
    gdal.UseExceptions()
    gt = info.geotransform
    bounds = (gt[0], gt[3] + info.rows * gt[5], gt[0] + info.cols * gt[1], gt[3])
    src = gdal.Open(path)
    nd = src.GetRasterBand(band).GetNoDataValue()
    opts = gdal.WarpOptions(format="MEM", outputBounds=bounds, width=info.cols,
                            height=info.rows, dstSRS=info.projection_wkt or None,
                            resampleAlg=resampling, srcBands=[band],
                            dstNodata=-3.4e38, srcNodata=nd)
    ds = gdal.Warp("", src, options=opts)
    a = ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
    a[a <= -3.0e38] = np.nan
    desc = f"{os.path.basename(path)} (resampled {resampling} to the DEM grid)"
    ds = None
    src = None
    return a, desc
