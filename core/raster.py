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
    overlapping valid cells and are not constant), their shares, and
    tied_share (share of E and S neighbour pairs with identical values).
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
    return {"duplicated_rows": dr, "duplicated_cols": dc,
            "duplicated_row_share": dr / max(nr, 1), "duplicated_col_share": dc / max(nc, 1),
            "tied_share": ties / max(pairs, 1)}


def audit_resampling(array, valid, threshold=0.002):
    """Warn when the DEM looks nearest-neighbour resampled (see resampling_stats).

    threshold : share of duplicated rows or columns that triggers the
        warning (0.2 %). Calibrated on real clips: nearest-neighbour
        reprojections showed 0.45-1.5 %, native and bilinear grids 0 %.
    """
    st = resampling_stats(array, valid)
    worst = max(st["duplicated_row_share"], st["duplicated_col_share"])
    if worst < threshold:
        return []
    return [
        f"{st['duplicated_rows']:,} rows ({100 * st['duplicated_row_share']:.2f}%) and "
        f"{st['duplicated_cols']:,} columns ({100 * st['duplicated_col_share']:.2f}%) exactly "
        f"repeat their neighbour, and {100 * st['tied_share']:.1f}% of neighbouring cells are "
        "exactly equal. The DEM was probably resampled or reprojected with nearest neighbour. "
        "The repeats create artificial flats and ties, so flow routing on them is arbitrary "
        "and a catchment can switch between crossings. Resample the source DEM with bilinear "
        "(or cubic) into the local UTM zone and run again."]


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
                 valid=None, compress=True):
    """Write a NumPy array as a GeoTIFF using the georeferencing in `info`."""
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

    drv = ogr.GetDriverByName("GPKG")
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

    drv = ogr.GetDriverByName("GPKG")
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
