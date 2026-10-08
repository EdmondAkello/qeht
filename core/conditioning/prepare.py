# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Prepare a DEM for hydrology (F9, v0.21).

Most bad runs came from DEM preparation outside the plugin: nearest-
neighbour reprojection, an export finer than the native grid (repeated
rows and columns) and undeclared NoData zeros. This module does the
preparation in-process with the GDAL bindings:

Target CRS: "auto" picks the WGS 84 / UTM zone of the extent centre
anywhere on Earth (UPS beyond 84 N / 80 S; Norway and Svalbard
exceptions), or any projected CRS the user gives.

1. gdal.BuildVRT over the tiles (in memory), with the NoData override as
   the source NoData when given.
2. Optional clip: the clip bounds (any CRS) are transformed to the target
   CRS (edges densified) and widened by the buffer.
3. gdal.Warp to the target CRS, cell size and resampling (bilinear or
   cubic; nearest is refused), float32, NoData -9999, target-aligned
   pixels (-tap).
4. audit_nodata and audit_resampling on the INPUT (the merged tiles) and
   on the OUTPUT. Repeats in the input mean the damage is upstream:
   resampling cannot undo it, get the native product.
5. GeoTIFF tags QEHT_DEM_SOURCE and QEHT_DEM_PREP; the design hydrology
   package reads dem_source from them when the user leaves it blank.

No QGIS imports; GDAL / OSR imports are function-local.
"""

import datetime
import json
import math
import os

import numpy as np

from ..raster import audit_nodata, audit_resampling, resampling_stats

RESAMPLING = ("bilinear", "cubic")
OUT_NODATA = -9999.0
AUDIT_MAX_CELLS = 16_000_000          # larger inputs are audited on a central window


def utm_zone(lon, lat):
    """UTM zone number at (lon, lat), with the Norway (32V) and Svalbard (31X-37X)
    exceptions of the UTM grid."""
    lon = ((lon + 180.0) % 360.0) - 180.0
    zone = min(int(math.floor((lon + 180.0) / 6.0)) + 1, 60)
    if 56.0 <= lat < 64.0 and 3.0 <= lon < 12.0:
        zone = 32
    elif 72.0 <= lat <= 84.0 and lon >= 0.0:
        if lon < 9.0:
            zone = 31
        elif lon < 21.0:
            zone = 33
        elif lon < 33.0:
            zone = 35
        elif lon < 42.0:
            zone = 37
    return zone


def utm_epsg(lon, lat):
    """EPSG code of the WGS 84 projected CRS for (lon, lat), anywhere on Earth:
    the UTM zone (north 326zz, south 327zz), or UPS North (32661) above 84 N and
    UPS South (32761) below 80 S, where UTM is not defined."""
    if lat > 84.0:
        return 32661
    if lat < -80.0:
        return 32761
    return (32700 if lat < 0 else 32600) + utm_zone(lon, lat)


def _crs_name(epsg):
    if epsg == 32661:
        return "WGS 84 / UPS North"
    if epsg == 32761:
        return "WGS 84 / UPS South"
    return f"WGS 84 / UTM zone {epsg % 100}{'S' if epsg > 32700 else 'N'}"


def metres_per_degree(lat):
    """(m per degree of longitude, m per degree of latitude) on WGS 84 at `lat`."""
    p = math.radians(lat)
    m_lat = 111132.92 - 559.82 * math.cos(2 * p) + 1.175 * math.cos(4 * p) - 0.0023 * math.cos(6 * p)
    m_lon = 111412.84 * math.cos(p) - 93.5 * math.cos(3 * p) + 0.118 * math.cos(5 * p)
    return m_lon, m_lat


def _srs(wkt_or_epsg):
    from osgeo import osr
    s = osr.SpatialReference()
    if isinstance(wkt_or_epsg, int):
        s.ImportFromEPSG(wkt_or_epsg)
    else:
        s.SetFromUserInput(str(wkt_or_epsg))
    s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return s


def transform_bounds(bounds, src, dst, densify=21):
    """(xmin, ymin, xmax, ymax) from CRS src to dst (WKT, 'EPSG:n' or int), edges densified."""
    from osgeo import osr
    t = osr.CoordinateTransformation(_srs(src), _srs(dst))
    x0, y0, x1, y1 = bounds
    u = np.linspace(0.0, 1.0, densify)
    pts = ([(x0 + (x1 - x0) * a, y0) for a in u] + [(x0 + (x1 - x0) * a, y1) for a in u]
           + [(x0, y0 + (y1 - y0) * a) for a in u] + [(x1, y0 + (y1 - y0) * a) for a in u])
    out = np.array([t.TransformPoint(float(x), float(y))[:2] for x, y in pts])
    return float(out[:, 0].min()), float(out[:, 1].min()), float(out[:, 0].max()), float(out[:, 1].max())


def native_cell_m(geotransform, wkt, centre_lat=None):
    """Input cell size in metres (x, y): degrees converted at the centre latitude
    for a geographic CRS, else the geotransform as it is."""
    s = _srs(wkt)
    dx, dy = abs(geotransform[1]), abs(geotransform[5])
    if s.IsGeographic():
        m_lon, m_lat = metres_per_degree(centre_lat or 0.0)
        return dx * m_lon, dy * m_lat
    f = s.GetLinearUnits() or 1.0
    return dx * f, dy * f


def _read_for_audit(ds):
    """Array, valid mask and declared NoData of band 1 (central window when large)."""
    b = ds.GetRasterBand(1)
    nd = b.GetNoDataValue()
    w, h = ds.RasterXSize, ds.RasterYSize
    if w * h > AUDIT_MAX_CELLS:
        side = int(math.sqrt(AUDIT_MAX_CELLS))
        xo, yo = max((w - side) // 2, 0), max((h - side) // 2, 0)
        a = b.ReadAsArray(xo, yo, min(side, w), min(side, h)).astype(np.float64)
        window = [xo, yo, min(side, w), min(side, h)]
    else:
        a = b.ReadAsArray().astype(np.float64)
        window = None
    valid = np.isfinite(a)
    if nd is not None:
        valid &= a != nd
    return a, valid, nd, window


def audit(ds):
    a, v, nd, window = _read_for_audit(ds)
    st = resampling_stats(a, v)
    return {"nodata_declared": nd, "window": window, "valid_cells": int(v.sum()),
            "warnings": audit_nodata(a, v, nd) + audit_resampling(a, v),
            "duplicated_row_share": st["duplicated_row_share"],
            "duplicated_col_share": st["duplicated_col_share"], "tied_share": st["tied_share"]}


def prepare_dem(paths, out_path, target="auto", cell_m=None, resampling="bilinear",
                clip_bounds=None, clip_crs=None, buffer_m=2000.0, nodata_override=None,
                source_text="", log_path=None, progress=None):
    """Merge, clip, reproject and audit DEM tiles -> the prepared DEM.

    target: "auto" (WGS 84 / UTM zone of the extent centre), an EPSG int or
    a WKT / 'EPSG:n' string. clip_bounds: (xmin, ymin, xmax, ymax) in
    clip_crs (WKT or 'EPSG:n'), or None for the whole merged extent.
    Returns the log dict (also written to log_path as JSON when given).
    """
    from osgeo import gdal
    gdal.UseExceptions()
    if resampling not in RESAMPLING:
        raise ValueError(f"Resampling must be one of {RESAMPLING}: nearest neighbour repeats "
                         "source rows and columns and creates artificial flats.")
    paths = [str(p) for p in paths]
    if not paths:
        raise ValueError("Give at least one DEM raster.")
    srcs = [gdal.Open(p) for p in paths]
    wkts = {ds.GetProjection() for ds in srcs}
    if any(not w for w in wkts):
        raise ValueError("A DEM tile has no CRS.")
    if len({_srs(w).ExportToProj4() for w in wkts}) > 1:
        raise ValueError("The DEM tiles are in different CRSs; give tiles of one product.")
    srcs = None
    vrt_path = f"/vsimem/qeht_prepare_{os.getpid()}_{id(paths)}.vrt"
    vopt = {}
    if nodata_override is not None:
        vopt.update(srcNodata=float(nodata_override), VRTNodata=float(nodata_override))
    vrt = gdal.BuildVRT(vrt_path, paths, options=gdal.BuildVRTOptions(**vopt))
    try:
        src_wkt = vrt.GetProjection()
        gt = vrt.GetGeoTransform()
        w, h = vrt.RasterXSize, vrt.RasterYSize
        src_bounds = (gt[0], gt[3] + h * gt[5], gt[0] + w * gt[1], gt[3])
        # centre for the auto UTM zone: of the clip when given, else of the tiles
        if clip_bounds is not None:
            cb = transform_bounds(clip_bounds, clip_crs or src_wkt, "EPSG:4326")
        else:
            cb = transform_bounds(src_bounds, src_wkt, "EPSG:4326")
        lon, lat = 0.5 * (cb[0] + cb[2]), 0.5 * (cb[1] + cb[3])
        zone_note = ""
        if target == "auto" or target is None:
            epsg = utm_epsg(lon, lat)
            dst = f"EPSG:{epsg}"
            target_desc = f"auto: {_crs_name(epsg)} (EPSG:{epsg}), centre {lon:.4f}, {lat:.4f}"
            corners = {utm_epsg(x, y) for x in (cb[0], cb[2]) for y in (cb[1], cb[3])}
            if len(corners) > 1:
                zone_note = (f"The extent spans {len(corners)} UTM zones; the zone of its centre "
                             "is used (distortion grows away from the central meridian - choose "
                             "a national or custom CRS for very wide areas).")
        else:
            dst = f"EPSG:{target}" if isinstance(target, int) else str(target)
            target_desc = dst if dst.upper().startswith("EPSG:") else _srs(dst).GetName()
        if _srs(dst).IsGeographic():
            raise ValueError("The target CRS is geographic; QEHT needs a projected, metric CRS.")
        ncx, ncy = native_cell_m(gt, src_wkt, lat)
        cell = float(cell_m) if cell_m else round(max(ncx, ncy), 1)
        if cell <= 0:
            raise ValueError("The cell size must be positive.")
        wopt = dict(format="GTiff", dstSRS=dst, xRes=cell, yRes=cell, targetAlignedPixels=True,
                    resampleAlg=resampling, outputType=gdal.GDT_Float32, dstNodata=OUT_NODATA,
                    creationOptions=["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"],
                    multithread=False)
        if nodata_override is not None:
            wopt["srcNodata"] = float(nodata_override)
        clip_out = None
        if clip_bounds is not None:
            tb = transform_bounds(clip_bounds, clip_crs or src_wkt, dst)
            b = float(buffer_m or 0.0)
            clip_out = (tb[0] - b, tb[1] - b, tb[2] + b, tb[3] + b)
            wopt["outputBounds"] = clip_out
        if progress is not None:
            wopt["callback"] = lambda frac, msg, data: progress(frac) or 1
        if os.path.exists(out_path):
            gdal.GetDriverByName("GTiff").Delete(out_path)
        out = gdal.Warp(out_path, vrt, options=gdal.WarpOptions(**wopt))
        a_in = audit(vrt)
        a_out = audit(out)
        when = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        prep = (f"QEHT prepare DEM {when}: {len(paths)} tile(s) -> {target_desc}, cell {cell:g} m, "
                f"{resampling}, float32, NoData {OUT_NODATA:g}, target-aligned"
                + (f", clipped with a {buffer_m:g} m buffer" if clip_out else "")
                + (f", input NoData override {nodata_override:g}" if nodata_override is not None else ""))
        md = {"QEHT_DEM_PREP": prep, "QEHT_TOOL": "Prepare DEM for hydrology"}
        if source_text:
            md["QEHT_DEM_SOURCE"] = source_text
        out.SetMetadata(md)
        ogt = out.GetGeoTransform()
        out_info = {"crs": target_desc, "cols": out.RasterXSize, "rows": out.RasterYSize,
                    "cell_m": [abs(ogt[1]), abs(ogt[5])],
                    "bounds": [ogt[0], ogt[3] + out.RasterYSize * ogt[5],
                               ogt[0] + out.RasterXSize * ogt[1], ogt[3]]}
        out.FlushCache()
        out = None
    finally:
        vrt = None
        gdal.Unlink(vrt_path)
    upstream = (a_in["duplicated_row_share"] >= 0.002 or a_in["duplicated_col_share"] >= 0.002)
    log = {
        "tool": "Prepare DEM for hydrology", "qeht_prep": prep, "inputs": paths,
        "dem_source": source_text or "", "input_crs": _srs(src_wkt).GetName(),
        "input_cell": [abs(gt[1]), abs(gt[5])], "input_cell_m": [ncx, ncy],
        "native_estimate_m": [ncx / max(1e-9, 1.0 - a_in["duplicated_col_share"]),
                              ncy / max(1e-9, 1.0 - a_in["duplicated_row_share"])],
        "target_crs": target_desc, "cell_m": cell, "resampling": resampling,
        "nodata_override": nodata_override, "clip_bounds_target": clip_out,
        "buffer_m": buffer_m if clip_out else None, "output": out_path, "output_info": out_info,
        "zone_note": zone_note, "audit_input": a_in, "audit_output": a_out, "input_damaged_upstream": bool(upstream),
        "advice": ("The INPUT already repeats rows or columns: the damage is upstream and "
                   "resampling cannot undo it. Get the native product (e.g. FABDEM / ALOS tiles "
                   "at their own grid) and prepare it here." if upstream else ""),
    }
    if log_path:
        with open(log_path, "w", encoding="utf-8") as f:
            json.dump(log, f, indent=1, default=str)
    return log
