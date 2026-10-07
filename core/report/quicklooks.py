# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Raster quicklooks (A6, v0.18): downsampled PNG + world file (.pgw) +
legend JSON for each raster, and a `rasters` index, so a designer working
outside a GIS can see the terrain and erosion grids under the crossings.

- Classified rasters (erosion classes) use their GDAL colour table and
  category names - the same colours as the .qml styles - and are
  downsampled by nearest neighbour, so class values are never blended.
- Continuous rasters use a stated ramp and stretch (2nd-98th percentile,
  or log for accumulation) and are block-averaged.
- The hillshade is computed from the raw DEM (Horn, azimuth 315°,
  altitude 45°) and blended over a hypsometric tint.

The PNG writer is pure Python (zlib); GDAL is imported only to read the
source rasters. No QGIS imports.
"""

import json
import math
import os
import struct
import zlib

import numpy as np

RAMPS = {   # name: list of (position 0-1, (r, g, b))
    "terrain": [(0.0, (46, 120, 92)), (0.3, (166, 196, 120)), (0.55, (236, 220, 160)),
                (0.8, (196, 150, 104)), (1.0, (150, 112, 92))],
    "blues": [(0.0, (240, 246, 251)), (0.5, (107, 174, 214)), (1.0, (8, 48, 107))],
    "viridis": [(0.0, (68, 1, 84)), (0.25, (59, 82, 139)), (0.5, (33, 145, 140)),
                (0.75, (94, 201, 98)), (1.0, (253, 231, 37))],
    "ylorrd": [(0.0, (255, 255, 204)), (0.5, (253, 141, 60)), (1.0, (189, 0, 38))],
}
INDEX_FIELDS = [("name", "text"), ("title", "text"), ("png", "text"), ("world_file", "text"),
                ("legend_json", "text"), ("crs_epsg", "text"), ("xmin", "real"), ("ymin", "real"),
                ("xmax", "real"), ("ymax", "real"), ("px_size_m", "real"),
                ("source_tool", "text"), ("source_raster", "text"), ("kind", "text")]


# -- PNG / world file ------------------------------------------------------------
def write_png(path, rgba):
    """8-bit RGBA PNG, pure Python. rgba: uint8 array (rows, cols, 4)."""
    a = np.ascontiguousarray(np.asarray(rgba, dtype=np.uint8))
    h, w = a.shape[:2]
    raw = b"".join(b"\x00" + a[r].tobytes() for r in range(h))

    def chunk(tag, data):
        return (struct.pack(">I", len(data)) + tag + data
                + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF))
    png = (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 6, 0, 0, 0))
           + chunk(b"IDAT", zlib.compress(raw, 6)) + chunk(b"IEND", b""))
    with open(path, "wb") as f:
        f.write(png)
    return path


def read_png_size(path):
    with open(path, "rb") as f:
        head = f.read(24)
    return struct.unpack(">II", head[16:24])           # (width, height)


def write_world_file(path, gt):
    """ESRI world file: pixel size x, rotation, rotation, pixel size y (negative),
    x and y of the CENTRE of the upper-left pixel."""
    lines = [gt[1], gt[4], gt[2], gt[5], gt[0] + gt[1] / 2.0 + gt[2] / 2.0,
             gt[3] + gt[4] / 2.0 + gt[5] / 2.0]
    with open(path, "w", encoding="ascii") as f:
        f.write("\n".join(f"{v:.10f}" for v in lines) + "\n")
    return path


# -- downsampling ------------------------------------------------------------------
def factor_for(shape, max_px):
    return max(1, int(math.ceil(max(shape) / float(max_px))))


def downsample(a, k, kind):
    """Integer factor k. classified: nearest (the cell nearest each block
    centre); continuous: mean of the finite cells in each block."""
    a = np.asarray(a)
    if k <= 1:
        return a.copy()
    rows, cols = a.shape
    nr, nc = int(math.ceil(rows / k)), int(math.ceil(cols / k))
    if kind == "classified":
        ri = np.minimum(np.arange(nr) * k + k // 2, rows - 1)
        ci = np.minimum(np.arange(nc) * k + k // 2, cols - 1)
        return a[np.ix_(ri, ci)]
    pad = np.full((nr * k, nc * k), np.nan)
    pad[:rows, :cols] = a
    blocks = pad.reshape(nr, k, nc, k)
    with np.errstate(invalid="ignore"):
        cnt = np.isfinite(blocks).sum(axis=(1, 3))
        s = np.nansum(blocks, axis=(1, 3))
    return np.where(cnt > 0, s / np.maximum(cnt, 1), np.nan)


# -- colour --------------------------------------------------------------------------
def ramp_rgb(t, ramp):
    stops = RAMPS[ramp] if isinstance(ramp, str) else ramp
    t = np.clip(np.asarray(t, float), 0.0, 1.0)
    pos = np.array([p for p, _ in stops])
    out = np.zeros(t.shape + (3,), float)
    for ch in range(3):
        out[..., ch] = np.interp(t, pos, [c[ch] for _, c in stops])
    return out


def render_continuous(v, ramp, lo, hi, log=False):
    v = np.asarray(v, float)
    ok = np.isfinite(v)
    if log:
        v = np.log10(np.maximum(v, 0.0) + 1.0)
        lo, hi = math.log10(max(lo, 0.0) + 1.0), math.log10(max(hi, 0.0) + 1.0)
    t = (v - lo) / (hi - lo) if hi > lo else np.zeros(v.shape)
    rgb = ramp_rgb(np.where(ok, t, 0.0), ramp)
    rgba = np.zeros(v.shape + (4,), np.uint8)
    rgba[..., :3] = np.round(rgb).astype(np.uint8)
    rgba[..., 3] = np.where(ok, 255, 0)
    return rgba


def render_classified(codes, colours, nodata=0):
    c = np.asarray(codes)
    rgba = np.zeros(c.shape + (4,), np.uint8)
    for code, rgb in colours.items():
        m = c == code
        rgba[m, 0], rgba[m, 1], rgba[m, 2] = rgb[:3]
        rgba[m, 3] = 255
    if nodata is not None:
        rgba[c == nodata, 3] = 0
    return rgba


def hillshade(z, cell_w, cell_h, azimuth=315.0, altitude=45.0):
    """Horn hillshade 0-1 (NaN where z is NaN)."""
    zp = np.pad(np.asarray(z, float), 1, mode="edge")
    zp = np.where(np.isnan(zp), np.nanmean(z) if np.isfinite(z).any() else 0.0, zp)
    a, b, c = zp[:-2, :-2], zp[:-2, 1:-1], zp[:-2, 2:]
    d, f = zp[1:-1, :-2], zp[1:-1, 2:]
    g, h, i = zp[2:, :-2], zp[2:, 1:-1], zp[2:, 2:]
    dzdx = ((c + 2 * f + i) - (a + 2 * d + g)) / (8.0 * cell_w)
    dzdy = ((g + 2 * h + i) - (a + 2 * b + c)) / (8.0 * cell_h)
    slope = np.arctan(np.hypot(dzdx, dzdy))
    aspect = np.arctan2(dzdy, -dzdx)
    az = math.radians(360.0 - azimuth + 90.0)
    alt = math.radians(altitude)
    hs = np.sin(alt) * np.cos(slope) + np.cos(alt) * np.sin(slope) * np.cos(az - aspect)
    return np.where(np.isfinite(z), np.clip(hs, 0.0, 1.0), np.nan)


def render_relief(z, cell_w, cell_h):
    """Hypsometric tint multiplied by the hillshade -> RGBA, plus the stretch."""
    ok = np.isfinite(z)
    lo, hi = (float(np.percentile(z[ok], 2)), float(np.percentile(z[ok], 98))) if ok.any() else (0, 1)
    tint = ramp_rgb(np.where(ok, (z - lo) / max(hi - lo, 1e-9), 0.0), "terrain")
    hs = hillshade(z, cell_w, cell_h)
    shade = 0.35 + 0.65 * np.nan_to_num(hs, nan=1.0)
    rgba = np.zeros(z.shape + (4,), np.uint8)
    rgba[..., :3] = np.clip(np.round(tint * shade[..., None]), 0, 255).astype(np.uint8)
    rgba[..., 3] = np.where(ok, 255, 0)
    return rgba, lo, hi


# -- driver ------------------------------------------------------------------------
def _read(path):
    from osgeo import gdal
    gdal.UseExceptions()
    ds = gdal.Open(path)
    b = ds.GetRasterBand(1)
    a = b.ReadAsArray().astype(np.float64)
    nd = b.GetNoDataValue()
    ct = b.GetColorTable()
    colours = {}
    if ct is not None:
        for k in range(ct.GetCount()):
            e = ct.GetColorEntry(k)
            if e[3] > 0:
                colours[k] = (int(e[0]), int(e[1]), int(e[2]))
    names = b.GetCategoryNames() or []
    gt = ds.GetGeoTransform()
    srs = ds.GetSpatialRef()
    epsg = ""
    if srs is not None:
        srs.AutoIdentifyEPSG()
        epsg = srs.GetAuthorityCode(None) or ""
    ds = None
    if nd is not None:
        a[a == nd] = np.nan
    return a, nd, gt, epsg, colours, names


def make_quicklook(item, out_dir, max_px=4096):
    """item: dict(name, title, path, kind = relief | continuous | classified,
    ramp, log, source_tool, units). Returns an index row (paths relative to
    out_dir's parent is the caller's choice; here they are file names)."""
    a, nd, gt, epsg, colours, names = _read(item["path"])
    if abs(gt[2]) > 0 or abs(gt[4]) > 0:
        raise ValueError(f"{item['path']}: rotated rasters are not supported")
    k = factor_for(a.shape, max_px)
    kind = item["kind"]
    legend = {"name": item["name"], "title": item["title"], "kind": kind,
              "source_raster": item["path"], "downsample_factor": k}
    if kind == "classified":
        codes = downsample(np.nan_to_num(a, nan=0).astype(np.int64), k, "classified")
        present = sorted(int(c) for c in np.unique(codes) if c != 0)
        cols = {c: colours.get(c, (128, 128, 128)) for c in present}
        rgba = render_classified(codes, cols, nodata=0)
        legend["classes"] = [{"value": c, "label": (names[c] if c < len(names) and names[c] else str(c)),
                              "colour": "#%02x%02x%02x" % cols[c]} for c in present]
        legend["resampling"] = "nearest (class values preserved)"
    elif kind == "relief":
        z = downsample(a, k, "continuous")
        rgba, lo, hi = render_relief(z, abs(gt[1]) * k, abs(gt[5]) * k)
        legend.update({"ramp": "terrain", "stretch": [lo, hi], "units": item.get("units", "m"),
                       "note": "hypsometric tint x hillshade (azimuth 315, altitude 45)",
                       "resampling": "block mean"})
    else:
        v = downsample(a, k, "continuous")
        fin = v[np.isfinite(v)]
        if item.get("log"):
            lo, hi = 0.0, float(fin.max()) if fin.size else 1.0
        else:
            lo, hi = ((float(np.percentile(fin, 2)), float(np.percentile(fin, 98)))
                      if fin.size else (0.0, 1.0))
        rgba = render_continuous(v, item.get("ramp", "viridis"), lo, hi, log=bool(item.get("log")))
        legend.update({"ramp": item.get("ramp", "viridis"), "stretch": [lo, hi],
                       "stretch_rule": "log10(1 + value), 0 to max" if item.get("log")
                       else "2nd to 98th percentile", "units": item.get("units", ""),
                       "resampling": "block mean"})
    os.makedirs(out_dir, exist_ok=True)
    png = os.path.join(out_dir, item["name"] + ".png")
    pgw = os.path.join(out_dir, item["name"] + ".pgw")
    leg = os.path.join(out_dir, item["name"] + "_legend.json")
    write_png(png, rgba)
    gt_s = (gt[0], gt[1] * k, 0.0, gt[3], 0.0, gt[5] * k)
    write_world_file(pgw, gt_s)
    h, w = rgba.shape[:2]
    xmin, xmax = gt[0], gt[0] + gt[1] * a.shape[1]
    ymax, ymin = gt[3], gt[3] + gt[5] * a.shape[0]
    legend.update({"png": os.path.basename(png), "world_file": os.path.basename(pgw),
                   "extent": [xmin, ymin, xmax, ymax], "crs_epsg": epsg,
                   "png_size": [w, h], "px_size_m": abs(gt[1]) * k})
    with open(leg, "w", encoding="utf-8") as f:
        json.dump(legend, f, indent=1)
    return {"name": item["name"], "title": item["title"], "png": png, "world_file": pgw,
            "legend_json": leg, "crs_epsg": epsg, "xmin": xmin, "ymin": ymin, "xmax": xmax,
            "ymax": ymax, "px_size_m": abs(gt[1]) * k, "source_tool": item.get("source_tool", ""),
            "source_raster": item["path"], "kind": "classified" if kind == "classified" else "continuous"}


def standard_items(raw_dem=None, accumulation=None, erosion_folder=None):
    """The A6 raster list from what a run has: relief, accumulation, erosion classes, STI."""
    items = []
    if raw_dem:
        items.append(dict(name="relief", title="Terrain (hillshade over elevation)", path=raw_dem,
                          kind="relief", source_tool="raw DEM", units="m"))
    if accumulation:
        items.append(dict(name="flow_accumulation", title="Flow accumulation (log stretch)",
                          path=accumulation, kind="continuous", ramp="blues", log=True,
                          source_tool="Flow accumulation", units="cells"))
    if erosion_folder:
        from ..erosion.io import FILES
        for key, title, kind, ramp in (
                ("spi_class", "SPI classes", "classified", None),
                ("rusle_class", "RUSLE soil loss classes", "classified", None),
                ("ls_class", "LS terrain potential classes (LS-only)", "classified", None),
                ("combined_class", "Combined erosion severity", "classified", None),
                ("sti", "Sediment transport index (overland)", "continuous", "ylorrd")):
            p = os.path.join(erosion_folder, FILES[key])
            if os.path.exists(p):
                items.append(dict(name=key, title=title, path=p, kind=kind, ramp=ramp,
                                  source_tool="Erosion indices and RUSLE"))
    return items
