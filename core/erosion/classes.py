# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Erosion severity classes (decision D7 + the A14 2025 scheme), WP-F v0.14.

A class scheme = increasing breaks, class names and a SEVERITY SCORE per
class on a common 1-5 scale (1 = lowest, 5 = most severe). Class rasters
store the class code (1..n, 0 = NoData); combinations and the crossing
MCDMA use the scores, so schemes with 4 and 5 classes can be combined.

Schemes shipped (all editable by passing custom breaks):
* spi_percentile (D7 default): ln(SPI) at the 50/75/90/97th percentiles of
  the analysed extent - RELATIVE classes (rank cells within this DEM).
* spi_a14: fixed ln(SPI) breaks 0 / 5 / 10 -> Low, Moderate, High, Severe
  (Akello & Omosa 2025, A14 corridor, AJERI; scores 1, 2, 3, 5).
* rusle_d7 (default): A t/ha/yr < 5 / 5-12 / 12-25 / 25-50 / > 50 ->
  Very low ... Very high. Context: tolerable soil loss in the central
  Kenyan highlands is 2.2-10 t/ha/yr (Angima et al. 2003).
* rusle_a14: 0-5 Slight, 5-10 Moderate, 10-20 High, 20-40 Very high,
  > 40 Severe (after Watene et al. 2021; Akello & Omosa 2025).
* ls_percentile: LS percentiles 50/75/90/97 - "terrain potential" for the
  LS-only fallback (not soil loss).
* volume_a14: sediment volume m3/yr < 1,000 Low, 1,000-5,000 Moderate,
  5,000-15,000 High, > 15,000 Severe (scores 1, 2, 3, 5).
Combined raster (D7): 5 x 5 matrix on [SPI score - 1][RUSLE-or-LS score - 1],
default = the higher of the two.
"""

import numpy as np

SCORE_COLOURS = {1: (26, 150, 65), 2: (166, 217, 106), 3: (255, 255, 191),
                 4: (253, 174, 97), 5: (215, 25, 28)}
FIVE = ["Very low", "Low", "Moderate", "High", "Very high"]
PERCENTILES = (50.0, 75.0, 90.0, 97.0)

SCHEMES = {
    "spi_percentile": {"breaks": None, "names": FIVE, "scores": [1, 2, 3, 4, 5],
                       "basis": "ln(SPI) percentiles 50/75/90/97 of this extent (relative)"},
    "spi_a14": {"breaks": (0.0, 5.0, 10.0), "names": ["Low", "Moderate", "High", "Severe"],
                "scores": [1, 2, 3, 5], "basis": "fixed ln(SPI) 0/5/10 (Akello & Omosa 2025)"},
    "rusle_d7": {"breaks": (5.0, 12.0, 25.0, 50.0), "names": FIVE, "scores": [1, 2, 3, 4, 5],
                 "basis": "RUSLE A t/ha/yr 5/12/25/50 (QEHT D7)"},
    "rusle_a14": {"breaks": (5.0, 10.0, 20.0, 40.0),
                  "names": ["Slight", "Moderate", "High", "Very high", "Severe"],
                  "scores": [1, 2, 3, 4, 5],
                  "basis": "RUSLE A t/ha/yr 5/10/20/40 (Watene et al. 2021; Akello & Omosa 2025)"},
    "ls_percentile": {"breaks": None, "names": FIVE, "scores": [1, 2, 3, 4, 5],
                      "basis": "LS percentiles 50/75/90/97 (terrain potential, not soil loss)"},
    "volume_a14": {"breaks": (1000.0, 5000.0, 15000.0), "names": ["Low", "Moderate", "High", "Severe"],
                   "scores": [1, 2, 3, 5], "basis": "sediment volume m3/yr 1,000/5,000/15,000 "
                   "(Akello & Omosa 2025)"},
}
DEFAULT_MATRIX = [[max(i, j) + 1 for j in range(5)] for i in range(5)]


class ClassScheme(object):
    def __init__(self, breaks, names, scores, basis):
        b = tuple(float(v) for v in breaks)
        if len(b) + 1 != len(names) or len(names) != len(scores):
            raise ValueError("a scheme needs n breaks, n + 1 names and n + 1 scores")
        if any(y <= x for x, y in zip(b, b[1:])):
            raise ValueError("class breaks must be strictly increasing")
        if min(scores) < 1 or max(scores) > 5:
            raise ValueError("severity scores must lie in 1-5")
        self.breaks, self.names, self.scores, self.basis = b, list(names), list(scores), basis

    def classify(self, x, valid):
        out = np.zeros(np.shape(x), dtype=np.uint8)
        ok = np.asarray(valid, dtype=bool) & np.isfinite(x)
        out[ok] = (np.searchsorted(np.asarray(self.breaks), np.asarray(x)[ok], side="right") + 1
                   ).astype(np.uint8)
        return out

    def score_of(self, codes):
        lut = np.zeros(len(self.names) + 1, dtype=np.uint8)
        lut[1:] = self.scores
        return lut[np.asarray(codes, dtype=np.int64)]

    def describe(self):
        return {"breaks": list(self.breaks), "names": self.names, "scores": self.scores,
                "basis": self.basis}


def percentile_breaks(x, valid, percentiles=PERCENTILES):
    vals = np.asarray(x)[np.asarray(valid, dtype=bool) & np.isfinite(x)]
    if vals.size == 0:
        raise ValueError("no valid cells to derive percentile breaks")
    b = np.percentile(vals, percentiles)
    for i in range(1, len(b)):                       # plateaus: keep breaks increasing
        if b[i] <= b[i - 1]:
            b[i] = np.nextafter(b[i - 1], np.inf)
    return tuple(float(v) for v in b)


def get_scheme(key, x=None, valid=None, breaks=None):
    """A ClassScheme by key; percentile schemes derive breaks from x; custom breaks override."""
    s = SCHEMES[key]
    if breaks is not None:
        return ClassScheme(breaks, s["names"], s["scores"], "user breaks (" + key + ")")
    if s["breaks"] is None:
        return ClassScheme(percentile_breaks(x, valid), s["names"], s["scores"], s["basis"])
    return ClassScheme(s["breaks"], s["names"], s["scores"], s["basis"])


def combine_scores(score_a, score_b, matrix=None):
    m = np.asarray(DEFAULT_MATRIX if matrix is None else matrix, dtype=np.int64)
    if m.shape != (5, 5) or m.min() < 1 or m.max() > 5:
        raise ValueError("combination matrix must be 5 x 5 with values 1-5")
    out = np.zeros(np.shape(score_a), dtype=np.uint8)
    ok = (score_a > 0) & (score_b > 0)
    out[ok] = m[score_a[ok].astype(np.int64) - 1, score_b[ok].astype(np.int64) - 1]
    return out


COMBINED = ClassScheme((1.5, 2.5, 3.5, 4.5), FIVE, [1, 2, 3, 4, 5],
                       "combined severity score (D7 matrix)")


def class_extents(codes, scheme, cell_area_m2, channel=None):
    """Area and share per class, split channel / hillslope when a mask is given."""
    total = int((codes > 0).sum())
    rows = []
    for k, name in enumerate(scheme.names, start=1):
        m = codes == k
        n = int(m.sum())
        row = {"class": k, "name": name, "score": scheme.scores[k - 1], "cells": n,
               "area_km2": n * cell_area_m2 / 1e6, "pct": 100.0 * n / total if total else 0.0}
        if channel is not None:
            row["channel_km2"] = int((m & channel).sum()) * cell_area_m2 / 1e6
            row["hillslope_km2"] = int((m & ~channel).sum()) * cell_area_m2 / 1e6
        rows.append(row)
    return rows


def _xml(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            .replace('"', "&quot;"))


def qml_paletted(scheme, title):
    """QGIS style (.qml): paletted renderer with the scheme's names and severity colours."""
    items = "\n".join(
        f'        <paletteEntry value="{k}" alpha="255" '
        f'color="#{"%02x%02x%02x" % SCORE_COLOURS[scheme.scores[k - 1]]}" label="{_xml(n)}"/>'
        for k, n in enumerate(scheme.names, start=1))
    return ("<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>\n"
            '<qgis styleCategories="Symbology" version="3.22">\n  <pipe>\n'
            '    <rasterrenderer type="paletted" band="1" opacity="1" alphaBand="-1" nodataColor="">\n'
            '      <rasterTransparency/>\n      <colorPalette>\n' + items +
            '\n      </colorPalette>\n    </rasterrenderer>\n  </pipe>\n'
            f'  <!-- {_xml(title)} -->\n</qgis>\n')


def write_class_raster(path, codes, scheme, geotransform, projection_wkt, title):
    """uint8 GeoTIFF, NoData 0, colour table, category names, RAT and a .qml sidecar."""
    from osgeo import gdal
    gdal.UseExceptions()
    rows, cols = codes.shape
    ds = gdal.GetDriverByName("GTiff").Create(path, cols, rows, 1, gdal.GDT_Byte,
                                              options=["COMPRESS=DEFLATE"])
    ds.SetGeoTransform(tuple(geotransform))
    if projection_wkt:
        ds.SetProjection(projection_wkt)
    band = ds.GetRasterBand(1)
    ct = gdal.ColorTable()                      # colour table before any pixel is written
    ct.SetColorEntry(0, (0, 0, 0, 0))
    for k, sc in enumerate(scheme.scores, start=1):
        ct.SetColorEntry(k, SCORE_COLOURS[sc] + (255,))
    band.SetRasterColorTable(ct)
    band.SetRasterColorInterpretation(gdal.GCI_PaletteIndex)
    band.SetNoDataValue(0)
    band.WriteArray(codes.astype(np.uint8))
    band.SetCategoryNames(["NoData"] + scheme.names)
    rat = gdal.RasterAttributeTable()
    rat.CreateColumn("Value", gdal.GFT_Integer, gdal.GFU_MinMax)
    rat.CreateColumn("Class", gdal.GFT_String, gdal.GFU_Name)
    rat.CreateColumn("Score", gdal.GFT_Integer, gdal.GFU_Generic)
    for i, name in enumerate(scheme.names):
        rat.SetValueAsInt(i, 0, i + 1)
        rat.SetValueAsString(i, 1, name)
        rat.SetValueAsInt(i, 2, scheme.scores[i])
    band.SetDefaultRAT(rat)
    band.SetDescription(title)
    ds.SetMetadataItem("QEHT_CLASS_BASIS", scheme.basis)
    ds.SetMetadataItem("QEHT_CLASS_BREAKS", ", ".join(f"{b:g}" for b in scheme.breaks))
    band = None
    ds = None
    with open(path.rsplit(".", 1)[0] + ".qml", "w", encoding="utf-8") as f:
        f.write(qml_paletted(scheme, title))
    return path
