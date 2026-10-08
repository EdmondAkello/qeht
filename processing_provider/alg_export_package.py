# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Export a design hydrology package to KMZ and XLSX (F11, v0.23)."""

import os

from qgis.core import (QgsProcessingParameterFile, QgsProcessingParameterFileDestination,
                       QgsProcessingParameterString, QgsProcessingParameterBoolean,
                       QgsProcessingOutputString, QgsProcessingException)

from .base import QehtAlgorithm


def wgs84_transformer(package_path):
    """points in the package CRS -> (lon, lat), with OSR (in-process GDAL)."""
    from osgeo import ogr, osr
    ogr.UseExceptions()
    osr.UseExceptions()
    ds = ogr.Open(package_path)
    src = None
    for name in ("crossings", "catchments", "flowpaths"):
        lyr = ds.GetLayerByName(name)
        if lyr is not None and lyr.GetSpatialRef() is not None:
            src = lyr.GetSpatialRef().Clone()
            break
    ds = None
    if src is None:
        raise QgsProcessingException("The package layers have no CRS.")
    dst = osr.SpatialReference()
    dst.ImportFromEPSG(4326)
    for s in (src, dst):
        s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tr = osr.CoordinateTransformation(src, dst)

    def to_wgs84(points):
        pts = [[float(p[0]), float(p[1])] for p in points]
        return [tuple(q[:2]) for q in tr.TransformPoints(pts)] if pts else []
    return to_wgs84


def relief_quicklook(package_path):
    """(png path, (xmin, ymin, xmax, ymax)) of the relief quicklook, or None."""
    from ..core.interop import gpkg
    if "rasters" not in dict(gpkg.list_tables(package_path)):
        return None
    for r in gpkg.read_table(package_path, "rasters", with_geometry=False):
        if r.get("name") == "relief" and r.get("png"):
            p = r["png"]
            if not os.path.isabs(p):
                p = os.path.join(os.path.dirname(os.path.abspath(package_path)), p)
            if os.path.exists(p):
                return p, (r["xmin"], r["ymin"], r["xmax"], r["ymax"])
    return None


class ExportPackageAlgorithm(QehtAlgorithm):

    def name(self): return "exportpackage"
    def displayName(self): return "Export package to KMZ / XLSX"
    def group(self): return "Workflow"
    def groupId(self): return "workflow"

    def shortHelpString(self):
        return (
            "Writes a design hydrology package as a <b>KMZ</b> for Google Earth and as an "
            "<b>XLSX</b> workbook, both in pure Python (no extra software).\n\n"
            "<b>KMZ:</b> crossings (existing blue, proposed orange) with a card per crossing "
            "(location, catchment, flow path, time of concentration by every method, runoff, "
            "channel section, erosion), translucent catchments with labels, longest flow paths, "
            "the road with kilometre posts, coverage findings, sags and flat stretches, channel "
            "sections and side-drain siltation (hidden by default), the relief quicklook as a "
            "terrain overlay when the package has one, a run record, a legend and a title block "
            "(the legend and title need matplotlib, which ships with QGIS on most installs). "
            "Coordinates are transformed to WGS 84.\n\n"
            "<b>XLSX:</b> sheets Catchment characteristics (design order, units in row 2, panes "
            "frozen at B3), Time of concentration, Coverage findings, Run metadata and Field "
            "dictionary. The column order is one config list (core/report/xlsx_layout.py).\n\n"
            "Nothing is recalculated: every value is read from the package.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile(
            "PACKAGE", "Design hydrology package", fileFilter="GeoPackage (*.gpkg)"))
        self.addParameter(QgsProcessingParameterString("TITLE", "Title (KMZ document and title "
                                                       "block)", optional=True))
        self.addParameter(QgsProcessingParameterBoolean(
            "RELIEF", "KMZ: include the relief quicklook as a terrain overlay (if present)",
            defaultValue=True))
        self.addParameter(QgsProcessingParameterFileDestination(
            "KMZ", "KMZ", fileFilter="KMZ (*.kmz)", optional=True))
        self.addParameter(QgsProcessingParameterFileDestination(
            "XLSX", "XLSX workbook", fileFilter="Excel workbook (*.xlsx)", optional=True))
        self.addOutput(QgsProcessingOutputString("SUMMARY", "Summary"))

    def processAlgorithm(self, parameters, context, feedback):
        from ..core.interop.heas_exchange import validate_exchange
        pkg = self.parameterAsFile(parameters, "PACKAGE", context)
        errors, _ = validate_exchange(pkg)
        if errors:
            raise QgsProcessingException("Not a valid design hydrology package: " + "; ".join(errors))
        title = (self.parameterAsString(parameters, "TITLE", context) or "").strip() or \
            os.path.splitext(os.path.basename(pkg))[0]
        kmz = self.parameterAsFileOutput(parameters, "KMZ", context)
        xlsx = self.parameterAsFileOutput(parameters, "XLSX", context)
        if not kmz and not xlsx:
            raise QgsProcessingException("Choose a KMZ, an XLSX or both.")
        res, done = {}, []
        if kmz:
            from ..core.report.kmz import write_kmz
            if not kmz.lower().endswith(".kmz"):
                kmz += ".kmz"
            relief = relief_quicklook(pkg) if self.parameterAsBool(parameters, "RELIEF", context) else None
            write_kmz(pkg, kmz, wgs84_transformer(pkg), title=title, relief=relief)
            res["KMZ"] = kmz
            done.append(f"KMZ {kmz}" + (" (with relief)" if relief else ""))
        if xlsx:
            from ..core.report.xlsx_layout import write_workbook
            if not xlsx.lower().endswith(".xlsx"):
                xlsx += ".xlsx"
            write_workbook(pkg, xlsx, title=title)
            res["XLSX"] = xlsx
            done.append(f"XLSX {xlsx}")
        res["SUMMARY"] = "; ".join(done)
        feedback.pushInfo(res["SUMMARY"])
        return res
