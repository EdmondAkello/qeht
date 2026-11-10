# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Review register and crossing sheets from a design hydrology package (1.2)."""

import os

from qgis.core import (QgsProcessingParameterFile, QgsProcessingParameterFileDestination,
                       QgsProcessingParameterRasterLayer, QgsProcessingParameterString,
                       QgsProcessingParameterVectorDestination, QgsProcessing,
                       QgsProcessingException)

from .base import QehtAlgorithm
from .alg_review_depressions import add_review_parameters, review_settings, write_depressions


def package_review(package, raw=None, filled=None, valid=None, info=None, settings=None,
                   flat=None, warnings=(), dem_audit=None, feedback=None):
    """Depressions (when raw + filled are given), register rows. Shared by
    the package review tool and the pipeline. -> (register, depressions, summary)."""
    from ..core.interop import gpkg
    from ..core.review.depressions import depression_table
    from ..core.review.register import build_register
    from ..core.review.sheets import alignment_from_package
    settings = dict(settings or {})
    deps, summ = [], None
    if raw is not None and filled is not None:
        al = alignment_from_package(package)
        xs = [(c["outlet_uid"], c["chainage_m"]) for c in
              gpkg.read_table(package, "crossings", with_geometry=False)
              if c.get("chainage_m") is not None]
        deps, _, summ = depression_table(raw, filled, valid, info.geotransform, alignment=al,
                                         crossings=xs or None, **settings)
    reg = build_register(package, {"warnings": list(warnings), "dem_audit": dem_audit,
                                   "flat_check": flat or {}},
                         depressions=deps, near_m=settings.get("near_m", 100.0))
    return reg, deps, summ


class ReviewPackageAlgorithm(QehtAlgorithm):

    def name(self): return "reviewpackage"
    def displayName(self): return "Review register and crossing sheets"
    def group(self): return "Review"
    def groupId(self): return "review"

    def shortHelpString(self):
        return (
            "Builds the <b>review register</b> of a design hydrology package: every finding in "
            "one table with an ID (R-001 ...), severity (high / medium / low / info), "
            "category, crossing, chainage, the evidence and a suggested check, plus empty "
            "status and response columns for the review. It collects truncated catchments, "
            "catchments on the clip edge, crossings whose area or channel depends on the flat "
            "method, DEM-uncertainty switching, outlets off mapped drainage, far snaps, Tc "
            "outside calibration ranges, rising flow paths, missing USLE K and proposed "
            "crossings.\n\n"
            "Also writes the <b>crossing review sheets</b>: one printable page per crossing "
            "with a plan, the characteristics, the checks and its register findings, and a "
            "sign-off box. Print to PDF from a browser.\n\n"
            "Give the <b>DEM the fill ran on</b> and the <b>filled DEM</b> to add the "
            "flat-method check and the fill-depression review (deep fills at or near the road "
            "become findings; see Review fill depressions). Without them the register uses "
            "what the package holds.\n\n"
            "The one-click pipeline writes both files automatically.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile(
            "PACKAGE", "Design hydrology package", fileFilter="GeoPackage (*.gpkg)"))
        self.addParameter(QgsProcessingParameterRasterLayer(
            "DEM", "DEM the fill ran on (optional)", optional=True))
        self.addParameter(QgsProcessingParameterRasterLayer(
            "FILLED", "Filled DEM (optional)", optional=True))
        add_review_parameters(self)
        self.addParameter(QgsProcessingParameterString("TITLE", "Title", optional=True))
        self.addParameter(QgsProcessingParameterFileDestination(
            "REGISTER", "Review register", fileFilter="CSV (*.csv)"))
        self.addParameter(QgsProcessingParameterFileDestination(
            "SHEETS", "Crossing review sheets", fileFilter="HTML (*.html)"))
        self.addParameter(QgsProcessingParameterVectorDestination(
            "DEPRESSIONS", "Fill depressions", type=QgsProcessing.SourceType.TypeVectorPolygon,
            optional=True, createByDefault=False))

    def processAlgorithm(self, parameters, context, feedback):
        from ..core.interop.heas_exchange import validate_exchange
        from ..core.raster import read_dem
        from ..core.review.register import write_register_csv, register_counts
        from ..core.review.sheets import write_sheets
        pkg = self.parameterAsFile(parameters, "PACKAGE", context)
        errors, warnings = validate_exchange(pkg)
        if errors:
            raise QgsProcessingException("Not a valid design hydrology package: " + "; ".join(errors))
        has_dem = self.parameterAsRasterLayer(parameters, "DEM", context) is not None
        has_fill = self.parameterAsRasterLayer(parameters, "FILLED", context) is not None
        raw = filled = valid = info = None
        flat = {}
        if has_dem != has_fill:
            feedback.pushWarning("Give both the DEM and the filled DEM for the depression "
                                 "review; using the package only.")
        elif has_dem:
            raw, valid, info = read_dem(self.raster_path(parameters, "DEM", context))
            filled, fvalid, finfo = read_dem(self.raster_path(parameters, "FILLED", context))
            if (finfo.rows, finfo.cols) != (info.rows, info.cols):
                raise QgsProcessingException("The filled DEM is not on the DEM's grid.")
            valid = valid & fvalid
            self.check_size(feedback, info.rows, info.cols)
            flat = self._flat_check(pkg, filled, valid, finfo, feedback)
        s = review_settings(self, parameters, context)
        reg, deps, summ = package_review(pkg, raw, filled, valid, info, s, flat, warnings)
        out = self.parameterAsFileOutput(parameters, "REGISTER", context)
        write_register_csv(reg, out)
        n = register_counts(reg)
        feedback.pushInfo(f"Review register: {len(reg)} finding(s) - " + ", ".join(
            f"{v} {k}" for k, v in n.items() if v) + f" -> {out}")
        title = self.parameterAsString(parameters, "TITLE", context) or \
            os.path.splitext(os.path.basename(pkg))[0]
        sheets = self.parameterAsFileOutput(parameters, "SHEETS", context)
        write_sheets(pkg, sheets, reg, {"flat_check": flat}, deps, title=title)
        feedback.pushInfo(f"Crossing sheets: {sheets}")
        result = {"REGISTER": out, "SHEETS": sheets}
        dpath = self.parameterAsOutputLayer(parameters, "DEPRESSIONS", context)
        if dpath and deps:
            write_depressions(self, dpath, deps, info.projection_wkt)
            result["DEPRESSIONS"] = dpath
        if summ:
            feedback.pushInfo(f"Fill depressions: {summ['reported']} listed, "
                              f"{summ['flag_at_road']} against the road, "
                              f"{summ['flag_near_road']} near it.")
        return result

    @staticmethod
    def _flat_check(pkg, filled, valid, info, feedback):
        from ..core.interop import gpkg
        from ..core.flow.sensitivity import flat_method_sensitivity
        rows = gpkg.read_table(pkg, "crossings", with_geometry=False)
        md = {r["key"]: r["value"] for r in gpkg.read_table(pkg, "qeht_run_metadata")}
        gt = info.geotransform
        rc = [(int((c["outlet_y"] - gt[3]) // gt[5]), int((c["outlet_x"] - gt[0]) // gt[1]))
              for c in rows]
        if not rc or not all(0 <= r < info.rows and 0 <= c < info.cols for r, c in rc):
            feedback.pushWarning("Crossings fall outside the filled DEM; no flat-method check.")
            return {}
        try:
            km2 = float(md.get("stream_threshold_km2") or "")
            thr = max(1.0, round(km2 * 1e6 / (info.cell_width * info.cell_height)))
        except ValueError:
            thr = None
        res = flat_method_sensitivity(filled, valid, info.cell_width, info.cell_height, rc,
                                      stream_threshold_cells=thr)
        return {c["outlet_uid"]: r for c, r in zip(rows, res)}
