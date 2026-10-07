# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Run report from a design hydrology package (F8, v0.17)."""

import os

from qgis.core import (QgsProcessingParameterFile, QgsProcessingParameterFileDestination,
                       QgsProcessingParameterString, QgsProcessingParameterBoolean,
                       QgsProcessingException)

from .base import QehtAlgorithm


class RunReportAlgorithm(QehtAlgorithm):

    def name(self): return "runreport"
    def displayName(self): return "Run report from a design hydrology package"
    def group(self): return "Workflow"
    def groupId(self): return "workflow"

    def shortHelpString(self):
        return ("Writes a one-page HTML report for a design hydrology package: summary cards, a "
                "schematic plan (catchments, flow paths, road, crossings), the crossing "
                "schedule, the drainage coverage findings, inputs and provenance (DEM source "
                "and hash, conditioning, flat method, thresholds, soil / curve-number / rainfall "
                "/ erosion sources with PROXY and ESTIMATE labels) and the curve number lookup. "
                "Everything is read from the package; nothing is recalculated. Optionally also "
                "writes the catchment characteristics table (one row per crossing) as CSV.\n\n"
                "The one-click pipeline writes this report automatically, with its DEM checks, "
                "flat-method check and warnings added. Print the page to PDF from a browser.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile(
            "PACKAGE", "Design hydrology package", fileFilter="GeoPackage (*.gpkg)"))
        self.addParameter(QgsProcessingParameterString("TITLE", "Report title", optional=True))
        self.addParameter(QgsProcessingParameterBoolean(
            "CSV", "Also write the catchment characteristics table (CSV, next to the report)",
            defaultValue=True))
        self.addParameter(QgsProcessingParameterFileDestination(
            "OUTPUT", "Run report", fileFilter="HTML (*.html)"))

    def processAlgorithm(self, parameters, context, feedback):
        from ..core.interop.heas_exchange import validate_exchange
        from ..core.report.run_report import write_report
        from ..core.report.characteristics import write_characteristics_csv
        pkg = self.parameterAsFile(parameters, "PACKAGE", context)
        errors, warnings = validate_exchange(pkg)
        if errors:
            raise QgsProcessingException("Not a valid design hydrology package: " + "; ".join(errors))
        out = self.parameterAsFileOutput(parameters, "OUTPUT", context)
        title = self.parameterAsString(parameters, "TITLE", context) or \
            os.path.splitext(os.path.basename(pkg))[0]
        write_report(pkg, out, title=title, extra={"warnings": warnings} if warnings else None)
        result = {"OUTPUT": out}
        if self.parameterAsBool(parameters, "CSV", context):
            csv_path, n = write_characteristics_csv(
                pkg, os.path.splitext(out)[0] + "_characteristics.csv")
            feedback.pushInfo(f"Characteristics table: {n} crossing(s) -> {csv_path}")
            result["CSV"] = csv_path
        feedback.pushInfo(f"Report: {out}")
        return result
