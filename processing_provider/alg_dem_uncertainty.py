# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""DEM uncertainty at crossings (F15, v0.27)."""

import os

import numpy as np
from qgis.core import (QgsProcessing, QgsProcessingException, QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterFeatureSource, QgsProcessingParameterField,
                       QgsProcessingParameterNumber, QgsProcessingParameterVectorDestination,
                       QgsProcessingParameterFileDestination, QgsProcessingOutputString)

from .base import QehtAlgorithm
from ..core.raster import read_dem
from ..core.watershed import uncertainty as un


class DemUncertaintyAlgorithm(QehtAlgorithm):

    def name(self): return "demuncertainty"
    def displayName(self): return "DEM uncertainty at crossings"
    def group(self): return "Road drainage"
    def groupId(self): return "roaddrainage"

    def shortHelpString(self):
        return (
            "Turns each crossing's catchment area, flow-path length, 10-85 slope and Kirpich Tc "
            "into a range under DEM error (Monte Carlo).\n\n"
            "For each realisation a spatially correlated Gaussian error (standard deviation "
            "sigma, correlation exp(-r^2/L^2)) is added to the raw DEM, which is filled, routed "
            "(D8, Barnes) and accumulated; each crossing's pour point is snapped to the stream "
            "within the snap radius, or recorded as lost.\n\n"
            "<b>Outputs per crossing:</b> P10, P50, P90 and the coefficient of variation of each "
            "value (unc_*), the share of realisations where the outlet lost its stream and where "
            "the area moved by more than 25 % (catchment switching).\n\n"
            "<b>Vertical error:</b> no default. Presets with their sources: AW3D30 4.4 m "
            "(Tadono et al. 2016, 4.40 m RMSE at 5,121 check points); FABDEM 2.5 m, an ESTIMATE "
            "from Hawker et al. 2022 (mean absolute error 1.1 m built-up to 2.9 m forest); or a "
            "custom value. Correlation length 90 m, 50 realisations and seed 2026 by default; "
            "the seed is recorded so the run repeats exactly.\n\n"
            "<b>Run time:</b> one fill, D8 and accumulation per realisation; the time per "
            "realisation is reported before the loop. Clip a large DEM to the road's "
            "contributing area first.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer("DEM", "Raw DEM (projected, metres)"))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "CROSSINGS", "Crossings (candidates, a package's crossings or pour points)",
            [QgsProcessing.SourceType.TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterField(
            "ID_FIELD", "ID field (e.g. outlet_uid)", parentLayerParameterName="CROSSINGS",
            optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            "THRESHOLD", "Stream threshold (cells)", QgsProcessingParameterNumber.Type.Double,
            defaultValue=200.0, minValue=1.0))
        self.addParameter(QgsProcessingParameterNumber(
            "SNAP", "Snap radius (cells)", QgsProcessingParameterNumber.Type.Integer,
            defaultValue=5, minValue=1))
        self.add_uncertainty_parameters(with_switch=False)
        self.addParameter(QgsProcessingParameterVectorDestination(
            "OUTPUT", "Crossings with DEM uncertainty", QgsProcessing.SourceType.TypeVectorPoint))
        self.addParameter(QgsProcessingParameterFileDestination(
            "LOG", "Settings and correlation-length sensitivity (JSON)", fileFilter="JSON (*.json)",
            optional=True))
        self.addOutput(QgsProcessingOutputString("SUMMARY", "Summary"))

    def processAlgorithm(self, parameters, context, feedback):
        from qgis.core import QgsCoordinateReferenceSystem, QgsCoordinateTransform, QgsProject
        sigma, preset, src, corr, n, seed, sens = self.uncertainty_settings(parameters, context)
        z, valid, info = read_dem(self.raster_path(parameters, "DEM", context))
        self.check_size(feedback, info.rows, info.cols)
        pts = self.parameterAsSource(parameters, "CROSSINGS", context)
        idf = self.field_parameter(parameters, "ID_FIELD", context) if parameters.get("ID_FIELD") else ""
        dcrs = QgsCoordinateReferenceSystem()
        dcrs.createFromWkt(info.projection_wkt)
        tr = QgsCoordinateTransform(pts.sourceCrs(), dcrs, QgsProject.instance()) \
            if dcrs.isValid() and pts.sourceCrs() != dcrs else None
        gt = info.geotransform
        rows_, outlets = [], []
        for f in pts.getFeatures():
            g = f.geometry()
            if g is None or g.isEmpty():
                continue
            if tr is not None:
                g.transform(tr)
            p = g.asPoint()
            c, r = int((p.x() - gt[0]) // gt[1]), int((p.y() - gt[3]) // gt[5])
            if not (0 <= r < info.rows and 0 <= c < info.cols):
                continue
            uid = str(f[idf]) if idf else str(f.id())
            rows_.append(((p.x(), p.y()), {"outlet_uid": uid}))
            outlets.append((r, c))
        if not outlets:
            raise QgsProcessingException("No crossing lies on the DEM.")
        thr = self.parameterAsDouble(parameters, "THRESHOLD", context)
        snap = self.parameterAsInt(parameters, "SNAP", context)
        dem = np.where(valid, z, 0.0)
        if sens:
            feedback.pushInfo("The correlation-length sensitivity adds 2 x N realisations "
                              "(about three times the run time); switch it off to save time.")
        res, inf = un.run(dem, valid, gt, outlets, sigma, corr, n, seed, thr, snap,
                          progress=lambda f_: feedback.setProgress(100.0 * f_), log=feedback.pushInfo)
        for (_, a), b in zip(rows_, res):
            a.update(b)
        sensitivity, uids = None, None
        if sens:
            areas = [(b["unc_area_p50"] or 0) for b in res]
            big = sorted(range(len(res)), key=lambda k: -areas[k])[:3]
            uids = [rows_[k][1]["outlet_uid"] for k in big]
            sensitivity = un.corr_sensitivity(dem, valid, gt, [outlets[k] for k in big], sigma, corr,
                                              n, seed, thr, snap)
        fields = [("outlet_uid", "text")] + [(k, "text" if k == "unc_note" else
                                              ("int" if k == "unc_n" else "real"))
                  for k in un.UNC_FIELDS]
        out = self.parameterAsOutputLayer(parameters, "OUTPUT", context)
        self.write_vector(out, "dem_uncertainty", info.projection_wkt, "point", fields, rows_)
        log = self.parameterAsFileOutput(parameters, "LOG", context) or \
            os.path.splitext(out)[0] + "_uncertainty.json"
        with open(log, "w", encoding="utf-8") as fh:
            fh.write(un.params_json(inf, preset, src, sensitivity, uids))
        n_sw = sum(1 for b in res if (b["unc_switch_pct"] or 0) > 0)
        summary = (f"{len(res)} crossing(s), sigma {sigma:g} m ({preset}), L {corr:g} m, N {n}, "
                   f"seed {seed}; {n_sw} with catchment switching")
        feedback.pushInfo(summary)
        return {"OUTPUT": out, "LOG": log, "SUMMARY": summary}
