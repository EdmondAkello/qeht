# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Land-cover scenario against the baseline (F14, v0.26).

A second land-cover raster with the same classes (for example a planned
development or a future land-use map) is run through the SAME lookups and
the SAME hydrologic soil groups as the baseline, so every difference comes
from land cover only:

    cn_ii_scn, cn_export_scn     curve number of the scenario
    d_cn                         cn_ii_scn - cn_ii
    rational_c_scn, d_rational_c Rational C (only with a user C lookup)
    lc_pct_<group>_scn           land-cover shares of the scenario
    ero_c_mean_scn               mean C from the WorldCover -> C lookup
    ero_a_mean_tha_scn           mean soil loss with that C: A_scn = A x C_scn / C per cell
                                 (A = R K LS C P, so R, K, LS and P are unchanged)
    d_ero_a_mean_tha             ero_a_mean_tha_scn - the baseline mean over the same cells

No QGIS imports, no GDAL.
"""

import json

import numpy as np

from .curve_number import RunoffInputs, LC_GROUPS

SCN_FIELDS = (["scenario_name", "cn_ii_scn", "cn_export_scn", "d_cn", "rational_c_scn",
               "d_rational_c"] + [f"lc_pct_{g}_scn" for g in LC_GROUPS]
              + ["ero_c_mean_scn", "ero_a_mean_tha_scn", "d_ero_a_mean_tha"])


class Scenario(object):
    def __init__(self, name, classes, baseline, soil_loss=None, c_grid=None, lc_dataset=""):
        """baseline: the RunoffInputs of the run (lookups and HSG are reused).
        soil_loss / c_grid: the erosion run's A and C grids (optional)."""
        self.name = name or "scenario"
        self.lc_dataset = lc_dataset
        self.runoff = RunoffInputs(classes, hsg=baseline.hsg, condition=baseline.condition,
                                   amc=baseline.amc, cn_lookup=baseline.cn_lookup,
                                   cn_lookup_path=baseline.cn_lookup_path,
                                   rc_lookup=baseline.rc_lookup_in,
                                   rc_lookup_path=baseline.rc_lookup_path, lc_dataset=lc_dataset)
        self.base = baseline
        self.a, self.c_scn, self.a_scn = None, None, None
        if soil_loss is not None and c_grid is not None:
            from ..erosion.rusle import c_from_worldcover
            self.c_scn, _ = c_from_worldcover(np.asarray(classes, float))
            a = np.asarray(soil_loss, float)
            c = np.asarray(c_grid, float)
            with np.errstate(divide="ignore", invalid="ignore"):
                self.a_scn = np.where(c > 0, a * self.c_scn / c, np.nan)
            self.a = np.where(np.isfinite(self.a_scn), a, np.nan)

    def block(self, mask):
        out = {k: None for k in SCN_FIELDS}
        out["scenario_name"] = self.name
        b = self.runoff.block(mask)
        base = self.base.block(mask)
        out["cn_ii_scn"], out["cn_export_scn"] = b.get("cn_ii"), b.get("cn_export")
        out["rational_c_scn"] = b.get("rational_c")
        if b.get("cn_ii") is not None and base.get("cn_ii") is not None:
            out["d_cn"] = b["cn_ii"] - base["cn_ii"]
        if b.get("rational_c") is not None and base.get("rational_c") is not None:
            out["d_rational_c"] = b["rational_c"] - base["rational_c"]
        for g in LC_GROUPS:
            out[f"lc_pct_{g}_scn"] = b.get(f"lc_pct_{g}")
        if self.c_scn is not None:
            m = np.asarray(mask, bool)
            c = self.c_scn[m]
            c = c[np.isfinite(c)]
            out["ero_c_mean_scn"] = float(c.mean()) if c.size else None
            s, a = self.a_scn[m], self.a[m]
            ok = np.isfinite(s) & np.isfinite(a)
            if ok.any():
                out["ero_a_mean_tha_scn"] = float(s[ok].mean())
                out["d_ero_a_mean_tha"] = float(s[ok].mean() - a[ok].mean())
        return out

    def meta_json(self, c_source=""):
        return json.dumps({"name": self.name, "lc_dataset": self.lc_dataset,
                           "cn_lookup_id": self.runoff.cn_lookup_name,
                           "rc_lookup_id": self.runoff.rc_lookup_name,
                           "hsg": "baseline hydrologic soil groups",
                           "erosion": ("A_scn = A x C_scn / C per cell; C_scn from the WorldCover -> "
                                       "C lookup" + (f"; baseline C: {c_source}" if c_source else ""))
                           if self.c_scn is not None else "",
                           "note": "the change comes from land cover only (same lookups and HSG)"},
                          sort_keys=True)
