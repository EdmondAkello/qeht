# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Build design hydrology package - one self-describing GeoPackage per run.

Shown as "Build design hydrology package" since 0.15 (was "Build HEAS exchange
package"); the algorithm id buildheasexchange and the schema qeht-heas-1 are
unchanged so saved models and HEAS imports keep working.
"""

import os

from qgis.core import (
    QgsProcessing, QgsProcessingException, QgsProcessingContext,
    QgsProcessingParameterRasterLayer, QgsProcessingParameterFeatureSource,
    QgsProcessingParameterField, QgsProcessingParameterString,
    QgsProcessingParameterEnum, QgsProcessingParameterNumber,
    QgsProcessingParameterBoolean, QgsProcessingParameterFileDestination,
    QgsProcessingOutputString,
    QgsCoordinateReferenceSystem,
)
import numpy as np

from .base import QehtAlgorithm
from ..core.raster import read_dem
from ..core.grid import decode_d8
from ..core.watershed.delineate import extract_streams
from ..core.network.crossings import select_crossings, CANDIDATE_FIELDS
from ..core.interop.heas_exchange import (build_exchange_records, write_exchange,
                                          validate_exchange, crs_check,
                                          file_fingerprint, plugin_version,
                                          ExchangeError, SCHEMA_VERSION)

FDR = "FDR"; DEM = "DEM"; RAW_DEM = "RAW_DEM"; FAC = "FAC"; ORDER = "ORDER"
ROAD = "ROAD"; START = "START"; REVERSE = "REVERSE"
POINTS = "POINTS"; ID_FIELD = "ID_FIELD"; ID_PREFIX = "ID_PREFIX"; ID_ORDER = "ID_ORDER"
SNAP = "SNAP"; SNAP_THRESHOLD = "SNAP_THRESHOLD"; LOCAL = "LOCAL"
DEM_SOURCE = "DEM_SOURCE"; CONDITIONING = "CONDITIONING"; FLAT_METHOD = "FLAT_METHOD"
CSV = "CSV"; OUTPUT = "OUTPUT"; SUMMARY = "SUMMARY"

ORDER_OPTIONS = ["Automatic: along the chainage when known, else downstream first",
                 "Downstream first (largest contributing area = 001)",
                 "Pour-point layer order"]
ORDER_KEYS = ["auto", "downstream", "input"]


class BuildHeasExchangeAlgorithm(QehtAlgorithm):

    def name(self): return "buildheasexchange"
    def displayName(self): return "Build design hydrology package"
    def group(self): return "Interoperability"
    def groupId(self): return "interop"

    def shortHelpString(self):
        return (
            "One self-describing GeoPackage of crossings, catchments and flow "
            "paths linked by outlet_uid, for spreadsheets, reports or design "
            "software. Runs snapping, catchment delineation, longest flow path and "
            "catchment characteristics for every pour point.\n\n"
            f"(Formerly 'Build HEAS exchange package'. HEAS, a separate design tool, "
            f"imports the package directly - schema {SCHEMA_VERSION} - but it is not "
            "required.)\n\n"
            "<b>Layers:</b> crossings (snapped outlets), catchments, flowpaths, "
            "plus qeht_run_metadata (DEM, method, thresholds, QEHT version) and "
            "qeht_field_dictionary (meaning, unit and method of every field).\n\n"
            "<b>outlet_uid</b> links the three layers. It is assigned once per "
            "run: taken from the ID attribute you choose (e.g. existing culvert "
            "numbers), or sequential with your prefix (X001, X002 ...), numbered "
            "downstream-first by default. Duplicate or empty IDs stop the run "
            "with a list. <i>outlet_id</i> (the pour-point feature id) is still "
            "written but is NOT stable - do not link on it.\n\n"
            "<b>A projected, metric CRS is required.</b> A DEM in degrees is "
            "refused; reproject it (e.g. to UTM) first.\n\n"
            "<b>Crossing candidates:</b> give the layer from 'Road crossing "
            "candidates' as the crossings. The accepted candidates are used (or the "
            "recommended ones if none is accepted); each keeps its outlet cell "
            "(no snapping) and its chainage, IDs are numbered along the chainage, "
            "and the full candidate layer is stored in the package as "
            "crossing_candidates. With hand-placed points, give the road alignment "
            "to get chainages.\n\n"
            "<b>Catchments:</b> full upstream area per crossing by default "
            "(overlapping - what a culvert design flow needs). 'Local' gives "
            "non-overlapping areas, computed upstream-first so the result does "
            "not depend on the order of the pour-point layer.\n\n"
            "<b>Slopes (four domains, never merged):</b> catch_slope_horn (Horn "
            "mean), catch_relief_ratio (relief / LFP length; a basin-steepness index, not a Tc "
            "input), lfp_slope (drop / "
            "length) and lfp_slope_1085 (10-85 along the LFP, measured from the "
            "OUTLET; lfp_L10_m, lfp_L85_m, lfp_z10_m, lfp_z85_m let you check it "
            "by hand).\n\n"
            "<b>Flow path segments:</b> each flow path is split where it reaches the "
            "snap-to-stream threshold (the channel head): lfp_overland_m / _slope, "
            "lfp_channel_m / _slope / _slope_1085, and the overland part as sheet flow "
            "(up to the cap, default 100 m) plus shallow concentrated flow - inputs for "
            "Kerby and segmental (TR-55) Tc. <b>Basin shape:</b> perimeter, form factor, "
            "elongation and circularity ratios, drainage density, stream frequency and "
            "highest Strahler order, for the report. <b>Channel slopes at each crossing:</b> "
            "ch_slope_us (main stem upstream) and ch_slope_ds (D8 path downstream) over "
            "200 m by default. With an erosion folder, each crossing also gets the "
            "deposition indicator (ero_dep_ratio / ero_dep_flag) and the local overland "
            "STI; with a road too, the side-drain siltation indicator per side along the "
            "corridor (layer corridor_sti; relative classes, screening only).\n\n"
            "<b>Curve number and Rational C (optional):</b> give a land-cover raster "
            "(ESA WorldCover classes) with a soil source: per cell land cover x hydrologic "
            "soil group -> CN from TR-55 Table 2-2 (a PROXY match to WorldCover, condition "
            "fair/good/poor; replace it with your own lookup CSV), averaged over the "
            "catchment (cn_ii) and exported at the chosen AMC (cn_export). Rational C needs "
            "your lookup CSV - none ships. Land-cover shares lc_pct_* are always given.\n\n"
            "<b>Floodplain width indicator (with a road):</b> at crossings of at least 10 km2 "
            "the width along the road where the ground is below the bed + 0.5, 1 and 2 m "
            "(fp_w_*), from the alignment profile, and the same by height above nearest "
            "drainage (fp_hand_w_*). A terrain indicator for the split of the check flood "
            "between the main structure and relief culverts - not a flood level.\n\n"
            "<b>Channel section at each crossing:</b> a transect of the raw DEM perpendicular to "
            "the flow, 30 m downstream of the outlet by default: bed level, bank-full width and "
            "depth by break of slope, widths below bed + 0.5, 1 and 2 m, side slopes H:V and the "
            "station-elevation list (xs_*; layer xs_transects), with a quality flag. Indicative "
            "only; on a 30 m DEM a small channel is below the grid resolution. Use survey where "
            "available.\n\n"
            "<b>Check against mapped drainage (optional):</b> give mapped waterways (OSM, "
            "HydroRIVERS or a national layer). DEM streams at the comparison threshold are "
            "compared with them within the tolerance: precision, recall and F1 overall "
            "(metadata) and per catchment (map_precision, map_recall, map_f1); per crossing the "
            "distance to the nearest mapped waterway and whether it agrees (map_agrees); mapped "
            "rivers crossing the road without a crossing within 50 m become coverage findings; "
            "DEM stream runs that leave the mapped course form the layer drainage_divergence.\n\n"
            "<b>Time of concentration:</b> Kirpich, Kerby + Kirpich, SCS lag, TR-55 segments "
            "and Bransby-Williams side by side (tc_*_min), each with a validity flag from its "
            "published calibration range and its inputs in tc_basis_json. QEHT does not pick a "
            "method. TR-55 needs the 2-yr 24-h rainfall P2 (no default); sheet-flow n and Kerby N "
            "come from the land cover (PROXY lookups from TR-55 Table 3-1 and Kerby 1959) unless "
            "you give your own CSV; the TR-55 channel velocity uses Manning on the channel "
            "section (n 0.035 by default).\n\n"
            "<b>Rainfall (optional):</b> rainfall zone polygons give each catchment its "
            "dominant zone (rain_zone, rain_zone_pct, all shares in rain_zones_json; a "
            "warning when the dominant zone covers less than 80 %); a mean annual rainfall "
            "raster (mm/yr) gives map_mm; optionally rusle_r is ESTIMATED from it with a "
            "published R-P relation (Renard & Freimund 1994 or Lo et al. 1985).\n\n"
            "<b>Coverage check (with a road alignment):</b> every place a stream of at least "
            "the stream-threshold area crosses the road with no crossing within 50 m becomes "
            "a PROPOSED crossing (status = proposed, IDs P001...), delineated and "
            "characterised like the others; low points of the ground profile (sags) with "
            "enough local area against the embankment become proposed crossings too; flat "
            "stretches where relief culverts for sheet flow may be needed are listed. "
            "Layers coverage_check, sag_points, flat_stretches. Adopt, move or delete the "
            "proposed crossings, then run 'Renumber and relink'.\n\n"
            "<b>Soils (optional), from any source:</b> soil map polygons (with the SOTWIS / "
            "SOTER SQLite database, a CSV table, SOTWIS or plain sand/silt/clay/oc fields, "
            "or your own field names under the advanced parameters); a soil unit raster "
            "+ CSV (e.g. HWSD v2); texture rasters or a SoilGrids folder (depth-weighted); "
            "plus an optional hydrologic soil group raster (e.g. HYSOGs250m) or field and "
            "K raster or field that override the texture-based values. Adds the soil "
            "block to every catchment: "
            "texture, organic carbon, coarse fragments, USLE K (Williams/EPIC), a "
            "texture-based hydrologic-group PROXY and the share of the catchment "
            "covered. Default depth 0-20 cm.\n\n"
            "Supply the RAW DEM for reported elevations and slopes; routing uses "
            "the flow-direction grid. Everything runs in-process.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(FDR, "Flow direction (D8-coded)"))
        self.addParameter(QgsProcessingParameterRasterLayer(FAC, "Flow accumulation"))
        self.addParameter(QgsProcessingParameterRasterLayer(
            RAW_DEM, "Raw DEM (for reported elevations and slopes)"))
        self.addParameter(QgsProcessingParameterRasterLayer(
            ORDER, "Stream order raster (optional, for crossing attributes)", optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource(
            POINTS, "Crossings / pour points", [QgsProcessing.SourceType.TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterFeatureSource(
            ROAD, "Road alignment (optional; chainage for hand-placed points)",
            [QgsProcessing.SourceType.TypeVectorLine], optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            START, "Start chainage (m)", QgsProcessingParameterNumber.Type.Double,
            defaultValue=0.0))
        self.addParameter(QgsProcessingParameterBoolean(
            REVERSE, "Reverse chainage direction", defaultValue=False))
        self.addParameter(QgsProcessingParameterField(
            ID_FIELD, "ID attribute (optional; blank = sequential IDs)",
            parentLayerParameterName=POINTS, optional=True))
        self.addParameter(QgsProcessingParameterString(
            ID_PREFIX, "ID prefix (sequential IDs only; ignored with an ID attribute)",
            defaultValue="X", optional=True))
        self.addParameter(QgsProcessingParameterEnum(
            ID_ORDER, "Sequential numbering order", options=ORDER_OPTIONS, defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber(
            SNAP, "Snap radius (cells)", QgsProcessingParameterNumber.Type.Integer,
            defaultValue=5, minValue=0, maxValue=100))
        self.addParameter(QgsProcessingParameterNumber(
            SNAP_THRESHOLD, "Snap-to-stream threshold (cells; 0 = max accumulation)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=200.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterBoolean(
            LOCAL, "Non-overlapping (local) catchments", defaultValue=False))
        self.addParameter(QgsProcessingParameterString(
            DEM_SOURCE, "DEM source (recorded in metadata, e.g. 'ALOS 30 m')", optional=True))
        self.addParameter(QgsProcessingParameterString(
            CONDITIONING, "Conditioning used (recorded, e.g. 'fill, min_slope 0')",
            optional=True))
        self.addParameter(QgsProcessingParameterString(
            FLAT_METHOD, "Flat resolution used for flow direction (recorded, e.g. 'barnes')",
            optional=True))
        self.addParameter(QgsProcessingParameterBoolean(
            "COVERAGE", "Coverage check with a road alignment: missing crossings, sags and flat "
            "stretches; delineate proposed crossings", defaultValue=True))
        for key, label, default in (
                ("COV_MIN_AREA", "Coverage: smallest stream area to need a crossing (km2; 0 = stream threshold)", 0.0),
                ("COV_SEARCH", "Coverage: an existing crossing within this chainage covers a stream (m)", 50.0),
                ("COV_MERGE", "Coverage: merge stream stations within (m)", 30.0),
                ("COV_SMALL", "Coverage: flag existing crossings smaller than (km2)", 0.01),
                ("SAG_SMOOTH", "Sags: smoothing window along the profile (m)", 30.0),
                ("SAG_DEPTH", "Sags: minimum depth (m)", 0.3),
                ("SAG_MIN_AREA", "Sags: area that makes a sag a proposed crossing (km2)", 0.05),
                ("FLAT_SLOPE", "Flat stretches: slope below (%)", 0.5),
                ("FLAT_CROSSFALL", "Flat stretches: cross-fall measured over (m)", 100.0),
                ("FLAT_MIN_LEN", "Flat stretches: minimum length (m)", 300.0)):
            self.addParameter(self._advanced(QgsProcessingParameterNumber(
                key, label, QgsProcessingParameterNumber.Type.Double, defaultValue=default,
                minValue=0.0)))
        self.addParameter(self._advanced(QgsProcessingParameterString(
            "PROPOSED_PREFIX", "Prefix for proposed crossings", defaultValue="P", optional=True)))
        self.addParameter(QgsProcessingParameterNumber(
            "CH_SLOPE_DIST", "Approach / exit channel length for the crossing channel slopes (m)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=200.0, minValue=1.0))
        self.addParameter(QgsProcessingParameterNumber(
            "SHEET_CAP", "Sheet-flow cap within the overland part of the flow path (m)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=100.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterRasterLayer(
            "FILLED", "Filled DEM (optional; ponding depth in the alignment profile)",
            optional=True))
        from qgis.core import QgsProcessingParameterFolderDestination as _FD
        self.addParameter(_FD("QUICKLOOKS", "Raster quicklooks folder (PNG + world file + legend; "
                              "optional)", optional=True, createByDefault=False))
        self.addParameter(self._advanced(QgsProcessingParameterNumber(
            "QL_MAX_PX", "Quicklook size, long side (pixels)",
            QgsProcessingParameterNumber.Type.Integer, defaultValue=4096, minValue=64)))
        self.addParameter(self._advanced(QgsProcessingParameterNumber(
            "FP_MIN_AREA", "Floodplain width indicator: crossings of at least (km2)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=10.0, minValue=0.0)))
        self.addParameter(self._advanced(QgsProcessingParameterBoolean(
            "FP_HAND", "Floodplain width indicator: also by height above nearest drainage (HAND)",
            defaultValue=True)))
        self.addParameter(self._advanced(QgsProcessingParameterBoolean(
            "XS", "Channel cross-section downstream of each crossing (indicative)",
            defaultValue=True)))
        for key, label, default in (
                ("XS_DIST", "Channel section: distance downstream of the outlet (m)", 30.0),
                ("XS_HALF", "Channel section: half width of the transect (m)", 150.0),
                ("XS_BANK_SLOPE", "Channel section: bank top where the side slope falls below (m/m)",
                 0.05)):
            self.addParameter(self._advanced(QgsProcessingParameterNumber(
                key, label, QgsProcessingParameterNumber.Type.Double, defaultValue=default,
                minValue=0.0)))
        self.addParameter(QgsProcessingParameterNumber(
            "PROFILE_STEP", "Alignment profile station spacing (m)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=10.0, minValue=0.5))
        self.add_soil_parameters("SOIL", optional=True)
        self.add_runoff_parameters()
        self.add_rainfall_parameters()
        from qgis.core import QgsProcessingParameterFeatureSource as _FS, QgsProcessing as _QP
        self.addParameter(_FS(
            "BURN_LOG", "Breach log from 'Burn crossings through embankments' (optional; "
            "recorded in the package)", [_QP.SourceType.TypeVectorLine], optional=True))
        from qgis.core import QgsProcessingParameterFile as _PF
        self.addParameter(_PF(
            "EROSION", "Erosion output folder (optional; adds the ero_* erosion block)",
            behavior=_PF.Behavior.Folder, optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "MAPPED", "Mapped waterways for the drainage check (lines, e.g. OSM or HydroRIVERS; "
            "optional)", [QgsProcessing.SourceType.TypeVectorLine], optional=True))
        self.addParameter(QgsProcessingParameterField(
            "MAPPED_NAME_FIELD", "Mapped waterways: name field", parentLayerParameterName="MAPPED",
            optional=True))
        self.addParameter(self._advanced(QgsProcessingParameterString(
            "MAPPED_SOURCE", "Mapped waterways: source for the record (e.g. OSM 2026-09)",
            optional=True)))
        for key, label, default in (
                ("MAPPED_KM2", "Mapped drainage: comparison threshold (km2; about 1 for OSM, 10 for "
                 "HydroRIVERS)", 1.0),
                ("MAPPED_TOL", "Mapped drainage: tolerance (m)", 60.0),
                ("MAPPED_MIN_CELLS", "Mapped drainage: shortest divergence reach (cells)", 10.0),
                ("MAPPED_FAR", "Mapped drainage: no comparison beyond (m from the outlet)", 1000.0),
                ("CSTI_HALF", "Side-drain siltation: corridor half-width (m)", 50.0),
                ("CSTI_OFFSET", "Side-drain siltation: sampling interval across the road (m)", 10.0),
                ("CSTI_ROAD_HALF", "Side-drain siltation: road strip half-width (m)", 5.0),
                ("CSTI_DRAIN_SLOPE", "Side-drain siltation: flag below this longitudinal slope (%)",
                 1.0)):
            self.addParameter(self._advanced(QgsProcessingParameterNumber(
                key, label, QgsProcessingParameterNumber.Type.Double, defaultValue=default,
                minValue=0.0)))
        self.add_uncertainty_parameters()
        self.addParameter(QgsProcessingParameterNumber(
            "TC_P2", "Time of concentration: 2-yr 24-h rainfall P2 for TR-55 (mm; no default)",
            QgsProcessingParameterNumber.Type.Double, optional=True, minValue=0.0))
        from qgis.core import QgsProcessingParameterRasterLayer as _RL2, QgsProcessingParameterFile as _PF2
        self.addParameter(_RL2("TC_P2_RASTER", "Time of concentration: P2 raster (mm; mean along "
                               "the flow path; optional)", optional=True))
        self.addParameter(self._advanced(QgsProcessingParameterBoolean(
            "TC", "Time of concentration by five methods", defaultValue=True)))
        self.addParameter(self._advanced(QgsProcessingParameterNumber(
            "TC_CHANNEL_N", "Time of concentration: channel Manning n (TR-55 channel part)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=0.035, minValue=0.001)))
        for key, label in (("TC_SHEET_N_CSV", "Time of concentration: sheet-flow n lookup CSV "
                            "(class, value; replaces the PROXY)"),
                           ("TC_KERBY_N_CSV", "Time of concentration: Kerby N lookup CSV "
                            "(class, value; replaces the PROXY)")):
            self.addParameter(self._advanced(_PF2(key, label, optional=True,
                                                  fileFilter="CSV (*.csv)")))
        self.addParameter(QgsProcessingParameterBoolean(
            CSV, "Also write one CSV per layer (for spreadsheets)", defaultValue=False))
        self.addParameter(QgsProcessingParameterFileDestination(
            OUTPUT, "Exchange GeoPackage", fileFilter="GeoPackage (*.gpkg)"))
        self.addOutput(QgsProcessingOutputString(SUMMARY, "Summary"))

    def processAlgorithm(self, parameters, context, feedback):
        fdr_path = self.raster_path(parameters, FDR, context)
        fac_path = self.raster_path(parameters, FAC, context)
        raw_path = self.raster_path(parameters, RAW_DEM, context)
        out_path = self.parameterAsFileOutput(parameters, OUTPUT, context)
        if not out_path.lower().endswith(".gpkg"):
            out_path += ".gpkg"

        d8, valid, info = read_dem(fdr_path)

        # CRS rule first - fail before any heavy work. The WKT check decides
        # (geographic refused, metres required); QGIS only supplies a clean
        # EPSG code and name when it recognises the CRS.
        ok, epsg, msg = crs_check(info.projection_wkt)
        dem_crs = QgsCoordinateReferenceSystem()
        dem_crs.createFromWkt(info.projection_wkt)
        if dem_crs.isValid() and dem_crs.isGeographic():
            ok, msg = False, ("The DEM is in a geographic CRS (degrees). Reproject it to "
                              "a projected, metric CRS (e.g. UTM) and re-run.")
        if not ok:
            raise QgsProcessingException(msg)
        if dem_crs.isValid() and dem_crs.authid().upper().startswith("EPSG:"):
            epsg = int(dem_crs.authid().split(":")[1])

        self.check_size(feedback, info.rows, info.cols)
        direction = decode_d8(d8.astype(np.int32))
        accum, _, _ = read_dem(fac_path)
        elevation, raw_valid, raw_info = read_dem(raw_path)
        if (raw_info.rows, raw_info.cols) != (info.rows, info.cols):
            raise QgsProcessingException(
                "The raw DEM and the flow-direction grid have different sizes. "
                "Use the raw DEM the flow direction was derived from.")
        # NoData in the raw DEM must not become an elevation (-9999 would
        # pass as the catchment minimum): mark it NaN, statistics skip NaN.
        elevation = np.where(raw_valid, elevation, np.nan)
        stream_order = None
        if self.parameterAsRasterLayer(parameters, ORDER, context) is not None:
            so, so_valid, _ = read_dem(self.raster_path(parameters, ORDER, context))
            stream_order = np.where(so_valid, so, np.nan)

        id_field = self.field_parameter(parameters, ID_FIELD, context)
        prefix = (self.parameterAsString(parameters, ID_PREFIX, context) or "").strip()
        if id_field:
            prefix = ""   # the prefix applies to sequential IDs only (0.13.1)
        order = ORDER_KEYS[self.parameterAsEnum(parameters, ID_ORDER, context)]
        snap_radius = self.parameterAsInt(parameters, SNAP, context)
        snap_threshold = self.parameterAsDouble(parameters, SNAP_THRESHOLD, context)
        local = self.parameterAsBool(parameters, LOCAL, context)
        sheet_cap = (self.parameterAsDouble(parameters, "SHEET_CAP", context)
                     if "SHEET_CAP" in parameters else 100.0)

        points = self.read_pour_points(parameters, POINTS, context, info, feedback,
                                       id_field=id_field,
                                       extra_fields=[f for f, _ in CANDIDATE_FIELDS])
        extra_layers = []
        points, candidate_mode, all_points = self.candidate_selection(points, feedback)
        if candidate_mode:
            extra_layers.append(("crossing_candidates", "POINT", CANDIDATE_FIELDS,
                                 [((p["x"], p["y"]), {f: p.get("attr_" + f)
                                                      for f, _ in CANDIDATE_FIELDS})
                                  for p in all_points],
                                 "All crossing candidates (audit trail)"))

        burn_summary = ""
        burn_src = self.parameterAsSource(parameters, "BURN_LOG", context) \
            if "BURN_LOG" in parameters else None
        if burn_src is not None:
            from .alg_burn_crossings import LOG_FIELDS as BURN_FIELDS
            from qgis.core import QgsCoordinateReferenceSystem as _CRS, QgsCoordinateTransform as _CT, QgsProject as _P
            dcrs = _CRS(); dcrs.createFromWkt(info.projection_wkt)
            tr = _CT(burn_src.sourceCrs(), dcrs, _P.instance()) \
                if dcrs.isValid() and burn_src.sourceCrs() != dcrs else None
            names = burn_src.fields().names()
            brows = []
            try:                                   # QGIS >= 3.36
                from qgis.core import Qgis as _Q
                _skip = _Q.ProcessingFeatureSourceFlag.SkipGeometryValidityChecks
            except AttributeError:
                from qgis.core import QgsProcessingFeatureSource as _PFS
                _skip = _PFS.Flag.FlagSkipGeometryValidityChecks
            from qgis.core import QgsFeatureRequest as _FR
            # skipped breaches are logged as zero-length lines: keep them (audit record)
            for f in burn_src.getFeatures(_FR(), _skip):
                g = f.geometry()
                if g is None or g.isEmpty():
                    continue
                if tr is not None:
                    g.transform(tr)
                line = g.asMultiPolyline()[0] if g.isMultipart() else g.asPolyline()
                if len(line) < 2:
                    continue
                attrs = {k: (f[k] if k in names else None) for k, _ in BURN_FIELDS}
                attrs = {k: (None if str(v) == "NULL" else v) for k, v in attrs.items()}
                brows.append(([(pt.x(), pt.y()) for pt in line], attrs))
            n_burned = sum(1 for _, a in brows if (a.get("cells") or 0) > 0)
            cut = sum(float(a.get("total_cut_m3") or 0) for _, a in brows)
            burn_summary = (f"burn crossings: {n_burned} of {len(brows)} breached, "
                            f"{cut:,.0f} m3 cut (see layer burn_log)")
            extra_layers.append(("burn_log", "LINESTRING", BURN_FIELDS, brows,
                                 "Breaches cut through embankments before routing"))
            feedback.pushInfo("Conditioning: " + burn_summary)

        alignment = self.read_alignment(
            parameters, ROAD, context, info, feedback,
            start_chainage=self.parameterAsDouble(parameters, START, context),
            reverse=self.parameterAsBool(parameters, REVERSE, context))
        if alignment is not None:
            missing = [p for p in points if p.get("chainage") is None]
            if missing:
                ch, off, _ = alignment.locate([p["x"] for p in missing], [p["y"] for p in missing])
                for p, c in zip(missing, ch):
                    p["chainage"] = float(c)
            extra_layers.append((
                "road_alignment", "LINESTRING",
                [("part", "int"), ("chainage_from_m", "real"), ("chainage_to_m", "real")],
                [([tuple(v) for v in part], {"part": k + 1,
                  "chainage_from_m": float(alignment.ch0[alignment.part == k].min()),
                  "chainage_to_m": float((alignment.ch0 + alignment.seg_len)[alignment.part == k].max())})
                 for k, part in enumerate(alignment.parts)],
                "Road alignment used for chainage"))
            from ..core.network.profile import profile_with_crossings, profile_summary
            from ..core.interop.field_dictionary import OPTIONAL_LAYERS
            filled = None
            if self.parameterAsRasterLayer(parameters, "FILLED", context) is not None:
                fz, fv, fi = read_dem(self.raster_path(parameters, "FILLED", context))
                if (fi.rows, fi.cols) != (info.rows, info.cols):
                    raise QgsProcessingException("The filled DEM is not on the flow-direction grid.")
                filled = np.where(fv, fz, np.nan)
            src = self.parameterAsSource(parameters, ROAD, context)
            align_name = src.sourceName() if src is not None else ""
            profile_step = self.parameterAsDouble(parameters, "PROFILE_STEP", context)
            prof, _ = profile_with_crossings(
                alignment, elevation, info.geotransform, direction=direction, valid=valid,
                accumulation=accum, filled=filled,
                stream_threshold_cells=snap_threshold if snap_threshold > 0 else 200.0,
                stream_order=stream_order, step=profile_step, align_name=align_name)
            pfields = [(f[0], f[1]) for f in OPTIONAL_LAYERS["alignment_profile"][1]]
            extra_layers.append(("alignment_profile", "POINT", pfields,
                                 [((r["x"], r["y"]), r) for r in prof],
                                 OPTIONAL_LAYERS["alignment_profile"][2]))
            feedback.pushInfo("Alignment profile: " + profile_summary(prof))
        if order == "auto":
            order = ("chainage" if all(p.get("chainage") is not None for p in points)
                     else "downstream")
        elif order == "chainage" and not all(p.get("chainage") is not None for p in points):
            raise QgsProcessingException("Chainage numbering needs a candidate layer or a road alignment.")
        stream_mask = None
        if snap_radius > 0 and snap_threshold > 0:
            stream_mask = extract_streams(accum, valid, threshold_cells=snap_threshold)

        feedback.pushInfo(f"{len(points)} pour points; "
                          f"{'local' if local else 'full upstream'} catchments; "
                          f"IDs {'from ' + id_field if id_field else 'sequential ' + (prefix or '') + '001...'}")
        soil = self.load_soil(parameters, context, info, feedback, "SOIL")
        runoff = self.load_runoff(parameters, context, info, soil, feedback)
        rainfall = self.load_rainfall(parameters, context, info, feedback)
        erosion, erosion_run = None, None
        ero_folder = self.parameterAsFile(parameters, "EROSION", context) if "EROSION" in parameters else ""
        if ero_folder:
            from ..core.erosion.io import load_erosion_inputs
            try:
                erosion, erosion_run = load_erosion_inputs(ero_folder, info, read_dem)
            except (ValueError, FileNotFoundError) as e:
                raise QgsProcessingException(str(e))
            feedback.pushInfo(f"Erosion block from {ero_folder} ({erosion.mode}).")
        scenario = self.load_scenario(parameters, context, info, runoff, feedback, erosion)
        try:
            crossings, catchments, flowpaths, issues, id_info = build_exchange_records(
                direction, valid, accum, elevation, info.geotransform, points,
                snap_radius_cells=snap_radius, stream_mask=stream_mask, local=local,
                stream_order=stream_order,
                id_scheme="attribute" if id_field else "sequential",
                id_prefix=prefix, id_order=order, soil=soil, erosion=erosion,
                progress=self.make_progress(feedback, weight=0.9),
                channel_threshold_cells=snap_threshold if snap_threshold > 0 else None,
                sheet_cap_m=sheet_cap, runoff=runoff, rainfall=rainfall, scenario=scenario,
                channel_slope_m=(self.parameterAsDouble(parameters, "CH_SLOPE_DIST", context)
                                 if "CH_SLOPE_DIST" in parameters else 200.0))
        except ExchangeError as e:
            raise QgsProcessingException(str(e))
        for msg in issues:
            feedback.pushWarning(msg)

        # -- drainage coverage check: missing crossings, sags, flat stretches (A2/A3)
        coverage_md, n_proposed = "", 0
        if alignment is not None and self.parameterAsBool(parameters, "COVERAGE", context):
            from ..core.network.coverage import run_coverage, COVERAGE_FIELDS, SAG_FIELDS, FLAT_FIELDS
            from ..core.interop.field_dictionary import OPTIONAL_LAYERS
            import json as _json
            fill_for_wall = filled
            if fill_for_wall is None:
                from ..core.conditioning.fill import fill_depressions
                feedback.pushInfo("Coverage check: no filled DEM given - filling the raw DEM for the "
                                  "sag areas.")
                fz, _, _ = fill_depressions(np.nan_to_num(elevation, nan=0.0), raw_valid & valid,
                                            cell_width=info.cell_width, cell_height=info.cell_height)
                fill_for_wall = np.where(raw_valid, fz, np.nan)
            cp = {k: self.parameterAsDouble(parameters, key, context) for k, key in (
                ("merge_m", "COV_MERGE"), ("search_m", "COV_SEARCH"), ("small_area_km2", "COV_SMALL"),
                ("sag_smooth_m", "SAG_SMOOTH"), ("sag_min_depth_m", "SAG_DEPTH"),
                ("sag_min_area_km2", "SAG_MIN_AREA"), ("flat_slope_pct", "FLAT_SLOPE"),
                ("crossfall_m", "FLAT_CROSSFALL"), ("flat_min_len_m", "FLAT_MIN_LEN"))
                if key in parameters}
            min_a = self.parameterAsDouble(parameters, "COV_MIN_AREA", context) \
                if "COV_MIN_AREA" in parameters else 0.0
            pprefix = (self.parameterAsString(parameters, "PROPOSED_PREFIX", context) or "P").strip() \
                if "PROPOSED_PREFIX" in parameters else "P"
            try:
                cov = run_coverage(
                    direction, valid, accum, elevation, fill_for_wall, info.geotransform, alignment,
                    prof, crossings,
                    build_kwargs=dict(local=False, stream_order=stream_order, soil=soil,
                                      erosion=erosion, runoff=runoff, rainfall=rainfall,
                                      scenario=scenario,
                                      channel_threshold_cells=snap_threshold if snap_threshold > 0 else None,
                                      sheet_cap_m=sheet_cap,
                                      channel_slope_m=(self.parameterAsDouble(parameters, "CH_SLOPE_DIST", context)
                                                       if "CH_SLOPE_DIST" in parameters else 200.0)),
                    min_area_km2=min_a if min_a > 0 else None, proposed_prefix=pprefix,
                    stream_threshold_cells=snap_threshold if snap_threshold > 0 else 200.0, **cp)
            except ExchangeError as e:
                raise QgsProcessingException(str(e))
            for msg in cov["issues"]:
                feedback.pushWarning(msg)
            crossings += cov["crossings"]; catchments += cov["catchments"]; flowpaths += cov["flowpaths"]
            n_proposed = len(cov["crossings"])
            extra_layers.append(("coverage_check", "POINT", COVERAGE_FIELDS, cov["coverage"],
                                 OPTIONAL_LAYERS["coverage_check"][2]))
            extra_layers.append(("sag_points", "POINT", SAG_FIELDS,
                                 [((s_["x"], s_["y"]), s_) for s_ in cov["sags"]],
                                 OPTIONAL_LAYERS["sag_points"][2]))
            flines = []
            for f_ in cov["flats"]:
                cs_ = [f_["chainage_m"]] + [r["chainage_m"] for r in prof
                                            if f_["chainage_m"] < r["chainage_m"] < f_["chainage_to_m"]] \
                    + [f_["chainage_to_m"]]
                flines.append(([alignment.point_at(c_) for c_ in cs_], f_))
            extra_layers.append(("flat_stretches", "LINESTRING", FLAT_FIELDS, flines,
                                 OPTIONAL_LAYERS["flat_stretches"][2]))
            coverage_md = _json.dumps(dict(cp, min_area_km2=min_a or "stream threshold",
                                           proposed_prefix=pprefix), sort_keys=True)
            feedback.pushInfo(
                f"Coverage check: {cov['n_missing']} stream crossing(s) without a culvert, "
                f"{len(cov['sags'])} sag point(s), {len(cov['flats'])} flat stretch(es); "
                f"{n_proposed} proposed crossing(s) delineated ({pprefix}001...).")
            for _, a in cov["coverage"]:
                feedback.pushInfo(f"  ch {a['chainage_m']:,.0f}: {a['issue']} - {a['note']}")

        # -- floodplain width indicator at large crossings (A4) ----------------
        fp_md = ""
        if alignment is not None and prof:
            from ..core.network.floodplain import floodplain_block, hand_grid
            import json as _jfp
            fp_min = self.parameterAsDouble(parameters, "FP_MIN_AREA", context) \
                if "FP_MIN_AREA" in parameters else 10.0
            use_hand = self.parameterAsBool(parameters, "FP_HAND", context) \
                if "FP_HAND" in parameters else True
            big = [a for _, a in crossings if (a.get("acc_at_outlet_km2") or 0) >= fp_min]
            hs = None
            if big and use_hand:
                hmask = extract_streams(accum, valid, threshold_cells=snap_threshold
                                        if snap_threshold > 0 else 200.0)
                hg = hand_grid(direction, valid, elevation, hmask)
                gt_ = info.geotransform
                rr = [min(max(int((r["y"] - gt_[3]) // gt_[5]), 0), info.rows - 1) for r in prof]
                cc = [min(max(int((r["x"] - gt_[0]) // gt_[1]), 0), info.cols - 1) for r in prof]
                hs = hg[rr, cc]
            for (_, a), b in zip(crossings, floodplain_block([a for _, a in crossings], prof,
                                                            min_area_km2=fp_min,
                                                            hand_at_stations=hs)):
                a.update(b)
            fp_md = _jfp.dumps({"min_area_km2": fp_min, "dz_m": [0.5, 1.0, 2.0],
                                "bed_window_m": 50.0, "hand": bool(use_hand),
                                "note": "terrain indicator only; not a flood level"})
            if big:
                feedback.pushInfo(f"Floodplain width indicator at {len(big)} crossing(s) of "
                                  f">= {fp_min:g} km2 (profile{' + HAND' if use_hand else ''}).")
                for _, a in crossings:
                    if a.get("fp_method"):
                        feedback.pushInfo(
                            f"  {a['outlet_uid']}: width at bed + 1 m "
                            f"{a['fp_w_1p0_m'] if a['fp_w_1p0_m'] is None else round(a['fp_w_1p0_m'])} m"
                            + (f" (HAND {round(a['fp_hand_w_1p0_m'])} m)"
                               if a.get("fp_hand_w_1p0_m") is not None else "")
                            + (f" - {a['fp_note']}" if a.get("fp_note") else ""))

        # -- channel cross-section downstream of each crossing (F5) ------------
        xs_md = ""
        if self.parameterAsBool(parameters, "XS", context) if "XS" in parameters else True:
            from ..core.watershed.section import section_block, params_json as _xs_params
            from ..core.interop.field_dictionary import OPTIONAL_LAYERS
            xp = {k: (self.parameterAsDouble(parameters, key, context) if key in parameters else dflt)
                  for k, key, dflt in (("dist_m", "XS_DIST", 30.0), ("half_m", "XS_HALF", 150.0),
                                       ("bank_slope", "XS_BANK_SLOPE", 0.05))}
            blocks, xlines = section_block([a for _, a in crossings], direction, valid, elevation,
                                           info.geotransform, **xp)
            for (_, a), b in zip(crossings, blocks):
                a.update(b)
            xfields = [(f[0], f[1]) for f in OPTIONAL_LAYERS["xs_transects"][1]]
            extra_layers.append(("xs_transects", "LINESTRING", xfields, xlines,
                                 OPTIONAL_LAYERS["xs_transects"][2]))
            xs_md = _xs_params(xp["dist_m"], xp["half_m"], None, xp["bank_slope"])
            nq = {q: sum(1 for b in blocks if b["xs_quality"] == q) for q in ("high", "medium", "low")}
            feedback.pushInfo(f"Channel sections {xp['dist_m']:g} m downstream: {nq['high']} high, "
                              f"{nq['medium']} medium, {nq['low']} low quality (indicative).")
            for _, a in crossings:
                if a.get("xs_quality") == "low":
                    feedback.pushWarning(f"  {a['outlet_uid']}: channel section low quality - "
                                         f"{a.get('xs_note')}")

        # -- catchments cut by an automatic clip (F13) -----------------------------
        from ..core.raster import raster_tags as _rtags
        from ..core.conditioning.autoclip import TAG as _CLIPTAG
        clip_md = _rtags(raw_path).get(_CLIPTAG, "")
        if clip_md:
            from ..core.geometry.rasterize import polygon_window
            n_cut = 0
            for g, a in catchments:
                pw = polygon_window(g, info.geotransform, (info.rows, info.cols))
                if pw is None:
                    continue
                r0, r1, c0, c1, lab = pw
                at_edge = (r0 == 0 and lab[0].any()) or (r1 == info.rows and lab[-1].any()) \
                    or (c0 == 0 and lab[:, 0].any()) or (c1 == info.cols and lab[:, -1].any())
                cut = int(bool(at_edge) or _touches_nodata(lab, raw_valid[r0:r1, c0:c1]))
                a["clip_edge"] = cut
                n_cut += cut
            if n_cut:
                feedback.pushWarning(f"{n_cut} catchment(s) touch the edge of the automatic clip "
                                     "(clip_edge = 1): enlarge the clip margin and run again.")

        # -- check against mapped drainage (F12) -----------------------------------
        mapped_md = ""
        if parameters.get("MAPPED") and self.parameterAsSource(parameters, "MAPPED", context) is not None:
            mapped_md = self._mapped_block(parameters, context, feedback, info, direction, valid,
                                           accum, alignment, crossings, catchments, extra_layers)

        # -- DEM-error sensitivity per crossing (F15; off by default, slow) ------
        unc_md = ""
        if parameters.get("UNC") and self.parameterAsBool(parameters, "UNC", context):
            unc_md = self._uncertainty_block(parameters, context, feedback, info, elevation,
                                             raw_valid, crossings, snap_threshold, snap_radius)

        # -- time of concentration by five methods (F10) -------------------------
        tc_md = ""
        if self.parameterAsBool(parameters, "TC", context) if "TC" in parameters else True:
            tc_md = self._tc_block(parameters, context, feedback, info, crossings, catchments,
                                   flowpaths, runoff)

        # -- side-drain siltation indicator along the corridor (STI R3) ---------
        csti_md = ""
        if alignment is not None and erosion is not None and getattr(erosion, "sti", None) is not None:
            from ..core.erosion import io as _eio
            from ..core.erosion.corridor import sample_corridor
            from ..core.erosion.side_drain import corridor_sti, params_json as _cs_params
            from ..core.interop.field_dictionary import OPTIONAL_LAYERS
            comb, cv, _ = read_dem(_eio.path_of(ero_folder, "combined_class"))
            _, reaches = sample_corridor(alignment, {}, np.where(cv, comb, 0).astype(np.uint8),
                                         info.geotransform, step=profile_step)
            pch = np.array([r_["chainage_m"] for r_ in prof])
            psl = np.array([np.nan if r_["slope_long_pct"] is None else r_["slope_long_pct"]
                            for r_ in prof])

            def _slope_at(ch_):
                return psl[np.abs(pch[None, :] - np.asarray(ch_)[:, None]).argmin(axis=1)] \
                    if pch.size else np.full(len(ch_), np.nan)
            gp_ = lambda k, d: self.parameterAsDouble(parameters, k, context) if k in parameters else d  # noqa: E731
            cs_half, cs_off = gp_("CSTI_HALF", 50.0), gp_("CSTI_OFFSET", 10.0)
            cs_road, cs_slope = gp_("CSTI_ROAD_HALF", 5.0), gp_("CSTI_DRAIN_SLOPE", 1.0)
            sst, srr, sbr = corridor_sti(alignment, erosion.sti, direction, valid, info.geotransform,
                                         reaches, channel=erosion.channel, slope_at=_slope_at,
                                         step=profile_step, half_width=cs_half, offset_step=cs_off,
                                         road_half_m=cs_road, drain_slope_pct=cs_slope)
            lay = OPTIONAL_LAYERS["corridor_sti"]
            crow = []
            for r_, b_ in zip(reaches, srr):
                cs_ = [r_["ch_start"]] + [s_["chainage"] for s_ in sst
                                          if r_["ch_start"] < s_["chainage"] < r_["ch_end"]] + [r_["ch_end"]]
                line_ = [alignment.point_at(c_) for c_ in cs_]
                if len(line_) >= 2:
                    crow.append((line_, dict(r_, **b_)))
            extra_layers.append(("corridor_sti", "LINESTRING", [(f[0], f[1]) for f in lay[1]], crow,
                                 lay[2]))
            csti_md = _cs_params(profile_step, cs_half, cs_off, cs_road, cs_slope, sbr,
                                 "alignment profile on the raw DEM (centred difference)")
            for side, lbl in (("lhs", "left"), ("rhs", "right")):
                flen = sum(b_.get(f"siltation_len_{side}_m") or 0.0 for b_ in srr)
                feedback.pushInfo(f"Side drains, {lbl}: {flen:,.0f} m flagged for siltation "
                                  "(screening; layer corridor_sti).")

        # -- raster quicklooks (A6) ---------------------------------------------
        extra_tables, ql_md = [], ""
        ql_dir = self.parameterAsString(parameters, "QUICKLOOKS", context) \
            if parameters.get("QUICKLOOKS") not in (None, "") else ""
        if ql_dir:
            from ..core.report.quicklooks import make_quicklook, standard_items, INDEX_FIELDS
            from ..core.interop.field_dictionary import OPTIONAL_LAYERS
            import json as _jql
            ql_dir = self.parameterAsFileOutput(parameters, "QUICKLOOKS", context) or ql_dir
            max_px = self.parameterAsInt(parameters, "QL_MAX_PX", context) \
                if "QL_MAX_PX" in parameters else 4096
            base_dir = os.path.dirname(os.path.abspath(out_path))
            qrows = []
            for item in standard_items(raw_path, self.raster_path(parameters, FAC, context),
                                       ero_folder or None):
                try:
                    row = make_quicklook(item, ql_dir, max_px)
                except (ValueError, RuntimeError) as e:
                    feedback.pushWarning(f"Quicklook {item['name']}: {e}")
                    continue
                for k in ("png", "world_file", "legend_json"):
                    row[k] = os.path.relpath(row[k], base_dir).replace(os.sep, "/")
                qrows.append(row)
            extra_tables.append(("rasters", INDEX_FIELDS, qrows, OPTIONAL_LAYERS["rasters"][2]))
            ql_md = _jql.dumps({"folder": ql_dir, "max_px": max_px,
                                "rasters": [r["name"] for r in qrows]})
            feedback.pushInfo(f"Quicklooks: {len(qrows)} raster(s) in {ql_dir}")

        from ..core.raster import raster_tags
        tags_fdr = raster_tags(fdr_path)
        tags_raw = raster_tags(raw_path)
        tags_fill = raster_tags(self.raster_path(parameters, "FILLED", context)) \
            if self.parameterAsRasterLayer(parameters, "FILLED", context) is not None else {}
        cell_area = info.cell_width * info.cell_height
        md = dict(id_info)
        md.update({
            "qeht_version": plugin_version(),
            "cell_size_m": f"{info.cell_width:g} x {info.cell_height:g}",
            "dem_path": fdr_path, "raw_dem_path": raw_path,
            "dem_sha256": file_fingerprint(raw_path),
            "dem_source": (self.parameterAsString(parameters, DEM_SOURCE, context)
                           or (f"{tags_raw['QEHT_DEM_SOURCE']} (from the DEM tags)"
                               if tags_raw.get("QEHT_DEM_SOURCE") else "")
                           or f"file {os.path.basename(raw_path)} (not described by the user)"),
            "dem_prep": tags_raw.get("QEHT_DEM_PREP", ""),
            "conditioning": "; ".join(x for x in (
                self.parameterAsString(parameters, CONDITIONING, context)
                or tags_fill.get("QEHT_CONDITIONING", ""),
                burn_summary) if x) or "not recorded",
            "conditioning_burn": burn_summary,
            "flat_method": (self.parameterAsString(parameters, FLAT_METHOD, context)
                            or tags_fdr.get("QEHT_FLAT_METHOD") or "not recorded"),
            "tie_rule": "QEHT D8: steepest drop/distance; ties by lowest internal index",
            "stream_threshold_cells": f"{snap_threshold:g}",
            "stream_threshold_km2": f"{snap_threshold * cell_area / 1e6:g}",
            "snap_radius_cells": str(snap_radius),
            "sheet_cap_m": f"{sheet_cap:g}",
            "channel_slope_m": f"{(self.parameterAsDouble(parameters, 'CH_SLOPE_DIST', context) if 'CH_SLOPE_DIST' in parameters else 200.0):g}",
            "snap_strategy": ("none" if snap_radius <= 0 else
                              "nearest_stream" if stream_mask is not None else "max_accumulation"),
            "catchment_mode": "local" if local else "full",
            "id_attribute": id_field or "",
            "soil_dataset": soil.info["soil_dataset"] if soil else "",
            "soil_depth_cm": soil.info["soil_depth_cm"] if soil else "",
            "soil_hsg_source": soil.info.get("hsg_source", "") if soil else "",
            "soil_k_source": soil.info.get("k_source", "") if soil else "",
            "runoff_json": runoff.meta_json() if runoff is not None else "",
            "rainfall_json": rainfall.meta_json() if rainfall is not None else "",
            "coverage_params_json": coverage_md,
            "fp_params_json": fp_md,
            "quicklook_params_json": ql_md,
            "xs_params_json": xs_md,
            "corridor_sti_params_json": csti_md,
            "tc_params_json": tc_md,
            "mapped_drainage_json": mapped_md,
            "autoclip_json": clip_md,
            "uncertainty_json": unc_md,
            "scenario_json": scenario.meta_json(
                next((f.get("source", "") for f in (erosion_run or {}).get("factors", [])
                      if f.get("factor") == "C"), "")) if scenario is not None else "",
            "n_proposed": str(n_proposed),
            "crossing_source": ("crossing candidates" if candidate_mode else "pour points"),
            "chainage_start_m": (f"{alignment.start_chainage:g}" if alignment is not None else ""),
            "alignment_source": (self.parameterAsSource(parameters, ROAD, context).sourceName()
                                 if alignment is not None else ""),
            "alignment_step_m": (f"{self.parameterAsDouble(parameters, 'PROFILE_STEP', context):g}"
                                 if alignment is not None else ""),
            "parameters_json": self.parameters_record(parameters, context),
            "erosion_json": ({k: erosion_run.get(k) for k in (
                "mode", "factors", "indices", "schemes", "mcdma_weights", "bulk_density_kgm3",
                "sdr_model", "channel_threshold_cells", "sti", "deposition")} if erosion_run else ""),
        })
        try:
            write_exchange(out_path, crossings, catchments, flowpaths, md,
                           info.projection_wkt, crs_epsg=epsg,
                           crs_name=dem_crs.description() if dem_crs.isValid() else "",
                           csv_dir=(os.path.splitext(out_path)[0] + "_csv")
                           if self.parameterAsBool(parameters, CSV, context) else None,
                           extra_layers=extra_layers, extra_tables=extra_tables)
        except ExchangeError as e:
            raise QgsProcessingException(str(e))

        feedback.pushInfo(
            f"QEHT {plugin_version()} \u00b7 flat method "
            f"{md['flat_method']} "
            f"\u00b7 10-85 reference outlet")
        errors, warnings = validate_exchange(out_path)
        for w in warnings:
            feedback.pushWarning(w)
        if errors:
            raise QgsProcessingException("The package failed validation: " + "; ".join(errors))

        for _, a in catchments:
            feedback.pushInfo(
                f"  {a['outlet_uid']}: A={a['area_km2']:.4f} km2  "
                f"Shorn={a.get('catch_slope_horn', float('nan')):.4f}  "
                f"RR={a.get('catch_relief_ratio', float('nan')):.4f}  "
                f"L={a.get('lfp_length_km', float('nan')):.3f} km")

        stem = os.path.splitext(os.path.basename(out_path))[0]
        for lname in ("catchments", "flowpaths", "crossings"):
            try:
                context.addLayerToLoadOnCompletion(
                    f"{out_path}|layername={lname}",
                    QgsProcessingContext.LayerDetails(f"{stem} - {lname}",
                                                      context.project(), lname))
            except Exception as e:      # loading is a convenience only
                feedback.pushInfo(f"(Could not auto-load {lname}: {e})")

        summary = (f"{len(crossings)} crossings, {len(catchments)} catchments, "
                   f"{len(flowpaths)} flow paths -> {out_path} ({SCHEMA_VERSION})")
        feedback.pushInfo(summary)
        return {OUTPUT: out_path, SUMMARY: summary}

    def _tc_block(self, parameters, context, feedback, info, crossings, catchments, flowpaths,
                  runoff):
        """F10: tc_* fields on every crossing; returns tc_params_json."""
        from ..core.runoff.tc import (tc_block, along_path_value, read_value_lookup, params_json,
                                      SHEET_N, KERBY_N, SHEET_N_ID, KERBY_N_ID, CHANNEL_N)
        gt = info.geotransform

        def cells(xs, ys, grid):
            c = [int((x - gt[0]) // gt[1]) for x in xs]
            r = [int((y - gt[3]) // gt[5]) for y in ys]
            return [grid[r_, c_] if 0 <= r_ < info.rows and 0 <= c_ < info.cols else None
                    for r_, c_ in zip(r, c)]
        sheet, kerby, sid, kid = SHEET_N, KERBY_N, SHEET_N_ID, KERBY_N_ID
        try:
            if parameters.get("TC_SHEET_N_CSV"):
                pth = self.parameterAsFile(parameters, "TC_SHEET_N_CSV", context)
                sheet, sid = read_value_lookup(pth), f"user lookup {pth}"
            if parameters.get("TC_KERBY_N_CSV"):
                pth = self.parameterAsFile(parameters, "TC_KERBY_N_CSV", context)
                kerby, kid = read_value_lookup(pth), f"user lookup {pth}"
        except (OSError, ValueError) as e:
            raise QgsProcessingException(str(e))
        p2, p2_src, p2_grid = None, "", None
        if parameters.get("TC_P2") not in (None, ""):
            p2 = self.parameterAsDouble(parameters, "TC_P2", context)
            p2_src = "user value"
        elif "TC_P2_RASTER" in parameters and \
                self.parameterAsRasterLayer(parameters, "TC_P2_RASTER", context) is not None:
            from ..core.raster import warp_to_grid
            p2_grid, desc = warp_to_grid(self.raster_path(parameters, "TC_P2_RASTER", context), info)
            p2_src = f"raster {desc}, mean along the longest flow path"
        chn = self.parameterAsDouble(parameters, "TC_CHANNEL_N", context) \
            if "TC_CHANNEL_N" in parameters else CHANNEL_N
        cat = {a["outlet_uid"]: a for _, a in catchments}
        paths = {a["outlet_uid"]: (g, a) for g, a in flowpaths}
        n_ok = 0
        for _, xa in crossings:
            uid = xa["outlet_uid"]
            if uid not in paths:
                xa["tc_note"] = "no flow path"
                continue
            line, fa = paths[uid]
            sn = kn = None
            if runoff is not None:
                at = lambda xs, ys: cells(xs, ys, runoff.classes)   # noqa: E731
                sn = along_path_value(line, min(fa.get("lfp_overland_m") or 0.0, 30.0), at, sheet)
                kn = along_path_value(line, fa.get("lfp_overland_m"), at, kerby)
            p2_here = p2
            if p2_grid is not None:
                import math as _m
                v = [x for x in cells([p[0] for p in line], [p[1] for p in line], p2_grid)
                     if x is not None and _m.isfinite(x)]
                p2_here = sum(v) / len(v) if v else None
            b = tc_block(fa, cat.get(uid, {}), xa, sheet_n=sn, kerby_n=kn, p2_mm=p2_here,
                         channel_n=chn)
            xa.update(b)
            n_ok += b["tc_kirpich_min"] is not None
        feedback.pushInfo(f"Time of concentration: {n_ok} crossing(s) with Kirpich; TR-55 "
                          + ("with P2 " + p2_src if p2_src else "empty (give P2)")
                          + ("" if runoff is not None else "; Kerby and TR-55 sheet need land cover")
                          + ". Every method is reported; none is chosen.")
        return params_json(p2, p2_src, chn, sid, kid, "cn_ii")

    def _mapped_block(self, parameters, context, feedback, info, direction, valid, accum,
                      alignment, crossings, catchments, extra_layers):
        """F12: map_* fields, optional layers and coverage findings; returns mapped_drainage_json."""
        from qgis.core import QgsCoordinateReferenceSystem, QgsCoordinateTransform, QgsProject
        from ..core.network import mapped as mp
        from ..core.geometry.rasterize import polygon_window
        from ..core.interop.field_dictionary import OPTIONAL_LAYERS
        src = self.parameterAsSource(parameters, "MAPPED", context)
        name_field = self.field_parameter(parameters, "MAPPED_NAME_FIELD", context) \
            if parameters.get("MAPPED_NAME_FIELD") else ""
        dcrs = QgsCoordinateReferenceSystem()
        dcrs.createFromWkt(info.projection_wkt)
        tr = QgsCoordinateTransform(src.sourceCrs(), dcrs, QgsProject.instance()) \
            if dcrs.isValid() and src.sourceCrs() != dcrs else None
        lines, names = [], []
        for f in src.getFeatures():
            g = f.geometry()
            if g is None or g.isEmpty():
                continue
            if tr is not None:
                g.transform(tr)
            parts = g.asMultiPolyline() if g.isMultipart() else [g.asPolyline()]
            nm = f[name_field] if name_field else None
            nm = None if nm is None or str(nm) == "NULL" else str(nm)
            for ln in parts:
                if len(ln) >= 2:
                    lines.append([(p.x(), p.y()) for p in ln])
                    names.append(nm)
        gp = lambda k, d: self.parameterAsDouble(parameters, k, context) if k in parameters else d  # noqa: E731
        km2, tol, mc = gp("MAPPED_KM2", 1.0), gp("MAPPED_TOL", 60.0), int(gp("MAPPED_MIN_CELLS", 10.0))
        far, search = gp("MAPPED_FAR", 1000.0), gp("COV_SEARCH", 50.0)
        gt = info.geotransform
        cell_area = info.cell_width * info.cell_height
        stream = extract_streams(accum, valid, threshold_cells=max(1.0, km2 * 1e6 / cell_area))
        ctx = mp.prepare(stream, direction, valid, lines, gt, tol)
        overall = mp.scores(ctx)
        shape = (info.rows, info.cols)
        for g, a in catchments:
            pw = polygon_window(g, gt, shape)
            if pw is None:
                continue
            r0, r1, c0, c1, lab = pw
            region = np.zeros(shape, bool)
            region[r0:r1, c0:c1] = lab
            p, r, f1 = mp.scores(ctx, region)
            a.update({"map_precision": p, "map_recall": r, "map_f1": f1})
        for (_, a), b in zip(crossings, mp.crossing_fields([a for _, a in crossings], lines, names,
                                                            tol, km2, far_m=far)):
            a.update(b)
        n_unc = 0
        if alignment is not None:
            found = mp.uncovered_rivers(lines, names, alignment,
                                        [a.get("chainage_m") for _, a in crossings],
                                        search_m=search)
            n_unc = len(found)
            if found:
                cov = [e for e in extra_layers if e[0] == "coverage_check"]
                if cov:
                    cov[0][3].extend(found)
                else:
                    from ..core.network.coverage import COVERAGE_FIELDS
                    extra_layers.append(("coverage_check", "POINT", COVERAGE_FIELDS, found,
                                         OPTIONAL_LAYERS["coverage_check"][2]))
                for (_, a_) in found:
                    feedback.pushWarning(f"  ch {a_['chainage_m']:,.0f}: {a_['note']}")
        div = mp.divergence_reaches(stream, ctx["stream"] & ctx["near_map"], direction, valid,
                                    accum, gt, min_cells=mc)
        for name in ("drainage_divergence", "mapped_rivers_used"):
            lay = OPTIONAL_LAYERS[name]
            rows = div if name == "drainage_divergence" else mp.clip_lines(lines, names, gt, shape)
            extra_layers.append((name, lay[0], [(f[0], f[1]) for f in lay[1]], rows, lay[2]))
        fmt = lambda v: "-" if v is None else f"{v:.2f}"  # noqa: E731
        feedback.pushInfo(f"Mapped drainage ({len(lines)} line(s), streams >= {km2:g} km2, "
                          f"tolerance {tol:g} m): precision {fmt(overall[0])}, recall "
                          f"{fmt(overall[1])}, F1 {fmt(overall[2])}; "
                          f"{sum(1 for _, a in crossings if a.get('map_agrees') == 0)} crossing(s) "
                          f"off the map, {n_unc} mapped river(s) without a crossing, {len(div)} "
                          "divergence reach(es).")
        return mp.params_json((self.parameterAsString(parameters, "MAPPED_SOURCE", context) or
                               src.sourceName()), km2, tol, mc, overall, len(lines), name_field,
                              far_m=far, road_search_m=search)


    def _uncertainty_block(self, parameters, context, feedback, info, elevation, raw_valid,
                           crossings, snap_threshold, snap_radius):
        """F15 in the package: unc_* on every crossing; returns uncertainty_json."""
        from ..core.watershed import uncertainty as un
        sigma, preset, src, corr, n, seed, sens = self.uncertainty_settings(parameters, context)
        gt = info.geotransform
        outlets, idx = [], []
        for i, (_, a) in enumerate(crossings):
            c = int((a["outlet_x"] - gt[0]) // gt[1])
            r = int((a["outlet_y"] - gt[3]) // gt[5])
            if 0 <= r < info.rows and 0 <= c < info.cols:
                outlets.append((r, c))
                idx.append(i)
        thr = snap_threshold if snap_threshold > 0 else 200.0
        snap = max(int(snap_radius), 1)
        dem = np.where(raw_valid, elevation, 0.0)
        if sens:
            feedback.pushInfo("The correlation-length sensitivity adds 2 x N realisations "
                              "(about three times the run time); switch it off to save time.")
        res, inf = un.run(dem, raw_valid, gt, outlets, sigma, corr, n, seed, thr, snap,
                          progress=lambda f: feedback.setProgress(100.0 * f),
                          log=feedback.pushInfo)
        for i, b in zip(idx, res):
            crossings[i][1].update(b)
        sensitivity, uids = None, None
        if sens and outlets:
            big = sorted(range(len(idx)), key=lambda k: -(crossings[idx[k]][1].get("acc_at_outlet_km2") or 0))[:3]
            uids = [crossings[idx[k]][1]["outlet_uid"] for k in big]
            feedback.pushInfo(f"Correlation-length sensitivity (0.5x, 2x) for {', '.join(uids)}")
            sensitivity = un.corr_sensitivity(dem, raw_valid, gt, [outlets[k] for k in big], sigma,
                                              corr, n, seed, thr, snap)
        n_sw = sum(1 for _, a in crossings if (a.get("unc_switch_pct") or 0) > 0)
        feedback.pushInfo(f"DEM uncertainty: sigma {sigma:g} m ({preset}), L {corr:g} m, N {n}, "
                          f"seed {seed}; {n_sw} crossing(s) with catchment switching.")
        return un.params_json(inf, preset, src, sensitivity, uids)


def _touches_nodata(mask, valid):
    """True when a mask cell is 8-adjacent to an invalid cell inside the window."""
    inv = np.pad(~np.asarray(valid, bool), 1, constant_values=False)
    m = np.asarray(mask, bool)
    rows, cols = m.shape
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if (m & inv[1 + dr:1 + dr + rows, 1 + dc:1 + dc + cols]).any():
                return True
    return False
