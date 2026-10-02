# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Road x drainage crossing candidates (WP-A)."""

from qgis.core import (
    QgsProcessingParameterRasterLayer, QgsProcessingParameterFeatureSource,
    QgsProcessingParameterNumber, QgsProcessingParameterBoolean,
    QgsProcessingParameterString, QgsProcessingParameterVectorDestination,
    QgsProcessing, QgsProcessingException,
)
import numpy as np

from .base import QehtAlgorithm
from ..core.raster import read_dem
from ..core.grid import decode_d8
from ..core.watershed.delineate import extract_streams
from ..core.flow.streamlines import vectorize_streams
from ..core.network.crossings import (find_crossing_candidates, CANDIDATE_FIELDS,
                                      PARALLEL_FIELDS)

FDR = "FDR"; FAC = "FAC"; STREAMS = "STREAMS"; THRESHOLD = "THRESHOLD"; ORDER = "ORDER"
ROAD = "ROAD"; START = "START"; REVERSE = "REVERSE"; MIN_AREA = "MIN_AREA"
HALFWIDTH = "HALFWIDTH"; MIN_PARALLEL = "MIN_PARALLEL"; MERGE = "MERGE"; PREFIX = "PREFIX"
CANDIDATES = "CANDIDATES"; PARALLEL = "PARALLEL"


class _StatusValueMap(object):
    """Post-processor: make 'status' a drop-down (candidate/accepted/rejected)."""

    def __init__(self):
        from qgis.core import QgsProcessingLayerPostProcessorInterface

        class _P(QgsProcessingLayerPostProcessorInterface):
            def postProcessLayer(self, layer, context, feedback):
                try:
                    from qgis.core import QgsEditorWidgetSetup
                    idx = layer.fields().indexOf("status")
                    if idx >= 0:
                        layer.setEditorWidgetSetup(idx, QgsEditorWidgetSetup(
                            "ValueMap", {"map": {"candidate": "candidate",
                                                 "accepted": "accepted",
                                                 "rejected": "rejected"}}))
                except Exception as e:          # cosmetic only
                    feedback.pushInfo(f"(status drop-down not set: {e})")
        self.instance = _P()


class CrossingCandidatesAlgorithm(QehtAlgorithm):

    def name(self): return "crossingcandidates"
    def displayName(self): return "Road crossing candidates"
    def group(self): return "Road drainage"
    def groupId(self): return "roaddrainage"

    def shortHelpString(self):
        return (
            "Finds every place the drainage network crosses a road alignment and "
            "proposes one crossing per drainage line, without deleting anything.\n\n"
            "<b>Method.</b> Each D8 flow link of the stream network near the road "
            "is intersected with the centreline - the same segments QEHT's stream "
            "polylines are made of - so a crossing cannot be missed through a "
            "digitising gap. Every intersection is a candidate with chainage, "
            "crossing angle (90 = square), contributing area, Strahler order, "
            "reach id and the side the flow comes from (L/R of the chainage "
            "direction).\n\n"
            "<b>Parallel flow.</b> A stream that runs alongside the road within "
            "the corridor half-width for at least the minimum parallel length "
            "crosses the centreline many times on a DEM. Those candidates share a "
            "cluster, and the stream is written to the parallel-reaches layer - "
            "that is where a side drain has to carry the water to the crossing. "
            "Candidates closer than the merge distance also share a cluster.\n\n"
            "<b>Recommendation, not deletion.</b> In each cluster the most "
            "downstream candidate (largest contributing area) gets recommended = 1. "
            "All candidates have status = candidate: set accepted / rejected, then "
            "run 'Build design hydrology package' on this layer. It uses the "
            "accepted candidates, or the recommended ones if none is accepted, "
            "numbers crossings along the chainage and takes each outlet cell "
            "exactly as found here (no snapping).\n\n"
            "Chainage runs along the road layer's feature order and digitising "
            "direction from the start value; tick 'Reverse' to run it the other way.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(FDR, "Flow direction (D8-coded)"))
        self.addParameter(QgsProcessingParameterRasterLayer(FAC, "Flow accumulation"))
        self.addParameter(QgsProcessingParameterRasterLayer(
            STREAMS, "Stream raster (optional; blank = accumulation threshold)", optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            THRESHOLD, "Stream threshold (cells) when no stream raster is given",
            QgsProcessingParameterNumber.Type.Double, defaultValue=200.0, minValue=1.0))
        self.addParameter(QgsProcessingParameterRasterLayer(
            ORDER, "Stream order raster (optional)", optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource(
            ROAD, "Road alignment (centreline)", [QgsProcessing.SourceType.TypeVectorLine]))
        self.addParameter(QgsProcessingParameterNumber(
            START, "Start chainage (m)", QgsProcessingParameterNumber.Type.Double, defaultValue=0.0))
        self.addParameter(QgsProcessingParameterBoolean(
            REVERSE, "Reverse chainage direction", defaultValue=False))
        self.addParameter(QgsProcessingParameterNumber(
            MIN_AREA, "Minimum contributing area (km2)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=0.01, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber(
            HALFWIDTH, "Corridor half-width for parallel flow (m)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=30.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber(
            MIN_PARALLEL, "Minimum parallel length (m)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=100.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber(
            MERGE, "Cluster merge distance along chainage (m; 0 = off)",
            QgsProcessingParameterNumber.Type.Double, defaultValue=0.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterString(
            PREFIX, "Candidate ID prefix", defaultValue="C", optional=True))
        self.addParameter(QgsProcessingParameterVectorDestination(
            CANDIDATES, "Crossing candidates", QgsProcessing.SourceType.TypeVectorPoint))
        self.addParameter(QgsProcessingParameterVectorDestination(
            PARALLEL, "Parallel reaches (side-drain hints)", QgsProcessing.SourceType.TypeVectorLine))

    def processAlgorithm(self, parameters, context, feedback):
        d8, valid, info = read_dem(self.raster_path(parameters, FDR, context))
        self.check_size(feedback, info.rows, info.cols)
        direction = decode_d8(d8.astype(np.int32))
        accum, _, _ = read_dem(self.raster_path(parameters, FAC, context))
        if self.parameterAsRasterLayer(parameters, STREAMS, context) is not None:
            st, sv, _ = read_dem(self.raster_path(parameters, STREAMS, context))
            stream_mask = sv & (st > 0)
            feedback.pushInfo(f"Stream network from raster: {int(stream_mask.sum()):,} cells")
        else:
            thr = self.parameterAsDouble(parameters, THRESHOLD, context)
            stream_mask = extract_streams(accum, valid, threshold_cells=thr)
            feedback.pushInfo(f"Stream network at >= {thr:,.0f} cells: "
                              f"{int(stream_mask.sum()):,} cells")
        order = None
        if self.parameterAsRasterLayer(parameters, ORDER, context) is not None:
            so, sov, _ = read_dem(self.raster_path(parameters, ORDER, context))
            order = np.where(sov, so, np.nan)
        feedback.setProgressText("Labelling stream reaches")
        reaches = vectorize_streams(direction, valid, stream_mask,
                                    info.cell_width, info.cell_height)
        reach_ids = np.zeros(direction.shape, dtype=np.float64)
        for rch in reaches:
            for r, c in rch["cells"][:-1] or rch["cells"]:
                if reach_ids[r, c] == 0:
                    reach_ids[r, c] = rch["reach_id"]
        feedback.setProgress(40)

        alignment = self.read_alignment(
            parameters, ROAD, context, info, feedback,
            start_chainage=self.parameterAsDouble(parameters, START, context),
            reverse=self.parameterAsBool(parameters, REVERSE, context))
        cand, par, summary = find_crossing_candidates(
            direction, valid, accum, stream_mask, info.geotransform, alignment,
            min_area_km2=self.parameterAsDouble(parameters, MIN_AREA, context),
            corridor_halfwidth_m=self.parameterAsDouble(parameters, HALFWIDTH, context),
            min_parallel_m=self.parameterAsDouble(parameters, MIN_PARALLEL, context),
            merge_distance_m=self.parameterAsDouble(parameters, MERGE, context),
            stream_order=order, reach_ids=reach_ids,
            prefix=(self.parameterAsString(parameters, PREFIX, context) or "").strip())
        feedback.setProgress(85)

        cand_out = self.parameterAsOutputLayer(parameters, CANDIDATES, context)
        par_out = self.parameterAsOutputLayer(parameters, PARALLEL, context)
        self.write_vector(cand_out, "crossing_candidates", info.projection_wkt, "point",
                          CANDIDATE_FIELDS, cand)
        self.write_vector(par_out, "parallel_reaches", info.projection_wkt, "line",
                          PARALLEL_FIELDS, par)
        for _, a in cand:
            feedback.pushInfo(
                f"  {a['cand_id']}  ch {a['chainage_m']:,.1f}  A={a['acc_km2']:.4f} km2  "
                f"angle {a['crossing_angle_deg']:.0f}  cluster {a['cluster_id']}"
                f"{'  RECOMMENDED' if a['recommended'] else ''}"
                f"{'  parallel' if a['parallel_reach'] else ''}")
        self.report_stats(feedback, "Crossing candidates", summary)
        if context.willLoadLayerOnCompletion(cand_out):
            self._post = _StatusValueMap()
            context.layerToLoadOnCompletionDetails(cand_out).setPostProcessor(self._post.instance)
        return {CANDIDATES: cand_out, PARALLEL: par_out}
