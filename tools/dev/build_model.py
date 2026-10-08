"""Build examples/qeht_corridor.model3 (F16) with the QGIS model API, so the example stays
in step with the parameter names. Run from the parent of the repo folder:

    QT_QPA_PLATFORM=offscreen python3.12 qeht/tools/dev/build_model.py

Fill -> D8 -> accumulation -> streams -> crossing candidates -> Build design hydrology
package; exposed inputs DEM, road and stream threshold (cells).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..")))

from qgis.core import QgsApplication  # noqa: E402

app = QgsApplication([], False)
app.initQgis()
from qgis.PyQt.QtCore import QPointF  # noqa: E402
from qgis.core import (QgsProcessingModelAlgorithm, QgsProcessingModelParameter,  # noqa: E402
                       QgsProcessingModelChildAlgorithm, QgsProcessingModelChildParameterSource as Src,
                       QgsProcessingModelOutput, QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterFeatureSource, QgsProcessingParameterNumber,
                       QgsProcessing)

OUT = os.path.abspath(os.path.join(HERE, "..", "..", "examples", "qeht_corridor.model3"))


def build():
    m = QgsProcessingModelAlgorithm("QEHT road corridor (example)", "QEHT examples")
    m.setHelpContent({"ALG_DESC": "Fill, D8 flow direction, flow accumulation, streams, road "
                      "crossing candidates and the design hydrology package for a road corridor. "
                      "An example of chaining QEHT tools in the Graphical Modeler; the one-click "
                      "pipeline does more."})
    for i, (p, x) in enumerate((
            (QgsProcessingParameterRasterLayer("dem", "DEM (raw, projected, metres)"), 120),
            (QgsProcessingParameterFeatureSource("road", "Road alignment",
                                                 [QgsProcessing.SourceType.TypeVectorLine]), 420),
            (QgsProcessingParameterNumber("threshold", "Stream threshold (cells)",
                                          QgsProcessingParameterNumber.Type.Double, 200.0,
                                          minValue=1.0), 720))):
        mp = QgsProcessingModelParameter(p.name())
        mp.setPosition(QPointF(x, 60))
        m.addModelParameter(p, mp)

    def child(cid, alg, desc, y, sources, outputs=None, x=300):
        c = QgsProcessingModelChildAlgorithm(alg)
        c.setChildId(cid)
        c.setDescription(desc)
        c.setPosition(QPointF(x, y))
        for k, v in sources.items():
            c.addParameterSources(k, [v])
        if outputs:
            c.setModelOutputs(outputs)          # before adding: the model stores a copy
        m.addChildAlgorithm(c)

    P = Src.fromModelParameter
    V = Src.fromStaticValue
    O = Src.fromChildOutput
    child("fill", "qeht:filldepressions", "Fill depressions", 160,
          {"DEM": P("dem"), "MIN_SLOPE": V(0.0)})
    child("fdr", "qeht:flowdirection", "D8 flow direction (Barnes)", 260,
          {"DEM": O("fill", "OUTPUT"), "RESOLVE_FLATS": V(True), "FLAT_METHOD": V(1)})
    child("fac", "qeht:flowaccumulation", "Flow accumulation", 360,
          {"FDR": O("fdr", "OUTPUT"), "QUANTITY": V(0)})
    child("streams", "qeht:streamnetwork", "Stream network + Strahler", 460,
          {"FDR": O("fdr", "OUTPUT"), "FAC": O("fac", "OUTPUT"), "MODE": V(0),
           "THRESHOLD": P("threshold")})
    child("cand", "qeht:crossingcandidates", "Road crossing candidates", 560,
          {"FDR": O("fdr", "OUTPUT"), "FAC": O("fac", "OUTPUT"), "STREAMS": O("streams", "STREAMS"),
           "THRESHOLD": P("threshold"), "ORDER": O("streams", "ORDER"), "ROAD": P("road"),
           "START": V(0.0), "REVERSE": V(False)})
    out = QgsProcessingModelOutput("Design hydrology package", "Design hydrology package")
    out.setChildId("package")
    out.setChildOutputName("OUTPUT")
    child("package", "qeht:buildheasexchange", "Build design hydrology package", 660,
          {"FDR": O("fdr", "OUTPUT"), "FAC": O("fac", "OUTPUT"), "RAW_DEM": P("dem"),
           "ORDER": O("streams", "ORDER"), "POINTS": O("cand", "CANDIDATES"),
           "ROAD": P("road"), "SNAP_THRESHOLD": P("threshold"), "FILLED": O("fill", "OUTPUT")},
          outputs={"Design hydrology package": out})
    m.updateDestinationParameters()
    return m


if __name__ == "__main__":
    from qeht.processing_provider.provider import QehtProvider
    provider = QehtProvider()                       # keep a reference: the registry does not own it
    QgsApplication.processingRegistry().addProvider(provider)
    model = build()
    ok, errs = model.validate()
    if not ok:
        print("model does not validate:", errs)
        sys.exit(1)
    model.toFile(OUT)
    print(OUT, [p.name() for p in model.parameterDefinitions()])
    app.exitQgis()
