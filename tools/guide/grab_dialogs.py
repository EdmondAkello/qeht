"""Screenshots of QEHT Processing dialogs with the demo layers loaded (offscreen Qt)."""
import sys, os
os.environ["QT_QPA_PLATFORM"] = "offscreen"
_REPO_PARENT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
G = os.path.abspath(os.environ.get("QEHT_GUIDE_BUILD", "guide_build"))  # holds demo/, run/, fig/
sys.path.insert(0, _REPO_PARENT); sys.path.append('/usr/share/qgis/python/plugins')
from qgis.core import QgsApplication, QgsProject, QgsRasterLayer, QgsVectorLayer
app = QgsApplication([], True); app.initQgis()
from qgis.gui import QgsGui
QgsGui.editorWidgetRegistry().initEditors()
from processing.core.Processing import Processing
Processing.initialize()
from qeht.processing_provider.provider import QehtProvider
prov = QehtProvider(); QgsApplication.processingRegistry().addProvider(prov)
from processing.gui.AlgorithmDialog import AlgorithmDialog
from qgis.PyQt.QtGui import QFont
f = app.font(); f.setPointSizeF(9.5); app.setFont(f)

D = f"{G}/demo"
FIG = f"{G}/fig"
L = {}
for name, path in [("demo_dem", "demo_dem.tif"), ("demo_landcover", "demo_landcover.tif"),
                   ("demo_landcover_2040", "demo_landcover_2040.tif"), ("demo_map", "demo_map.tif"),
                   ("demo_r", "demo_r.tif")]:
    L[name] = QgsRasterLayer(f"{D}/{path}", name)
for name, path in [("demo_road", "demo_road.gpkg"), ("demo_soils", "demo_soils.gpkg"),
                   ("demo_rain_zones", "demo_rain_zones.gpkg"), ("demo_rivers", "demo_rivers.gpkg")]:
    L[name] = QgsVectorLayer(f"{D}/{path}", name, "ogr")
R = f"{G}/run"
for name, path in [("filled", "rasters/filled.tif"), ("flow_direction", "rasters/flow_direction.tif"),
                   ("flow_accumulation", "rasters/flow_accumulation.tif"), ("strahler", "rasters/strahler.tif"),
                   ("streams", "rasters/streams.tif")]:
    L[name] = QgsRasterLayer(f"{R}/{path}", name)
for name, path in [("crossing_candidates", "layers/crossing_candidates.gpkg"), ("crossings", "layers/crossings.gpkg")]:
    L[name] = QgsVectorLayer(f"{R}/{path}", name, "ogr")
for lyr in L.values():
    QgsProject.instance().addMapLayer(lyr)

def grab(alg_id, out, w, h, params=None, scroll=0, crop=None):
    alg = QgsApplication.processingRegistry().createAlgorithmById(alg_id)
    dlg = AlgorithmDialog(alg)
    if params:
        dlg.setParameters({k: (L[v].id() if isinstance(v, str) and v in L else v) for k, v in params.items()})
    dlg.setFixedSize(w, h); dlg.show(); app.processEvents()
    if scroll:
        from qgis.PyQt.QtWidgets import QScrollArea
        for sa in dlg.findChildren(QScrollArea):
            sa.verticalScrollBar().setValue(scroll)
        app.processEvents()
    img = dlg.grab()
    if crop:
        img = img.copy(0, 0, crop, img.height())
    img.save(f"{FIG}/{out}.png")
    print("saved", out)
    dlg.close()

PIPE = {"DEM": "demo_dem", "ROAD": "demo_road", "STREAM_KM2": 0.5, "FLAT_METHOD": 1,
        "NODATA_OVERRIDE": 0, "DEM_SOURCE": "synthetic demonstration DEM, 30 m", "RUN_NAME": "Demo Road",
        "SOIL_POLYGONS": "demo_soils", "LANDCOVER": "demo_landcover", "LANDCOVER_SCN": "demo_landcover_2040",
        "SCENARIO_NAME": "2040 build-out", "RAIN_ZONES": "demo_rain_zones", "RAIN_ZONE_FIELD": "zone",
        "RAIN_MAP": "demo_map", "MAPPED": "demo_rivers", "MAPPED_NAME_FIELD": "name", "MAPPED_KM2": 8.0,
        "TC_P2": 60.0, "R_RASTER": "demo_r", "OUTPUT": "C:/Projects/Demo Road/qeht_run"}
grab("qeht:hydrologypipeline", "dlg_pipeline_1", 1000, 1500, PIPE, crop=720)
grab("qeht:hydrologypipeline", "dlg_pipeline_2", 1000, 1500, PIPE, scroll=1180, crop=720)
grab("qeht:preparedem", "dlg_prepare", 1000, 900, {"DEMS": ["demo_dem"], "SOURCE": "synthetic demonstration DEM"})
grab("qeht:filldepressions", "dlg_fill", 900, 560, {"INPUT": "demo_dem", "NODATA_OVERRIDE": 0})
grab("qeht:flowdirection", "dlg_d8", 900, 520, {"DEM": "filled"})
grab("qeht:crossingcandidates", "dlg_candidates", 1000, 1000, {"FDR": "flow_direction", "FAC": "flow_accumulation", "STREAMS": "streams", "ORDER": "strahler", "ROAD": "demo_road", "THRESHOLD": 556})
grab("qeht:buildheasexchange", "dlg_package", 1000, 1050, {"FDR": "flow_direction", "FAC": "flow_accumulation", "RAW_DEM": "demo_dem", "ORDER": "strahler", "POINTS": "crossing_candidates", "ROAD": "demo_road", "FILLED": "filled", "DEM_SOURCE": "synthetic demonstration DEM, 30 m"})
grab("qeht:soilparameters", "dlg_soils", 1000, 1000, {"SOIL_POLYGONS": "demo_soils", "LANDCOVER": "demo_landcover", "RAIN_ZONES": "demo_rain_zones", "RAIN_ZONE_FIELD": "zone", "RAIN_MAP": "demo_map"})
grab("qeht:erosionindices", "dlg_erosion", 1000, 1000, {"RAW_DEM": "demo_dem", "FAC": "flow_accumulation", "R_RASTER": "demo_r", "SOIL_POLYGONS": "demo_soils", "WORLDCOVER": "demo_landcover"})
grab("qeht:demuncertainty", "dlg_uncertainty", 1000, 950, {"DEM": "demo_dem", "CROSSINGS": "crossings", "THRESHOLD": 556, "UNC_PRESET": 1})
grab("qeht:exportpackage", "dlg_export", 900, 560)
grab("qeht:autoclip", "dlg_autoclip", 900, 700, {"DEM": "demo_dem", "ROAD": "demo_road"})
grab("qeht:drainagecoverage", "dlg_coverage", 1000, 1000, {"RAW_DEM": "demo_dem", "FILLED": "filled", "FDR": "flow_direction", "FAC": "flow_accumulation", "ROAD": "demo_road", "CROSSINGS": "crossings", "THRESHOLD": 556})
