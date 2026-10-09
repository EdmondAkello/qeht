import os, sys, json
_REPO_PARENT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
G = os.path.abspath(os.environ.get("QEHT_GUIDE_BUILD", "guide_build"))  # holds demo/, run/, fig/
sys.path.insert(0, _REPO_PARENT)
from qgis.core import QgsApplication
app = QgsApplication([], False); app.initQgis()
sys.path.append('/usr/share/qgis/python/plugins')
import processing
from processing.core.Processing import Processing
Processing.initialize()
from qeht.processing_provider.provider import QehtProvider
_prov = QehtProvider(); print("added", QgsApplication.processingRegistry().addProvider(_prov), len(_prov.algorithms()))
D = f"{G}/demo"
OUT = f"{G}/run"
p = {"DEM": f"{D}/demo_dem.tif", "ROAD": f"{D}/demo_road.gpkg", "START": 0.0, "STREAM_KM2": 0.5,
     "FLAT_METHOD": 1, "NODATA_OVERRIDE": 0, "DEM_SOURCE": "synthetic demonstration DEM, 30 m",
     "RUN_NAME": "Demo Road", "MAPPED": f"{D}/demo_rivers.gpkg", "MAPPED_NAME_FIELD": "name", "MAPPED_KM2": 8.0, "TC_P2": 60.0,
     "SOIL_POLYGONS": f"{D}/demo_soils.gpkg", "LANDCOVER": f"{D}/demo_landcover.tif",
     "LANDCOVER_SCN": f"{D}/demo_landcover_2040.tif",
     "RAIN_ZONES": f"{D}/demo_rain_zones.gpkg", "RAIN_ZONE_FIELD": "zone", "RAIN_MAP": f"{D}/demo_map.tif",
     "RAIN_MAP_DATASET": "synthetic climatology", "R_RASTER": f"{D}/demo_r.tif",
     "EROSION": True, "UNC": True, "UNC_PRESET": 1, "UNC_N": 30, "SCENARIO_NAME": "2040 build-out", "COVERAGE": True, "QUICKLOOKS": True, "EXPORTS": True, "PACKAGE": True,
     "LOAD": False, "OUTPUT": OUT}
if len(sys.argv) > 1 and sys.argv[1] == "review":
    p.update({"STAGE": 1, "OUTPUT": OUT + "_review"})
r = processing.run("qeht:hydrologypipeline", p)
print(json.dumps({k: str(v) for k, v in r.items()}, indent=1))
