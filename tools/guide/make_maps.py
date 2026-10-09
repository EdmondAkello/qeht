"""Map figures for the user guide, drawn by QGIS print layouts (offscreen)."""
import os, json
os.environ["QT_QPA_PLATFORM"] = "offscreen"
from qgis.core import (QgsApplication, QgsProject, QgsRasterLayer, QgsVectorLayer, QgsPrintLayout,
                       QgsLayoutItemMap, QgsLayoutItemScaleBar, QgsLayoutItemLegend, QgsLayoutItemLabel,
                       QgsLayoutPoint, QgsLayoutSize, QgsUnitTypes, QgsLayoutExporter, QgsRectangle,
                       QgsCoordinateReferenceSystem, QgsSingleSymbolRenderer, QgsCategorizedSymbolRenderer,
                       QgsRendererCategory, QgsGraduatedSymbolRenderer, QgsRendererRange, QgsPalLayerSettings,
                       QgsTextFormat, QgsTextBufferSettings, QgsVectorLayerSimpleLabeling, QgsMarkerSymbol,
                       QgsLineSymbol, QgsFillSymbol, QgsLegendStyle, QgsRuleBasedRenderer,
                       QgsPalettedRasterRenderer, QgsFeatureRequest)
from qgis.PyQt.QtGui import QColor, QFont
from qgis.PyQt.QtCore import QRectF
app = QgsApplication([], True); app.initQgis()

_REPO_PARENT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
G = os.path.abspath(os.environ.get("QEHT_GUIDE_BUILD", "guide_build"))  # holds demo/, run/, fig/
R = f"{G}/run"; D = f"{G}/demo"; FIG = f"{G}/fig"
PKG = f"{R}/package/design_hydrology.gpkg"
CRS = QgsCoordinateReferenceSystem("EPSG:32737")
P = QgsProject.instance(); P.setCrs(CRS)
FONT = "Inter"


def raster(path, name, opacity=1.0):
    l = QgsRasterLayer(path, name); l.setCrs(CRS); l.renderer().setOpacity(opacity)
    P.addMapLayer(l, False); return l


def vector(path, name, layer=None, subset=None):
    src = f"{path}|layername={layer}" if layer else path
    l = QgsVectorLayer(src, name, "ogr")
    if not l.isValid():
        raise RuntimeError(f"cannot open {src}")
    if subset: l.setSubsetString(subset)
    P.addMapLayer(l, False); return l


def label(l, field, size=7.0, color="#1d2330", expr=False, dist=1.6):
    s = QgsPalLayerSettings(); s.fieldName = field; s.isExpression = expr
    f = QgsTextFormat(); f.setFont(QFont(FONT)); f.setSize(size); f.setColor(QColor(color))
    b = QgsTextBufferSettings(); b.setEnabled(True); b.setSize(0.9); b.setColor(QColor("white")); f.setBuffer(b)
    s.setFormat(f); s.dist = dist
    l.setLabeling(QgsVectorLayerSimpleLabeling(s)); l.setLabelsEnabled(True)


def marker(shape, color, size, outline="#ffffff", ow=0.5):
    return QgsMarkerSymbol.createSimple({"name": shape, "color": color, "size": str(size),
                                         "outline_color": outline, "outline_width": str(ow)})


def line(color, width, style="solid"):
    return QgsLineSymbol.createSimple({"color": color, "width": str(width), "line_style": style,
                                       "capstyle": "round", "joinstyle": "round"})


def single(l, sym):
    l.setRenderer(QgsSingleSymbolRenderer(sym))


def road_style(l, w=1.5):
    sym = line("#ffffff", w + 1.0)
    top = line("#1d2330", w)
    sym.appendSymbolLayer(top.symbolLayer(0).clone()); single(l, sym)


def streams_style(l, color="#2f7fbf", scale=1.0):
    ranges = []
    for o, w in ((1, 0.22), (2, 0.38), (3, 0.6), (4, 0.9), (5, 1.25), (6, 1.6)):
        ranges.append(QgsRendererRange(o - 0.5, o + 0.5, line(color, w * scale), f"Strahler {o}"))
    l.setRenderer(QgsGraduatedSymbolRenderer('"order"', ranges))


def crossings_style(l, size=3.0, labels=True, lsize=6.6):
    cats = [QgsRendererCategory("existing", marker("circle", "#1f6f8b", size), "Existing crossing"),
            QgsRendererCategory("proposed", marker("triangle", "#d9822b", size + 0.7), "Proposed crossing")]
    l.setRenderer(QgsCategorizedSymbolRenderer("status", cats))
    if labels: label(l, "outlet_uid", lsize)


def poly(l, color="#24324f", width=0.3, fill="0,0,0,0", style="solid"):
    single(l, QgsFillSymbol.createSimple({"color": fill, "outline_color": color,
                                          "outline_width": str(width), "outline_style": style}))


FULL = None


def clamp(extent):
    if FULL is None:
        return extent
    dx = dy = 0.0
    if extent.xMaximum() > FULL.xMaximum(): dx = FULL.xMaximum() - extent.xMaximum()
    if extent.xMinimum() + dx < FULL.xMinimum(): dx = FULL.xMinimum() - extent.xMinimum()
    if extent.yMaximum() > FULL.yMaximum(): dy = FULL.yMaximum() - extent.yMaximum()
    if extent.yMinimum() + dy < FULL.yMinimum(): dy = FULL.yMinimum() - extent.yMinimum()
    return QgsRectangle(extent.xMinimum() + dx, extent.yMinimum() + dy, extent.xMaximum() + dx, extent.yMaximum() + dy)


def layout(name, layers, extent, w_mm=180.0, h_mm=None, legend=None, scalebar=True, north=True,
           dpi=220, legend_pos="below", legend_cols=3):
    if h_mm is None:
        h_mm = w_mm * extent.height() / extent.width()
    else:  # widen the extent to the frame
        want = extent.width() / extent.height(); have = w_mm / h_mm
        c = extent.center()
        if have > want:
            hw = extent.height() * have / 2; extent = QgsRectangle(c.x() - hw, extent.yMinimum(), c.x() + hw, extent.yMaximum())
        else:
            hh = extent.width() / have / 2; extent = QgsRectangle(extent.xMinimum(), c.y() - hh, extent.xMaximum(), c.y() + hh)
    extent = clamp(extent)
    below = bool(legend) and legend_pos == "below"
    band = 0.0
    if below:
        rows = sum(len(l.renderer().legendSymbolItems()) + (1 if len(l.renderer().legendSymbolItems()) > 1 else 0)
                   for l in legend)
        band = 7.0 + 5.4 * -(-rows // legend_cols)
    lay = QgsPrintLayout(P); lay.initializeDefaults()
    lay.pageCollection().page(0).setPageSize(QgsLayoutSize(w_mm, h_mm + band, QgsUnitTypes.LayoutMillimeters))
    m = QgsLayoutItemMap(lay); m.attemptSetSceneRect(QRectF(0, 0, w_mm, h_mm))
    m.setFrameEnabled(False); m.setCrs(CRS); m.setLayers(layers); m.setExtent(extent)
    m.setBackgroundColor(QColor("#eef1f5")); lay.addLayoutItem(m)
    tf7 = QgsTextFormat(); tf7.setFont(QFont(FONT)); tf7.setSize(7)
    if scalebar:
        sb = QgsLayoutItemScaleBar(lay); sb.setStyle("Single Box"); sb.setLinkedMap(m)
        sb.applyDefaultSize(); sb.setUnits(QgsUnitTypes.DistanceKilometers); sb.setUnitLabel("km")
        target = extent.width() / 1000.0 / 9.0
        nice = min([0.1, 0.2, 0.25, 0.5, 1, 2, 2.5, 5], key=lambda v: abs(v - target))
        sb.setUnitsPerSegment(nice); sb.setNumberOfSegments(2); sb.setNumberOfSegmentsLeft(0)
        sb.setHeight(1.5); sb.setTextFormat(tf7)
        sb.setBackgroundEnabled(True); sb.setBackgroundColor(QColor(255, 255, 255, 220))
        sb.setBoxContentSpace(1.2)
        lay.addLayoutItem(sb); sb.refresh(); sb.resizeToMinimumWidth()
        sb.attemptMove(QgsLayoutPoint(3, h_mm - sb.sizeWithUnits().height() - 3, QgsUnitTypes.LayoutMillimeters))
    if north:
        lab = QgsLayoutItemLabel(lay); lab.setText("▲ N")
        tf = QgsTextFormat(); tf.setFont(QFont(FONT, 8, QFont.Bold)); tf.setSize(8); lab.setTextFormat(tf)
        lab.setBackgroundEnabled(True); lab.setBackgroundColor(QColor(255, 255, 255, 220))
        lab.setMarginX(1.2); lab.setMarginY(0.8); lab.adjustSizeToText(); lay.addLayoutItem(lab)
        lab.attemptMove(QgsLayoutPoint(w_mm - lab.sizeWithUnits().width() - 3, 3, QgsUnitTypes.LayoutMillimeters))
    if legend:
        lg = QgsLayoutItemLegend(lay); lg.setLinkedMap(m); lg.setAutoUpdateModel(False)
        root = lg.model().rootGroup(); root.removeAllChildren()
        for l in legend:
            root.addLayer(l)
        lg.setTitle("")
        for st in (QgsLegendStyle.Subgroup, QgsLegendStyle.SymbolLabel, QgsLegendStyle.Group):
            lg.rstyle(st).setTextFormat(tf7)
        lg.setSymbolWidth(5); lg.setSymbolHeight(2.8)
        lg.setColumnCount(legend_cols if below else 1); lg.setSplitLayer(True); lg.setEqualColumnWidth(True)
        lg.setBackgroundEnabled(True); lg.setBackgroundColor(QColor(255, 255, 255, 255 if below else 230))
        lg.setFrameEnabled(not below); lg.setFrameStrokeColor(QColor("#c9d0db"))
        from qgis.core import QgsLayoutItem
        ref = {"tl": QgsLayoutItem.UpperLeft, "tr": QgsLayoutItem.UpperRight, "below": QgsLayoutItem.UpperLeft,
               "bl": QgsLayoutItem.LowerLeft, "br": QgsLayoutItem.LowerRight}[legend_pos]
        lay.addLayoutItem(lg); lg.setResizeToContents(True); lg.refresh(); lg.adjustBoxSize()
        lg.setReferencePoint(ref)
        x = w_mm - 3 if legend_pos.endswith("r") else (0 if below else 3)
        y = h_mm + 1 if below else (h_mm - 3 if legend_pos.startswith("b") else (11 if legend_pos == "tr" and north else 3))
        if below:
            lg.setResizeToContents(False)
            lg.attemptResize(QgsLayoutSize(w_mm, band - 1, QgsUnitTypes.LayoutMillimeters))
        lg.attemptMove(QgsLayoutPoint(x, y, QgsUnitTypes.LayoutMillimeters))
    st = QgsLayoutExporter.ImageExportSettings(); st.dpi = dpi
    print("map", name, QgsLayoutExporter(lay).exportToImage(f"{FIG}/{name}.png", st))


# chainage -> extent from the alignment profile
prof = QgsVectorLayer(f"{R}/layers/alignment_profile.gpkg", "p", "ogr")
PTS = sorted((f["chainage_m"], f.geometry().asPoint()) for f in prof.getFeatures())


def around(ch0, ch1, pad):
    xs = [p.x() for c, p in PTS if ch0 <= c <= ch1]; ys = [p.y() for c, p in PTS if ch0 <= c <= ch1]
    return QgsRectangle(min(xs) - pad, min(ys) - pad, max(xs) + pad, max(ys) + pad)


def feature_extent(layer, expr, pad):
    e = QgsRectangle(); e.setMinimal()
    for f in layer.getFeatures(QgsFeatureRequest().setFilterExpression(expr)):
        e.combineExtentWith(f.geometry().boundingBox())
    return QgsRectangle(e.xMinimum() - pad, e.yMinimum() - pad, e.xMaximum() + pad, e.yMaximum() + pad)


# ---------------------------------------------------------------- layers
relief = raster(f"{R}/quicklooks/relief.png", "Relief")
relief_soft = raster(f"{R}/quicklooks/relief.png", "Relief", 0.45)
streams = vector(f"{R}/layers/streams.gpkg", "Streams"); streams_style(streams)
road = vector(f"{R}/layers/road_alignment.gpkg", "Road centreline"); road_style(road)
road_thin = vector(f"{R}/layers/road_alignment.gpkg", "Road centreline"); road_style(road_thin, 1.0)
cr = vector(f"{R}/layers/crossings.gpkg", "Crossings"); crossings_style(cr)
cr_small = vector(f"{R}/layers/crossings.gpkg", "Crossings"); crossings_style(cr_small, 2.3, lsize=5.8)
cr_nolab = vector(f"{R}/layers/crossings.gpkg", "Crossings"); crossings_style(cr_nolab, 2.6, labels=False)
ca = vector(f"{R}/layers/catchments.gpkg", "Catchments"); poly(ca, "#24324f", 0.28)
fp = vector(f"{R}/layers/flowpaths.gpkg", "Longest flow paths"); single(fp, line("#c0392b", 0.4, "dash"))
kp = vector(f"{R}/layers/alignment_profile.gpkg", "Kilometre posts", subset='"chainage_m" % 2000 = 0')
single(kp, marker("square", "#ffffff", 1.4, "#1d2330", 0.35))
label(kp, "'km ' || to_string(\"chainage_m\" / 1000)", 6.0, "#1d2330", expr=True, dist=1.2)
full = relief.extent()
FULL = full
st_leg = QgsVectorLayer("LineString?crs=EPSG:32737", "Streams (width by Strahler order)", "memory")
single(st_leg, line("#2f7fbf", 0.6)); P.addMapLayer(st_leg, False)
st_leg_g = QgsVectorLayer("LineString?crs=EPSG:32737", "DEM streams (width by Strahler order)", "memory")
single(st_leg_g, line("#5b6575", 0.5)); P.addMapLayer(st_leg_g, False)

# 1. overview and cover
layout("map_overview", [kp, cr_small, road, fp, ca, streams, relief], full, 180,
       legend=[cr, road, st_leg, ca, fp, kp])
layout("cover_map", [cr_nolab, road, streams, relief], full, 210, legend=None, scalebar=False,
       north=False, dpi=200)

# 2. parallel valley: candidates and clusters, then the crossings and coverage findings
cand = vector(f"{R}/layers/crossing_candidates.gpkg", "Crossing candidates")
rr = QgsRuleBasedRenderer(marker("diamond", "#3a9e5f", 3.0)); root = rr.rootRule(); root.removeChildAt(0)
root.appendChild(QgsRuleBasedRenderer.Rule(marker("diamond", "#3a9e5f", 3.4), 0, 0, '"recommended" = 1', "Recommended candidate"))
root.appendChild(QgsRuleBasedRenderer.Rule(marker("diamond", "#ffffff", 2.6, "#3a9e5f", 0.6), 0, 0, '"recommended" = 0', "Other candidate in the cluster"))
cand.setRenderer(rr)
label(cand, "\"cand_id\" || ' (' || \"cluster_id\" || ')'", 6.0, "#24553a", expr=True, dist=2.0)
par = vector(PKG, "Parallel reach", "crossing_candidates", subset="0")
try:
    par = vector(f"{R}/layers/parallel_reaches.gpkg", "Parallel reach (side drain)")
    single(par, line("#7a4fb3", 2.2)); par.renderer().symbol().setOpacity(0.7)
except RuntimeError:
    pass
ex_par = around(5300, 10900, 450)
layout("map_candidates", [cand, road, par, streams, relief_soft], ex_par, 180, h_mm=120,
       legend=[cand, par, road, st_leg])

cov = vector(f"{R}/layers/coverage_check.gpkg", "Coverage findings")
cats = [QgsRendererCategory("missing_crossing", marker("cross2", "#c0392b", 3.4, "#c0392b", 0.7), "Missing crossing"),
        QgsRendererCategory("sag_point", marker("circle", "#ffffff", 2.0, "#7a4fb3", 0.7), "Sag point"),
        QgsRendererCategory("flat_stretch", marker("square", "#e0b400", 2.0, "#8a6d00", 0.4), "Flat stretch (start)"),
        QgsRendererCategory("small_area", marker("circle", "#9aa7bd", 1.8), "Small area"),
        QgsRendererCategory("mapped_river_uncovered", marker("cross2", "#1f6f8b", 3.2, "#1f6f8b", 0.7), "Mapped river, no crossing")]
cov.setRenderer(QgsCategorizedSymbolRenderer("issue", cats))
flats = vector(f"{R}/layers/flat_stretches.gpkg", "Flat stretch")
single(flats, line("#e0b400", 3.4)); flats.renderer().symbol().setOpacity(0.8)
cov_pts = vector(f"{R}/layers/coverage_check.gpkg", "Coverage findings", subset="\"issue\" <> 'flat_stretch'")
cov_pts.setRenderer(cov.renderer().clone())
layout("map_coverage", [cr, cov_pts, road, flats, streams, relief_soft], ex_par, 180, h_mm=120,
       legend=[cr, cov_pts, road, st_leg])

# 3. the lake bed: flat stretches, sags, crossings
ex_lake = around(15600, 22000, 600)
layout("map_lakebed", [cr, cov_pts, road, flats, streams, ca, relief_soft], ex_lake, 180, h_mm=105,
       legend=[cr, cov_pts, flats, ca, st_leg])

# 4. one catchment: flow path, segments, 10-85 points, channel section
uid = "X014"
ca1 = vector(f"{R}/layers/catchments.gpkg", f"Catchment {uid}", subset=f"\"outlet_uid\" = '{uid}'")
poly(ca1, "#1f6f8b", 0.6, "31,111,139,40")
fp1 = vector(f"{R}/layers/flowpaths.gpkg", "Longest flow path", subset=f"\"outlet_uid\" = '{uid}'")
single(fp1, line("#c0392b", 0.8))
xs1 = vector(PKG, "Channel section", "xs_transects", subset=f"\"outlet_uid\" = '{uid}'")
single(xs1, line("#d9822b", 1.0))
cr1 = vector(f"{R}/layers/crossings.gpkg", "Crossing", subset=f"\"outlet_uid\" = '{uid}'"); crossings_style(cr1, 3.4)
pts = vector(f"{G}/meta/lfp_points.gpkg", "Flow path") if os.path.exists(f"{G}/meta/lfp_points.gpkg") else None
lyrs = [cr1] + ([pts] if pts else []) + [xs1, fp1, road, streams, ca1, relief_soft]
if pts:
    rr = QgsRuleBasedRenderer(marker("circle", "#ffffff", 2.4)); root = rr.rootRule(); root.removeChildAt(0)
    root.appendChild(QgsRuleBasedRenderer.Rule(marker("circle", "#ffffff", 2.4, "#c0392b", 0.8), 0, 0, "\"kind\" IN ('10', '85')", "10 % and 85 % points of the path"))
    root.appendChild(QgsRuleBasedRenderer.Rule(marker("star", "#e0b400", 3.6, "#5b4500", 0.4), 0, 0, "\"kind\" = 'head'", "Channel head"))
    pts.setRenderer(rr); label(pts, "lab", 6.4, "#7a1f16")
e1 = feature_extent(ca1, "1=1", 500)
layout("map_catchment", lyrs, e1, 180, h_mm=125,
       legend=[cr1, ca1, fp1] + ([pts] if pts else []) + [xs1, st_leg])

# 5. mapped drainage: mapped lines, DEM streams, divergence, agreement at crossings
mapped = vector(PKG, "Mapped waterways", "mapped_rivers_used"); single(mapped, line("#1f6f8b", 1.6))
mapped.renderer().symbol().setOpacity(0.55)
div = vector(PKG, "Divergence reach", "drainage_divergence"); single(div, line("#c0392b", 1.2))
agree = vector(f"{R}/layers/crossings.gpkg", "Crossing on the mapped waterway", subset='"status" IS NOT NULL')
P.addMapLayer(agree, False)
# map_agrees lives in the package crossings
agree = vector(PKG, "Crossings", "crossings")
rr = QgsRuleBasedRenderer(marker("circle", "#3a9e5f", 3.0)); root = rr.rootRule(); root.removeChildAt(0)
root.appendChild(QgsRuleBasedRenderer.Rule(marker("circle", "#3a9e5f", 3.2), 0, 0, '"map_agrees" = 1', "On the mapped waterway"))
root.appendChild(QgsRuleBasedRenderer.Rule(marker("circle", "#c0392b", 3.2), 0, 0, '"map_agrees" = 0', "Off the mapped waterway"))
root.appendChild(QgsRuleBasedRenderer.Rule(marker("circle", "#ffffff", 2.2, "#5b6575", 0.5), 0, 0, '"map_agrees" IS NULL', "Not compared (small or far)"))
agree.setRenderer(rr); label(agree, "outlet_uid", 6.0)
streams_dem = vector(f"{R}/layers/streams.gpkg", "DEM streams"); streams_style(streams_dem, "#5b6575", 0.8)
ex_map = QgsRectangle(full.xMinimum() + 4500, full.yMinimum() + 1000, full.xMaximum(), full.yMinimum() + 13500)
layout("map_mapped", [agree, road_thin, div, mapped, streams_dem, relief_soft], ex_map, 180, h_mm=110,
       legend=[agree, mapped, st_leg_g, div])

# 6. erosion: combined class quicklook with the road and side-drain siltation
comb = raster(f"{R}/quicklooks/combined_class.png", "Combined erosion class")
sti = vector(PKG, "Side-drain siltation", "corridor_sti")
layout("map_erosion", [cr_small, road, comb, relief_soft], around(1200, 9000, 900), 180, h_mm=110,
       legend=None)

# 7. land cover and scenario
for nm, path in (("map_landcover", f"{D}/demo_landcover.tif"), ("map_landcover_2040", f"{D}/demo_landcover_2040.tif")):
    lc = raster(path, "Land cover (WorldCover classes)")
    classes = [(10, "#006400", "Tree cover"), (20, "#ffbb22", "Shrubland"), (30, "#ffff4c", "Grassland"),
               (40, "#f096ff", "Cropland"), (50, "#fa0000", "Built-up"), (60, "#b4b4b4", "Bare / sparse")]
    cl = [QgsPalettedRasterRenderer.Class(v, QColor(c), t) for v, c, t in classes]
    rd = QgsPalettedRasterRenderer(lc.dataProvider(), 1, cl); rd.setOpacity(0.85); lc.setRenderer(rd)
    layout(nm, [cr_nolab, road_thin, ca, lc], full, 88, legend=[lc] if nm == "map_landcover" else None,
           legend_pos="bl", scalebar=nm == "map_landcover", dpi=240)
print("done")
