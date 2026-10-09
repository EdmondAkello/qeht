"""Charts for the user guide (matplotlib), from the demonstration run."""
import json, sqlite3, csv, sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from osgeo import gdal, ogr
gdal.UseExceptions(); ogr.UseExceptions()
import os
_REPO_PARENT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
G = os.path.abspath(os.environ.get("QEHT_GUIDE_BUILD", "guide_build"))  # holds demo/, run/, fig/
sys.path.insert(0, _REPO_PARENT)

R = f"{G}/run"; D = f"{G}/demo"; FIG = f"{G}/fig"
PKG = f"{R}/package/design_hydrology.gpkg"
INK, MUTED, LINE, NAVY, ACC, ORANGE, RED, GREEN = "#1d2330", "#5b6575", "#d9dee7", "#1f2a44", "#1f6f8b", "#d9822b", "#c0392b", "#3a9e5f"
plt.rcParams.update({"font.family": "Inter", "font.size": 8.5, "axes.edgecolor": "#9aa7bd", "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlesize": 9.5, "axes.titleweight": "bold",
                     "axes.titlecolor": NAVY, "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "figure.dpi": 100, "savefig.dpi": 220})
db = sqlite3.connect(PKG); db.row_factory = sqlite3.Row
CR = [dict(r) for r in db.execute("select * from crossings order by chainage_m")]


def save(fig, name):
    fig.savefig(f"{FIG}/{name}.png", bbox_inches="tight", pad_inches=0.06, facecolor="white"); plt.close(fig)
    print("chart", name)


# 1 alignment ground profile -------------------------------------------------------------
prof = [dict(r) for r in db.execute("select * from alignment_profile order by chainage_m")]
ch = np.array([p["chainage_m"] for p in prof]) / 1000
z = np.array([p["z_dem_m"] if p["z_dem_m"] is not None else np.nan for p in prof])
zf = np.array([p["z_fill_m"] if p["z_fill_m"] is not None else np.nan for p in prof])
acc = np.array([p["acc_km2"] or np.nan for p in prof])
fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.2, 4.2), sharex=True, gridspec_kw={"height_ratios": [2.3, 1]})
for f in db.execute("select chainage_m, chainage_to_m from flat_stretches"):
    for a in (a1, a2):
        a.axvspan(f[0] / 1000, f[1] / 1000, color="#f3d65b", alpha=0.35, lw=0)
a1.fill_between(ch, z, zf, where=zf > z + 0.05, color=ACC, alpha=0.35, lw=0, label="Ponding removed by the fill")
a1.plot(ch, z, color=INK, lw=0.9, label="Ground (raw DEM, bilinear)")
for c in CR:
    if c["chainage_m"] is None: continue
    zc = np.interp(c["chainage_m"] / 1000, ch, z)
    ex = c["status"] == "existing"
    a1.plot(c["chainage_m"] / 1000, zc, "o" if ex else "^", ms=4.2 if ex else 4.8, color=ACC if ex else ORANGE,
            mec="white", mew=0.6, zorder=5)
sag = [dict(r) for r in db.execute("select chainage_m, z_dem_m from sag_points")]
a1.plot([s["chainage_m"] / 1000 for s in sag], [s["z_dem_m"] for s in sag], "o", ms=3.0, mfc="white", mec="#7a4fb3", mew=0.8, zorder=4)
a1.set_ylabel("Elevation (m)")
h = [Line2D([], [], color=INK, lw=0.9, label="Ground (raw DEM)"), Patch(color=ACC, alpha=0.35, label="Ponding the fill removed"),
     Line2D([], [], marker="o", ls="", color=ACC, mec="white", label="Existing crossing"),
     Line2D([], [], marker="^", ls="", color=ORANGE, mec="white", label="Proposed crossing"),
     Line2D([], [], marker="o", ls="", mfc="white", mec="#7a4fb3", label="Sag point"),
     Patch(color="#f3d65b", alpha=0.5, label="Flat stretch")]
a1.legend(handles=h, ncol=3, loc="upper right", fontsize=7.4)
a2.semilogy(ch, acc, color=ACC, lw=0.8)
st = np.array([p["stream"] == 1 for p in prof])
a2.semilogy(ch[st], acc[st], "|", color=NAVY, ms=6)
a2.set_ylabel("Area (km²)"); a2.set_xlabel("Chainage (km)")
a2.set_ylim(0.0009, 200)
save(fig, "chart_profile")

# 2 channel sections ---------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.6), sharey=False)
for a, uid in zip(axs, ("P012", "X014")):
    c = next(x for x in CR if x["outlet_uid"] == uid)
    se = np.array(json.loads(c["xs_station_elev_json"]))
    a.plot(se[:, 0], se[:, 1], color=INK, lw=1.0)
    zb = c["xs_bed_m"]
    a.axhline(zb, color=MUTED, lw=0.6, ls=":")
    if c["xs_bankfull_d_m"]:
        a.fill_between(se[:, 0], se[:, 1], zb + c["xs_bankfull_d_m"], where=se[:, 1] < zb + c["xs_bankfull_d_m"], color=ACC, alpha=0.3, lw=0)
        a.axhline(zb + c["xs_bankfull_d_m"], color=ACC, lw=0.8)
    for dz, col in ((1.0, ORANGE),):
        a.axhline(zb + dz, color=col, lw=0.7, ls="--")
    a.set_title(f"{uid}: {c['acc_at_outlet_km2']:.1f} km², quality {c['xs_quality']}", loc="left")
    bw = c["xs_bankfull_w_m"]; bd = c["xs_bankfull_d_m"]
    a.text(0.02, 0.95, f"bank-full {bw:.0f} m wide, {bd:.1f} m deep\nwidth at bed + 1 m: {c['xs_w_1p0_m']:.0f} m",
           transform=a.transAxes, va="top", fontsize=7.3, color=INK)
    a.set_xlabel("Offset, right looking downstream (m)")
axs[0].set_ylabel("Elevation (m)")
fig.legend(handles=[Line2D([], [], color=INK, label="Raw DEM, bilinear"), Patch(color=ACC, alpha=0.3, label="Below bank-full level"),
                    Line2D([], [], color=ACC, label="Bank-full level"), Line2D([], [], color=ORANGE, ls="--", label="Bed + 1 m"),
                    Line2D([], [], color=MUTED, ls=":", label="Bed")], ncol=5, loc="lower center", fontsize=7, bbox_to_anchor=(0.5, -0.1))
save(fig, "chart_sections")

# 3 time of concentration ----------------------------------------------------------------
M = [("tc_kirpich", "Kirpich"), ("tc_kerby_kirpich", "Kerby + Kirpich"), ("tc_scs_lag", "SCS lag"),
     ("tc_tr55", "TR-55"), ("tc_bransby_williams", "Bransby-Williams")]
ex = [c for c in CR if c["status"] == "existing"]
fig, a = plt.subplots(figsize=(7.2, 3.0))
cols = [ACC, GREEN, ORANGE, "#7a4fb3", RED]
xs = np.arange(len(ex))
for k, ((f, lab), col) in enumerate(zip(M, cols)):
    for i, c in enumerate(ex):
        v = c[f + "_min"]
        if v is None: continue
        within = (c[f + "_flag"] or "").startswith("within") or (c[f + "_flag"] or "") == "rural catchments"
        a.plot(i + (k - 2) * 0.13, v, "o", ms=4, color=col if within else "white", mec=col, mew=1.0)
a.set_yscale("log"); a.set_ylabel("Tc (min)")
a.set_xticks(xs); a.set_xticklabels([c["outlet_uid"] for c in ex], rotation=90, fontsize=7)
h = [Line2D([], [], marker="o", ls="", color=col, label=lab) for (f, lab), col in zip(M, cols)] + \
    [Line2D([], [], marker="o", ls="", mfc="white", mec=MUTED, label="outside the method's range")]
a.legend(handles=h, ncol=6, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, 1.13))
a.grid(axis="y", color=LINE, lw=0.5)
save(fig, "chart_tc")

# 4 DEM uncertainty -----------------------------------------------------------------------
fig, a = plt.subplots(figsize=(7.2, 2.9))
for i, c in enumerate(ex):
    det = c["acc_at_outlet_km2"]
    if c["unc_area_p10"] is None:
        a.plot(i, det, "x", color=MUTED); continue
    sw = (c["unc_switch_pct"] or 0) >= 25
    a.plot([i, i], [c["unc_area_p10"], c["unc_area_p90"]], color=RED if sw else ACC, lw=3.2, solid_capstyle="butt", alpha=0.75)
    a.plot(i, c["unc_area_p50"], "_", color=INK, ms=9, mew=1.2)
    a.plot(i, det, "o", ms=3.6, color="white", mec=INK, mew=0.9, zorder=5)
a.set_yscale("log"); a.set_ylabel("Catchment area (km²)")
a.set_xticks(range(len(ex))); a.set_xticklabels([c["outlet_uid"] for c in ex], rotation=90, fontsize=7)
a.grid(axis="y", color=LINE, lw=0.5)
a.legend(handles=[Line2D([], [], color=ACC, lw=3.2, alpha=0.75, label="P10–P90"),
                  Line2D([], [], color=RED, lw=3.2, alpha=0.75, label="P10–P90, area moves > 25 % in ≥ 25 % of runs"),
                  Line2D([], [], marker="_", ls="", color=INK, ms=9, mew=1.2, label="P50"),
                  Line2D([], [], marker="o", ls="", mfc="white", mec=INK, label="Deterministic area")],
         ncol=4, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, 1.14))
save(fig, "chart_uncertainty")

# 5 flat methods on the lake bed -----------------------------------------------------------
from qeht.core.conditioning.fill import fill_depressions
from qeht.core.flow.direction import d8_direction
from qeht.core.flow.accumulation import flow_accumulation
ds = gdal.Open(f"{D}/demo_dem.tif"); zz = ds.ReadAsArray().astype(float); gt = ds.GetGeoTransform()
r0, r1, c0, c1 = 440, 700, 380, 714
sub = zz[r0:r1, c0:c1]; valid = sub > 0
fz, _, _ = fill_depressions(np.where(valid, sub, 0.0), valid, cell_width=30.0, cell_height=30.0)
ls = matplotlib.colors.LightSource(315, 45)
shade = ls.hillshade(np.where(valid, sub, np.nan), vert_exag=4, dx=30, dy=30)
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.95))
ext = [gt[0] + c0 * 30, gt[0] + c1 * 30, gt[3] - r1 * 30, gt[3] - r0 * 30]
nflat = None
streams = {}
for a, (meth, lab) in zip(axs, (("barnes", "Barnes 2014 (default)"), ("toward", "Toward lower terrain"))):
    d, info = d8_direction(fz, valid, 30.0, 30.0, flat_method=meth)
    acc_, _ = flow_accumulation(d, valid)
    s = acc_ >= 556
    streams[meth] = s
    a.imshow(shade, cmap="gray", extent=ext, vmin=0.2, vmax=1.0)
    rgba = np.zeros(s.shape + (4,)); rgba[s] = matplotlib.colors.to_rgba(ACC if meth == "barnes" else RED)
    a.imshow(rgba, extent=ext, interpolation="nearest")
    a.set_title(lab, loc="left"); a.set_xticks([]); a.set_yticks([])
    for sp in a.spines.values(): sp.set_visible(True); sp.set_color(LINE)
    rd = ogr.Open(f"{R}/layers/road_alignment.gpkg"); rl = rd.GetLayer()
    for fe in rl:
        p = np.array(fe.GetGeometryRef().GetPoints()); a.plot(p[:, 0], p[:, 1], color=INK, lw=1.6)
    a.set_xlim(ext[0], ext[1]); a.set_ylim(ext[2], ext[3])
both = streams["barnes"] & streams["toward"]; anyv = streams["barnes"] | streams["toward"]
iou = both.sum() / anyv.sum()
flat_share = float((np.isclose(fz, np.roll(fz, 1, 0)) & valid).mean())
fig.text(0.01, -0.02, f"Streams ≥ 0.5 km² on the same filled DEM. Cells in a stream under both methods: {100 * iou:.0f} % of the cells in a stream under either.", fontsize=7.3, color=MUTED)
save(fig, "chart_flats")
print("stream IoU", iou)

# 6 land cover baseline and scenario --------------------------------------------------------
classes = [(10, "#006400", "Tree cover"), (20, "#ffbb22", "Shrubland"), (30, "#ffff4c", "Grassland"),
           (40, "#f096ff", "Cropland"), (50, "#fa0000", "Built-up"), (60, "#b4b4b4", "Bare / sparse")]
cmap = matplotlib.colors.ListedColormap(["#ffffff"] + [c for _, c, _ in classes])
norm = matplotlib.colors.BoundaryNorm([-5, 5, 15, 25, 35, 45, 55, 65], cmap.N)
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.6))
for a, (fn, tit) in zip(axs, (("demo_landcover.tif", "Baseline"), ("demo_landcover_2040.tif", "Scenario: 2040 build-out"))):
    d = gdal.Open(f"{D}/{fn}"); lc = d.ReadAsArray(); g2 = d.GetGeoTransform()
    e2 = [g2[0], g2[0] + g2[1] * lc.shape[1], g2[3] + g2[5] * lc.shape[0], g2[3]]
    a.imshow(lc, cmap=cmap, norm=norm, extent=e2, interpolation="nearest")
    rd = ogr.Open(f"{R}/layers/road_alignment.gpkg"); rl = rd.GetLayer()
    for fe in rl:
        p = np.array(fe.GetGeometryRef().GetPoints()); a.plot(p[:, 0], p[:, 1], color=INK, lw=1.4)
    a.set_title(tit, loc="left"); a.set_xticks([]); a.set_yticks([])
    a.set_xlim(e2[0] + 1000, e2[1] - 200); a.set_ylim(e2[2] + 3500, e2[3] - 4500)
fig.subplots_adjust(bottom=0.1, wspace=0.05)
fig.legend(handles=[Patch(color=c, label=t) for _, c, t in classes], ncol=6, loc="lower center", fontsize=7.2, bbox_to_anchor=(0.5, 0.0))
save(fig, "chart_landcover")

# 7 erosion quicklooks ------------------------------------------------------------------------
from PIL import Image
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.2))
for a, nm in zip(axs, ("spi_class", "rusle_class")):
    im = np.array(Image.open(f"{R}/quicklooks/{nm}.png"))
    leg = json.load(open(f"{R}/quicklooks/{nm}_legend.json"))
    pgw = [float(x) for x in open(f"{R}/quicklooks/{nm}.pgw").read().split()]
    e3 = [pgw[4] - pgw[0] / 2, pgw[4] + pgw[0] * (im.shape[1] - 0.5), pgw[5] + pgw[3] * (im.shape[0] - 0.5), pgw[5] - pgw[3] / 2]
    a.imshow(shade_full := ls.hillshade(np.where(zz > 0, zz, np.nan), vert_exag=3, dx=30, dy=30), cmap="gray", extent=[gt[0], gt[0] + 30 * zz.shape[1], gt[3] - 30 * zz.shape[0], gt[3]], vmin=0, vmax=1.2)
    a.imshow(im, extent=e3, alpha=0.8, interpolation="nearest")
    rd = ogr.Open(f"{R}/layers/road_alignment.gpkg"); rl = rd.GetLayer()
    for fe in rl:
        p = np.array(fe.GetGeometryRef().GetPoints()); a.plot(p[:, 0], p[:, 1], color=INK, lw=1.4)
    a.set_xlim(gt[0] + 500, gt[0] + 10500); a.set_ylim(gt[3] - 14500, gt[3] - 7500)
    a.set_title(leg["title"], loc="left"); a.set_xticks([]); a.set_yticks([])
    a.legend(handles=[Patch(color=c["colour"], label=c["label"]) for c in leg["classes"]], fontsize=6.6, loc="lower left",
             frameon=True, facecolor="white", framealpha=0.9, edgecolor=LINE)
save(fig, "chart_erosion")

# 8 side-drain siltation along the corridor -------------------------------------------------
cs_ = [dict(r) for r in db.execute("select * from corridor_sti order by ch_start")]
fig, a = plt.subplots(figsize=(7.2, 2.2))
for side, col, sgn in (("lhs", ACC, 1), ("rhs", ORANGE, -1)):
    for r in cs_:
        v = r[f"sti_p90_{side}"]
        if v is None: continue
        a.bar((r["ch_start"] + r["ch_end"]) / 2000, sgn * v, width=(r["ch_end"] - r["ch_start"]) / 1000, color=col, alpha=0.75, lw=0)
        if (r[f"siltation_len_{side}_m"] or 0) > 0:
            a.plot((r["ch_start"] + r["ch_end"]) / 2000, sgn * (v + 2), "v" if sgn < 0 else "^", color=RED, ms=3.5)
a.axhline(0, color=INK, lw=0.6)
a.set_xlabel("Chainage (km)"); a.set_ylabel("STI p90  (LHS up, RHS down)", fontsize=7.5)
a.legend(handles=[Patch(color=ACC, alpha=0.75, label="Left-hand side"), Patch(color=ORANGE, alpha=0.75, label="Right-hand side"),
                  Line2D([], [], marker="^", ls="", color=RED, label="Siltation flag (high class, slope < 1 %)")], ncol=3, fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 1.0))
yl = max(abs(x) for x in a.get_ylim()); a.set_ylim(-yl, yl)
save(fig, "chart_siltation")

# 9 DEM checks: nearest-neighbour regrid and the undeclared NoData strip ----------------------
from qeht.core.raster import resampling_stats, audit_nodata
src = gdal.Open(f"{D}/demo_dem.tif")
nn = gdal.Warp("/vsimem/nn.tif", src, xRes=27.0, yRes=27.0, resampleAlg="near", srcNodata=0, dstNodata=-9999)
bl = gdal.Warp("/vsimem/bl.tif", src, xRes=27.0, yRes=27.0, resampleAlg="bilinear", srcNodata=0, dstNodata=-9999)
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6))
stats = {}
for a, (dsx, tit) in zip(axs[:2], ((nn, "Nearest neighbour to 27 m"), (bl, "Bilinear to 27 m"))):
    arr = dsx.ReadAsArray().astype(float); v = arr > -9000
    stt = resampling_stats(arr, v); stats[tit] = stt
    win = arr[560:640, 520:600]
    a.imshow(ls.hillshade(win, vert_exag=6, dx=27, dy=27), cmap="gray")
    a.set_title(tit, loc="left", fontsize=8.5); a.set_xticks([]); a.set_yticks([])
    a.text(0.02, 0.03, f"repeated rows {100 * stt['duplicated_row_share']:.1f} %, columns {100 * stt['duplicated_col_share']:.1f} %",
           transform=a.transAxes, fontsize=6.8, color="white", bbox=dict(fc=INK, ec="none", alpha=0.75, pad=1.5))
raw = src.ReadAsArray().astype(float)
a = axs[2]
row = raw[380, -40:]
a.plot(np.arange(-40, 0), row, color=INK, lw=1.0, drawstyle="steps-mid")
a.set_title("Undeclared NoData strip", loc="left", fontsize=8.5)
a.set_xlabel("Columns from the east edge"); a.set_ylabel("Elevation (m)")
a.annotate("0 m, no NoData\nin the header", xy=(-3, 0), xytext=(-36, 500), fontsize=6.8, color=RED,
           arrowprops=dict(arrowstyle="->", color=RED, lw=0.8))
msg = audit_nodata(raw, np.ones(raw.shape, bool), None)
fig.subplots_adjust(wspace=0.35)
save(fig, "chart_dem_checks")
print(json.dumps({k: {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items()} for k, v in stats.items()}))
print(msg)
