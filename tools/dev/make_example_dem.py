"""Make examples/example_dem.tif: a synthetic DEM with no real-world correspondence.

    python3.12 qeht/tools/dev/make_example_dem.py [output.tif]

760 x 720 cells at 30 m in EPSG:32737 (WGS 84 / UTM 37S) at an invented origin. The terrain is
built to exercise every QEHT tool the way a real DEM does:
  - a regional tilt (higher to the north) and a volcanic cone with a summit crater (a closed
    depression for Fill to resolve);
  - fractal relief (power-law spectrum), then dendritic valleys carved by QEHT's own flow
    accumulation over a few iterations, so the stream network looks like real drainage;
  - a low-relief plain in the south;
  - values stored in whole metres, like ALOS AW3D30, so the plain has real flats;
  - an undeclared NoData strip (value 0, no NoData in the header) along the east edge, which the
    DEM checks are expected to catch.
Seeded, so the file is reproducible.
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
from qeht.core.conditioning.fill import fill_depressions  # noqa: E402
from qeht.core.flow.direction import d8_direction  # noqa: E402
from qeht.core.flow.accumulation import flow_accumulation  # noqa: E402

ROWS, COLS, CS = 760, 720, 30.0
X0, Y0 = 410000.0, 9890000.0          # invented upper-left corner, UTM 37S
SEED = 2026


def fractal(rows, cols, beta, rng):
    ky = np.fft.fftfreq(rows)[:, None]
    kx = np.fft.rfftfreq(cols)[None, :]
    k = np.hypot(kx, ky)
    k[0, 0] = 1.0
    amp = k ** (-beta / 2.0)
    amp[0, 0] = 0.0
    phase = rng.uniform(0, 2 * np.pi, amp.shape)
    f = np.fft.irfft2(amp * np.exp(1j * phase), s=(rows, cols))
    return (f - f.mean()) / f.std()


def blur(a, sigma_cells):
    """Gaussian blur by FFT (valley width)."""
    ky = np.fft.fftfreq(a.shape[0])[:, None]
    kx = np.fft.rfftfreq(a.shape[1])[None, :]
    g = np.exp(-2 * (np.pi * sigma_cells) ** 2 * (kx ** 2 + ky ** 2))
    return np.fft.irfft2(np.fft.rfft2(a) * g, s=a.shape)


def build():
    rng = np.random.default_rng(SEED)
    r = np.arange(ROWS)[:, None] * CS
    c = np.arange(COLS)[None, :] * CS
    H, W = ROWS * CS, COLS * CS
    z = (1650.0 + 520.0 * (1.0 - r / H) ** 1.3) * np.ones((1, COLS))  # tilt, higher to the north
    vy, vx = 0.30 * H, 0.66 * W                                       # volcano
    d = np.hypot(r - vy, c - vx)
    z += 900.0 * np.exp(-(d / 4200.0) ** 1.6)
    z -= 160.0 * np.exp(-(d / 650.0) ** 4)                            # crater
    z += 420.0 * np.exp(-((r - 0.18 * H) ** 2 / (2 * 2600.0 ** 2) + (c - 0.2 * W) ** 2 / (2 * 3400.0 ** 2)))
    relief = 1.0 - 0.75 * np.clip((r - 0.62 * H) / (0.3 * H), 0, 1)  # plain in the south
    z += relief * (70.0 * fractal(ROWS, COLS, 3.8, rng) + 4.0 * fractal(ROWS, COLS, 2.4, rng))
    valid = np.ones(z.shape, bool)
    for it in range(5):                                               # carve valleys
        filled, _, _ = fill_depressions(z, valid, cell_width=CS, cell_height=CS)
        dirs, _ = d8_direction(filled, valid, CS, CS, flat_method="barnes")
        acc, _ = flow_accumulation(dirs, valid)
        area_km2 = (acc + 1) * CS * CS / 1e6
        cut = 16.0 * np.clip(np.log10(area_km2 / 0.15), 0, None)        # only channels > 0.15 km2
        z = z - relief * blur(cut, 1.6 + 0.4 * it)
    keep = np.exp(-(d / 700.0) ** 4) > 0.5                            # keep the crater closed
    z[keep] = np.minimum(z[keep], np.percentile(z[keep], 60))
    z = np.round(z)                                                    # whole metres, like ALOS
    z[:, -6:] = 0.0                                                    # undeclared NoData strip
    return z.astype(np.float32)


def write(path, z):
    from osgeo import gdal, osr
    gdal.UseExceptions()
    ds = gdal.GetDriverByName("GTiff").Create(path, COLS, ROWS, 1, gdal.GDT_Float32,
                                              options=["COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"])
    ds.SetGeoTransform((X0, CS, 0.0, Y0, 0.0, -CS))
    srs = osr.SpatialReference(); srs.ImportFromEPSG(32737)
    ds.SetProjection(srs.ExportToWkt())
    ds.SetMetadataItem("QEHT_DEM_SOURCE", "synthetic example (tools/dev/make_example_dem.py); "
                       "no real-world correspondence")
    ds.GetRasterBand(1).WriteArray(z)
    ds = None


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)),
                                                              "examples", "example_dem.tif")
    z = build()
    write(out, z)
    v = z[z > 0]
    print(f"wrote {out}: {ROWS} x {COLS} at {CS:g} m, elevation {v.min():.0f}-{v.max():.0f} m")
