# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""The erosion output folder: file names, run record, and reading it back.

"Erosion indices and RUSLE" writes one folder; "Sample erosion along
alignment" and "Build HEAS exchange package" read it, so every downstream
product uses the same grids, factors and class schemes, recorded once in
erosion_run.json.
"""

import json
import os

import numpy as np

FILES = {
    "slope_deg": "slope_deg.tif", "tan_beta": "tan_beta.tif",
    "spec_catch_area": "spec_catch_area_m.tif", "spi": "spi.tif", "ln_spi": "ln_spi.tif",
    "twi": "twi.tif", "ls": "ls_factor.tif", "r": "r_factor.tif", "k": "k_factor.tif",
    "c": "c_factor.tif", "p": "p_factor.tif", "soil_loss": "rusle_soil_loss_t_ha_yr.tif",
    "spi_class": "spi_class.tif", "rusle_class": "rusle_class.tif", "ls_class": "ls_class.tif",
    "combined_class": "erosion_combined_class.tif", "channel": "channel_mask.tif",
    "sti": "sti_overland.tif",
}
RUN_JSON = "erosion_run.json"
EXTENTS_CSV = "class_extents.csv"


def path_of(folder, key):
    return os.path.join(folder, FILES[key])


def write_run(folder, record):
    with open(os.path.join(folder, RUN_JSON), "w", encoding="utf-8") as f:
        json.dump(record, f, indent=1, default=float)


def read_run(folder):
    p = os.path.join(folder, RUN_JSON)
    if not os.path.exists(p):
        raise FileNotFoundError(f"{p} not found - is this a folder written by "
                                "'Erosion indices and RUSLE'?")
    with open(p, encoding="utf-8") as f:
        return json.load(f)


def load_erosion_inputs(folder, info, read_grid):
    """ErosionInputs from an erosion folder, checked against the DEM grid `info`.

    read_grid(path) -> (array, valid, info) (e.g. core.raster.read_dem).
    Raises ValueError when the folder was computed on a different grid.
    """
    from .classes import ClassScheme
    from .summary import ErosionInputs
    run = read_run(folder)
    grids = {}
    for key in ("ln_spi", "ls", "twi", "tan_beta", "soil_loss", "k", "c", "p", "channel", "sti"):
        p = path_of(folder, key)
        if not os.path.exists(p):
            grids[key] = None
            continue
        a, v, gi = read_grid(p)
        if (gi.rows, gi.cols) != (info.rows, info.cols) or \
                any(abs(x - y) > 1e-6 * max(1.0, abs(y)) for x, y in zip(gi.geotransform, info.geotransform)):
            raise ValueError(f"{os.path.basename(p)} is on a different grid from the flow "
                             "direction raster. Run 'Erosion indices and RUSLE' on the same DEM.")
        grids[key] = np.where(v, a, np.nan)
    sch = run["schemes"]
    spi = ClassScheme(**{k: sch["spi"][k] for k in ("breaks", "names", "scores", "basis")})
    ru_key = "rusle" if grids["soil_loss"] is not None else "ls"
    ru = ClassScheme(**{k: sch[ru_key][k] for k in ("breaks", "names", "scores", "basis")})
    vol = ClassScheme(**{k: sch["volume"][k] for k in ("breaks", "names", "scores", "basis")})
    ch = grids["channel"]
    channel = np.nan_to_num(ch, nan=0.0) > 0 if ch is not None else np.zeros((info.rows, info.cols), bool)
    inp = ErosionInputs(grids["ln_spi"], grids["ls"], grids["twi"], grids["tan_beta"], channel,
                        spi, ru, soil_loss=grids["soil_loss"], k=grids["k"], c=grids["c"],
                        p=grids["p"], volume_scheme=vol, weights=run.get("mcdma_weights"),
                        bulk_kgm3=run.get("bulk_density_kgm3", 1400.0),
                        matrix=run.get("combination_matrix"))
    inp.sti = grids.get("sti")
    dep = run.get("deposition") or {}
    inp.dep_near_m = float(dep.get("near_m", 100.0))
    inp.dep_far_m = float(dep.get("far_m", 500.0))
    inp.dep_breaks = tuple(dep.get("breaks", (0.7, 1.3)))
    return inp, run
