# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Soil erodibility (USLE/RUSLE K), texture class and hydrologic-group proxy.

Williams (EPIC) K  -  default
-----------------------------
Sharpley & Williams (1990), as used in SWAT (Neitsch et al. 2011, eq. 4:1.1.5):

    K = f_csand * f_cl-si * f_orgc * f_hisand            (US customary units)

    f_csand  = 0.2 + 0.3 exp[-0.0256 SAN (1 - SIL/100)]
    f_cl-si  = (SIL / (CLA + SIL)) ^ 0.3
    f_orgc   = 1 - 0.25 C / (C + exp(3.72 - 2.95 C))
    f_hisand = 1 - 0.7 SN1 / (SN1 + exp(-5.51 + 22.9 SN1)),   SN1 = 1 - SAN/100

SAN, SIL, CLA in %, C = organic carbon in %. The coefficient in f_csand is
0.0256; the SWAT 2009 theory PDF misprints it as 0.256, which drives f_csand
to 0.2 for almost any sandy soil. K_SI = 0.1317 x K_US, in
t.ha.h/(ha.MJ.mm) - the unit of the ESDAC global K rasters.

RUSLE geometric-mean-diameter K  -  optional
--------------------------------------------
Renard et al. (1997), for soils without structure/permeability data:

    K_SI = 0.0034 + 0.0405 exp[-0.5 ((log10 Dg + 1.659) / 0.7101)^2]
    Dg (mm) = exp(0.01 sum f_i ln m_i),  m = 0.001 (clay), 0.026 (silt),
              1.025 (sand) mm

Coarse fragments
----------------
CFRG = exp(-0.053 x rock %) (SWAT, MUSLE coarse-fragment factor), reported
separately - it multiplies the MUSLE product, it is not folded into K.

Hydrologic soil group - PROXY
-----------------------------
From texture, then drainage: A = sand, loamy sand; B = sandy loam, loam,
silt loam, silt; C = sandy clay loam; D = clay loam, silty clay loam, sandy
clay, silty clay, clay. Poorly / very poorly drained (P, V) -> D;
imperfectly drained (I) moves A/B to C. This is a texture-based proxy, not
a measured infiltration class, and every output labels it as such.

No QGIS imports, no GDAL.
"""

import math

K_US_TO_SI = 0.1317
NAN = float("nan")


def _ok(*vals):
    return all(v is not None and math.isfinite(v) and v >= 0 for v in vals)


def williams_k(sand, silt, clay, oc_pct):
    """Williams/EPIC K. Returns (K_US, K_SI, factors dict); NaN if inputs missing."""
    if not _ok(sand, silt, clay, oc_pct) or (silt + clay) <= 0:
        return NAN, NAN, {}
    f_csand = 0.2 + 0.3 * math.exp(-0.0256 * sand * (1.0 - silt / 100.0))
    f_clsi = (silt / (clay + silt)) ** 0.3
    c = oc_pct
    f_orgc = 1.0 - 0.25 * c / (c + math.exp(3.72 - 2.95 * c))
    sn1 = 1.0 - sand / 100.0
    f_hisand = 1.0 - 0.7 * sn1 / (sn1 + math.exp(-5.51 + 22.9 * sn1))
    k_us = f_csand * f_clsi * f_orgc * f_hisand
    return k_us, k_us * K_US_TO_SI, {"f_csand": f_csand, "f_clsi": f_clsi,
                                     "f_orgc": f_orgc, "f_hisand": f_hisand}


def dg_k(sand, silt, clay):
    """RUSLE K from geometric mean particle diameter (SI units). NaN if missing."""
    if not _ok(sand, silt, clay):
        return NAN
    tot = sand + silt + clay
    if tot <= 0:
        return NAN
    f = [clay / tot * 100.0, silt / tot * 100.0, sand / tot * 100.0]
    m = [0.001, 0.026, 1.025]
    dg = math.exp(0.01 * sum(fi * math.log(mi) for fi, mi in zip(f, m)))
    return 0.0034 + 0.0405 * math.exp(-0.5 * ((math.log10(dg) + 1.659) / 0.7101) ** 2)


def cfrg(rock_pct):
    """Coarse-fragment factor exp(-0.053 rock%). 1.0 when rock is unknown/0."""
    if rock_pct is None or not math.isfinite(rock_pct) or rock_pct < 0:
        return NAN
    return math.exp(-0.053 * rock_pct)


def usda_texture(sand, silt, clay):
    """USDA texture class name from %sand, %silt, %clay (normalised to 100)."""
    if not _ok(sand, silt, clay) or (sand + silt + clay) <= 0:
        return None
    t = sand + silt + clay
    sa, si, cl = 100.0 * sand / t, 100.0 * silt / t, 100.0 * clay / t
    if si + 1.5 * cl < 15:
        return "sand"
    if si + 1.5 * cl < 30:
        return "loamy sand"
    if (7 <= cl < 20 and sa > 52 and si + 2 * cl >= 30) or (cl < 7 and si < 50 and si + 2 * cl >= 30):
        return "sandy loam"
    if 7 <= cl < 27 and 28 <= si < 50 and sa <= 52:
        return "loam"
    if (si >= 50 and 12 <= cl < 27) or (50 <= si < 80 and cl < 12):
        return "silt loam"
    if si >= 80 and cl < 12:
        return "silt"
    if 20 <= cl < 35 and si < 28 and sa > 45:
        return "sandy clay loam"
    if 27 <= cl < 40 and 20 < sa <= 45:
        return "clay loam"
    if 27 <= cl < 40 and sa <= 20:
        return "silty clay loam"
    if cl >= 35 and sa > 45:
        return "sandy clay"
    if cl >= 40 and si >= 40:
        return "silty clay"
    if cl >= 40 and sa <= 45 and si < 40:
        return "clay"
    return "loam"   # boundary rounding


_HSG_TEXTURE = {
    "sand": "A", "loamy sand": "A",
    "sandy loam": "B", "loam": "B", "silt loam": "B", "silt": "B",
    "sandy clay loam": "C",
    "clay loam": "D", "silty clay loam": "D", "sandy clay": "D",
    "silty clay": "D", "clay": "D",
}


def hsg_proxy(texture, drain=None):
    """Hydrologic soil group PROXY from texture and FAO drainage class."""
    g = _HSG_TEXTURE.get(texture)
    if g is None:
        return None
    d = (drain or "").strip().upper()[:1]
    if d in ("P", "V"):
        return "D"
    if d == "I" and g in ("A", "B"):
        return "C"
    return g


# -- vectorised versions (v0.15, soil rasters) -------------------------------
# Same formulas as the scalar functions above, on NumPy arrays; NaN where an
# input is missing. Tested against the scalar versions in tests/test_soils_any.

import numpy as _np

TEXTURES = ("sand", "loamy sand", "sandy loam", "loam", "silt loam", "silt",
            "sandy clay loam", "clay loam", "silty clay loam", "sandy clay",
            "silty clay", "clay")
HSG_CODES = ("A", "B", "C", "D")        # code 1..4; 0 = unknown


def _arr(*vals):
    out = [_np.asarray(v, dtype=float) for v in vals]
    ok = _np.ones(_np.broadcast(*out).shape, bool)
    for v in out:
        ok &= _np.isfinite(v) & (v >= 0)
    return out, ok


def williams_k_arr(sand, silt, clay, oc_pct):
    """K_SI per cell (Williams/EPIC x 0.1317)."""
    (sa, si, cl, c), ok = _arr(sand, silt, clay, oc_pct)
    ok &= (si + cl) > 0
    with _np.errstate(all="ignore"):
        f_csand = 0.2 + 0.3 * _np.exp(-0.0256 * sa * (1.0 - si / 100.0))
        f_clsi = (si / (cl + si)) ** 0.3
        f_orgc = 1.0 - 0.25 * c / (c + _np.exp(3.72 - 2.95 * c))
        sn1 = 1.0 - sa / 100.0
        f_hisand = 1.0 - 0.7 * sn1 / (sn1 + _np.exp(-5.51 + 22.9 * sn1))
        k = f_csand * f_clsi * f_orgc * f_hisand * K_US_TO_SI
    return _np.where(ok, k, _np.nan)


def dg_k_arr(sand, silt, clay):
    (sa, si, cl), ok = _arr(sand, silt, clay)
    tot = sa + si + cl
    ok &= tot > 0
    with _np.errstate(all="ignore"):
        s = (cl / tot * 100.0 * math.log(0.001) + si / tot * 100.0 * math.log(0.026)
             + sa / tot * 100.0 * math.log(1.025))
        dg = _np.exp(0.01 * s)
        k = 0.0034 + 0.0405 * _np.exp(-0.5 * ((_np.log10(dg) + 1.659) / 0.7101) ** 2)
    return _np.where(ok, k, _np.nan)


def usda_texture_arr(sand, silt, clay):
    """Texture code per cell: 1..12 = TEXTURES index + 1, 0 = unknown."""
    (sand, silt, clay), ok = _arr(sand, silt, clay)
    t = sand + silt + clay
    ok &= t > 0
    with _np.errstate(all="ignore"):
        sa, si, cl = 100.0 * sand / t, 100.0 * silt / t, 100.0 * clay / t
    conds = [
        si + 1.5 * cl < 15,
        si + 1.5 * cl < 30,
        ((7 <= cl) & (cl < 20) & (sa > 52) & (si + 2 * cl >= 30)) | ((cl < 7) & (si < 50) & (si + 2 * cl >= 30)),
        (7 <= cl) & (cl < 27) & (28 <= si) & (si < 50) & (sa <= 52),
        ((si >= 50) & (12 <= cl) & (cl < 27)) | ((50 <= si) & (si < 80) & (cl < 12)),
        (si >= 80) & (cl < 12),
        (20 <= cl) & (cl < 35) & (si < 28) & (sa > 45),
        (27 <= cl) & (cl < 40) & (20 < sa) & (sa <= 45),
        (27 <= cl) & (cl < 40) & (sa <= 20),
        (cl >= 35) & (sa > 45),
        (cl >= 40) & (si >= 40),
        (cl >= 40) & (sa <= 45) & (si < 40),
    ]
    code = _np.select(conds, list(range(1, 13)), default=4)     # "loam" (boundary rounding)
    return _np.where(ok, code, 0).astype(_np.int8)


_TEX_HSG = {"sand": 1, "loamy sand": 1, "sandy loam": 2, "loam": 2, "silt loam": 2, "silt": 2,
            "sandy clay loam": 3, "clay loam": 4, "silty clay loam": 4, "sandy clay": 4,
            "silty clay": 4, "clay": 4}


def hsg_proxy_arr(texture_code, drain=None):
    """HSG proxy code per cell (1..4 = A..D, 0 unknown) from texture code and
    an optional array of FAO drainage letters (as single-character strings)."""
    lut = _np.array([0] + [_TEX_HSG[t] for t in TEXTURES], dtype=_np.int8)
    g = lut[_np.asarray(texture_code, dtype=_np.int64)]
    if drain is not None:
        d = _np.char.upper(_np.char.strip(_np.asarray(drain, dtype=str)))
        first = _np.array([x[:1] for x in d.ravel()]).reshape(d.shape) if d.size else d
        g = _np.where((g > 0) & _np.isin(first, ("P", "V")), 4, g)
        g = _np.where((g > 0) & (g <= 2) & (first == "I"), 3, g)
    return g.astype(_np.int8)
