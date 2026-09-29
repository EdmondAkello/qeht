# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Peak-memory estimate for a run, from measured per-step peaks (v0.12).

QEHT holds the whole grid in memory. The bytes-per-cell figures below were
measured with tracemalloc on a 2.25-million-cell surface (30 % filled
cells) for the v0.12 vectorised core; RESIDENT covers the arrays a tool
keeps while a step runs (input DEM, validity mask, outputs). They are
estimates, rounded up - use them to decide whether to clip, not as a
guarantee. No QGIS imports.
"""

import os

STEP_BYTES_PER_CELL = {
    "fill": 44,
    "flow_direction": 119,           # toward-lower flats
    "flow_direction_barnes": 165,    # Barnes / hybrid flats
    "accumulation": 82,
    "strahler": 83,
    "catchment": 82,
    "longest_flow_path": 82,
    "slope": 33,
}
RESIDENT_BYTES_PER_CELL = 24


def estimate_peak_memory(rows, cols, steps=None):
    """{'steps_gb': {step: GB}, 'peak_gb': max, 'megapixels': ...}."""
    n = float(rows) * float(cols)
    steps = steps or list(STEP_BYTES_PER_CELL)
    per = {s: n * (STEP_BYTES_PER_CELL[s] + RESIDENT_BYTES_PER_CELL) / 1e9 for s in steps}
    return {"steps_gb": per, "peak_gb": max(per.values()), "megapixels": n / 1e6}


def available_memory_gb():
    """Free physical memory in GB (best effort; None when unknown)."""
    try:
        if os.name == "nt":
            import ctypes

            class _MS(ctypes.Structure):
                _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                            ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                            ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                            ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                            ("sullAvailExtendedVirtual", ctypes.c_ulonglong)]
            ms = _MS(); ms.dwLength = ctypes.sizeof(_MS)
            if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(ms)):
                return ms.ullAvailPhys / 1e9
            return None
        pages = os.sysconf("SC_AVPHYS_PAGES")
        return pages * os.sysconf("SC_PAGE_SIZE") / 1e9
    except (AttributeError, ValueError, OSError):
        return None
