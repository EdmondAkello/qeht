# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Stable outlet identifiers (`outlet_uid`).

Why this exists
---------------
Up to v0.8.3 every output carried `outlet_id` = the pour-point layer's
feature id. That integer changes when the layer is edited, re-saved or
filtered, so a downstream tool (HEAS) had to guess which catchment belonged
to which crossing. `outlet_uid` is a text identifier assigned ONCE per run
and written, unchanged, on every crossing, catchment and flow path.

ID schemes (decision D2, 2026-09-29)
------------------------------------
"sequential" (default)
    prefix + zero-padded number, e.g. X001, X002 ... Numbering order:
      order="downstream"  - largest contributing area first (the most
                            downstream crossing gets 001), ties broken by
                            input order. Used when no alignment is given.
      order="input"       - the order of the pour-point layer.
      order="chainage"    - ascending chainage values (supplied by the
                            caller; the crossing tool, WP-A, provides them).
    Chainage is never encoded in the ID.
"attribute"
    Taken verbatim from a pour-point attribute the user chooses (e.g.
    existing culvert numbers). Values are stripped; empty values are an
    error.

Uniqueness is enforced in every scheme; duplicates stop the run with a
list of the offending values rather than being renumbered silently.

No QGIS imports.
"""

ID_SCHEMES = ("sequential", "attribute")
ORDERS = ("downstream", "input", "chainage")


class OutletIdError(ValueError):
    """Raised when identifiers cannot be assigned or are not unique."""


def find_duplicates(uids):
    """Return a sorted list of values that occur more than once."""
    seen, dup = set(), set()
    for u in uids:
        if u in seen:
            dup.add(u)
        seen.add(u)
    return sorted(dup)


def validate_unique(uids, what="outlet_uid"):
    """Raise OutletIdError listing every duplicate/empty value."""
    empty = [i for i, u in enumerate(uids) if u is None or str(u).strip() == ""]
    if empty:
        raise OutletIdError(
            f"{len(empty)} {what} value(s) are empty (feature positions "
            f"{', '.join(str(i + 1) for i in empty[:20])}"
            f"{' ...' if len(empty) > 20 else ''}).")
    dup = find_duplicates(uids)
    if dup:
        raise OutletIdError(
            f"{what} must be unique; duplicated: {', '.join(map(str, dup[:50]))}"
            f"{' ...' if len(dup) > 50 else ''}. Fix the source values and run again.")
    return True


def numbering_order(n, order="downstream", accumulation=None, chainage=None):
    """Positions 0..n-1 in the order they should receive numbers."""
    if order not in ORDERS:
        raise OutletIdError(f"Unknown order '{order}'. Use one of {ORDERS}.")
    idx = list(range(n))
    if order == "input":
        return idx
    if order == "downstream":
        if accumulation is None or len(accumulation) != n:
            raise OutletIdError("order='downstream' needs the accumulation at each outlet.")
        # Largest area first; stable on input order for ties.
        return sorted(idx, key=lambda i: (-float(accumulation[i]), i))
    if chainage is None or len(chainage) != n:
        raise OutletIdError("order='chainage' needs a chainage for each outlet.")
    return sorted(idx, key=lambda i: (float(chainage[i]), i))


def sequential_uids(n, prefix="X", width=3, start=1, order="downstream",
                    accumulation=None, chainage=None):
    """Sequential IDs, returned in INPUT order (uids[i] belongs to outlet i).

    `width` is the minimum digit count; it grows automatically so that
    1000+ outlets never produce ambiguous IDs (X1000 after X999).
    """
    prefix = "" if prefix is None else str(prefix).strip()
    if n <= 0:
        return []
    last = int(start) + n - 1
    width = max(int(width), len(str(last)))
    uids = [None] * n
    for k, pos in enumerate(numbering_order(n, order, accumulation, chainage)):
        uids[pos] = f"{prefix}{int(start) + k:0{width}d}"
    validate_unique(uids)
    return uids


def attribute_uids(values, prefix=""):
    """IDs from a user attribute (one value per outlet, input order)."""
    prefix = "" if prefix is None else str(prefix).strip()
    uids = []
    for v in values:
        s = "" if v is None else str(v).strip()
        # QGIS NULL arrives as a QVariant-like object whose str() is 'NULL'
        if s.upper() == "NULL":
            s = ""
        uids.append(f"{prefix}{s}" if s else "")
    validate_unique(uids)
    return uids


def assign_uids(n, scheme="sequential", prefix="X", width=3, start=1,
                order="downstream", accumulation=None, chainage=None,
                attribute_values=None):
    """Single entry point used by the Processing tools. Returns (uids, info).

    `info` is a dict suitable for run metadata (id_scheme, id_prefix, ...).
    """
    if scheme == "attribute":
        if attribute_values is None:
            raise OutletIdError("scheme='attribute' needs attribute_values.")
        uids = attribute_uids(attribute_values, prefix=prefix)
        info = {"id_scheme": "attribute", "id_prefix": prefix or ""}
    elif scheme == "sequential":
        uids = sequential_uids(n, prefix=prefix, width=width, start=start,
                               order=order, accumulation=accumulation,
                               chainage=chainage)
        info = {"id_scheme": "sequential", "id_prefix": prefix or "",
                "id_order": order, "id_start": int(start),
                "id_width": max(int(width), len(str(int(start) + max(n, 1) - 1)))}
    else:
        raise OutletIdError(f"Unknown id scheme '{scheme}'. Use one of {ID_SCHEMES}.")
    return uids, info
