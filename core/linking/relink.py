# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Renumber and relink (decision D3).

After an exchange package has been produced, the engineer may delete, move
or add crossing points in QGIS. Re-running with this tool:
  * re-sorts the remaining crossings (by chainage when every crossing has
    one, else downstream-first),
  * re-issues gapless sequential outlet_uid values,
  * recomputes every link (catchment, flow path, all attributes),
  * writes `renumber_log` (old -> new) so external references can be
    updated.
IDs are only ever re-issued when the user runs this tool.

No QGIS imports.
"""

LOG_FIELDS = [("old_uid", "text"), ("new_uid", "text"), ("change", "text"),
              ("chainage_m", "real"), ("moved_m", "real")]


def renumber_log(crossings_new, old_by_fid, old_all_uids):
    """Build renumber_log rows.

    crossings_new : list of (geom, attrs) from build_exchange_records; attrs
        carry `outlet_id` = the fid of the edited crossing feature
    old_by_fid : {fid: {"uid": old outlet_uid or None, "x": .., "y": ..}}
        taken from the edited crossings layer (None uid = point added)
    old_all_uids : every outlet_uid in the package before editing (from its
        catchments layer), so deleted crossings can be listed

    change values: unchanged | renumbered | new | deleted; `moved_m` is the
    distance between the old outlet and the new outlet location.
    """
    rows, seen = [], set()
    for (x, y), a in crossings_new:
        old = old_by_fid.get(a.get("outlet_id"), {})
        ou = old.get("uid")
        ou = None if ou is None or str(ou).strip() in ("", "NULL") else str(ou)
        moved = None
        if old.get("x") is not None and old.get("y") is not None:
            moved = ((x - float(old["x"])) ** 2 + (y - float(old["y"])) ** 2) ** 0.5
        if ou is None:
            change = "new"
        else:
            seen.add(ou)
            change = "unchanged" if ou == a["outlet_uid"] else "renumbered"
        rows.append({"old_uid": ou, "new_uid": a["outlet_uid"], "change": change,
                     "chainage_m": a.get("chainage_m"), "moved_m": moved})
    for u in sorted(set(u for u in old_all_uids if u) - seen):
        rows.append({"old_uid": u, "new_uid": None, "change": "deleted",
                     "chainage_m": None, "moved_m": None})
    return rows
