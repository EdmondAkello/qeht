# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.


def classFactory(iface):
    from .plugin import QehtPlugin
    return QehtPlugin(iface)
