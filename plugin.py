# -*- coding: utf-8 -*-
# QEHT - Licensed under the GNU General Public License v2 or later.
"""Plugin shell. Deliberately thin: it registers the Processing provider
and does nothing else. All functionality lives in the algorithms."""

from qgis.core import QgsApplication

from .processing_provider.provider import QehtProvider


class QehtPlugin(object):

    def __init__(self, iface):
        self.iface = iface
        self.provider = None

    def initGui(self):
        self.initProcessing()

    def initProcessing(self):
        self.provider = QehtProvider()
        QgsApplication.processingRegistry().addProvider(self.provider)

    def unload(self):
        if self.provider is not None:
            QgsApplication.processingRegistry().removeProvider(self.provider)
            self.provider = None
