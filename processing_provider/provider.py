# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Processing provider registration.

Everything is exposed as a QgsProcessingAlgorithm rather than a bespoke
dialog. That buys parameter widgets, validation, batch mode, the
Graphical Modeler, the history log and processing.run() scriptability
without writing any of it - and the Modeler is a big part of why similar
commercial hydrology toolsets are useful in the first place.
"""

import os

from qgis.core import QgsProcessingProvider
from qgis.PyQt.QtGui import QIcon

from .alg_fill import FillDepressionsAlgorithm
from .alg_flow_direction import FlowDirectionAlgorithm
from .alg_flow_accumulation import FlowAccumulationAlgorithm
from .alg_streams import StreamNetworkAlgorithm
from .alg_watershed import DelineateCatchmentAlgorithm
from .alg_longest_flowpath import LongestFlowPathAlgorithm
from .alg_streamlines import StreamLinesAlgorithm
from .alg_characteristics import CatchmentCharacteristicsAlgorithm
from .alg_build_exchange import BuildHeasExchangeAlgorithm


class QehtProvider(QgsProcessingProvider):

    def loadAlgorithms(self):
        for alg in (
            FillDepressionsAlgorithm(),
            FlowDirectionAlgorithm(),
            FlowAccumulationAlgorithm(),
            StreamNetworkAlgorithm(),
            DelineateCatchmentAlgorithm(),
            StreamLinesAlgorithm(),
            LongestFlowPathAlgorithm(),
            CatchmentCharacteristicsAlgorithm(),
            BuildHeasExchangeAlgorithm(),
        ):
            self.addAlgorithm(alg)

    def id(self):
        return "qeht"

    def name(self):
        return "Engineering Hydrology"

    def longName(self):
        return "QGIS Engineering Hydrology Toolkit"

    def icon(self):
        path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "icon.svg")
        return QIcon(path) if os.path.exists(path) else QgsProcessingProvider.icon(self)
