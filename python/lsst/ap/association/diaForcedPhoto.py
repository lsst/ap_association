#
# LSST Data Management System
# Copyright 2017 LSST/AURA.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

"""Class for performing forced photometry on calexps using DIASource positions.
"""


from lsst.meas.base import ForcedMeasurementTask
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase


class DiaForcedPhotoCcdConfig(pexConfig.Config):
    forcedMeasurement = pexConfig.ConfigurableField(
        target=ForcedMeasurementTask,
        doc='Specify where and how to load and store DIAObjects and '
        'DIASources.',
    )

    def setDefaults(self):
        self.forcedMeasurement.config.


class DiaForcedPhotoCcdTask(pipeBase.Task):
    
    onfigClass = DiaForcedPhotoCcdConfig
    _DefaultName = "diaForcedPhoto"

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.makeSubtask('forcedMeasurement')

    def run(self):
        pass

    
