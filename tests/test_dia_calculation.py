# This file is part of ap_association.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import pandas as pd
import unittest

from lsst.ap.association import (
    DiaObjectCalculationTask,
    DiaObjectCalculationConfig,
    DiaObjectCalculationPlugin)
from lsst.meas.base.pluginRegistry import register
import lsst.utils.tests


@register("testDiaPlugin")
class TestDiaPlugin(DiaObjectCalculationPlugin):

    def __init__(self, config, name, metadata):
        DiaObjectCalculationPlugin.__init__(self, config, name, metadata)

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self, diaObject, psFluxes, **kwargs):
        """
        """
        diaObject["meanFlux"] = np.mean(psFluxes)
        diaObject["stdFlux"] = np.std(psFluxes, ddof=1)


@register("testDependentDiaPlugin")
class TestDependentDiaPlugin(DiaObjectCalculationPlugin):

    def __init__(self, config, name, metadata):
        DiaObjectCalculationPlugin.__init__(self, config, name, metadata)

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_MOMMENTS_CALCULATED

    def calculate(self, diaObject, psFluxes, psFluxErrs, **kwargs):
        """
        """
        diaObject["chiFlux"] = np.sum(
            ((psFluxes - diaObject["meanFlux"]) / psFluxErrs) ** 2)


class TestMeanPosition(unittest.TestCase):

    def setUp(self):
        self.diaObjects = pd.DataFrame(
            data=[{"diaObjectId": objId} for objId in range(5)])
        self.diaObjects.set_index("diaObjectId", inplace=True)

        diaSources = [{"diaSourceId": objId, "diaObjectId": objId,
                       "psFlux": 0., "psFluxErr": 1.,
                       "totFlux": 0., "totFluxErr": 1.,
                       "midPointTai": 0, "filterName": "g"}
                      for objId in range(5)]
        diaSources.extend([{"diaSourceId": 5 + objId, "diaObjectId": objId,
                            "psFlux": 0., "psFluxErr": 1.,
                            "totFlux": 0., "totFluxErr": 1.,
                            "midPointTai": 0, "filterName": "r"}
                           for objId in range(5)])
        diaSources.extend([{"diaSourceId": 10, "diaObjectId": 0,
                            "psFlux": 1., "psFluxErr": 1.,
                            "totFlux": 0., "totFluxErr": 0.,
                            "midPointTai": 0, "filterName": "g"},
                           {"diaSourceId": 11, "diaObjectId": 1,
                            "psFlux": 1., "psFluxErr": 1.,
                            "totFlux": 0., "totFluxErr": 0.,
                            "midPointTai": 0, "filterName": "g"},
                           {"diaSourceId": 12, "diaObjectId": 2,
                            "psFlux": np.nan, "psFluxErr": 1.,
                            "totFlux": 0., "totFluxErr": 0.,
                            "midPointTai": 0, "filterName": "g"},
                           {"diaSourceId": 13, "diaObjectId": 13,
                            "psFlux": 1., "psFluxErr": 1.,
                            "totFlux": 0., "totFluxErr": 0.,
                            "midPointTai": 0, "filterName": "g"}])
        self.diaSources = pd.DataFrame(data=diaSources)
        self.diaSources.set_index(["diaObjectId", "filterName", "diaSourceId"],
                                  inplace=True)

        self.newDiaObjectId = 13
        self.updatedDiaObjectIds = np.array([0, 1, 2, 13], dtype=np.int)

        conf = DiaObjectCalculationConfig()
        conf.plugins = ["testDiaPlugin",
                        "testDependentDiaPlugin"]
        self.diaObjCalTask = DiaObjectCalculationTask(config=conf)

    def testRun(self):
        """Test the run method and that diaObjects are updated correctly.
        """
        results = self.diaObjCalTask.run(self.diaObjects,
                                         self.diaSources,
                                         self.updatedDiaObjectIds,
                                         "g")
        diaObjectCat = results.diaObjectCat
        updatedDiaObjects = results.updatedDiaObjects
        updatedDiaObjects.set_index("diaObjectId", inplace=True)

        self.assertEqual(len(diaObjectCat), len(self.diaObjects) + 1)
        self.assertEqual(len(updatedDiaObjects),
                         len(self.updatedDiaObjectIds))

        for objId, diaObject in updatedDiaObjects.iterrows():
            if objId == self.newDiaObjectId:
                self.assertEqual(diaObject["meanFlux"], 1.)
                self.assertTrue(np.isnan(diaObject["stdFlux"]))
                self.assertAlmostEqual(diaObject["chiFlux"], 0.0)
            elif objId == 2:
                self.assertAlmostEqual(diaObject["meanFlux"], 0.0)
                self.assertTrue(np.isnan(diaObject["stdFlux"]))
                self.assertAlmostEqual(diaObject["chiFlux"], 0.0)
            else:
                self.assertAlmostEqual(diaObject["meanFlux"], 0.5)
                self.assertAlmostEqual(diaObject["stdFlux"],
                                       0.7071067811865476)
                self.assertAlmostEqual(diaObject["chiFlux"], 0.5)


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
