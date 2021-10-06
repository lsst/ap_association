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

import unittest
import numpy as np
import pandas as pd

import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.pipe.base as pipeBase
import lsst.utils.tests
from unittest.mock import patch, Mock, DEFAULT

from lsst.ap.association import DiaPipelineTask


class TestDiaPipelineTask(unittest.TestCase):

    @classmethod
    def _makeDefaultConfig(cls, doPackageAlerts=False):
        config = DiaPipelineTask.ConfigClass()
        config.apdb.db_url = "sqlite://"
        config.doPackageAlerts = doPackageAlerts
        return config

    def setUp(self):
        # schemas are persisted in both Gen 2 and Gen 3 butler as prototypical catalogs
        srcSchema = afwTable.SourceTable.makeMinimalSchema()
        srcSchema.addField("base_PixelFlags_flag", type="Flag")
        srcSchema.addField("base_PixelFlags_flag_offimage", type="Flag")
        self.srcSchema = afwTable.SourceCatalog(srcSchema)

    def tearDown(self):
        pass

    def testRunWithAlerts(self):
        """Test running while creating and packaging alerts.
        """
        self._testRun(True)

    def testRunWithoutAlerts(self):
        """Test running without creating and packaging alerts.
        """
        self._testRun(False)

    def _testRun(self, doPackageAlerts=False):
        """Test the normal workflow of each ap_pipe step.
        """
        config = self._makeDefaultConfig(doPackageAlerts=doPackageAlerts)
        task = DiaPipelineTask(config=config)
        # Set DataFrame index testing to always return False. Mocks return
        # true for this check otherwise.
        task.testDataFrameIndex = lambda x: False
        diffIm = Mock(spec=afwImage.ExposureF)
        exposure = Mock(spec=afwImage.ExposureF)
        template = Mock(spec=afwImage.ExposureF)
        diaSrc = Mock(sepc=pd.DataFrame)
        ccdExposureIdBits = 32

        # Each of these subtasks should be called once during diaPipe
        # execution. We use mocks here to check they are being executed
        # appropriately.
        subtasksToMock = [
            "diaCatalogLoader",
            "associator",
            "diaCalculation",
            "diaForcedSource",
        ]
        if doPackageAlerts:
            subtasksToMock.append("alertPackager")
        else:
            self.assertFalse(hasattr(task, "alertPackager"))

        # apdb isn't a subtask, but still needs to be mocked out for correct
        # execution in the test environment.
        with patch.multiple(
            task, **{task: DEFAULT for task in subtasksToMock + ["apdb"]}
        ):
            result = task.run(diaSrc,
                              diffIm,
                              exposure,
                              template,
                              ccdExposureIdBits,
                              "g")
            for subtaskName in subtasksToMock:
                getattr(task, subtaskName).run.assert_called_once()
            pipeBase.testUtils.assertValidOutput(task, result)
            self.assertEqual(result.apdbMarker.db_url, "sqlite://")

    def test_createDiaObjects(self):
        """Test that creating new DiaObjects works as expected.
        """
        nSources = 5
        diaSources = pd.DataFrame(data=[
            {"ra": 0.04*idx, "decl": 0.04*idx,
             "diaSourceId": idx + 1 + nSources, "diaObjectId": 0,
             "ssObjectId": 0}
            for idx in range(nSources)])

        config = self._makeDefaultConfig(doPackageAlerts=False)
        task = DiaPipelineTask(config=config)
        result = task.createNewDiaObjects(diaSources)
        self.assertEqual(nSources, len(result.newDiaObjects))
        self.assertTrue(np.all(np.equal(
            result.diaSources["diaObjectId"].to_numpy(),
            result.diaSources["diaSourceId"].to_numpy())))
        self.assertTrue(np.all(np.equal(
            result.newDiaObjects["diaObjectId"].to_numpy(),
            result.diaSources["diaSourceId"].to_numpy())))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
