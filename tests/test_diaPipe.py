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
from lsst.pipe.base.testUtils import assertValidOutput
import lsst.utils.tests
import lsst.utils.timer
from unittest.mock import patch, Mock, MagicMock, DEFAULT

from lsst.ap.association import DiaPipelineTask


class TestDiaPipelineTask(unittest.TestCase):

    @classmethod
    def _makeDefaultConfig(cls,
                           doPackageAlerts=False,
                           doSolarSystemAssociation=False):
        config = DiaPipelineTask.ConfigClass()
        config.apdb.db_url = "sqlite://"
        config.doPackageAlerts = doPackageAlerts
        config.doSolarSystemAssociation = doSolarSystemAssociation
        return config

    def setUp(self):
        # schemas are persisted in both Gen 2 and Gen 3 butler as prototypical catalogs
        srcSchema = afwTable.SourceTable.makeMinimalSchema()
        srcSchema.addField("base_PixelFlags_flag", type="Flag")
        srcSchema.addField("base_PixelFlags_flag_offimage", type="Flag")
        self.srcSchema = afwTable.SourceCatalog(srcSchema)

    def tearDown(self):
        pass

    def testRun(self):
        """Test running while creating and packaging alerts.
        """
        self._testRun(doPackageAlerts=True, doSolarSystemAssociation=True)

    def testRunWithSolarSystemAssociation(self):
        """Test running while creating and packaging alerts.
        """
        self._testRun(doPackageAlerts=False, doSolarSystemAssociation=True)

    def testRunWithAlerts(self):
        """Test running while creating and packaging alerts.
        """
        self._testRun(doPackageAlerts=True, doSolarSystemAssociation=False)

    def testRunWithoutAlertsOrSolarSystem(self):
        """Test running without creating and packaging alerts.
        """
        self._testRun(doPackageAlerts=False, doSolarSystemAssociation=False)

    def _testRun(self, doPackageAlerts=False, doSolarSystemAssociation=False):
        """Test the normal workflow of each ap_pipe step.
        """
        config = self._makeDefaultConfig(
            doPackageAlerts=doPackageAlerts,
            doSolarSystemAssociation=doSolarSystemAssociation)
        task = DiaPipelineTask(config=config)
        # Set DataFrame index testing to always return False. Mocks return
        # true for this check otherwise.
        task.testDataFrameIndex = lambda x: False
        diffIm = Mock(spec=afwImage.ExposureF)
        exposure = Mock(spec=afwImage.ExposureF)
        template = Mock(spec=afwImage.ExposureF)
        diaSrc = MagicMock(spec=pd.DataFrame())
        ssObjects = MagicMock(spec=pd.DataFrame())
        ccdExposureIdBits = 32

        # Each of these subtasks should be called once during diaPipe
        # execution. We use mocks here to check they are being executed
        # appropriately.
        subtasksToMock = [
            "diaCatalogLoader",
            "diaCalculation",
            "diaForcedSource",
        ]
        if doPackageAlerts:
            subtasksToMock.append("alertPackager")
        else:
            self.assertFalse(hasattr(task, "alertPackager"))

        if not doSolarSystemAssociation:
            self.assertFalse(hasattr(task, "solarSystemAssociator"))

        def concatMock(_data, **_kwargs):
            return MagicMock(spec=pd.DataFrame)

        # Mock out the run() methods of these two Tasks to ensure they
        # return data in the correct form.
        @lsst.utils.timer.timeMethod
        def solarSystemAssociator_run(self, unAssocDiaSources, solarSystemObjectTable, diffIm):
            return lsst.pipe.base.Struct(nTotalSsObjects=42,
                                         nAssociatedSsObjects=30,
                                         ssoAssocDiaSources=MagicMock(spec=pd.DataFrame()),
                                         unAssocDiaSources=MagicMock(spec=pd.DataFrame()))

        @lsst.utils.timer.timeMethod
        def associator_run(self, table, diaObjects, exposure_time=None):
            return lsst.pipe.base.Struct(nUpdatedDiaObjects=2, nUnassociatedDiaObjects=3,
                                         matchedDiaSources=MagicMock(spec=pd.DataFrame()),
                                         unAssocDiaSources=MagicMock(spec=pd.DataFrame()))

        # apdb isn't a subtask, but still needs to be mocked out for correct
        # execution in the test environment.
        with patch.multiple(
            task, **{task: DEFAULT for task in subtasksToMock + ["apdb"]}
        ):
            with patch('lsst.ap.association.diaPipe.pd.concat', new=concatMock), \
                patch('lsst.ap.association.association.AssociationTask.run', new=associator_run), \
                patch('lsst.ap.association.ssoAssociation.SolarSystemAssociationTask.run',
                      new=solarSystemAssociator_run):

                result = task.run(diaSrc,
                                  ssObjects,
                                  diffIm,
                                  exposure,
                                  template,
                                  ccdExposureIdBits,
                                  "g")
                for subtaskName in subtasksToMock:
                    getattr(task, subtaskName).run.assert_called_once()
                assertValidOutput(task, result)
                self.assertEqual(result.apdbMarker.db_url, "sqlite://")
                meta = task.getFullMetadata()
                # Check that the expected metadata has been set.
                self.assertEqual(meta["diaPipe.numUpdatedDiaObjects"], 2)
                self.assertEqual(meta["diaPipe.numUnassociatedDiaObjects"], 3)
                # and that associators ran once or not at all.
                self.assertEqual(len(meta.getArray("diaPipe:associator.associator_runEndUtc")), 1)
                if doSolarSystemAssociation:
                    self.assertEqual(len(meta.getArray("diaPipe:solarSystemAssociator."
                                                       "solarSystemAssociator_runEndUtc")), 1)
                else:
                    self.assertNotIn("diaPipe:solarSystemAssociator", meta)

    def test_createDiaObjects(self):
        """Test that creating new DiaObjects works as expected.
        """
        nSources = 5
        diaSources = pd.DataFrame(data=[
            {"ra": 0.04*idx, "dec": 0.04*idx,
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
