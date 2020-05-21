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

import contextlib
import os
import unittest

import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.pipe.base as pipeBase
from lsst.utils import getPackageDir
import lsst.utils.tests
from unittest.mock import patch, Mock

from lsst.ap.association import DiaPipelineTask


class TestDiaPipelineTask(unittest.TestCase):

    @classmethod
    def _makeDefaultConfig(cls, doPackageAlerts=False):
        config = DiaPipelineTask.ConfigClass()
        config.apdb.db_url = "sqlite://"
        config.apdb.isolation_level = "READ_UNCOMMITTED"
        config.diaSourceDpddifier.copyColumns = {"id": "id",
                                                 "parent": "parent",
                                                 "coord_ra": "coord_ra",
                                                 "coord_dec": "coord_dec"}
        config.diaSourceDpddifier.flagMap = os.path.join(
            getPackageDir("ap_association"),
            "tests",
            "test-flag-map.yaml")
        config.doPackageAlerts = doPackageAlerts
        return config

    @contextlib.contextmanager
    def mockPatchSubtasks(self, task, doPackageAlerts=False):
        """Make mocks for all the ap_pipe subtasks.

        This is needed because the task itself cannot be a mock.
        The task's subtasks do not exist until the task is created, so
        this allows us to mock them instead.

        Parameters
        ----------
        task : `lsst.ap.association.DiaPipelineTask`
            The task whose subtasks will be mocked.

        Yields
        ------
        subtasks : `lsst.pipe.base.Struct`
            All mocks created by this context manager, including:

            ``diaCatalogLoader``
            ``dpddifier``
            ``associator``
            ``forcedSource``
                a mock for the corresponding subtask. Mocks do not return any
                particular value, but have mocked methods that can be queried
                for calls by ApPipeTask
        """
        if doPackageAlerts:
            with patch.object(task, "diaCatalogLoader") as mockDiaCatLoader, \
                    patch.object(task, "diaSourceDpddifier") as mockDpddifier, \
                    patch.object(task, "associator") as mockAssociator, \
                    patch.object(task, "diaForcedSource") as mockForcedSource, \
                    patch.object(task, "apdb") as mockApdb, \
                    patch.object(task, "alertPackager") as mockAlertPackager:
                yield pipeBase.Struct(diaCatalogLoader=mockDiaCatLoader,
                                      dpddifier=mockDpddifier,
                                      associator=mockAssociator,
                                      diaForcedSource=mockForcedSource,
                                      apdb=mockApdb,
                                      alertPackager=mockAlertPackager)
        else:
            with patch.object(task, "diaCatalogLoader") as mockDiaCatLoader, \
                    patch.object(task, "diaSourceDpddifier") as mockDpddifier, \
                    patch.object(task, "associator") as mockAssociator, \
                    patch.object(task, "diaForcedSource") as mockForcedSource, \
                    patch.object(task, "apdb") as mockApdb:
                yield pipeBase.Struct(diaCatalogLoader=mockDiaCatLoader,
                                      dpddifier=mockDpddifier,
                                      associator=mockAssociator,
                                      diaForcedSource=mockForcedSource,
                                      apdb=mockApdb)

    def setUp(self):
        self.srcSchema = afwTable.SourceTable.makeMinimalSchema()
        self.srcSchema.addField("base_PixelFlags_flag", type="Flag")
        self.srcSchema.addField("base_PixelFlags_flag_offimage", type="Flag")

    def tearDown(self):
        pass

    def testRunQuantum(self):
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
        task = DiaPipelineTask(
            config=config,
            initInputs={"diaSourceSchema": self.srcSchema})
        diffIm = Mock(spec=afwImage.ExposureF)
        exposure = Mock(spec=afwImage.ExposureF)
        diaSrc = Mock(sepc=afwTable.SourceCatalog)
        ccdExposureIdBits = 32
        with self.mockPatchSubtasks(task, doPackageAlerts) as subtasks:
            result = task.run(diaSrc, diffIm, exposure, ccdExposureIdBits)
            subtasks.dpddifier.run.assert_called_once()
            subtasks.dpddifier.run.assert_called_once()
            subtasks.associator.run.assert_called_once()
            subtasks.diaForcedSource.run.assert_called_once()
            if doPackageAlerts:
                subtasks.alertPackager.run.assert_called_once()
            else:
                self.assertFalse(hasattr(subtasks, 'alertPackager'))
            self.assertEqual(result.apdb_marker.db_url, "sqlite://")
            self.assertEqual(result.apdb_marker.isolation_level,
                             "READ_UNCOMMITTED")


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
