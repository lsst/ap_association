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

import tempfile
import unittest
from unittest.mock import patch, Mock, MagicMock, DEFAULT
import warnings

import numpy as np
import pandas as pd

import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.dax.apdb as daxApdb
from lsst.meas.base import IdGenerator
import lsst.pex.config as pexConfig
import lsst.utils.tests
from lsst.pipe.base.testUtils import assertValidOutput

from lsst.ap.association import DiaPipelineTask
from utils_tests import makeExposure, makeDiaObjects, makeDiaSources, makeDiaForcedSources, \
    makeSolarSystemSources


def _makeMockDataFrame():
    """Create a new mock of a DataFrame.

    Returns
    -------
    mock : `unittest.mock.Mock`
        A mock guaranteed to accept all operations used by `pandas.DataFrame`.
    """
    with warnings.catch_warnings():
        # spec triggers deprecation warnings on DataFrame, but will
        # automatically adapt to any removals.
        warnings.simplefilter("ignore", category=DeprecationWarning)
        return MagicMock(spec=pd.DataFrame())


class TestDiaPipelineTask(unittest.TestCase):

    @classmethod
    def _makeDefaultConfig(cls, config_file, **kwargs):
        config = DiaPipelineTask.ConfigClass()
        config.doConfigureApdb = False
        config.apdb_config_url = config_file
        config.update(**kwargs)
        return config

    def setUp(self):
        # Create an instance of random generator with fixed seed.
        rng = np.random.default_rng(1234)
        self.rng = rng

        # schemas are persisted in both Gen 2 and Gen 3 butler as prototypical catalogs
        srcSchema = afwTable.SourceTable.makeMinimalSchema()
        srcSchema.addField("base_PixelFlags_flag", type="Flag")
        srcSchema.addField("base_PixelFlags_flag_offimage", type="Flag")
        self.srcSchema = afwTable.SourceCatalog(srcSchema)
        self.exposure = makeExposure(False, False)
        self.diaObjects = makeDiaObjects(20, self.exposure, rng)
        self.diaSources = makeDiaSources(
            100, self.diaObjects["diaObjectId"].to_numpy(), self.exposure, rng)
        self.diaForcedSources = makeDiaForcedSources(
            200, self.diaObjects["diaObjectId"].to_numpy(), self.exposure, rng)
        self.ssSources = makeSolarSystemSources(
            20, self.diaObjects["diaObjectId"].to_numpy(), self.exposure, rng)

        apdb_config = daxApdb.ApdbSql.init_database(db_url="sqlite://")
        self.config_file = tempfile.NamedTemporaryFile()
        self.addCleanup(self.config_file.close)
        apdb_config.save(self.config_file.name)

    # TODO: remove on DM-43419
    def testConfigApdbNestedOk(self):
        config = DiaPipelineTask.ConfigClass()
        config.doConfigureApdb = True
        with self.assertWarns(FutureWarning):
            config.apdb.db_url = "sqlite://"
            config.freeze()
            config.validate()

    # TODO: remove on DM-43419
    def testConfigApdbNestedInvalid(self):
        config = DiaPipelineTask.ConfigClass()
        config.doConfigureApdb = True
        # Don't set db_url
        config.freeze()
        with self.assertRaises(pexConfig.FieldValidationError):
            config.validate()

    # TODO: remove on DM-43419
    def testConfigApdbFileOk(self):
        config = DiaPipelineTask.ConfigClass()
        config.doConfigureApdb = False
        config.apdb_config_url = "some/file/path.yaml"
        config.freeze()
        config.validate()

    # TODO: remove on DM-43419
    def testConfigApdbFileInvalid(self):
        config = DiaPipelineTask.ConfigClass()
        config.doConfigureApdb = False
        # Don't set apdb_config_url
        config.freeze()
        with self.assertRaises(pexConfig.FieldValidationError):
            config.validate()

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

    def _testRun(self, doPackageAlerts=False, doSolarSystemAssociation=False, doRunForcedMeasurement=False):
        """Test the normal workflow of each ap_pipe step.
        """
        config = self._makeDefaultConfig(
            config_file=self.config_file.name,
            doPackageAlerts=doPackageAlerts,
            doSolarSystemAssociation=doSolarSystemAssociation,
            doRunForcedMeasurement=doRunForcedMeasurement,
        )
        task = DiaPipelineTask(config=config)
        # Set DataFrame index testing to always return False. Mocks return
        # true for this check otherwise.
        task.testDataFrameIndex = lambda x: False
        diffIm = Mock(spec=afwImage.ExposureF)
        exposure = Mock(spec=afwImage.ExposureF)
        template = Mock(spec=afwImage.ExposureF)
        diaSrc = _makeMockDataFrame()
        ssObjects = _makeMockDataFrame()

        # Each of these subtasks should be called once during diaPipe
        # execution. We use mocks here to check they are being executed
        # appropriately.
        subtasksToMock = [
            "diaCalculation",
        ]
        if doPackageAlerts:
            subtasksToMock.append("alertPackager")
        else:
            self.assertFalse(hasattr(task, "alertPackager"))

        if not doSolarSystemAssociation:
            self.assertFalse(hasattr(task, "solarSystemAssociator"))

        def concatMock(_data, **_kwargs):
            return _makeMockDataFrame()

        # Mock out the run() methods of these two Tasks to ensure they
        # return data in the correct form.
        def solarSystemAssociator_run(unAssocDiaSources, solarSystemObjectTable, diffIm):
            return lsst.pipe.base.Struct(nTotalSsObjects=42,
                                         nAssociatedSsObjects=30,
                                         ssoAssocDiaSources=_makeMockDataFrame(),
                                         unAssocDiaSources=_makeMockDataFrame())

        def associator_run(table, diaObjects):
            return lsst.pipe.base.Struct(nUpdatedDiaObjects=2, nUnassociatedDiaObjects=3,
                                         matchedDiaSources=_makeMockDataFrame(),
                                         unAssocDiaSources=_makeMockDataFrame())

        def updateObjectTableMock(diaObjects, diaSources):
            pass

        # apdb isn't a subtask, but still needs to be mocked out for correct
        # execution in the test environment.
        with patch.multiple(task, **{task: DEFAULT for task in subtasksToMock + ["apdb"]}), \
            patch('lsst.ap.association.diaPipe.pd.concat', side_effect=concatMock), \
            patch('lsst.ap.association.diaPipe.DiaPipelineTask.updateObjectTable',
                  side_effect=updateObjectTableMock), \
            patch('lsst.ap.association.association.AssociationTask.run',
                  side_effect=associator_run) as mainRun, \
            patch('lsst.ap.association.ssoAssociation.SolarSystemAssociationTask.run',
                  side_effect=solarSystemAssociator_run) as ssRun:

            result = task.run(diaSrc,
                              None,
                              diffIm,
                              exposure,
                              template,
                              self.diaObjects,
                              self.diaSources,
                              self.diaForcedSources,
                              "g",
                              IdGenerator(),
                              solarSystemObjectTable=ssObjects)
            for subtaskName in subtasksToMock:
                getattr(task, subtaskName).run.assert_called_once()
            assertValidOutput(task, result)
            # Exact type and contents of apdbMarker are undefined.
            self.assertIsInstance(result.apdbMarker, pexConfig.Config)
            meta = task.getFullMetadata()
            # Check that the expected metadata has been set.
            self.assertEqual(meta["diaPipe.numUpdatedDiaObjects"], 2)
            self.assertEqual(meta["diaPipe.numUnassociatedDiaObjects"], 3)
            # and that associators ran once or not at all.
            mainRun.assert_called_once()
            if doSolarSystemAssociation:
                ssRun.assert_called_once()
            else:
                ssRun.assert_not_called()

    def test_createDiaObjects(self):
        """Test that creating new DiaObjects works as expected.
        """
        nSources = 5
        diaSources = pd.DataFrame(data=[
            {"ra": 0.04*idx, "dec": 0.04*idx,
             "diaSourceId": idx + 1 + nSources, "diaObjectId": 0,
             "ssObjectId": 0}
            for idx in range(nSources)])

        config = self._makeDefaultConfig(config_file=self.config_file.name, doPackageAlerts=False)
        task = DiaPipelineTask(config=config)
        result = task.createNewDiaObjects(diaSources)
        self.assertEqual(nSources, len(result.newDiaObjects))
        self.assertTrue(np.all(np.equal(
            result.diaSources["diaObjectId"].to_numpy(),
            result.diaSources["diaSourceId"].to_numpy())))
        self.assertTrue(np.all(np.equal(
            result.newDiaObjects["diaObjectId"].to_numpy(),
            result.diaSources["diaSourceId"].to_numpy())))

    def test_purgeDiaObjects(self):
        """Remove diaOjects that are outside an image's bounding box.
        """

        config = self._makeDefaultConfig(config_file=self.config_file.name, doPackageAlerts=False)
        task = DiaPipelineTask(config=config)
        exposure = makeExposure(False, False)
        nObj0 = 20

        # Create diaObjects
        diaObjects = makeDiaObjects(nObj0, exposure, self.rng)
        # Shrink the bounding box so that some of the diaObjects will be outside
        bbox = exposure.getBBox()
        size = np.minimum(bbox.getHeight(), bbox.getWidth())
        bbox.grow(-size//4)
        exposureCut = exposure[bbox]
        sizeCut = np.minimum(bbox.getHeight(), bbox.getWidth())
        buffer = 10
        bbox.grow(buffer)

        def check_diaObjects(bbox, wcs, diaObjects):
            raVals = diaObjects.ra.to_numpy()
            decVals = diaObjects.dec.to_numpy()
            xVals, yVals = wcs.skyToPixelArray(raVals, decVals, degrees=True)
            selector = bbox.contains(xVals, yVals)
            return selector

        selector0 = check_diaObjects(bbox, exposureCut.getWcs(), diaObjects)
        nIn0 = np.count_nonzero(selector0)
        nOut0 = np.count_nonzero(~selector0)
        self.assertEqual(nObj0, nIn0 + nOut0)

        # Add an ID that is not in the diaObject table. It should not get removed.
        diaObjectIds0 = diaObjects["diaObjectId"].copy(deep=True)
        diaObjectIds0[max(diaObjectIds0.index) + 1] = 999
        diaObjects1, objIds = task.purgeDiaObjects(exposureCut.getBBox(), exposureCut.getWcs(), diaObjects,
                                                   diaObjectIds=diaObjectIds0, buffer=buffer)
        diaObjectIds1 = diaObjects1["diaObjectId"]
        # Verify that the bounding box was not changed
        sizeCheck = np.minimum(exposureCut.getBBox().getHeight(), exposureCut.getBBox().getWidth())
        self.assertEqual(sizeCut, sizeCheck)
        selector1 = check_diaObjects(bbox, exposureCut.getWcs(), diaObjects1)
        nIn1 = np.count_nonzero(selector1)
        nOut1 = np.count_nonzero(~selector1)
        nObj1 = len(diaObjects1)
        self.assertEqual(nObj1, nIn0)
        # Verify that not all diaObjects were removed
        self.assertGreater(nObj1, 0)
        # Check that some diaObjects were removed
        self.assertLess(nObj1, nObj0)
        # Verify that no objects outside the bounding box remain
        self.assertEqual(nOut1, 0)
        # Verify that no objects inside the bounding box were removed
        self.assertEqual(nIn1, nIn0)
        # The length of the updated object IDs should equal the number of objects
        # plus one, since we added an extra ID.
        self.assertEqual(nObj1 + 1, len(objIds))
        # All of the object IDs extracted from the catalog should be in the pruned object IDs
        self.assertTrue(set(objIds).issuperset(diaObjectIds1))
        # The pruned object IDs should contain entries that are not in the catalog
        self.assertFalse(set(diaObjectIds1).issuperset(objIds))
        # Some IDs should have been removed
        self.assertLess(len(objIds), len(diaObjectIds0))

    def test_mergeEmptyCatalog(self):
        """Test that a catalog is unchanged if it is merged with an empty
        catalog.
        """
        diaSourcesBase = self.diaSources

        config = self._makeDefaultConfig(config_file=self.config_file.name, doPackageAlerts=False)
        task = DiaPipelineTask(config=config)
        # Include some but not all columns that should be in diaSourcesBase, and some that are mis-matched
        diaSourcesEmpty = pd.DataFrame(columns=["ra", "dec", "foo"])
        diaSourcesTest = task.mergeCatalogs(diaSourcesBase, diaSourcesEmpty, "diaSourcesTest")
        self.assertTrue(diaSourcesBase.equals(diaSourcesTest))

    def test_mergeCatalogs(self):
        """Test that a merged catalog is concatenated correctly.
        """
        diaSourcesBase = self.diaSources
        nBase = len(diaSourcesBase)
        nNew = int(nBase/2)

        diaSourcesNew = makeDiaSources(nNew, self.diaObjects["diaObjectId"].to_numpy(), self.exposure,
                                       self.rng)
        config = self._makeDefaultConfig(config_file=self.config_file.name, doPackageAlerts=False)
        task = DiaPipelineTask(config=config)
        diaSourcesTest = task.mergeCatalogs(diaSourcesBase, diaSourcesNew, "diaSourcesTest")
        self.assertEqual(len(diaSourcesTest), nBase + nNew)
        diaSourcesExtract1 = diaSourcesTest.iloc[:nBase]
        diaSourcesExtract2 = diaSourcesTest.iloc[nBase:]

        self.assertTrue(diaSourcesBase.equals(diaSourcesExtract1))
        self.assertTrue(diaSourcesNew.equals(diaSourcesExtract2))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
