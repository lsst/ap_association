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
import unittest.mock

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

import lsst.utils.tests
from lsst.pex.config import Config
from lsst.daf.base import PropertySet
from lsst.dax.ppdb import Ppdb
from lsst.pipe.base import Task, Struct
from lsst.verify import Name, Measurement
from lsst.verify.gen2tasks.testUtils import MetricTaskTestCase
from lsst.verify.tasks import MetricComputationError
from lsst.verify.tasks.testUtils import MetadataMetricTestCase, PpdbMetricTestCase

from lsst.ap.association.metrics import \
    NumberNewDiaObjectsMetricTask, \
    NumberUnassociatedDiaObjectsMetricTask, \
    FractionUpdatedDiaObjectsMetricTask, \
    TotalUnassociatedDiaObjectsMetricTask


def _makeAssociationMetadata(numUpdated=27, numNew=4, numUnassociated=15):
    metadata = PropertySet()
    metadata.add("association.numUpdatedDiaObjects", numUpdated)
    metadata.add("association.numNewDiaObjects", numNew)
    metadata.add("association.numUnassociatedDiaObjects", numUnassociated)
    return metadata


class TestNewDiaObjects(MetadataMetricTestCase):

    @classmethod
    def makeTask(cls):
        return NumberNewDiaObjectsMetricTask()

    def testValid(self):
        metadata = _makeAssociationMetadata()
        result = self.task.run([metadata])
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.numNewDiaObjects"))
        self.assertEqual(meas.quantity, metadata.getAsDouble("association.numNewDiaObjects") * u.count)

    def testNoNew(self):
        metadata = _makeAssociationMetadata(numNew=0)
        result = self.task.run([metadata])
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.numNewDiaObjects"))
        self.assertEqual(meas.quantity, 0.0 * u.count)

    def testMissingData(self):
        result = self.task.run([None])
        meas = result.measurement
        self.assertIsNone(meas)

    def testNoDataExpected(self):
        result = self.task.run([])
        meas = result.measurement
        self.assertIsNone(meas)

    def testAssociationFailed(self):
        result = self.task.run([PropertySet()])
        meas = result.measurement
        self.assertIsNone(meas)

    def testBadlyTypedKeys(self):
        metadata = _makeAssociationMetadata()
        metadata.set("association.numNewDiaObjects", "Ultimate Answer")

        with self.assertRaises(MetricComputationError):
            self.task.run([metadata])

    def testGetInputDatasetTypes(self):
        config = self.taskClass.ConfigClass()
        config.connections.taskName = "test"
        types = self.taskClass.getInputDatasetTypes(config)
        # dict.keys() is a collections.abc.Set, which has a narrower interface than __builtins__.set...
        self.assertSetEqual(set(types.keys()), {"metadata"})
        self.assertEqual(types["metadata"], "test_metadata")

    def testFineGrainedMetric(self):
        metadata = _makeAssociationMetadata()
        inputData = {"metadata": [metadata]}
        inputDataIds = {"metadata": [{"visit": 42, "ccd": 1}]}
        outputDataId = {"measurement": {"visit": 42, "ccd": 1}}
        measDirect = self.task.run([metadata]).measurement
        measIndirect = self.task.adaptArgsAndRun(inputData, inputDataIds, outputDataId).measurement

        assert_quantity_allclose(measIndirect.quantity, measDirect.quantity)

    def testCoarseGrainedMetric(self):
        metadata = _makeAssociationMetadata()
        nCcds = 3
        inputData = {"metadata": [metadata] * nCcds}
        inputDataIds = {"metadata": [{"visit": 42, "ccd": x} for x in range(nCcds)]}
        outputDataId = {"measurement": {"visit": 42}}
        measDirect = self.task.run([metadata]).measurement
        measMany = self.task.adaptArgsAndRun(inputData, inputDataIds, outputDataId).measurement

        assert_quantity_allclose(measMany.quantity, nCcds * measDirect.quantity)


class TestUnassociatedDiaObjects(MetadataMetricTestCase):

    @classmethod
    def makeTask(cls):
        return NumberUnassociatedDiaObjectsMetricTask()

    def testValid(self):
        metadata = _makeAssociationMetadata()
        result = self.task.run([metadata])
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.numUnassociatedDiaObjects"))
        self.assertEqual(meas.quantity,
                         metadata.getAsDouble("association.numUnassociatedDiaObjects") * u.count)

    def testAllUpdated(self):
        metadata = _makeAssociationMetadata(numUnassociated=0)
        result = self.task.run([metadata])
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.numUnassociatedDiaObjects"))
        self.assertEqual(meas.quantity, 0.0 * u.count)

    def testMissingData(self):
        result = self.task.run([None])
        meas = result.measurement
        self.assertIsNone(meas)

    def testNoDataExpected(self):
        result = self.task.run([])
        meas = result.measurement
        self.assertIsNone(meas)

    def testAssociationFailed(self):
        result = self.task.run([PropertySet()])
        meas = result.measurement
        self.assertIsNone(meas)

    def testBadlyTypedKeys(self):
        metadata = _makeAssociationMetadata()
        metadata.set("association.numUnassociatedDiaObjects", "Ultimate Answer")

        with self.assertRaises(MetricComputationError):
            self.task.run([metadata])

    def testGetInputDatasetTypes(self):
        config = self.taskClass.ConfigClass()
        config.connections.taskName = "test"
        types = self.taskClass.getInputDatasetTypes(config)
        # dict.keys() is a collections.abc.Set, which has a narrower interface than __builtins__.set...
        self.assertSetEqual(set(types.keys()), {"metadata"})
        self.assertEqual(types["metadata"], "test_metadata")

    def testFineGrainedMetric(self):
        metadata = _makeAssociationMetadata()
        inputData = {"metadata": [metadata]}
        inputDataIds = {"metadata": [{"visit": 42, "ccd": 1}]}
        outputDataId = {"measurement": {"visit": 42, "ccd": 1}}
        measDirect = self.task.run([metadata]).measurement
        measIndirect = self.task.adaptArgsAndRun(inputData, inputDataIds, outputDataId).measurement

        assert_quantity_allclose(measIndirect.quantity, measDirect.quantity)

    def testCoarseGrainedMetric(self):
        metadata = _makeAssociationMetadata()
        nCcds = 3
        inputData = {"metadata": [metadata] * nCcds}
        inputDataIds = {"metadata": [{"visit": 42, "ccd": x} for x in range(nCcds)]}
        outputDataId = {"measurement": {"visit": 42}}
        measDirect = self.task.run([metadata]).measurement
        measMany = self.task.adaptArgsAndRun(inputData, inputDataIds, outputDataId).measurement

        assert_quantity_allclose(measMany.quantity, nCcds * measDirect.quantity)


class TestFracUpdatedDiaObjects(MetadataMetricTestCase):

    @classmethod
    def makeTask(cls):
        return FractionUpdatedDiaObjectsMetricTask()

    def testValid(self):
        metadata = _makeAssociationMetadata()
        result = self.task.run([metadata])
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.fracUpdatedDiaObjects"))
        nUpdated = metadata.getAsDouble("association.numUpdatedDiaObjects")
        nTotal = metadata.getAsDouble("association.numUnassociatedDiaObjects") + nUpdated
        self.assertEqual(meas.quantity, nUpdated / nTotal * u.dimensionless_unscaled)

    def testNoUpdated(self):
        metadata = _makeAssociationMetadata(numUpdated=0)
        result = self.task.run([metadata])
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.fracUpdatedDiaObjects"))
        self.assertEqual(meas.quantity, 0.0 * u.dimensionless_unscaled)

    def testAllUpdated(self):
        metadata = _makeAssociationMetadata(numUnassociated=0)
        result = self.task.run([metadata])
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.fracUpdatedDiaObjects"))
        self.assertEqual(meas.quantity, 1.0 * u.dimensionless_unscaled)

    def testNoObjects(self):
        metadata = _makeAssociationMetadata(numUpdated=0, numUnassociated=0)
        with self.assertRaises(MetricComputationError):
            self.task.run([metadata])

    def testMissingData(self):
        result = self.task.run([None])
        meas = result.measurement
        self.assertIsNone(meas)

    def testNoDataExpected(self):
        result = self.task.run([])
        meas = result.measurement
        self.assertIsNone(meas)

    def testAssociationFailed(self):
        result = self.task.run([PropertySet()])
        meas = result.measurement
        self.assertIsNone(meas)

    def testBadlyTypedKeys(self):
        metadata = _makeAssociationMetadata()
        metadata.set("association.numUnassociatedDiaObjects", "Ultimate Answer")

        with self.assertRaises(MetricComputationError):
            self.task.run([metadata])

    def testGetInputDatasetTypes(self):
        config = self.taskClass.ConfigClass()
        config.connections.taskName = "test"
        types = self.taskClass.getInputDatasetTypes(config)
        # dict.keys() is a collections.abc.Set, which has a narrower interface than __builtins__.set...
        self.assertSetEqual(set(types.keys()), {"metadata"})
        self.assertEqual(types["metadata"], "test_metadata")

    def testFineGrainedMetric(self):
        metadata = _makeAssociationMetadata()
        inputData = {"metadata": [metadata]}
        inputDataIds = {"metadata": [{"visit": 42, "ccd": 1}]}
        outputDataId = {"measurement": {"visit": 42, "ccd": 1}}
        measDirect = self.task.run([metadata]).measurement
        measIndirect = self.task.adaptArgsAndRun(inputData, inputDataIds, outputDataId).measurement

        assert_quantity_allclose(measIndirect.quantity, measDirect.quantity)

    def testCoarseGrainedMetric(self):
        metadata = _makeAssociationMetadata()
        nCcds = 3
        inputData = {"metadata": [metadata] * nCcds}
        inputDataIds = {"metadata": [{"visit": 42, "ccd": x} for x in range(nCcds)]}
        outputDataId = {"measurement": {"visit": 42}}
        measDirect = self.task.run([metadata]).measurement
        measMany = self.task.adaptArgsAndRun(inputData, inputDataIds, outputDataId).measurement

        assert_quantity_allclose(measMany.quantity, measDirect.quantity)


class TestTotalUnassociatedObjects(PpdbMetricTestCase):

    @staticmethod
    def _makePpdb(dummy_dbInfo):
        """Create a dummy ppdb.

        We don't have access to the ppdb in the task directly so mocking
        return values is difficult. We thus make use of the dummy dbInfo
        that is passed to the init task to pass values to the ppdb object
        instantiated.
        """
        ppdb = unittest.mock.Mock(Ppdb)
        try:
            test_value = dummy_dbInfo["test_value"]
            ppdb.countUnassociatedObjects = unittest.mock.MagicMock(
                return_value=test_value)
        except TypeError:
            ppdb.countUnassociatedObjects = unittest.mock.MagicMock(
                return_value=42)
        return ppdb

    @classmethod
    def makeTask(cls):
        class SimpleDbLoader(Task):
            ConfigClass = Config

            def run(self, dummy):
                if dummy is not None:
                    return Struct(ppdb=cls._makePpdb(dummy))
                else:
                    return Struct(ppdb=None)

        config = TotalUnassociatedDiaObjectsMetricTask.ConfigClass()
        config.dbLoader.retarget(SimpleDbLoader)
        return TotalUnassociatedDiaObjectsMetricTask(config=config)

    def setUp(self):
        super().setUp()
        # Do the patch here to avoid passing extra arguments to superclass tests

    def testValid(self):
        result = self.task.adaptArgsAndRun({"dbInfo": [{"test_value": 42}]},
                                           {"dbInfo": {}},
                                           {"measurement": {}})
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.totalUnassociatedDiaObjects"))
        self.assertEqual(meas.quantity, 42 * u.count)

    def testAllAssociated(self):
        result = self.task.adaptArgsAndRun({"dbInfo": [{"test_value": 0}]},
                                           {"dbInfo": {}},
                                           {"measurement": {}})
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.totalUnassociatedDiaObjects"))
        self.assertEqual(meas.quantity, 0.0 * u.count)

    def testMissingData(self):
        result = self.task.adaptArgsAndRun({"dbInfo": None}, {"dbInfo": {}}, {"measurement": {}})
        meas = result.measurement
        self.assertIsNone(meas)

    def testFineGrainedMetric(self):
        with self.assertRaises(ValueError):
            self.task.adaptArgsAndRun({"dbInfo": "DB source"}, {"dbInfo": {}}, {"measurement": {"visit": 42}})

    def testGetInputDatasetTypes(self):
        config = self.taskClass.ConfigClass()
        config.connections.dbInfo = "absolutely anything"
        types = self.taskClass.getInputDatasetTypes(config)
        # dict.keys() is a collections.abc.Set, which has a narrower interface than __builtins__.set...
        self.assertSetEqual(set(types.keys()), {"dbInfo"})
        self.assertEqual(types["dbInfo"], "absolutely anything")

    # Override because of non-standard adaptArgsAndRun
    def testCallAddStandardMetadata(self):
        dummy = Measurement('foo.bar', 0.0)
        with unittest.mock.patch.multiple(
                self.taskClass, autospec=True,
                makeMeasurement=unittest.mock.DEFAULT,
                addStandardMetadata=unittest.mock.DEFAULT) as mockDict:
            mockDict['makeMeasurement'].return_value = Struct(measurement=dummy)

            dataId = {}
            result = self.task.adaptArgsAndRun(
                {"dbInfo": "DB Source"},
                {"dbInfo": dataId},
                {'measurement': dataId})
            mockDict['addStandardMetadata'].assert_called_once_with(
                self.task, result.measurement, dataId)


# Hack around unittest's hacky test setup system
del MetricTaskTestCase
del MetadataMetricTestCase
del PpdbMetricTestCase


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
