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

import lsst.utils.tests
from lsst.pex.config import Config
from lsst.daf.base import PropertySet
from lsst.dax.apdb import Apdb
from lsst.pipe.base import Task, Struct
import lsst.pipe.base.testUtils
from lsst.verify import Name
from lsst.verify.tasks import MetricComputationError
from lsst.verify.tasks.testUtils import MetricTaskTestCase, MetadataMetricTestCase, ApdbMetricTestCase

from lsst.ap.association.metrics import \
    NumberNewDiaObjectsMetricTask, \
    NumberUnassociatedDiaObjectsMetricTask, \
    FractionUpdatedDiaObjectsMetricTask, \
    NumberSolarSystemObjectsMetricTask, \
    NumberAssociatedSolarSystemObjectsMetricTask, \
    TotalUnassociatedDiaObjectsMetricTask


def _makeAssociationMetadata(numUpdated=27, numNew=4, numUnassociated=15, numSso=5, numAssocSso=1):
    metadata = PropertySet()
    metadata.add("association.numUpdatedDiaObjects", numUpdated)
    metadata.add("association.numNewDiaObjects", numNew)
    metadata.add("association.numTotalSolarSystemObjects", numSso)
    metadata.add("association.numAssociatedSsObjects", numAssocSso)
    metadata.add("association.numUnassociatedDiaObjects", numUnassociated)
    return metadata


class TestNewDiaObjects(MetadataMetricTestCase):

    @classmethod
    def makeTask(cls):
        return NumberNewDiaObjectsMetricTask()

    def testValid(self):
        metadata = _makeAssociationMetadata()
        result = self.task.run(metadata)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.numNewDiaObjects"))
        self.assertEqual(meas.quantity, metadata.getAsDouble("association.numNewDiaObjects") * u.count)

    def testNoNew(self):
        metadata = _makeAssociationMetadata(numNew=0)
        result = self.task.run(metadata)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.numNewDiaObjects"))
        self.assertEqual(meas.quantity, 0.0 * u.count)

    def testMissingData(self):
        result = self.task.run(None)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement
        self.assertIsNone(meas)

    def testAssociationFailed(self):
        result = self.task.run(PropertySet())
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement
        self.assertIsNone(meas)

    def testBadlyTypedKeys(self):
        metadata = _makeAssociationMetadata()
        metadata.set("association.numNewDiaObjects", "Ultimate Answer")

        with self.assertRaises(MetricComputationError):
            self.task.run(metadata)


class TestUnassociatedDiaObjects(MetadataMetricTestCase):

    @classmethod
    def makeTask(cls):
        return NumberUnassociatedDiaObjectsMetricTask()

    def testValid(self):
        metadata = _makeAssociationMetadata()
        result = self.task.run(metadata)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.numUnassociatedDiaObjects"))
        self.assertEqual(meas.quantity,
                         metadata.getAsDouble("association.numUnassociatedDiaObjects") * u.count)

    def testAllUpdated(self):
        metadata = _makeAssociationMetadata(numUnassociated=0)
        result = self.task.run(metadata)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.numUnassociatedDiaObjects"))
        self.assertEqual(meas.quantity, 0.0 * u.count)

    def testMissingData(self):
        result = self.task.run(None)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement
        self.assertIsNone(meas)

    def testAssociationFailed(self):
        result = self.task.run(PropertySet())
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement
        self.assertIsNone(meas)

    def testBadlyTypedKeys(self):
        metadata = _makeAssociationMetadata()
        metadata.set("association.numUnassociatedDiaObjects", "Ultimate Answer")

        with self.assertRaises(MetricComputationError):
            self.task.run(metadata)


class TestFracUpdatedDiaObjects(MetadataMetricTestCase):

    @classmethod
    def makeTask(cls):
        return FractionUpdatedDiaObjectsMetricTask()

    def testValid(self):
        metadata = _makeAssociationMetadata()
        result = self.task.run(metadata)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.fracUpdatedDiaObjects"))
        nUpdated = metadata.getAsDouble("association.numUpdatedDiaObjects")
        nTotal = metadata.getAsDouble("association.numUnassociatedDiaObjects") + nUpdated
        self.assertEqual(meas.quantity, nUpdated / nTotal * u.dimensionless_unscaled)

    def testNoUpdated(self):
        metadata = _makeAssociationMetadata(numUpdated=0)
        result = self.task.run(metadata)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.fracUpdatedDiaObjects"))
        self.assertEqual(meas.quantity, 0.0 * u.dimensionless_unscaled)

    def testAllUpdated(self):
        metadata = _makeAssociationMetadata(numUnassociated=0)
        result = self.task.run(metadata)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.fracUpdatedDiaObjects"))
        self.assertEqual(meas.quantity, 1.0 * u.dimensionless_unscaled)

    def testNoObjects(self):
        metadata = _makeAssociationMetadata(numUpdated=0, numUnassociated=0)
        result = self.task.run(metadata)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertIsNone(meas)

    def testMissingData(self):
        result = self.task.run(None)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement
        self.assertIsNone(meas)

    def testAssociationFailed(self):
        result = self.task.run(PropertySet())
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement
        self.assertIsNone(meas)

    def testBadlyTypedKeys(self):
        metadata = _makeAssociationMetadata()
        metadata.set("association.numUnassociatedDiaObjects", "Ultimate Answer")

        with self.assertRaises(MetricComputationError):
            self.task.run(metadata)


class TestNumberSolarSystemObjects(MetadataMetricTestCase):

    @classmethod
    def makeTask(cls):
        return NumberSolarSystemObjectsMetricTask()

    def testValid(self):
        metadata = _makeAssociationMetadata()
        result = self.task.run(metadata)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.numTotalSolarSystemObjects"))
        self.assertEqual(meas.quantity,
                         metadata.getAsDouble("association.numTotalSolarSystemObjects") * u.count)

    def testAllUpdated(self):
        metadata = _makeAssociationMetadata(numSso=0)
        result = self.task.run(metadata)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.numTotalSolarSystemObjects"))
        self.assertEqual(meas.quantity, 0.0 * u.count)

    def testMissingData(self):
        result = self.task.run(None)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement
        self.assertIsNone(meas)

    def testAssociationFailed(self):
        result = self.task.run(PropertySet())
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement
        self.assertIsNone(meas)

    def testBadlyTypedKeys(self):
        metadata = _makeAssociationMetadata()
        metadata.set("association.numTotalSolarSystemObjects", "Ultimate Answer")

        with self.assertRaises(MetricComputationError):
            self.task.run(metadata)


class TestNumberAssociatedSolarSystemObjects(MetadataMetricTestCase):

    @classmethod
    def makeTask(cls):
        return NumberAssociatedSolarSystemObjectsMetricTask()

    def testValid(self):
        metadata = _makeAssociationMetadata()
        result = self.task.run(metadata)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.numAssociatedSsObjects"))
        self.assertEqual(meas.quantity,
                         metadata.getAsDouble("association.numAssociatedSsObjects") * u.count)

    def testAllUpdated(self):
        metadata = _makeAssociationMetadata(numAssocSso=0)
        result = self.task.run(metadata)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.numAssociatedSsObjects"))
        self.assertEqual(meas.quantity, 0.0 * u.count)

    def testMissingData(self):
        result = self.task.run(None)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement
        self.assertIsNone(meas)

    def testAssociationFailed(self):
        result = self.task.run(PropertySet())
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement
        self.assertIsNone(meas)

    def testBadlyTypedKeys(self):
        metadata = _makeAssociationMetadata()
        metadata.set("association.numAssociatedSsObjects", "Ultimate Answer")

        with self.assertRaises(MetricComputationError):
            self.task.run(metadata)


class TestTotalUnassociatedObjects(ApdbMetricTestCase):

    @staticmethod
    def _makeApdb(dummy_dbInfo):
        """Create a dummy apdb.

        We don't have access to the apdb in the task directly so mocking
        return values is difficult. We thus make use of the dummy dbInfo
        that is passed to the init task to pass values to the apdb object
        instantiated.
        """
        apdb = unittest.mock.Mock(Apdb)
        test_value = dummy_dbInfo["test_value"]
        apdb.countUnassociatedObjects = unittest.mock.MagicMock(
            return_value=test_value)
        return apdb

    @classmethod
    def makeTask(cls):
        class SimpleDbLoader(Task):
            ConfigClass = Config

            def run(self, dummy):
                if dummy is not None:
                    return Struct(apdb=cls._makeApdb(dummy))
                else:
                    return Struct(apdb=None)

        config = TotalUnassociatedDiaObjectsMetricTask.ConfigClass()
        config.dbLoader.retarget(SimpleDbLoader)
        return TotalUnassociatedDiaObjectsMetricTask(config=config)

    @classmethod
    def makeDbInfo(cls):
        return {"test_value": "whatever"}

    def setUp(self):
        super().setUp()

    def testValid(self):
        result = self.task.run([{"test_value": 42}])
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.totalUnassociatedDiaObjects"))
        self.assertEqual(meas.quantity, 42 * u.count)

    def testAllAssociated(self):
        result = self.task.run([{"test_value": 0}])
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ap_association.totalUnassociatedDiaObjects"))
        self.assertEqual(meas.quantity, 0.0 * u.count)

    def testMissingData(self):
        result = self.task.run(None)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement
        self.assertIsNone(meas)

    def testFineGrainedMetric(self):
        with self.assertRaises(ValueError):
            self.task.run([self.makeDbInfo()], outputDataId={"visit": 42})


# Hack around unittest's hacky test setup system
del MetricTaskTestCase
del MetadataMetricTestCase
del ApdbMetricTestCase


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
