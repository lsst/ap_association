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

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

import lsst.utils.tests
from lsst.daf.base import PropertySet
from lsst.verify import Name
from lsst.verify.gen2tasks.testUtils import MetricTaskTestCase
from lsst.verify.tasks import MetricComputationError
from lsst.verify.tasks.testUtils import MetadataMetricTestCase

from lsst.ap.association.metrics import \
    NumberNewDiaObjectsMetricTask, \
    NumberUnassociatedDiaObjectsMetricTask, \
    FractionUpdatedDiaObjectsMetricTask


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
        config.metadata.name = "test_metadata"
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
        config.metadata.name = "test_metadata"
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
        config.metadata.name = "test_metadata"
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


# Hack around unittest's hacky test setup system
del MetricTaskTestCase
del MetadataMetricTestCase


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
