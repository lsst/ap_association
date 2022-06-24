# This file is part of ap_association
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

import os
import unittest

import numpy as np

from lsst.ap.association.transformDiaSourceCatalog import (TransformDiaSourceCatalogConfig,
                                                           TransformDiaSourceCatalogTask)
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.meas.base.tests as measTests
from lsst.pipe.base import Struct
from lsst.utils import getPackageDir
import lsst.utils.tests

from lsst.ap.association.transformDiaSourceCatalog import UnpackApdbFlags


class TestTransformDiaSourceCatalogTask(unittest.TestCase):
    def setUp(self):
        self.nSources = 10
        self.bboxSize = 18
        self.xyLoc = 100
        self.bbox = geom.Box2I(geom.Point2I(0, 0),
                               geom.Extent2I(1024, 1153))
        dataset = measTests.TestDataset(self.bbox)
        for srcIdx in range(self.nSources-1):
            dataset.addSource(100000.0, geom.Point2D(self.xyLoc, self.xyLoc))
        # Ensure the last source has no peak `significance` field.
        dataset.addSource(100000.0, geom.Point2D(self.xyLoc, self.xyLoc), setPeakSignificance=False)
        schema = dataset.makeMinimalSchema()
        schema.addField("base_PixelFlags_flag", type="Flag")
        schema.addField("base_PixelFlags_flag_offimage", type="Flag")
        self.exposure, self.inputCatalog = dataset.realize(10.0, schema, randomSeed=1234)
        # Create schemas for use in initializing the TransformDiaSourceCatalog
        # task.
        self.initInputs = {"diaSourceSchema": Struct(schema=schema)}
        self.initInputsBadFlags = {"diaSourceSchema": Struct(schema=dataset.makeMinimalSchema())}

        self.expId = 4321
        self.date = dafBase.DateTime(nsecs=1400000000 * 10**9)
        detector = DetectorWrapper(id=23, bbox=self.exposure.getBBox()).detector
        visit = afwImage.VisitInfo(
            exposureId=self.expId,
            exposureTime=200.,
            date=self.date)
        self.exposure.info.id = self.expId
        self.exposure.setDetector(detector)
        self.exposure.getInfo().setVisitInfo(visit)
        self.filterName = 'g'
        self.exposure.setFilter(afwImage.FilterLabel(band=self.filterName, physical='g.MP9401'))
        scale = 2
        scaleErr = 1
        self.photoCalib = afwImage.PhotoCalib(scale, scaleErr)
        self.exposure.setPhotoCalib(self.photoCalib)

    def test_run(self):
        """Test output dataFrame is created and values are correctly inserted
        from the exposure.
        """
        transConfig = TransformDiaSourceCatalogConfig()
        transConfig.flagMap = os.path.join(
            getPackageDir("ap_association"),
            "tests",
            "data",
            "test-flag-map.yaml")
        transConfig.functorFile = os.path.join(getPackageDir("ap_association"),
                                               "tests",
                                               "data",
                                               "testDiaSource.yaml")
        transformTask = TransformDiaSourceCatalogTask(initInputs=self.initInputs,
                                                      config=transConfig)
        result = transformTask.run(self.inputCatalog,
                                   self.exposure,
                                   self.filterName,
                                   ccdVisitId=self.expId)

        self.assertEqual(len(result.diaSourceTable), len(self.inputCatalog))
        np.testing.assert_array_equal(result.diaSourceTable["bboxSize"], [self.bboxSize]*self.nSources)
        np.testing.assert_array_equal(result.diaSourceTable["ccdVisitId"], [self.expId]*self.nSources)
        np.testing.assert_array_equal(result.diaSourceTable["filterName"], [self.filterName]*self.nSources)
        np.testing.assert_array_equal(result.diaSourceTable["midPointTai"],
                                      [self.date.get(system=dafBase.DateTime.MJD)]*self.nSources)
        np.testing.assert_array_equal(result.diaSourceTable["diaObjectId"], [0]*self.nSources)
        # The final snr value should be NaN because it doesn't have a peak significance field.
        expect_snr = [397.887353515625]*9
        expect_snr.append(np.nan)
        self.assertTrue(np.array_equal(result.diaSourceTable["snr"], expect_snr, equal_nan=True))

    def test_run_dia_source_wrong_flags(self):
        """Test that the proper errors are thrown when requesting flag columns
        that are not in the input schema.
        """
        with self.assertRaises(KeyError):
            TransformDiaSourceCatalogTask(initInputs=self.initInputsBadFlags)

    def test_computeBBoxSize(self):
        """Test the values created for diaSourceBBox.
        """
        transConfig = TransformDiaSourceCatalogConfig()
        transConfig.flagMap = os.path.join(
            getPackageDir("ap_association"),
            "tests",
            "data",
            "test-flag-map.yaml")
        transform = TransformDiaSourceCatalogTask(initInputs=self.initInputs,
                                                  config=transConfig)
        bboxArray = transform.computeBBoxSizes(self.inputCatalog)

        # Default in catalog is 18.
        self.assertEqual(bboxArray[0], self.bboxSize)

    def test_bit_unpacker(self):
        """Test that the integer bit packer is functioning correctly.
        """
        transConfig = TransformDiaSourceCatalogConfig()
        transConfig.flagMap = os.path.join(
            getPackageDir("ap_association"),
            "tests",
            "data",
            "test-flag-map.yaml")
        transConfig.functorFile = os.path.join(getPackageDir("ap_association"),
                                               "tests",
                                               "data",
                                               "testDiaSource.yaml")
        transform = TransformDiaSourceCatalogTask(initInputs=self.initInputs,
                                                  config=transConfig)
        for idx, obj in enumerate(self.inputCatalog):
            if idx in [1, 3, 5]:
                obj.set("base_PixelFlags_flag", 1)
            if idx in [1, 4, 6]:
                obj.set("base_PixelFlags_flag_offimage", 1)
        outputCatalog = transform.run(self.inputCatalog,
                                      self.exposure,
                                      self.filterName,
                                      ccdVisitId=self.expId).diaSourceTable

        unpacker = UnpackApdbFlags(transConfig.flagMap, "DiaSource")
        flag_values = unpacker.unpack(outputCatalog["flags"], "flags")

        for idx, flag in enumerate(flag_values):
            if idx in [1, 3, 5]:
                self.assertTrue(flag['base_PixelFlags_flag'])
            else:
                self.assertFalse(flag['base_PixelFlags_flag'])

            if idx in [1, 4, 6]:
                self.assertTrue(flag['base_PixelFlags_flag_offimage'])
            else:
                self.assertFalse(flag['base_PixelFlags_flag_offimage'])


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
