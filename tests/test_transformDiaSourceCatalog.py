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
import lsst.utils.tests

from lsst.ap.association.transformDiaSourceCatalog import UnpackApdbFlags

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class TestTransformDiaSourceCatalogTask(unittest.TestCase):
    def setUp(self):
        # Create an instance of random generator with fixed seed.
        rng = np.random.default_rng(1234)

        self.nSources = 10
        # Default PSF size (psfDim in makeEmptyExposure) in TestDataset results
        # in an 18 pixel wide source box.
        self.bboxSize = 18
        self.yLoc = 100
        self.bbox = geom.Box2I(geom.Point2I(0, 0),
                               geom.Extent2I(1024, 1153))
        dataset = measTests.TestDataset(self.bbox)
        for srcIdx in range(self.nSources-1):
            # Place sources at (index, yLoc), so we can distinguish them later.
            dataset.addSource(100000.0, geom.Point2D(srcIdx, self.yLoc))
        # Ensure the last source has no peak `significance` field.
        dataset.addSource(100000.0, geom.Point2D(srcIdx+1, self.yLoc), setPeakSignificance=False)
        schema = dataset.makeMinimalSchema()
        schema.addField("base_PixelFlags_flag", type="Flag")
        schema.addField("base_PixelFlags_flag_offimage", type="Flag")
        schema.addField("sky_source", type="Flag", doc="Sky objects.")
        schema.addField("reliability", doc="real/bogus score of this source", type=float)
        self.exposure, self.inputCatalog = dataset.realize(10.0, schema, randomSeed=1234)
        # Create schemas for use in initializing the TransformDiaSourceCatalog task.
        self.initInputs = {"diaSourceSchema": Struct(schema=schema)}
        self.initInputsBadFlags = {"diaSourceSchema": Struct(schema=dataset.makeMinimalSchema())}

        self.inputCatalog["reliability"] = rng.random(self.nSources, dtype=np.float32)

        self.expId = 4321
        self.date = dafBase.DateTime(nsecs=1400000000 * 10**9)
        detector = DetectorWrapper(id=23, bbox=self.exposure.getBBox()).detector
        visit = afwImage.VisitInfo(
            id=self.expId,
            exposureTime=200.,
            date=self.date)
        self.exposure.info.id = self.expId + (10000*detector.getId())
        self.exposure.setDetector(detector)
        self.exposure.info.setVisitInfo(visit)
        self.band = 'g'
        self.exposure.setFilter(afwImage.FilterLabel(band=self.band, physical='g.MP9401'))
        scale = 2
        scaleErr = 1
        self.photoCalib = afwImage.PhotoCalib(scale, scaleErr)
        self.exposure.setPhotoCalib(self.photoCalib)

        self.config = TransformDiaSourceCatalogConfig()
        self.config.flagMap = os.path.join(TESTDIR, "data", "test-flag-map.yaml")
        self.config.functorFile = os.path.join(TESTDIR,
                                               "data",
                                               "testDiaSource.yaml")

    def test_run(self):
        """Test output dataFrame is created and values are correctly inserted
        from the exposure.
        """
        transformTask = TransformDiaSourceCatalogTask(initInputs=self.initInputs,
                                                      config=self.config)
        result = transformTask.run(self.inputCatalog,
                                   self.exposure,
                                   self.band)

        self.assertEqual(len(result.diaSourceTable), len(self.inputCatalog))
        np.testing.assert_array_equal(result.diaSourceTable["bboxSize"], [self.bboxSize]*self.nSources)
        np.testing.assert_array_equal(result.diaSourceTable["visit"],
                                      [self.exposure.visitInfo.id]*self.nSources)
        np.testing.assert_array_equal(result.diaSourceTable["detector"],
                                      [self.exposure.detector.getId()]*self.nSources)
        np.testing.assert_array_equal(result.diaSourceTable["band"], [self.band]*self.nSources)
        np.testing.assert_array_equal(result.diaSourceTable["midpointMjdTai"],
                                      [self.date.get(system=dafBase.DateTime.MJD)]*self.nSources)
        np.testing.assert_array_equal(result.diaSourceTable["diaObjectId"], [0]*self.nSources)
        np.testing.assert_array_equal(result.diaSourceTable["x"], np.arange(self.nSources))
        # The final snr value should be NaN because it doesn't have a peak significance field.
        expect_snr = [397.887353515625]*9
        expect_snr.append(np.nan)
        # Have to use allclose because assert_array_equal doesn't support equal_nan.
        np.testing.assert_allclose(result.diaSourceTable["snr"], expect_snr, equal_nan=True, rtol=0)

    def test_run_with_apdb_schema(self):
        """Test output dataFrame is created and values are correctly inserted
        from the exposure.
        """
        self.config.doUseApdbSchema = True
        transformTask = TransformDiaSourceCatalogTask(initInputs=self.initInputs,
                                                      config=self.config)
        result = transformTask.run(self.inputCatalog,
                                   self.exposure,
                                   self.band)

        self.assertEqual(len(result.diaSourceTable), len(self.inputCatalog))
        np.testing.assert_array_equal(result.diaSourceTable["bboxSize"], [self.bboxSize]*self.nSources)
        np.testing.assert_array_equal(result.diaSourceTable["visit"],
                                      [self.exposure.visitInfo.id]*self.nSources)
        np.testing.assert_array_equal(result.diaSourceTable["detector"],
                                      [self.exposure.detector.getId()]*self.nSources)
        np.testing.assert_array_equal(result.diaSourceTable["band"], [self.band]*self.nSources)
        np.testing.assert_array_equal(result.diaSourceTable["midpointMjdTai"],
                                      [self.date.get(system=dafBase.DateTime.MJD)]*self.nSources)
        np.testing.assert_array_equal(result.diaSourceTable["diaObjectId"], [0]*self.nSources)
        np.testing.assert_array_equal(result.diaSourceTable["x"], np.arange(self.nSources))
        # The final snr value should be NaN because it doesn't have a peak significance field.
        expect_snr = [397.887353515625]*9
        expect_snr.append(np.nan)
        # Downcast to single precision to match `sdm_schemas`
        expect_snr = np.array(expect_snr, dtype=np.float32)
        # Have to use allclose because assert_array_equal doesn't support equal_nan.
        np.testing.assert_allclose(result.diaSourceTable["snr"], expect_snr, equal_nan=True, rtol=0)

    def test_run_dia_source_wrong_flags(self):
        """Test that the proper errors are thrown when requesting flag columns
        that are not in the input schema.
        """
        with self.assertRaises(KeyError):
            TransformDiaSourceCatalogTask(initInputs=self.initInputsBadFlags)

    def test_computeBBoxSize(self):
        transform = TransformDiaSourceCatalogTask(initInputs=self.initInputs,
                                                  config=self.config)
        boxSizes = transform.computeBBoxSizes(self.inputCatalog)

        for size in boxSizes:
            self.assertEqual(size, self.bboxSize)
        self.assertEqual(len(boxSizes), self.nSources)

    def test_flag_existence_check(self):
        unpacker = UnpackApdbFlags(self.config.flagMap, "DiaSource")

        self.assertTrue(unpacker.flagExists('base_PixelFlags_flag'))
        self.assertFalse(unpacker.flagExists(''))
        with self.assertRaisesRegex(ValueError, 'column doesNotExist not in flag map'):
            unpacker.flagExists('base_PixelFlags_flag', columnName='doesNotExist')

    def test_flag_bitmask(self):
        """Test that we get the expected bitmask back from supplied flag names.
        """
        unpacker = UnpackApdbFlags(self.config.flagMap, "DiaSource")

        with self.assertRaisesRegex(ValueError, "flag '' not included"):
            unpacker.makeFlagBitMask([''])
        with self.assertRaisesRegex(ValueError, 'column doesNotExist not in flag map'):
            unpacker.makeFlagBitMask(['base_PixelFlags_flag'], columnName='doesNotExist')
        self.assertEqual(unpacker.makeFlagBitMask(['base_PixelFlags_flag']), np.uint64(1))
        self.assertEqual(unpacker.makeFlagBitMask(['base_PixelFlags_flag_offimage']), np.uint64(4))
        self.assertEqual(unpacker.makeFlagBitMask(['base_PixelFlags_flag',
                                                   'base_PixelFlags_flag_offimage']),
                         np.uint64(5))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
