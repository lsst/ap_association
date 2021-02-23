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

from lsst.ap.association.transformDiaSourceCatalog import (TransformDiaSourceCatalogConfig,
                                                           TransformDiaSourceCatalogTask)
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.meas.base.tests as measTests
from lsst.utils import getPackageDir
import lsst.utils.tests


class TestTransformDiaSourceCatalogTask(unittest.TestCase):

    def setUp(self):
        nSources = 10
        self.bboxSize = 18
        self.xyLoc = 100
        self.bbox = geom.Box2I(geom.Point2I(0, 0),
                               geom.Extent2I(1024, 1153))
        dataset = measTests.TestDataset(self.bbox)
        for srcIdx in range(nSources):
            dataset.addSource(100000.0, geom.Point2D(self.xyLoc, self.xyLoc))
        schema = dataset.makeMinimalSchema()
        self.exposure, self.inputCatalog = dataset.realize(10.0, schema, randomSeed=1234)

        self.expId = 4321
        self.date = dafBase.DateTime(nsecs=1400000000 * 10**9)
        detector = DetectorWrapper(id=23, bbox=self.exposure.getBBox()).detector
        visit = afwImage.VisitInfo(
            exposureId=self.expId,
            exposureTime=200.,
            date=self.date)
        self.exposure.setDetector(detector)
        self.exposure.getInfo().setVisitInfo(visit)
        self.filterName = 'g'
        self.exposure.setFilterLabel(afwImage.FilterLabel(band=self.filterName, physical='g.MP9401'))
        scale = 2
        scaleErr = 1
        self.photoCalib = afwImage.PhotoCalib(scale, scaleErr)
        self.exposure.setPhotoCalib(self.photoCalib)

    def test_run(self):
        """Test output dataFrame is created and values are correctly inserted
        from the exposure.
        """
        transformConfig = TransformDiaSourceCatalogConfig()
        transformConfig.functorFile = os.path.join(getPackageDir("ap_association"),
                                                   "tests/data/",
                                                   "testDiaSource.yaml")
        transformTask = TransformDiaSourceCatalogTask(config=transformConfig)
        result = transformTask.run(self.inputCatalog,
                                   self.exposure,
                                   self.filterName,
                                   ccdVisitId=self.expId)

        self.assertEqual(len(result.diaSourceTable), len(self.inputCatalog))
        for idx, src in result.diaSourceTable.iterrows():
            self.assertEqual(src["bboxSize"], self.bboxSize)
            self.assertEqual(src["ccdVisitId"], self.expId)
            self.assertEqual(src["filterName"], self.filterName)
            self.assertEqual(src["midPointTai"],
                             self.date.get(system=dafBase.DateTime.MJD))
            self.assertEqual(src["pixelId"], 0)
            self.assertEqual(src["diaObjectId"], 0)

    def test_computeBBoxSize(self):
        """Test the values created for diaSourceBBox.
        """
        transform = TransformDiaSourceCatalogTask()
        bboxArray = transform.computeBBoxSizes(self.inputCatalog)

        # Default in catalog is 18.
        self.assertEqual(bboxArray[0], self.bboxSize)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
