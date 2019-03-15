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

import numpy as np
import unittest

from lsst.afw.cameraGeom.testUtils import DetectorWrapper
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.image.utils as afwImageUtils
import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
import lsst.meas.algorithms as measAlg
from lsst.ap.association import \
    DiaForcedSourceTask, \
    make_dia_object_schema
import lsst.utils.tests


def create_test_dia_objects(n_points, wcs, startPos=100):
    """Create dummy DIASources or DIAObjects for use in our tests.
    Parameters
    ----------
    n_points : `int`
        Number of DiaObject test points to create.
    wcs : `lsst.afw.geom.SkyWcs`
        Wcs to convert RA/Dec to pixel x/y.
    startPos : `int`
        Start position to iterate from when creating test DiaObjects

    Returns
    -------
    test_points : `lsst.afw.table.SourceCatalog`
        Catalog of points to test.
    """
    objects = afwTable.SourceCatalog(make_dia_object_schema())

    for src_idx in range(n_points):
        src = objects.addNew()
        src['id'] = src_idx
        src.setCoord(wcs.pixelToSky(startPos + src_idx,
                                    startPos + src_idx))
    return objects


class TestDiaForcedSource(unittest.TestCase):

    def setUp(self):
        # CFHT Filters from the camera mapper.
        self.filter_names = ["u", "g", "r", "i", "z"]
        afwImageUtils.resetFilters()
        afwImageUtils.defineFilter('u', lambdaEff=374, alias="u.MP9301")
        afwImageUtils.defineFilter('g', lambdaEff=487, alias="g.MP9401")
        afwImageUtils.defineFilter('r', lambdaEff=628, alias="r.MP9601")
        afwImageUtils.defineFilter('i', lambdaEff=778, alias="i.MP9701")
        afwImageUtils.defineFilter('z', lambdaEff=1170, alias="z.MP9801")

        # metadata taken from CFHT data
        # v695856-e0/v695856-e0-c000-a00.sci_img.fits
        self.metadata = dafBase.PropertySet()

        self.metadata.set("SIMPLE", "T")
        self.metadata.set("BITPIX", -32)
        self.metadata.set("NAXIS", 2)
        self.metadata.set("NAXIS1", 1024)
        self.metadata.set("NAXIS2", 1153)
        self.metadata.set("RADECSYS", 'FK5')
        self.metadata.set("EQUINOX", 2000.)

        self.metadata.setDouble("CRVAL1", 215.604025685476)
        self.metadata.setDouble("CRVAL2", 53.1595451514076)
        self.metadata.setDouble("CRPIX1", 1109.99981456774)
        self.metadata.setDouble("CRPIX2", 560.018167811613)
        self.metadata.set("CTYPE1", 'RA---SIN')
        self.metadata.set("CTYPE2", 'DEC--SIN')

        self.metadata.setDouble("CD1_1", 5.10808596133527E-05)
        self.metadata.setDouble("CD1_2", 1.85579539217196E-07)
        self.metadata.setDouble("CD2_2", -5.10281493481982E-05)
        self.metadata.setDouble("CD2_1", -8.27440751733828E-07)

        self.wcs = afwGeom.makeSkyWcs(self.metadata)

        self.calibration = 10000
        self.calibrationErr = 100
        self.exposureId = 1234
        self.exposureTime = 200.
        self.imageSize = [1024, 1153]
        self.dateTime = "2014-05-13T17:00:00.000000000"

        # Make images  with one source in them and distinct values and
        # variance for each image.
        # Direct Image
        source_image = afwImage.MaskedImageF(
            lsst.geom.ExtentI(self.imageSize[0] + 1, self.imageSize[1] + 1))
        source_image.image[100, 100, afwImage.LOCAL] = 10
        source_image.getVariance().set(1)
        bbox = lsst.geom.BoxI(
            lsst.geom.PointI(1, 1),
            lsst.geom.ExtentI(self.imageSize[0],
                              self.imageSize[1]))
        masked_image = afwImage.MaskedImageF(source_image, bbox, afwImage.LOCAL)
        self.exposure = afwImage.makeExposure(masked_image, self.wcs)

        detector = DetectorWrapper(
            id=23, bbox=self.exposure.getBBox()).detector
        visit = afwImage.VisitInfo(
            exposureId=self.exposureId,
            exposureTime=self.exposureTime,
            date=dafBase.DateTime(self.dateTime,
                                  dafBase.DateTime.Timescale.TAI))
        self.exposure.setDetector(detector)
        self.exposure.getInfo().setVisitInfo(visit)
        self.exposure.setFilter(afwImage.Filter('g'))
        self.exposure.setPhotoCalib(afwImage.PhotoCalib(self.calibration, self.calibrationErr))

        # Difference Image
        source_image = afwImage.MaskedImageF(
            lsst.geom.ExtentI(self.imageSize[0] + 1, self.imageSize[1] + 1))
        source_image.image[100, 100, afwImage.LOCAL] = 20
        source_image.getVariance().set(2)
        bbox = lsst.geom.BoxI(
            lsst.geom.PointI(1, 1),
            lsst.geom.ExtentI(self.imageSize[0],
                              self.imageSize[1]))
        masked_image = afwImage.MaskedImageF(source_image, bbox, afwImage.LOCAL)
        self.diffim = afwImage.makeExposure(masked_image, self.wcs)
        self.diffim.setDetector(detector)
        self.diffim.getInfo().setVisitInfo(visit)
        self.diffim.setFilter(afwImage.Filter('g'))
        self.diffim.setPhotoCalib(afwImage.PhotoCalib(self.calibration, self.calibrationErr))

        self.expIdBits = 16

        FWHM = 5
        psf = measAlg.DoubleGaussianPsf(15, 15, FWHM/(2*np.sqrt(2*np.log(2))))
        self.exposure.setPsf(psf)
        self.diffim.setPsf(psf)

        self.testDiaObjects = create_test_dia_objects(5, self.wcs)

    def testRun(self):
        """Test that forced source catalogs are successfully created and have
        sensible values.
        """
        dfs = DiaForcedSourceTask()
        result = dfs.run(
            self.testDiaObjects, self.expIdBits, self.exposure, self.diffim)

        direct_values = [19.98544842, 16.00974072, 8.2299179, 2.71486044, 0.57469884]
        direct_var = 7.52143925
        diff_values = [39.97089683, 32.01948144, 16.45983579, 5.42972089, 1.14939768]
        diff_var = 10.6369214

        self.assertEqual(len(result.directForcedSources), len(self.testDiaObjects))
        self.assertEqual(len(result.diffForcedSources), len(self.testDiaObjects))

        for dirRec, diffRec, testObj, dirVal, diffVal in zip(result.directForcedSources,
                                                             result.diffForcedSources,
                                                             self.testDiaObjects,
                                                             direct_values,
                                                             diff_values):
            self.assertEqual(dirRec["id"], diffRec["id"])
            self.assertEqual(dirRec["coord_ra"], testObj["coord_ra"])
            self.assertEqual(dirRec["coord_dec"], testObj["coord_dec"])
            self.assertEqual(diffRec["coord_ra"], testObj["coord_ra"])
            self.assertEqual(diffRec["coord_dec"], testObj["coord_dec"])

            self.assertAlmostEqual(dirRec["slot_PsfFlux_instFlux"], dirVal)
            self.assertAlmostEqual(dirRec["slot_PsfFlux_instFluxErr"],
                                   direct_var)

            self.assertAlmostEqual(diffRec["slot_PsfFlux_instFlux"], diffVal)
            self.assertAlmostEqual(diffRec["slot_PsfFlux_instFluxErr"],
                                   diff_var)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
