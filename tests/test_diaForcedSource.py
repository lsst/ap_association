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
import unittest.mock

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

    # Add points outside of the CCD.
    # Fully outside.
    src = objects.addNew()
    src['id'] = 10000000
    src.setCoord(wcs.pixelToSky(-100000,
                                -100000))
    # y outside
    src = objects.addNew()
    src['id'] = 10000001
    src.setCoord(wcs.pixelToSky(100,
                                -100000))
    # x outside
    src = objects.addNew()
    src['id'] = 10000001
    src.setCoord(wcs.pixelToSky(-100000,
                                100))

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
        self.expected_n_columns = 11

    def testRun(self):
        """Test that forced source catalogs are successfully created and have
        sensible values.
        """
        test_objects = self._convert_to_pandas(self.testDiaObjects)
        test_objects.rename(columns={"id": "diaObjectId"},
                            inplace=True)
        test_objects.set_index("diaObjectId", drop=False, inplace=True)
        testIds = test_objects.loc[:, "diaObjectId"].to_numpy()[:6]

        dfs = DiaForcedSourceTask()
        dia_forced_sources = dfs.run(
            test_objects, testIds, self.expIdBits, self.exposure, self.diffim)

        direct_values = [199854.48417094944, 160097.40719241602,
                         82299.17897267535, 27148.604434624354,
                         5746.988388215507]
        direct_var = [75240.939811168, 75231.42933749466,
                      75218.89495113207, 75214.88248249644,
                      75214.41447602339]
        diff_values = [399708.9683418989, 320194.81438483205,
                       164598.3579453507, 54297.20886924871,
                       11493.976776431015]
        diff_var = [106444.28782374493, 106417.39592887461,
                    106381.94840437356, 106370.59980584883,
                    106369.27608815048]

        # Should be number of test objects minus one as one object is purposely
        # outside of the ccd area.
        self.assertEqual(len(dia_forced_sources), len(self.testDiaObjects) - 2)
        self.assertEqual(len(dia_forced_sources.columns),
                         self.expected_n_columns)

        for (diaFS_id, diaFS), testObj, dirVal, diffVal, dirVar, diffVar in zip(dia_forced_sources.iterrows(),
                                                                                self.testDiaObjects,
                                                                                direct_values,
                                                                                diff_values,
                                                                                direct_var,
                                                                                diff_var):
            self.assertAlmostEqual(diaFS["psFlux"] / diffVal, 1.)
            self.assertAlmostEqual(diaFS["psFluxErr"] / diffVar, 1.)

            self.assertAlmostEqual(diaFS["totFlux"] / dirVal, 1.)
            self.assertAlmostEqual(diaFS["totFluxErr"] / dirVar, 1.)

            self.assertEqual(diaFS["ccdVisitId"], self.exposureId)

    def _convert_to_pandas(self, dia_objects):
        """Convert input afw table to pandas.

        Parameters
        ----------
        dia_objects : `lsst.afw.table.SourceCatalog`
            Catalog to convert

        Returns
        -------
        output_catalog : `pandas.DataFrame`
            Converted catalog
        """
        output_catalog = dia_objects.asAstropy().to_pandas()
        output_catalog.rename(columns={"id": "diaObjectId",
                                       "coord_ra": "ra",
                                       "coord_dec": "decl"},
                              inplace=True)

        output_catalog.loc[:, "ra"] = np.degrees(output_catalog["ra"])
        output_catalog.loc[:, "decl"] = np.degrees(output_catalog["decl"])

        return output_catalog


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
