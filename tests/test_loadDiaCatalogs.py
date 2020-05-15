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

import os
import numpy as np
import pandas as pd
import tempfile
import unittest

from lsst.afw.cameraGeom.testUtils import DetectorWrapper
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.image.utils as afwImageUtils
from lsst.ap.association import (LoadDiaCatalogsTask,
                                 LoadDiaCatalogsConfig,
                                 make_dia_source_schema,
                                 make_dia_object_schema)
import lsst.daf.base as dafBase
from lsst.dax.apdb import Apdb, ApdbConfig
import lsst.geom as geom
import lsst.sphgeom as sphgeom
from lsst.utils import getPackageDir
import lsst.utils.tests


def _data_file_name(basename, module_name):
    """Return path name of a data file.

    Parameters
    ----------
    basename : `str`
        Name of the file to add to the path string.
    module_name : `str`
        Name of lsst stack package environment variable.

    Returns
    -------
    data_file_path : `str`
       Full path of the file to load from the "data" directory in a given
       repository.
    """
    return os.path.join(getPackageDir(module_name), "data", basename)


def makeExposure(flipX=False, flipY=False):
    """Create an exposure and flip the x or y (or both) coordinates.

    Returns bounding boxes that are right or left handed around the bounding
    polygon.

    Parameters
    ----------
    flipX : `bool`
        Flip the x coordinate in the WCS.
    flipY : `bool`
        Flip the y coordinate in the WCS.

    Returns
    -------
    exposure : `lsst.afw.image.Exposure`
        Exposure with a valid bounding box and wcs.
    """
    metadata = dafBase.PropertySet()

    metadata.set("SIMPLE", "T")
    metadata.set("BITPIX", -32)
    metadata.set("NAXIS", 2)
    metadata.set("NAXIS1", 1024)
    metadata.set("NAXIS2", 1153)
    metadata.set("RADECSYS", 'FK5')
    metadata.set("EQUINOX", 2000.)

    metadata.setDouble("CRVAL1", 215.604025685476)
    metadata.setDouble("CRVAL2", 53.1595451514076)
    metadata.setDouble("CRPIX1", 1109.99981456774)
    metadata.setDouble("CRPIX2", 560.018167811613)
    metadata.set("CTYPE1", 'RA---SIN')
    metadata.set("CTYPE2", 'DEC--SIN')

    xFlip = 1
    if flipX:
        xFlip = -1
    yFlip = 1
    if flipY:
        yFlip = -1
    metadata.setDouble("CD1_1", xFlip * 5.10808596133527E-05)
    metadata.setDouble("CD1_2", yFlip * 1.85579539217196E-07)
    metadata.setDouble("CD2_2", yFlip * -5.10281493481982E-05)
    metadata.setDouble("CD2_1", xFlip * -8.27440751733828E-07)

    wcs = afwGeom.makeSkyWcs(metadata)
    exposure = afwImage.makeExposure(
        afwImage.makeMaskedImageFromArrays(np.ones((1024, 1153))), wcs)
    detector = DetectorWrapper(id=23, bbox=exposure.getBBox()).detector
    visit = afwImage.VisitInfo(
        exposureId=1234,
        exposureTime=200.,
        date=dafBase.DateTime("2014-05-13T17:00:00.000000000",
                              dafBase.DateTime.Timescale.TAI))
    exposure.setDetector(detector)
    exposure.getInfo().setVisitInfo(visit)
    exposure.setFilter(afwImage.Filter('g'))

    return exposure


def makeDiaObjects(nObjects, exposure, pixelator):
    """Make a test set of DiaObjects.

    Parameters
    ----------
    nObjects : `int`
        Number of objects to create.
    exposure : `lsst.afw.image.Exposure`
        Exposure to create objects over.
    pixelator : `lsst.sphgeom.HtmPixelization`
        Object to compute spatial indicies from.

    Returns
    -------
    diaObjects : `pandas.DataFrame`
        DiaObjects generated across the exposure.
    """
    bbox = geom.Box2D(exposure.getBBox())
    rand_x = np.random.uniform(bbox.getMinX(), bbox.getMaxX(), size=nObjects)
    rand_y = np.random.uniform(bbox.getMinY(), bbox.getMaxY(), size=nObjects)

    midPointTaiMJD = exposure.getInfo().getVisitInfo().getDate().get(
        system=dafBase.DateTime.MJD)

    wcs = exposure.getWcs()

    data = []
    for idx, (x, y) in enumerate(zip(rand_x, rand_y)):
        coord = wcs.pixelToSky(x, y)
        htmIdx = pixelator.index(coord.getVector())
        newObject = {"ra": coord.getRa().asDegrees(),
                     "decl": coord.getRa().asDegrees(),
                     "radecTai": midPointTaiMJD,
                     "diaObjectId": idx,
                     "pixelId": htmIdx,
                     "pmParallaxNdata": 0,
                     "nearbyObj1": 0,
                     "nearbyObj2": 0,
                     "nearbyObj3": 0}
        for f in ["u", "g", "r", "i", "z", "y"]:
            newObject["%sPSFluxNdata" % f] = 0
        data.append(newObject)

    return pd.DataFrame(data=data)


def makeDiaSources(nSources, diaObjectIds, exposure, pixelator):
    """Make a test set of DiaSources.

    Parameters
    ----------
    nSources : `int`
        Number of sources to create.
    diaObjectIds : `numpy.ndarray`
        Integer Ids of diaobjects to "associate" with the DiaSources.
    exposure : `lsst.afw.image.Exposure`
        Exposure to create sources over.
    pixelator : `lsst.sphgeom.HtmPixelization`
        Object to compute spatial indicies from.

    Returns
    -------
    diaSources : `pandas.DataFrame`
        DiaSources generated across the exposure.
    """
    bbox = geom.Box2D(exposure.getBBox())
    rand_x = np.random.uniform(bbox.getMinX(), bbox.getMaxX(), size=nSources)
    rand_y = np.random.uniform(bbox.getMinY(), bbox.getMaxY(), size=nSources)
    rand_ids = diaObjectIds[np.random.randint(len(diaObjectIds), size=nSources)]

    midPointTaiMJD = exposure.getInfo().getVisitInfo().getDate().get(
        system=dafBase.DateTime.MJD)

    wcs = exposure.getWcs()

    data = []
    for idx, (x, y, objId) in enumerate(zip(rand_x, rand_y, rand_ids)):
        coord = wcs.pixelToSky(x, y)
        htmIdx = pixelator.index(coord.getVector())
        data.append({"ra": coord.getRa().asDegrees(),
                     "decl": coord.getRa().asDegrees(),
                     "diaObjectId": objId,
                     "diaSourceId": idx,
                     "pixelId": htmIdx,
                     "midPointTai": midPointTaiMJD})

    return pd.DataFrame(data=data)


class TestLoadDiaCatalogs(unittest.TestCase):

    def setUp(self):
        np.random.seed(1234)

        # CFHT Filters from the camera mapper.
        self.filter_names = ["g"]
        afwImageUtils.resetFilters()
        afwImageUtils.defineFilter('g', lambdaEff=487, alias="g.MP9401")

        self.db_file_fd, self.db_file = tempfile.mkstemp(
            dir=os.path.dirname(__file__))

        self.apdbConfig = ApdbConfig()
        self.apdbConfig.db_url = "sqlite:///" + self.db_file
        self.apdbConfig.isolation_level = "READ_UNCOMMITTED"
        self.apdbConfig.dia_object_index = "baseline"
        self.apdbConfig.dia_object_columns = []
        self.apdbConfig.schema_file = _data_file_name(
            "apdb-schema.yaml", "dax_apdb")
        self.apdbConfig.column_map = _data_file_name(
            "apdb-ap-pipe-afw-map.yaml", "ap_association")
        self.apdbConfig.extra_schema_file = _data_file_name(
            "apdb-ap-pipe-schema-extra.yaml", "ap_association")

        self.apdb = Apdb(config=self.apdbConfig,
                         afw_schemas=dict(DiaObject=make_dia_object_schema(),
                                          DiaSource=make_dia_source_schema()))
        self.apdb.makeSchema()

        # Expected HTM pixel ranges for max range=4 and level = 20. This
        # set of pixels should be same for the WCS created by default in
        # makeExposure and for one with a flipped y axis.
        self.ranges = np.sort(np.array([15154776375296, 15154779521024,
                                        15154788958208, 15154792103936]))

        self.pixelator = sphgeom.HtmPixelization(20)
        self.exposure = makeExposure(False, False)

        self.diaObjects = makeDiaObjects(20, self.exposure, self.pixelator)
        self.diaSources = makeDiaSources(
            100,
            self.diaObjects["diaObjectId"].to_numpy(),
            self.exposure,
            self.pixelator)

        self.apdb.storeDiaSources(self.diaSources)
        self.dateTime = \
            self.exposure.getInfo().getVisitInfo().getDate().toPython()
        self.apdb.storeDiaObjects(self.diaObjects,
                                  self.dateTime)

    def tearDown(self):
        os.close(self.db_file_fd)
        os.remove(self.db_file)

    def testRun(self):
        """Test the full run method for the loader.
        """
        diaLoader = LoadDiaCatalogsTask()
        result = diaLoader.run(self.exposure, self.apdb)

        self.assertEqual(len(result.diaObjects), len(self.diaObjects))
        self.assertEqual(len(result.diaSources), len(self.diaSources))

    def testLoadDiaObjects(self):
        """Test that the correct number of diaObjects are loaded.
        """
        diaLoader = LoadDiaCatalogsTask()
        normPixels = diaLoader._getPixelRanges(self.exposure)
        diaObjects = diaLoader.loadDiaObjects(normPixels,
                                              self.apdb)
        self.assertEqual(len(diaObjects), len(self.diaObjects))

    def testLoadDiaSourcesByPixelId(self):
        """Test that the correct number of diaSources are loaded.

        Also check that they can be properly loaded both by location and
        ``diaObjectId``.
        """
        self._testLoadDiaSources(True)

    def testLoadDiaSourcesByDiaObjectId(self):
        """Test that the correct number of diaSources are loaded.

        Also check that they can be properly loaded both by location and
        ``diaObjectId``.
        """
        self._testLoadDiaSources(False)

    def _testLoadDiaSources(self, loadByPixelId):
        """Test that DiaSources are loaded correctly.

        Parameters
        ----------
        loadByPixelId : `bool`
            Load DiaSources by ``pixelId`` if ``True`` and by ``diaObjectId``
            if ``False``.
        """
        diaConfig = LoadDiaCatalogsConfig()
        diaConfig.loadDiaSourcesByPixelId = loadByPixelId
        diaLoader = LoadDiaCatalogsTask(config=diaConfig)

        normPixels = diaLoader._getPixelRanges(self.exposure)
        diaSources = diaLoader.loadDiaSources(self.diaObjects,
                                              normPixels,
                                              self.dateTime,
                                              self.apdb)
        self.assertEqual(len(diaSources), len(self.diaSources))

    def testGetPixelRanges(self):
        """Test the same pixels are returned for flips/ordering changes in
        the WCS.
        """
        diaConfig = LoadDiaCatalogsConfig()
        diaConfig.htmMaxRanges = 4
        diaLoader = LoadDiaCatalogsTask(config=diaConfig)

        # Make two exposures, one with a flipped y axis to get left vs. right
        # handed.
        exposure = makeExposure(False, False)
        exposureFlip = makeExposure(False, True)

        normPixels = diaLoader._getPixelRanges(exposure)
        flipPixels = diaLoader._getPixelRanges(exposureFlip)

        for normPix, flipPix, testPix in zip(
                np.sort(np.array(normPixels).flatten()),
                np.sort(np.array(flipPixels).flatten()),
                self.ranges):
            self.assertEqual(normPix, flipPix)
            self.assertEqual(normPix, testPix)
            self.assertEqual(flipPix, testPix)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
