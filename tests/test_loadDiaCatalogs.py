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
import tempfile
import unittest
import yaml

from lsst.afw.cameraGeom.testUtils import DetectorWrapper
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
from lsst.ap.association import LoadDiaCatalogsTask, LoadDiaCatalogsConfig
import lsst.daf.base as dafBase
from lsst.dax.apdb import ApdbSql, ApdbSqlConfig, ApdbTables
from lsst.utils import getPackageDir
import lsst.utils.tests
import utils_tests


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
    exposure.info.id = 1234
    exposure.setDetector(detector)
    exposure.getInfo().setVisitInfo(visit)
    exposure.setFilter(afwImage.FilterLabel(band='g'))

    return exposure


class TestLoadDiaCatalogs(unittest.TestCase):

    def setUp(self):
        np.random.seed(1234)

        self.db_file_fd, self.db_file = tempfile.mkstemp(
            dir=os.path.dirname(__file__))

        self.apdbConfig = ApdbSqlConfig()
        self.apdbConfig.db_url = "sqlite:///" + self.db_file
        self.apdbConfig.dia_object_index = "baseline"
        self.apdbConfig.dia_object_columns = []

        self.apdb = ApdbSql(config=self.apdbConfig)
        self.apdb.makeSchema()

        self.exposure = makeExposure(False, False)

        self.diaObjects = utils_tests.makeDiaObjects(20, self.exposure)
        self.diaSources = utils_tests.makeDiaSources(
            100,
            self.diaObjects["diaObjectId"].to_numpy(),
            self.exposure)
        self.diaForcedSources = utils_tests.makeDiaForcedSources(
            200,
            self.diaObjects["diaObjectId"].to_numpy(),
            self.exposure)

        self.dateTime = self.exposure.getInfo().getVisitInfo().getDate()
        self.apdb.store(self.dateTime,
                        self.diaObjects,
                        self.diaSources,
                        self.diaForcedSources)

        # These columns are not in the DPDD, yet do appear in DiaSource.yaml.
        # We don't need to check them against the default APDB schema.
        self.ignoreColumns = ["band", "bboxSize", "isDipole"]

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
        self.assertEqual(len(result.diaForcedSources),
                         len(self.diaForcedSources))

    def testLoadDiaObjects(self):
        """Test that the correct number of diaObjects are loaded.
        """
        diaLoader = LoadDiaCatalogsTask()
        region = diaLoader._getRegion(self.exposure)
        diaObjects = diaLoader.loadDiaObjects(region,
                                              self.apdb)
        self.assertEqual(len(diaObjects), len(self.diaObjects))

    def testLoadDiaForcedSources(self):
        """Test that the correct number of diaForcedSources are loaded.
        """
        diaLoader = LoadDiaCatalogsTask()
        region = diaLoader._getRegion(self.exposure)
        diaForcedSources = diaLoader.loadDiaForcedSources(
            self.diaObjects,
            region,
            self.dateTime,
            self.apdb)
        self.assertEqual(len(diaForcedSources), len(self.diaForcedSources))

    def testLoadDiaSources(self):
        """Test that the correct number of diaSources are loaded.

        Also check that they can be properly loaded both by location and
        ``diaObjectId``.
        """
        diaConfig = LoadDiaCatalogsConfig()
        diaLoader = LoadDiaCatalogsTask(config=diaConfig)

        region = diaLoader._getRegion(self.exposure)
        diaSources = diaLoader.loadDiaSources(self.diaObjects,
                                              region,
                                              self.dateTime,
                                              self.apdb)
        self.assertEqual(len(diaSources), len(self.diaSources))

    def test_apdbSchema(self):
        """Test that the default DiaSource schema from dax_apdb agrees with the
        column names defined here in ap_association/data/DiaSource.yaml.
        """
        tableDef = self.apdb.tableDef(ApdbTables.DiaSource)
        apdbSchemaColumns = [column.name for column in tableDef.columns]

        functorFile = _data_file_name("DiaSource.yaml", "ap_association")
        with open(functorFile) as yaml_stream:
            diaSourceFunctor = yaml.safe_load_all(yaml_stream)
            for functor in diaSourceFunctor:
                diaSourceColumns = [column for column in list(functor['funcs'].keys())
                                    if column not in self.ignoreColumns]
            self.assertLess(set(diaSourceColumns), set(apdbSchemaColumns))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
