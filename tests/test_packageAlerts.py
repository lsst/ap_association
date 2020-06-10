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
import shutil
import tempfile
import unittest

from lsst.ap.association import (PackageAlertsConfig,
                                 PackageAlertsTask,
                                 make_dia_source_schema,
                                 make_dia_object_schema)
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
import lsst.afw.fits as afwFits
import lsst.afw.image as afwImage
import lsst.afw.image.utils as afwImageUtils
import lsst.daf.base as dafBase
from lsst.dax.apdb import Apdb, ApdbConfig
import lsst.geom as geom
import lsst.meas.base.tests
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


def makeDiaObjects(nObjects, exposure):
    """Make a test set of DiaObjects.

    Parameters
    ----------
    nObjects : `int`
        Number of objects to create.
    exposure : `lsst.afw.image.Exposure`
        Exposure to create objects over.

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
        htmIdx = 1
        newObject = {"ra": coord.getRa().asDegrees(),
                     "decl": coord.getDec().asDegrees(),
                     "radecTai": midPointTaiMJD,
                     "diaObjectId": idx,
                     "pixelId": htmIdx,
                     "pmParallaxNdata": 0,
                     "nearbyObj1": 0,
                     "nearbyObj2": 0,
                     "nearbyObj3": 0,
                     "flags": 1,
                     "nDiaSources": 5}
        for f in ["u", "g", "r", "i", "z", "y"]:
            newObject["%sPSFluxNdata" % f] = 0
        data.append(newObject)

    return pd.DataFrame(data=data)


def makeDiaSources(nSources, diaObjectIds, exposure):
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

    midPointTaiMJD = exposure.getInfo().getVisitInfo().getDate().get(
        system=dafBase.DateTime.MJD)

    wcs = exposure.getWcs()
    ccdVisitId = exposure.getInfo().getVisitInfo().getExposureId()

    data = []
    for idx, (x, y) in enumerate(zip(rand_x, rand_y)):
        coord = wcs.pixelToSky(x, y)
        htmIdx = 1
        objId = diaObjectIds[idx % len(diaObjectIds)]
        # Put together the minimum values for the alert.
        data.append({"ra": coord.getRa().asDegrees(),
                     "decl": coord.getDec().asDegrees(),
                     "x": x,
                     "y": y,
                     "ccdVisitId": ccdVisitId,
                     "diaObjectId": objId,
                     "ssObjectId": 0,
                     "parentDiaSourceId": 0,
                     "prv_procOrder": 0,
                     "diaSourceId": idx,
                     "pixelId": htmIdx,
                     "midPointTai": midPointTaiMJD + 1.0 * idx,
                     "filterName": exposure.getFilter().getCanonicalName(),
                     "filterId": 0,
                     "psNdata": 0,
                     "trailNdata": 0,
                     "dipNdata": 0,
                     "flags": 1})

    return pd.DataFrame(data=data)


def _roundTripThroughApdb(objects, sources, dateTime):
    """Run object and source catalogs through the Apdb to get the correct
    table schemas.

    Parameters
    ----------
    objects : `pandas.DataFrame`
        Set of test DiaObjects to round trip.
    sources : `pandas.DataFrame`
        Set of test DiaSources to round trip.
    dateTime : `datetime.datetime`
        Time for the Apdb.

    Returns
    -------
    objects : `pandas.DataFrame`
        Round tripped objects.
    sources : `pandas.DataFrame`
        Round tripped sources.
    """
    tmpFile = tempfile.NamedTemporaryFile()

    apdbConfig = ApdbConfig()
    apdbConfig.db_url = "sqlite:///" + tmpFile.name
    apdbConfig.isolation_level = "READ_UNCOMMITTED"
    apdbConfig.dia_object_index = "baseline"
    apdbConfig.dia_object_columns = []
    apdbConfig.schema_file = _data_file_name(
        "apdb-schema.yaml", "dax_apdb")
    apdbConfig.column_map = _data_file_name(
        "apdb-ap-pipe-afw-map.yaml", "ap_association")
    apdbConfig.extra_schema_file = _data_file_name(
        "apdb-ap-pipe-schema-extra.yaml", "ap_association")

    apdb = Apdb(config=apdbConfig,
                afw_schemas=dict(DiaObject=make_dia_object_schema(),
                                 DiaSource=make_dia_source_schema()))
    apdb.makeSchema()

    minId = objects["pixelId"].min()
    maxId = objects["pixelId"].max()
    diaObjects = apdb.getDiaObjects([[minId, maxId + 1]], return_pandas=True).append(objects)
    diaSources = apdb.getDiaSources(np.unique(objects["diaObjectId"]),
                                    dateTime,
                                    return_pandas=True).append(sources)

    apdb.storeDiaSources(diaSources)
    apdb.storeDiaObjects(diaObjects, dateTime)

    diaObjects = apdb.getDiaObjects([[minId, maxId + 1]], return_pandas=True)
    diaSources = apdb.getDiaSources(np.unique(diaObjects["diaObjectId"]),
                                    dateTime,
                                    return_pandas=True)
    diaObjects.set_index("diaObjectId", drop=False, inplace=True)
    diaSources.set_index(["diaObjectId", "filterName", "diaSourceId"],
                         drop=False,
                         inplace=True)

    return (diaObjects, diaSources)


class TestPackageAlerts(unittest.TestCase):

    def setUp(self):
        np.random.seed(1234)
        self.cutoutSize = 35
        self.center = lsst.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -30),
                                    lsst.geom.Extent2I(140, 160))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(100000.0, self.center)
        exposure, catalog = self.dataset.realize(
            10.0,
            self.dataset.makeMinimalSchema(),
            randomSeed=0)
        self.exposure = exposure
        detector = DetectorWrapper(id=23, bbox=exposure.getBBox()).detector
        self.exposure.setDetector(detector)

        visit = afwImage.VisitInfo(
            exposureId=1234,
            exposureTime=200.,
            date=dafBase.DateTime("2014-05-13T17:00:00.000000000",
                                  dafBase.DateTime.Timescale.TAI))
        self.exposure.getInfo().setVisitInfo(visit)

        self.filter_names = ["g"]
        afwImageUtils.resetFilters()
        afwImageUtils.defineFilter('g', lambdaEff=487, alias="g.MP9401")
        self.exposure.setFilter(afwImage.Filter('g'))

        diaObjects = makeDiaObjects(2, self.exposure)
        diaSourceHistory = makeDiaSources(10,
                                          diaObjects["diaObjectId"],
                                          self.exposure)
        self.diaObjects, diaSourceHistory = _roundTripThroughApdb(
            diaObjects,
            diaSourceHistory,
            self.exposure.getInfo().getVisitInfo().getDate().toPython())
        self.diaObjects.replace(to_replace=[None], value=np.nan, inplace=True)
        diaSourceHistory.replace(to_replace=[None], value=np.nan, inplace=True)
        diaSourceHistory["programId"] = 0

        self.diaSources = diaSourceHistory.loc[
            [(0, "g", 8), (1, "g", 9)], :]
        self.diaSources["bboxSize"] = self.cutoutSize
        self.diaSourceHistory = diaSourceHistory.drop(labels=[(0, "g", 8),
                                                              (1, "g", 9)])

    def testCreateBBox(self):
        """Test the bbox creation
        """
        packConfig = PackageAlertsConfig()
        # Just create a minimum less than the default cutout.
        packConfig.minCutoutSize = self.cutoutSize - 5
        packageAlerts = PackageAlertsTask(config=packConfig)
        bbox = packageAlerts.createDiaSourceBBox(packConfig.minCutoutSize - 5)
        self.assertTrue(bbox == geom.Extent2I(packConfig.minCutoutSize,
                                              packConfig.minCutoutSize))
        # Test that the cutout size is correct.
        bbox = packageAlerts.createDiaSourceBBox(self.cutoutSize)
        self.assertTrue(bbox == geom.Extent2I(self.cutoutSize,
                                              self.cutoutSize))

    def testMakeCutoutBytes(self):
        """Test round tripping an exposure/cutout to bytes and back.
        """
        packageAlerts = PackageAlertsTask()

        sphPoint = self.exposure.getWcs().pixelToSky(self.center)
        cutout = self.exposure.getCutout(sphPoint,
                                         geom.Extent2I(self.cutoutSize,
                                                       self.cutoutSize))

        cutoutBytes = packageAlerts.makeCutoutBytes(cutout)
        tempMemFile = afwFits.MemFileManager(len(cutoutBytes))
        tempMemFile.setData(cutoutBytes, len(cutoutBytes))
        cutoutFromBytes = afwImage.ExposureF(tempMemFile)
        self.assertTrue(
            np.all(cutout.getImage().array == cutoutFromBytes.getImage().array))

    def testMakeAlertDict(self):
        """Test stripping data from the various data products and into a
        dictionary "alert".
        """
        packageAlerts = PackageAlertsTask()
        alertId = 1234

        for srcIdx, diaSource in self.diaSources.iterrows():
            sphPoint = geom.SpherePoint(diaSource["ra"],
                                        diaSource["decl"],
                                        geom.degrees)
            cutout = self.exposure.getCutout(sphPoint,
                                             geom.Extent2I(self.cutoutSize,
                                                           self.cutoutSize))
            cutputBytes = packageAlerts.makeCutoutBytes(cutout)
            objSources = self.diaSourceHistory.loc[srcIdx[0]]
            alert = packageAlerts.makeAlertDict(
                alertId,
                diaSource,
                self.diaObjects.loc[srcIdx[0]],
                objSources,
                cutout,
                None)
            self.assertEqual(len(alert), 9)

            self.assertEqual(alert["alertId"], alertId)
            self.assertEqual(alert["diaSource"], diaSource.to_dict())
            self.assertEqual(alert["cutoutDifference"],
                             cutputBytes)
            self.assertEqual(alert["cutoutTemplate"],
                             None)

    def testRun(self):
        """Test the run method of package alerts.
        """
        packConfig = PackageAlertsConfig()
        tempdir = tempfile.mkdtemp(prefix='alerts')
        packConfig.alertWriteLocation = tempdir
        packageAlerts = PackageAlertsTask(config=packConfig)

        packageAlerts.run(self.diaSources,
                          self.diaObjects,
                          self.diaSourceHistory,
                          self.exposure,
                          None,
                          None)

        ccdVisitId = self.exposure.getInfo().getVisitInfo().getExposureId()
        with open(os.path.join(tempdir, f"{ccdVisitId}.avro"), 'rb') as f:
            writer_schema, data = \
                packageAlerts.alertSchema.retrieve_alerts(f)
        self.assertEqual(len(data), len(self.diaSources))
        for idx, alert in enumerate(data):
            for key, value in alert["diaSource"].items():
                if isinstance(value, float):
                    if np.isnan(self.diaSources.iloc[idx][key]):
                        self.assertTrue(np.isnan(value))
                    else:
                        self.assertAlmostEqual(
                            1 - value / self.diaSources.iloc[idx][key],
                            0.)
                else:
                    self.assertEqual(value, self.diaSources.iloc[idx][key])
            sphPoint = geom.SpherePoint(alert["diaSource"]["ra"],
                                        alert["diaSource"]["decl"],
                                        geom.degrees)
            cutout = self.exposure.getCutout(sphPoint,
                                             geom.Extent2I(self.cutoutSize,
                                                           self.cutoutSize))
            self.assertEqual(alert["cutoutDifference"],
                             packageAlerts.makeCutoutBytes(cutout))

        shutil.rmtree(tempdir)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
