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

import io
import os
import numpy as np
import pandas as pd
import shutil
import tempfile
import unittest

from astropy import wcs
from astropy.nddata import CCDData

from lsst.ap.association import PackageAlertsConfig, PackageAlertsTask
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
from lsst.dax.apdb import ApdbSql, ApdbSqlConfig
import lsst.geom as geom
import lsst.meas.base.tests
from lsst.sphgeom import Box
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
        newObject = {"ra": coord.getRa().asDegrees(),
                     "decl": coord.getDec().asDegrees(),
                     "radecTai": midPointTaiMJD,
                     "diaObjectId": idx + 1,
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
    ccdVisitId = exposure.info.id

    data = []
    for idx, (x, y) in enumerate(zip(rand_x, rand_y)):
        coord = wcs.pixelToSky(x, y)
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
                     "diaSourceId": idx + 1,
                     "midPointTai": midPointTaiMJD + 1.0 * idx,
                     "filterName": exposure.getFilter().bandLabel,
                     "psNdata": 0,
                     "trailNdata": 0,
                     "dipNdata": 0,
                     "flags": 1})

    return pd.DataFrame(data=data)


def makeDiaForcedSources(nSources, diaObjectIds, exposure):
    """Make a test set of DiaSources.

    Parameters
    ----------
    nSources : `int`
        Number of sources to create.
    diaObjectIds : `numpy.ndarray`
        Integer Ids of diaobjects to "associate" with the DiaSources.
    exposure : `lsst.afw.image.Exposure`
        Exposure to create sources over.

    Returns
    -------
    diaSources : `pandas.DataFrame`
        DiaSources generated across the exposure.
    """
    midPointTaiMJD = exposure.getInfo().getVisitInfo().getDate().get(
        system=dafBase.DateTime.MJD)

    ccdVisitId = exposure.info.id

    data = []
    for idx in range(nSources):
        objId = diaObjectIds[idx % len(diaObjectIds)]
        # Put together the minimum values for the alert.
        data.append({"diaForcedSourceId": idx + 1,
                     "ccdVisitId": ccdVisitId + idx,
                     "diaObjectId": objId,
                     "midPointTai": midPointTaiMJD + 1.0 * idx,
                     "filterName": exposure.getFilter().bandLabel,
                     "flags": 0})

    return pd.DataFrame(data=data)


def _roundTripThroughApdb(objects, sources, forcedSources, dateTime):
    """Run object and source catalogs through the Apdb to get the correct
    table schemas.

    Parameters
    ----------
    objects : `pandas.DataFrame`
        Set of test DiaObjects to round trip.
    sources : `pandas.DataFrame`
        Set of test DiaSources to round trip.
    forcedSources : `pandas.DataFrame`
        Set of test DiaForcedSources to round trip.
    dateTime : `lsst.daf.base.DateTime`
        Time for the Apdb.

    Returns
    -------
    objects : `pandas.DataFrame`
        Round tripped objects.
    sources : `pandas.DataFrame`
        Round tripped sources.
    """
    tmpFile = tempfile.NamedTemporaryFile()

    apdbConfig = ApdbSqlConfig()
    apdbConfig.db_url = "sqlite:///" + tmpFile.name
    apdbConfig.dia_object_index = "baseline"
    apdbConfig.dia_object_columns = []

    apdb = ApdbSql(config=apdbConfig)
    apdb.makeSchema()

    wholeSky = Box.full()
    diaObjects = pd.concat([apdb.getDiaObjects(wholeSky), objects])
    diaSources = pd.concat([apdb.getDiaSources(wholeSky, [], dateTime), sources])
    diaForcedSources = pd.concat([apdb.getDiaForcedSources(wholeSky, [], dateTime), forcedSources])

    apdb.store(dateTime, diaObjects, diaSources, diaForcedSources)

    diaObjects = apdb.getDiaObjects(wholeSky)
    diaSources = apdb.getDiaSources(wholeSky, np.unique(diaObjects["diaObjectId"]), dateTime)
    diaForcedSources = apdb.getDiaForcedSources(
        wholeSky, np.unique(diaObjects["diaObjectId"]), dateTime)

    diaObjects.set_index("diaObjectId", drop=False, inplace=True)
    diaSources.set_index(["diaObjectId", "filterName", "diaSourceId"],
                         drop=False,
                         inplace=True)
    diaForcedSources.set_index(["diaObjectId"], drop=False, inplace=True)

    return (diaObjects, diaSources, diaForcedSources)


class TestPackageAlerts(lsst.utils.tests.TestCase):

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
        self.exposure.info.id = 1234
        self.exposure.getInfo().setVisitInfo(visit)

        self.exposure.setFilter(afwImage.FilterLabel(band='g', physical="g.MP9401"))

        diaObjects = makeDiaObjects(2, self.exposure)
        diaSourceHistory = makeDiaSources(10,
                                          diaObjects["diaObjectId"],
                                          self.exposure)
        diaForcedSources = makeDiaForcedSources(10,
                                                diaObjects["diaObjectId"],
                                                self.exposure)
        self.diaObjects, diaSourceHistory, self.diaForcedSources = _roundTripThroughApdb(
            diaObjects,
            diaSourceHistory,
            diaForcedSources,
            self.exposure.getInfo().getVisitInfo().getDate())
        self.diaObjects.replace(to_replace=[None], value=np.nan, inplace=True)
        diaSourceHistory.replace(to_replace=[None], value=np.nan, inplace=True)
        self.diaForcedSources.replace(to_replace=[None], value=np.nan, inplace=True)
        diaSourceHistory["programId"] = 0

        self.diaSources = diaSourceHistory.loc[
            [(1, "g", 9), (2, "g", 10)], :]
        self.diaSources["bboxSize"] = self.cutoutSize
        self.diaSourceHistory = diaSourceHistory.drop(labels=[(1, "g", 9),
                                                              (2, "g", 10)])

        self.cutoutWcs = wcs.WCS(naxis=2)
        self.cutoutWcs.wcs.crpix = [self.center[0], self.center[1]]
        self.cutoutWcs.wcs.crval = [
            self.exposure.getWcs().getSkyOrigin().getRa().asDegrees(),
            self.exposure.getWcs().getSkyOrigin().getDec().asDegrees()]
        self.cutoutWcs.wcs.cd = self.exposure.getWcs().getCdMatrix()
        self.cutoutWcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    def testCreateExtent(self):
        """Test the extent creation for the cutout bbox.
        """
        packConfig = PackageAlertsConfig()
        # Just create a minimum less than the default cutout.
        packConfig.minCutoutSize = self.cutoutSize - 5
        packageAlerts = PackageAlertsTask(config=packConfig)
        extent = packageAlerts.createDiaSourceExtent(
            packConfig.minCutoutSize - 5)
        self.assertTrue(extent == geom.Extent2I(packConfig.minCutoutSize,
                                                packConfig.minCutoutSize))
        # Test that the cutout size is correct.
        extent = packageAlerts.createDiaSourceExtent(self.cutoutSize)
        self.assertTrue(extent == geom.Extent2I(self.cutoutSize,
                                                self.cutoutSize))

    def testCreateCcdDataCutout(self):
        """Test that the data is being extracted into the CCDData cutout
        correctly.
        """
        packageAlerts = PackageAlertsTask()

        diaSrcId = 1234
        ccdData = packageAlerts.createCcdDataCutout(
            self.exposure,
            self.exposure.getWcs().getSkyOrigin(),
            self.exposure.getBBox().getDimensions(),
            self.exposure.getPhotoCalib(),
            diaSrcId)
        calibExposure = self.exposure.getPhotoCalib().calibrateImage(
            self.exposure.getMaskedImage())

        self.assertFloatsAlmostEqual(ccdData.wcs.wcs.cd,
                                     self.cutoutWcs.wcs.cd)
        self.assertFloatsAlmostEqual(ccdData.data,
                                     calibExposure.getImage().array)

        ccdData = packageAlerts.createCcdDataCutout(
            self.exposure,
            geom.SpherePoint(0, 0, geom.degrees),
            self.exposure.getBBox().getDimensions(),
            self.exposure.getPhotoCalib(),
            diaSrcId)
        self.assertTrue(ccdData is None)

    def testMakeLocalTransformMatrix(self):
        """Test that the local WCS approximation is correct.
        """
        packageAlerts = PackageAlertsTask()

        sphPoint = self.exposure.getWcs().pixelToSky(self.center)
        cutout = self.exposure.getCutout(sphPoint,
                                         geom.Extent2I(self.cutoutSize,
                                                       self.cutoutSize))
        cd = packageAlerts.makeLocalTransformMatrix(
            cutout.getWcs(), self.center, sphPoint)
        self.assertFloatsAlmostEqual(
            cd,
            cutout.getWcs().getCdMatrix(),
            rtol=1e-11,
            atol=1e-11)

    def testStreamCcdDataToBytes(self):
        """Test round tripping an CCDData cutout to bytes and back.
        """
        packageAlerts = PackageAlertsTask()

        sphPoint = self.exposure.getWcs().pixelToSky(self.center)
        cutout = self.exposure.getCutout(sphPoint,
                                         geom.Extent2I(self.cutoutSize,
                                                       self.cutoutSize))
        cutoutCcdData = CCDData(
            data=cutout.getImage().array,
            wcs=self.cutoutWcs,
            unit="adu")

        cutoutBytes = packageAlerts.streamCcdDataToBytes(cutoutCcdData)
        with io.BytesIO(cutoutBytes) as bytesIO:
            cutoutFromBytes = CCDData.read(bytesIO, format="fits")
        self.assertFloatsAlmostEqual(cutoutCcdData.data, cutoutFromBytes.data)

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
            ccdCutout = packageAlerts.createCcdDataCutout(
                cutout,
                sphPoint,
                geom.Extent2I(self.cutoutSize, self.cutoutSize),
                cutout.getPhotoCalib(),
                1234)
            cutoutBytes = packageAlerts.streamCcdDataToBytes(
                ccdCutout)
            objSources = self.diaSourceHistory.loc[srcIdx[0]]
            objForcedSources = self.diaForcedSources.loc[srcIdx[0]]
            alert = packageAlerts.makeAlertDict(
                alertId,
                diaSource,
                self.diaObjects.loc[srcIdx[0]],
                objSources,
                objForcedSources,
                ccdCutout,
                ccdCutout)
            self.assertEqual(len(alert), 9)

            self.assertEqual(alert["alertId"], alertId)
            self.assertEqual(alert["diaSource"], diaSource.to_dict())
            self.assertEqual(alert["cutoutDifference"],
                             cutoutBytes)
            self.assertEqual(alert["cutoutTemplate"],
                             cutoutBytes)

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
                          self.diaForcedSources,
                          self.exposure,
                          self.exposure,
                          None)

        ccdVisitId = self.exposure.info.id
        with open(os.path.join(tempdir, f"{ccdVisitId}.avro"), 'rb') as f:
            writer_schema, data_stream = \
                packageAlerts.alertSchema.retrieve_alerts(f)
            data = list(data_stream)
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
            ccdCutout = packageAlerts.createCcdDataCutout(
                cutout,
                sphPoint,
                geom.Extent2I(self.cutoutSize, self.cutoutSize),
                cutout.getPhotoCalib(),
                1234)
            self.assertEqual(alert["cutoutDifference"],
                             packageAlerts.streamCcdDataToBytes(ccdCutout))

        shutil.rmtree(tempdir)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
