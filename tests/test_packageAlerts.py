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
from unittest.mock import patch, Mock
import sys
from astropy import wcs
from astropy.nddata import CCDData
import logging

from lsst.ap.association import PackageAlertsConfig, PackageAlertsTask
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
from lsst.dax.apdb import Apdb, ApdbSql, ApdbSqlConfig
import lsst.geom as geom
import lsst.meas.base.tests
from lsst.sphgeom import Box
import lsst.utils.tests
import utils_tests

_log = logging.getLogger("lsst." + __name__)
_log.setLevel(logging.DEBUG)

try:
    import confluent_kafka  # noqa: F401
    from confluent_kafka import KafkaException
except ModuleNotFoundError as e:
    _log.error('Kafka module not found: {}'.format(e))


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

    Apdb.makeSchema(apdbConfig)
    apdb = ApdbSql(config=apdbConfig)

    wholeSky = Box.full()
    diaObjects = pd.concat([apdb.getDiaObjects(wholeSky), objects])
    diaSources = pd.concat(
        [apdb.getDiaSources(wholeSky, [], dateTime), sources])
    diaForcedSources = pd.concat(
        [apdb.getDiaForcedSources(wholeSky, [], dateTime), forcedSources])

    apdb.store(dateTime, diaObjects, diaSources, diaForcedSources)

    diaObjects = apdb.getDiaObjects(wholeSky)
    diaSources = apdb.getDiaSources(wholeSky,
                                    np.unique(diaObjects["diaObjectId"]),
                                    dateTime)
    diaForcedSources = apdb.getDiaForcedSources(
        wholeSky, np.unique(diaObjects["diaObjectId"]), dateTime)

    diaObjects.set_index("diaObjectId", drop=False, inplace=True)
    diaSources.set_index(["diaObjectId", "band", "diaSourceId"],
                         drop=False,
                         inplace=True)
    diaForcedSources.set_index(["diaObjectId"], drop=False, inplace=True)

    return (diaObjects, diaSources, diaForcedSources)


def mock_alert(alert_id):
    """Generate a minimal mock alert.
    """
    return {
        "alertId": alert_id,
        "diaSource": {
            # Below are all the required fields containing random values.
            "midpointMjdTai": 5,
            "diaSourceId": 4,
            "ccdVisitId": 2,
            "band": 'g',
            "ra": 12.5,
            "dec": -16.9,
            # These types are 32-bit floats in the avro schema, so we have to
            # make them that type here, so that they round trip appropriately.
            "x": np.float32(15.7),
            "y": np.float32(89.8),
            "apFlux": np.float32(54.85),
            "apFluxErr": np.float32(70.0),
            "snr": np.float32(6.7),
            "psfFlux": np.float32(700.0),
            "psfFluxErr": np.float32(90.0),
            "flags": 12345,
        }
    }


class TestPackageAlerts(lsst.utils.tests.TestCase):
    kafka_enabled = "confluent_kafka" in sys.modules

    def __init__(self, *args, **kwargs):
        TestPackageAlerts.kafka_enabled = "confluent_kafka" in sys.modules
        _log.debug('TestPackageAlerts: kafka_enabled={}'.format(self.kafka_enabled))
        super(TestPackageAlerts, self).__init__(*args, **kwargs)

    def setUp(self):
        patcher = patch.dict(os.environ, {"AP_KAFKA_PRODUCER_PASSWORD": "fake_password",
                                          "AP_KAFKA_PRODUCER_USERNAME": "fake_username",
                                          "AP_KAFKA_SERVER": "fake_server",
                                          "AP_KAFKA_TOPIC": "fake_topic"})
        self.environ = patcher.start()
        self.addCleanup(patcher.stop)
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
            exposureTime=200.,
            date=dafBase.DateTime("2014-05-13T17:00:00.000000000",
                                  dafBase.DateTime.Timescale.TAI))
        self.exposure.info.id = 1234
        self.exposure.info.setVisitInfo(visit)

        self.exposure.setFilter(
            afwImage.FilterLabel(band='g', physical="g.MP9401"))

        diaObjects = utils_tests.makeDiaObjects(2, self.exposure)
        diaSourceHistory = utils_tests.makeDiaSources(10,
                                                      diaObjects[
                                                          "diaObjectId"],
                                                      self.exposure)
        diaForcedSources = utils_tests.makeDiaForcedSources(10,
                                                            diaObjects[
                                                                "diaObjectId"],
                                                            self.exposure)
        self.diaObjects, diaSourceHistory, self.diaForcedSources = _roundTripThroughApdb(
            diaObjects,
            diaSourceHistory,
            diaForcedSources,
            self.exposure.visitInfo.date)
        self.diaObjects.replace(to_replace=[None], value=np.nan, inplace=True)
        diaSourceHistory.replace(to_replace=[None], value=np.nan, inplace=True)
        self.diaForcedSources.replace(to_replace=[None], value=np.nan,
                                      inplace=True)
        diaSourceHistory["programId"] = 0

        self.diaSources = diaSourceHistory.loc[[(1, "g", 9), (2, "g", 10)], :]
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
                                        diaSource["dec"],
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

    def test_produceAlerts_empty_password(self):
        """ Test that produceAlerts raises if the password is empty or missing.
        """
        self.environ['AP_KAFKA_PRODUCER_PASSWORD'] = ""
        with self.assertRaisesRegex(ValueError, "Kafka password"):
            packConfig = PackageAlertsConfig(doProduceAlerts=True)
            packageAlerts = PackageAlertsTask(config=packConfig)  # noqa: F841

        del self.environ['AP_KAFKA_PRODUCER_PASSWORD']
        with self.assertRaisesRegex(ValueError, "Kafka password"):
            packConfig = PackageAlertsConfig(doProduceAlerts=True)
            packageAlerts = PackageAlertsTask(config=packConfig)  # noqa: F841

    def test_produceAlerts_empty_username(self):
        """ Test that produceAlerts raises if the username is empty or missing.
        """
        self.environ['AP_KAFKA_PRODUCER_USERNAME'] = ""
        with self.assertRaisesRegex(ValueError, "Kafka username"):
            packConfig = PackageAlertsConfig(doProduceAlerts=True)
            packageAlerts = PackageAlertsTask(config=packConfig)  # noqa: F841

        del self.environ['AP_KAFKA_PRODUCER_USERNAME']
        with self.assertRaisesRegex(ValueError, "Kafka username"):
            packConfig = PackageAlertsConfig(doProduceAlerts=True)
            packageAlerts = PackageAlertsTask(config=packConfig)  # noqa: F841

    def test_produceAlerts_empty_server(self):
        """ Test that produceAlerts raises if the server is empty or missing.
        """
        self.environ['AP_KAFKA_SERVER'] = ""
        with self.assertRaisesRegex(ValueError, "Kafka server"):
            packConfig = PackageAlertsConfig(doProduceAlerts=True)
            packageAlerts = PackageAlertsTask(config=packConfig)  # noqa: F841

        del self.environ['AP_KAFKA_SERVER']
        with self.assertRaisesRegex(ValueError, "Kafka server"):
            packConfig = PackageAlertsConfig(doProduceAlerts=True)
            packageAlerts = PackageAlertsTask(config=packConfig)  # noqa: F841

    def test_produceAlerts_empty_topic(self):
        """ Test that produceAlerts raises if the topic is empty or missing.
        """
        self.environ['AP_KAFKA_TOPIC'] = ""
        with self.assertRaisesRegex(ValueError, "Kafka topic"):
            packConfig = PackageAlertsConfig(doProduceAlerts=True)
            packageAlerts = PackageAlertsTask(config=packConfig)  # noqa: F841

        del self.environ['AP_KAFKA_TOPIC']
        with self.assertRaisesRegex(ValueError, "Kafka topic"):
            packConfig = PackageAlertsConfig(doProduceAlerts=True)
            packageAlerts = PackageAlertsTask(config=packConfig)  # noqa: F841

    @patch('confluent_kafka.Producer')
    @unittest.skipIf("confluent_kafka" not in sys.modules, 'Kafka is not enabled')
    def test_produceAlerts_success(self, mock_producer):
        """ Test that produceAlerts calls the producer on all provided alerts
        when the alerts are all under the batch size limit.
        """
        packConfig = PackageAlertsConfig(doProduceAlerts=True)
        packageAlerts = PackageAlertsTask(config=packConfig)
        alerts = [mock_alert(1), mock_alert(2)]
        ccdVisitId = 123

        # Create a variable and assign it an instance of the patched kafka producer
        producer_instance = mock_producer.return_value
        producer_instance.produce = Mock()
        producer_instance.flush = Mock()
        packageAlerts.produceAlerts(alerts, ccdVisitId)

        self.assertEqual(producer_instance.produce.call_count, len(alerts))
        self.assertEqual(producer_instance.flush.call_count, len(alerts)+1)

    @patch('confluent_kafka.Producer')
    @unittest.skipIf("confluent_kafka" not in sys.modules, 'Kafka is not enabled')
    def test_produceAlerts_one_failure(self, mock_producer):
        """ Test that produceAlerts correctly fails on one alert
        and is writing the failure to disk.
        """
        counter = 0

        # confluent_kafka is not visible to mock_producer and needs to be
        # re-imported here.
        def mock_produce(*args, **kwargs):
            nonlocal counter
            counter += 1
            if counter == 2:
                raise KafkaException
            else:
                return

        packConfig = PackageAlertsConfig(doProduceAlerts=True)
        packageAlerts = PackageAlertsTask(config=packConfig)

        patcher = patch("builtins.open")
        patch_open = patcher.start()
        alerts = [mock_alert(1), mock_alert(2), mock_alert(3)]
        ccdVisitId = 123

        producer_instance = mock_producer.return_value
        producer_instance.produce = Mock(side_effect=mock_produce)
        producer_instance.flush = Mock()

        packageAlerts.produceAlerts(alerts, ccdVisitId)

        self.assertEqual(producer_instance.produce.call_count, len(alerts))
        self.assertEqual(patch_open.call_count, 1)
        self.assertIn("123_2.avro", patch_open.call_args.args[0])
        # Because one produce raises, we call flush one fewer times than in the success
        # test above.
        self.assertEqual(producer_instance.flush.call_count, len(alerts))
        patcher.stop()

    def testRun_without_produce(self):
        """Test the run method of package alerts with produce set to False.
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
                          self.exposure)

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
                                        alert["diaSource"]["dec"],
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

    @patch.object(PackageAlertsTask, 'produceAlerts')
    @patch('confluent_kafka.Producer')
    @unittest.skipIf("confluent_kafka" not in sys.modules, 'Kafka is not enabled')
    def testRun_with_produce(self, mock_produceAlerts, mock_producer):
        """Test that packageAlerts calls produceAlerts when doProduceAlerts
        is set to True.
        """
        packConfig = PackageAlertsConfig(doProduceAlerts=True)
        packageAlerts = PackageAlertsTask(config=packConfig)

        packageAlerts.run(self.diaSources,
                          self.diaObjects,
                          self.diaSourceHistory,
                          self.diaForcedSources,
                          self.exposure,
                          self.exposure)

        self.assertEqual(mock_produceAlerts.call_count, 1)

    def testRun_without_produce_or_write(self):
        """Test that packageAlerts calls produceAlerts when doProduceAlerts
        is set to True.
        """
        packConfig = PackageAlertsConfig(doProduceAlerts=False,
                                         doWriteAlerts=False)
        packageAlerts = PackageAlertsTask(config=packConfig)

        with self.assertRaisesRegex(Exception, "Neither produce alerts"):
            packageAlerts.run(self.diaSources,
                              self.diaObjects,
                              self.diaSourceHistory,
                              self.diaForcedSources,
                              self.exposure,
                              self.exposure)

    def test_serialize_alert_round_trip(self):
        """Test that values in the alert packet exactly round trip.
        """
        ConfigClass = PackageAlertsConfig()
        packageAlerts = PackageAlertsTask(config=ConfigClass)

        alert = mock_alert(1)
        serialized = PackageAlertsTask._serialize_alert(packageAlerts, alert)
        deserialized = PackageAlertsTask._deserialize_alert(packageAlerts, serialized)

        for field in alert['diaSource']:
            self.assertEqual(alert['diaSource'][field], deserialized['diaSource'][field])
        self.assertEqual(1, deserialized["alertId"])


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
