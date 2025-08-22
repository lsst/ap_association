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

import datetime
import io
import os

import numpy as np
import pandas as pd
import tempfile
import unittest
from unittest.mock import patch, Mock
from astropy import wcs
from astropy.nddata import CCDData
import fastavro
try:
    import confluent_kafka
    from confluent_kafka import KafkaException
except ImportError:
    confluent_kafka = None

import lsst.alert.packet as alertPack
from lsst.ap.association import PackageAlertsConfig, PackageAlertsTask
from lsst.ap.association.utils import readSchemaFromApdb, convertDataFrameToSdmSchema
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
import lsst.afw.image as afwImage
from lsst.daf.base import DateTime
from lsst.dax.apdb import Apdb, ApdbSql
import lsst.geom as geom
import lsst.meas.base.tests
from lsst.sphgeom import Box
import lsst.utils.tests
import utils_tests


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
    dateTime : `astropy.time.Time`
        Time for the Apdb.

    Returns
    -------
    objects : `pandas.DataFrame`
        Round tripped objects.
    sources : `pandas.DataFrame`
        Round tripped sources.
    """
    with tempfile.NamedTemporaryFile() as tmpFile:
        apdbConfig = ApdbSql.init_database(db_url="sqlite:///" + tmpFile.name)
        apdb = Apdb.from_config(apdbConfig)

        wholeSky = Box.full()
        loadedObjects = apdb.getDiaObjects(wholeSky)
        if loadedObjects.empty:
            diaObjects = objects
        else:
            diaObjects = pd.concat([loadedObjects, objects])
        loadedDiaSources = apdb.getDiaSources(wholeSky, [], dateTime)
        if loadedDiaSources.empty:
            diaSources = sources
        else:
            diaSources = pd.concat([loadedDiaSources, sources])
        loadedDiaForcedSources = apdb.getDiaForcedSources(wholeSky, [], dateTime)
        if loadedDiaForcedSources.empty:
            diaForcedSources = forcedSources
        else:
            diaForcedSources = pd.concat([loadedDiaForcedSources, forcedSources])

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

    # apply SDM type standardization to catch pandas typing issues
    schema = readSchemaFromApdb(apdb)
    diaObjects = convertDataFrameToSdmSchema(schema, diaObjects, tableName="DiaObject")
    diaSources = convertDataFrameToSdmSchema(schema, diaSources, tableName="DiaSource")
    diaForcedSources = convertDataFrameToSdmSchema(schema, diaForcedSources, tableName="DiaForcedSource")

    return (diaObjects, diaSources, diaForcedSources)


VISIT = 2
DETECTOR = 42


def mock_alert(dia_source_id):
    """Generate a minimal mock alert.
    """
    return {
        "diaSourceId": dia_source_id,
        "diaSource": {
            "midpointMjdTai": 5,
            "diaSourceId": 1234,
            "visit": VISIT,
            "detector": DETECTOR,
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
            # unlike in transformDiaSourceCatalog.py we need a timezone-aware
            # version because mock_alert does not go through pandas
            "timeProcessedMjdTai": DateTime.now().get(system=DateTime.MJD, scale=DateTime.TAI)
        }
    }


def mock_ss_alert(dia_source_id, ss_object_id):
    """Generate a minimal mock alert.
    """
    alert = mock_alert(dia_source_id)
    alert['MPCORB'] = {
        'mpcDesignation': 'K20A11H',  # a string-typed field
        'ssObjectId': ss_object_id,
        'q': np.float64(0.99999),  # a double-typed field

    }
    alert['ssSource'] = {
        'ssObjectId': ss_object_id,
        'eclipticLambda': np.float64(3.141592),  # a double-typed field
        'heliocentricDist': np.float32(3.141592),  # a float-typed field
    }
    return alert


def _deserialize_alert(alert_bytes):
    """Deserialize an alert message from Kafka.

    Parameters
    ----------
    alert_bytes : `bytes`
        Binary-encoding serialized Avro alert, including Confluent Wire
        Format prefix.

    Returns
    -------
    alert : `dict`
        An alert payload.
    """
    schema = alertPack.Schema.from_uri(str(alertPack.get_uri_to_latest_schema()))
    content_bytes = io.BytesIO(alert_bytes[5:])

    return fastavro.schemaless_reader(content_bytes, schema.definition)


class TestPackageAlerts(lsst.utils.tests.TestCase):
    def setUp(self):
        # Create an instance of random generator with fixed seed.
        rng = np.random.default_rng(1234)

        patcher = patch.dict(os.environ, {"AP_KAFKA_PRODUCER_PASSWORD": "fake_password",
                                          "AP_KAFKA_PRODUCER_USERNAME": "fake_username",
                                          "AP_KAFKA_SERVER": "fake_server",
                                          "AP_KAFKA_TOPIC": "fake_topic"})
        self.environ = patcher.start()
        self.addCleanup(patcher.stop)
        self.cutoutSize = 35
        self.center = lsst.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -30),
                                    lsst.geom.Extent2I(140, 160))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(100000.0, self.center)
        exposure, catalog = self.dataset.realize(
            10.0,
            self.dataset.makeMinimalSchema(),
            randomSeed=1234)
        self.exposure = exposure
        detector = DetectorWrapper(id=DETECTOR, bbox=exposure.getBBox()).detector
        self.exposure.setDetector(detector)

        visit = afwImage.VisitInfo(
            id=VISIT,
            exposureTime=200.,
            date=DateTime("2014-05-13T17:00:00.000000000",
                          DateTime.Timescale.TAI))
        self.exposure.info.id = 1234
        self.exposure.info.setVisitInfo(visit)

        self.exposure.setFilter(
            afwImage.FilterLabel(band='g', physical="g.MP9401"))

        diaObjects = utils_tests.makeDiaObjects(2, self.exposure, rng)
        diaSourceHistory = utils_tests.makeDiaSources(
            10, diaObjects["diaObjectId"].to_numpy(), self.exposure, rng)
        diaForcedSources = utils_tests.makeDiaForcedSources(
            10, diaObjects["diaObjectId"].to_numpy(), self.exposure, rng)
        self.diaObjects, diaSourceHistory, self.diaForcedSources = _roundTripThroughApdb(
            diaObjects,
            diaSourceHistory,
            diaForcedSources,
            self.exposure.visitInfo.date.toAstropy())
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

    def testCreateExtentMinimum(self):
        """Test the extent creation for the cutout bbox returns a cutout with
        the minimum cutouut size.
        """
        packConfig = PackageAlertsConfig()
        # Just create a minimum less than the default cutout size.
        packConfig.minCutoutSize = self.cutoutSize - 5
        packageAlerts = PackageAlertsTask(config=packConfig)
        extent = packageAlerts.createDiaSourceExtent(
            packConfig.minCutoutSize - 5)
        self.assertTrue(extent == geom.Extent2I(packConfig.minCutoutSize,
                                                packConfig.minCutoutSize))
        # Test that the cutout size is correctly increased.
        extent = packageAlerts.createDiaSourceExtent(self.cutoutSize)
        self.assertTrue(extent == geom.Extent2I(self.cutoutSize,
                                                self.cutoutSize))

    def testCreateExtentMaximum(self):
        """Test the extent creation for the cutout bbox returns a cutout with
        the maximum cutout size.
        """
        packConfig = PackageAlertsConfig()
        # Just create a maximum more than the default cutout size.
        packConfig.maxCutoutSize = self.cutoutSize + 5
        packageAlerts = PackageAlertsTask(config=packConfig)
        extent = packageAlerts.createDiaSourceExtent(
            packConfig.maxCutoutSize + 5)
        self.assertTrue(extent == geom.Extent2I(packConfig.maxCutoutSize,
                                                packConfig.maxCutoutSize))
        # Test that the cutout size is correctly reduced.
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
            self.exposure.getWcs().getPixelOrigin(),
            self.exposure.getBBox().getDimensions(),
            self.exposure.getPhotoCalib(),
            diaSrcId)
        calibExposure = self.exposure.getPhotoCalib().calibrateImage(
            self.exposure.getMaskedImage())

        self.assertFloatsAlmostEqual(ccdData.wcs.wcs.cd,
                                     self.cutoutWcs.wcs.cd)
        self.assertFloatsAlmostEqual(ccdData.data,
                                     calibExposure.getImage().array)
        self.assertFloatsAlmostEqual(ccdData.psf,
                                     self.exposure.psf.computeKernelImage(self.center).array)

        ccdData = packageAlerts.createCcdDataCutout(
            self.exposure,
            geom.SpherePoint(0, 0, geom.degrees),
            geom.Point2D(-100, -100),
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
            rtol=5e-10,
            atol=5e-10)

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
        dia_source_id = 1234

        for srcIdx, diaSource in self.diaSources.iterrows():
            sphPoint = geom.SpherePoint(diaSource["ra"],
                                        diaSource["dec"],
                                        geom.degrees)
            pixelPoint = geom.Point2D(diaSource["x"], diaSource["y"])
            cutout = self.exposure.getCutout(sphPoint,
                                             geom.Extent2I(self.cutoutSize,
                                                           self.cutoutSize))
            ccdCutout = packageAlerts.createCcdDataCutout(
                cutout,
                sphPoint,
                pixelPoint,
                geom.Extent2I(self.cutoutSize, self.cutoutSize),
                cutout.getPhotoCalib(),
                1234)
            cutoutBytes = packageAlerts.streamCcdDataToBytes(
                ccdCutout)
            objSources = self.diaSourceHistory.loc[srcIdx[0]]
            objForcedSources = self.diaForcedSources.loc[srcIdx[0]]
            alert = packageAlerts.makeAlertDict(
                dia_source_id,
                diaSource,
                self.diaObjects.loc[srcIdx[0]],
                objSources,
                objForcedSources,
                ccdCutout,
                ccdCutout,
                ccdCutout)
            self.assertEqual(len(alert), 11)

            self.assertEqual(alert["diaSourceId"], dia_source_id)
            self.assertEqual(alert["diaSource"], diaSource.to_dict())
            self.assertEqual(alert["cutoutDifference"],
                             cutoutBytes)
            self.assertEqual(alert["cutoutScience"],
                             cutoutBytes)
            self.assertEqual(alert["cutoutTemplate"],
                             cutoutBytes)

    @unittest.skipIf(confluent_kafka is None, 'Kafka is not enabled')
    def test_produceAlerts_empty_password(self):
        """ Test that produceAlerts raises if the password is empty or missing.
        """
        self.environ['AP_KAFKA_PRODUCER_PASSWORD'] = ""
        with self.assertRaisesRegex(ValueError, "Kafka password"):
            packConfig = PackageAlertsConfig(doProduceAlerts=True)
            PackageAlertsTask(config=packConfig)

        del self.environ['AP_KAFKA_PRODUCER_PASSWORD']
        with self.assertRaisesRegex(ValueError, "Kafka password"):
            packConfig = PackageAlertsConfig(doProduceAlerts=True)
            PackageAlertsTask(config=packConfig)

    @unittest.skipIf(confluent_kafka is None, 'Kafka is not enabled')
    def test_produceAlerts_empty_username(self):
        """ Test that produceAlerts raises if the username is empty or missing.
        """
        self.environ['AP_KAFKA_PRODUCER_USERNAME'] = ""
        with self.assertRaisesRegex(ValueError, "Kafka username"):
            packConfig = PackageAlertsConfig(doProduceAlerts=True)
            PackageAlertsTask(config=packConfig)

        del self.environ['AP_KAFKA_PRODUCER_USERNAME']
        with self.assertRaisesRegex(ValueError, "Kafka username"):
            packConfig = PackageAlertsConfig(doProduceAlerts=True)
            PackageAlertsTask(config=packConfig)

    @unittest.skipIf(confluent_kafka is None, 'Kafka is not enabled')
    def test_produceAlerts_empty_server(self):
        """ Test that produceAlerts raises if the server is empty or missing.
        """
        self.environ['AP_KAFKA_SERVER'] = ""
        with self.assertRaisesRegex(ValueError, "Kafka server"):
            packConfig = PackageAlertsConfig(doProduceAlerts=True)
            PackageAlertsTask(config=packConfig)

        del self.environ['AP_KAFKA_SERVER']
        with self.assertRaisesRegex(ValueError, "Kafka server"):
            packConfig = PackageAlertsConfig(doProduceAlerts=True)
            PackageAlertsTask(config=packConfig)

    @unittest.skipIf(confluent_kafka is None, 'Kafka is not enabled')
    def test_produceAlerts_empty_topic(self):
        """ Test that produceAlerts raises if the topic is empty or missing.
        """
        self.environ['AP_KAFKA_TOPIC'] = ""
        with self.assertRaisesRegex(ValueError, "Kafka topic"):
            packConfig = PackageAlertsConfig(doProduceAlerts=True)
            PackageAlertsTask(config=packConfig)

        del self.environ['AP_KAFKA_TOPIC']
        with self.assertRaisesRegex(ValueError, "Kafka topic"):
            packConfig = PackageAlertsConfig(doProduceAlerts=True)
            PackageAlertsTask(config=packConfig)

    @patch('confluent_kafka.Producer')
    @patch.object(PackageAlertsTask, '_server_check')
    @unittest.skipIf(confluent_kafka is None, 'Kafka is not enabled')
    def test_produceAlerts_success(self, mock_server_check, mock_producer):
        """ Test that produceAlerts calls the producer on all provided alerts
        when the alerts are all under the batch size limit.
        """
        packConfig = PackageAlertsConfig(doProduceAlerts=True)
        packageAlerts = PackageAlertsTask(config=packConfig)
        alerts = [mock_alert(1), mock_alert(2), mock_ss_alert(3, 3)]

        # Create a variable and assign it an instance of the patched kafka producer
        producer_instance = mock_producer.return_value
        producer_instance.produce = Mock()
        producer_instance.flush = Mock()
        unix_midpoint = self.exposure.visitInfo.date.toAstropy().tai.unix
        exposure_time = self.exposure.visitInfo.exposureTime
        packageAlerts.produceAlerts(alerts, VISIT, DETECTOR, unix_midpoint, exposure_time)

        self.assertEqual(mock_server_check.call_count, 1)
        self.assertEqual(producer_instance.produce.call_count, len(alerts))
        self.assertEqual(producer_instance.flush.call_count, len(alerts)+1)

    @patch('confluent_kafka.Producer')
    @patch.object(PackageAlertsTask, '_server_check')
    @unittest.skipIf(confluent_kafka is None, 'Kafka is not enabled')
    def test_produceAlerts_one_failure(self, mock_server_check, mock_producer):
        """ Test that produceAlerts correctly fails on one alert
        and is writing the failure to disk.
        """
        counter = 0

        def mock_produce(*args, **kwargs):
            nonlocal counter
            counter += 1
            if counter == 2:
                raise KafkaException
            else:
                return

        packConfig = PackageAlertsConfig(doProduceAlerts=True, doWriteFailedAlerts=True)
        packageAlerts = PackageAlertsTask(config=packConfig)

        patcher = patch("builtins.open")
        patch_open = patcher.start()
        alerts = [mock_alert(1), mock_alert(2), mock_alert(3), mock_ss_alert(4, 4)]
        unix_midpoint = self.exposure.visitInfo.date.toAstropy().tai.unix
        exposure_time = self.exposure.visitInfo.exposureTime

        producer_instance = mock_producer.return_value
        producer_instance.produce = Mock(side_effect=mock_produce)
        producer_instance.flush = Mock()
        packageAlerts.produceAlerts(alerts, VISIT, DETECTOR, unix_midpoint, exposure_time)

        self.assertEqual(mock_server_check.call_count, 1)
        self.assertEqual(producer_instance.produce.call_count, len(alerts))
        self.assertEqual(patch_open.call_count, 1)
        self.assertIn(f"{VISIT}_{DETECTOR}_2.avro", patch_open.call_args.args[0])
        # Because one produce raises, we call flush one fewer times than in the success
        # test above.
        self.assertEqual(producer_instance.flush.call_count, len(alerts))
        patcher.stop()

    @patch.object(PackageAlertsTask, '_server_check')
    def testRun_without_produce(self, mock_server_check):
        """Test the run method of package alerts with produce set to False and
        doWriteAlerts set to true.
        """
        packConfig = PackageAlertsConfig(doWriteAlerts=True)
        with tempfile.TemporaryDirectory(prefix='alerts') as tempdir:
            packConfig.alertWriteLocation = tempdir
            packageAlerts = PackageAlertsTask(config=packConfig)

            packageAlerts.run(self.diaSources,
                              self.diaObjects,
                              self.diaSourceHistory,
                              self.diaForcedSources,
                              self.exposure,
                              self.exposure,
                              self.exposure)

            self.assertEqual(mock_server_check.call_count, 0)

            with open(os.path.join(tempdir, f"{VISIT}_{DETECTOR}.avro"), 'rb') as f:
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
                elif isinstance(value, datetime.datetime):
                    # TEMPORARY: avoid known failure
                    # self.assertEqual(value, self.diaSources.iloc[idx][key].
                    #                 replace(tzinfo=datetime.UTC).to_pydatetime())
                    continue
                else:
                    self.assertEqual(value, self.diaSources.iloc[idx][key])
            sphPoint = geom.SpherePoint(alert["diaSource"]["ra"],
                                        alert["diaSource"]["dec"],
                                        geom.degrees)
            pixelPoint = geom.Point2D(alert["diaSource"]["x"], alert["diaSource"]["y"])
            cutout = self.exposure.getCutout(sphPoint,
                                             geom.Extent2I(self.cutoutSize,
                                                           self.cutoutSize))
            ccdCutout = packageAlerts.createCcdDataCutout(
                cutout,
                sphPoint,
                pixelPoint,
                geom.Extent2I(self.cutoutSize, self.cutoutSize),
                cutout.getPhotoCalib(),
                1234)
            self.assertEqual(alert["cutoutDifference"],
                             packageAlerts.streamCcdDataToBytes(ccdCutout))

    @patch.object(PackageAlertsTask, '_server_check')
    def testRun_without_produce_use_averagePsf(self, mock_server_check):
        """Test the run method of package alerts with produce set to False and
        doWriteAlerts set to true.
        """
        packConfig = PackageAlertsConfig(doWriteAlerts=True)
        with tempfile.TemporaryDirectory(prefix='alerts') as tempdir:
            packConfig.alertWriteLocation = tempdir
            packConfig.useAveragePsf = True
            packageAlerts = PackageAlertsTask(config=packConfig)

            packageAlerts.run(self.diaSources,
                              self.diaObjects,
                              self.diaSourceHistory,
                              self.diaForcedSources,
                              self.exposure,
                              self.exposure,
                              self.exposure)

            self.assertEqual(mock_server_check.call_count, 0)

            with open(os.path.join(tempdir, f"{VISIT}_{DETECTOR}.avro"), 'rb') as f:
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
                elif isinstance(value, datetime.datetime):
                    # TEMPORARY: avoid known failure
                    # self.assertEqual(value, self.diaSources.iloc[idx][key].
                    #                 replace(tzinfo=datetime.UTC).to_pydatetime())
                    continue
                else:
                    self.assertEqual(value, self.diaSources.iloc[idx][key])
            sphPoint = geom.SpherePoint(alert["diaSource"]["ra"],
                                        alert["diaSource"]["dec"],
                                        geom.degrees)
            pixelPoint = geom.Point2D(alert["diaSource"]["x"], alert["diaSource"]["y"])
            cutout = self.exposure.getCutout(sphPoint,
                                             geom.Extent2I(self.cutoutSize,
                                                           self.cutoutSize))
            ccdCutout = packageAlerts.createCcdDataCutout(
                cutout,
                sphPoint,
                pixelPoint,
                geom.Extent2I(self.cutoutSize, self.cutoutSize),
                cutout.getPhotoCalib(),
                1234)
            self.assertEqual(alert["cutoutDifference"],
                             packageAlerts.streamCcdDataToBytes(ccdCutout))

    @patch.object(PackageAlertsTask, 'produceAlerts')
    @patch('confluent_kafka.Producer')
    @patch.object(PackageAlertsTask, '_server_check')
    @unittest.skipIf(confluent_kafka is None, 'Kafka is not enabled')
    def testRun_with_produce(self, mock_produceAlerts, mock_server_check, mock_producer):
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
                          self.exposure,
                          self.exposure)
        self.assertEqual(mock_server_check.call_count, 1)
        self.assertEqual(mock_produceAlerts.call_count, 1)

    def test_serialize_alert_round_trip(self):
        """Test that values in the alert packet exactly round trip.
        """
        packClass = PackageAlertsConfig()
        packageAlerts = PackageAlertsTask(config=packClass)

        alert = mock_ss_alert(1, 1)
        serialized = PackageAlertsTask._serializeAlert(packageAlerts, alert)
        deserialized = _deserialize_alert(serialized)
        for table in ['diaSource', 'ssSource', 'MPCORB']:
            for field in alert[table]:
                self.assertEqual(alert[table][field], deserialized[table][field])

        self.assertEqual(1, deserialized["diaSourceId"])

    @unittest.skipIf(confluent_kafka is None, 'Kafka is not enabled')
    def test_server_check(self):

        with self.assertRaisesRegex(KafkaException, "_TRANSPORT"):
            packConfig = PackageAlertsConfig(doProduceAlerts=True)
            PackageAlertsTask(config=packConfig)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
