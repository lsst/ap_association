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

__all__ = ("PackageAlertsConfig", "PackageAlertsTask")

import io
import os
import sys

from astropy import wcs
import astropy.units as u
from astropy.nddata import CCDData, VarianceUncertainty
import pandas as pd
import struct
import fastavro
# confluent_kafka is not in the standard Rubin environment as it is a third
# party package and is only needed when producing alerts.
try:
    import confluent_kafka
    from confluent_kafka import KafkaException
except ImportError:
    confluent_kafka = None

import lsst.alert.packet as alertPack
import lsst.afw.geom as afwGeom
import lsst.geom as geom
import lsst.pex.config as pexConfig
from lsst.pex.exceptions import InvalidParameterError
import lsst.pipe.base as pipeBase
from lsst.utils.timer import timeMethod


class PackageAlertsConfig(pexConfig.Config):
    """Config class for AssociationTask.
    """
    schemaFile = pexConfig.Field(
        dtype=str,
        doc="Schema definition file URI for the avro alerts.",
        default=str(alertPack.get_uri_to_latest_schema())
    )
    minCutoutSize = pexConfig.RangeField(
        dtype=int,
        min=0,
        max=1000,
        default=30,
        doc="Dimension of the square image cutouts to package in the alert."
    )
    alertWriteLocation = pexConfig.Field(
        dtype=str,
        doc="Location to write alerts to.",
        default=os.path.join(os.getcwd(), "alerts"),
    )

    doProduceAlerts = pexConfig.Field(
        dtype=bool,
        doc="Turn on alert production to kafka if true and if confluent_kafka is in the environment.",
        default=False,
    )

    doWriteAlerts = pexConfig.Field(
        dtype=bool,
        doc="Write alerts to disk if true.",
        default=False,
    )


class PackageAlertsTask(pipeBase.Task):
    """Tasks for packaging Dia and Pipelines data into Avro alert packages.
    """
    ConfigClass = PackageAlertsConfig
    _DefaultName = "packageAlerts"

    _scale = (1.0 * geom.arcseconds).asDegrees()

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.alertSchema = alertPack.Schema.from_uri(self.config.schemaFile)
        os.makedirs(self.config.alertWriteLocation, exist_ok=True)

        if self.config.doProduceAlerts:
            if confluent_kafka is not None:
                self.password = os.getenv("AP_KAFKA_PRODUCER_PASSWORD")
                if not self.password:
                    raise ValueError("Kafka password environment variable was not set.")
                self.username = os.getenv("AP_KAFKA_PRODUCER_USERNAME")
                if not self.username:
                    raise ValueError("Kafka username environment variable was not set.")
                self.server = os.getenv("AP_KAFKA_SERVER")
                if not self.server:
                    raise ValueError("Kafka server environment variable was not set.")
                self.kafkaTopic = os.getenv("AP_KAFKA_TOPIC")
                if not self.kafkaTopic:
                    raise ValueError("Kafka topic environment variable was not set.")

                # confluent_kafka configures all of its classes with dictionaries. This one
                # sets up the bare minimum that is needed.
                self.kafkaConfig = {
                    # This is the URL to use to connect to the Kafka cluster.
                    "bootstrap.servers": self.server,
                    # These next two properties tell the Kafka client about the specific
                    # authentication and authorization protocols that should be used when
                    # connecting.
                    "security.protocol": "SASL_PLAINTEXT",
                    "sasl.mechanisms": "SCRAM-SHA-512",
                    # The sasl.username and sasl.password are passed through over
                    # SCRAM-SHA-512 auth to connect to the cluster. The username is not
                    # sensitive, but the password is (of course) a secret value which
                    # should never be committed to source code.
                    "sasl.username": self.username,
                    "sasl.password": self.password,
                    # Batch size limits the largest size of a kafka alert that can be sent.
                    # We set the batch size to 2 Mb.
                    "batch.size": 2097152,
                    "linger.ms": 5,
                }
                self.producer = confluent_kafka.Producer(**self.kafkaConfig)

            else:
                raise RuntimeError("Produce alerts is set but confluent_kafka is not present in "
                                   "the environment. Alerts will not be sent to the alert stream.")

    @timeMethod
    def run(self,
            diaSourceCat,
            diaObjectCat,
            diaSrcHistory,
            diaForcedSources,
            diffIm,
            calexp,
            template,
            ):
        """Package DiaSources/Object and exposure data into Avro alerts.

        Alerts can be sent to the alert stream if ``doProduceAlerts`` is set
        and written to disk if ``doWriteAlerts`` is set. Both can be set at the
        same time, and are independent of one another.

        Writes Avro alerts to a location determined by the
        ``alertWriteLocation`` configurable.

        Parameters
        ----------
        diaSourceCat : `pandas.DataFrame`
            New DiaSources to package. DataFrame should be indexed on
            ``["diaObjectId", "band", "diaSourceId"]``
        diaObjectCat : `pandas.DataFrame`
            New and updated DiaObjects matched to the new DiaSources. DataFrame
            is indexed on ``["diaObjectId"]``
        diaSrcHistory : `pandas.DataFrame`
            12 month history of DiaSources matched to the DiaObjects. Excludes
            the newest DiaSource and is indexed on
            ``["diaObjectId", "band", "diaSourceId"]``
        diaForcedSources : `pandas.DataFrame`
            12 month history of DiaForcedSources matched to the DiaObjects.
            ``["diaObjectId"]``
        diffIm : `lsst.afw.image.ExposureF`
            Difference image the sources in ``diaSourceCat`` were detected in.
        calexp : `lsst.afw.image.ExposureF`
            Calexp used to create the ``diffIm``.
        template : `lsst.afw.image.ExposureF` or `None`
            Template image used to create the ``diffIm``.
        """
        alerts = []
        self._patchDiaSources(diaSourceCat)
        self._patchDiaSources(diaSrcHistory)
        ccdVisitId = diffIm.info.id
        diffImPhotoCalib = diffIm.getPhotoCalib()
        calexpPhotoCalib = calexp.getPhotoCalib()
        templatePhotoCalib = template.getPhotoCalib()
        for srcIndex, diaSource in diaSourceCat.iterrows():
            # Get all diaSources for the associated diaObject.
            # TODO: DM-31992 skip DiaSources associated with Solar System
            #       Objects for now.
            if srcIndex[0] == 0:
                continue
            diaObject = diaObjectCat.loc[srcIndex[0]]
            if diaObject["nDiaSources"] > 1:
                objSourceHistory = diaSrcHistory.loc[srcIndex[0]]
            else:
                objSourceHistory = None
            objDiaForcedSources = diaForcedSources.loc[srcIndex[0]]
            sphPoint = geom.SpherePoint(diaSource["ra"],
                                        diaSource["dec"],
                                        geom.degrees)

            cutoutExtent = self.createDiaSourceExtent(diaSource["bboxSize"])
            diffImCutout = self.createCcdDataCutout(
                diffIm,
                sphPoint,
                cutoutExtent,
                diffImPhotoCalib,
                diaSource["diaSourceId"])
            calexpCutout = self.createCcdDataCutout(
                calexp,
                sphPoint,
                cutoutExtent,
                calexpPhotoCalib,
                diaSource["diaSourceId"])
            templateCutout = self.createCcdDataCutout(
                template,
                sphPoint,
                cutoutExtent,
                templatePhotoCalib,
                diaSource["diaSourceId"])

            # TODO: Create alertIds DM-24858
            alertId = diaSource["diaSourceId"]
            alerts.append(
                self.makeAlertDict(alertId,
                                   diaSource,
                                   diaObject,
                                   objSourceHistory,
                                   objDiaForcedSources,
                                   diffImCutout,
                                   calexpCutout,
                                   templateCutout))

        if self.config.doProduceAlerts:
            self.produceAlerts(alerts, ccdVisitId)

        if self.config.doWriteAlerts:
            with open(os.path.join(self.config.alertWriteLocation, f"{ccdVisitId}.avro"), "wb") as f:
                self.alertSchema.store_alerts(f, alerts)

    def _patchDiaSources(self, diaSources):
        """Add the ``programId`` column to the data.

        Parameters
        ----------
        diaSources : `pandas.DataFrame`
            DataFrame of DiaSources to patch.
        """
        diaSources["programId"] = 0

    def createDiaSourceExtent(self, bboxSize):
        """Create an extent for a box for the cutouts given the size of the
        square BBox that covers the source footprint.

        Parameters
        ----------
        bboxSize : `int`
            Size of a side of the square bounding box in pixels.

        Returns
        -------
        extent : `lsst.geom.Extent2I`
            Geom object representing the size of the bounding box.
        """
        if bboxSize < self.config.minCutoutSize:
            extent = geom.Extent2I(self.config.minCutoutSize,
                                   self.config.minCutoutSize)
        else:
            extent = geom.Extent2I(bboxSize, bboxSize)
        return extent

    def produceAlerts(self, alerts, ccdVisitId):
        """Serialize alerts and send them to the alert stream using
        confluent_kafka's producer.

        Parameters
        ----------
        alerts : `dict`
            Dictionary of alerts to be sent to the alert stream.
        ccdVisitId : `int`
            ccdVisitId of the alerts sent to the alert stream. Used to write
            out alerts which fail to be sent to the alert stream.
        """
        for alert in alerts:
            alertBytes = self._serializeAlert(alert, schema=self.alertSchema.definition, schema_id=1)
            try:
                self.producer.produce(self.kafkaTopic, alertBytes, callback=self._delivery_callback)
                self.producer.flush()

            except KafkaException as e:
                self.log.warning('Kafka error: {}, message was {} bytes'.format(e, sys.getsizeof(alertBytes)))

                with open(os.path.join(self.config.alertWriteLocation,
                                       f"{ccdVisitId}_{alert['alertId']}.avro"), "wb") as f:
                    f.write(alertBytes)

        self.producer.flush()

    def createCcdDataCutout(self, image, skyCenter, extent, photoCalib, srcId):
        """Grab an image as a cutout and return a calibrated CCDData image.

        Parameters
        ----------
        image : `lsst.afw.image.ExposureF`
            Image to pull cutout from.
        skyCenter : `lsst.geom.SpherePoint`
            Center point of DiaSource on the sky.
        extent : `lsst.geom.Extent2I`
            Bounding box to cutout from the image.
        photoCalib : `lsst.afw.image.PhotoCalib`
            Calibrate object of the image the cutout is cut from.
        srcId : `int`
            Unique id of DiaSource. Used for when an error occurs extracting
            a cutout.

        Returns
        -------
        ccdData : `astropy.nddata.CCDData` or `None`
            CCDData object storing the calibrate information from the input
            difference or template image.
        """
        # Catch errors in retrieving the cutout.
        try:
            cutout = image.getCutout(skyCenter, extent)
        except InvalidParameterError:
            point = image.getWcs().skyToPixel(skyCenter)
            imBBox = image.getBBox()
            if not geom.Box2D(image.getBBox()).contains(point):
                self.log.warning(
                    "DiaSource id=%i centroid lies at pixel (%.2f, %.2f) "
                    "which is outside the Exposure with bounding box "
                    "((%i, %i), (%i, %i)). Returning None for cutout...",
                    srcId, point.x, point.y,
                    imBBox.minX, imBBox.maxX, imBBox.minY, imBBox.maxY)
            else:
                raise InvalidParameterError(
                    "Failed to retrieve cutout from image for DiaSource with "
                    "id=%i. InvalidParameterError thrown during cutout "
                    "creation. Exiting."
                    % srcId)
            return None

        # Find the value of the bottom corner of our cutout's BBox and
        # subtract 1 so that the CCDData cutout position value will be
        # [1, 1].
        cutOutMinX = cutout.getBBox().minX - 1
        cutOutMinY = cutout.getBBox().minY - 1
        center = cutout.getWcs().skyToPixel(skyCenter)
        calibCutout = photoCalib.calibrateImage(cutout.getMaskedImage())

        cutoutWcs = wcs.WCS(naxis=2)
        cutoutWcs.array_shape = (cutout.getBBox().getWidth(),
                                 cutout.getBBox().getWidth())
        cutoutWcs.wcs.crpix = [center.x - cutOutMinX, center.y - cutOutMinY]
        cutoutWcs.wcs.crval = [skyCenter.getRa().asDegrees(),
                               skyCenter.getDec().asDegrees()]
        cutoutWcs.wcs.cd = self.makeLocalTransformMatrix(cutout.getWcs(),
                                                         center,
                                                         skyCenter)

        return CCDData(
            data=calibCutout.getImage().array,
            uncertainty=VarianceUncertainty(calibCutout.getVariance().array),
            flags=calibCutout.getMask().array,
            wcs=cutoutWcs,
            meta={"cutMinX": cutOutMinX,
                  "cutMinY": cutOutMinY},
            unit=u.nJy)

    def makeLocalTransformMatrix(self, wcs, center, skyCenter):
        """Create a local, linear approximation of the wcs transformation
        matrix.

        The approximation is created as if the center is at RA=0, DEC=0. All
        comparing x,y coordinate are relative to the position of center. Matrix
        is initially calculated with units arcseconds and then converted to
        degrees. This yields higher precision results due to quirks in AST.

        Parameters
        ----------
        wcs : `lsst.afw.geom.SkyWcs`
            Wcs to approximate
        center : `lsst.geom.Point2D`
            Point at which to evaluate the LocalWcs.
        skyCenter : `lsst.geom.SpherePoint`
            Point on sky to approximate the Wcs.

        Returns
        -------
        localMatrix : `numpy.ndarray`
            Matrix representation the local wcs approximation with units
            degrees.
        """
        blankCDMatrix = [[self._scale, 0], [0, self._scale]]
        localGnomonicWcs = afwGeom.makeSkyWcs(
            center, skyCenter, blankCDMatrix)
        measurementToLocalGnomonic = wcs.getTransform().then(
            localGnomonicWcs.getTransform().inverted()
        )
        localMatrix = measurementToLocalGnomonic.getJacobian(center)
        return localMatrix / 3600

    def makeAlertDict(self,
                      alertId,
                      diaSource,
                      diaObject,
                      objDiaSrcHistory,
                      objDiaForcedSources,
                      diffImCutout,
                      calexpCutout,
                      templateCutout):
        """Convert data and package into a dictionary alert.

        Parameters
        ----------
        diaSource : `pandas.DataFrame`
            New single DiaSource to package.
        diaObject : `pandas.DataFrame`
            DiaObject that ``diaSource`` is matched to.
        objDiaSrcHistory : `pandas.DataFrame`
            12 month history of ``diaObject`` excluding the latest DiaSource.
        objDiaForcedSources : `pandas.DataFrame`
            12 month history of ``diaObject`` forced measurements.
        diffImCutout : `astropy.nddata.CCDData` or `None`
            Cutout of the difference image around the location of ``diaSource``
            with a min size set by the ``cutoutSize`` configurable.
        calexpCutout : `astropy.nddata.CCDData` or `None`
            Cutout of the calexp around the location of ``diaSource``
            with a min size set by the ``cutoutSize`` configurable.
        templateCutout : `astropy.nddata.CCDData` or `None`
            Cutout of the template image around the location of ``diaSource``
            with a min size set by the ``cutoutSize`` configurable.
        """
        alert = dict()
        alert['alertId'] = alertId
        alert['diaSource'] = diaSource.to_dict()

        if objDiaSrcHistory is None:
            alert['prvDiaSources'] = objDiaSrcHistory
        else:
            alert['prvDiaSources'] = objDiaSrcHistory.to_dict("records")

        if isinstance(objDiaForcedSources, pd.Series):
            alert['prvDiaForcedSources'] = [objDiaForcedSources.to_dict()]
        else:
            alert['prvDiaForcedSources'] = objDiaForcedSources.to_dict("records")
        alert['prvDiaNondetectionLimits'] = None

        alert['diaObject'] = diaObject.to_dict()

        alert['ssObject'] = None

        if diffImCutout is None:
            alert['cutoutDifference'] = None
        else:
            alert['cutoutDifference'] = self.streamCcdDataToBytes(diffImCutout)

        if calexpCutout is None:
            alert['cutoutScience'] = None
        else:
            alert['cutoutScience'] = self.streamCcdDataToBytes(calexpCutout)

        if templateCutout is None:
            alert["cutoutTemplate"] = None
        else:
            alert["cutoutTemplate"] = self.streamCcdDataToBytes(templateCutout)

        return alert

    def streamCcdDataToBytes(self, cutout):
        """Serialize a cutout into bytes.

        Parameters
        ----------
        cutout : `astropy.nddata.CCDData`
            Cutout to serialize.

        Returns
        -------
        coutputBytes : `bytes`
            Input cutout serialized into byte data.
        """
        with io.BytesIO() as streamer:
            cutout.write(streamer, format="fits")
            cutoutBytes = streamer.getvalue()
        return cutoutBytes

    def _serializeAlert(self, alert, schema=None, schema_id=0):
        """Serialize an alert to a byte sequence for sending to Kafka.

        Parameters
        ----------
        alert : `dict`
            An alert payload to be serialized.
        schema : `dict`, optional
            An Avro schema definition describing how to encode `alert`. By default,
            the schema is None, which sets it to the latest schema available.
        schema_id : `int`, optional
            The Confluent Schema Registry ID of the schema. By default, 0 (an
            invalid ID) is used, indicating that the schema is not registered.

        Returns
        -------
        serialized : `bytes`
            The byte sequence describing the alert, including the Confluent Wire
            Format prefix.
        """
        if schema is None:
            schema = self.alertSchema.definition

        buf = io.BytesIO()
        # TODO: Use a proper schema versioning system (DM-42606)
        buf.write(self._serializeConfluentWireHeader(schema_id))
        fastavro.schemaless_writer(buf, schema, alert)
        return buf.getvalue()

    @staticmethod
    def _serializeConfluentWireHeader(schema_version):
        """Returns the byte prefix for Confluent Wire Format-style Kafka messages.

        Parameters
        ----------
        schema_version : `int`
            A version number which indicates the Confluent Schema Registry ID
            number of the Avro schema used to encode the message that follows this
            header.

        Returns
        -------
        header : `bytes`
            The 5-byte encoded message prefix.

        Notes
        -----
        The Confluent Wire Format is described more fully here:
        https://docs.confluent.io/current/schema-registry/serdes-develop/index.html#wire-format
        """
        ConfluentWireFormatHeader = struct.Struct(">bi")
        return ConfluentWireFormatHeader.pack(0, schema_version)

    def _delivery_callback(self, err, msg):
        if err:
            self.log.warning('Message failed delivery: %s\n' % err)
        else:
            self.log.debug('Message delivered to %s [%d] @ %d', msg.topic(), msg.partition(), msg.offset())
