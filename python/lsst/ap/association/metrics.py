# This file is part of ap_association.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__all__ = ["NumberNewDiaObjectsMetricTask",
           "NumberUnassociatedDiaObjectsMetricTask",
           "FractionUpdatedDiaObjectsMetricTask",
           "NumberSolarSystemObjectsMetricTask",
           "NumberAssociatedSolarSystemObjectsMetricTask",
           "TotalUnassociatedDiaObjectsMetricTask",
           ]


import astropy.units as u

from lsst.verify import Measurement
from lsst.verify.tasks import MetadataMetricTask, MetadataMetricConfig, \
    ApdbMetricTask, ApdbMetricConfig, MetricComputationError


class NumberNewDiaObjectsMetricConfig(MetadataMetricConfig):
    def setDefaults(self):
        self.connections.package = "ap_association"
        self.connections.metric = "numNewDiaObjects"


class NumberNewDiaObjectsMetricTask(MetadataMetricTask):
    """Task that computes the number of DIASources that create new DIAObjects
    in an image, visit, etc..
    """
    _DefaultName = "numNewDiaObjects"
    ConfigClass = NumberNewDiaObjectsMetricConfig

    def makeMeasurement(self, values):
        """Compute the number of new DIAObjects.

        Parameters
        ----------
        values : `dict` [`str`, `int` or `None`]
            A `dict` representation of the metadata. Each `dict` has the
            following key:

            ``"newObjects"``
                The number of new objects created for this image (`int`
                or `None`). May be `None` if the image was not successfully
                associated.

        Returns
        -------
        measurement : `lsst.verify.Measurement` or `None`
            The total number of new objects.
        """
        if values["newObjects"] is not None:
            try:
                nNew = int(values["newObjects"])
            except (ValueError, TypeError) as e:
                raise MetricComputationError("Corrupted value of numNewDiaObjects") from e
            else:
                return Measurement(self.config.metricName, nNew * u.count)
        else:
            self.log.info("Nothing to do: no association results found.")
            return None

    @classmethod
    def getInputMetadataKeys(cls, config):
        return {"newObjects": ".numNewDiaObjects"}


class NumberUnassociatedDiaObjectsMetricConfig(MetadataMetricConfig):
    def setDefaults(self):
        self.connections.package = "ap_association"
        self.connections.metric = "numUnassociatedDiaObjects"


class NumberUnassociatedDiaObjectsMetricTask(MetadataMetricTask):
    """Task that computes the number of previously-known DIAObjects that do
    not have detected DIASources in an image, visit, etc..
    """
    _DefaultName = "numUnassociatedDiaObjects"
    ConfigClass = NumberUnassociatedDiaObjectsMetricConfig

    def makeMeasurement(self, values):
        """Compute the number of non-updated DIAObjects.

        Parameters
        ----------
        values : `dict` [`str`, `int` or `None`]
            A `dict` representation of the metadata. Each `dict` has the
            following key:

            ``"unassociatedObjects"``
                The number of DIAObjects not associated with a DiaSource in
                this image (`int` or `None`). May be `None` if the image was
                not successfully associated.

        Returns
        -------
        measurement : `lsst.verify.Measurement` or `None`
            The total number of unassociated objects.
        """
        if values["unassociatedObjects"] is not None:
            try:
                nNew = int(values["unassociatedObjects"])
            except (ValueError, TypeError) as e:
                raise MetricComputationError("Corrupted value of numUnassociatedDiaObjects") from e
            else:
                return Measurement(self.config.metricName, nNew * u.count)
        else:
            self.log.info("Nothing to do: no association results found.")
            return None

    @classmethod
    def getInputMetadataKeys(cls, config):
        return {"unassociatedObjects": ".numUnassociatedDiaObjects"}


class FractionUpdatedDiaObjectsMetricConfig(MetadataMetricConfig):
    def setDefaults(self):
        self.connections.package = "ap_association"
        self.connections.metric = "fracUpdatedDiaObjects"


class FractionUpdatedDiaObjectsMetricTask(MetadataMetricTask):
    """Task that computes the fraction of previously created DIAObjects that
    have a new association in this image, visit, etc..
    """
    _DefaultName = "fracUpdatedDiaObjects"
    ConfigClass = FractionUpdatedDiaObjectsMetricConfig

    def makeMeasurement(self, values):
        """Compute the number of non-updated DIAObjects.

        AssociationTask reports each pre-existing DIAObject as either updated
        (associated with a new DIASource) or unassociated.

        Parameters
        ----------
        values : `dict` [`str`, `int` or `None`]
            A `dict` representation of the metadata. Each `dict` has the
            following keys:

            ``"updatedObjects"``
                The number of DIAObjects updated for this image (`int` or
                `None`). May be `None` if the image was not
                successfully associated.
            ``"unassociatedObjects"``
                The number of DIAObjects not associated with a DiaSource in
                this image (`int` or `None`). May be `None` if the image was
                not successfully associated.

        Returns
        -------
        measurement : `lsst.verify.Measurement` or `None`
            The total number of unassociated objects.
        """
        if values["updatedObjects"] is not None \
                and values["unassociatedObjects"] is not None:
            try:
                nUpdated = int(values["updatedObjects"])
                nUnassociated = int(values["unassociatedObjects"])
            except (ValueError, TypeError) as e:
                raise MetricComputationError("Corrupted value of numUpdatedDiaObjects "
                                             "or numUnassociatedDiaObjects") from e
            else:
                if nUpdated <= 0 and nUnassociated <= 0:
                    return None  # No pre-existing DIAObjects; no fraction to compute
                else:
                    fraction = nUpdated / (nUpdated + nUnassociated)
                    return Measurement(self.config.metricName, fraction * u.dimensionless_unscaled)
        else:
            self.log.info("Nothing to do: no association results found.")
            return None

    @classmethod
    def getInputMetadataKeys(cls, config):
        return {"updatedObjects": ".numUpdatedDiaObjects",
                "unassociatedObjects": ".numUnassociatedDiaObjects"}


class NumberSolarSystemObjectsMetricConfig(MetadataMetricConfig):
    def setDefaults(self):
        self.connections.package = "ap_association"
        self.connections.metric = "numTotalSolarSystemObjects"


class NumberSolarSystemObjectsMetricTask(MetadataMetricTask):
    """Task that computes the number of SolarSystemObjects that are
    observable within this detectorVisit.
    """
    _DefaultName = "numTotalSolarSystemObjects"
    ConfigClass = NumberSolarSystemObjectsMetricConfig

    def makeMeasurement(self, values):
        """Compute the total number of SolarSystemObjects within a
        detectorVisit.

        Parameters
        ----------
        values : `dict` [`str`, `int` or `None`]
            A `dict` representation of the metadata. Each `dict` has the
            following key:

            ``"numTotalSolarSystemObjects"``
                The number of SolarSystemObjects within the observable detector
                area (`int` or `None`). May be `None` if solar system
                association was not attempted or the image was not
                successfully associated.

        Returns
        -------
        measurement : `lsst.verify.Measurement` or `None`
            The total number of Solar System objects.
        """
        if values["numTotalSolarSystemObjects"] is not None:
            try:
                nNew = int(values["numTotalSolarSystemObjects"])
            except (ValueError, TypeError) as e:
                raise MetricComputationError("Corrupted value of numTotalSolarSystemObjects") from e
            else:
                return Measurement(self.config.metricName, nNew * u.count)
        else:
            self.log.info("Nothing to do: no solar system results found.")
            return None

    @classmethod
    def getInputMetadataKeys(cls, config):
        return {"numTotalSolarSystemObjects": ".numTotalSolarSystemObjects"}


class NumberAssociatedSolarSystemObjectsMetricConfig(MetadataMetricConfig):
    def setDefaults(self):
        self.connections.package = "ap_association"
        self.connections.metric = "numAssociatedSsObjects"


class NumberAssociatedSolarSystemObjectsMetricTask(MetadataMetricTask):
    """Number of SolarSystemObjects that were associated with new DiaSources
    for this detectorVisit.
    """
    _DefaultName = "numAssociatedSsObjects"
    ConfigClass = NumberAssociatedSolarSystemObjectsMetricConfig

    def makeMeasurement(self, values):
        """Compute the number of associated SolarSystemObjects.

        Parameters
        ----------
        values : `dict` [`str`, `int` or `None`]
            A `dict` representation of the metadata. Each `dict` has the
            following key:

            ``"numAssociatedSsObjects"``
                The number of successfully associated SolarSystem Objects
                (`int` or `None`). May be `None` if solar system association
                was not attempted or the image was not successfully associated.

        Returns
        -------
        measurement : `lsst.verify.Measurement` or `None`
            The total number of associated SolarSystemObjects.
        """
        if values["numAssociatedSsObjects"] is not None:
            try:
                nNew = int(values["numAssociatedSsObjects"])
            except (ValueError, TypeError) as e:
                raise MetricComputationError("Corrupted value of numAssociatedSsObjects") from e
            else:
                return Measurement(self.config.metricName, nNew * u.count)
        else:
            self.log.info("Nothing to do: no solar system results found.")
            return None

    @classmethod
    def getInputMetadataKeys(cls, config):
        return {"numAssociatedSsObjects": ".numAssociatedSsObjects"}


class TotalUnassociatedDiaObjectsMetricConfig(ApdbMetricConfig):
    def setDefaults(self):
        self.connections.package = "ap_association"
        self.connections.metric = "totalUnassociatedDiaObjects"


class TotalUnassociatedDiaObjectsMetricTask(ApdbMetricTask):
    """Task that computes the number of DIAObjects with only one
    associated DIASource.
    """
    _DefaultName = "totalUnassociatedDiaObjects"
    ConfigClass = TotalUnassociatedDiaObjectsMetricConfig

    def makeMeasurement(self, dbHandle, outputDataId):
        """Compute the number of unassociated DIAObjects.

        Parameters
        ----------
        dbHandle : `lsst.dax.apdb.Apdb`
            A database instance.
        outputDataId : any data ID type
            The subset of the database to which this measurement applies.
            Must be empty, as the number of unassociated sources is
            ill-defined for subsets of the dataset.

        Returns
        -------
        measurement : `lsst.verify.Measurement`
            The total number of unassociated objects.

        Raises
        ------
        MetricComputationError
            Raised on any failure to query the database.
        ValueError
            Raised if outputDataId is not empty
        """
        # All data ID types define keys()
        if outputDataId.keys() - {'instrument'}:
            raise ValueError("%s must not be associated with specific data IDs (gave %s)."
                             % (self.config.metricName, outputDataId))

        try:
            nUnassociatedDiaObjects = dbHandle.countUnassociatedObjects()
        except Exception as e:
            raise MetricComputationError("Could not get unassociated objects from database") from e

        meas = Measurement(self.config.metricName, nUnassociatedDiaObjects * u.count)
        return meas
