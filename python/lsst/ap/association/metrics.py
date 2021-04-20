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
           "FractionAssociatedSourcesMetricTask",
           "TotalUnassociatedDiaObjectsMetricTask",
           ]


import astropy.units as u

from lsst.pipe.base import connectionTypes, Struct
from lsst.verify import Measurement
from lsst.verify.gen2tasks import register
from lsst.verify.tasks import MetadataMetricTask, MetadataMetricConfig, \
    ApdbMetricTask, ApdbMetricConfig, MetricComputationError, \
    AbstractMetadataMetricTask, SingleMetadataMetricConnections


class NumberNewDiaObjectsMetricConfig(MetadataMetricConfig):
    def setDefaults(self):
        self.connections.package = "ap_association"
        self.connections.metric = "numNewDiaObjects"


@register("numNewDiaObjects")
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


@register("numUnassociatedDiaObjects")
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


@register("fracUpdatedDiaObjects")
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


class TotalUnassociatedDiaObjectsMetricConfig(ApdbMetricConfig):
    def setDefaults(self):
        self.connections.package = "ap_association"
        self.connections.metric = "totalUnassociatedDiaObjects"


@register("totalUnassociatedDiaObjects")
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


class FractionAssociatedSourcesConnections(
        SingleMetadataMetricConnections,
        dimensions={"instrument", "visit", "detector"},
        defaultTemplates={"labelName": "diaPipe",
                          "package": "ap_association",
                          "metric": "fracAssociatedDiaSources",
                          "coaddName": "deep",
                          "fakesType": ""}):
    # metadata input inherited from SingleMetadataMetricConnections
    diaSources = connectionTypes.Input(
        doc="The catalog of DIA sources.",
        name="{fakesType}{coaddName}Diff_diaSrc",
        storageClass="SourceCatalog",
        dimensions={"instrument", "visit", "detector"},
    )


class FractionAssociatedSourcesConfig(MetadataMetricConfig,
                                      pipelineConnections=FractionAssociatedSourcesConnections):
    pass


# Can't use MetadataMetricTask, because it assumes metadata is the only input
class FractionAssociatedSourcesMetricTask(AbstractMetadataMetricTask):
    """Task that computes the ratio of newly associated sources to DIA sources.
    """
    _DefaultName = "fracAssociatedDiaSources"
    ConfigClass = FractionAssociatedSourcesConfig

    def run(self, metadata, diaSources):
        """Compute the ratio of associated sources to input sources.

        Parameters
        ----------
        metadata : `lsst.daf.base.PropertySet` or `None`
            A metadata object produced by `DiaPipelineTask`.
        diaSources : `lsst.afw.table.SourceCatalog` or `None`
            The source catalog input to `DiaPipelineTask`.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            A `~lsst.pipe.base.Struct` containing the following component:
            - ``measurement``: the value of the metric
              (`lsst.verify.Measurement` or `None`)
        """
        if diaSources is None:
            self.log.info("Nothing to do: no DIA results found.")
            return Struct(measurement=None)
        nSources = len(diaSources)
        if nSources <= 0:
            raise MetricComputationError(
                "No DIA sources found; ratio of associations to DIA sources ill-defined.")

        metadataKeys = self.getInputMetadataKeys(self.config)
        if metadata is not None:
            data = self.extractMetadata(metadata, metadataKeys)
        else:
            data = {dataName: None for dataName in metadataKeys}

        if data["updatedObjects"] is not None and data["newObjects"] is not None:
            try:
                nUpdated = int(data["updatedObjects"])
                nNew = int(data["newObjects"])
            except (ValueError, TypeError) as e:
                raise MetricComputationError("Corrupted value of numUpdatedDiaObjects "
                                             "or numNewDiaObjects") from e
            else:
                fraction = (nUpdated + nNew) / nSources
                meas = Measurement(self.config.metricName, fraction * u.dimensionless_unscaled)
        else:
            self.log.info("Nothing to do: no association results found.")
            meas = None
        return Struct(measurement=meas)

    @classmethod
    def getInputMetadataKeys(cls, config):
        return {"updatedObjects": ".numUpdatedDiaObjects",
                "newObjects": ".numNewDiaObjects"}
