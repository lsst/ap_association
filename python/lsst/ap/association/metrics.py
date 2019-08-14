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
           "TotalUnassociatedDiaObjectsMetricTask",
           ]


import astropy.units as u

from lsst.verify import Measurement
from lsst.verify.gen2tasks import register
from lsst.verify.tasks import MetadataMetricTask, PpdbMetricTask, MetricComputationError


@register("numNewDiaObjects")
class NumberNewDiaObjectsMetricTask(MetadataMetricTask):
    """Task that computes the number of DIASources that create new DIAObjects
    in an image, visit, etc..
    """
    _DefaultName = "numNewDiaObjects"

    def makeMeasurement(self, values):
        """Compute the number of new DIAObjects.

        Parameters
        ----------
        values : sequence [`dict` [`str`, `int` or `None`]]
            A list where each element corresponds to a metadata object passed
            to `run`. Each `dict` has the following key:

            ``"newObjects"``
                The number of new objects created for this image (`int`
                or `None`). May be `None` if the image was not successfully
                associated.

        Returns
        -------
        measurement : `lsst.verify.Measurement` or `None`
            The total number of new objects.
        """
        nNew = 0
        associated = False
        for value in values:
            if value["newObjects"] is not None:
                try:
                    nNew += value["newObjects"]
                except TypeError as e:
                    raise MetricComputationError("Corrupted value of numNewDiaObjects") from e
                associated = True

        if associated:
            return Measurement(self.getOutputMetricName(self.config),
                               nNew * u.count)
        else:
            self.log.info("Nothing to do: no association results found.")
            return None

    @classmethod
    def getInputMetadataKeys(cls, config):
        return {"newObjects": ".numNewDiaObjects"}

    @classmethod
    def getOutputMetricName(cls, config):
        return "ap_association.numNewDiaObjects"


@register("numUnassociatedDiaObjects")
class NumberUnassociatedDiaObjectsMetricTask(MetadataMetricTask):
    """Task that computes the number of previously-known DIAObjects that do
    not have detected DIASources in an image, visit, etc..
    """
    _DefaultName = "numUnassociatedDiaObjects"

    def makeMeasurement(self, values):
        """Compute the number of non-updated DIAObjects.

        Parameters
        ----------
        values : sequence [`dict` [`str`, `int` or `None`]]
            A list where each element corresponds to a metadata object passed
            to `run`. Each `dict` has the following key:

            ``"unassociatedObjects"``
                The number of DIAObjects not associated with a DiaSource in
                this image (`int` or `None`). May be `None` if the image was
                not successfully associated.

        Returns
        -------
        measurement : `lsst.verify.Measurement` or `None`
            The total number of unassociated objects.
        """
        nNew = 0
        associated = False
        for value in values:
            if value["unassociatedObjects"] is not None:
                try:
                    nNew += value["unassociatedObjects"]
                except TypeError as e:
                    raise MetricComputationError("Corrupted value of numUnassociatedDiaObjects") from e
                associated = True

        if associated:
            return Measurement(self.getOutputMetricName(self.config),
                               nNew * u.count)
        else:
            self.log.info("Nothing to do: no association results found.")
            return None

    @classmethod
    def getInputMetadataKeys(cls, config):
        return {"unassociatedObjects": ".numUnassociatedDiaObjects"}

    @classmethod
    def getOutputMetricName(cls, config):
        return "ap_association.numUnassociatedDiaObjects"


@register("fracUpdatedDiaObjects")
class FractionUpdatedDiaObjectsMetricTask(MetadataMetricTask):
    """Task that computes the fraction of previously created DIAObjects that
    have a new association in this image, visit, etc..
    """
    _DefaultName = "fracUpdatedDiaObjects"

    def makeMeasurement(self, values):
        """Compute the number of non-updated DIAObjects.

        Parameters
        ----------
        values : sequence [`dict` [`str`, `int` or `None`]]
            A list where each element corresponds to a metadata object passed
            to `run`. Each `dict` has the following keys:

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
        nUpdated = 0
        nUnassociated = 0
        associated = False
        for value in values:
            if value["updatedObjects"] is not None \
                    and value["unassociatedObjects"] is not None:
                try:
                    nUpdated += value["updatedObjects"]
                    nUnassociated += value["unassociatedObjects"]
                except TypeError as e:
                    raise MetricComputationError("Corrupted value of numUpdatedDiaObjects "
                                                 "or numUnassociatedDiaObjects") from e
                associated = True

        if associated:
            if nUpdated <= 0 and nUnassociated <= 0:
                raise MetricComputationError("No pre-existing DIAObjects; can't compute updated fraction.")
            else:
                fraction = nUpdated / (nUpdated + nUnassociated)
                return Measurement(self.getOutputMetricName(self.config),
                                   fraction * u.dimensionless_unscaled)
        else:
            self.log.info("Nothing to do: no association results found.")
            return None

    @classmethod
    def getInputMetadataKeys(cls, config):
        return {"updatedObjects": ".numUpdatedDiaObjects",
                "unassociatedObjects": ".numUnassociatedDiaObjects"}

    @classmethod
    def getOutputMetricName(cls, config):
        return "ap_association.fracUpdatedDiaObjects"


@register("totalUnassociatedDiaObjects")
class TotalUnassociatedDiaObjectsMetricTask(PpdbMetricTask):
    """Task that computes the number of DIAObjects with only one
    associated DIASource.
    """
    _DefaultName = "totalUnassociatedDiaObjects"

    def makeMeasurement(self, dbHandle, outputDataId):
        """Compute the number of unassociated DIAObjects.

        Parameters
        ----------
        dbHandle : `lsst.dax.ppdb.Ppdb`
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
        if outputDataId.keys():
            raise ValueError("%s must not be associated with specific data IDs."
                             % self.getOutputMetricName(self.config))

        try:
            nUnassociatedDiaObjects = dbHandle.countUnassociatedObjects()
        except Exception as e:
            raise MetricComputationError("Could not get unassociated objects from database") from e

        meas = Measurement(self.getOutputMetricName(self.config), nUnassociatedDiaObjects * u.count)
        return meas

    @classmethod
    def getOutputMetricName(cls, config):
        return "ap_association.totalUnassociatedDiaObjects"
