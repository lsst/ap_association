#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#

"""PipelineTask for associating single-frame sources with solar system objects.

"""

__all__ = ("SsSingleFrameAssociationConfig",
           "SsSingleFrameAssociationTask",
           "SsSingleFrameAssociationConnections")

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as connTypes
from lsst.ap.association.ssoAssociation import SolarSystemAssociationTask
from lsst.utils.timer import timeMethod


class SsSingleFrameAssociationConnections(
        pipeBase.PipelineTaskConnections,
        dimensions=("instrument", "visit", "detector"),
        defaultTemplates={"fakesType": ""}):
    """Butler connections for SsSingleFrameAssociationTask.
    """
    sourceTable = connTypes.Input(
        doc="Catalog of calibrated DiaSources.",
        name="initial_stars_footprints_detector",
        storageClass="DataFrame",
        dimensions=("instrument", "visit", "detector"),
    )
    solarSystemObjectTable = connTypes.Input(
        doc="Catalog of SolarSolarSystem objects expected to be observable in "
            "this detectorVisit.",
        name="preloaded_SsObjects",
        storageClass="DataFrame",
        dimensions=("instrument", "group", "detector"),
        minimum=0,
    )
    exposure = connTypes.Input(
        doc="Calibrated exposure differenced with a template image during "
            "image differencing.",
        name="calexp",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector"),
    )
    associatedSsSources = connTypes.Output(
        doc="ssSource record columns which can be computed as of association",
        name="Diff_associatedSsSources",
        storageClass="DataFrame",
        dimensions=("instrument", "visit", "detector"),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)


class SsSingleFrameAssociationConfig(pipeBase.PipelineTaskConfig,
                                     pipelineConnections=SsSingleFrameAssociationConnections):
    """Config for DiaPipelineTask.
    """
    validBands = pexConfig.ListField(
        dtype=str,
        default=["u", "g", "r", "i", "z", "y"],
        doc="List of bands that are valid for AP processing. To process a "
            "band not on this list, the appropriate band specific columns "
            "must be added to the Apdb schema in dax_apdb.",
    )
    solarSystemAssociator = pexConfig.ConfigurableField(
        target=SolarSystemAssociationTask,
        doc="Task used to associate DiaSources with SolarSystemObjects.",
    )
    imagePixelMargin = pexConfig.RangeField(
        dtype=int,
        default=10,
        min=0,
        doc="Pad the image by this many pixels before removing off-image "
            "diaObjects for association.",
    )


class SsSingleFrameAssociationTask(pipeBase.PipelineTask):
    """Task for loading, associating and storing Difference Image Analysis
    (DIA) Objects and Sources.
    """
    ConfigClass = SsSingleFrameAssociationConfig
    _DefaultName = "ssSingleFrameAssociation"

    def __init__(self, initInputs=None, **kwargs):
        super().__init__(**kwargs)
        self.makeSubtask("solarSystemAssociator")

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        inputs["idGenerator"] = self.config.idGenerator.apply(butlerQC.quantum.dataId)
        inputs["band"] = butlerQC.quantum.dataId["band"]
        inputs["legacySolarSystemTable"] = None

        outputs = self.run(**inputs)

        butlerQC.put(outputs, outputRefs)

    @timeMethod
    def run(self,
            sourceTable,
            exposure,
            band,
            solarSystemObjectTable):
        """Process DiaSources and DiaObjects.

        Load previous DiaObjects and their DiaSource history. Calibrate the
        values in the diaSourceCat. Associate new DiaSources with previous
        DiaObjects. Run forced photometry at the updated DiaObject locations.
        Store the results in the Alert Production Database (Apdb).

        Parameters
        ----------
        exposure : `lsst.afw.image.ExposureF`
            Calibrated exposure
        band : `str`
            The band in which the new DiaSources were detected.
        solarSystemObjectTable : `pandas.DataFrame`
            Preloaded Solar System objects expected to be visible in the image.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            Results struct with components.

            - ``associatedSsSources`` : Catalog of ssSource records.
              (`pandas.DataFrame`)

        Raises
        ------
        RuntimeError
            Raised if duplicate DiaObjects or duplicate DiaSources are found.
        """
        # Associate DiaSources with DiaObjects
        associatedSsSources = self.associateSsSources(sourceTable, solarSystemObjectTable)

        return pipeBase.Struct(associatedSsSources=associatedSsSources)

    @timeMethod
    def associateSources(self, sourceTable, solarSystemObjectTable):
        """Associate single-image sources with ssObjects.

        Parameters
        ----------
        sourceTable : `pandas.DataFrame`
            Newly detected sources.
        solarSystemObjectTable : `pandas.DataFrame`
            Preloaded Solar System objects expected to be visible in the image.

        Returns
        -------
        associatedSsSources : `pandas.DataFrame`
            Table of new ssSources after association.
        """
        ssoAssocResult = self.solarSystemAssociator.run(sourceTable, solarSystemObjectTable)
        associatedSsSources = ssoAssocResult.ssSourceData
        return associatedSsSources
