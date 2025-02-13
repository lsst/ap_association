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
from astropy.units import deg

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
    exposure = connTypes.Input(
        doc="Exposure from which source table was generated",
        name="initial_pvi",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector"),
    )
    sourceTable = connTypes.Input(
        doc="Catalog of calibrated Sources.",
        name="initial_stars_footprints_detector",
        storageClass="SourceCatalog",
        dimensions=("instrument", "visit", "detector"),
    )
    solarSystemObjectTable = connTypes.Input(
        doc="Optional catalog of SolarSolarSystem objects expected to be"
            "observable in this detectorVisit.",
        name="preloaded_SsObjects",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "group", "detector"),
        minimum=0,
    )
    associatedSsSources = connTypes.Output(
        doc="ssSource record columns which can be computed as of association",
        name="ssSingleFrameAssociatedSources",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "visit", "detector"),
    )
    unassociatedSsObjects = connTypes.Output(
        doc="Expected locations of an ssObject with no source",
        name="ssSingleFrameUnassociatedObjects",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "visit", "detector"),
    )


class SsSingleFrameAssociationConfig(pipeBase.PipelineTaskConfig,
                                     pipelineConnections=SsSingleFrameAssociationConnections):
    """Config for DiaPipelineTask.
    """
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
        inputs["band"] = butlerQC.quantum.dataId["band"]
        outputs = self.run(**inputs)

        butlerQC.put(outputs, outputRefs)

    @timeMethod
    def run(self,
            exposure,
            sourceTable,
            band,
            solarSystemObjectTable=None):
        """Process DiaSources and DiaObjects.

        Load previous DiaObjects and their DiaSource history. Calibrate the
        values in the diaSourceCat. Associate new DiaSources with previous
        DiaObjects. Run forced photometry at the updated DiaObject locations.
        Store the results in the Alert Production Database (Apdb).

        Parameters
        ----------
        exposure : `lsst.afw.image.ExposureF`
            Calibrated exposure with wcs and midpoint time.
        sourceTable : `lsst.afw.table.SourceCatalog`
            Newly detected sources.
        band : `str`
            The band in which the new DiaSources were detected.
        solarSystemObjectTable : `astropy.table.Table` or `None`
            Preloaded Solar System objects expected to be visible in the image.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            Results struct with components.
            - ``ssoAssocDiaSources`` : DiaSources that were associated with
              solar system objects in this visit. (`Astropy.table.Table`)
            - ``unAssocDiaSources`` : Set of DiaSources that were not
              associated with any solar system object. (`astropy.table.Table`)
            - ``nTotalSsObjects`` : Total number of SolarSystemObjects
              contained in the CCD footprint. (`int`)
            - ``nAssociatedSsObjects`` : Number of SolarSystemObjects
              that were associated with DiaSources.
            - ``ssSourceData`` : ssSource table data. (`Astropy.table.Table`)


        Raises
        ------
        RuntimeError
            Raised if duplicate DiaObjects or duplicate DiaSources are found.
        """
        if solarSystemObjectTable is None:
            raise pipeBase.NoWorkFound("No ephemerides to associate. Skipping ssSingleFrameAssociation.")

        # Associate DiaSources with DiaObjects
        sourceTable = sourceTable.asAstropy()
        sourceTable['ra'] = sourceTable['coord_ra'].to(deg).value
        sourceTable['dec'] = sourceTable['coord_dec'].to(deg).value
        return self.solarSystemAssociator.run(sourceTable, solarSystemObjectTable,
                                              exposure.visitInfo, exposure.getBBox(), exposure.wcs)
