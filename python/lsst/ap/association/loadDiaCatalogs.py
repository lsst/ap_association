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

"""Task for pre-loading DiaSources and DiaObjects within ap_pipe.
"""
import pandas as pd

import lsst.dax.apdb as daxApdb
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as connTypes
import lsst.sphgeom

from lsst.utils.timer import timeMethod

from lsst.ap.association.utils import convertTableToSdmSchema, readSchemaFromApdb, getMidpointFromTimespan

__all__ = ("LoadDiaCatalogsTask", "LoadDiaCatalogsConfig")


class LoadDiaCatalogsConnections(pipeBase.PipelineTaskConnections,
                                 dimensions=("instrument", "group", "detector")):
    regionTime = connTypes.Input(
        doc="The predicted exposure region and time",
        name="regionTimeInfo",
        storageClass="RegionTimeInfo",
        dimensions=("instrument", "group", "detector"),
    )
    diaObjects = connTypes.Output(
        doc="DiaObjects preloaded from the APDB.",
        name="preloaded_diaObjects",
        storageClass="DataFrame",
        dimensions=("instrument", "group", "detector"),
    )
    diaSources = connTypes.Output(
        doc="DiaSources preloaded from the APDB.",
        name="preloaded_diaSources",
        storageClass="DataFrame",
        dimensions=("instrument", "group", "detector"),
    )
    diaForcedSources = connTypes.Output(
        doc="DiaForcedSources preloaded from the APDB.",
        name="preloaded_diaForcedSources",
        storageClass="DataFrame",
        dimensions=("instrument", "group", "detector"),
    )


class LoadDiaCatalogsConfig(pipeBase.PipelineTaskConfig,
                            pipelineConnections=LoadDiaCatalogsConnections):
    """Config class for LoadDiaCatalogsConfig.
    """
    apdb_config_url = pexConfig.Field(
        dtype=str,
        default=None,
        optional=False,
        doc="A config file specifying the APDB and its connection parameters, "
            "typically written by the apdb-cli command-line utility. "
            "The database must already be initialized.",
    )

    pixelMargin = pexConfig.RangeField(
        doc="Padding to add to 4 all edges of the bounding box (pixels)",
        dtype=int,
        default=250,
        min=0,
        deprecated="This config has been replaced by `angleMargin`"
                   "Will be removed after v28.",
    )
    angleMargin = pexConfig.RangeField(
        doc="Padding to add to the radius of the bounding circle (arcseconds)",
        dtype=float,
        default=20,
        min=0,
    )
    doLoadForcedSources = pexConfig.Field(
        dtype=bool,
        default=True,
        deprecated="Added to allow disabling forced sources for performance "
                   "reasons during the ops rehearsal. "
                   "It is expected to be removed.",
        doc="Load forced DiaSource history from the APDB? "
            "This should only be turned off for debugging purposes.",
    )


class LoadDiaCatalogsTask(pipeBase.PipelineTask):
    """Retrieve DiaObjects and associated DiaSources from the Apdb given an
    input exposure.
    """
    ConfigClass = LoadDiaCatalogsConfig
    _DefaultName = "loadDiaCatalogs"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.apdb = daxApdb.Apdb.from_uri(self.config.apdb_config_url)

    @timeMethod
    def run(self, regionTime):
        """Preload all DiaObjects and DiaSources from the Apdb given the
        current exposure.

        Parameters
        ----------
        regionTime : `lsst.pipe.base.utils.RegionTimeInfo`
            A serializable container for a sky region and timespan.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components.

            - ``diaObjects`` : Complete set of DiaObjects covering the input
              exposure padded by ``pixelMargin``. DataFrame is indexed by
              the ``diaObjectId`` column. (`pandas.DataFrame`)
            - ``diaSources`` : Complete set of DiaSources covering the input
              exposure padded by ``pixelMargin``. DataFrame is indexed by
              ``diaObjectId``, ``band``, ``diaSourceId`` columns.
              (`pandas.DataFrame`)
            - ``diaForcedSources`` : Complete set of forced photometered fluxes
            on the past 12 months of difference images at DiaObject locations.

        Raises
        ------
        RuntimeError
            Raised if the Database query failed to load DiaObjects.
        """
        baseRegion = regionTime.region
        region = baseRegion.getBoundingCircle()
        region.dilateBy(lsst.sphgeom.Angle.fromDegrees(self.config.angleMargin/3600.))

        schema = readSchemaFromApdb(self.apdb)

        # This is the first database query.
        diaObjects = self.loadDiaObjects(region, schema)

        # Load diaSources and forced sources up to the time of the exposure
        # The timespan may include significant padding, so use the midpoint to
        #  avoid missing valid recent diaSources.
        visitTime = getMidpointFromTimespan(regionTime.timespan).tai

        diaSources = self.loadDiaSources(diaObjects, region, visitTime, schema)

        if self.config.doLoadForcedSources:
            diaForcedSources = self.loadDiaForcedSources(diaObjects, region, visitTime, schema)
        else:
            diaForcedSources = pd.DataFrame(columns=["diaObjectId", "diaForcedSourceId"])

        return pipeBase.Struct(
            diaObjects=diaObjects,
            diaSources=diaSources,
            diaForcedSources=diaForcedSources)

    @timeMethod
    def loadDiaObjects(self, region, schema):
        """Load DiaObjects from the Apdb based on their HTM location.

        Parameters
        ----------
        region : `sphgeom.Region`
            Region of interest.
        schema : 'dict' of `lsst.dax.apdb.apdbSchema.ApdbSchema`
            A dict of the schemas in the apdb.

        Returns
        -------
        diaObjects : `pandas.DataFrame`
            DiaObjects loaded from the Apdb that are within the area defined
            by ``pixelRanges``.
        """
        self.log.info("Loading DiaObjects")
        diaObjects = self.apdb.getDiaObjects(region)

        diaObjects.set_index("diaObjectId", drop=False, inplace=True)
        if diaObjects.index.has_duplicates:
            self.log.warning(
                "Duplicate DiaObjects loaded from the Apdb. This may cause "
                "downstream pipeline issues. Dropping duplicated rows")
            # Drop duplicates via index and keep the first appearance.
            diaObjects = diaObjects.groupby(diaObjects.index).first()

        return convertTableToSdmSchema(schema, diaObjects, tableName="DiaObject")

    @timeMethod
    def loadDiaSources(self, diaObjects, region, dateTime, schema):
        """Load DiaSources from the Apdb based on their diaObjectId or
        location.

        Variable used to load sources is set in config.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            DiaObjects loaded from the Apdb that are within the area defined
            by ``pixelRanges``.
        region : `sphgeom.Region`
            Region of interest.
        dateTime : `astropy.time.Time`
            Time of the current visit
        schema : 'dict' of `lsst.dax.apdb.apdbSchema.ApdbSchema`
            A dict of the schemas in the apdb.

        Returns
        -------
        DiaSources : `pandas.DataFrame`
            DiaSources loaded from the Apdb that are within the area defined
            by ``pixelRange`` and associated with ``diaObjects``.
        """
        self.log.info("Loading DiaSources")

        diaSources = self.apdb.getDiaSources(region, diaObjects.loc[:, "diaObjectId"], dateTime)

        diaSources.set_index(["diaObjectId", "band", "diaSourceId"],
                             drop=False,
                             inplace=True)
        if diaSources.index.has_duplicates:
            self.log.warning(
                "Duplicate DiaSources loaded from the Apdb. This may cause "
                "downstream pipeline issues. Dropping duplicated rows")
            # Drop duplicates via index and keep the first appearance. Reset
            # due to the index shape being slight different thatn expected.
            diaSources = diaSources.groupby(diaSources.index).first().reset_index(drop=True)
            diaSources.set_index(["diaObjectId", "band", "diaSourceId"],
                                 drop=False,
                                 inplace=True)

        return convertTableToSdmSchema(schema, diaSources, tableName="DiaSource")

    @timeMethod
    def loadDiaForcedSources(self, diaObjects, region, dateTime, schema):
        """Load DiaObjects from the Apdb based on their HTM location.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            DiaObjects loaded from the Apdb.
        region : `sphgeom.Region`
            Region of interest.
        dateTime : `astropy.time.Time`
            Time of the current visit
        schema : 'dict' of `lsst.dax.apdb.apdbSchema.ApdbSchema`
            A dict of the schemas in the apdb.

        Returns
        -------
        diaObjects : `pandas.DataFrame`
            DiaObjects loaded from the Apdb that are within the area defined
            by ``pixelRanges``.
        """
        self.log.info("Loading DiaForcedSources")

        if len(diaObjects) == 0:
            # If no diaObjects are available return an empty DataFrame with
            # the the column used for indexing later in AssociationTask.
            diaForcedSources = pd.DataFrame(columns=["diaObjectId",
                                                     "diaForcedSourceId"])
        else:
            diaForcedSources = self.apdb.getDiaForcedSources(
                region,
                diaObjects.loc[:, "diaObjectId"],
                dateTime)

        diaForcedSources.set_index(["diaObjectId", "diaForcedSourceId"],
                                   drop=False,
                                   inplace=True)
        if diaForcedSources.index.has_duplicates:
            self.log.warning(
                "Duplicate DiaForcedSources loaded from the Apdb. This may "
                "cause downstream pipeline issues. Dropping duplicated rows.")
            # Drop duplicates via index and keep the first appearance. Reset
            # due to the index shape being slightly different than expected.
            diaForcedSources = diaForcedSources.groupby(diaForcedSources.index).first()
            diaForcedSources.reset_index(drop=True, inplace=True)
            diaForcedSources.set_index(["diaObjectId", "diaForcedSourceId"],
                                       drop=False,
                                       inplace=True)

        return convertTableToSdmSchema(schema, diaForcedSources, tableName="DiaForcedSource")
