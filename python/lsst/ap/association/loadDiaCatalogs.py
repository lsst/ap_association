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
from sqlalchemy.exc import OperationalError, ProgrammingError

import lsst.dax.apdb as daxApdb
import lsst.geom as geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as connTypes
import lsst.sphgeom as sphgeom
from lsst.utils.timer import timeMethod

from lsst.ap.association.utils import convertTableToSdmSchema, readSchemaFromApdb

__all__ = ("LoadDiaCatalogsTask", "LoadDiaCatalogsConfig")


class LoadDiaCatalogsConnections(pipeBase.PipelineTaskConnections,
                                 dimensions=("instrument", "group", "detector")):
    regionTime = connTypes.Input(
        doc="The predicted exposure region and time",
        name="regionTimeInfo",
        storageClass="RegionTimeInfo",
        dimensions={"instrument", "group", "detector"},
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
        region = regionTime.region
        schema = readSchemaFromApdb(self.apdb)

        # This is the first database query.
        try:
            diaObjects = self.loadDiaObjects(region, self.apdb, schema)
        except (OperationalError, ProgrammingError) as e:
            raise RuntimeError(
                "Database query failed to load DiaObjects; did you call "
                "apdb_cli.py first? If you did, some other error occurred "
                "during database access of the DiaObject table.") from e

        # Load diaSources and forced sources up to the time of the beginning of the exposure
        visitTime = regionTime.timespan.begin.tai

        diaSources = self.loadDiaSources(diaObjects, region, visitTime, self.apdb, schema)

        diaForcedSources = self.loadDiaForcedSources(diaObjects, region, visitTime, self.apdb, schema)

        return pipeBase.Struct(
            diaObjects=diaObjects,
            diaSources=diaSources,
            diaForcedSources=diaForcedSources)

    @timeMethod
    def loadDiaObjects(self, region, apdb, schema):
        """Load DiaObjects from the Apdb based on their HTM location.

        Parameters
        ----------
        region : `sphgeom.Region`
            Region of interest.
        apdb : `lsst.dax.apdb.Apdb`
            Database connection object to load from.
        schema : 'dict' of `lsst.dax.apdb.apdbSchema.ApdbSchema`
            A dict of the schemas in the apdb.

        Returns
        -------
        diaObjects : `pandas.DataFrame`
            DiaObjects loaded from the Apdb that are within the area defined
            by ``pixelRanges``.
        """
        self.log.info("Loading DiaObjects")
        diaObjects = apdb.getDiaObjects(region)

        diaObjects.set_index("diaObjectId", drop=False, inplace=True)
        if diaObjects.index.has_duplicates:
            self.log.warning(
                "Duplicate DiaObjects loaded from the Apdb. This may cause "
                "downstream pipeline issues. Dropping duplicated rows")
            # Drop duplicates via index and keep the first appearance.
            diaObjects = diaObjects.groupby(diaObjects.index).first()

        return convertTableToSdmSchema(schema, diaObjects, tableName="DiaObject")

    @timeMethod
    def loadDiaSources(self, diaObjects, region, dateTime, apdb, schema):
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
        dateTime : `lsst.daf.base.DateTime`
            Time of the current visit
        apdb : `lsst.dax.apdb.Apdb`
            Database connection object to load from.
        schema : 'dict' of `lsst.dax.apdb.apdbSchema.ApdbSchema`
            A dict of the schemas in the apdb.

        Returns
        -------
        DiaSources : `pandas.DataFrame`
            DiaSources loaded from the Apdb that are within the area defined
            by ``pixelRange`` and associated with ``diaObjects``.
        """
        self.log.info("Loading DiaSources")

        diaSources = apdb.getDiaSources(region, diaObjects.loc[:, "diaObjectId"], dateTime)

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
    def loadDiaForcedSources(self, diaObjects, region, dateTime, apdb, schema):
        """Load DiaObjects from the Apdb based on their HTM location.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            DiaObjects loaded from the Apdb.
        region : `sphgeom.Region`
            Region of interest.
        dateTime : `lsst.daf.base.DateTime`
            Time of the current visit
        apdb : `lsst.dax.apdb.Apdb`
            Database connection object to load from.
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
            diaForcedSources = apdb.getDiaForcedSources(
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

    @timeMethod
    def _getRegion(self, exposure):
        """Calculate an enveloping region for an exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure object with calibrated WCS.

        Returns
        -------
        region : `sphgeom.Region`
            Region enveloping an exposure.
        """
        bbox = geom.Box2D(exposure.getBBox())
        bbox.grow(self.config.pixelMargin)
        wcs = exposure.getWcs()

        region = sphgeom.ConvexPolygon([wcs.pixelToSky(pp).getVector()
                                        for pp in bbox.getCorners()])

        return region
