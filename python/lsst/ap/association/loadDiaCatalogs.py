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
import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as connTypes
import lsst.sphgeom

from lsst.utils.timer import timeMethod, duration_from_timeMethod

from lsst.ap.association.utils import getMidpointFromTimespan, paddedRegion, readSchemaFromApdb
from lsst.pipe.tasks.schemaUtils import convertDataFrameToSdmSchema

__all__ = ("LoadDiaCatalogsTask", "LoadDiaCatalogsConfig", "dropDuplicateRows",
           "loadDiaObjectsFromApdb", "loadDiaSourcesFromApdb", "loadDiaForcedSourcesFromApdb")


def dropDuplicateRows(catalog, index, name, log):
    """Index a catalog loaded from the Apdb and drop any duplicate rows.

    Where several rows share an index, the first is kept exactly as the Apdb
    returned it and the others are discarded. The rows stay in the order the
    Apdb returned them.

    Parameters
    ----------
    catalog : `pandas.DataFrame`
        Catalog loaded from the Apdb. Left unchanged.
    index : `str` or `list` [`str`]
        Column or columns to index the catalog on.
    name : `str`
        Name of the catalog, for logging.
    log : `logging.Logger`
        Log to report duplicates to.

    Returns
    -------
    catalog : `pandas.DataFrame`
        A new catalog, indexed by ``index`` and free of duplicates.
    """
    catalog = catalog.set_index(index, drop=False)
    if catalog.index.has_duplicates:
        log.warning("Duplicate %s loaded from the Apdb. This may cause "
                    "downstream pipeline issues. Dropping duplicated rows.", name)
        catalog = catalog[~catalog.index.duplicated(keep="first")]
    return catalog


def loadDiaObjectsFromApdb(apdb, region, schema, log):
    """Load DiaObjects from the Apdb based on their HTM location.

    Parameters
    ----------
    apdb : `lsst.dax.apdb.Apdb`
        Database to load the DiaObjects from.
    region : `sphgeom.Region`
        Region of interest, including any padding.
    schema : `dict` of `lsst.dax.apdb.apdbSchema.ApdbSchema`
        A dict of the schemas in the apdb.
    log : `logging.Logger`
        Log to report the loaded catalog to.

    Returns
    -------
    diaObjects : `pandas.DataFrame`
        DiaObjects within ``region``, indexed by ``diaObjectId``.
    """
    diaObjects = apdb.getDiaObjects(region)
    diaObjects = dropDuplicateRows(diaObjects, "diaObjectId", "DiaObjects", log)
    log.info("Loaded %i DiaObjects", len(diaObjects))
    return convertDataFrameToSdmSchema(schema, diaObjects, tableName="DiaObject", skipIndex=True)


def loadDiaSourcesFromApdb(apdb, region, diaObjectIds, dateTime, schema, log):
    """Load DiaSources from the Apdb based on their diaObjectId or location.

    Parameters
    ----------
    apdb : `lsst.dax.apdb.Apdb`
        Database to load the DiaSources from.
    region : `sphgeom.Region`
        Region of interest, including any padding.
    diaObjectIds : `pandas.Series`
        Ids of the DiaObjects to load the history for.
    dateTime : `astropy.time.Time`
        Time of the current visit.
    schema : `dict` of `lsst.dax.apdb.apdbSchema.ApdbSchema`
        A dict of the schemas in the apdb.
    log : `logging.Logger`
        Log to report the loaded catalog to.

    Returns
    -------
    diaSources : `pandas.DataFrame`
        DiaSource history, indexed by ``diaObjectId``, ``band``, and
        ``diaSourceId``.
    """
    diaSources = apdb.getDiaSources(region, diaObjectIds, dateTime)
    diaSources = dropDuplicateRows(diaSources, ["diaObjectId", "band", "diaSourceId"], "DiaSources", log)
    log.info("Loaded %i DiaSources", len(diaSources))
    return convertDataFrameToSdmSchema(schema, diaSources, tableName="DiaSource", skipIndex=True)


def loadDiaForcedSourcesFromApdb(apdb, region, diaObjectIds, dateTime, schema, log):
    """Load DiaForcedSources from the Apdb based on their diaObjectId.

    Parameters
    ----------
    apdb : `lsst.dax.apdb.Apdb`
        Database to load the DiaForcedSources from.
    region : `sphgeom.Region`
        Region of interest, including any padding.
    diaObjectIds : `pandas.Series`
        Ids of the DiaObjects to load the history for.
    dateTime : `astropy.time.Time`
        Time of the current visit.
    schema : `dict` of `lsst.dax.apdb.apdbSchema.ApdbSchema`
        A dict of the schemas in the apdb.
    log : `logging.Logger`
        Log to report the loaded catalog to.

    Returns
    -------
    diaForcedSources : `pandas.DataFrame`
        DiaForcedSource history, indexed by ``diaObjectId`` and
        ``diaForcedSourceId``.
    """
    if len(diaObjectIds) == 0:
        # If no diaObjects are available return an empty DataFrame with
        # the minimal set of columns.
        diaForcedSources = pd.DataFrame(columns=["diaObjectId", "diaForcedSourceId"])
    else:
        diaForcedSources = apdb.getDiaForcedSources(region, diaObjectIds, dateTime)
    diaForcedSources = dropDuplicateRows(diaForcedSources, ["diaObjectId", "diaForcedSourceId"],
                                         "DiaForcedSources", log)
    nVisits = 0 if diaForcedSources.empty else len(set(diaForcedSources["visit"]))
    log.info("Loaded %i DiaForcedSources from %i visits", len(diaForcedSources), nVisits)
    return convertDataFrameToSdmSchema(schema, diaForcedSources, tableName="DiaForcedSource",
                                       skipIndex=True)


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
        storageClass="ArrowAstropy",
        dimensions=("instrument", "group", "detector"),
    )
    diaSources = connTypes.Output(
        doc="DiaSources preloaded from the APDB.",
        name="preloaded_diaSources",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "group", "detector"),
    )
    diaForcedSources = connTypes.Output(
        doc="DiaForcedSources preloaded from the APDB.",
        name="preloaded_diaForcedSources",
        storageClass="ArrowAstropy",
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
              exposure padded by ``angleMargin``. DataFrame is indexed by
              the ``diaObjectId`` column. (`pandas.DataFrame`)
            - ``diaSources`` : Complete set of DiaSources covering the input
              exposure padded by ``angleMargin``. DataFrame is indexed by
              ``diaObjectId``, ``band``, ``diaSourceId`` columns.
              (`pandas.DataFrame`)
            - ``diaForcedSources`` : Complete set of forced photometered
              fluxes on the past 12 months of difference images at DiaObject
              locations, indexed by ``diaObjectId`` and
              ``diaForcedSourceId``. (`pandas.DataFrame`)

        Raises
        ------
        RuntimeError
            Raised if the Database query failed to load DiaObjects.
        """
        region = paddedRegion(regionTime.region,
                              lsst.sphgeom.Angle.fromDegrees(self.config.angleMargin/3600.))
        schema = readSchemaFromApdb(self.apdb)

        try:
            # This is the first database query.
            try:
                diaObjects = self.loadDiaObjects(region, schema)
            finally:
                self.metadata["loadDiaObjectsDuration"] = duration_from_timeMethod(
                    self.metadata, "loadDiaObjects", clock="Utc")
                self.log.verbose("DiaObjects: Took %.4f seconds", self.metadata["loadDiaObjectsDuration"])

            # Load diaSources and forced sources up to the time of the exposure
            # The timespan may include significant padding, so use the midpoint to
            #  avoid missing valid recent diaSources.
            visitTime = getMidpointFromTimespan(regionTime.timespan)

            try:
                diaSources = self.loadDiaSources(diaObjects, region, visitTime, schema)
            finally:
                self.metadata["loadDiaSourcesDuration"] = duration_from_timeMethod(
                    self.metadata, "loadDiaSources", clock="Utc")
                self.log.verbose("DiaSources: Took %.4f seconds", self.metadata["loadDiaSourcesDuration"])

            if self.config.doLoadForcedSources:
                try:
                    diaForcedSources = self.loadDiaForcedSources(diaObjects, region, visitTime, schema)
                finally:
                    self.metadata["loadDiaForcedSourcesDuration"] = duration_from_timeMethod(
                        self.metadata, "loadDiaForcedSources", clock="Utc")
                    self.log.verbose("DiaForcedSources: Took %.4f seconds",
                                     self.metadata["loadDiaForcedSourcesDuration"])
            else:
                diaForcedSources = pd.DataFrame(columns=["diaObjectId", "diaForcedSourceId"])
                self.metadata["loadDiaForcedSourcesDuration"] = -1
        finally:
            # Loki can add up the three individual times, but a combined log puts less load on the server.
            self.log.verbose("All catalogs: Took %.4f seconds",
                             self.metadata.get("loadDiaObjectsDuration", 0)
                             + self.metadata.get("loadDiaSourcesDuration", 0)
                             + max(0, self.metadata.get("loadDiaForcedSourcesDuration", 0))
                             )

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
            DiaObjects loaded from the Apdb that are within ``region``,
            indexed by ``diaObjectId``.
        """
        return loadDiaObjectsFromApdb(self.apdb, region, schema, self.log)

    @timeMethod
    def loadDiaSources(self, diaObjects, region, dateTime, schema):
        """Load DiaSources from the Apdb based on their diaObjectId or
        location.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            DiaObjects to load the history for, indexed by ``diaObjectId``.
        region : `sphgeom.Region`
            Region of interest.
        dateTime : `astropy.time.Time`
            Time of the current visit
        schema : 'dict' of `lsst.dax.apdb.apdbSchema.ApdbSchema`
            A dict of the schemas in the apdb.

        Returns
        -------
        diaSources : `pandas.DataFrame`
            DiaSources loaded from the Apdb that are within ``region`` and
            associated with ``diaObjects``, indexed by ``diaObjectId``,
            ``band``, and ``diaSourceId``.
        """
        return loadDiaSourcesFromApdb(self.apdb, region, diaObjects.loc[:, "diaObjectId"], dateTime,
                                      schema, self.log)

    @timeMethod
    def loadDiaForcedSources(self, diaObjects, region, dateTime, schema):
        """Load DiaForcedSources from the Apdb based on their diaObjectId.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            DiaObjects to load the history for, indexed by ``diaObjectId``.
        region : `sphgeom.Region`
            Region of interest.
        dateTime : `astropy.time.Time`
            Time of the current visit
        schema : 'dict' of `lsst.dax.apdb.apdbSchema.ApdbSchema`
            A dict of the schemas in the apdb.

        Returns
        -------
        diaForcedSources : `pandas.DataFrame`
            DiaForcedSources loaded from the Apdb that are associated with
            ``diaObjects``, indexed by ``diaObjectId`` and
            ``diaForcedSourceId``.
        """
        return loadDiaForcedSourcesFromApdb(self.apdb, region, diaObjects.loc[:, "diaObjectId"],
                                            dateTime, schema, self.log)
