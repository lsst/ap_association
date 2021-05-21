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
import numpy as np
import pandas as pd
from sqlalchemy.exc import OperationalError, ProgrammingError

import lsst.geom as geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.sphgeom as sphgeom

__all__ = ("LoadDiaCatalogsTask", "LoadDiaCatalogsConfig")


class LoadDiaCatalogsConfig(pexConfig.Config):
    """Config class for LoadDiaCatalogsConfig.
    """
    htmLevel = pexConfig.RangeField(
        dtype=int,
        doc="Level of the HTM pixelization.",
        default=20,
        min=1,
    )
    htmMaxRanges = pexConfig.RangeField(
        dtype=int,
        doc="Maximum number of HTM (min, max) ranges to return.",
        default=128,
        min=2,
    )
    pixelMargin = pexConfig.RangeField(
        doc="Padding to add to 4 all edges of the bounding box (pixels)",
        dtype=int,
        default=250,
        min=0,
    )
    loadDiaSourcesByPixelId = pexConfig.Field(
        doc="Load DiaSources by their HTM pixelId instead of by their "
            "associated diaObjectId",
        dtype=bool,
        default=False,
    )


class LoadDiaCatalogsTask(pipeBase.Task):
    """Retrieve DiaObjects and associated DiaSources from the Apdb given an
    input exposure.
    """
    ConfigClass = LoadDiaCatalogsConfig
    _DefaultName = "loadDiaCatalogs"

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.pixelator = sphgeom.HtmPixelization(self.config.htmLevel)

    @pipeBase.timeMethod
    def run(self, exposure, apdb):
        """Preload all DiaObjects and DiaSources from the Apdb given the
        current exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            An exposure with a bounding box.
        apdb : `lsst.dax.apdb.Apdb`
            AP database connection object.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components.

            - ``diaObjects`` : Complete set of DiaObjects covering the input
              exposure padded by ``pixelMargin``. DataFrame is indexed by
              the ``diaObjectId`` column. (`pandas.DataFrame`)
            - ``diaSources`` : Complete set of DiaSources covering the input
              exposure padded by ``pixelMargin``. DataFrame is indexed by
              ``diaObjectId``, ``filterName``, ``diaSourceId`` columns.
              (`pandas.DataFrame`)
        """
        visiInfo = exposure.getInfo().getVisitInfo()
        pixelRanges = self._getPixelRanges(exposure)

        # This is the first database query
        try:
            diaObjects = self.loadDiaObjects(pixelRanges, apdb)
        except (OperationalError, ProgrammingError) as e:
            raise RuntimeError(
                "Database query failed to load DiaObjects; did you call "
                "make_apdb.py first? If you did, some other error occurred "
                "during database access of the DiaObject table.") from e

        dateTime = visiInfo.getDate().toPython()

        diaSources = self.loadDiaSources(diaObjects,
                                         dateTime,
                                         pixelRanges,
                                         apdb)

        diaForcedSources = self.loadDiaForcedSources(diaObjects,
                                                     dateTime,
                                                     apdb)

        return pipeBase.Struct(
            diaObjects=diaObjects,
            diaSources=diaSources,
            diaForcedSources=diaForcedSources)

    @pipeBase.timeMethod
    def loadDiaObjects(self, pixelRanges, apdb):
        """Load DiaObjects from the Apdb based on their HTM location.

        Parameters
        ----------
        pixelRanges : `tuple` [`int`]
            Ranges of pixel values that cover region of interest.
        apdb : `lsst.dax.apdb.Apdb`
            Database connection object to load from.

        Returns
        -------
        diaObjects : `pandas.DataFrame`
            DiaObjects loaded from the Apdb that are within the area defined
            by ``pixelRanges``.
        """
        if len(pixelRanges) == 0:
            # If no area is specified return an empty DataFrame with the
            # the column used for indexing later in AssociationTask.
            diaObjects = pd.DataFrame(columns=["diaObjectId"])
        else:
            diaObjects = apdb.getDiaObjects(pixelRanges, return_pandas=True)

        diaObjects.set_index("diaObjectId", drop=False, inplace=True)
        if diaObjects.index.has_duplicates:
            self.log.warn(
                "Duplicate DiaObjects loaded from the Apdb. This may cause "
                "downstream pipeline issues. Dropping duplicated rows")
            # Drop duplicates via index and keep the first appearance.
            diaObjects = diaObjects.groupby(diaObjects.index).first()

        return diaObjects.replace(to_replace=[None], value=np.nan)

    @pipeBase.timeMethod
    def loadDiaSources(self, diaObjects, pixelRanges, dateTime, apdb):
        """Load DiaSources from the Apdb based on their diaObjectId or
        pixelId location.

        Variable used to load sources is set in config.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            DiaObjects loaded from the Apdb that are within the area defined
            by ``pixelRanges``.
        pixelRanges : `list` of `tuples`
            Ranges of pixelIds that cover region of interest.
        dateTime : `datetime.datetime`
            Time of the current visit
        apdb : `lsst.dax.apdb.Apdb`
            Database connection object to load from.

        Returns
        -------
        DiaSources : `pandas.DataFrame`
            DiaSources loaded from the Apdb that are within the area defined
            by ``pixelRange`` and associated with ``diaObjects``.
        """
        if self.config.loadDiaSourcesByPixelId:
            if len(pixelRanges) == 0:
                # If no area is specified return an empty DataFrame with the
                # the column used for indexing later in AssociationTask.
                diaSources = pd.DataFrame(columns=["diaObjectId",
                                                   "filterName",
                                                   "diaSourceId"])
            else:
                diaSources = apdb.getDiaSourcesInRegion(pixelRanges,
                                                        dateTime,
                                                        return_pandas=True)
        else:
            if len(diaObjects) == 0:
                # If no diaObjects are available return an empty DataFrame with
                # the the column used for indexing later in AssociationTask.
                diaSources = pd.DataFrame(columns=["diaObjectId",
                                                   "filterName",
                                                   "diaSourceId"])
            else:
                diaSources = apdb.getDiaSources(
                    diaObjects.loc[:, "diaObjectId"],
                    dateTime,
                    return_pandas=True)

        diaSources.set_index(["diaObjectId", "filterName", "diaSourceId"],
                             drop=False,
                             inplace=True)
        if diaSources.index.has_duplicates:
            self.log.warn(
                "Duplicate DiaSources loaded from the Apdb. This may cause "
                "downstream pipeline issues. Dropping duplicated rows")
            # Drop duplicates via index and keep the first appearance. Reset
            # due to the index shape being slight different thatn expected.
            diaSources = diaSources.groupby(diaSources.index).first().reset_index(drop=True)
            diaSources.set_index(["diaObjectId", "filterName", "diaSourceId"],
                                 drop=False,
                                 inplace=True)

        return diaSources.replace(to_replace=[None], value=np.nan)

    @pipeBase.timeMethod
    def loadDiaForcedSources(self, diaObjects, dateTime, apdb):
        """Load DiaObjects from the Apdb based on their HTM location.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            DiaObjects loaded from the Apdb.
        dateTime : `datetime.datetime`
            Time of the current visit
        apdb : `lsst.dax.apdb.Apdb`
            Database connection object to load from.

        Returns
        -------
        diaObjects : `pandas.DataFrame`
            DiaObjects loaded from the Apdb that are within the area defined
            by ``pixelRanges``.
        """
        if len(diaObjects) == 0:
            # If no diaObjects are available return an empty DataFrame with
            # the the column used for indexing later in AssociationTask.
            diaForcedSources = pd.DataFrame(columns=["diaObjectId",
                                                     "diaForcedSourceId"])
        else:
            diaForcedSources = apdb.getDiaForcedSources(
                diaObjects.loc[:, "diaObjectId"],
                dateTime,
                return_pandas=True)

        diaForcedSources.set_index(["diaObjectId", "diaForcedSourceId"],
                                   drop=False,
                                   inplace=True)
        if diaForcedSources.index.has_duplicates:
            self.log.warn(
                "Duplicate DiaForcedSources loaded from the Apdb. This may "
                "cause downstream pipeline issues. Dropping duplicated rows.")
            # Drop duplicates via index and keep the first appearance. Reset
            # due to the index shape being slight different thatn expected.
            diaForcedSources = diaForcedSources.groupby(diaForcedSources.index).first()
            diaForcedSources.reset_index(drop=True, inplace=True)
            diaForcedSources.set_index(["diaObjectId", "diaForcedSourceId"],
                                       drop=False,
                                       inplace=True)

        return diaForcedSources.replace(to_replace=[None], value=np.nan)

    @pipeBase.timeMethod
    def _getPixelRanges(self, exposure):
        """Calculate covering HTM pixels for the current exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure object with calibrated WCS.

        Returns
        -------
        htmRanges : `list` of `tuples`
            A list of tuples containing `int` values.
        """
        bbox = geom.Box2D(exposure.getBBox())
        bbox.grow(self.config.pixelMargin)
        wcs = exposure.getWcs()

        region = sphgeom.ConvexPolygon([wcs.pixelToSky(pp).getVector()
                                        for pp in bbox.getCorners()])

        indices = self.pixelator.envelope(region, self.config.htmMaxRanges)

        return indices.ranges()
