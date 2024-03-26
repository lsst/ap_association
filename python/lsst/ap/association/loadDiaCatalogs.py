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
from lsst.utils.timer import timeMethod

__all__ = ("LoadDiaCatalogsTask", "LoadDiaCatalogsConfig")


class LoadDiaCatalogsConfig(pexConfig.Config):
    """Config class for LoadDiaCatalogsConfig.
    """
    pixelMargin = pexConfig.RangeField(
        doc="Padding to add to 4 all edges of the bounding box (pixels)",
        dtype=int,
        default=250,
        min=0,
    )


class LoadDiaCatalogsTask(pipeBase.Task):
    """Retrieve DiaObjects and associated DiaSources from the Apdb given an
    input exposure.
    """
    ConfigClass = LoadDiaCatalogsConfig
    _DefaultName = "loadDiaCatalogs"

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)

    @timeMethod
    def run(self, exposure, apdb, doLoadForcedSources=True):
        """Preload all DiaObjects and DiaSources from the Apdb given the
        current exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            An exposure with a bounding box.
        apdb : `lsst.dax.apdb.Apdb`
            AP database connection object.
        doLoadForcedSources : `bool`, optional
            Load forced DiaSource history from the APDB?
            This should only be turned off for debugging purposes.
            Added to allow disabling forced sources for performance
            reasons during the ops rehearsal.

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
        region = self._getRegion(exposure)

        # This is the first database query.
        try:
            diaObjects = self.loadDiaObjects(region, apdb)
        except (OperationalError, ProgrammingError) as e:
            raise RuntimeError(
                "Database query failed to load DiaObjects; did you call "
                "make_apdb.py first? If you did, some other error occurred "
                "during database access of the DiaObject table.") from e

        dateTime = exposure.visitInfo.date

        diaSources = self.loadDiaSources(diaObjects, region, dateTime, apdb)

        if doLoadForcedSources:
            diaForcedSources = self.loadDiaForcedSources(diaObjects, region, dateTime, apdb)
        else:
            diaForcedSources = pd.DataFrame(columns=["diaObjectId", "diaForcedSourceId"])

        return pipeBase.Struct(
            diaObjects=diaObjects,
            diaSources=diaSources,
            diaForcedSources=diaForcedSources)

    @timeMethod
    def loadDiaObjects(self, region, apdb):
        """Load DiaObjects from the Apdb based on their HTM location.

        Parameters
        ----------
        region : `sphgeom.Region`
            Region of interest.
        apdb : `lsst.dax.apdb.Apdb`
            Database connection object to load from.

        Returns
        -------
        diaObjects : `pandas.DataFrame`
            DiaObjects loaded from the Apdb that are within the area defined
            by ``pixelRanges``.
        """
        if region is None:
            # If no area is specified return an empty DataFrame with the
            # the column used for indexing later in AssociationTask.
            diaObjects = pd.DataFrame(columns=["diaObjectId"])
        else:
            diaObjects = apdb.getDiaObjects(region)

        diaObjects.set_index("diaObjectId", drop=False, inplace=True)
        if diaObjects.index.has_duplicates:
            self.log.warning(
                "Duplicate DiaObjects loaded from the Apdb. This may cause "
                "downstream pipeline issues. Dropping duplicated rows")
            # Drop duplicates via index and keep the first appearance.
            diaObjects = diaObjects.groupby(diaObjects.index).first()

        return diaObjects.replace(to_replace=[None], value=np.nan)

    @timeMethod
    def loadDiaSources(self, diaObjects, region, dateTime, apdb):
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

        Returns
        -------
        DiaSources : `pandas.DataFrame`
            DiaSources loaded from the Apdb that are within the area defined
            by ``pixelRange`` and associated with ``diaObjects``.
        """
        if region is None:
            # If no area is specified return an empty DataFrame with the
            # the column used for indexing later in AssociationTask.
            diaSources = pd.DataFrame(columns=["diaObjectId",
                                               "band",
                                               "diaSourceId"])
        else:
            diaSources = apdb.getDiaSources(region, diaObjects.loc[:, "diaObjectId"], dateTime.toAstropy())

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

        return diaSources.replace(to_replace=[None], value=np.nan)

    @timeMethod
    def loadDiaForcedSources(self, diaObjects, region, dateTime, apdb):
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
                region,
                diaObjects.loc[:, "diaObjectId"],
                dateTime.toAstropy())

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

        return diaForcedSources.replace(to_replace=[None], value=np.nan)

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
