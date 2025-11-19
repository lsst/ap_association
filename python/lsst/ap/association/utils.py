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

"""Utilities for working with the APDB.
"""
__all__ = ("getMidpointFromTimespan", "makeEmptyForcedSourceTable", "getRegion", "paddedRegion",
           "readSchemaFromApdb")


import math
import pandas as pd

from lsst.daf.butler import Timespan
from lsst.dax.apdb import Apdb, ApdbTables, schema_model
import lsst.geom
import lsst.sphgeom

from lsst.pipe.tasks.schemaUtils import convertDataFrameToSdmSchema


def readSchemaFromApdb(apdb: Apdb) -> dict[str, schema_model.Table | None]:
    """Extract the schema from an APDB instance.

    Parameters
    ----------
    apdb : `lsst.dax.apdb.Apdb`
        AP database connection object.

    Returns
    -------
    schemaTable : dict[str, schema_model.Table | None]
        A dict of the schemas in the given table defined in the specified file.
    """
    return {table.table_name(): apdb.tableDef(table) for table in ApdbTables}


def getMidpointFromTimespan(timespan, allowUnbounded=True):
    """Safely retrieve the midpoint in TAI from a Timespan.

    Parameters
    ----------
    timespan : `lsst.daf.butler.Timespan` of `astropy.time.Time`
        A Timespan centered on the midpoint of a visit.
    allowUnbounded : `bool`, optional
        If set, return the start or end of an unbounded timespan.

    Returns
    -------
    midpoint : `astropy.time.Time`
        The midpoint of the timespan.

    Raises
    ------
    ValueError
        Raised if either the start or end of the timespan is None, and
        ``allowUnbounded`` is not set.
    ValueError
        Raised if the timespan is empty.
    """
    if (timespan.begin == Timespan.EMPTY) or (timespan.begin == Timespan.EMPTY):
        raise ValueError("Cannot compute midpoint: EMPTY Timespan.")

    try:
        interval = timespan.end - timespan.begin
        return (timespan.begin + interval/2).tai
    except TypeError as e:
        if allowUnbounded:
            if timespan.end is not None:
                return timespan.end.tai
            elif timespan.begin is not None:
                return timespan.begin.tai
            else:
                raise ValueError("Cannot compute midpoint: unbounded timespan.") from e
        else:
            raise ValueError("Cannot compute midpoint: unbounded timespan.") from e


def makeEmptyForcedSourceTable(schema):
    """Return a dataframe with the correct columns for diaForcedSources table.

    Returns
    -------
    diaForcedSources : `pandas.DataFrame`
        Empty dataframe.
    """
    diaForcedSources = convertDataFrameToSdmSchema(schema, pd.DataFrame(), tableName="DiaForcedSource")
    return diaForcedSources


def getRegion(exposure):
    """Calculate an enveloping region for an exposure.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure object with calibrated WCS.

    Returns
    -------
    region : `lsst.sphgeom.Region`
        Region enveloping an exposure.
    """
    # Bounding box needs to be a `Box2D` not a `Box2I` for `wcs.pixelToSky()`
    bbox = lsst.geom.Box2D(exposure.getBBox())
    wcs = exposure.getWcs()

    region = lsst.sphgeom.ConvexPolygon([pp.getVector() for pp in wcs.pixelToSky(bbox.getCorners())])

    return region


def paddedRegion(region, margin):
    """Return a region that has been expanded by a buffer.

    Parameters
    ----------
    region : `lsst.sphgeom.Region`
        The region to pad.
    margin : `lsst.sphgeom.Angle`
        The amount by which to increase the region.

    Returns
    -------
    padded : `lsst.sphgeom.Region`
        An enlarged copy of ``region``.
    """
    # region is almost certainly a (padded) detector bounding box.
    if isinstance(region, lsst.sphgeom.ConvexPolygon):
        # This is an ad-hoc, approximate implementation. It should be good
        # enough for catalog loading, but is not a general-purpose solution.
        center = lsst.geom.SpherePoint(region.getCentroid())
        corners = [lsst.geom.SpherePoint(c) for c in region.getVertices()]
        # Approximate the region as a Euclidian square
        # geom.Angle(sphgeom.Angle) converter not pybind-wrapped???
        diagonal_margin = lsst.geom.Angle(margin.asRadians() * math.sqrt(2.0))
        padded = [c.offset(center.bearingTo(c), diagonal_margin) for c in corners]
        return lsst.sphgeom.ConvexPolygon.convexHull([c.getVector() for c in padded])
    elif hasattr(region, "dilatedBy"):
        return region.dilatedBy(margin)
    else:
        return region.getBoundingCircle().dilatedBy(margin)
