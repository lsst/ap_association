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

"""Helper functions for tests of DIA catalogs, including generating mock
catalogs for simulated APDB access.
"""
import datetime

import astropy.units
import pandas as pd
import numpy as np

from lsst.afw.cameraGeom.testUtils import DetectorWrapper
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import lsst.daf.butler as dafButler
import lsst.geom
from lsst.pipe.base.utils import RegionTimeInfo
import lsst.sphgeom


def makeDiaObjects(nObjects, exposure, rng):
    """Make a test set of DiaObjects.

    Parameters
    ----------
    nObjects : `int`
        Number of objects to create.
    exposure : `lsst.afw.image.Exposure`
        Exposure to create objects over.

    Returns
    -------
    diaObjects : `pandas.DataFrame`
        DiaObjects generated across the exposure.
    """
    bbox = lsst.geom.Box2D(exposure.getBBox())
    rand_x = rng.uniform(bbox.getMinX(), bbox.getMaxX(), size=nObjects)
    rand_y = rng.uniform(bbox.getMinY(), bbox.getMaxY(), size=nObjects)

    midpointMjdTai = exposure.visitInfo.date.get(system=dafBase.DateTime.MJD)

    data = []
    for idx, (x, y) in enumerate(zip(rand_x, rand_y)):
        coord = exposure.wcs.pixelToSky(x, y)
        newObject = {"ra": coord.getRa().asDegrees(),
                     "dec": coord.getDec().asDegrees(),
                     "radecMjdTai": midpointMjdTai,
                     "diaObjectId": idx + 1,
                     "pmParallaxNdata": 0,
                     "nearbyObj1": 0,
                     "nearbyObj2": 0,
                     "nearbyObj3": 0,
                     "nDiaSources": 5}
        for f in ["u", "g", "r", "i", "z", "y"]:
            newObject["%s_psfFluxNdata" % f] = 0
        data.append(newObject)

    return pd.DataFrame(data=data)


def makeDiaSources(nSources, diaObjectIds, exposure, rng, randomizeObjects=False):
    """Make a test set of DiaSources.

    Parameters
    ----------
    nSources : `int`
        Number of sources to create.
    diaObjectIds : `numpy.ndarray`
        Integer Ids of diaobjects to "associate" with the DiaSources.
    exposure : `lsst.afw.image.Exposure`
        Exposure to create sources over.
    randomizeObjects : `bool`, optional
        If True, randomly draw from `diaObjectIds` to generate the ids in the
        output catalog, otherwise just iterate through them, repeating as
        necessary to get nSources objectIds.

    Returns
    -------
    diaSources : `pandas.DataFrame`
        DiaSources generated across the exposure.
    """
    bbox = lsst.geom.Box2D(exposure.getBBox())
    rand_x = rng.uniform(bbox.getMinX(), bbox.getMaxX(), size=nSources)
    rand_y = rng.uniform(bbox.getMinY(), bbox.getMaxY(), size=nSources)
    if randomizeObjects:
        objectIds = diaObjectIds[rng.randint(len(diaObjectIds), size=nSources)]
    else:
        objectIds = diaObjectIds[[i % len(diaObjectIds) for i in range(nSources)]]

    midpointMjdTai = exposure.visitInfo.date.get(system=dafBase.DateTime.MJD)

    data = []
    for idx, (x, y, objId) in enumerate(zip(rand_x, rand_y, objectIds)):
        coord = exposure.wcs.pixelToSky(x, y)
        # Put together the minimum values for the alert.
        data.append({"ra": coord.getRa().asDegrees(),
                     "dec": coord.getDec().asDegrees(),
                     "x": x,
                     "y": y,
                     "visit": exposure.visitInfo.id,
                     "detector": exposure.detector.getId(),
                     "time_processed": datetime.datetime.now(),
                     "diaObjectId": objId,
                     "ssObjectId": 0,
                     "parentDiaSourceId": 0,
                     "diaSourceId": idx + 1,
                     "midpointMjdTai": midpointMjdTai + 1.0 * idx,
                     "band": exposure.getFilter().bandLabel,
                     "psfNdata": 0,
                     "trailNdata": 0,
                     "dipoleNdata": 0})

    return pd.DataFrame(data=data)


def makeDiaForcedSources(nForcedSources, diaObjectIds, exposure, rng, randomizeObjects=False):
    """Make a test set of DiaSources.

    Parameters
    ----------
    nForcedSources : `int`
        Number of sources to create.
    diaObjectIds : `numpy.ndarray`
        Integer Ids of diaobjects to "associate" with the DiaSources.
    exposure : `lsst.afw.image.Exposure`
        Exposure to create sources over.
    randomizeObjects : `bool`, optional
        If True, randomly draw from `diaObjectIds` to generate the ids in the
        output catalog, otherwise just iterate through them.

    Returns
    -------
    diaForcedSources : `pandas.DataFrame`
        DiaForcedSources generated across the exposure.
    """
    midpointMjdTai = exposure.visitInfo.date.get(system=dafBase.DateTime.MJD)
    visit = exposure.visitInfo.id
    detector = exposure.detector.getId()
    if randomizeObjects:
        objectIds = diaObjectIds[rng.randint(len(diaObjectIds), size=nForcedSources)]
    else:
        objectIds = diaObjectIds[[i % len(diaObjectIds) for i in range(nForcedSources)]]

    data = []
    bbox = exposure.getBBox()

    for i, objId in enumerate(objectIds):
        # Put together the minimum values for the alert.
        x = rng.uniform(bbox.minX, bbox.maxX)
        y = rng.uniform(bbox.minY, bbox.maxY)
        coord = exposure.wcs.pixelToSky(x, y)
        data.append({"diaForcedSourceId": i + 1,
                     "visit": visit + i,
                     "detector": detector,
                     "diaObjectId": objId,
                     "ra": coord.getRa().asDegrees(),
                     "dec": coord.getDec().asDegrees(),
                     "midpointMjdTai": midpointMjdTai + 1.0 * i,
                     "time_processed": datetime.datetime.now(),
                     "band": exposure.getFilter().bandLabel})

    return pd.DataFrame(data=data)


def makeExposure(flipX=False, flipY=False):
    """Create an exposure and flip the x or y (or both) coordinates.

    Returns bounding boxes that are right or left handed around the bounding
    polygon.

    Parameters
    ----------
    flipX : `bool`
        Flip the x coordinate in the WCS.
    flipY : `bool`
        Flip the y coordinate in the WCS.

    Returns
    -------
    exposure : `lsst.afw.image.Exposure`
        Exposure with a valid bounding box and wcs.
    """
    metadata = dafBase.PropertySet()

    metadata.set("SIMPLE", "T")
    metadata.set("BITPIX", -32)
    metadata.set("NAXIS", 2)
    metadata.set("NAXIS1", 1024)
    metadata.set("NAXIS2", 1153)
    metadata.set("RADECSYS", 'FK5')
    metadata.set("EQUINOX", 2000.)

    metadata.setDouble("CRVAL1", 215.604025685476)
    metadata.setDouble("CRVAL2", 53.1595451514076)
    metadata.setDouble("CRPIX1", 1109.99981456774)
    metadata.setDouble("CRPIX2", 560.018167811613)
    metadata.set("CTYPE1", 'RA---SIN')
    metadata.set("CTYPE2", 'DEC--SIN')

    xFlip = 1
    if flipX:
        xFlip = -1
    yFlip = 1
    if flipY:
        yFlip = -1
    metadata.setDouble("CD1_1", xFlip * 5.10808596133527E-05)
    metadata.setDouble("CD1_2", yFlip * 1.85579539217196E-07)
    metadata.setDouble("CD2_2", yFlip * -5.10281493481982E-05)
    metadata.setDouble("CD2_1", xFlip * -8.27440751733828E-07)

    wcs = afwGeom.makeSkyWcs(metadata)
    exposure = afwImage.makeExposure(
        afwImage.makeMaskedImageFromArrays(np.ones((1024, 1153))), wcs)
    detector = DetectorWrapper(id=23, bbox=exposure.getBBox()).detector
    visit = afwImage.VisitInfo(
        exposureTime=200.,
        date=dafBase.DateTime("2014-05-13T17:00:00.000000000",
                              dafBase.DateTime.Timescale.TAI))
    exposure.info.id = 1234
    exposure.setDetector(detector)
    exposure.info.setVisitInfo(visit)
    exposure.setFilter(afwImage.FilterLabel(band='g'))

    return exposure


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


def makeRegionTime(exposure=None, begin=None, end=None):
    """Make a `RegionTimeInfo` for testing

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`, optional
        Exposure to construct a ``RegionTimeInfo`` for.
        If None, a default exposure will be created.
    begin : `astropy.time.Time`, optional
        The start time of the interval.
        If `None`, calculate the start time from the midpoint of the exposure
        and the exposure length.
    end : `astropy.time.Time`, optional
        The end time of the interval.
        If `None`, calculate the end time from the midpoint of the exposure
        and the exposure length.

    Returns
    -------
    regionTime : `lsst.pipe.base.utils.RegionTimeInfo`
        Object containing the spatial region and temporal timespan for an exposure.
    """
    if exposure is None:
        exposure = makeExposure()
    region = getRegion(exposure)
    expTime = exposure.visitInfo.exposureTime*astropy.units.second
    # visitInfo time is the midpoint of the exposure.
    if begin is None:
        begin = exposure.visitInfo.date.toAstropy() - expTime/2
    if end is None:
        end = exposure.visitInfo.date.toAstropy() + expTime/2
    timespan = dafButler.Timespan(begin=begin, end=end)
    return RegionTimeInfo(region=region, timespan=timespan)
