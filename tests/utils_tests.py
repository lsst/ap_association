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
import pandas as pd
import numpy as np

import lsst.daf.base as dafBase
import lsst.geom


def makeDiaObjects(nObjects, exposure):
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
    rand_x = np.random.uniform(bbox.getMinX(), bbox.getMaxX(), size=nObjects)
    rand_y = np.random.uniform(bbox.getMinY(), bbox.getMaxY(), size=nObjects)

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
                     "flags": 1,
                     "nDiaSources": 5}
        for f in ["u", "g", "r", "i", "z", "y"]:
            newObject["%s_psfFluxNdata" % f] = 0
        data.append(newObject)

    return pd.DataFrame(data=data)


def makeDiaSources(nSources, diaObjectIds, exposure, randomizeObjects=False):
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
    rand_x = np.random.uniform(bbox.getMinX(), bbox.getMaxX(), size=nSources)
    rand_y = np.random.uniform(bbox.getMinY(), bbox.getMaxY(), size=nSources)
    if randomizeObjects:
        objectIds = diaObjectIds[np.random.randint(len(diaObjectIds), size=nSources)]
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
                     "ccdVisitId": exposure.info.id,
                     "time_processed": datetime.datetime.now(),
                     "diaObjectId": objId,
                     "ssObjectId": 0,
                     "parentDiaSourceId": 0,
                     "prv_procOrder": 0,
                     "diaSourceId": idx + 1,
                     "midpointMjdTai": midpointMjdTai + 1.0 * idx,
                     "band": exposure.getFilter().bandLabel,
                     "psfNdata": 0,
                     "trailNdata": 0,
                     "dipoleNdata": 0,
                     "flags": 1})

    return pd.DataFrame(data=data)


def makeDiaForcedSources(nForcedSources, diaObjectIds, exposure, randomizeObjects=False):
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
    ccdVisitId = exposure.info.id
    if randomizeObjects:
        objectIds = diaObjectIds[np.random.randint(len(diaObjectIds), size=nForcedSources)]
    else:
        objectIds = diaObjectIds[[i % len(diaObjectIds) for i in range(nForcedSources)]]

    data = []
    for i, objId in enumerate(objectIds):
        # Put together the minimum values for the alert.
        data.append({"diaForcedSourceId": i + 1,
                     "ccdVisitId": ccdVisitId + i,
                     "diaObjectId": objId,
                     "midpointMjdTai": midpointMjdTai + 1.0 * i,
                     "time_processed": datetime.datetime.now(),
                     "band": exposure.getFilter().bandLabel,
                     "flags": 0})

    return pd.DataFrame(data=data)
