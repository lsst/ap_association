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

"""Solar System Object Query to MPSky in place of a internal Rubin solar
system object caching/retrieval code.

Will compute the location for of known SSObjects within a visit-detector combination.
"""

__all__ = ["MPSkyEphemerisQueryConfig", "MPSkyEphemerisQueryTask"]


import os
import pandas as pd
import pyarrow as pa
import requests

from lsst.ap.association.utils import getMidpointFromTimespan
from lsst.geom import SpherePoint
import lsst.pex.config as pexConfig
from lsst.utils.timer import timeMethod

from lsst.pipe.base import connectionTypes, NoWorkFound, PipelineTask, \
    PipelineTaskConfig, PipelineTaskConnections, Struct


class MPSkyEphemerisQueryConnections(PipelineTaskConnections,
                                     dimensions=("instrument",
                                                 "group", "detector")):

    predictedRegionTime = connectionTypes.Input(
        doc="The predicted exposure region and time",
        name="regionTimeInfo",
        storageClass="RegionTimeInfo",
        dimensions={"instrument", "group", "detector"},
    )

    ssObjects = connectionTypes.Output(
        doc="MPSky-provided Solar System objects observable in this detector-visit",
        name="preloaded_SsObjects",
        storageClass="DataFrame",
        dimensions=("instrument", "group", "detector"),
    )


class MPSkyEphemerisQueryConfig(
        PipelineTaskConfig,
        pipelineConnections=MPSkyEphemerisQueryConnections):
    observerCode = pexConfig.Field(
        dtype=str,
        doc="IAU Minor Planet Center observer code for queries "
            "(default is X05 for Rubin Obs./LSST)",
        default='X05'
    )
    queryBufferRadiusDegrees = pexConfig.Field(
        dtype=float,
        doc="Buffer radius in degrees added to detector bounding circle for ephemeris "
            "cone search. Defaults to 10 deg/day * 30 minutes",
        default=0.208
    )
    mpSkyRequestTimeoutSeconds = pexConfig.Field(
        dtype=float,
        doc="Time in seconds to wait for mpSky request before failing ",
        default=1.0
    )


class MPSkyEphemerisQueryTask(PipelineTask):
    """Task to query the MPSky service and retrieve the solar system objects
    that are observable within the input visit.
    """
    ConfigClass = MPSkyEphemerisQueryConfig
    _DefaultName = "mpSkyEphemerisQuery"

    @timeMethod
    def run(self, predictedRegionTime):
        """Parse the information on the current visit and retrieve the
        observable solar system objects from MPSky.

        Parameters
        ----------
        predictedRegionTime : `lsst.pipe.base.utils.RegionTimeInfo`
            Predicted footprint and timespan of the exposure

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components:

            - ``ssObjects`` : `pandas.DataFrame`
                DataFrame containing Solar System Objects near the detector
                footprint as retrieved by MPSky. The columns are as follows:

                ``Name``
                    object name (`str`)
                ``ra``
                    RA in decimal degrees (`float`)
                ``dec``
                    DEC in decimal degrees (`float`)
                ``obj_poly``
                    DO NOT USE until DM-46069 is resolved
                ``obs_poly``
                    DO NOT USE until DM-46069 is resolved
        """
        # Get detector center and timespan midpoint from predictedRegionTime.
        region = predictedRegionTime.region
        timespan = predictedRegionTime.timespan
        expCenter = SpherePoint(region.getBoundingCircle().getCenter())
        expRadius = region.getBoundingCircle().getOpeningAngle().asDegrees()
        expMidPointEPOCH = getMidpointFromTimespan(timespan, allowUnbounded=False).mjd

        # MPSky service query
        mpSkyURL = os.environ.get('MP_SKY_URL', '')
        mpSkySsObjects = self._mpSkyConeSearch(expCenter, expMidPointEPOCH,
                                               expRadius + self.config.queryBufferRadiusDegrees, mpSkyURL)

        return Struct(
            ssObjects=mpSkySsObjects,
        )

    def read_mp_sky_response(self, response):
        """Extract ephemerides from an MPSky web request

        Parameters
        ----------
        response : `requests.Response`
            MPSky message

        Returns
        -------
        objID : `np.ndarray`
            Designations of nearby objects
        ra : `np.ndarray`
            Array of object right ascensions
        dec : `np.ndarray`
            Array of object declinations
        object_polynomial : `np.ndarray`, (N,M)
            Array of object cartesian position polynomials
        observer_polynomial : `np.ndarray`, (N,M)
            Array of observer cartesian position polynomials
        """
        with pa.input_stream(memoryview(response.content)) as stream:
            stream.seek(0)
            object_polynomial = pa.ipc.read_tensor(stream)
            observer_polynomial = pa.ipc.read_tensor(stream)
            with pa.ipc.open_stream(stream) as reader:
                columns = next(reader)
        objID = columns["name"].to_numpy(zero_copy_only=False)
        ra = columns["ra"].to_numpy()
        dec = columns["dec"].to_numpy()
        return objID, ra, dec, object_polynomial, observer_polynomial

    def _mpSkyConeSearch(self, expCenter, epochMJD, queryRadius, mpSkyURL):
        """Query MPSky ephemeris service for objects near the expected detector position

        Parameters
        ----------
        expCenter : `lsst.geom.SpherePoint`
            Center of search cone
        epochMJD : `float`
            Epoch of cone search, (MJD in UTC).
        queryRadius : `float`
            Radius of the cone search in degrees.
        mpSkyURL : `str`
            URL to query for MPSky.

        Returns
        -------
        mpSkySsObjects : `pandas.DataFrame`
            DataFrame with Solar System Object information and RA/DEC position
            within the visit.
        """
        fieldRA = expCenter.getRa().asDegrees()
        fieldDec = expCenter.getDec().asDegrees()

        params = {
            "t": epochMJD,
            "ra": fieldRA,
            "dec": fieldDec,
            "radius": queryRadius
        }

        try:
            response = requests.get(mpSkyURL, params=params, timeout=self.config.mpSkyRequestTimeoutSeconds)
            response.raise_for_status()
            objID, ra, dec, object_polynomial, observer_polynomial = self.read_mp_sky_response(response)

            mpSkySsObjects = pd.DataFrame()
            mpSkySsObjects['ObjID'] = objID
            mpSkySsObjects['ra'] = ra
            mpSkySsObjects['obj_poly'] = object_polynomial
            mpSkySsObjects['obs_poly'] = observer_polynomial
            mpSkySsObjects['dec'] = dec
            mpSkySsObjects['Err(arcsec)'] = 2
            mpSkySsObjects['ssObjectId'] = [abs(hash(v)) for v in mpSkySsObjects['ObjID'].values]
            nFound = len(mpSkySsObjects)

            if nFound == 0:
                self.log.info("No Solar System objects found for visit.")
            else:
                self.log.info("%d Solar System Objects in visit", nFound)
        except requests.RequestException as e:
            raise NoWorkFound(f"Failed to connect to the remote ephemerides service: {e}") from e

        return mpSkySsObjects
