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


import numpy as np
import os
import pandas as pd
import pyarrow as pa
import requests

from lsst.ap.association.utils import getMidpointFromTimespan
from lsst.geom import SpherePoint
import lsst.pex.config as pexConfig
from lsst.pipe.tasks.associationUtils import obj_id_to_ss_object_id
from lsst.utils.timer import timeMethod

from lsst.pipe.base import connectionTypes, NoWorkFound, PipelineTask, \
    PipelineTaskConfig, PipelineTaskConnections, Struct

ELEMENTS_LIST = "q,e,inc,node,argPeri,t_p_MJD_TDB,epochMJD_TDB".split(',')
ELEMENTS_FORMAT = " ".join(["{:> 11f}"] * len(ELEMENTS_LIST))
NAMING_MAP = {'inc': 'incl', 'argPeri': 'peri', 't_p_MJD_TDB': 't_p', 'epochMJD_TDB': 'epoch',
              'q': 'q', 'e': 'e', 'node': 'node'}

ELEMENTS_LIST_MPC_ORBITS = ['designation', 'id', 'packed_primary_provisional_designation',
                            'unpacked_primary_provisional_designation', 'mpc_orb_jsonb',
                            'created_at', 'updated_at', 'orbit_type_int', 'u_param', 'nopp',
                            'arc_length_total', 'arc_length_sel', 'nobs_total', 'nobs_total_sel',
                            'a', 'q', 'e', 'i', 'node', 'argperi', 'peri_time', 'yarkovsky', 'srp',
                            'a1', 'a2', 'a3', 'dt', 'mean_anomaly', 'period', 'mean_motion', 'a_unc',
                            'q_unc', 'e_unc', 'i_unc', 'node_unc', 'argperi_unc', 'peri_time_unc',
                            'yarkovsky_unc', 'srp_unc', 'a1_unc', 'a2_unc', 'a3_unc', 'dt_unc',
                            'mean_anomaly_unc', 'period_unc', 'mean_motion_unc', 'epoch_mjd', 'h', 'g',
                            'not_normalized_rms', 'normalized_rms', 'earth_moid', 'fitting_datetime']


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
        storageClass="ArrowAstropy",
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
        default=10.0
    )
    mpSkyFallbackURL = pexConfig.Field(
        dtype=str,
        doc="mpSky default URL if MP_SKY_URL environment variable unset",
        default="http://sdfiana014.sdf.slac.stanford.edu:3666/ephemerides/"
    )
    mpcorb_prefix = pexConfig.Field(
        dtype=str,
        doc="Column prefix to indicate MPCORB columns for downstream tasks.",
        default="MPCORB_"
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
                ``obj_X_poly``, ``obj_Y_poly``, ``obj_Z_poly``
                    Chebyshev coefficients for object path
                ``obs_X_poly``, ``obs_Y_poly``, ``obs_Z_poly``
                    Chebyshev coefficients for observer path
                ``t_min``
                    Lower time bound for polynomials
                ``t_max``
                    Upper time bound for polynomials
        """
        # Get detector center and timespan midpoint from predictedRegionTime.
        region = predictedRegionTime.region
        timespan = predictedRegionTime.timespan
        expCenter = SpherePoint(region.getBoundingCircle().getCenter())
        expRadius = region.getBoundingCircle().getOpeningAngle().asDegrees()
        expMidPointEPOCH = getMidpointFromTimespan(timespan, allowUnbounded=False).mjd

        # MPSky service query
        try:
            mpSkyURL = os.environ['MP_SKY_URL']
        except KeyError:
            self.log.warning("MP_SKY_URL is not defined. Falling back to %s.", self.config.mpSkyFallbackURL)
            mpSkyURL = self.config.mpSkyFallbackURL
        mpSkySsObjects, elements = self._mpSkyConeSearch(expCenter, expMidPointEPOCH,
                                                         expRadius + self.config.queryBufferRadiusDegrees,
                                                         mpSkyURL)
        if elements is not None:
            if 'mpc_orb_jsonb' in elements.columns:
                for element in ELEMENTS_LIST_MPC_ORBITS:
                    mpSkySsObjects[self.config.mpcorb_prefix + element] = elements[element]
            else:
                for element in ELEMENTS_LIST:
                    mpSkySsObjects[self.config.mpcorb_prefix + NAMING_MAP[element]] = elements[element]
        mpSkySsObjects = mpSkySsObjects.fillna(value=np.nan)
        # ssObjectId and mpcDesignation are also in the MPCORB schema so let's add them here.
        return Struct(
            ssObjects=mpSkySsObjects,
        )

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
        elements : `pandas.DataFrame` or `None`
            Orbital elements of ssObjects (to pass to alerts)
        """
        fieldRA = expCenter.getRa().asDegrees()
        fieldDec = expCenter.getDec().asDegrees()

        params = {
            "t": epochMJD,
            "ra": fieldRA,
            "dec": fieldDec,
            "radius": queryRadius,
            "return_elements": "extended"
        }

        response = None
        try:
            response = requests.get(mpSkyURL, params=params, timeout=self.config.mpSkyRequestTimeoutSeconds)
            response.raise_for_status()
            with pa.input_stream(memoryview(response.content)) as fp:
                fp.seek(0)
                p = pa.ipc.read_tensor(fp)
                op = pa.ipc.read_tensor(fp)
                with pa.ipc.open_stream(fp) as reader:
                    schema = reader.schema
                    r = next(reader)

            # construct an elements pandas dataframe
            if "q" in schema.names:
                if "mpc_orb_jsonb" in schema.names:
                    colnames = ELEMENTS_LIST_MPC_ORBITS
                else:
                    colnames = ELEMENTS_LIST
                cols = {col: r[col].to_numpy(zero_copy_only=False) for col in colnames}
                elements = pd.DataFrame(cols)
                for c in cols:
                    if c in ['created_at', 'fitting_datetime', 'updated_at']:
                        elements[c] = None
                        # TODO (DM-53952): Handle properly, probably as a datetime
            else:
                elements = None

            mpSkySsObjects = pd.DataFrame()
            mpSkySsObjects['ObjID'] = r["name"].to_numpy(zero_copy_only=False)
            mpSkySsObjects['ra'] = r["ra"].to_numpy()
            mpSkySsObjects['dec'] = r["dec"].to_numpy()
            mpSkySsObjects['tmin'] = r["tmin"].to_numpy()
            mpSkySsObjects['tmax'] = r["tmax"].to_numpy()

            mpSkySsObjects['trailedSourceMagTrue'] = np.nan
            if 'Vmag' in r:
                mpSkySsObjects['trailedSourceMagTrue'] = r["Vmag"].to_numpy()

            object_polynomial = p.to_numpy()
            observer_polynomial = op.to_numpy()

            mpSkySsObjects['obj_x_poly'] = [poly[0] for poly in object_polynomial.T]
            mpSkySsObjects['obj_y_poly'] = [poly[1] for poly in object_polynomial.T]
            mpSkySsObjects['obj_z_poly'] = [poly[2] for poly in object_polynomial.T]
            mpSkySsObjects['obs_x_poly'] = [observer_polynomial.T[0] for
                                            i in range(len(mpSkySsObjects))]
            mpSkySsObjects['obs_y_poly'] = [observer_polynomial.T[1] for
                                            i in range(len(mpSkySsObjects))]
            mpSkySsObjects['obs_z_poly'] = [observer_polynomial.T[2] for
                                            i in range(len(mpSkySsObjects))]
            mpSkySsObjects['Err(arcsec)'] = 2
            if 'packed_primary_provisional_designation' in elements.columns:
                desigs = elements['packed_primary_provisional_designation'].values
            else:
                desigs = mpSkySsObjects['ObjID'].values
            mpSkySsObjects['ssObjectId'] = [obj_id_to_ss_object_id(v) for v in desigs]

            nFound = len(mpSkySsObjects)

            if nFound == 0:
                self.log.info("No Solar System objects found for visit.")
            else:
                self.log.info("%d Solar System Objects in visit", nFound)
        except requests.RequestException as e:
            body = getattr(response, "text", "<No response>")
            message = f"Query to the remote ephemerides service failed. Got response {body}"
            self.log.error(message)
            raise NoWorkFound(f"{message}: {e}") from e

        return mpSkySsObjects, elements
