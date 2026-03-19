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
import pyarrow as pa
import requests

from astropy.table import Table, MaskedColumn
import astropy.units as u

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

# Units for orbital element columns (MPC orbits schema)
ELEMENTS_MPC_ORBITS_UNITS = {
    'a': u.AU, 'q': u.AU, 'e': u.dimensionless_unscaled,
    'i': u.deg, 'node': u.deg, 'argperi': u.deg,
    'peri_time': u.d, 'mean_anomaly': u.deg,
    'period': u.d, 'mean_motion': u.deg/u.d,
    'a_unc': u.AU, 'q_unc': u.AU, 'e_unc': u.dimensionless_unscaled,
    'i_unc': u.deg, 'node_unc': u.deg, 'argperi_unc': u.deg,
    'peri_time_unc': u.d, 'mean_anomaly_unc': u.deg,
    'period_unc': u.d, 'mean_motion_unc': u.deg/u.d,
    'earth_moid': u.AU, 'epoch_mjd': u.d,
    'arc_length_total': u.d, 'arc_length_sel': u.d,
    'h': u.mag, 'not_normalized_rms': u.arcsec, 'normalized_rms': u.dimensionless_unscaled,
}

# Units for classic orbital elements schema
ELEMENTS_UNITS = {
    'q': u.AU, 'e': u.dimensionless_unscaled,
    'incl': u.deg, 'node': u.deg, 'peri': u.deg,
    't_p': u.d, 'epoch': u.d,
}


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

            - ``ssObjects`` : `astropy.table.Table`
                Table containing Solar System Objects near the detector
                footprint as retrieved by MPSky. The columns are as follows:

                ``Name``
                    object name (`str`)
                ``ra``
                    RA in decimal degrees (`float`, degrees)
                ``dec``
                    DEC in decimal degrees (`float`, degrees)
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
            if 'mpc_orb_jsonb' in elements.colnames:
                for element in ELEMENTS_LIST_MPC_ORBITS:
                    mpSkySsObjects[self.config.mpcorb_prefix + element] = elements[element]
            else:
                for element in ELEMENTS_LIST:
                    mpSkySsObjects[self.config.mpcorb_prefix + NAMING_MAP[element]] = elements[element]

        # Fill masked values with NaN for float columns (equivalent to DataFrame.fillna)
        for colname in mpSkySsObjects.colnames:
            col = mpSkySsObjects[colname]
            if hasattr(col, 'filled') and col.dtype.kind == 'f':
                mpSkySsObjects[colname] = col.filled(np.nan)

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
        mpSkySsObjects : `astropy.table.Table`
            Table with Solar System Object information and RA/DEC position
            within the visit.
        elements : `astropy.table.Table` or `None`
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

            # Construct an elements astropy Table
            if "q" in schema.names:
                if "mpc_orb_jsonb" in schema.names:
                    colnames = ELEMENTS_LIST_MPC_ORBITS
                    units_map = ELEMENTS_MPC_ORBITS_UNITS
                else:
                    colnames = ELEMENTS_LIST
                    units_map = {NAMING_MAP[k]: v for k, v in ELEMENTS_UNITS.items()}

                # Columns that are timestamps/datetimes — store as masked strings
                DATETIME_COLS = {'created_at', 'fitting_datetime', 'updated_at'}
                # Columns that should be stored as plain strings (object dtype in arrow)
                STRING_COLS = {'designation', 'id', 'packed_primary_provisional_designation',
                               'unpacked_primary_provisional_designation', 'mpc_orb_jsonb'}

                elements = Table()
                for col in colnames:
                    data = r[col].to_numpy(zero_copy_only=False)
                    if col in DATETIME_COLS:
                        # Store as masked empty strings — TODO (DM-53952): handle as datetime
                        elements[col] = MaskedColumn(
                            np.full(len(data), "", dtype=str), mask=np.ones(len(data), dtype=bool)
                        )
                    elif col in STRING_COLS:
                        # Ensure clean string dtype, not object
                        elements[col] = np.array(data, dtype=str)
                    else:
                        unit = units_map.get(col, None)
                        if unit is None:
                            elements[col] = MaskedColumn(data)
                        else:
                            elements[col] = MaskedColumn(data, unit=unit)
            else:
                elements = None

            object_polynomial = p.to_numpy()
            observer_polynomial = op.to_numpy()

            n = len(r["name"])
            vmag_data = r["Vmag"].to_numpy() if 'Vmag' in r else np.full(n, np.nan)

            mpSkySsObjects = Table({
                'ObjID': r["name"].to_numpy(zero_copy_only=False).astype(str),
                'ra': MaskedColumn(r["ra"].to_numpy(), unit=u.deg),
                'dec': MaskedColumn(r["dec"].to_numpy(), unit=u.deg),
                'tmin': MaskedColumn(r["tmin"].to_numpy(), unit=u.d),
                'tmax': MaskedColumn(r["tmax"].to_numpy(), unit=u.d),
                'trailedSourceMagTrue': MaskedColumn(vmag_data, unit=u.mag),
                'obj_x_poly': [poly[0] for poly in object_polynomial.T],
                'obj_y_poly': [poly[1] for poly in object_polynomial.T],
                'obj_z_poly': [poly[2] for poly in object_polynomial.T],
                'obs_x_poly': [observer_polynomial.T[0] for _ in range(n)],
                'obs_y_poly': [observer_polynomial.T[1] for _ in range(n)],
                'obs_z_poly': [observer_polynomial.T[2] for _ in range(n)],
                'Err(arcsec)': MaskedColumn(np.full(n, 2.0), unit=u.arcsec),
            })

            if 'packed_primary_provisional_designation' in elements.colnames:
                desigs = elements['packed_primary_provisional_designation'].data
            else:
                desigs = mpSkySsObjects['ObjID'].data
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
