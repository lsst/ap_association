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

Will compute the location for of known SSObjects within a visit. This code
blocks on web requests, so should not be used as part of any real-time or
time-sensitive system. Use in a larger pipeline at your own risk.
"""

__all__ = ["MPSkyEphemerisQueryConfig", "MPSkyEphemerisQueryTask"]


import pandas as pd
import mpsky
import requests
from io import StringIO

from astropy.coordinates import Angle
from astropy import units as u

from lsst.daf.base import DateTime
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.utils.timer import timeMethod

from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections
import lsst.pipe.base.connectionTypes as connTypes

# Enforce an error for unsafe column/array value setting in pandas.
pd.options.mode.chained_assignment = 'raise'


class MPSkyEphemerisQueryConnections(PipelineTaskConnections,
                                      dimensions=("instrument",
                                                  "visit")):
    visitInfos = connTypes.Input(
        doc="Information defining the visit on a per detector basis.",
        name="raw.visitInfo",
        storageClass="VisitInfo",
        dimensions=("instrument", "exposure", "detector"),
        deferLoad=True,
        multiple=True
    )
    ssObjects = connTypes.Output(
        doc="Solar System objects observable in this visit retrieved from "
            "MPSky",
        name="visitSsObjects",
        storageClass="DataFrame",
        dimensions=("instrument", "visit"),
    )


class MPSkyEphemerisQueryConfig(
        PipelineTaskConfig,
        pipelineConnections=MPSkyEphemerisQueryConnections):
    observerCode = pexConfig.Field(
        dtype=str,
        doc="IAU Minor Planet Center observer code for queries " 
            "(Rubin Obs./LSST default is X05)",#[]X05? I11? 
        default='X05'
    )
    queryRadiusDegrees = pexConfig.Field(
        dtype=float,
        doc="On sky radius for Ephemeris cone search. Defaults "
            "to the radius of Rubin Obs FoV in degrees",
        default=1.75)


class MPSkyEphemerisQueryTask(PipelineTask):
    """Tasks to query the MPSky service and retrieve the solar system objects
    that are observable within the input visit.
    """
    ConfigClass = MPSkyEphemerisQueryConfig
    _DefaultName = "MPSkyEphemerisQuery"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        inputs["visit"] = butlerQC.quantum.dataId["visit"]

        outputs = self.run(**inputs)

        butlerQC.put(outputs, outputRefs)

    @timeMethod
    def run(self, visitInfos, visit):
        """Parse the information on the current visit and retrieve the
        observable solar system objects from MPSky.

        Parameters
        ----------
        visitInfos : `list` of `lsst.daf.butler.DeferredDatasetHandle`
            Set of visitInfos for raws covered by this visit/exposure. We
            only use the first instance to retrieve the exposure boresight.
        visit : `int`
            Id number of the visit being run.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components:

            - ``ssObjects``: `pandas.DataFrame`
                DataFrame containing Solar System Objects in field of view as
                retrieved by MPSky. The columns are as follows:
                ``Name``
                    object name (`str`)
                ``ra``
                    RA in decimal degrees (`float`)
                ``dec``
                    DEC in decimal degrees (`float`)
                ``obj_poly``
                    DO NOT USE until t_min issue is resolved
                ``obs_poly``
                    DO NOT USE until t_min issue is resolved



        """
        # Grab the visitInfo from the raw to get the information needed on the
        # full visit.
        visitInfo = visitInfos[0].get()

        # Midpoint time of the exposure in JD
        expMidPointEPOCH = visitInfo.date.get(system=DateTime.JD, scale=DateTime.UTC)

        # Boresight of the exposure on sky.
        expCenter = visitInfo.boresightRaDec

        # MPSky service query
        MPSkySsObjects= self._MPSkyConeSearch(self, expCenter, expMidPointEPOCH, self.config.queryRadiusDegrees) 
        # Add the visit as an extra column.
        MPSkySsObjects['visitId'] = visit

        return pipeBase.Struct(
            ssObjects=MPSkySsObjects,
        )

    def _MPSkyConeSearch(self, expCenter, epochJD, queryRadius):
        """Query MPSky ephemeris service using the exposure boresight.

        Parameters
        ----------
        expCenter : `lsst.geom.SpherePoint`
            Center of Exposure RADEC [deg]
        epochJD : `float`
            Mid point JD of exposure, in UTC [EPOCH].
        queryRadius : `float`
            Radius of the cone search in degrees.

        Returns
        -------
        MPSkySsObjects : `pandas.DataFrame`
            DataFrame with Solar System Object information and RA/DEC position
            within the visit.
        """

        fieldRA = expCenter.getRa().asDegrees()
        fieldDec = expCenter.getDec().asDegrees()
        observerMPCId = self.config.observerCode

        ObjID, ra, dec, obj_poly, obs_poly  = mpsky.query_service('https://sky.dirac.dev/ephemerides/', expMidPointEPOCH, expCenter, queryRadius)#query MPSky[]
        MPSkySsObjects = pd.DataFrame()
        MPSkySsObjects['ObjID'] = ObjID
        MPSkySsObjects['ra'] = ra
        MPSkySsObjects['dec'] = dec
        MPSkySsObjects['obj_poly'] = obj_poly
        MPSkySsObjects['obs_poly'] = obs_poly

        nFound = len(MPSkySsObjects)

        if nFound == 0:
            self.log.info("No Solar System objects found for visit.")

        self.log.info("%d Solar System Objects in visit", nFound)

        return MPSkySsObjects
