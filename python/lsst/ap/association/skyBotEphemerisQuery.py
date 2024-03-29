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

"""Solar System Object Query to Skybot in place of a internal Rubin solar
system object caching/retrieval code.

Will compute the location for of known SSObjects within a visit. This code
blocks on web requests, so should not be used as part of any real-time or
time-sensitive system. Use in a larger pipeline at your own risk.
"""

__all__ = ["SkyBotEphemerisQueryConfig", "SkyBotEphemerisQueryTask"]


from hashlib import blake2b
import pandas as pd
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


class SkyBotEphemerisQueryConnections(PipelineTaskConnections,
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
            "SkyBoy",
        name="visitSsObjects",
        storageClass="DataFrame",
        dimensions=("instrument", "visit"),
    )


class SkyBotEphemerisQueryConfig(
        PipelineTaskConfig,
        pipelineConnections=SkyBotEphemerisQueryConnections):
    observerCode = pexConfig.Field(
        dtype=str,
        doc="IAU Minor Planet Center observer code for queries "
            "(Rubin Obs./LSST default is I11)",
        default='I11'
    )
    queryRadiusDegrees = pexConfig.Field(
        dtype=float,
        doc="On sky radius for Ephemeris cone search. Also limits sky "
            "position error in ephemeris query. Defaults to the radius of "
            "Rubin Obs FoV in degrees",
        default=1.75)


class SkyBotEphemerisQueryTask(PipelineTask):
    """Tasks to query the SkyBot service and retrieve the solar system objects
    that are observable within the input visit.
    """
    ConfigClass = SkyBotEphemerisQueryConfig
    _DefaultName = "SkyBotEphemerisQuery"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        inputs["visit"] = butlerQC.quantum.dataId["visit"]

        outputs = self.run(**inputs)

        butlerQC.put(outputs, outputRefs)

    @timeMethod
    def run(self, visitInfos, visit):
        """Parse the information on the current visit and retrieve the
        observable solar system objects from SkyBot.

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
                retrieved by SkyBot. The columns are as follows; for more
                details see
                https://ssp.imcce.fr/webservices/skybot/api/conesearch/#output-results

                ``Num``
                    object number (`int`, optional)
                ``Name``
                    object name (`str`)
                ``RA(h)``
                    RA in HMS (`str`)
                ``DE(deg)``
                    DEC in DMS (`str`)
                ``Class``
                    Minor planet classification (`str`)
                ``Mv``
                    visual magnitude (`float`)
                ``Err(arcsec)``
                    position error (`float`)
                ``d(arcsec)``
                    distance from exposure boresight (`float`)?
                ``dRA(arcsec/h)``
                    proper motion in RA (`float`)
                ``dDEC(arcsec/h)``
                    proper motion in DEC (`float`)
                ``Dg(ua)``
                    geocentric distance (`float`)
                ``Dh(ua)``
                    heliocentric distance (`float`)
                ``Phase(deg)``
                    phase angle (`float`)
                ``SunElong(deg)``
                    solar elongation (`float`)
                ``ra``
                    RA in decimal degrees (`float`)
                ``dec``
                    DEC in decimal degrees (`float`)
                ``ssObjectId``
                    unique minor planet ID for internal use (`int`). Shared
                    across catalogs; the pair ``(ssObjectId, visitId)`` is
                    globally unique.
                ``visitId``
                    a copy of ``visit`` (`int`)
        """
        # Grab the visitInfo from the raw to get the information needed on the
        # full visit.
        visitInfo = visitInfos[0].get()

        # Midpoint time of the exposure in JD
        expMidPointEPOCH = visitInfo.date.get(system=DateTime.JD, scale=DateTime.UTC)

        # Boresight of the exposure on sky.
        expCenter = visitInfo.boresightRaDec

        # Skybot service query
        skybotSsObjects = self._skybotConeSearch(
            expCenter,
            expMidPointEPOCH,
            self.config.queryRadiusDegrees)

        # Add the visit as an extra column.
        skybotSsObjects['visitId'] = visit

        return pipeBase.Struct(
            ssObjects=skybotSsObjects,
        )

    def _skybotConeSearch(self, expCenter, epochJD, queryRadius):
        """Query IMCCE SkyBot ephemeris service for cone search using the
        exposure boresight.

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
        dfSSO : `pandas.DataFrame`
            DataFrame with Solar System Object information and RA/DEC position
            within the visit.
        """

        fieldRA = expCenter.getRa().asDegrees()
        fieldDec = expCenter.getDec().asDegrees()
        observerMPCId = self.config.observerCode
        radius = queryRadius
        orbitUncertaintyFilter = queryRadius

        # TODO: DM-31866
        q = ['http://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php?']
        q.append('-ra=' + str(fieldRA))
        q.append('&-dec=' + str(fieldDec))
        q.append('&-rd=' + str(radius))
        q.append('&-ep=' + str(epochJD))
        q.append('&-loc=' + observerMPCId)
        q.append('&-filter=' + str(orbitUncertaintyFilter))
        q.append('&-objFilter=111&-refsys=EQJ2000&-output=obs&-mime=text')
        query = ''.join(q)

        result = requests.request("GET", query)
        dfSSO = pd.read_csv(StringIO(result.text), sep='|', skiprows=2)
        if len(dfSSO) > 0:
            # Data has leading and trailing spaces hence the strip.
            columns = [col.strip() for col in dfSSO.columns]
            coldict = dict(zip(dfSSO.columns, columns))
            dfSSO.rename(columns=coldict, inplace=True)
            # Data returned in hourangle format.
            dfSSO["ra"] = Angle(dfSSO["RA(h)"], unit=u.hourangle).deg
            dfSSO["dec"] = Angle(dfSSO["DE(deg)"], unit=u.deg).deg
            # SkyBot returns a string name for the object. To store the id in
            # the Apdb we convert this string to an int by hashing the object
            # name. This is a stop gap until such a time as the Rubin
            # Ephemeris system exists and we create our own Ids. Use blake2b
            # is it can produce hashes that can fit in a 64bit int.
            dfSSO["ssObjectId"] = [
                int(blake2b(bytes(name, "utf-8"), digest_size=7).hexdigest(),
                    base=16)
                for name in dfSSO["Name"]
            ]
        else:
            self.log.info("No Solar System objects found for visit.")
            return pd.DataFrame()

        nFound = len(dfSSO)
        self.log.info("%d Solar System Objects in visit", nFound)

        return dfSSO
