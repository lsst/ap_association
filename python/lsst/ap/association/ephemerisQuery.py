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

"""Solar System Object Query to Skybot in place of the Ephemeris task and
internal Rubin data.

Will compute the location for of SSObjects for a know visit. Currently uses
a full external service across the web to fill data.
"""

__all__ = ["EphemerisQueryConfig", "EphemerisQueryTask"]


from hashlib import blake2b
import numpy as np
import pandas as pd
import requests
from io import StringIO

import lsst.geom as geom
from lsst.daf.base import DateTime
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections
import lsst.pipe.base.connectionTypes as connTypes

# Enforce an error for unsafe column/array value setting in pandas.
pd.options.mode.chained_assignment = 'raise'


class EphemerisQueryConnections(PipelineTaskConnections,
                                dimensions=("instrument",
                                            "visit")):
    visitInfos = connTypes.Input(
        doc="Information defining the visit on a per detector basis..",
        name="raw.visitInfo",
        storageClass="VisitInfo",
        dimensions=("instrument", "exposure", "detector"),
        deferLoad=True,
        multiple=True
    )
    ssObjects = connTypes.Output(
        doc="Solar System Objects for all difference images in diffIm.",
        name="visitSsObjects",
        storageClass="DataFrame",
        dimensions=("instrument", "visit"),
    )


class EphemerisQueryConfig(PipelineTaskConfig,
                           pipelineConnections=EphemerisQueryConnections):
    observerCode = pexConfig.Field(
        dtype=str,
        doc='IAU Minor Planet Center observer code for LSST.',
        default='I11'
    )
    queryRadiusDegrees = pexConfig.Field(
        dtype=float,
        doc='On sky radius for Ephemeris cone search.'
        'Also limits sky position error in ephemeris query. '
        'Defaults to the radius of Rubin Obs FoV in degrees',
        default=1.75)


class EphemerisQueryTask(PipelineTask):
    """
    """
    ConfigClass = EphemerisQueryConfig
    _DefaultName = "EphemerisQuery"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        inputs["visit"] = butlerQC.quantum.dataId["visit"]

        outputs = self.run(**inputs)

        butlerQC.put(outputs, outputRefs)

    @pipeBase.timeMethod
    def run(self, visitInfos, visit):
        """

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

            - ``ssObjects``: `pandas DataFrame`
                Contains pandas DataFrame with name, RADEC position, position
                uncertainty, and unique Id of on sky of Solar System Objects in
                Field of View as retrieved by SkyBot.
        """
        # Grab the visitInfo from the raw to get the information needed on the
        # full visit.
        visitInfo = visitInfos[0].get(
            datasetType=self.config.connections.visitInfos,
            immediate=True)

        # Midpoint time of the exposure in MJD
        expMidPointMJD = visitInfo.date.get(system=DateTime.MJD)

        # Boresight of the exposure on sky.
        expCenter = visitInfo.boresightRaDec
        ra = expCenter.getRa().asDegrees()
        decl =  expCenter.getDec().asDegrees()

        # Skybot service query
        skybotSsObjects = self._skybotConeSearch(
            expCenter,
            expMidPointMJD,
            self.config.queryRadiusDegrees)
        # Add the visit as an extra column.
        skybotSsObjects['visitId'] = visit

        return pipeBase.Struct(
            ssObjects=skybotSsObjects,
        )

    def _skybotConeSearch(self, expCenter, expMidPointMJD, queryRadius):
        """Query IMCCE SkyBot ephemeris service for cone search to get RADEC
        positions of Solar System Objects in FOV.

        Parameters
        ----------
        expCenter : `lsst.geom.SpherePoint`
            Center of Exposure RADEC [deg]
        expMidPointMJD : `float`
            Mid point time of exposure [MJD].

        Returns
        -------
        dfSSO : `pandas.DataFrame`
            DataFrame with Solar System Object information and RADEC position
            within the exposure.
        """

        fieldRA = expCenter.getRa().asDegrees()
        fieldDec = expCenter.getDec().asDegrees()
        epochJD = expMidPointMJD + 2400000.5
        observerMPCId = self.config.observerCode
        radius = queryRadius
        orbitUncertaintyFilter = queryRadius

        q = ['http://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php?']
        q.append('-ra=' + str(fieldRA))
        q.append('&-dec=' + str(fieldDec))
        q.append('&-rd=' + str(radius))
        q.append('&-ep=' + str(epochJD))
        q.append('&-loc=' + observerMPCId)
        q.append('&-filter=' + str(orbitUncertaintyFilter))
        q.append('&-objFilter=111&-refsys=EQJ2000&-output=obs&-mime=text')
        query = ''.join(q)

        conedf = pd.DataFrame()
        result = requests.request("GET", query)
        dfSSO = pd.read_csv(StringIO(result.text), sep='|', skiprows=2)
        if len(dfSSO) > 0:
            columns = [col.strip() for col in dfSSO.columns]
            coldict = dict(zip(dfSSO.columns, columns))
            dfSSO.rename(columns=coldict, inplace=True)
            dfSSO["ra"] = self._rahms2radeg(dfSSO["RA(h)"])
            dfSSO["decl"] = self._decdms2decdeg(dfSSO["DE(deg)"])
            # SkyBot returns a string name for the object. To store the id in
            # the Apdb we convert this string to an int by hashing the object
            # name. This is a stop gap until such a time as the Rubin
            # Emphemeris system exists and we create our own Ids.
            dfSSO["ssObjectId"] = [
                int(blake2b(bytes(name, "utf-8"), digest_size=8).hexdigest(),
                    base=16)
                for name in dfSSO["Name"]
            ]

        else:
            self.log.info("No Solar System objects found for visit.")
            return pd.DataFrame()

        nFound = len(dfSSO)
        self.log.info(f"{nFound} Solar System Objects in visit")

        return dfSSO

    def _decdms2decdeg(self, decdms):
        """Convert Declination from degrees minutes seconds to decimal degrees.

        Parameters
        ----------
        decdms : `list` of `str`
            Declination string "degrees minutes seconds"

        Returns
        -------
        decdeg : `numpy.ndarray`, (N,)
            Declination in degrees.
        """
        decdeg = np.empty(len(decdms))
        for idx, dec in enumerate(decdms):
            deglist = [float(d) for d in dec.split()]
            decdeg[idx] = deglist[0] + deglist[1]/60 + deglist[2]/3600
        return decdeg

    def _rahms2radeg(self, rahms):
        """Convert Right Ascension from hours minutes seconds to decimal
        degrees.

        Parameters
        ----------
        rahms : `list` of `str`
            Declination string "hours minutes seconds"
        Returns
        -------
        radeg : `numpy.ndarray`, (N,)
            Right Ascension in degrees
        """
        radeg = np.empty(len(rahms))
        for idx, ra in enumerate(rahms):
            ralist = [float(r) for r in ra.split()]
            radeg[idx] = (ralist[0]/24 + ralist[1]/1440 + ralist[2]/86400)*360
        return radeg
