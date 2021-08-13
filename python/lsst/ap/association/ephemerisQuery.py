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

"""Solar System Object Query to Ephemeris Service.

Will currently compute the location for of SSObjects for a known set of PVIs.
NOT TO BE RUN AS PART OF A PIPELINE! This is an Ephemeris Query step for
integration only.
"""

__all__ = ["EphemerisQueryConfig", "EphemerisQueryTask"]


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
                                dimensions=("skymap", "instrument"),
                                defaultTemplates={"coaddName": "deep",
                                                  "fakesType": ""}):
    diffIms = connTypes.Input(
        doc="Difference image on which the DiaSources were detected.",
        name="{fakesType}{coaddName}Diff_differenceExp",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector"),
        multiple=True,
        deferLoad=True,
    )
    ssObjects = connTypes.Output(
        doc="Solar System Objects for all difference images in diffIm.",
        name="{fakesType}{coaddName}_ccdVisitSsObjects",
        storageClass="DataFrame",
        dimensions=("skymap", "instrument"),
    )
    # ssoEphemerisDB = cT.Output(
    #     doc="Solar System Object ephemeris database.",
    #     name="ssoEphemerisDB",
    #     storageClass="DataFrame",
    #     dimensions=("tract", "skymap")
    # )


class EphemerisQueryConfig(PipelineTaskConfig,
                           pipelineConnections=EphemerisQueryConnections):

    observerCode = pexConfig.Field(
        dtype=str,
        doc='IAU Minor Planet Center observer code for LSST.',
        default='I11'
    )
    queryRadius = pexConfig.Field(
        dtype=str,
        doc='On sky radius for Ephemeris cone search.'
        'Also limits sky position error in ephemeris query.',
        default=600/3600)


class EphemerisQueryTask(PipelineTask):
    """Fetch latest MPCORB and Solar System Object ephemeris data and compute
    expected SSObjects within a set of difference images.

    Create interpolants in ICRF xyz on the observer unit sphere
    rather than RADEC to facilitate calculations of objects in field
    of view through dot products. Interpolated xyz values can be turned
    into RADEC.
    """
    ConfigClass = EphemerisQueryConfig
    _DefaultName = "EphemerisQuery"

    @pipeBase.timeMethod
    def run(self, diffIms):
        """Load MPCORB.DAT file containing orbits of known Solar System Objects.
        Query Ephemeris service for RADEC positions of Solar System Objects
        in Field of View.

        Parameters
        ----------
        diffIm : `list` of `lst.daf.butler.DeferredDatasetHandle`
            Full set of difference images in a repo to pre-compute SSObject
            locations within.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components:

            - ``ssObjects``: `pandas DataFrame`
                Contains pandas DataFrame with name, RADEC position and position
                uncertainty on sky of Solar System Objects in Field of View.
        """

        # Loop over our list of difference image exposures. For the output
        # DataFrame, my guess is you should try to follow the columns in the
        # DPDD (http://ls.st/dpdd) as closely as possible. If you need extra
        # columns or don't fill some, that's fine but the should follow those
        # naming conventions.
        ssODF = []
        for diffImRef in diffIms:
            # Load the difference image.
            diffIm = diffImRef.get(
                datasetType=self.config.connections.diffIms,
                immediate=True)

            # mid exposure MJD
            expMidPointMJD = diffIm.getInfo().getVisitInfo().getDate().get(system=DateTime.MJD)

            # expTime in seconds.
            # expTime = diffIm.getInfo().getVisitInfo().getExposureTime()

            # The center here is a SpherePoint object. It has a method called
            # `getVector` if you need the unit sphere 3 vector. You can access
            # those values as x, y, z from the object returned by getVector.
            expCenter = diffIm.getWcs().pixelToSky(geom.Box2D(diffIm.getBBox()).getCenter())

            # Make sure to add this value as a column to the output dataFrame.
            exposureId = diffIm.getInfo().getVisitInfo().getExposureId()

            # Ephemeris service query
            skybot = self._skybotConeSearch(expCenter, expMidPointMJD)
            skybot['ccdVisitId'] = exposureId
            ssODF.append(skybot)

            self.log.info(f"Finished ccdVisit {exposureId}")

        ssObjectDataFrame = pd.concat(ssODF)

        return pipeBase.Struct(
            ssObjects=ssObjectDataFrame,
        )

    def _skybotConeSearch(self, expCenter, expMidPointMJD):
        """Query IMCCE SkyBot ephemeris service for cone search to get RADEC
        positions of Solar System Objects in FOV.

        Parameters
        ----------
        expCenter : diffIm.getWcs().pixelToSky(geom.Box2D(diffIm.getBBox()).getCenter()) object
            Center of Exposure RADEC [deg]
        expMidPointMJD : diffIm.getInfo().getVisitInfo().getDate().get(system=DateTime.MJD)
            Mid point time of exposure [MJD]

        Returns
        -------
        dfSSO ... Pandas DataFrame
            DataFrame with Solar System Object information and RADEC positin in Field of View
        """

        fieldRA = expCenter.getRa().asDegrees()
        fieldDec = expCenter.getDec().asDegrees()
        epochJD = expMidPointMJD + 2400000.5
        observerMPCId = self.config.observerCode
        radius = self.config.queryRadius
        orbitUncertaintyFilter = self.config.queryRadius

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
        try:
            r = requests.request("GET", query)
            dfSSO = pd.read_csv(StringIO(r.text), sep='|', skiprows=2)
            columns = [col.strip() for col in conedf.columns]
            coldict = dict(zip(conedf.columns, columns))
            dfSSO.rename(columns=coldict, inplace=True)
            dfSSO['RA(deg)'] = self._rahms2radeg(conedf['RA(h)'])
            dfSSO['DEC(deg)'] = self._decdms2decdeg(conedf['DE(deg)'])
            dfSSO.drop(columns=['RA(h)', 'DE(deg)', 'Dg(ua)', 'Dh(ua)', 'd(arcsec)'], inplace=True)

        except RuntimeError:
            self.log.info("No Solar System objects found for ccdVisit.")
            return pd.DataFrame()

        nFound = len(dfSSO)
        self.log.info(f"{nFound} Solar System Objects in ccdVisit")

        return dfSSO

    def _decdms2decdeg(self, decdms):
        """Convert Declination from degrees minutes seconds to decimal degrees.

        Parameters
        ----------
        decdms ... str
            Declination string "degrees minutes seconds"

        Returns
        -------
        decdeg ... float / array of floats
            Declination in degrees
        """
        decdeg = np.zeros(len(decdms))
        i = 0
        for dec in decdms:
            deglist = [float(d) for d in dec.split()]
            decdeg[i] = deglist[0] + deglist[1]/60 + deglist[2]/3600
            i = i+1
        return decdeg

    def _rahms2radeg(rahms):
        """Convert Right Ascension from hours minutes seconds to decimal degrees.

        Parameters
        ----------
        rahms ... str
            Declination string "hours minutes seconds"

        Returns
        -------
        radeg ... float / array of floats
            Declination in degrees
        """
        radeg = np.zeros(len(rahms))
        i = 0
        for ra in rahms:
            ralist = [float(r) for r in ra.split()]
            radeg[i] = (ralist[0]/24 + ralist[1]/1440 + ralist[2]/86400)*360
            i = i+1
        return radeg
