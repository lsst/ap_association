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

"""Solar System Object Ephemeris Precomputation.

Will currently compute the location for of SSObjects for a known set of PVIs.
NOT TO BE RUN AS PART OF A PIPELINE! This is a pre-computation step for
integration only.
"""

__all__ = ["PreComputeEphemerisConfig", "PreComputeEphemerisTask"]

### for pipeline testing only

from astroquery.imcce import Skybot
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
###


import numpy as np
import pandas as pd
import scipy.interpolate as spi

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


class sso:
    """Class for storing Solar System Object data.
    """
    def __init__(self):
        self.ssoId = []
        self.name = []
        self.interpolant = []


class PreComputeEphemerisConnections(PipelineTaskConnections,
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

class PreComputeEphemerisConfig(PipelineTaskConfig,
                pipelineConnections=PreComputeEphemerisConnections):
    ssoEphemerisTstart = pexConfig.Field(
        dtype=float,
        doc='First/Start epoch for Solar System Object ephemeris generation [JD].',
        default=2459580.5
    )
    ssoEphemerisTstop = pexConfig.Field(
        dtype=float,
        doc='Last/Stop epoch for Solar System Object ephemeris generation [JD].',
        # default=2463232.5
        default=2459590.5
    )
    mpcorbIdxMin = pexConfig.Field(
        dtype=int,
        doc="TODO: Siegfried writes doc",
        default=0,
    )
    mpcorbIdxMax = pexConfig.Field(
        dtype=int,
        doc="TODO: Siegfried writes doc",
        default=-1,
    )
    path2mpcorb = pexConfig.Field(
        dtype=str,
        doc='Path to MPCORB.DAT file containing the latest list of '
        'Solar System objects and their orbits. Can be local or url.',
        default='https://minorplanetcenter.net/iau/MPCORB/MPCORB.DAT'
    )
    observerCode = pexConfig.Field(
        dtype=str,
        doc='IAU Minor Planet Center observer code for LSST.',
        default='I11'
    )
    queryStep =  pexConfig.Field(
        dtype=str,
        doc='Time step for JPL Horizons ephemeris web API query.'
        'The unit has to be separated from the number by one '
        'whitespace encoded as %20. Units are h for hours and d '
        'for days.',

        default='2%20h')


class PreComputeEphemerisTask(PipelineTask):
#class PreComputeEphemerisTask(pipeBase.Task):
    """Fetch latest MPCORB and Solar System Object ephemeris data and compute
    expected SSObjects within a set of difference images.

    Create interpolants in ICRF xyz on the observer unit sphere
    rather than RADEC to facilitate calculations of objects in field
    of view through dot products. Interpolated xyz values can be turned
    into RADEC.
    """
    ConfigClass = PreComputeEphemerisConfig
    _DefaultName = "preComputeEphemeris"

    @pipeBase.timeMethod
    def run(self, diffIms):
        """Load MPCORB.DAT file containing orbits of known Solar System Objects.
        Query JPL Horizons for ephemeris between startJD and endJD epochs
        and create interpolants for unit sphere xyz coordinates in ICRF.
        Object IDs and interpolants are stored in ssoList object.

        Parameters
        ----------
        diffIm : `list` of `lst.daf.butler.DeferredDatasetHandle`
            Full set of difference images in a repo to pre-compute SSObject
            locations within.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components:

            - ``ssoList``: ``sso`` class object
                   Contains lists of name (`str`), ssoId (`int`) and
                   interpolant (`scipy.interpolate.CubicSpline`) for SSOs.

            - ``mpcorbdf``: `pandas.DataFrame`
                    MPCORB.DAT DataFrame containing identifiers and orbit data for SSOs.
        """

        """
        mpcorbIdxRange = [self.config.mpcorbIdxMin, self.config.mpcorbIdxMax]

        # fetch MPCORB file containing orbits of known SSOs and put them into a pandas DataFrame
        mpcorbdf = self.fetchMPCORB(url=self.config.path2mpcorb)

        # fetch a rough grid of ephemerides from JPL Horizons
        # and create ICRF xyz interpolants that can be turned into RADEC values
        startJD = self.config.ssoEphemerisTstart
        endJD = self.config.ssoEphemerisTstop

        ssoList=sso()

        for index, row in mpcorbdf[self.config.mpcorbIdxMin:self.config.mpcorbIdxMax].iterrows():
            name, number = self._convertSSOname(row['ssoName'])
            try:
                [horizonsdf, ssoInterpolant] = self.createSSOInterpolant(name, number)
            except:
                raise RuntimeError(
                    'Could not download ephemeris for'
                    ' the following object:', name)
            ssoList.name.append(name)
            ssoList.ssoId.append(row['ssoId'])
            ssoList.interpolant.append(ssoInterpolant)
        #else:
        #    selmpcorbdf = mpcorbdf[mpcorbIdxRange[0]:mpcorbIdxRange[1]]
        #    newmpcorbdf = sel[~selmpcorbdf['ssoName'].isin(ssoList.name)]
        #    for index, row in newmpcorbdf.iterrows():
        #        name, number = self._convertSSOname(row['ssoName'])
        #        try:
        #            [horizonsdf, ssoInterpolant] = self.createSSOInterpolant(name,number)
        #
        #        except:
        #            raise RuntimeError('Could not download ephemeris for'+
        #           ' the following object:',name)
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

            # MJD in nanoseconds.
            expMidPointMJD = diffIm.getInfo().getVisitInfo().getDate().get(system=DateTime.MJD)

            # expTime in seconds.
            expTime = diffIm.getInfo().getVisitInfo().getExposureTime()

            # The center here is a SpherePoint object. It has a method called
            # `getVector` if you need the unit sphere 3 vector. You can access
            # those values as x, y, z from the object returned by getVector.
            expCenter = diffIm.getWcs().pixelToSky(geom.Box2D(diffIm.getBBox()).getCenter())

            # Make sure to add this value as a column to the output dataFrame.
            exposureId = diffIm.getInfo().getVisitInfo().getExposureId()

            # TODO: Siegfried, test against the splines. Append to output pandas DataFrame if
            # inside exposure.
            skybot = self._testSkybot(expCenter, expMidPointMJD, 600/3600)
            skybot['ccdVisitId'] = exposureId
            ssODF.append(skybot)

            self.log.info(f"Finished ccdVisit {exposureId}")

        ssObjectDataFrame = pd.concat(ssODF)

        return pipeBase.Struct(
            ssObjects=ssObjectDataFrame,
        )


    def _testSkybot(self, expCenter, expMidPointMJD, radius):
        """Extracts SSOs per field in cone search.
        """

        field = SkyCoord(expCenter.getRa().asDegrees()*u.deg, expCenter.getDec().asDegrees()*u.deg)
        epoch = Time(expMidPointMJD, format='mjd',scale='utc')
        try:
            res = Skybot.cone_search(field, radius*u.deg, epoch)
        except RuntimeError:
            self.log.info("No Solar System objects found for ccdVisit.")
            return pd.DataFrame()

        dfSSO = pd.DataFrame()

        dfSSO['name'] = res['Name']
        dfSSO['ra'] = res['RA']
        dfSSO['decl'] = res['DEC']
        nFound = len(dfSSO)
        self.log.info(f"{nFound} Solar System Objects in ccdVisit")

        return dfSSO

    def _convertSSOname(self, name):
        """Convert MPCORB name to JPL Horizons query compatible name/number.

        Parameters:
        -----------
        name: `str`
            MPCORB readable name

        Results:
        --------
        ssoName: `str`
            JPL Horizons query compatible name
        ssoNumber: `str`
               JPL Horizons query compatible number for numbered objects
        """

        ssoName=None
        objectNumber=None

        splt = name.split(') ')

        if (len(splt)>1):
           ssoNumber=(splt[0].split('('))[1]
           ssoName=splt[1]
        else:
           ssoName=splt[0]

        return ssoName, ssoNumber

    #@pipeBase.timeMethod
    def fetchMPCORB(self, url=None, skiprows=None, filternan=None):
        """Fetch IAU Minor planet center MPCORB file from url and converts it to Pandas Data Frame.

        Parameters
        ----------
        url : `str`
            Path/URL for MPCORB.DAT file
        skiprows : `int`
            Number of header lines to skip when reading MPCORB.DAT file.

        Returns
        -------
        mpcorbdf:  `pandas.DataFrame`
            Pandas DataFrame containing MPCORB orbit data.

        """
        if (url is None):
            url='https://minorplanetcenter.net/iau/MPCORB/MPCORB.DAT'
        if (skiprows is None):
            skiprows = 43
        if (filternan is None):
            filternan = True

        mpcorb_col_numbers=[(0,7),(8,13),(14,19),(20,25),(26,35),(37,46),
                (48,57),(59,68),(70,79),(80,91),(92,103),(106,116),
                (117,122),(123,126),(127,136),(137,141),
                (142, 145),(146,149),(150,160),(166,194),(194,202)]
        col_names=['Designation','H','G','epoch','M','argperi','node','i',
                   'e','n','a','reference','N_Obs', 'N_Opp',
                   'yr_of_first_and_last_Obs', 'r.m.s',
                   'coarsePerts', 'precisePerts', 'computer',
                   'ssoName', 'lastObs']

        dtp = [str,float,float,str,float,float,float,float,float,float,float]
        dtypes = dict(zip(col_names,dtp))

        mpcorbdf = pd.read_fwf(url,skiprows=skiprows,colspecs=mpcorb_col_numbers,
                       names=col_names,dytpe=dtypes,index_col=False)
        mpcorbdf['ssoId'] = mpcorbdf.index

        if (filternan):
          mpcorbdf.dropna(subset=['a', 'e','i','node','argperi','M','epoch', 'H'],inplace=True)

        return mpcorbdf

    @pipeBase.timeMethod
    def createSSOInterpolant(self, ssoName, ssoNumber):

        """Create scipy interpolant for position on the unit sphere as seen from an observer.

        Parameters
        -----------
        ssoName : `str`
            JPL Horizons web api compatible name of Solar Sytem Object (e.g. 'Pallas')
        ssoNumber : `str`
            JPL Horizons web api compatible number of Solar System Object (e.g. '1' for Ceres)

        Returns
        -------
        ephemeris   : JPL Horizons ephemeris for a
                                given solar system object (`pandas.DataFrame`)
        interpolant : ICRF xyz [au] unit sphere interpolant (function of JD)
                                Velocities can be interpolated via ipos.derivative(nu=1)
                               (`scipy.interpolate.CubicSpline`)
        """

        tminJD = self.config.ssoEphemerisTstart
        tmaxJD = self.config.ssoEphemerisTstop
        observer = self.config.observerCode
        queryStep = self.config.queryStep

        hheader=['timeJD','solar_presence','lunar_presence',
                     'RA','DEC','dRAcosDEC',
                     'dDEC','APmag','surface_brigthness',
                     'range','range_rate','RA_3sigma','DEC_3sigma','empty1']
        hquant='1,3,9,20,36'

        #Start and stop times of the survey
        tstart = 'JD'+str(tminJD-1.)
        tstop = 'JD'+str(tmaxJD+1.)

        epoch_list={'start':tstart, 'stop':tstop, 'step':queryStep}

        try:
            query = self._makeJPLquery(ssoName,observer,tstart,tstop,queryStep,hquant)
            r = requests.request("GET", query)
            qstring = r.text
            if(not "$$SOE" in r.text):
        # try again with the ssoNumber
                query = self._makeJPLquery(ssoNumber,observer,tstart,tstop,queryStep,hquant)
                r1 = requests.request("GET", query)
                qstring = r1.text
                if(not "$$SOE" in r1.text):
                    raise ValueError('ERROR: JPL Horizons query failed.')
            selstring=(qstring.split("$$SOE",1)[1]).split("$$EOE",1)[0]
            # Query JPL Horizons web api
            ephemeris = pd.read_csv(StringIO(selstring),sep=',',names=hheader).drop(columns=['empty1'])
        # convert RADEC to unit ICRF
            xyz = self._radec2icrf(ephemeris['RA'].values,ephemeris['DEC'].values)
            # create a scipy interpolant for the xyz coordiantes
            interpolant = spi.CubicSpline(ephemeris['timeJD'].values,xyz,axis=1, extrapolate=None)

        except:
            raise ValueError('ERROR: Could not generate ephemeris from JPL Horizons')

        return ephemeris, interpolant

    def _makeJPLquery(self, ssoName, observer, tstart, tstop, queryStep, hquant):
        """
        """
        q = []
        q.append("https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1")
    # the semicolon %3B after the query excludes spacecraft and major bodies with the same name
        q.append("&COMMAND='%20"+str(ssoName)+"%3B'")
        q.append("&CENTER='"+str(observer)+"'")
        q.append("&START_TIME='"+str(tstart)+"'")
        q.append("&STOP_TIME='"+str(tstop)+"'")
        q.append("&STEP_SIZE='"+str(queryStep)+"'")
        q.append("&OBJ_DATA='NO'")
        q.append("&CAL_FORMAT='JD'&ANG_FORMAT='DEG'&MAKE_EPHEM='YES'&TABLE_TYPE='OBSERVER'")
        q.append("&QUANTITIES='"+hquant+"'&CSV_FORMAT='YES'")

        query = ''.join(q)

        return query

    def _icrf2radec(self, pos, deg=None):
        """Convert ICRF xyz to Right Ascension and Declination.
        Geometric states on unit sphere, no light travel time/aberration correction.

        Parameters:
        -----------
        pos : `float` or `np.array` of `float` dimension n x 3 vector
            position on unit sphere (ICRF)
        deg : `boolean`
            True: angles in degrees, False: angles in radians

        Returns:
        --------
        ra : `float` or `np.array` of `float`
            Right Ascension [deg]
        dec : `float` or `np.array` of `float`
            Declination [deg]
        """
        if deg is None:
           deg = True

        norm = np.linalg.norm
        array = np.array
        arctan2 = np.arctan2
        arcsin = np.arcsin
        rad2deg = np.rad2deg
        modulo = np.mod
        pix2 = 2.*np.pi

        if(pos.ndim>1):
            r = norm(pos,axis=1)

            xu = pos[:,0]/r
            yu = pos[:,1]/r
            zu = pos[:,2]/r
        else:
            r = norm(pos)
            xu = pos[0]/r
            yu = pos[1]/r
            zu = pos[2]/r


        phi = arctan2(yu,xu)
        delta = arcsin(zu)

        if(deg):
            ra = modulo(rad2deg(phi)+360,360)
            dec = rad2deg(delta)
        else:
            ra = modulo(phi+pix2,pix2)
            dec = delta

        return ra, dec


    def _radec2icrf(self, ra, dec, deg=None):
        """Convert Right Ascension and Declination to ICRF xyz unit vector.
        Geometric states on unit sphere, no light travel time/aberration correction.

        Parameters:
        -----------
        ra : `float` or `np.array` of `float`
            Right Ascension [deg]
        dec : `float` or `np.array` of `float`
            Declination [deg]
        deg : `boolean`
            True: angles in degrees, False: angles in radians

        Returns:
        --------
        result : `np.array` of `float` dimension n x 3
            3D vector of unit length (ICRF)
        """

        if (deg == None):
            deg = True

        deg2rad = np.deg2rad
        array = np.array
        cos = np.cos
        sin = np.sin

        if(deg):
            a = deg2rad(ra)
            d = deg2rad(dec)
        else:
            a = array(ra)
            d = array(dec)

        cosd = cos(d)
        x = cosd*cos(a)
        y = cosd*sin(a)
        z = sin(d)

        return array([x, y, z])
