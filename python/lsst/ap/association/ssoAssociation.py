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

"""Spatial association for Solar System Objects."""

__all__ = ["SolarSystemAssociationConfig", "SolarSystemAssociationTask"]

from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
import healpy as hp
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.utils.timer import timeMethod


class SolarSystemAssociationConfig(pexConfig.Config):
    """Config class for SolarSystemAssociationTask.
    """
    maxDistArcSeconds = pexConfig.Field(
        dtype=float,
        doc='Maximum distance in arcseconds to test for a DIASource to be a '
            'match to a SSObject.',
        default=2.0,
    )
    maxPixelMargin = pexConfig.RangeField(
        doc="Maximum padding to add to the ccd bounding box before masking "
            "SolarSystem objects to the ccd footprint. The bounding box will "
            "be padded by the minimum of this number or the max uncertainty "
            "of the SolarSystemObjects in pixels.",
        dtype=int,
        default=100,
        min=0,
    )


class SolarSystemAssociationTask(pipeBase.Task):
    """Associate DIASources into existing SolarSystem Objects.

    This task performs the association of detected DIASources in a visit
    with known solar system objects.
    """
    ConfigClass = SolarSystemAssociationConfig
    _DefaultName = "ssoAssociation"

    @timeMethod
    def run(self, diaSourceCatalog, solarSystemObjects, exposure):
        """Create a searchable tree of unassociated DiaSources and match
        to the nearest ssoObject.

        Parameters
        ----------
        diaSourceCatalog : `pandas.DataFrame`
            Catalog of DiaSources. Modified in place to add ssObjectId to
            successfully associated DiaSources.
        solarSystemObjects : `pandas.DataFrame`
            Set of solar system objects that should be within the footprint
            of the current visit.
        exposure : `lsst.afw.image.ExposureF`
            Exposure where the DiaSources in ``diaSourceCatalog`` were
            detected in.

        Returns
        -------
        resultsStruct : `lsst.pipe.base.Struct`

            - ``ssoAssocDiaSources`` : DiaSources that were associated with
              solar system objects in this visit. (`pandas.DataFrame`)
            - ``unAssocDiaSources`` : Set of DiaSources that were not
              associated with any solar system object. (`pandas.DataFrame`)
            - ``nTotalSsObjects`` : Total number of SolarSystemObjects
              contained in the CCD footprint. (`int`)
            - ``nAssociatedSsObjects`` : Number of SolarSystemObjects
              that were associated with DiaSources. (`int`)
            - ``ssSourceData`` : ssSource table data. (`Astropy.table.Table`)
        """
        nSolarSystemObjects = len(solarSystemObjects)
        if nSolarSystemObjects <= 0:
            return self._return_empty(diaSourceCatalog, solarSystemObjects)

        mjd_midpoint = exposure.visitInfo.date.toAstropy().tai.mjd
        ref_time = mjd_midpoint - solarSystemObjects["tmin"].values[0]

        solarSystemObjects['obs_position'] = solarSystemObjects.apply(lambda row: np.array([
            np.polynomial.chebyshev.chebval(ref_time, row['obs_x_poly']),
            np.polynomial.chebyshev.chebval(ref_time, row['obs_y_poly']),
            np.polynomial.chebyshev.chebval(ref_time, row['obs_z_poly'])
        ]), axis=1)

        solarSystemObjects['obj_position'] = solarSystemObjects.apply(lambda row: np.array([
            np.polynomial.chebyshev.chebval(ref_time, row['obj_x_poly']),
            np.polynomial.chebyshev.chebval(ref_time, row['obj_y_poly']),
            np.polynomial.chebyshev.chebval(ref_time, row['obj_z_poly'])
        ]), axis=1)

        vector = np.vstack(solarSystemObjects['obj_position'].values
                           - solarSystemObjects['obs_position'].values)
        solarSystemObjects[['ra', 'dec']] = np.vstack(hp.vec2ang(vector, lonlat=True)).T
        solarSystemObjects['obs_position_x'], solarSystemObjects['obs_position_y'], \
            solarSystemObjects['obs_position_z'] = np.array(list(solarSystemObjects['obs_position'].values)).T
        solarSystemObjects['obj_position_x'], solarSystemObjects['obj_position_y'], \
            solarSystemObjects['obj_position_z'] = np.array(list(solarSystemObjects['obj_position'].values)).T

        maskedObjects = self._maskToCcdRegion(
            solarSystemObjects,
            exposure,
            solarSystemObjects["Err(arcsec)"].max()).copy()
        nSolarSystemObjects = len(maskedObjects)
        if nSolarSystemObjects <= 0:
            return self._return_empty(diaSourceCatalog, maskedObjects)

        maxRadius = np.deg2rad(self.config.maxDistArcSeconds / 3600)

        # Transform DIA RADEC coordinates to unit sphere xyz for tree building.
        vectors = self._radec_to_xyz(diaSourceCatalog["ra"],
                                     diaSourceCatalog["dec"])

        # Create KDTree of DIA sources
        tree = cKDTree(vectors)

        nFound = 0
        # Query the KDtree for DIA nearest neighbors to SSOs. Currently only
        # picks the DiaSource with the shortest distance. We can do something
        # fancier later.
        ssSourceData = []
        ras, decs, residual_ras, residual_decs, dia_ids = [], [], [], [], []
        diaSourceCatalog["ssObjectId"] = 0
        for index, ssObject in maskedObjects.iterrows():
            ssoVect = self._radec_to_xyz(ssObject["ra"], ssObject["dec"])
            # Which DIA Sources fall within r?
            dist, idx = tree.query(ssoVect, distance_upper_bound=maxRadius)
            if len(idx) == 1 and np.isfinite(dist[0]):
                nFound += 1
                diaSourceCatalog.loc[diaSourceCatalog.index[idx[0]], "ssObjectId"] = ssObject["ssObjectId"]
                ssSourceData.append(ssObject[["ssObjectId", "obs_position_x", "obs_position_y",
                                              "obs_position_z", "obj_position_x", "obj_position_y",
                                              "obj_position_z"]].values)
                dia_ra = diaSourceCatalog.loc[diaSourceCatalog.index[idx[0]], "ra"]
                dia_dec = diaSourceCatalog.loc[diaSourceCatalog.index[idx[0]], "dec"]
                dia_id = diaSourceCatalog.loc[diaSourceCatalog.index[idx[0]], "diaSourceId"]
                ras.append(dia_ra)
                decs.append(dia_dec)
                dia_ids.append(dia_id)
                residual_ras.append(dia_ra - ssObject["ra"])
                residual_decs.append(dia_dec - ssObject["dec"])
                maskedObjects.loc[index, 'associated'] = True
            else:
                maskedObjects.loc[index, 'associated'] = False

        self.log.info("Successfully associated %d / %d SolarSystemObjects.", nFound, nSolarSystemObjects)
        self.metadata['nAssociatedSsObjects'] = nFound
        self.metadata['nExpectedSsObjects'] = nSolarSystemObjects
        assocSourceMask = diaSourceCatalog["ssObjectId"] != 0
        unAssocObjectMask = np.logical_not(maskedObjects['associated'].values)
        ssSourceData = pd.DataFrame(ssSourceData, columns=["ssObjectId", "obs_position_x", "obs_position_y",
                                                           "obs_position_z", "obj_position_x",
                                                           "obj_position_y", "obj_position_z"])
        ssSourceData["ra"] = ras
        ssSourceData["dec"] = decs
        ssSourceData["residualRa"] = residual_ras
        ssSourceData["residualDec"] = residual_decs
        ssSourceData["diaSourceId"] = dia_ids
        coords = SkyCoord(ra=ssSourceData['ra'].values * u.deg, dec=ssSourceData['dec'].values * u.deg)
        ssSourceData['galacitcL'] = coords.galactic.l.deg
        ssSourceData['galacticB'] = coords.galactic.b.deg
        ssSourceData['eclipticLambda'] = coords.barycentrictrueecliptic.lon.deg
        ssSourceData['eclipticBeta'] = coords.barycentrictrueecliptic.lat.deg
        unassociatedObjects = maskedObjects[unAssocObjectMask].drop(columns=['obs_x_poly', 'obs_y_poly',
                                                                             'obs_z_poly', 'obj_x_poly',
                                                                             'obj_y_poly', 'obj_z_poly',
                                                                             'obs_position', 'obj_position',
                                                                             'associated'])
        return pipeBase.Struct(
            ssoAssocDiaSources=diaSourceCatalog[assocSourceMask].reset_index(drop=True),
            unAssocDiaSources=diaSourceCatalog[~assocSourceMask].reset_index(drop=True),
            nTotalSsObjects=nSolarSystemObjects,
            nAssociatedSsObjects=nFound,
            associatedSsSources=Table.from_pandas(ssSourceData),
            unassociatedSsObjects=Table.from_pandas(unassociatedObjects))

    def _maskToCcdRegion(self, solarSystemObjects, exposure, marginArcsec):
        """Mask the input SolarSystemObjects to only those in the exposure
        bounding box.

        Parameters
        ----------
        solarSystemObjects : `pandas.DataFrame`
            SolarSystemObjects to mask to ``exposure``.
        exposure : `lsst.afw.image.ExposureF`
            Exposure to mask to.
        marginArcsec : `float`
            Maximum possible matching radius to pad onto the exposure bounding
            box. If greater than ``maxPixelMargin``, ``maxPixelMargin`` will
            be used.

        Returns
        -------
        maskedSolarSystemObjects : `pandas.DataFrame`
            Set of SolarSystemObjects contained within the exposure bounds.
        """
        if len(solarSystemObjects) == 0:
            return solarSystemObjects
        wcs = exposure.getWcs()
        padding = min(
            int(np.ceil(marginArcsec / wcs.getPixelScale(exposure.getBBox().getCenter()).asArcseconds())),
            self.config.maxPixelMargin)

        return solarSystemObjects[exposure.containsSkyCoords(
            solarSystemObjects['ra'].to_numpy() * u.degree,
            solarSystemObjects['dec'].to_numpy() * u.degree,
            padding)]

    def _radec_to_xyz(self, ras, decs):
        """Convert input ra/dec coordinates to spherical unit-vectors.

        Parameters
        ----------
        ras : `array-like`
            RA coordinates of objects in degrees.
        decs : `array-like`
            DEC coordinates of objects in degrees.

        Returns
        -------
        vectors : `numpy.ndarray`, (N, 3)
            Output unit-vectors
        """
        ras = np.radians(ras)
        decs = np.radians(decs)
        try:
            vectors = np.empty((len(ras), 3))
        except TypeError:
            vectors = np.empty((1, 3))

        sin_dec = np.sin(np.pi / 2 - decs)
        vectors[:, 0] = sin_dec * np.cos(ras)
        vectors[:, 1] = sin_dec * np.sin(ras)
        vectors[:, 2] = np.cos(np.pi / 2 - decs)

        return vectors

    def _return_empty(self, diaSourceCatalog, emptySolarSystemObjects):
        """Return a struct with all appropriate empty values for no SSO associations.

        Parameters
        ----------
        diaSourceCatalog : `pandas.DataFrame`
            Used for column names
        emptySolarSystemObjects : `pandas.DataFrame`
            Used for column names.
        Returns
        -------
        results : `lsst.pipe.base.Struct`
            Results struct with components.
            - ``ssoAssocDiaSources`` : Empty. (`pandas.DataFrame`)
            - ``unAssocDiaSources`` : Input DiaSources. (`pandas.DataFrame`)
            - ``nTotalSsObjects`` : Zero. (`int`)
            - ``nAssociatedSsObjects`` : Zero.
            - ``associatedSsSources`` : Empty. (`Astropy.table.Table`)
            - ``unassociatedSsObjects`` : Empty. (`Astropy.table.Table`)


        Raises
        ------
        RuntimeError
            Raised if duplicate DiaObjects or duplicate DiaSources are found.
        """
        self.log.info("No SolarSystemObjects found in detector bounding box.")
        return pipeBase.Struct(
            ssoAssocDiaSources=pd.DataFrame(columns=diaSourceCatalog.columns),
            unAssocDiaSources=diaSourceCatalog,
            nTotalSsObjects=0,
            nAssociatedSsObjects=0,
            associatedSsSources=Table(names=["ssObjectId", "ra", "dec", "obs_position", "obj_position",
                                      "residual_ras", "residual_decs"]),
            unassociatedSsObjects=Table(names=emptySolarSystemObjects.columns)
        )
