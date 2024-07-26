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

__all__ = ("SatelliteFilterTask", "SatelliteFilterConfig")

import numpy as np
import math

import lsst.sphgeom as sphgeom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.utils.timer import timeMethod


class SatelliteFilterConfig(pexConfig.Config):
    """Config class for TrailedSourceFilterTask.
    """

    psf_multiplier = pexConfig.Field(
        dtype=float,
        doc="Multiply the psf by this value.",
        default=2.0,
    )


class SatelliteFilterTask(pipeBase.Task):

    ConfigClass = SatelliteFilterConfig
    _DefaultName = "satelliteSourceFilter"

    @timeMethod
    def run(self, dia_sources, psf, sat_coords):

        angles = self._angle_between_points(sat_coords)
        sph_coords = self.sph_sat_coords(dia_sources)
        tracts = self.satellite_tracts(psf, angles, sat_coords)


        trail_mask = self._check_tracts(sph_coords, tracts)

        return pipeBase.Struct(
            diaSources=dia_sources[~trail_mask].reset_index(drop=True))

    def sph_sat_coords(self, dia_sources):

        sphere_coords = []
        for source in dia_sources.iterrows():
            sphere_coords.append(sphgeom.UnitVector3d(sphgeom.LonLat.fromDegrees(source[1]['ra'], source[1]['dec'])))

        return np.array(sphere_coords)

    def _angle_between_points(self, sat_coords):
        dx = sat_coords[:, 1, 0].flatten() - sat_coords[:, 0, 0].flatten()
        dy = sat_coords[:, 1, 1].flatten() - sat_coords[:, 0, 1].flatten()

        # Angle in radians
        angle_array = np.arctan2(dy, dx)

        return angle_array

    def satellite_tracts(self, psf, theta, sat_coords):
        """ Calculate the satellite tracts
        """
        tracts = []

        perp_slopes = -1.0/np.tan(theta)

        corner1 = [sat_coords[:,0,0] + psf * perp_slopes, sat_coords[:,0,1] + psf * perp_slopes]
        corner2 = [sat_coords[:,0,0] - psf * perp_slopes, sat_coords[:,0,1] - psf * perp_slopes]
        corner3 = [sat_coords[:,1,0] + psf * perp_slopes, sat_coords[:,1,1] + psf * perp_slopes]
        corner4 = [sat_coords[:,1,0] - psf * perp_slopes, sat_coords[:,1,1] - psf * perp_slopes]

        for i in range(len(theta)):
            if (np.isfinite(corner1[0][i]) and np.isfinite(corner1[1][i]) and np.isfinite(corner2[0][i])
                    and np.isfinite(corner2[1][i]) and np.isfinite(corner3[0][i])
                    and np.isfinite(corner3[1][i]) and np.isfinite(corner4[0][i]) and np.isfinite(corner4[1][i])):
                print(corner1[0][i], corner1[1][i])
                tract = sphgeom.ConvexPolygon([sphgeom.UnitVector3d(sphgeom.LonLat.fromDegrees(corner1[0][i], corner1[1][i])),
                                               sphgeom.UnitVector3d(
                                                   sphgeom.LonLat.fromDegrees(
                                                       corner2[0][i],
                                                       corner2[1][i])),
                                               sphgeom.UnitVector3d(
                                                   sphgeom.LonLat.fromDegrees(
                                                       corner3[0][i],
                                                       corner3[1][i])),
                                               sphgeom.UnitVector3d(
                                                   sphgeom.LonLat.fromDegrees(
                                                       corner4[0][i],
                                                       corner4[1][i]))
                    ])
                tracts.append(tract)

        return tracts

    def _check_tracts(self, sphere_coords, tracts):
        """ Check if sources in the catalog fall within the calculated
        satellite boundaries. If so, add them to a mask of sources which will
        be dropped.
        """
        sat_mask = []
        for tract in tracts:
            if tract.contains(sphere_coords[:,0],sphere_coords[:,1],sphere_coords[:,2]).any():
                if sat_mask == []:

                    sat_mask = tract.contains(sphere_coords[:,0],sphere_coords[:,1],sphere_coords[:,2])
                else:
                    sat_mask |= tract.contains(sphere_coords[:,0],sphere_coords[:,1],sphere_coords[:,2])

        return sat_mask