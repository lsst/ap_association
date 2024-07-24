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


        midpoints = self._midpoint(sat_coords)
        angles = self._angle_between_points(sat_coords)

        ellipses = self.satellite_ellipse(midpoints, psf, angles, sat_coords)

        trail_mask = self._check_satellites(dia_sources, ellipses)

        return pipeBase.Struct(
            diaSources=dia_sources[~trail_mask].reset_index(drop=True))

    def _midpoint(sat_coords):
        xm = (sat_coords[:, 1, 0].flatten() + sat_coords[:, 0, 0].flatten()) / 2.0
        ym = (sat_coords[:, 1, 1].flatten() + sat_coords[:, 0, 1].flatten()) / 2.0
        return np.array([xm, ym])

    def _angle_between_points(sat_coords):
        dx = sat_coords[:, 1, 0].flatten() - sat_coords[:, 0, 0].flatten()
        dy = sat_coords[:, 1, 1].flatten() - sat_coords[:, 0, 1].flatten()

        # Angle in radians
        angle_array = np.arctan2(dy, dx)

        return angle_array

    def satellite_ellipse(self, center, theta):
        ellipses=[]
        for i in range(len(center)):
            ellipse = sphgeom.Ellipse(
                center=sphgeom.UnitVector3d(
                    lsst.sphgeom._sphgeom.LonLat.fromDegrees(center[i, 0],
                                                             center[i, 1])
                ),
                alpha=sphgeom.Angle.fromDegrees(0),
                beta=sphgeom.Angle.fromDegrees(0),
                orientation=sphgeom.Angle.fromDegrees(theta[i]),
            )

            ellipses.append(ellipse)

        return ellipses

    def _in_ellipse(self, catalog_radec, ellipses):
        # Check if the



        for ellipse in ellipses:
            if ellipse.contains(sat_radec):
                return point.within(polygon)
            else:
                continue

    def _check_satellites(self, catalog, x, y):
        """ Check if sources in the catalog fall within the calculated
        satellite boundaries. If so, add them to a mask of sources which will
        be dropped.
        """
        mask = []

        for k in range(len(catalog)):
            ra = catalog[k]["base_SdssCentroid_x"]
            dec = catalog[k]["base_SdssCentroid_y"]
            check = self._in_ellipse(np.array([ra, dec]), x, y)
            print(check)

            mask.append[k]

        return mask
