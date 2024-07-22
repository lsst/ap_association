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
from lsst.geom import Point, Polygon

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

        ra, dec = self.satellite_ellipse(midpoints, psf, angles, sat_coords)

        trail_mask = self._check_satellites(dia_sources, ra, dec)

        return pipeBase.Struct(
            diaSources=dia_sources[~trail_mask].reset_index(drop=True))

    def _midpoint(sat_coords):
        xm = (sat_coords[:, 1, 0].flatten() + sat_coords[:, 0, 0].flatten()) / 2
        ym = (sat_coords[:, 1, 1].flatten() + sat_coords[:, 0, 1].flatten()) / 2
        return np.array([xm, ym])

    def _angle_between_points(sat_coords):
        # Compute differences
        dx = sat_coords[:, 1, 0].flatten() - sat_coords[:, 0, 0].flatten()
        dy = sat_coords[:, 1, 1].flatten() - sat_coords[:, 0, 1].flatten()

        # Angle in radians
        angle_array = np.arctan2(dy, dx)

        return angle_array

    def satellite_ellipse(center, psf, theta, array):
        # Generate points for ellipse
        t = np.linspace(0, 2 * np.pi, 100)
        t_large = np.repeat([t], len(theta), axis=0)

        # The semi major axis with the psf buffer times two so the ellipse
        # doesn't terminate exactly at the endpoints
        a = np.linalg.norm(array[:, 0] - array[:, 1], axis=1) / 2.0 + 2.0*psf
        # Semi minor axis that is the width of the psf
        b = psf


        # Rotate to match with the angle. Need to ask about
        # how much satellites curve in images.
        x = center[0, :, np.newaxis] + a[:, np.newaxis] * np.cos(
            theta[:, np.newaxis]) * np.cos(t_large) - b * np.sin(
            theta[:, np.newaxis]) * np.sin(t_large)
        y = center[1, :, np.newaxis] + a[:, np.newaxis] * np.sin(
            theta[:, np.newaxis]) * np.cos(t_large) + b * np.cos(
            theta[:, np.newaxis]) * np.sin(t_large)

        return x, y

    def _in_ellipse(self, sat_radec, ra, dec):
        point = Point(sat_radec[0], sat_radec[1])
        for i, k in enumerate(ra):
            polygon = Polygon(list(zip(ra[i], dec[i])))
            if point.within(polygon):
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
