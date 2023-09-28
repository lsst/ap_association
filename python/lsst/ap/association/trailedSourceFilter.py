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

__all__ = ("TrailedSourceFilterTask", "TrailedSourceFilterConfig")

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.utils.timer import timeMethod


class TrailedSourceFilterConfig(pexConfig.Config):
    """Config class for TrailedSourceFilterTask.
    """

    max_trail_length = pexConfig.Field(
        dtype=float,
        doc="Length of long trailed sources to remove from the input catalog, "
            "in arcseconds per second. Default comes from DMTN-199, which "
            "requires removal of sources with trails longer than 10 "
            "degrees/day, which is 36000/3600/24 arcsec/second, or roughly"
            "0.416 arcseconds per second.",
        default=36000/3600.0/24.0,
    )


class TrailedSourceFilterTask(pipeBase.Task):
    """Find trailed sources in DIASources and filter them as per DMTN-199
    guidelines.

    This task checks the length of trailLength in the DIASource catalog using
    a given arcsecond/second rate from max_trail_length and the exposure time.
    The two values are used to calculate the maximum allowed trail length and
    filters out any trail longer than the maximum. The max_trail_length is
    outlined in DMTN-199 and determines the default value.
    """

    ConfigClass = TrailedSourceFilterConfig
    _DefaultName = "trailedSourceFilter"

    @timeMethod
    def run(self, dia_sources, exposure_time):
        """Remove trailed sources longer than ``config.max_trail_length`` from
        the input catalog.

        Parameters
        ----------
        dia_sources : `pandas.DataFrame`
            New DIASources to be checked for trailed sources.
        exposure_time : `float`
            Exposure time from difference image.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components.

            - ``dia_sources`` : DIASource table that is free from unwanted
             trailed sources. (`pandas.DataFrame`)

            - ``trailed_dia_sources`` : DIASources that have trails which
            exceed max_trail_length/second*exposure_time.
            (`pandas.DataFrame`)
        """
        trail_mask = self._check_dia_source_trail(dia_sources, exposure_time)

        return pipeBase.Struct(
            diaSources=dia_sources[~trail_mask].reset_index(drop=True),
            trailedDiaSources=dia_sources[trail_mask].reset_index(drop=True))

    def _check_dia_source_trail(self, dia_sources, exposure_time):
        """Find DiaSources that have long trails.

        Return a mask of sources with lengths greater than
        ``config.max_trail_length``  multiplied by the exposure time.

        Parameters
        ----------
        dia_sources : `pandas.DataFrame`
            Input DIASources to check for trail lengths.
        exposure_time : `float`
            Exposure time from difference image.

        Returns
        -------
        trail_mask : `pandas.DataFrame`
            Boolean mask for DIASources which are greater than the
            cutoff length.
        """
        trail_mask = (dia_sources.loc[:, "trailLength"].values[:]
                      >= (self.config.max_trail_length*exposure_time))

        return trail_mask
