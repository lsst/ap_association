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

import numpy as np

import lsst.geom as geom

from .diaCalculation import (
    DiaObjectCalculationPluginConfig,
    DiaObjectCalculationPlugin)
from lsst.meas.base.pluginRegistry import register

__all__ = ("MeanDiaPositionConfig", "MeanDiaPosition",
           "WeightedMeanDiaPsFluxConfig", "WeightedMeanDiaPsFlux")


class MeanDiaPositionConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_meanPosition")
class MeanDiaPosition(DiaObjectCalculationPlugin):
    """Compute the mean position of a DiaObject given a set of DiaSources.
    """

    ConfigClass = MeanDiaPositionConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self, diaObject, diaSources, **kwargs):
        """Compute the mean ra/dec position of the diaObject given the
        diaSource locations.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        diaSources : `pandas.DataFrame`
            Catalog of DiaSources summarized by this DiaObject.
        """
        aveCoord = geom.averageSpherePoint(
            list(geom.SpherePoint(src["ra"], src["decl"], geom.degrees)
                 for idx, src in diaSources.iterrows()))
        if not (np.isfinite(aveCoord.getRa().asDegrees()) and
                np.isfinite(aveCoord.getDec().asDegrees())):
            self.fail(diaObject)
        else:
            diaObject["ra"] = aveCoord.getRa().asDegrees()
            diaObject["decl"] = aveCoord.getDec().asDegrees()

    def fail(self, diaObject, error=None):
        """Set diaObject position values to nan.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        error : `BaseException`
            Error to pass.
        """
        diaObject["ra"] = np.nan
        diaObject["decl"] = np.nan


class WeightedMeanDiaPsFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_meanFlux")
class WeightedMeanDiaPsFlux(DiaObjectCalculationPlugin):
    """Compute the weighted mean and mean error on the point source fluxes
    of the DiaSource measured on the difference image.
    """

    ConfigClass = WeightedMeanDiaPsFluxConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObject,
                  diaSources,
                  filterDiaFluxes,
                  filterName,
                  **kwargs):
        """Compute the weighted mean and mean error of the point source flux.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        diaSources : `pandas.DataFrame`
            DataFrame representing all diaSources associated with this
            diaObject.
        filterDiaFluxes : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``filterName``.
        filterName : `str`
            Simple name of the filter for the flux being calculated.
        """
        if len(filterDiaFluxes) > 0:
            tot_weight = np.nansum(1 / filterDiaFluxes["psFluxErr"] ** 2)
            fluxMean = np.nansum(filterDiaFluxes["psFlux"] /
                                 filterDiaFluxes["psFluxErr"] ** 2)
            fluxMean /= tot_weight
            fluxMeanErr = np.sqrt(1 / tot_weight)
        else:
            fluxMean = np.nan
            fluxMeanErr = np.nan

        if np.isfinite(fluxMean) and np.isfinite(fluxMeanErr):
            diaObject["%sPSFluxMean" % filterName] = fluxMean
            diaObject["%sPSFluxMeanErr" % filterName] = fluxMeanErr
        else:
            self.fail(diaObject, filterName)

    def fail(self, diaObject, filterName, error=None):
        """Set diaObject position values to nan.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        filterName : `str`
            Simple name of the filter for the flux being calculated.
        error : `BaseException`
            Error to pass.
        """
        diaObject["%sPSFluxMean" % filterName] = np.nan
        diaObject["%sPSFluxMeanErr" % filterName] = np.nan
