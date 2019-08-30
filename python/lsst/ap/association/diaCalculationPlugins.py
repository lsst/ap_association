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

    def __init__(self, config, name, metadata):
        DiaObjectCalculationPlugin.__init__(self, config, name, metadata)

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
        coordList = [geom.SpherePoint(src["ra"], src["decl"], geom.degrees)
                     for idx, src in diaSources.iterrows()]
        aveCoord = geom.averageSpherePoint(coordList)
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

    def __init__(self, config, name, metadata):
        DiaObjectCalculationPlugin.__init__(self, config, name, metadata)

    def calculate(self, diaObject, psFluxes, psFluxErrs, filterName, **kwargs):
        """Compute the mean ra/dec position of the diaObject given the
        diaSource locations.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        psFluxes : `numpy.ndarray`
            Point source fluxes on the difference image trimmed of nans and in
            one filter.
        psFluxErrs : `numpy.ndarray`
            Point source flux errors on the difference image trimmed of nans
            and in one filter.
        filterName : `str`
            Simple name of the filter for the flux being calculated.
        """
        if len(psFluxes) > 0:
            psFluxMean = np.average(psFluxes, weights=1 / psFluxErrs ** 2)
            psFluxMeanErr = np.sqrt(1 / np.sum(1 / psFluxErrs ** 2))
        else:
            psFluxMean = np.nan
            psFluxMeanErr = np.nan

        if np.isfinite(psFluxMean) and np.isfinite(psFluxMeanErr):
            diaObject["%sPSFluxMean" % filterName] = psFluxMean
            diaObject["%sPSFluxMeanErr" % filterName] = psFluxMeanErr
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
