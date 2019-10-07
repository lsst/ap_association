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

from astropy.stats import median_absolute_deviation
import numpy as np
from scipy.optimize import lsq_linear
from scipy.stats import skew

import lsst.geom as geom
import lsst.pex.config as pexConfig

from .diaCalculation import (
    DiaObjectCalculationPluginConfig,
    DiaObjectCalculationPlugin)
from lsst.meas.base.pluginRegistry import register

__all__ = ("MeanDiaPositionConfig", "MeanDiaPosition",
           "WeightedMeanDiaPsFluxConfig", "WeightedMeanDiaPsFlux",
           "PercentileDiaPsFlux", "PercentileDiaPsFluxConfig",
           "SigmaDiaPsFlux", "SigmaDiaPsFluxConfig",
           "Chi2DiaPsFlux", "Chi2DiaPsFluxConfig",
           "MadDiaPsFlux", "MadDiaPsFluxConfig",
           "SkewDiaPsFlux", "SkewDiaPsFluxConfig",
           "MinMaxDiaPsFlux", "MinMaxDiaPsFluxConfig",
           "MaxSlopeDiaPsFlux", "MaxSlopeDiaPsFluxConfig",
           "ErrMeanDiaPsFlux", "ErrMeanDiaPsFluxConfig",
           "LinearFitDiaPsFlux", "LinearFitDiaPsFluxConfig",
           "StetsonJDiaPsFlux", "StetsonJDiaPsFluxConfig",
           "WeightedMeanDiaTotFlux", "WeightedMeanDiaTotFluxConfig",
           "SigmaDiaTotFlux", "SigmaDiaTotFluxConfig")


class MeanDiaPositionConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_meanPosition")
class MeanDiaPosition(DiaObjectCalculationPlugin):
    """Compute the mean position of a DiaObject given a set of DiaSources.
    """

    ConfigClass = MeanDiaPositionConfig
    outputCols = ["ra", "decl"]

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
        error : `BaseException` or `None`
            Error to pass. Kept for consistency with CatologCalculationPlugin.
            Unused.
        """
        diaObject["ra"] = np.nan
        diaObject["decl"] = np.nan


class WeightedMeanDiaPsFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_meanFlux")
class WeightedMeanDiaPsFlux(DiaObjectCalculationPlugin):
    """Compute the weighted mean and mean error on the point source fluxes
    of the DiaSource measured on the difference image.

    Additionally store number of usable data points.
    """

    ConfigClass = WeightedMeanDiaPsFluxConfig
    outputCols = ["PSFluxMean", "PSFluxMeanErr", "PSFFluxNdata"]

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
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaFluxes) > 0:
            tot_weight = np.nansum(1 / filterDiaFluxes["psFluxErr"] ** 2)
            fluxMean = np.nansum(filterDiaFluxes["psFlux"] /
                                 filterDiaFluxes["psFluxErr"] ** 2)
            fluxMean /= tot_weight
            fluxMeanErr = np.sqrt(1 / tot_weight)
            nFluxData = np.sum(np.isfinite(filterDiaFluxes["psFlux"]))
        else:
            fluxMean = np.nan
            fluxMeanErr = np.nan
            nFluxData = 0
        if np.isfinite(fluxMean) and np.isfinite(fluxMeanErr):
            diaObject["{}PSFluxMean".format(filterName)] = fluxMean
            diaObject["{}PSFluxMeanErr".format(filterName)] = fluxMeanErr
            diaObject["{}PSFluxNdata".format(filterName)] = nFluxData
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
        error : `BaseException` or `None`
            Error to pass.
        """
        diaObject["{}PSFluxMean".format(filterName)] = np.nan
        diaObject["{}PSFluxMeanErr".format(filterName)] = np.nan
        diaObject["{}PSFluxNdata".format(filterName)] = 0


class PercentileDiaPsFluxConfig(DiaObjectCalculationPluginConfig):
    percentiles = pexConfig.ListField(
        dtype=int,
        default=[5, 25, 50, 75, 95],
        doc="Percentiles to calculate to compute values for. Should be "
            "integer values."
    )


@register("ap_percentileFlux")
class PercentileDiaPsFlux(DiaObjectCalculationPlugin):
    """Compute percentiles of diaSource fluxes.
    """

    ConfigClass = PercentileDiaPsFluxConfig
    # Output columns are created upon instantiation of the class.
    outputCols = []

    def __init__(self, config, name, metadata, **kwargs):
        DiaObjectCalculationPlugin.__init__(self,
                                            config,
                                            name,
                                            metadata,
                                            **kwargs)
        self.outputCols = ["PSFluxPercentile{:02d}".format(percent)
                           for percent in self.config.percentiles]

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObject,
                  diaSources,
                  filterDiaFluxes,
                  filterName,
                  **kwargs):
        """Compute the percentile fluxes of the point source flux.

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
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaFluxes) > 0:
            pTiles = np.nanpercentile(filterDiaFluxes["psFlux"],
                                      self.config.percentiles)
            for pTile, tilePercent in zip(pTiles, self.config.percentiles):
                diaObject[
                    "{}PSFluxPercentile{:02d}".format(filterName,
                                                      tilePercent)] = pTile
        else:
            self.fail(diaObject, filterName)

    def fail(self, diaObject, filterName, error=None):
        """Set diaObject values to nan.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        filterName : `str`
            Simple name of the filter for the flux being calculated.
        error : `BaseException` or `None`
            Error to pass.
        """
        for pTile in self.config.percentiles:
            diaObject["{}PSFluxPercentile{:02d}".format(filterName,
                                                        pTile)] = np.nan


class SigmaDiaPsFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_sigmaFlux")
class SigmaDiaPsFlux(DiaObjectCalculationPlugin):
    """Compute scatter of diaSource fluxes.
    """

    ConfigClass = SigmaDiaPsFluxConfig
    # Output columns are created upon instantiation of the class.
    outputCols = ["PSFluxSigma"]

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObject,
                  diaSources,
                  filterDiaFluxes,
                  filterName,
                  **kwargs):
        """Compute the sigma fluxes of the point source flux.

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
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaFluxes) > 1:
            diaObject["{}PSFluxSigma".format(filterName)] = np.nanstd(
                filterDiaFluxes["psFlux"])
        else:
            self.fail(diaObject, filterName)

    def fail(self, diaObject, filterName, error=None):
        """Set diaObject values to nan.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        filterName : `str`
            Simple name of the filter for the flux being calculated.
        error : `BaseException` or `None`
            Error to pass.
        """
        diaObject["{}PSFluxSigma".format(filterName)] = np.nan


class Chi2DiaPsFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_chi2Flux")
class Chi2DiaPsFlux(DiaObjectCalculationPlugin):
    """Compute chi2 of diaSource fluxes.
    """

    ConfigClass = Chi2DiaPsFluxConfig

    # Required input Cols
    inputCols = ["PSFluxMean"]
    # Output columns are created upon instantiation of the class.
    outputCols = ["PSFluxChi2"]

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_MOMENTS_CALCULATED

    def calculate(self,
                  diaObject,
                  diaSources,
                  filterDiaFluxes,
                  filterName,
                  **kwargs):
        """Compute the chi2 of the point source fluxes.

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
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaFluxes) > 0:
            diaObject["{}PSFluxChi2".format(filterName)] = np.nansum(
                ((filterDiaFluxes["psFlux"] -
                  diaObject["{}PSFluxMean".format(filterName)]) /
                 filterDiaFluxes["psFluxErr"]) ** 2)
        else:
            self.fail(diaObject, filterName)

    def fail(self, diaObject, filterName, error=None):
        """Set diaObject values to nan.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        filterName : `str`
            Simple name of the filter for the flux being calculated.
        error : `BaseException` or `None`
            Error to pass.
        """
        diaObject["{}PSFluxChi2".format(filterName)] = np.nan


class MadDiaPsFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_madFlux")
class MadDiaPsFlux(DiaObjectCalculationPlugin):
    """Compute median absolute deviation of diaSource fluxes.
    """

    ConfigClass = MadDiaPsFluxConfig

    # Required input Cols
    # Output columns are created upon instantiation of the class.
    outputCols = ["PSFluxMAD"]

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObject,
                  diaSources,
                  filterDiaFluxes,
                  filterName,
                  **kwargs):
        """Compute the median absolute deviation of the point source fluxes.

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
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaFluxes) > 0:
            diaObject["{}PSFluxMAD".format(filterName)] = (
                median_absolute_deviation(filterDiaFluxes["psFlux"],
                                          ignore_nan=True)
            )
        else:
            self.fail(diaObject, filterName)

    def fail(self, diaObject, filterName, error=None):
        """Set diaObject values to nan.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        filterName : `str`
            Simple name of the filter for the flux being calculated.
        error : `BaseException` or `None`
            Error to pass.
        """
        diaObject["{}PSFluxMAD".format(filterName)] = np.nan


class SkewDiaPsFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_skewFlux")
class SkewDiaPsFlux(DiaObjectCalculationPlugin):
    """Compute the skew of diaSource fluxes.
    """

    ConfigClass = SkewDiaPsFluxConfig

    # Required input Cols
    # Output columns are created upon instantiation of the class.
    outputCols = ["PSFluxSkew"]

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObject,
                  diaSources,
                  filterDiaFluxes,
                  filterName,
                  **kwargs):
        """Compute the skew of the point source fluxes.

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
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaFluxes) > 0:
            fluxes = filterDiaFluxes["psFlux"]
            diaObject["{}PSFluxSkew".format(filterName)] = (
                skew(fluxes[~np.isnan(fluxes)])
            )
        else:
            self.fail(diaObject, filterName)

    def fail(self, diaObject, filterName, error=None):
        """Set diaObject values to nan.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        filterName : `str`
            Simple name of the filter for the flux being calculated.
        error : `BaseException` or `None`
            Error to pass.
        """
        diaObject["{}PSFluxSkew".format(filterName)] = np.nan


class MinMaxDiaPsFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_minMaxFlux")
class MinMaxDiaPsFlux(DiaObjectCalculationPlugin):
    """Compute min/max of diaSource fluxes.
    """

    ConfigClass = MinMaxDiaPsFluxConfig

    # Required input Cols
    # Output columns are created upon instantiation of the class.
    outputCols = ["PSFluxMin", "PSFluxMax"]

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObject,
                  diaSources,
                  filterDiaFluxes,
                  filterName,
                  **kwargs):
        """Compute min/max of the point source fluxes.

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
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaFluxes) > 0:
            fluxes = filterDiaFluxes["psFlux"]
            diaObject["{}PSFluxMin".format(filterName)] = np.min(fluxes)
            diaObject["{}PSFluxMax".format(filterName)] = np.max(fluxes)
        else:
            self.fail(diaObject, filterName)

    def fail(self, diaObject, filterName, error=None):
        """Set diaObject values to nan.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        filterName : `str`
            Simple name of the filter for the flux being calculated.
        error : `BaseException` or `None`
            Error to pass.
        """
        diaObject["{}PSFluxMin".format(filterName)] = np.nan
        diaObject["{}PSFluxMax".format(filterName)] = np.nan


class MaxSlopeDiaPsFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_maxSlopeFlux")
class MaxSlopeDiaPsFlux(DiaObjectCalculationPlugin):
    """Compute the maximum ratio time ordered deltaFlux / deltaTime.
    """

    ConfigClass = MinMaxDiaPsFluxConfig

    # Required input Cols
    # Output columns are created upon instantiation of the class.
    outputCols = ["PSFluxMaxSlope"]

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObject,
                  diaSources,
                  filterDiaFluxes,
                  filterName,
                  **kwargs):
        """Compute the maximum ratio time ordered deltaFlux / deltaTime.

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
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaFluxes) > 1:
            tmpDiaSources = filterDiaFluxes[~np.isnan(filterDiaFluxes["psFlux"])]
            fluxes = tmpDiaSources["psFlux"].to_numpy()
            times = tmpDiaSources["midPointTai"].to_numpy()
            diaObject["{}PSFluxMaxSlope".format(filterName)] = np.max(
                (fluxes[1:] - fluxes[:-1]) / (times[1:] - times[:-1]))
        else:
            self.fail(diaObject, filterName)

    def fail(self, diaObject, filterName, error=None):
        """Set diaObject values to nan.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        filterName : `str`
            Simple name of the filter for the flux being calculated.
        error : `BaseException` or `None`
            Error to pass.
        """
        diaObject["{}PSFluxMaxSlope".format(filterName)] = np.nan


class ErrMeanDiaPsFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_meanErrFlux")
class ErrMeanDiaPsFlux(DiaObjectCalculationPlugin):
    """Compute the mean of the dia source errors.
    """

    ConfigClass = ErrMeanDiaPsFluxConfig

    # Required input Cols
    # Output columns are created upon instantiation of the class.
    outputCols = ["PSFluxErrMean"]

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObject,
                  diaSources,
                  filterDiaFluxes,
                  filterName,
                  **kwargs):
        """Compute the mean of the dia source errors.

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
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaFluxes) > 0:
            diaObject["{}PSFluxErrMean".format(filterName)] = np.nanmean(
                filterDiaFluxes["psFluxErr"])
        else:
            self.fail(diaObject, filterName)

    def fail(self, diaObject, filterName, error=None):
        """Set diaObject values to nan.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        filterName : `str`
            Simple name of the filter for the flux being calculated.
        error : `BaseException` or `None`
            Error to pass.
        """
        diaObject["{}PSFluxErrMean".format(filterName)] = np.nan


class LinearFitDiaPsFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_linearFit")
class LinearFitDiaPsFlux(DiaObjectCalculationPlugin):
    """Compute fit a linear model to flux vs time.
    """

    ConfigClass = LinearFitDiaPsFluxConfig

    # Required input Cols
    # Output columns are created upon instantiation of the class.
    outputCols = ["PSFluxLinearSlope", "PSFluxLinearIntercept"]

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObject,
                  diaSources,
                  filterDiaFluxes,
                  filterName,
                  **kwargs):
        """Compute fit a linear model to flux vs time.

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
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaFluxes) > 1:
            tmpDiaSources = filterDiaFluxes[
                ~np.logical_or(np.isnan(filterDiaFluxes["psFlux"]),
                               np.isnan(filterDiaFluxes["psFluxErr"]))]
            fluxes = tmpDiaSources["psFlux"].to_numpy()
            errors = tmpDiaSources["psFluxErr"].to_numpy()
            times = tmpDiaSources["midPointTai"].to_numpy()
            A = np.array([times / errors, 1 / errors]).transpose()
            m, b = lsq_linear(A, fluxes / errors).x
            diaObject["{}PSFluxLinearSlope".format(filterName)] = m
            diaObject["{}PSFluxLinearIntercept".format(filterName)] = b
        else:
            self.fail(diaObject, filterName)

    def fail(self, diaObject, filterName, error=None):
        """Set diaObject values to nan.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        filterName : `str`
            Simple name of the filter for the flux being calculated.
        error : `BaseException` or `None`
            Error to pass.
        """
        diaObject["{}PSFluxLinearSlope".format(filterName)] = np.nan
        diaObject["{}PSFluxLinearIntercept".format(filterName)] = np.nan


class StetsonJDiaPsFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_stetsonJ")
class StetsonJDiaPsFlux(DiaObjectCalculationPlugin):
    """Compute the StetsonJ statistic on the DIA point source fluxes.
    """

    ConfigClass = LinearFitDiaPsFluxConfig

    # Required input Cols
    inputCols = ["PSFluxMean"]
    # Output columns are created upon instantiation of the class.
    outputCols = ["PSFluxStetsonJ"]

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_MOMENTS_CALCULATED

    def calculate(self,
                  diaObject,
                  diaSources,
                  filterDiaFluxes,
                  filterName,
                  **kwargs):
        """Compute the StetsonJ statistic on the DIA point source fluxes.

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
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaFluxes) > 1:
            tmpDiaSources = filterDiaFluxes[
                ~np.logical_or(np.isnan(filterDiaFluxes["psFlux"]),
                               np.isnan(filterDiaFluxes["psFluxErr"]))]
            fluxes = tmpDiaSources["psFlux"].to_numpy()
            errors = tmpDiaSources["psFluxErr"].to_numpy()

            diaObject["{}PSFluxStetsonJ".format(filterName)] = self._stetson_J(
                fluxes,
                errors,
                diaObject["{}PSFluxMean".format(filterName)])
        else:
            self.fail(diaObject, filterName)

    def _stetson_J(self, fluxes, errors, mean=None):
        """Compute the single band stetsonJ statistic.

        Parameters
        ----------
        fluxes : `numpy.ndarray` (N,)
            Calibrated lightcurve flux values.
        errors : `numpy.ndarray` (N,)
            Errors on the calibrated lightcurve fluxes.
        mean : `float`
            Starting mean from previous plugin.

        Returns
        -------
        stetsonJ : `float`
            stetsonJ statistic for the input fluxes and errors.

        References
        ----------
        .. [1] Stetson, P. B., "On the Automatic Determination of Light-Curve
           Parameters for Cepheid Variables", PASP, 108, 851S, 1996
        """
        n_points = len(fluxes)
        flux_mean = self._stetson_mean(fluxes, errors, mean)
        delta_val = (
            np.sqrt(n_points / (n_points - 1)) * (fluxes - flux_mean) / errors)
        p_k = delta_val ** 2 - 1

        return np.mean(np.sign(p_k) * np.sqrt(np.fabs(p_k)))

    def _stetson_mean(self,
                      values,
                      errors,
                      mean=None,
                      alpha=2.,
                      beta=2.,
                      n_iter=20,
                      tol=1e-6):
        """Compute the stetson mean of the fluxes which down-weights outliers.

        Weighted biased on an error weighted difference scaled by a constant
        (1/``a``) and raised to the power beta. Higher betas more harshly
        penalize outliers and ``a`` sets the number of sigma where a weighted
        difference of 1 occurs.

        Parameters
        ----------
        values : `numpy.dnarray`, (N,)
            Input values to compute the mean of.
        errors : `numpy.ndarray`, (N,)
            Errors on the input values.
        mean : `float`
            Starting mean value or None.
        alpha : `float`
            Scalar down-weighting of the fractional difference. lower->more
            clipping
        beta : `float`
            Power law slope of the used to down-weight outliers. higher->more
            clipping
        n_iter : `int`
            Number of iterations of clipping.
        tol : `float`
            Fractional and absolute tolerance goal on the change in the mean
            before exiting early.

        Returns
        -------
        mean : `float`
            Weighted stetson mean result.

        References
        ----------
        .. [1] Stetson, P. B., "On the Automatic Determination of Light-Curve
           Parameters for Cepheid Variables", PASP, 108, 851S, 1996
        """
        n_points = len(values)
        n_factor = np.sqrt(n_points / (n_points - 1))
        inv_var = 1 / errors ** 2

        if mean is None:
            mean = np.average(values, weights=inv_var)
        for iter_idx in range(n_iter):
            chi = np.fabs(n_factor * (values - mean) / errors)
            tmp_mean = np.average(
                values,
                weights=inv_var / (1 + (chi / alpha) ** beta))
            diff = np.fabs(tmp_mean - mean)
            mean = tmp_mean
            if diff / mean < tol and diff < tol:
                break
        return mean

    def fail(self, diaObject, filterName, error=None):
        """Set diaObject values to nan.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        filterName : `str`
            Simple name of the filter for the flux being calculated.
        error : `BaseException` or `None`
            Error to pass.
        """
        diaObject["{}PSFluxStetsonJ".format(filterName)] = np.nan


class WeightedMeanDiaTotFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_meanTotFlux")
class WeightedMeanDiaTotFlux(DiaObjectCalculationPlugin):
    """Compute the weighted mean and mean error on the point source fluxes
    forced photometered at the DiaSource location in the calibrated image.

    Additionally store number of usable data points.
    """

    ConfigClass = WeightedMeanDiaPsFluxConfig
    outputCols = ["TOTFluxMean", "TOTFluxMeanErr"]

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
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaFluxes) > 0:
            tot_weight = np.nansum(1 / filterDiaFluxes["totFluxErr"] ** 2)
            fluxMean = np.nansum(filterDiaFluxes["totFlux"] /
                                 filterDiaFluxes["totFluxErr"] ** 2)
            fluxMean /= tot_weight
            fluxMeanErr = np.sqrt(1 / tot_weight)
        else:
            fluxMean = np.nan
            fluxMeanErr = np.nan
        if np.isfinite(fluxMean) and np.isfinite(fluxMeanErr):
            diaObject["{}TOTFluxMean".format(filterName)] = fluxMean
            diaObject["{}TOTFluxMeanErr".format(filterName)] = fluxMeanErr
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
        error : `BaseException` or `None`
            Error to pass.
        """
        diaObject["{}TOTFluxMean".format(filterName)] = np.nan
        diaObject["{}TOTFluxMeanErr".format(filterName)] = np.nan


class SigmaDiaTotFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_sigmaTotFlux")
class SigmaDiaTotFlux(DiaObjectCalculationPlugin):
    """Compute scatter of diaSource fluxes.
    """

    ConfigClass = SigmaDiaPsFluxConfig
    # Output columns are created upon instantiation of the class.
    outputCols = ["TOTFluxSigma"]

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObject,
                  diaSources,
                  filterDiaFluxes,
                  filterName,
                  **kwargs):
        """Compute the sigma fluxes of the point source flux measured on the
        calibrated image.

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
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaFluxes) > 1:
            diaObject["{}TOTFluxSigma".format(filterName)] = np.nanstd(
                filterDiaFluxes["totFlux"])
        else:
            self.fail(diaObject, filterName)

    def fail(self, diaObject, filterName, error=None):
        """Set diaObject values to nan.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in.
        filterName : `str`
            Simple name of the filter for the flux being calculated.
        error : `BaseException` or `None`
            Error to pass.
        """
        diaObject["{}TOTFluxSigma".format(filterName)] = np.nan
