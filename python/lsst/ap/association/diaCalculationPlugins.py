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

"""Plugins for use in DiaSource summary statistics.

Output columns must be
as defined in the schema of the Ppdb both in name and units.
"""

from astropy.stats import median_absolute_deviation
import numpy as np
from scipy.optimize import lsq_linear
from scipy.stats import skew

import lsst.geom as geom
from lsst.meas.algorithms.indexerRegistry import IndexerRegistry
import lsst.pex.config as pexConfig

from .diaCalculation import (
    DiaObjectCalculationPluginConfig,
    DiaObjectCalculationPlugin)
from lsst.meas.base.pluginRegistry import register

__all__ = ("MeanDiaPositionConfig", "MeanDiaPosition",
           "HTMIndexDiaPosition", "HTMIndexDiaPositionConfig",
           "NumDiaSourcesDiaPlugin", "NumDiaSourcesDiaPluginConfig",
           "SimpleSourceFlagDiaPlugin", "SimpleSourceFlagDiaPluginConfig",
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
    outputCols = ["ra", "decl", "radecTai"]

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
            self.fail(diaObject, self.outputCols)
        else:
            diaObject["ra"] = aveCoord.getRa().asDegrees()
            diaObject["decl"] = aveCoord.getDec().asDegrees()
            diaObject["radecTai"] = np.max(diaSources["midPointTai"])


class HTMIndexDiaPositionConfig(DiaObjectCalculationPluginConfig):

    indexer = IndexerRegistry.makeField(
        doc='Select the spatial indexer to use within the database. For this '
            'plugin we enforce HTM pixelization. The configuration as is '
            'allows for different resolutions to be used.',
        default="HTM",
    )

    def validate(self):
        self.indexer == "HTM"


@register("ap_HTMIndex")
class HTMIndexDiaPosition(DiaObjectCalculationPlugin):
    """Compute the mean position of a DiaObject given a set of DiaSources.
    """

    ConfigClass = HTMIndexDiaPositionConfig
    inputCols = ["ra", "decl"]
    outputCols = ["pixelId"]

    def __init__(self, config, name, metadata):
        DiaObjectCalculationPlugin.__init__(self, config, name, metadata)
        self.indexer = IndexerRegistry[self.config.indexer.name](
            self.config.indexer.active)

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_MOMENTS_CALCULATED

    def calculate(self, diaObject, **kwargs):
        """Compute the mean position of a DiaObject given a set of DiaSources

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in and read ra/decl from.
        """
        diaObject["pixelId"] = self.indexer.indexPoints([diaObject["ra"]],
                                                        [diaObject["decl"]])[0]


class NumDiaSourcesDiaPluginConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_nDiaSources")
class NumDiaSourcesDiaPlugin(DiaObjectCalculationPlugin):
    """Compute the total number of DiaSources associated with this DiaObject.
    """

    ConfigClass = NumDiaSourcesDiaPluginConfig
    outputCols = ["nDiaSources"]

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self, diaObject, diaSources, **kwargs):
        """Compute the total number of DiaSources associated with this DiaObject.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in and read ra/decl from.
        """
        diaObject["nDiaSources"] = len(diaSources)


class SimpleSourceFlagDiaPluginConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_diaObjectFlag")
class SimpleSourceFlagDiaPlugin(DiaObjectCalculationPlugin):
    """Find if any DiaSource is flagged.

    Set the DiaObject flag if any DiaSource is flagged.
    """

    ConfigClass = NumDiaSourcesDiaPluginConfig
    outputCols = ["flags"]

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self, diaObject, diaSources, **kwargs):
        """Find if any DiaSource is flagged.

        Set the DiaObject flag if any DiaSource is flagged.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in and read ra/decl from.
        """
        if np.any(diaSources["flags"] > 0):
            diaObject["flags"] = 1
        else:
            diaObject["flags"] = 0


class WeightedMeanDiaPsFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_meanFlux")
class WeightedMeanDiaPsFlux(DiaObjectCalculationPlugin):
    """Compute the weighted mean and mean error on the point source fluxes
    of the DiaSource measured on the difference image.

    Additionally store number of usable data points.
    """

    ConfigClass = WeightedMeanDiaPsFluxConfig
    outputCols = ["PSFluxMean", "PSFluxMeanErr", "PSFluxNdata"]

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObject,
                  diaSources,
                  filterDiaSources,
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
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``filterName``.
        filterName : `str`
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaSources) > 0:
            tot_weight = np.nansum(1 / filterDiaSources["psFluxErr"] ** 2)
            fluxMean = np.nansum(filterDiaSources["psFlux"] /
                                 filterDiaSources["psFluxErr"] ** 2)
            fluxMean /= tot_weight
            fluxMeanErr = np.sqrt(1 / tot_weight)
            nFluxData = np.sum(np.isfinite(filterDiaSources["psFlux"]))
        else:
            fluxMean = np.nan
            fluxMeanErr = np.nan
            nFluxData = 0
        if np.isfinite(fluxMean) and np.isfinite(fluxMeanErr):
            diaObject["{}PSFluxMean".format(filterName)] = fluxMean
            diaObject["{}PSFluxMeanErr".format(filterName)] = fluxMeanErr
            diaObject["{}PSFluxNdata".format(filterName)] = nFluxData
        else:
            self.fail(diaObject,
                      ["{}{}".format(filterName, colName)
                       for colName in self.outputCols])

    def fail(self, diaObject, columns, error=None):
        """Set diaObject position values to nan.

        Since we set an explicit value instead of nan for all, we override
        the fail method.
        """
        for colName in columns:
            if colName.endswith("Ndata"):
                diaObject[colName] = 0
            else:
                diaObject[colName] = np.nan


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
                  filterDiaSources,
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
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``filterName``.
        filterName : `str`
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaSources) > 0:
            pTiles = np.nanpercentile(filterDiaSources["psFlux"],
                                      self.config.percentiles)
            for pTile, tilePercent in zip(pTiles, self.config.percentiles):
                diaObject[
                    "{}PSFluxPercentile{:02d}".format(filterName,
                                                      tilePercent)] = pTile
        else:
            self.fail(diaObject,
                      ["{}{}".format(filterName, colName)
                       for colName in self.outputCols])


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
                  filterDiaSources,
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
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``filterName``.
        filterName : `str`
            Simple, string name of the filter for the flux being calculated.
        """
        # Set "delta degrees of freedom (ddf)" to 1 to calculate the unbiased
        # estimator of scatter (i.e. 'N - 1' instead of 'N').
        if len(filterDiaSources) > 1:
            diaObject["{}PSFluxSigma".format(filterName)] = np.nanstd(
                filterDiaSources["psFlux"],
                ddof=1)
        else:
            self.fail(diaObject,
                      ["{}{}".format(filterName, colName)
                       for colName in self.outputCols])


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
                  filterDiaSources,
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
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``filterName``.
        filterName : `str`
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaSources) > 0:
            delta = (filterDiaSources["psFlux"] -
                     diaObject["{}PSFluxMean".format(filterName)])
            diaObject["{}PSFluxChi2".format(filterName)] = np.nansum(
                (delta / filterDiaSources["psFluxErr"]) ** 2)
        else:
            self.fail(diaObject,
                      ["{}{}".format(filterName, colName)
                       for colName in self.outputCols])


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
                  filterDiaSources,
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
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``filterName``.
        filterName : `str`
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaSources) > 0:
            diaObject["{}PSFluxMAD".format(filterName)] = (
                median_absolute_deviation(filterDiaSources["psFlux"],
                                          ignore_nan=True)
            )
        else:
            self.fail(diaObject,
                      ["{}{}".format(filterName, colName)
                       for colName in self.outputCols])


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
                  filterDiaSources,
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
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``filterName``.
        filterName : `str`
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaSources) > 0:
            fluxes = filterDiaSources["psFlux"]
            diaObject["{}PSFluxSkew".format(filterName)] = (
                skew(fluxes[~np.isnan(fluxes)])
            )
        else:
            self.fail(diaObject,
                      ["{}{}".format(filterName, colName)
                       for colName in self.outputCols])


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
                  filterDiaSources,
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
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``filterName``.
        filterName : `str`
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaSources) > 0:
            fluxes = filterDiaSources["psFlux"]
            diaObject["{}PSFluxMin".format(filterName)] = np.min(fluxes)
            diaObject["{}PSFluxMax".format(filterName)] = np.max(fluxes)
        else:
            self.fail(diaObject,
                      ["{}{}".format(filterName, colName)
                       for colName in self.outputCols])


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
                  filterDiaSources,
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
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``filterName``.
        filterName : `str`
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaSources) > 1:
            tmpDiaSources = filterDiaSources[~np.isnan(filterDiaSources["psFlux"])]
            fluxes = tmpDiaSources["psFlux"].to_numpy()
            times = tmpDiaSources["midPointTai"].to_numpy()
            diaObject["{}PSFluxMaxSlope".format(filterName)] = np.max(
                (fluxes[1:] - fluxes[:-1]) / (times[1:] - times[:-1]))
        else:
            self.fail(diaObject,
                      ["{}{}".format(filterName, colName)
                       for colName in self.outputCols])


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
                  filterDiaSources,
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
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``filterName``.
        filterName : `str`
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaSources) > 0:
            diaObject["{}PSFluxErrMean".format(filterName)] = np.nanmean(
                filterDiaSources["psFluxErr"])
        else:
            self.fail(diaObject,
                      ["{}{}".format(filterName, colName)
                       for colName in self.outputCols])


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
                  filterDiaSources,
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
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``filterName``.
        filterName : `str`
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaSources) > 1:
            tmpDiaSources = filterDiaSources[
                ~np.logical_or(np.isnan(filterDiaSources["psFlux"]),
                               np.isnan(filterDiaSources["psFluxErr"]))]
            fluxes = tmpDiaSources["psFlux"].to_numpy()
            errors = tmpDiaSources["psFluxErr"].to_numpy()
            times = tmpDiaSources["midPointTai"].to_numpy()
            A = np.array([times / errors, 1 / errors]).transpose()
            m, b = lsq_linear(A, fluxes / errors).x
            diaObject["{}PSFluxLinearSlope".format(filterName)] = m
            diaObject["{}PSFluxLinearIntercept".format(filterName)] = b
        else:
            self.fail(diaObject,
                      ["{}{}".format(filterName, colName)
                       for colName in self.outputCols])


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
                  filterDiaSources,
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
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``filterName``.
        filterName : `str`
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaSources) > 1:
            tmpDiaSources = filterDiaSources[
                ~np.logical_or(np.isnan(filterDiaSources["psFlux"]),
                               np.isnan(filterDiaSources["psFluxErr"]))]
            fluxes = tmpDiaSources["psFlux"].to_numpy()
            errors = tmpDiaSources["psFluxErr"].to_numpy()

            diaObject["{}PSFluxStetsonJ".format(filterName)] = self._stetson_J(
                fluxes,
                errors,
                diaObject["{}PSFluxMean".format(filterName)])
        else:
            self.fail(diaObject,
                      ["{}{}".format(filterName, colName)
                       for colName in self.outputCols])

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
            clipping. (Default value is 2.)
        beta : `float`
            Power law slope of the used to down-weight outliers. higher->more
            clipping. (Default value is 2.)
        n_iter : `int`
            Number of iterations of clipping.
        tol : `float`
            Fractional and absolute tolerance goal on the change in the mean
            before exiting early. (Default value is 1e-6)

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
                  filterDiaSources,
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
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``filterName``.
        filterName : `str`
            Simple, string name of the filter for the flux being calculated.
        """
        if len(filterDiaSources) > 0:
            tot_weight = np.nansum(1 / filterDiaSources["totFluxErr"] ** 2)
            fluxMean = np.nansum(filterDiaSources["totFlux"] /
                                 filterDiaSources["totFluxErr"] ** 2)
            fluxMean /= tot_weight
            fluxMeanErr = np.sqrt(1 / tot_weight)
        else:
            fluxMean = np.nan
            fluxMeanErr = np.nan
        if np.isfinite(fluxMean) and np.isfinite(fluxMeanErr):
            diaObject["{}TOTFluxMean".format(filterName)] = fluxMean
            diaObject["{}TOTFluxMeanErr".format(filterName)] = fluxMeanErr
        else:
            self.fail(diaObject,
                      ["{}{}".format(filterName, colName)
                       for colName in self.outputCols])


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
                  filterDiaSources,
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
        filterDiaSources : `pandas.DataFrame`
            DataFrame representing diaSources associated with this
            diaObject that are observed in the band pass ``filterName``.
        filterName : `str`
            Simple, string name of the filter for the flux being calculated.
        """
        # Set "delta degrees of freedom (ddf)" to 1 to calculate the unbiased
        # estimator of scatter (i.e. 'N - 1' instead of 'N').
        if len(filterDiaSources) > 1:
            diaObject["{}TOTFluxSigma".format(filterName)] = np.nanstd(
                filterDiaSources["totFlux"],
                ddof=1)
        else:
            self.fail(diaObject,
                      ["{}{}".format(filterName, colName)
                       for colName in self.outputCols])
