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
as defined in the schema of the Apdb both in name and units.
"""

import functools
import warnings

from astropy.stats import median_absolute_deviation
import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear

import lsst.geom as geom
import lsst.pex.config as pexConfig
import lsst.sphgeom as sphgeom

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


def catchWarnings(_func=None, *, warns=[]):
    """Decorator for generically catching numpy warnings.
    """
    def decoratorCatchWarnings(func):
        @functools.wraps(func)
        def wrapperCatchWarnings(*args, **kwargs):
            with warnings.catch_warnings():
                for val in warns:
                    warnings.filterwarnings("ignore", val)
                return func(*args, **kwargs)
        return wrapperCatchWarnings

    if _func is None:
        return decoratorCatchWarnings
    else:
        return decoratorCatchWarnings(_func)


class MeanDiaPositionConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_meanPosition")
class MeanDiaPosition(DiaObjectCalculationPlugin):
    """Compute the mean position of a DiaObject given a set of DiaSources.
    """

    ConfigClass = MeanDiaPositionConfig

    plugType = 'multi'

    outputCols = ["ra", "decl", "radecTai"]
    needsFilter = False

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self, diaObjects, diaSources, **kwargs):
        """Compute the mean ra/dec position of the diaObject given the
        diaSource locations.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            Summary objects to store values in.
        diaSources : `pandas.DataFrame` or `pandas.DataFrameGroupBy`
            Catalog of DiaSources summarized by this DiaObject.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        for outCol in self.outputCols:
            if outCol not in diaObjects.columns:
                diaObjects[outCol] = np.nan

        def _computeMeanPos(df):
            aveCoord = geom.averageSpherePoint(
                list(geom.SpherePoint(src["ra"], src["decl"], geom.degrees)
                     for idx, src in df.iterrows()))
            ra = aveCoord.getRa().asDegrees()
            decl = aveCoord.getDec().asDegrees()
            if np.isnan(ra) or np.isnan(decl):
                radecTai = np.nan
            else:
                radecTai = df["midPointTai"].max()

            return pd.Series({"ra": aveCoord.getRa().asDegrees(),
                              "decl": aveCoord.getDec().asDegrees(),
                              "radecTai": radecTai})

        ans = diaSources.apply(_computeMeanPos)
        diaObjects.loc[:, ["ra", "decl", "radecTai"]] = ans


class HTMIndexDiaPositionConfig(DiaObjectCalculationPluginConfig):

    htmLevel = pexConfig.Field(
        dtype=int,
        doc="Level of the HTM pixelization.",
        default=20,
    )


@register("ap_HTMIndex")
class HTMIndexDiaPosition(DiaObjectCalculationPlugin):
    """Compute the mean position of a DiaObject given a set of DiaSources.
    """
    ConfigClass = HTMIndexDiaPositionConfig

    plugType = 'single'

    inputCols = ["ra", "decl"]
    outputCols = ["pixelId"]
    needsFilter = False

    def __init__(self, config, name, metadata):
        DiaObjectCalculationPlugin.__init__(self, config, name, metadata)
        self.pixelator = sphgeom.HtmPixelization(self.config.htmLevel)

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_MOMENTS_CALCULATED

    def calculate(self, diaObjects, diaObjectId, **kwargs):
        """Compute the mean position of a DiaObject given a set of DiaSources

        Parameters
        ----------
        diaObjects : `pandas.dataFrame`
            Summary objects to store values in and read ra/decl from.
        diaObjectId : `int`
            Id of the diaObject to update.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        sphPoint = geom.SpherePoint(
            diaObjects.at[diaObjectId, "ra"] * geom.degrees,
            diaObjects.at[diaObjectId, "decl"] * geom.degrees)
        diaObjects.at[diaObjectId, "pixelId"] = self.pixelator.index(
            sphPoint.getVector())


class NumDiaSourcesDiaPluginConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_nDiaSources")
class NumDiaSourcesDiaPlugin(DiaObjectCalculationPlugin):
    """Compute the total number of DiaSources associated with this DiaObject.
    """

    ConfigClass = NumDiaSourcesDiaPluginConfig
    outputCols = ["nDiaSources"]
    plugType = "multi"
    needsFilter = False

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self, diaObjects, diaSources, **kwargs):
        """Compute the total number of DiaSources associated with this DiaObject.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in and read ra/decl from.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        diaObjects.loc[:, "nDiaSources"] = diaSources.diaObjectId.count()


class SimpleSourceFlagDiaPluginConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_diaObjectFlag")
class SimpleSourceFlagDiaPlugin(DiaObjectCalculationPlugin):
    """Find if any DiaSource is flagged.

    Set the DiaObject flag if any DiaSource is flagged.
    """

    ConfigClass = NumDiaSourcesDiaPluginConfig
    outputCols = ["flags"]
    plugType = "multi"
    needsFilter = False

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self, diaObjects, diaSources, **kwargs):
        """Find if any DiaSource is flagged.

        Set the DiaObject flag if any DiaSource is flagged.

        Parameters
        ----------
        diaObject : `dict`
            Summary object to store values in and read ra/decl from.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        diaObjects.loc[:, "flags"] = diaSources.flags.any().astype(np.uint64)


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
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    @catchWarnings(warns=["invalid value encountered",
                          "divide by zero"])
    def calculate(self,
                  diaObjects,
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
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        meanName = "{}PSFluxMean".format(filterName)
        errName = "{}PSFluxMeanErr".format(filterName)
        nDataName = "{}PSFluxNdata".format(filterName)
        if meanName not in diaObjects.columns:
            diaObjects[meanName] = np.nan
        if errName not in diaObjects.columns:
            diaObjects[errName] = np.nan
        if nDataName not in diaObjects.columns:
            diaObjects[nDataName] = 0

        def _weightedMean(df):
            tmpDf = df[~np.logical_or(np.isnan(df["psFlux"]),
                                      np.isnan(df["psFluxErr"]))]
            tot_weight = np.nansum(1 / tmpDf["psFluxErr"] ** 2)
            fluxMean = np.nansum(tmpDf["psFlux"]
                                 / tmpDf["psFluxErr"] ** 2)
            fluxMean /= tot_weight
            if tot_weight > 0:
                fluxMeanErr = np.sqrt(1 / tot_weight)
            else:
                fluxMeanErr = np.nan
            nFluxData = len(tmpDf)

            return pd.Series({meanName: fluxMean,
                              errName: fluxMeanErr,
                              nDataName: nFluxData},
                             dtype="object")

        diaObjects.loc[:, [meanName, errName, nDataName]] = \
            filterDiaSources.apply(_weightedMean)


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
    plugType = "multi"
    needsFilter = True

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

    @catchWarnings(warns=["All-NaN slice encountered"])
    def calculate(self,
                  diaObjects,
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
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        pTileNames = []
        for tilePercent in self.config.percentiles:
            pTileName = "{}PSFluxPercentile{:02d}".format(filterName,
                                                          tilePercent)
            pTileNames.append(pTileName)
            if pTileName not in diaObjects.columns:
                diaObjects[pTileName] = np.nan

        def _fluxPercentiles(df):
            pTiles = np.nanpercentile(df["psFlux"], self.config.percentiles)
            return pd.Series(
                dict((tileName, pTile)
                     for tileName, pTile in zip(pTileNames, pTiles)))

        diaObjects.loc[:, pTileNames] = filterDiaSources.apply(_fluxPercentiles)


class SigmaDiaPsFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_sigmaFlux")
class SigmaDiaPsFlux(DiaObjectCalculationPlugin):
    """Compute scatter of diaSource fluxes.
    """

    ConfigClass = SigmaDiaPsFluxConfig
    # Output columns are created upon instantiation of the class.
    outputCols = ["PSFluxSigma"]
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObjects,
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
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        # Set "delta degrees of freedom (ddf)" to 1 to calculate the unbiased
        # estimator of scatter (i.e. 'N - 1' instead of 'N').
        diaObjects.loc[:, "{}PSFluxSigma".format(filterName)] = \
            filterDiaSources.psFlux.std()


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
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_MOMENTS_CALCULATED

    @catchWarnings(warns=["All-NaN slice encountered"])
    def calculate(self,
                  diaObjects,
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
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        meanName = "{}PSFluxMean".format(filterName)

        def _chi2(df):
            delta = (df["psFlux"]
                     - diaObjects.at[df.diaObjectId.iat[0], meanName])
            return np.nansum((delta / df["psFluxErr"]) ** 2)

        diaObjects.loc[:, "{}PSFluxChi2".format(filterName)] = \
            filterDiaSources.apply(_chi2)


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
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    @catchWarnings(warns=["All-NaN slice encountered"])
    def calculate(self,
                  diaObjects,
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
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        diaObjects.loc[:, "{}PSFluxMAD".format(filterName)] = \
            filterDiaSources.psFlux.apply(median_absolute_deviation,
                                          ignore_nan=True)


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
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObjects,
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
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        diaObjects.loc[:, "{}PSFluxSkew".format(filterName)] = \
            filterDiaSources.psFlux.skew()


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
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObjects,
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
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        minName = "{}PSFluxMin".format(filterName)
        if minName not in diaObjects.columns:
            diaObjects[minName] = np.nan
        maxName = "{}PSFluxMax".format(filterName)
        if maxName not in diaObjects.columns:
            diaObjects[maxName] = np.nan

        diaObjects.loc[:, minName] = filterDiaSources.psFlux.min()
        diaObjects.loc[:, maxName] = filterDiaSources.psFlux.max()


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
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObjects,
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
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """

        def _maxSlope(df):
            tmpDf = df[~np.logical_or(np.isnan(df["psFlux"]),
                                      np.isnan(df["midPointTai"]))]
            if len(tmpDf) < 2:
                return np.nan
            times = tmpDf["midPointTai"].to_numpy()
            timeArgs = times.argsort()
            times = times[timeArgs]
            fluxes = tmpDf["psFlux"].to_numpy()[timeArgs]
            return (np.diff(fluxes) / np.diff(times)).max()

        diaObjects.loc[:, "{}PSFluxMaxSlope".format(filterName)] = \
            filterDiaSources.apply(_maxSlope)


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
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObjects,
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
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        diaObjects.loc[:, "{}PSFluxErrMean".format(filterName)] = \
            filterDiaSources.psFluxErr.mean()


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
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObjects,
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
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """

        mName = "{}PSFluxLinearSlope".format(filterName)
        if mName not in diaObjects.columns:
            diaObjects[mName] = np.nan
        bName = "{}PSFluxLinearIntercept".format(filterName)
        if bName not in diaObjects.columns:
            diaObjects[bName] = np.nan

        def _linearFit(df):
            tmpDf = df[~np.logical_or(
                np.isnan(df["psFlux"]),
                np.logical_or(np.isnan(df["psFluxErr"]),
                              np.isnan(df["midPointTai"])))]
            if len(tmpDf) < 2:
                return pd.Series({mName: np.nan, bName: np.nan})
            fluxes = tmpDf["psFlux"].to_numpy()
            errors = tmpDf["psFluxErr"].to_numpy()
            times = tmpDf["midPointTai"].to_numpy()
            A = np.array([times / errors, 1 / errors]).transpose()
            m, b = lsq_linear(A, fluxes / errors).x
            return pd.Series({mName: m, bName: b})

        diaObjects.loc[:, [mName, bName]] = filterDiaSources.apply(_linearFit)


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
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_MOMENTS_CALCULATED

    def calculate(self,
                  diaObjects,
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
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        meanName = "{}PSFluxMean".format(filterName)

        def _stetsonJ(df):
            tmpDf = df[~np.logical_or(np.isnan(df["psFlux"]),
                                      np.isnan(df["psFluxErr"]))]
            if len(tmpDf) < 2:
                return np.nan
            fluxes = tmpDf["psFlux"].to_numpy()
            errors = tmpDf["psFluxErr"].to_numpy()

            return self._stetson_J(
                fluxes,
                errors,
                diaObjects.at[tmpDf.diaObjectId.iat[0], meanName])

        diaObjects.loc[:, "{}PSFluxStetsonJ".format(filterName)] = \
            filterDiaSources.apply(_stetsonJ)

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
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    @catchWarnings(warns=["invalid value encountered",
                          "divide by zero"])
    def calculate(self,
                  diaObjects,
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
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        totMeanName = "{}TOTFluxMean".format(filterName)
        if totMeanName not in diaObjects.columns:
            diaObjects[totMeanName] = np.nan
        totErrName = "{}TOTFluxMeanErr".format(filterName)
        if totErrName not in diaObjects.columns:
            diaObjects[totErrName] = np.nan

        def _meanFlux(df):
            tmpDf = df[~np.logical_or(np.isnan(df["totFlux"]),
                                      np.isnan(df["totFluxErr"]))]
            tot_weight = np.nansum(1 / tmpDf["totFluxErr"] ** 2)
            fluxMean = np.nansum(tmpDf["totFlux"]
                                 / tmpDf["totFluxErr"] ** 2)
            fluxMean /= tot_weight
            fluxMeanErr = np.sqrt(1 / tot_weight)

            return pd.Series({totMeanName: fluxMean,
                              totErrName: fluxMeanErr})

        diaObjects.loc[:, [totMeanName, totErrName]] = \
            filterDiaSources.apply(_meanFlux)


class SigmaDiaTotFluxConfig(DiaObjectCalculationPluginConfig):
    pass


@register("ap_sigmaTotFlux")
class SigmaDiaTotFlux(DiaObjectCalculationPlugin):
    """Compute scatter of diaSource fluxes.
    """

    ConfigClass = SigmaDiaPsFluxConfig
    # Output columns are created upon instantiation of the class.
    outputCols = ["TOTFluxSigma"]
    plugType = "multi"
    needsFilter = True

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObjects,
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
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        # Set "delta degrees of freedom (ddf)" to 1 to calculate the unbiased
        # estimator of scatter (i.e. 'N - 1' instead of 'N').
        diaObjects.loc[:, "{}TOTFluxSigma".format(filterName)] = \
            filterDiaSources.totFlux.std()
