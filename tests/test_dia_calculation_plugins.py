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
import pandas as pd
from scipy.stats import skew
import unittest

from lsst.ap.association import (
    MeanDiaPosition, MeanDiaPositionConfig,
    WeightedMeanDiaPsFlux, WeightedMeanDiaPsFluxConfig,
    PercentileDiaPsFlux, PercentileDiaPsFluxConfig,
    SigmaDiaPsFlux, SigmaDiaPsFluxConfig,
    Chi2DiaPsFlux, Chi2DiaPsFluxConfig,
    MadDiaPsFlux, MadDiaPsFluxConfig,
    SkewDiaPsFlux, SkewDiaPsFluxConfig,
    MinMaxDiaPsFlux, MinMaxDiaPsFluxConfig,
    MaxSlopeDiaPsFlux, MaxSlopeDiaPsFluxConfig,
    ErrMeanDiaPsFlux, ErrMeanDiaPsFluxConfig,
    LinearFitDiaPsFlux, LinearFitDiaPsFluxConfig,
    StetsonJDiaPsFlux, StetsonJDiaPsFlux)
import lsst.utils.tests


class TestMeanPosition(unittest.TestCase):

    def testCalculate(self):
        """Test mean position calculation.
        """
        n_sources = 10

        # Test expected means.
        diaObject = dict()
        diaSources = pd.DataFrame(data={"ra": np.linspace(-1, 1, n_sources),
                                        "decl": np.zeros(n_sources)})
        plug = MeanDiaPosition(MeanDiaPositionConfig(),
                               "ap_meanPosition",
                               None)
        plug.calculate(diaObject, diaSources)

        self.assertAlmostEqual(diaObject["ra"], 0.0)
        self.assertAlmostEqual(diaObject["decl"], 0.0)

        diaObject = dict()
        diaSources = pd.DataFrame(data={"ra": np.zeros(n_sources),
                                        "decl": np.linspace(-1, 1, n_sources)})
        plug.calculate(diaObject, diaSources)

        self.assertAlmostEqual(diaObject["ra"], 0.0)
        self.assertAlmostEqual(diaObject["decl"], 0.0)

        # Test failure modes.
        diaObject = dict()
        diaSources = pd.DataFrame(data={"ra": np.full(n_sources, np.nan),
                                        "decl": np.zeros(n_sources)})
        plug.calculate(diaObject, diaSources)

        self.assertTrue(np.isnan(diaObject["ra"]))
        self.assertTrue(np.isnan(diaObject["decl"]))

        diaObject = dict()
        diaSources = pd.DataFrame(data={"ra": np.zeros(n_sources),
                                        "decl": np.full(n_sources, np.nan)})
        plug.calculate(diaObject, diaSources)

        self.assertTrue(np.isnan(diaObject["ra"]))
        self.assertTrue(np.isnan(diaObject["decl"]))


class TestWeightedMeanDiaPsFlux(unittest.TestCase):

    def testCalculate(self):
        """Test mean value calculation.
        """
        n_sources = 10
        diaObject = dict()
        diaSources = pd.DataFrame(data={"psFlux": np.linspace(-1, 1, n_sources),
                                        "psFluxErr": np.ones(n_sources)})

        plug = WeightedMeanDiaPsFlux(WeightedMeanDiaPsFluxConfig(),
                                     "ap_meanFlux",
                                     None)
        plug.calculate(diaObject, diaSources, diaSources, "u")

        self.assertAlmostEqual(diaObject["uPSFluxMean"], 0.0)
        self.assertAlmostEqual(diaObject["uPSFluxMeanErr"], np.sqrt(1 / n_sources))
        self.assertEqual(diaObject["uPSFluxNdata"], n_sources)

        diaObject = dict()
        plug.calculate(diaObject, [], [], "g")

        self.assertTrue(np.isnan(diaObject["gPSFluxMean"]))
        self.assertTrue(np.isnan(diaObject["gPSFluxMeanErr"]))
        self.assertEqual(diaObject["gPSFluxNdata"], 0)

        diaObject = dict()
        diaSources.loc[4, "psFlux"] = np.nan
        plug.calculate(diaObject, diaSources, diaSources, "r")

        self.assertTrue(~np.isnan(diaObject["rPSFluxMean"]))
        self.assertTrue(~np.isnan(diaObject["rPSFluxMeanErr"]))
        self.assertEqual(diaObject["rPSFluxNdata"], n_sources - 1)


class TestPercentileDiaPsFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux percentile calculation.
        """
        n_sources = 10
        diaObject = dict()
        fluxes = np.linspace(-1, 1, n_sources)
        diaSources = pd.DataFrame(data={"psFlux": fluxes,
                                        "psFluxErr": np.ones(n_sources)})

        plug = PercentileDiaPsFlux(PercentileDiaPsFluxConfig(),
                                   "ap_percentileFlux",
                                   None)
        plug.calculate(diaObject, diaSources, diaSources, "u")
        for pTile, testVal in zip(plug.config.percentiles,
                                  np.nanpercentile(
                                      fluxes,
                                      plug.config.percentiles)):
            self.assertAlmostEqual(
                diaObject["uPSFluxPercentile{:02d}".format(pTile)],
                testVal)

        diaObject = dict()
        plug.calculate(diaObject, [], [], "g")
        for pTile in plug.config.percentiles:
            self.assertTrue(np.isnan(
                diaObject["gPSFluxPercentile{:02d}".format(pTile)]))

        diaObject = dict()
        diaSources.loc[4, "psFlux"] = np.nan
        fluxes[4] = np.nan
        plug.calculate(diaObject, diaSources, diaSources, "r")
        for pTile, testVal in zip(plug.config.percentiles,
                                  np.nanpercentile(
                                      fluxes,
                                      plug.config.percentiles)):
            self.assertAlmostEqual(
                diaObject["rPSFluxPercentile{:02d}".format(pTile)],
                testVal)


class TestSigmaDiaPsFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux scatter calculation.
        """
        n_sources = 10
        diaObject = dict()
        fluxes = np.linspace(-1, 1, n_sources)
        diaSources = pd.DataFrame(data={"psFlux": fluxes,
                                        "psFluxErr": np.ones(n_sources)})

        plug = SigmaDiaPsFlux(SigmaDiaPsFluxConfig(),
                              "ap_sigmaFlux",
                              None)
        plug.calculate(diaObject, diaSources, diaSources, "u")
        self.assertAlmostEqual(diaObject["uPSFluxSigma"],
                               np.nanstd(fluxes))

        # test no inputs
        diaObject = dict()
        plug.calculate(diaObject, [], [], "g")
        self.assertTrue(np.isnan(diaObject["gPSFluxSigma"]))

        # test one input
        diaObject = dict()
        diaSources = pd.DataFrame(data={"psFlux": [fluxes[0]],
                                        "psFluxErr": [np.ones(n_sources)[0]]})
        plug.calculate(diaObject, diaSources, diaSources, "g")
        self.assertTrue(np.isnan(diaObject["gPSFluxSigma"]))

        diaObject = dict()
        fluxes[4] = np.nan
        diaSources = pd.DataFrame(data={"psFlux": fluxes,
                                        "psFluxErr": np.ones(n_sources)})
        plug.calculate(diaObject, diaSources, diaSources, "r")
        self.assertAlmostEqual(diaObject["rPSFluxSigma"],
                               np.nanstd(fluxes))


class TestChi2DiaPsFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux chi2 calculation.
        """
        n_sources = 10
        diaObject = dict()
        diaObject["uPSFluxMean"] = 0.0
        fluxes = np.linspace(-1, 1, n_sources)
        diaSources = pd.DataFrame(data={"psFlux": fluxes,
                                        "psFluxErr": np.ones(n_sources)})

        plug = Chi2DiaPsFlux(Chi2DiaPsFluxConfig(),
                             "ap_chi2Flux",
                             None)
        plug.calculate(diaObject, diaSources, diaSources, "u")
        self.assertAlmostEqual(
            diaObject["uPSFluxChi2"],
            np.nansum(((diaSources["psFlux"] -
                        np.nanmean(diaSources["psFlux"])) /
                       diaSources["psFluxErr"]) ** 2))

        # test no inputs
        diaObject = dict()
        plug.calculate(diaObject, [], [], "g")
        self.assertTrue(np.isnan(diaObject["gPSFluxChi2"]))

        diaObject = dict()
        diaSources.loc[4, "psFlux"] = np.nan
        fluxes[4] = np.nan
        diaObject["rPSFluxMean"] = np.nanmean(fluxes)
        plug.calculate(diaObject, diaSources, diaSources, "r")
        self.assertAlmostEqual(
            diaObject["rPSFluxChi2"],
            np.nansum(((diaSources["psFlux"] -
                        np.nanmean(diaSources["psFlux"])) /
                       diaSources["psFluxErr"]) ** 2))


class TestMadDiaPsFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux median absolute deviation calculation.
        """
        n_sources = 10
        diaObject = dict()
        fluxes = np.linspace(-1, 1, n_sources)
        diaSources = pd.DataFrame(data={"psFlux": fluxes,
                                        "psFluxErr": np.ones(n_sources)})

        plug = MadDiaPsFlux(MadDiaPsFluxConfig(),
                            "ap_madFlux",
                            None)
        plug.calculate(diaObject, diaSources, diaSources, "u")
        self.assertAlmostEqual(diaObject["uPSFluxMAD"],
                               median_absolute_deviation(fluxes,
                                                         ignore_nan=True))

        # test no inputs
        diaObject = dict()
        plug.calculate(diaObject, [], [], "g")
        self.assertTrue(np.isnan(diaObject["gPSFluxMAD"]))

        diaObject = dict()
        fluxes[4] = np.nan
        diaSources = pd.DataFrame(data={"psFlux": fluxes,
                                        "psFluxErr": np.ones(n_sources)})
        plug.calculate(diaObject, diaSources, diaSources, "r")
        self.assertAlmostEqual(diaObject["rPSFluxMAD"],
                               median_absolute_deviation(fluxes,
                                                         ignore_nan=True))


class TestSkewDiaPsFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux skew calculation.
        """
        n_sources = 10
        diaObject = dict()
        fluxes = np.linspace(-1, 1, n_sources)
        diaSources = pd.DataFrame(data={"psFlux": fluxes,
                                        "psFluxErr": np.ones(n_sources)})

        plug = SkewDiaPsFlux(SkewDiaPsFluxConfig(),
                             "ap_skewFlux",
                             None)
        plug.calculate(diaObject, diaSources, diaSources, "u")
        self.assertAlmostEqual(diaObject["uPSFluxSkew"], skew(fluxes))

        # test no inputs
        diaObject = dict()
        plug.calculate(diaObject, [], [], "g")
        self.assertTrue(np.isnan(diaObject["gPSFluxSkew"]))

        diaObject = dict()
        fluxes[4] = np.nan
        diaSources = pd.DataFrame(data={"psFlux": fluxes,
                                        "psFluxErr": np.ones(n_sources)})
        plug.calculate(diaObject, diaSources, diaSources, "r")
        cutFluxes = fluxes[~np.isnan(fluxes)]
        self.assertAlmostEqual(diaObject["rPSFluxSkew"], skew(cutFluxes))


class TestMinMaxDiaPsFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux min/max calculation.
        """
        n_sources = 10
        diaObject = dict()
        fluxes = np.linspace(-1, 1, n_sources)
        diaSources = pd.DataFrame(data={"psFlux": fluxes,
                                        "psFluxErr": np.ones(n_sources)})

        plug = MinMaxDiaPsFlux(MinMaxDiaPsFluxConfig(),
                               "ap_minMaxFlux",
                               None)
        plug.calculate(diaObject, diaSources, diaSources, "u")
        self.assertEqual(diaObject["uPSFluxMin"], -1)
        self.assertEqual(diaObject["uPSFluxMax"], 1)

        # test no inputs
        diaObject = dict()
        plug.calculate(diaObject, [], [], "g")
        self.assertTrue(np.isnan(diaObject["gPSFluxMin"]))
        self.assertTrue(np.isnan(diaObject["gPSFluxMax"]))

        diaObject = dict()
        fluxes[4] = np.nan
        diaSources = pd.DataFrame(data={"psFlux": fluxes,
                                        "psFluxErr": np.ones(n_sources)})
        plug.calculate(diaObject, diaSources, diaSources, "r")
        self.assertEqual(diaObject["rPSFluxMin"], -1)
        self.assertEqual(diaObject["rPSFluxMax"], 1)


class TestMaxSlopeDiaPsFlux(unittest.TestCase):

    def testCalculate(self):
        """Test flux maximum slope.
        """
        n_sources = 10
        diaObject = dict()
        fluxes = np.linspace(-1, 1, n_sources)
        times = np.linspace(0, 1, n_sources)
        diaSources = pd.DataFrame(data={"psFlux": fluxes,
                                        "psFluxErr": np.ones(n_sources),
                                        "midPointTai": times})

        plug = MaxSlopeDiaPsFlux(MaxSlopeDiaPsFluxConfig(),
                                 "ap_maxSlopeFlux",
                                 None)
        plug.calculate(diaObject, diaSources, diaSources, "u")
        self.assertEqual(diaObject["uPSFluxMaxSlope"], 2.)

        # test no inputs
        diaObject = dict()
        plug.calculate(diaObject, [], [], "g")
        self.assertTrue(np.isnan(diaObject["gPSFluxMaxSlope"]))

        # test one input
        diaObject = dict()
        diaSources = pd.DataFrame(data={"psFlux": [fluxes[0]],
                                        "psFluxErr": [np.ones(n_sources)[0]],
                                        "midPointTai": [times[0]]})
        plug.calculate(diaObject, diaSources, diaSources, "g")
        self.assertTrue(np.isnan(diaObject["gPSFluxMaxSlope"]))

        diaObject = dict()
        fluxes[4] = np.nan
        diaSources = pd.DataFrame(data={"psFlux": fluxes,
                                        "psFluxErr": np.ones(n_sources),
                                        "midPointTai": times})
        plug.calculate(diaObject, diaSources, diaSources, "r")
        self.assertEqual(diaObject["rPSFluxMaxSlope"], 2.)


class TestErrMeanDiaPsFlux(unittest.TestCase):

    def testCalculate(self):
        """Test error mean calculation.
        """
        n_sources = 10
        diaObject = dict()
        fluxes = np.linspace(-1, 1, n_sources)
        errors = np.linspace(1, 2, n_sources)
        diaSources = pd.DataFrame(data={"psFlux": fluxes,
                                        "psFluxErr": errors})

        plug = ErrMeanDiaPsFlux(ErrMeanDiaPsFluxConfig(),
                                "ap_errMeanFlux",
                                None)
        plug.calculate(diaObject, diaSources, diaSources, "u")
        self.assertEqual(diaObject["uPSFluxErrMean"], np.nanmean(errors))

        # test no inputs
        diaObject = dict()
        plug.calculate(diaObject, [], [], "g")
        self.assertTrue(np.isnan(diaObject["gPSFluxErrMean"]))

        diaObject = dict()
        errors[4] = np.nan
        diaSources = pd.DataFrame(data={"psFlux": fluxes,
                                        "psFluxErr": errors})
        plug.calculate(diaObject, diaSources, diaSources, "r")
        self.assertEqual(diaObject["rPSFluxErrMean"], np.nanmean(errors))


class TestLinearFitDiaPsFlux(unittest.TestCase):

    def testCalculate(self):
        """Test a linear fit to flux vs time.
        """
        n_sources = 10
        diaObject = dict()
        fluxes = np.linspace(-1, 1, n_sources)
        errors = np.linspace(1, 2, n_sources)
        times = np.linspace(0, 1, n_sources)
        diaSources = pd.DataFrame(data={"psFlux": fluxes,
                                        "psFluxErr": errors,
                                        "midPointTai": times})

        plug = LinearFitDiaPsFlux(LinearFitDiaPsFluxConfig(),
                                  "ap_LinearFit",
                                  None)
        plug.calculate(diaObject, diaSources, diaSources, "u")
        self.assertAlmostEqual(diaObject["uPSFluxLinearSlope"], 2.)
        self.assertAlmostEqual(diaObject["uPSFluxLinearIntercept"], -1.)

        # test no inputs
        diaObject = dict()
        plug.calculate(diaObject, [], [], "g")
        self.assertTrue(np.isnan(diaObject["gPSFluxLinearSlope"]))
        self.assertTrue(np.isnan(diaObject["gPSFluxLinearIntercept"]))

        # test no inputs
        diaObject = dict()
        diaSources = pd.DataFrame(data={"psFlux": [fluxes[0]],
                                        "psFluxErr": [errors[0]],
                                        "midPointTai": [times[0]]})
        plug.calculate(diaObject, diaSources, diaSources, "g")
        self.assertTrue(np.isnan(diaObject["gPSFluxLinearSlope"]))
        self.assertTrue(np.isnan(diaObject["gPSFluxLinearIntercept"]))

        diaObject = dict()
        fluxes[7] = np.nan
        errors[4] = np.nan
        diaSources = pd.DataFrame(data={"psFlux": fluxes,
                                        "psFluxErr": errors,
                                        "midPointTai": times})
        plug.calculate(diaObject, diaSources, diaSources, "r")
        self.assertAlmostEqual(diaObject["rPSFluxLinearSlope"], 2.)
        self.assertAlmostEqual(diaObject["rPSFluxLinearIntercept"], -1.)


class TestStetsonJDiaPsFlux(unittest.TestCase):

    def testCalculate(self):
        """Test the stetsonJ statistic.
        """
        n_sources = 10
        diaObject = dict(uPSFluxMean=0)
        fluxes = np.linspace(-1, 1, n_sources)
        errors = np.ones(n_sources)
        times = np.linspace(0, 1, n_sources)
        diaSources = pd.DataFrame(data={"psFlux": fluxes,
                                        "psFluxErr": errors,
                                        "midPointTai": times})

        plug = StetsonJDiaPsFlux(LinearFitDiaPsFluxConfig(),
                                 "ap_LinearFit",
                                 None)
        plug.calculate(diaObject, diaSources, diaSources, "u")
        # Expected StetsonJ for the values created. Confirmed using Cesimum's
        # implementation. http://github.com/cesium-ml/cesium
        self.assertAlmostEqual(diaObject["uPSFluxStetsonJ"],
                               -0.5958393936080928)

        # test no inputs
        diaObject = dict(gPSFluxMean=0)
        plug.calculate(diaObject, [], [], "g")
        self.assertTrue(np.isnan(diaObject["gPSFluxStetsonJ"]))

        # test no inputs
        diaObject = dict(gPSFluxMean=0)
        diaSources = pd.DataFrame(data={"psFlux": [fluxes[0]],
                                        "psFluxErr": [errors[0]],
                                        "midPointTai": [times[0]]})
        plug.calculate(diaObject, diaSources, diaSources, "g")
        self.assertTrue(np.isnan(diaObject["gPSFluxStetsonJ"]))

        fluxes[7] = np.nan
        errors[4] = np.nan
        nonNanMask = np.logical_and(~np.isnan(fluxes),
                                    ~np.isnan(errors))
        diaObject = dict(rPSFluxMean=np.average(fluxes[nonNanMask],
                                                weights=errors[nonNanMask]))
        diaSources = pd.DataFrame(data={"psFlux": fluxes,
                                        "psFluxErr": errors,
                                        "midPointTai": times})
        plug.calculate(diaObject, diaSources, diaSources, "r")
        self.assertAlmostEqual(diaObject["rPSFluxStetsonJ"],
                               -0.5412797916187173)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
