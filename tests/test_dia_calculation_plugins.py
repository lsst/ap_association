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
    SkewDiaPsFlux, SkewDiaPsFluxConfig)
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


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
