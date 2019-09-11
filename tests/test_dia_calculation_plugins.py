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
import pandas as pd
import unittest

from lsst.ap.association import (
    MeanDiaPosition,
    MeanDiaPositionConfig,
    WeightedMeanDiaPsFlux,
    WeightedMeanDiaPsFluxConfig)
import lsst.utils.tests


class TestMeanPosition(unittest.TestCase):

    def testCalculate(self):
        """Test that forced source catalogs are successfully created and have
        sensible values.
        """
        n_sources = 10

        # Test expected means.
        diaObject = dict()
        diaSources = pd.DataFrame(data={"ra": np.linspace(-1, 1, n_sources),
                                        "decl": np.zeros(n_sources)})
        mean_pos = MeanDiaPosition(MeanDiaPositionConfig(),
                                   "ap_meanPosition",
                                   None)
        mean_pos.calculate(diaObject, diaSources)

        self.assertAlmostEqual(diaObject["ra"], 0.0)
        self.assertAlmostEqual(diaObject["decl"], 0.0)

        diaObject = dict()
        diaSources = pd.DataFrame(data={"ra": np.zeros(n_sources),
                                        "decl": np.linspace(-1, 1, n_sources)})
        mean_pos.calculate(diaObject, diaSources)

        self.assertAlmostEqual(diaObject["ra"], 0.0)
        self.assertAlmostEqual(diaObject["decl"], 0.0)

        # Test failure modes.
        diaObject = dict()
        diaSources = pd.DataFrame(data={"ra": np.full(n_sources, np.nan),
                                        "decl": np.zeros(n_sources)})
        mean_pos.calculate(diaObject, diaSources)

        self.assertTrue(np.isnan(diaObject["ra"]))
        self.assertTrue(np.isnan(diaObject["decl"]))

        diaObject = dict()
        diaSources = pd.DataFrame(data={"ra": np.zeros(n_sources),
                                        "decl": np.full(n_sources, np.nan)})
        mean_pos.calculate(diaObject, diaSources)

        self.assertTrue(np.isnan(diaObject["ra"]))
        self.assertTrue(np.isnan(diaObject["decl"]))


class TestWeightedMeanDiaPsFlux(unittest.TestCase):

    def testCalculate(self):
        """Test that forced source catalogs are successfully created and have
        sensible values.
        """
        n_sources = 10
        diaObject = dict()
        diaSources = pd.DataFrame(data={"psFlux": np.linspace(-1, 1, n_sources),
                                        "psFluxErr": np.ones(n_sources)})

        mean_flux = WeightedMeanDiaPsFlux(WeightedMeanDiaPsFluxConfig(),
                                          "ap_meanFlux",
                                          None)
        mean_flux.calculate(diaObject, diaSources, diaSources, "u")

        self.assertAlmostEqual(diaObject["uPSFluxMean"], 0.0)
        self.assertAlmostEqual(diaObject["uPSFluxMeanErr"], np.sqrt(1 / n_sources))

        diaObject = dict()
        mean_flux.calculate(diaObject, [], [], "g")

        self.assertTrue(np.isnan(diaObject["gPSFluxMean"]))
        self.assertTrue(np.isnan(diaObject["gPSFluxMeanErr"]))

        diaObject = dict()
        diaSources.loc[4, "psFlux"] = np.nan
        mean_flux.calculate(diaObject, diaSources, diaSources, "r")

        self.assertTrue(~np.isnan(diaObject["rPSFluxMean"]))
        self.assertTrue(~np.isnan(diaObject["rPSFluxMeanErr"]))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
