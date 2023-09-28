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

import unittest
from lsst.ap.association import TrailedSourceFilterTask
import numpy as np
import pandas as pd
import lsst.utils.tests


class TestTrailedSourceFilterTask(unittest.TestCase):

    def setUp(self):
        """Create sets of diaSources.

        The trail lengths of the dia sources are 0, 5.5, 11, 16.5, 21.5
        arcseconds.
        """
        rng = np.random.default_rng(1234)
        scatter = 0.1 / 3600
        self.nSources = 5
        self.diaSources = pd.DataFrame(data=[
            {"ra": 0.04*idx + scatter*rng.uniform(-1, 1),
             "dec": 0.04*idx + scatter*rng.uniform(-1, 1),
             "diaSourceId": idx, "diaObjectId": 0, "trailLength": 5.5*idx}
            for idx in range(self.nSources)])
        self.exposure_time = 30.0

    def test_run(self):
        """Run trailedSourceFilterTask with the default max distance.

        With the default settings and an exposure of 30 seconds, the max trail
        length is 12.5 arcseconds. Two out of five of the diaSources will be
        filtered out of the final results and put into results.trailedSources.
        """
        trailedSourceFilterTask = TrailedSourceFilterTask()
        results = trailedSourceFilterTask.run(self.diaSources, self.exposure_time)

        self.assertEqual(len(results.diaSources), 3)
        np.testing.assert_array_equal(results.diaSources['diaSourceId'].values, [0, 1, 2])
        np.testing.assert_array_equal(results.trailedDiaSources['diaSourceId'].values, [3, 4])

    def test_run_short_max_trail(self):
        """Run trailedSourceFilterTask with aggressive trail length cutoff

        With a max_trail_length config of 0.01 arcseconds/second and an
        exposure of 30 seconds,the max trail length is 0.3 arcseconds. Only the
        source with a trail of 0 stays in the catalog and the rest are filtered
        out and put into results.trailedSources.
        """
        config = TrailedSourceFilterTask.ConfigClass()
        config.max_trail_length = 0.01
        trailedSourceFilterTask = TrailedSourceFilterTask(config=config)
        results = trailedSourceFilterTask.run(self.diaSources, self.exposure_time)

        self.assertEqual(len(results.diaSources), 1)
        np.testing.assert_array_equal(results.diaSources['diaSourceId'].values, [0])
        np.testing.assert_array_equal(results.trailedDiaSources['diaSourceId'].values, [1, 2, 3, 4])

    def test_run_no_trails(self):
        """Run trailedSourceFilterTask with a long trail length so that
        every source in the catalog is in the final diaSource catalog.

        With a max_trail_length config of 10 arcseconds/second and an
        exposure of 30 seconds,the max trail length is 300 arcseconds. All
        sources in the initial catalog should be in the final diaSource
        catalog.
        """
        config = TrailedSourceFilterTask.ConfigClass()
        config.max_trail_length = 10.00
        trailedSourceFilterTask = TrailedSourceFilterTask(config=config)
        results = trailedSourceFilterTask.run(self.diaSources, self.exposure_time)

        self.assertEqual(len(results.diaSources), 5)
        self.assertEqual(len(results.trailedDiaSources), 0)
        np.testing.assert_array_equal(results.diaSources["diaSourceId"].values, [0, 1, 2, 3, 4])
        np.testing.assert_array_equal(results.trailedDiaSources["diaSourceId"].values, [])

    def test_check_dia_source_trail(self):
        """Test the source trail mask filter.

        Test that the mask filter returns the expected mask array.
        """
        trailedSourceFilterTask = TrailedSourceFilterTask()
        mask = trailedSourceFilterTask._check_dia_source_trail(self.diaSources, self.exposure_time)
        np.testing.assert_array_equal(mask, [False, False, False, True, True])


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
