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

        With a maxTrailLength config of 0.01 arcseconds/second and an
        exposure of 30 seconds,the max trail length is 0.3 arcseconds. Only the
        source with a trail of 0 stays in the catalog and the rest are filtered
        out and put into results.trailedSources.
        """
        config = TrailedSourceFilterTask.ConfigClass()
        config.maxTrailLength = 0.01
        trailedSourceFilterTask = TrailedSourceFilterTask(config=config)
        results = trailedSourceFilterTask.run(self.diaSources, self.exposure_time)

        self.assertEqual(len(results.diaSources), 1)
        np.testing.assert_array_equal(results.diaSources['diaSourceId'].values, [0])
        np.testing.assert_array_equal(results.trailedDiaSources['diaSourceId'].values, [1, 2, 3, 4])

    def test_run_no_trails(self):
        """Run trailedSourceFilterTask with a long trail length so that
        every source in the catalog is in the final diaSource catalog.

        With a maxTrailLength config of 10 arcseconds/second and an
        exposure of 30 seconds,the max trail length is 300 arcseconds. All
        sources in the initial catalog should be in the final diaSource
        catalog.
        """
        config = TrailedSourceFilterTask.ConfigClass()
        config.maxTrailLength = 10.00
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
