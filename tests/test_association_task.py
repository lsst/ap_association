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
import lsst.geom as geom
import lsst.utils.tests

from lsst.ap.association import AssociationTask


class TestAssociationTask(unittest.TestCase):

    def setUp(self):
        """Create sets of diaSources and diaObjects.
        """
        rng = np.random.default_rng(1234)
        self.nObjects = 5
        scatter = 0.1/3600
        self.diaObjects = pd.DataFrame(data=[
            {"ra": 0.04*(idx + 1), "dec": 0.04*(idx + 1),
             "diaObjectId": idx + 1}
            for idx in range(self.nObjects)])
        self.diaObjects.set_index("diaObjectId", drop=False, inplace=True)
        self.nSources = 5
        self.diaSources = pd.DataFrame(data=[
            {"ra": 0.04*idx + scatter*rng.uniform(-1, 1),
             "dec": 0.04*idx + scatter*rng.uniform(-1, 1),
             "diaSourceId": idx + 1 + self.nObjects, "diaObjectId": 0, "trailLength": 5.5*idx}
            for idx in range(self.nSources)])
        self.diaSourceZeroScatter = pd.DataFrame(data=[
            {"ra": 0.04*idx,
             "dec": 0.04*idx,
             "diaSourceId": idx + 1 + self.nObjects, "diaObjectId": 0, "trailLength": 5.5*idx}
            for idx in range(self.nSources)])
        self.exposure_time = 30.0

    def test_run(self):
        """Test the full task by associating a set of diaSources to
        existing diaObjects.
        """
        config = AssociationTask.ConfigClass()
        config.doTrailedSourceFilter = False
        assocTask = AssociationTask(config=config)
        results = assocTask.run(self.diaSources, self.diaObjects, exposure_time=self.exposure_time)

        self.assertEqual(results.nUpdatedDiaObjects, len(self.diaObjects) - 1)
        self.assertEqual(results.nUnassociatedDiaObjects, 1)
        self.assertEqual(len(results.matchedDiaSources),
                         len(self.diaObjects) - 1)
        self.assertEqual(len(results.unAssocDiaSources), 1)
        for test_obj_id, expected_obj_id in zip(
                results.matchedDiaSources["diaObjectId"].to_numpy(),
                [1, 2, 3, 4]):
            self.assertEqual(test_obj_id, expected_obj_id)
        for test_obj_id, expected_obj_id in zip(
                results.unAssocDiaSources["diaObjectId"].to_numpy(),
                [0]):
            self.assertEqual(test_obj_id, expected_obj_id)

    def test_run_trailed_sources(self):
        """Test the full task by associating a set of diaSources to
        existing diaObjects when trailed sources are filtered.

        This should filter out two of the five sources based on trail length,
        leaving one unassociated diaSource and two associated diaSources.
        """
        assocTask = AssociationTask()
        results = assocTask.run(self.diaSources, self.diaObjects, exposure_time=self.exposure_time)

        self.assertEqual(results.nUpdatedDiaObjects, len(self.diaObjects) - 3)
        self.assertEqual(results.nUnassociatedDiaObjects, 3)
        self.assertEqual(len(results.matchedDiaSources), len(self.diaObjects) - 3)
        self.assertEqual(len(results.unAssocDiaSources), 1)
        np.testing.assert_array_equal(results.matchedDiaSources["diaObjectId"].values, [1, 2])
        np.testing.assert_array_equal(results.unAssocDiaSources["diaObjectId"].values, [0])

    def test_run_no_existing_objects(self):
        """Test the run method with a completely empty database.
        """
        assocTask = AssociationTask()
        results = assocTask.run(
            self.diaSources,
            pd.DataFrame(columns=["ra", "dec", "diaObjectId", "trailLength"]),
            exposure_time=self.exposure_time)
        self.assertEqual(results.nUpdatedDiaObjects, 0)
        self.assertEqual(results.nUnassociatedDiaObjects, 0)
        self.assertEqual(len(results.matchedDiaSources), 0)
        self.assertTrue(np.all(results.unAssocDiaSources["diaObjectId"] == 0))

    def test_associate_sources(self):
        """Test the performance of the associate_sources method in
        AssociationTask.
        """
        assoc_task = AssociationTask()
        assoc_result = assoc_task.associate_sources(
            self.diaObjects, self.diaSources)

        for test_obj_id, expected_obj_id in zip(
                assoc_result.diaSources["diaObjectId"].to_numpy(),
                [0, 1, 2, 3, 4]):
            self.assertEqual(test_obj_id, expected_obj_id)

    def test_score_and_match(self):
        """Test association between a set of sources and an existing
        DIAObjectCollection.
        """

        assoc_task = AssociationTask()
        score_struct = assoc_task.score(self.diaObjects,
                                        self.diaSourceZeroScatter,
                                        1.0 * geom.arcseconds)
        self.assertFalse(np.isfinite(score_struct.scores[0]))
        for src_idx in range(1, len(self.diaSources)):
            # Our scores should be extremely close to 0 but not exactly so due
            # to machine noise.
            self.assertAlmostEqual(score_struct.scores[src_idx], 0.0,
                                   places=16)

        # After matching each DIAObject should now contain 2 DIASources
        # except the last DIAObject in this collection which should be
        # newly created during the matching step and contain only one
        # DIASource.
        match_result = assoc_task.match(
            self.diaObjects, self.diaSources, score_struct)
        self.assertEqual(match_result.nUpdatedDiaObjects, 4)
        self.assertEqual(match_result.nUnassociatedDiaObjects, 1)

    def test_remove_nan_dia_sources(self):
        """Test removing DiaSources with NaN locations.
        """
        self.diaSources.loc[2, "ra"] = np.nan
        self.diaSources.loc[3, "dec"] = np.nan
        self.diaSources.loc[4, "ra"] = np.nan
        self.diaSources.loc[4, "dec"] = np.nan
        assoc_task = AssociationTask()
        out_dia_sources = assoc_task.check_dia_source_radec(self.diaSources)
        self.assertEqual(len(out_dia_sources), len(self.diaSources) - 3)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
