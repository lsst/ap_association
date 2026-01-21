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

import os
import astropy.units
import numpy as np
import pandas as pd
import tempfile
import unittest

from lsst.ap.association import DeduplicateAllSkyDiaObjectsTask
from lsst.ap.association.utils import getMidpointFromTimespan, readSchemaFromApdb
from lsst.dax.apdb import Apdb, ApdbSql
from lsst.utils import getPackageDir
import lsst.utils.tests
from utils_tests import makeExposure, makeDiaObjects, makeDiaSources, makeDiaForcedSources, makeRegionTime


def _data_file_name(basename, module_name):
    """Return path name of a data file.

    Parameters
    ----------
    basename : `str`
        Name of the file to add to the path string.
    module_name : `str`
        Name of lsst stack package environment variable.

    Returns
    -------
    data_file_path : `str`
       Full path of the file to load from the "data" directory in a given
       repository.
    """
    return os.path.join(getPackageDir(module_name), "data", basename)


class TestDeduplicateAllSkyDiaObjects(unittest.TestCase):

    def setUp(self):
        # Create an instance of random generator with fixed seed.
        rng = np.random.default_rng(1234)

        self.db_file_fd, self.db_file = tempfile.mkstemp(
            dir=os.path.dirname(__file__))
        self.addCleanup(os.remove, self.db_file)
        self.addCleanup(os.close, self.db_file_fd)

        self.apdbConfig = ApdbSql.init_database(db_url="sqlite:///" + self.db_file)
        self.config_file = tempfile.NamedTemporaryFile()
        self.addCleanup(self.config_file.close)
        self.apdbConfig.save(self.config_file.name)
        self.apdb = Apdb.from_config(self.apdbConfig)
        self.schema = readSchemaFromApdb(self.apdb)

        self.exposure = makeExposure(False, False)
        self.regionTime = makeRegionTime(exposure=self.exposure)
        self.dateTime = getMidpointFromTimespan(self.regionTime.timespan)

        original_diaObjects = makeDiaObjects(50, self.exposure, rng)

        # create duplicates
        dup_0_0 = {"ra": original_diaObjects.iloc[0]['ra'],
                   "dec": original_diaObjects.iloc[0]['dec'] + 0.2/3600,
                   "diaObjectId": 1000,
                   "nDiaSources": 1}
        for f in ["u", "g", "r", "i", "z", "y"]:
            dup_0_0["%s_psfFluxNdata" % f] = 0
        dup_0_1 = {"ra": original_diaObjects.iloc[0]['ra'] + 0.1/3600,
                   "dec": original_diaObjects.iloc[0]['dec'] + 0.2/3600,
                   "diaObjectId": 1001,
                   "nDiaSources": 1}
        for f in ["u", "g", "r", "i", "z", "y"]:
            dup_0_1["%s_psfFluxNdata" % f] = 0
        dup_0_2 = {"ra": original_diaObjects.iloc[0]['ra'],
                   "dec": original_diaObjects.iloc[0]['dec'] - 0.1/3600,
                   "diaObjectId": 1002,
                   "nDiaSources": 1}
        for f in ["u", "g", "r", "i", "z", "y"]:
            dup_0_2["%s_psfFluxNdata" % f] = 0
        dup_0_3 = {"ra": original_diaObjects.iloc[0]['ra'] + 0.2/3600,
                   "dec": original_diaObjects.iloc[0]['dec'] - 0.1/3600,
                   "diaObjectId": 1003,
                   "nDiaSources": 1}
        for f in ["u", "g", "r", "i", "z", "y"]:
            dup_0_3["%s_psfFluxNdata" % f] = 0

        dup_1_0 = {"ra": original_diaObjects.iloc[1]['ra'],
                   "dec": original_diaObjects.iloc[1]['dec'] + 0.2/3600,
                   "diaObjectId": 2000,
                   "nDiaSources": 1}
        for f in ["u", "g", "r", "i", "z", "y"]:
            dup_1_0["%s_psfFluxNdata" % f] = 0
        dup_1_1 = {"ra": original_diaObjects.iloc[1]['ra'] + 0.1/3600,
                   "dec": original_diaObjects.iloc[1]['dec'] + 0.2/3600,
                   "diaObjectId": 2001,
                   "nDiaSources": 1}
        for f in ["u", "g", "r", "i", "z", "y"]:
            dup_1_1["%s_psfFluxNdata" % f] = 0
        dup_1_2 = {"ra": original_diaObjects.iloc[1]['ra'],
                   "dec": original_diaObjects.iloc[1]['dec'] - 0.1/3600,
                   "diaObjectId": 2002,
                   "nDiaSources": 1}
        for f in ["u", "g", "r", "i", "z", "y"]:
            dup_1_2["%s_psfFluxNdata" % f] = 0

        dup_2_0 = {"ra": original_diaObjects.iloc[2]['ra'],
                   "dec": original_diaObjects.iloc[2]['dec'] + 0.2/3600,
                   "diaObjectId": 3000,
                   "nDiaSources": 1}
        for f in ["u", "g", "r", "i", "z", "y"]:
            dup_2_0["%s_psfFluxNdata" % f] = 0

        data = [dup_0_0, dup_0_1, dup_0_2, dup_0_3,
                dup_1_0, dup_1_1, dup_1_2,
                dup_2_0]
        duplicate_diaObjects = pd.DataFrame(data=data).set_index("diaObjectId", drop=False)

        self.diaObjects = pd.concat([original_diaObjects, duplicate_diaObjects])

        self.diaSources = makeDiaSources(
            100, self.diaObjects["diaObjectId"].to_numpy(), self.exposure, rng)
        self.diaForcedSources = makeDiaForcedSources(
            200, self.diaObjects["diaObjectId"].to_numpy(), self.exposure, rng)

        # Store the test diaSources as though they were observed a month before
        # the current exposure.
        dateTime = self.regionTime.timespan.begin.tai - 30 * astropy.units.day
        self.apdb.store(dateTime,
                        self.diaObjects,
                        self.diaSources,
                        self.diaForcedSources)

        # These columns are not in the DPDD, yet do appear in DiaSource.yaml.
        # We don't need to check them against the default APDB schema.
        self.ignoreColumns = ["band", "bboxSize", "isDipole", "flags"]

    def _makeConfig(self, **kwargs):
        config = DeduplicateAllSkyDiaObjectsTask.ConfigClass()
        config.apdb_config_url = self.config_file.name
        config.earliestMidpointMjdTai = self.dateTime.mjd
        config.update(**kwargs)
        return config

    def testRun(self):
        """Test the full run method for deduplication.
        """
        diaConfig = self._makeConfig()
        task = DeduplicateAllSkyDiaObjectsTask(config=diaConfig)
        result = task.run()

        self.assertIsNotNone(result.diaObjectDeduplicationMap)

        self.assertEqual(len(result.diaObjectDeduplicationMap), 8)

        self.assertIn('removedDiaObjectId', result.diaObjectDeduplicationMap.columns)
        self.assertIn('keptDiaObjectId', result.diaObjectDeduplicationMap.columns)

        removed_ids = set(result.diaObjectDeduplicationMap['removedDiaObjectId'])
        expected_duplicate_ids = {1000, 1001, 1002, 1003, 2000, 2001, 2002, 3000}
        self.assertEqual(removed_ids, expected_duplicate_ids)

        kept_ids = set(result.diaObjectDeduplicationMap['keptDiaObjectId'])
        self.assertTrue(kept_ids.issubset({1, 2, 3}))

        # Check that diaSources have been remapped from removed diaObjectIds
        # (We can't easily compare to the kept ids because non-duplicated ids
        # are not included in the deduplication map.)
        updated_diaSources = self.apdb.getDiaSources(region=self.regionTime.region, object_ids=None,
                                                     visit_time=self.dateTime, start_time=None)
        source_object_ids = set(updated_diaSources['diaObjectId'])
        self.assertTrue(source_object_ids.isdisjoint(removed_ids),
                        "DiaSources should not reference removed diaObjectIds")

    def test_remap_clusters(self):
        """Test the cluster remapping logic.
        """
        diaConfig = self._makeConfig()
        task = DeduplicateAllSkyDiaObjectsTask(config=diaConfig)

        # Create a simple test DataFrame with cluster labels
        test_data = pd.DataFrame({
            'diaObjectId': [1, 2, 3, 4, 5, 6],
            'validityStartMjdTai': [100.0, 101.0, 102.0, 200.0, 201.0, 300.0],
            'cluster_label': [0, 0, 0, 1, 1, 2],
            'ra': [10.0, 10.0, 10.0, 20.0, 20.0, 30.0],
            'dec': [5.0, 5.0, 5.0, 10.0, 10.0, 15.0]
        })

        dedup_map = task.remap_clusters(test_data)

        self.assertEqual(len(dedup_map), 3)  # 2 + 1 remappings

        # Check that objects 2 and 3 map to object 1
        obj2_mapping = dedup_map[dedup_map['removedDiaObjectId'] == 2]
        obj3_mapping = dedup_map[dedup_map['removedDiaObjectId'] == 3]
        self.assertEqual(len(obj2_mapping), 1)
        self.assertEqual(len(obj3_mapping), 1)
        self.assertEqual(obj2_mapping.iloc[0]['keptDiaObjectId'], 1)
        self.assertEqual(obj3_mapping.iloc[0]['keptDiaObjectId'], 1)

        # Check that object 5 maps to object 4
        obj5_mapping = dedup_map[dedup_map['removedDiaObjectId'] == 5]
        self.assertEqual(len(obj5_mapping), 1)
        self.assertEqual(obj5_mapping.iloc[0]['keptDiaObjectId'], 4)

    def test_cluster(self):
        """Test the clustering algorithm.
        """
        diaConfig = self._makeConfig(maxClusteringDistance=1.5)
        task = DeduplicateAllSkyDiaObjectsTask(config=diaConfig)

        cluster_labels = task.cluster(self.diaObjects)

        self.assertEqual(len(cluster_labels), len(self.diaObjects))

        group1_ids = [1, 1000, 1001, 1002, 1003]
        group1_mask = self.diaObjects['diaObjectId'].isin(group1_ids)
        group1_labels = cluster_labels[group1_mask]

        self.assertEqual(len(group1_labels.unique()), 1,
                         "Objects in duplicate group 1 should share a cluster label")

        group2_ids = [2, 2000, 2001, 2002]
        group2_mask = self.diaObjects['diaObjectId'].isin(group2_ids)
        group2_labels = cluster_labels[group2_mask]

        self.assertEqual(len(group2_labels.unique()), 1,
                         "Objects in duplicate group 2 should share a cluster label")

        group3_ids = [3, 3000]
        group3_mask = self.diaObjects['diaObjectId'].isin(group3_ids)
        group3_labels = cluster_labels[group3_mask]

        self.assertEqual(len(group3_labels.unique()), 1,
                         "Objects in duplicate group 3 should share a cluster label")

        self.assertNotEqual(group1_labels.iloc[0], group2_labels.iloc[0],
                            "Group 1 and Group 2 should have different cluster labels")
        self.assertNotEqual(group1_labels.iloc[0], group3_labels.iloc[0],
                            "Group 1 and Group 3 should have different cluster labels")
        self.assertNotEqual(group2_labels.iloc[0], group3_labels.iloc[0],
                            "Group 2 and Group 3 should have different cluster labels")

    def test_count_duplicates(self):
        """Test the duplicate counting method.
        """
        diaConfig = self._makeConfig(maxClusteringDistance=1.5)
        task = DeduplicateAllSkyDiaObjectsTask(config=diaConfig)

        test_catalog = pd.DataFrame({
            'diaObjectId': [1, 2, 3, 4, 5],
            'ra': [10.0, 10.0001,  # Pair 1
                   10.1, 10.1001,  # Pair 2
                   10.3],          # Single
            'dec': [5.0, 5.0001,
                    5.1, 5.1001,
                    5.3]
        })

        duplicate_count = task.count_duplicates(test_catalog)

        self.assertEqual(duplicate_count, 4)

        # Test with custom radius
        duplicate_count_small = task.count_duplicates(test_catalog, radius_arcsec=0.1)
        self.assertEqual(duplicate_count_small, 0)

        # Test with large radius
        duplicate_count_large = task.count_duplicates(test_catalog, radius_arcsec=3600.0)
        self.assertEqual(duplicate_count_large, 5)


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
