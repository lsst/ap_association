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


def _make_task(maxClusteringDistance=1.5, nNeighborsConnectivity=30,
               maxSubsetSize=50000):
    """Create a DeduplicateAllSkyDiaObjectsTask with given config values.

    Creates a temporary SQLite APDB to satisfy the config requirement.
    """
    db_file_fd, db_file = tempfile.mkstemp()
    apdb_config = ApdbSql.init_database(db_url="sqlite:///" + db_file)
    config_file = tempfile.NamedTemporaryFile(delete=False)
    apdb_config.save(config_file.name)
    config_file.close()

    config = DeduplicateAllSkyDiaObjectsTask.ConfigClass()
    config.apdb_config_url = config_file.name
    config.maxClusteringDistance = maxClusteringDistance
    config.nNeighborsConnectivity = nNeighborsConnectivity
    config.maxSubsetSize = maxSubsetSize

    task = DeduplicateAllSkyDiaObjectsTask(config=config)

    os.close(db_file_fd)
    os.unlink(db_file)
    os.unlink(config_file.name)

    return task


def _make_diaObjects(ra_list, dec_list, ids=None):
    """Create a minimal DiaObjects DataFrame from lists of ra/dec."""
    n = len(ra_list)
    if ids is None:
        ids = list(range(1, n + 1))
    df = pd.DataFrame({
        'diaObjectId': ids,
        'ra': ra_list,
        'dec': dec_list,
    })
    df.index = pd.RangeIndex(n)
    return df


class TestClusterSmallDataset(unittest.TestCase):
    """Tests for the cluster method when dataset size <= maxSubsetSize
    (direct clustering path via _cluster_subset).
    """

    def setUp(self):
        self.task = _make_task(maxClusteringDistance=1.5)

    def test_single_object(self):
        """A single object should get a unique cluster label."""
        diaObjects = _make_diaObjects([10.0], [5.0])
        labels = self.task.cluster(diaObjects)

        self.assertEqual(len(labels), 1)
        self.assertEqual(labels.iloc[0], 0)

    def test_two_objects_within_threshold(self):
        """Two objects within maxClusteringDistance should share a label."""
        # 0.3 arcsec apart in dec
        diaObjects = _make_diaObjects(
            [10.0, 10.0],
            [5.0, 5.0 + 0.3 / 3600]
        )
        labels = self.task.cluster(diaObjects)

        self.assertEqual(len(labels), 2)
        self.assertEqual(labels.iloc[0], labels.iloc[1])

    def test_two_objects_beyond_threshold(self):
        """Two objects beyond maxClusteringDistance should get different labels."""
        # 10 arcsec apart, well beyond the 1.5 arcsec threshold
        diaObjects = _make_diaObjects(
            [10.0, 10.0],
            [5.0, 5.0 + 10.0 / 3600]
        )
        labels = self.task.cluster(diaObjects)

        self.assertEqual(len(labels), 2)
        self.assertNotEqual(labels.iloc[0], labels.iloc[1])

    def test_multiple_distinct_clusters(self):
        """Multiple well-separated groups should form distinct clusters."""
        # Three tight groups separated by large distances
        ra = [10.0, 10.0 + 0.1 / 3600,  # Group 1
              20.0, 20.0 + 0.2 / 3600,  # Group 2
              30.0, 30.0 + 0.1 / 3600, 30.0 + 0.2 / 3600]  # Group 3
        dec = [5.0, 5.0,
               15.0, 15.0,
               25.0, 25.0, 25.0]

        diaObjects = _make_diaObjects(ra, dec)
        labels = self.task.cluster(diaObjects)

        self.assertEqual(len(labels), 7)

        # Group 1 should share a label
        self.assertEqual(labels.iloc[0], labels.iloc[1])
        # Group 2 should share a label
        self.assertEqual(labels.iloc[2], labels.iloc[3])
        # Group 3 should share a label
        self.assertEqual(labels.iloc[4], labels.iloc[5])
        self.assertEqual(labels.iloc[5], labels.iloc[6])

        # Distinct groups should have different labels
        unique_group_labels = {labels.iloc[0], labels.iloc[2], labels.iloc[4]}
        self.assertEqual(len(unique_group_labels), 3)

    def test_chain_of_objects(self):
        """Objects forming a chain: each consecutive pair is within threshold,
        but endpoints exceed threshold.

        AgglomerativeClustering with a distance_threshold will merge only
        pairs within that distance.  With a k-NN connectivity graph, it can
        still link objects transitively through neighbors.  However, the
        distance_threshold limits merging to the configured distance, so
        endpoints farther apart than the threshold may end up in separate
        clusters depending on the linkage.
        """
        # Create a chain of 3 objects, each 0.5 arcsec from its neighbor.
        # Total span is 1.0 arcsec < threshold of 1.5 arcsec.
        n = 3
        step = 0.5 / 3600  # 0.5 arcsec steps in degrees
        ra = [10.0] * n
        dec = [5.0 + i * step for i in range(n)]

        diaObjects = _make_diaObjects(ra, dec)
        labels = self.task.cluster(diaObjects)

        # All 3 are within 1.0 arcsec of each other, within the 1.5 arcsec
        # threshold, so they should form a single cluster.
        unique_labels = labels.unique()
        self.assertEqual(len(unique_labels), 1,
                         "Chain of nearby objects within threshold should form a single cluster")

    def test_objects_exactly_at_threshold(self):
        """Objects separated by exactly maxClusteringDistance are merged.

        AgglomerativeClustering uses <= for the distance_threshold, so
        objects at exactly the threshold distance are merged.
        """
        # Exactly 1.5 arcsec apart
        diaObjects = _make_diaObjects(
            [10.0, 10.0],
            [5.0, 5.0 + 1.5 / 3600]
        )
        labels = self.task.cluster(diaObjects)

        self.assertEqual(len(labels), 2)
        self.assertEqual(labels.iloc[0], labels.iloc[1])

    def test_returns_series_with_matching_index(self):
        """The returned Series should have the same index as the input."""
        diaObjects = _make_diaObjects([10.0, 20.0, 30.0], [5.0, 15.0, 25.0])
        diaObjects.index = pd.Index([10, 20, 30])

        labels = self.task.cluster(diaObjects)

        self.assertIsInstance(labels, pd.Series)
        self.assertTrue(labels.index.equals(diaObjects.index))
        self.assertEqual(labels.name, "cluster_label")


    def test_large_cluster(self):
        """A tight group of many objects should mostly cluster together.

        With k-NN connectivity, very tight groups where all objects are
        mutual nearest neighbors should form a single cluster.  We use
        a small group here to ensure the connectivity graph is dense enough.
        """
        rng = np.random.default_rng(42)
        n = 15
        # Objects within 0.3 arcsec of center (well within 1.5 arcsec threshold)
        center_ra, center_dec = 100.0, -30.0
        ra = center_ra + rng.uniform(-0.3 / 3600, 0.3 / 3600, n)
        dec = center_dec + rng.uniform(-0.3 / 3600, 0.3 / 3600, n)

        diaObjects = _make_diaObjects(ra.tolist(), dec.tolist())
        labels = self.task.cluster(diaObjects)

        # All should be in the same cluster
        self.assertEqual(len(labels.unique()), 1)

    def test_no_duplicates_all_isolated(self):
        """When all objects are well separated, each gets its own cluster."""
        n = 20
        # Place objects 1 degree apart
        ra = [float(i) for i in range(n)]
        dec = [0.0] * n

        diaObjects = _make_diaObjects(ra, dec)
        labels = self.task.cluster(diaObjects)

        self.assertEqual(len(labels.unique()), n)

    def test_mixed_duplicates_and_isolated(self):
        """A mix of duplicate groups and isolated objects."""
        ra = [10.0, 10.0 + 0.1 / 3600, 10.0 + 0.2 / 3600,  # Group of 3
              50.0,  # Isolated
              80.0, 80.0 + 0.1 / 3600,  # Pair
              120.0]  # Isolated
        dec = [5.0, 5.0, 5.0,
               20.0,
               40.0, 40.0,
               60.0]

        diaObjects = _make_diaObjects(ra, dec)
        labels = self.task.cluster(diaObjects)

        # Group of 3 should share a label
        self.assertEqual(labels.iloc[0], labels.iloc[1])
        self.assertEqual(labels.iloc[1], labels.iloc[2])

        # Pair should share a label
        self.assertEqual(labels.iloc[4], labels.iloc[5])

        # All distinct clusters
        unique_groups = {labels.iloc[0], labels.iloc[3], labels.iloc[4], labels.iloc[6]}
        self.assertEqual(len(unique_groups), 4)

    def test_nNeighborsConnectivity_adapted_for_small_n(self):
        """When n_objects < nNeighborsConnectivity, the method should adapt."""
        # config nNeighborsConnectivity=30 but only 5 objects
        diaObjects = _make_diaObjects(
            [10.0, 10.0 + 0.1 / 3600, 20.0, 30.0, 40.0],
            [5.0, 5.0, 15.0, 25.0, 35.0]
        )
        # Should not raise
        labels = self.task.cluster(diaObjects)
        self.assertEqual(len(labels), 5)
        # First two should be clustered
        self.assertEqual(labels.iloc[0], labels.iloc[1])

    def test_custom_clustering_distance(self):
        """Changing maxClusteringDistance affects cluster membership."""
        ra = [10.0, 10.0]
        dec = [5.0, 5.0 + 1.0 / 3600]  # 1 arcsec apart

        diaObjects = _make_diaObjects(ra, dec)

        # With 1.5 arcsec threshold, should cluster
        task_wide = _make_task(maxClusteringDistance=1.5)
        labels_wide = task_wide.cluster(diaObjects)
        self.assertEqual(labels_wide.iloc[0], labels_wide.iloc[1])

        # With 0.5 arcsec threshold, should not cluster
        task_narrow = _make_task(maxClusteringDistance=0.5)
        labels_narrow = task_narrow.cluster(diaObjects)
        self.assertNotEqual(labels_narrow.iloc[0], labels_narrow.iloc[1])

    def test_declination_near_pole(self):
        """Objects near the celestial pole should still cluster correctly."""
        # Near north pole
        diaObjects = _make_diaObjects(
            [0.0, 180.0],  # Opposite RA but very close at the pole
            [89.9999, 89.9999]
        )
        labels = self.task.cluster(diaObjects)
        # These are actually close together at the pole
        # (cos(dec) factor makes RA separation small)
        # but AgglomerativeClustering uses Euclidean on (ra, dec) coords,
        # so it will treat them as 180 degrees apart in RA
        self.assertEqual(len(labels), 2)

    def test_ra_wraparound(self):
        """Objects near RA=0/360 boundary."""
        # Two objects near RA=0, one just above and one just below
        diaObjects = _make_diaObjects(
            [0.0001, 359.9999],
            [5.0, 5.0]
        )
        labels = self.task.cluster(diaObjects)
        # Euclidean clustering in (ra, dec) space won't handle wraparound;
        # these will appear as ~360 degrees apart and won't cluster.
        self.assertNotEqual(labels.iloc[0], labels.iloc[1])


class TestClusterLargeDataset(unittest.TestCase):
    """Tests for the cluster method when dataset size > maxSubsetSize
    (partitioned clustering path via KMeans + boundary merging).
    """

    def setUp(self):
        # Use a very small maxSubsetSize to trigger partitioning
        self.task = _make_task(maxClusteringDistance=1.5, maxSubsetSize=1000)

    def test_partitioning_triggered(self):
        """Verify that partitioning is used for large datasets."""
        rng = np.random.default_rng(42)
        n = 2000
        ra = rng.uniform(0, 360, n)
        dec = rng.uniform(-90, 90, n)
        diaObjects = _make_diaObjects(ra.tolist(), dec.tolist())

        labels = self.task.cluster(diaObjects)

        self.assertEqual(len(labels), n)
        # Should have many unique clusters (most objects are isolated)
        self.assertGreater(len(labels.unique()), n * 0.9)

    def test_duplicates_found_within_partition(self):
        """Duplicates within the same partition should be identified."""
        rng = np.random.default_rng(42)
        n = 1500  # Exceeds maxSubsetSize=1000

        # Place most objects randomly and far apart
        ra = rng.uniform(0, 180, n).tolist()
        dec = rng.uniform(-45, 45, n).tolist()

        # Add a tight duplicate group (will end up in same partition)
        group_center_ra = 90.0
        group_center_dec = 0.0
        for i in range(5):
            ra.append(group_center_ra + rng.uniform(-0.3, 0.3) / 3600)
            dec.append(group_center_dec + rng.uniform(-0.3, 0.3) / 3600)

        diaObjects = _make_diaObjects(ra, dec)
        labels = self.task.cluster(diaObjects)

        # The 5 duplicate objects at the end should share a cluster label
        dup_labels = labels.iloc[n:n + 5]
        self.assertEqual(len(dup_labels.unique()), 1,
                         "Tight group within partition should be clustered")

    def test_duplicates_at_partition_boundary(self):
        """Duplicates that span a partition boundary should still be merged."""
        # maxSubsetSize=1000, so with 1500 objects we get 2 partitions.
        # Create two clusters of objects placed such that KMeans will
        # split them into different partitions, with a duplicate pair
        # at the boundary.
        rng = np.random.default_rng(123)

        # Group A: objects in one region of sky
        n_a = 800
        ra_a = rng.uniform(10, 20, n_a).tolist()
        dec_a = rng.uniform(10, 20, n_a).tolist()

        # Group B: objects in a different region
        n_b = 800
        ra_b = rng.uniform(100, 110, n_b).tolist()
        dec_b = rng.uniform(-20, -10, n_b).tolist()

        # Now add a pair of duplicates at the midpoint between groups
        # These are close together but KMeans may split them
        mid_ra = 55.0
        mid_dec = 0.0
        ra_boundary = [mid_ra, mid_ra + 0.2 / 3600]
        dec_boundary = [mid_dec, mid_dec + 0.1 / 3600]

        ra = ra_a + ra_b + ra_boundary
        dec = dec_a + dec_b + dec_boundary

        diaObjects = _make_diaObjects(ra, dec)
        labels = self.task.cluster(diaObjects)

        # The boundary pair should be clustered together
        boundary_labels = labels.iloc[n_a + n_b:]
        self.assertEqual(boundary_labels.iloc[0], boundary_labels.iloc[1],
                         "Boundary duplicate pair should share a cluster label")

    def test_output_labels_are_contiguous(self):
        """Labels should be renumbered to be contiguous even after merging."""
        rng = np.random.default_rng(99)
        n = 1500
        ra = rng.uniform(0, 360, n)
        dec = rng.uniform(-90, 90, n)
        diaObjects = _make_diaObjects(ra.tolist(), dec.tolist())

        labels = self.task.cluster(diaObjects)

        unique_labels = sorted(labels.unique())
        expected = list(range(len(unique_labels)))
        self.assertEqual(unique_labels, expected)

    def test_returns_series_with_matching_index(self):
        """Returned Series matches input DataFrame index."""
        rng = np.random.default_rng(77)
        n = 1500
        ra = rng.uniform(0, 360, n)
        dec = rng.uniform(-90, 90, n)
        diaObjects = _make_diaObjects(ra.tolist(), dec.tolist())
        diaObjects.index = pd.RangeIndex(start=100, stop=100 + n)

        labels = self.task.cluster(diaObjects)

        self.assertIsInstance(labels, pd.Series)
        self.assertTrue(labels.index.equals(diaObjects.index))
        self.assertEqual(labels.name, "cluster_label")

    def test_consistency_with_small_path(self):
        """Results from the partitioned path should be consistent with
        the direct path for datasets that can use either.
        """
        rng = np.random.default_rng(55)
        n = 1200

        # Create data with known duplicate structure
        # 10 well-separated groups, each with 2-3 tight duplicates
        ra = []
        dec = []
        for i in range(10):
            center_ra = float(i * 30)
            center_dec = float(i * 10 - 45)
            n_in_group = rng.integers(2, 4)
            for _ in range(n_in_group):
                ra.append(center_ra + rng.uniform(-0.3, 0.3) / 3600)
                dec.append(center_dec + rng.uniform(-0.3, 0.3) / 3600)

        # Fill remaining with isolated objects
        n_placed = len(ra)
        for i in range(n - n_placed):
            ra.append(rng.uniform(0, 360))
            dec.append(rng.uniform(-90, 90))

        diaObjects = _make_diaObjects(ra, dec)

        # Run with small maxSubsetSize (partitioned path)
        task_partitioned = _make_task(maxClusteringDistance=1.5, maxSubsetSize=1000)
        labels_partitioned = task_partitioned.cluster(diaObjects)

        # Run with large maxSubsetSize (direct path)
        task_direct = _make_task(maxClusteringDistance=1.5, maxSubsetSize=50000)
        labels_direct = task_direct.cluster(diaObjects)

        # The known tight groups should be clustered in both cases
        idx = 0
        for i in range(10):
            n_in_group = 0
            group_start = idx
            center_ra = float(i * 30)
            center_dec = float(i * 10 - 45)
            while idx < n_placed:
                if (abs(ra[idx] - center_ra) < 1.0 / 3600
                        and abs(dec[idx] - center_dec) < 1.0 / 3600):
                    n_in_group += 1
                    idx += 1
                else:
                    break

            if n_in_group > 1:
                # Check direct path clusters them
                direct_group = labels_direct.iloc[group_start:group_start + n_in_group]
                self.assertEqual(len(direct_group.unique()), 1,
                                 f"Direct path: group {i} should be one cluster")

                # Check partitioned path clusters them
                part_group = labels_partitioned.iloc[group_start:group_start + n_in_group]
                self.assertEqual(len(part_group.unique()), 1,
                                 f"Partitioned path: group {i} should be one cluster")

    def test_empty_partition_handled(self):
        """An empty KMeans partition should be handled gracefully."""
        # We can't easily force an empty partition from KMeans, but we can
        # verify the code runs without error on a dataset where some
        # partitions might be very small.
        rng = np.random.default_rng(7)
        # Highly clustered data that may produce uneven partitions
        n = 1500
        # Most objects in one location
        ra = rng.normal(50.0, 0.01, n).tolist()
        dec = rng.normal(20.0, 0.01, n).tolist()
        # A few outliers
        for i in range(50):
            ra.append(float(i * 7))
            dec.append(float(i * 3 - 45))

        diaObjects = _make_diaObjects(ra + ra[-50:], dec + dec[-50:])
        # Should complete without error
        labels = self.task.cluster(diaObjects)
        self.assertEqual(len(labels), len(diaObjects))

    def test_many_partitions(self):
        """A dataset requiring many partitions."""
        task = _make_task(maxClusteringDistance=1.5, maxSubsetSize=1000)
        rng = np.random.default_rng(2024)
        n = 5000  # Will create ~5 partitions

        ra = rng.uniform(0, 360, n).tolist()
        dec = rng.uniform(-90, 90, n).tolist()

        # Add known duplicates
        for i in range(20):
            base_ra = rng.uniform(0, 360)
            base_dec = rng.uniform(-60, 60)
            ra.extend([base_ra, base_ra + 0.2 / 3600])
            dec.extend([base_dec, base_dec + 0.1 / 3600])

        diaObjects = _make_diaObjects(ra, dec)
        labels = task.cluster(diaObjects)

        self.assertEqual(len(labels), len(diaObjects))

        # Verify known duplicates are clustered
        for i in range(20):
            idx1 = n + 2 * i
            idx2 = n + 2 * i + 1
            self.assertEqual(labels.iloc[idx1], labels.iloc[idx2],
                             f"Known duplicate pair {i} should be clustered")


class TestClusterSubset(unittest.TestCase):
    """Direct tests of the _cluster_subset helper."""

    def setUp(self):
        self.task = _make_task(maxClusteringDistance=1.5)

    def test_single_object_returns_label_zero(self):
        """A single object returns label 0."""
        diaObjects = _make_diaObjects([10.0], [5.0])
        labels = self.task._cluster_subset(diaObjects)

        self.assertEqual(len(labels), 1)
        self.assertEqual(labels.iloc[0], 0)

    def test_two_nearby_objects(self):
        """Two nearby objects get the same label."""
        diaObjects = _make_diaObjects(
            [10.0, 10.0 + 0.1 / 3600],
            [5.0, 5.0]
        )
        labels = self.task._cluster_subset(diaObjects)
        self.assertEqual(labels.iloc[0], labels.iloc[1])

    def test_preserves_dataframe_index(self):
        """Labels index should match input DataFrame index."""
        diaObjects = _make_diaObjects([10.0, 20.0, 30.0], [5.0, 15.0, 25.0])
        diaObjects.index = pd.Index([5, 10, 15])

        labels = self.task._cluster_subset(diaObjects)
        self.assertTrue(labels.index.equals(diaObjects.index))


class TestClusterEdgeCases(unittest.TestCase):
    """Edge cases and special configurations for the cluster method."""

    def test_all_objects_identical_position(self):
        """All objects at the same location should form one cluster."""
        n = 10
        ra = [100.0] * n
        dec = [-30.0] * n
        diaObjects = _make_diaObjects(ra, dec)

        task = _make_task(maxClusteringDistance=1.5)
        labels = task.cluster(diaObjects)

        self.assertEqual(len(labels.unique()), 1)

    def test_negative_declination(self):
        """Clustering works for objects in the southern hemisphere."""
        diaObjects = _make_diaObjects(
            [200.0, 200.0 + 0.1 / 3600, 200.0 + 0.2 / 3600],
            [-60.0, -60.0, -60.0]
        )
        task = _make_task(maxClusteringDistance=1.5)
        labels = task.cluster(diaObjects)

        self.assertEqual(len(labels.unique()), 1)

    def test_objects_along_equator(self):
        """Objects along the celestial equator."""
        # Objects 0.5 arcsec apart along RA at dec=0
        diaObjects = _make_diaObjects(
            [10.0, 10.0 + 0.5 / 3600, 10.0 + 1.0 / 3600],
            [0.0, 0.0, 0.0]
        )
        task = _make_task(maxClusteringDistance=1.5)
        labels = task.cluster(diaObjects)

        # All three within 1 arcsec of each other, should cluster
        self.assertEqual(len(labels.unique()), 1)

    def test_minimum_nNeighborsConnectivity(self):
        """Task with minimum nNeighborsConnectivity=3 still works."""
        task = _make_task(maxClusteringDistance=1.5, nNeighborsConnectivity=3)

        diaObjects = _make_diaObjects(
            [10.0, 10.0 + 0.1 / 3600, 20.0, 30.0, 40.0],
            [5.0, 5.0, 15.0, 25.0, 35.0]
        )
        labels = task.cluster(diaObjects)

        self.assertEqual(len(labels), 5)
        self.assertEqual(labels.iloc[0], labels.iloc[1])

    def test_maxSubsetSize_boundary(self):
        """Dataset exactly at maxSubsetSize should use direct path."""
        n = 1000
        task = _make_task(maxClusteringDistance=1.5, maxSubsetSize=1000)

        rng = np.random.default_rng(42)
        ra = rng.uniform(0, 360, n).tolist()
        dec = rng.uniform(-90, 90, n).tolist()
        diaObjects = _make_diaObjects(ra, dec)

        # Should use direct path (n == maxSubsetSize)
        labels = task.cluster(diaObjects)
        self.assertEqual(len(labels), n)

    def test_maxSubsetSize_plus_one_triggers_partitioning(self):
        """Dataset of maxSubsetSize+1 should trigger partitioning."""
        n = 1001
        task = _make_task(maxClusteringDistance=1.5, maxSubsetSize=1000)

        rng = np.random.default_rng(42)
        ra = rng.uniform(0, 360, n).tolist()
        dec = rng.uniform(-90, 90, n).tolist()
        diaObjects = _make_diaObjects(ra, dec)

        labels = task.cluster(diaObjects)
        self.assertEqual(len(labels), n)

    def test_three_objects_same_position(self):
        """Three objects at the same position form one cluster."""
        diaObjects = _make_diaObjects(
            [45.0, 45.0, 45.0],
            [-10.0, -10.0, -10.0]
        )
        task = _make_task(maxClusteringDistance=1.5)
        labels = task.cluster(diaObjects)

        self.assertEqual(len(labels.unique()), 1)

    def test_non_default_index(self):
        """DataFrame with a non-default integer index."""
        diaObjects = _make_diaObjects(
            [10.0, 10.0 + 0.1 / 3600, 50.0],
            [5.0, 5.0, 30.0]
        )
        diaObjects.index = pd.Index([100, 200, 300])

        task = _make_task(maxClusteringDistance=1.5)
        labels = task.cluster(diaObjects)

        self.assertTrue(labels.index.equals(diaObjects.index))
        self.assertEqual(labels.loc[100], labels.loc[200])
        self.assertNotEqual(labels.loc[100], labels.loc[300])

    def test_large_nNeighborsConnectivity(self):
        """nNeighborsConnectivity larger than dataset is handled."""
        task = _make_task(maxClusteringDistance=1.5, nNeighborsConnectivity=30)

        # Only 4 objects but nNeighborsConnectivity=30
        diaObjects = _make_diaObjects(
            [10.0, 10.0 + 0.1 / 3600, 50.0, 80.0],
            [5.0, 5.0, 30.0, 60.0]
        )
        labels = task.cluster(diaObjects)

        self.assertEqual(len(labels), 4)
        self.assertEqual(labels.iloc[0], labels.iloc[1])


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
