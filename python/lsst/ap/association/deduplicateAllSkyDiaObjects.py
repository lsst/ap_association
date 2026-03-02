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

from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering, MiniBatchKMeans
from sklearn.neighbors import kneighbors_graph, BallTree

import lsst.dax.apdb as daxApdb
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

__all__ = ("DeduplicateAllSkyDiaObjectsTask", "DeduplicateAllSkyDiaObjectsConfig")


class DeduplicateAllSkyDiaObjectsConfig(pexConfig.Config):
    """Configuration for DeduplicateAllSkyDiaObjectsTask.

    Parameters
    ----------
    maxClusteringDistance : `float`
        Maximum distance to merge clusters of duplicate DiaObjects, in
        arcseconds. Default is 1.5 arcseconds.
    earliestMidpointMjdTai : `float`
        MidpointMjdTai (Modified Julian Date in TAI scale) of the
        earliest DiaSource to be reassigned during deduplication.
        Default corresponds to 2025-09-01.
    maxDiaObjects : `int`
        Maximum number of DiaObjects to process. An exception is raised
        if more DiaObjects are found to prevent excessive memory/computation.
    """
    apdb_config_url = pexConfig.Field(
        dtype=str,
        default=None,
        optional=False,
        doc="A config file specifying the APDB and its connection parameters, "
            "typically written by the apdb-cli command-line utility. "
            "The database must already be initialized.",
    )
    maxClusteringDistance = pexConfig.RangeField(
        doc="Maximum distence to merge clusters of duplicate diaObjects (arcseconds)",
        dtype=float,
        default=1.5,
        min=0,
    )
    earliestMidpointMjdTai = pexConfig.Field(
        dtype=float,
        default=Time('2025-09-01', scale='tai').mjd,
        doc="MidpointMjdTai of earliest DIASource to be reassigned during "
            "deduplication.",
    )
    maxDiaObjects = pexConfig.RangeField(
        doc="Maximum number of DiaObjects to process. An exception is raised "
            "if more DiaObjects are found to prevent excessive memory/computation.",
        dtype=int,
        default=1000000,
        min=1,
    )
    nNeighborsConnectivity = pexConfig.RangeField(
        doc="Number of neighbors to use for clustering.  Larger values are more"
            "accurate but more expensive computationally.",
        dtype=int,
        default=30,
        min=3,
    )
    maxSubsetSize = pexConfig.RangeField(
        doc="Maximum number of DiaObjects to process in a single subset during "
            "clustering. Larger values use more memory but may produce better "
            "results at cluster boundaries.",
        dtype=int,
        default=50000,
        min=1000,
    )


class DeduplicateAllSkyDiaObjectsTask(pipeBase.Task):
    """Deduplicate DiaObjects by clustering spatially and update the APDB.

    This task identifies and merges duplicate DiaObjects in the Alert
    Production Database (APDB) that are spatially coincident. It uses
    agglomerative clustering to group nearby objects, then reassigns
    DiaSources from duplicate objects to a single kept DiaObject (the
    oldest one in each cluster).
    """
    ConfigClass = DeduplicateAllSkyDiaObjectsConfig
    _DefaultName = "deduplicateAllSkyDiaObjects"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.apdb = daxApdb.Apdb.from_uri(self.config.apdb_config_url)

    def run(self):
        """Load DiaObjects and cluster them to remove duplicates.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components:

            ``diaObjectDeduplicationMap`` : `astropy.table.Table`
                Table mapping duplicate DiaObject IDs to the IDs they
                should be merged into. Contains columns
                'removedDiaObjectId' and 'keptDiaObjectId'.

        Raises
        ------
        lsst.pipe.base.NoWorkFound
            If no duplicate DiaObjects are found in the APDB.
        """

        # Return all DiaObject versions created since the last deduplication
        # run--typically this will be since yesterday morning.
        # TODO: consider adding a config for the "since" argument
        all_diaObjects = self.apdb.getDiaObjectsForDedup()
        # fix getDiaObjectsForDedup returning a non-unique pandas index
        all_diaObjects = all_diaObjects.reset_index().drop('index', axis=1)

        if len(all_diaObjects) > self.config.maxDiaObjects:
            raise RuntimeError(f"Found {len(all_diaObjects)} DiaObjects, which exceeds the "
                               f"configured maximum of {self.config.maxDiaObjects}. "
                               "Aborting to avoid excessive memory/computation.")
        self.log.info(f"Loaded {len(all_diaObjects)} DiaObjects.")

        # only keep the latest versions of the DIAObjects
        all_diaObjects.sort_values(['diaObjectId', 'validityStartMjdTai'], inplace=True)
        wdup = all_diaObjects.duplicated(subset='diaObjectId', keep='last')
        # copy to avoid setting on copy errors
        diaObjects = all_diaObjects.loc[~wdup, :].copy()
        self.log.info(f"Filtered to {len(all_diaObjects)} latest DiaObjects.")

        if len(diaObjects) == 0:
            raise pipeBase.NoWorkFound("No DiaObjects available for deduplication.")

        duplicate_count = self.count_duplicates(diaObjects)

        if duplicate_count == 0:
            raise pipeBase.NoWorkFound("No duplicate DiaObjects found.")
        else:
            self.log.info(f"Found {duplicate_count} duplicates in {len(diaObjects)} DiaObjects.")

        cluster_labels = self.cluster(diaObjects)

        diaObjects.loc[:, 'cluster_label'] = cluster_labels

        deduplication_map = self.remap_clusters(diaObjects)

        self.log.info(f"Reassigned {len(deduplication_map)} duplicates to "
                      f"{len(np.unique(deduplication_map['keptDiaObjectId']))} DiaObjects.")

        removed_ids = deduplication_map['removedDiaObjectId'].to_list()
        wremoved = diaObjects['diaObjectId'].apply(lambda x: x in removed_ids)
        kept_diaObjects = diaObjects.loc[~wremoved, :]
        updated_duplicate_count = self.count_duplicates(kept_diaObjects)
        self.log.info(f"After deduplication, {updated_duplicate_count} duplicates "
                      f"remain of {len(kept_diaObjects)} DiaObjects.")

        # TODO: consider a "dry run" configuration
        self.remove_apdb_duplicates(diaObjects, deduplication_map)

        return pipeBase.Struct(
            diaObjectDeduplicationMap=Table.from_df(deduplication_map))

    def remove_apdb_duplicates(self, diaObjects, deduplication_map):
        """Reassign DiaSources and remove DiaObjects per provided map.

        This method performs the actual database modifications to
        deduplicate DiaObjects. It reassigns DiaSources from duplicate
        DiaObjects to their corresponding kept objects, then marks the
        duplicate DiaObjects as invalid by setting their validity end
        time.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            DataFrame containing DiaObjects, including both kept and
            removed objects. Must include columns: 'diaObjectId',
            'validityStartMjdTai', 'ra', 'dec'.
        deduplication_map : `pandas.DataFrame`
            Table specifying how to reassociate duplicates. Must contain
            columns:
            - 'removedDiaObjectId' : IDs of DiaObjects to be removed
            - 'keptDiaObjectId' : IDs of DiaObjects to keep
        """

        start_time = Time(self.config.earliestMidpointMjdTai,
                          format='mjd', scale='tai')

        # make the records needed by the apdb methods
        DiaObjectsToRemove = []
        for idx, idi in deduplication_map['removedDiaObjectId'].items():
            w = diaObjects['diaObjectId'] == idi
            assert (np.sum(w) == 1)
            DiaObjectsToRemove.extend([daxApdb.recordIds.DiaObjectId.from_named_tuple(row)
                                       for row in diaObjects.loc[w].itertuples()])

        diaSourcesToReassign = \
            self.apdb.getDiaSourcesForDiaObjects(DiaObjectsToRemove,
                                                 start_time=start_time,
                                                 max_dist_arcsec=self.config.maxClusteringDistance)

        id_map = {}
        for old_id, new_id in zip(deduplication_map['removedDiaObjectId'],
                                  deduplication_map['keptDiaObjectId']):
            w = diaSourcesToReassign['diaObjectId'] == old_id
            reassign_count = np.sum(w)

            if reassign_count:
                wnew = diaObjects['diaObjectId'] == new_id
                assert (np.sum(wnew) == 1)
                newDiaObjectId = [daxApdb.recordIds.DiaObjectId.from_named_tuple(row)
                                  for row in diaObjects.loc[wnew].itertuples()][0]
                for row in diaSourcesToReassign.loc[w].itertuples():
                    ds = daxApdb.recordIds.DiaSourceId.from_named_tuple(row)
                    id_map[ds] = newDiaObjectId
                self.log.verbose('Reassigned %d diaSources from diaObject %d to %d' %
                                 (reassign_count, old_id, new_id))

            else:
                self.log.verbose('No diaSources found for diaObject %d' % old_id)

        self.apdb.reassignDiaSourcesToDiaObjects(id_map)

        self.apdb.setValidityEnd(DiaObjectsToRemove, Time.now())

        # Clear the table staging DiaObjects for deduplication
        self.apdb.resetDedup()

    def remap_clusters(self, diaObjects):
        """Use cluster labels to determine which DiaObjects to remap.

        For each cluster containing multiple DiaObjects, this method
        identifies the oldest object (by validityStartMjdTai) to keep
        and marks all others for removal.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            DataFrame of DiaObjects with cluster assignments. Must
            include columns:

            - 'cluster_label' : Cluster ID from clustering algorithm
            - 'diaObjectId' : Unique identifier for each DiaObject
            - 'validityStartMjdTai' : Start time of object validity

        Returns
        -------
        deduplication_map : `pandas.DataFrame`
            Mapping of removed DiaObject IDs to kept IDs. Contains
            columns:

            - 'removedDiaObjectId' : IDs to be removed
            - 'keptDiaObjectId' : IDs to keep

            Only includes entries for clusters with multiple objects.
        """

        grp = diaObjects.groupby('cluster_label')
        dedup_map = []
        for drp_id_i, group in grp:
            if len(group) > 1:
                ids = group['diaObjectId'].tolist()
                # keep the oldest one
                min_age_idx = group['validityStartMjdTai'].argmin()
                id_keep = ids.pop(min_age_idx)
                for idi in ids:
                    dedup_map.append((idi, id_keep))

        return pd.DataFrame(dedup_map, columns=['removedDiaObjectId', 'keptDiaObjectId'])

    def cluster(self, diaObjects):
        """Use agglomerative clustering to identify duplicate groups.

        This method applies hierarchical agglomerative clustering to
        DiaObjects based on their sky positions (RA, Dec). Objects
        within the configured maximum distance are grouped into
        clusters representing duplicate sets.

        For large datasets, the data is first partitioned using KMeans
        to create manageable subsets, then agglomerative clustering is
        applied within each subset. Objects near partition boundaries
        are processed with overlap to avoid missing cross-boundary
        duplicates.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            DataFrame of DiaObjects to cluster. Must include columns:

            - 'ra' : Right ascension in degrees
            - 'dec' : Declination in degrees

        Returns
        -------
        cluster_labels : `pandas.Series`
            Series of cluster labels (integers) indexed to match the
            input DataFrame. Objects with the same label belong to the
            same cluster and are considered duplicates.

        Notes
        -----
        A k-nearest neighbors connectivity graph is used to improve
        performance. However, it prevents objects from being grouped
        together if they are not among each other's k nearest
        neighbors, even if they fall within the desired distance
        threshold. `radius_neighbors_graph with a distance threshold
        would avoid this issue but has proven computationally
        infeasible.
        """
        n_objects = len(diaObjects)

        # If the dataset is small enough, use direct clustering
        if n_objects <= self.config.maxSubsetSize:
            return self._cluster_subset(diaObjects)

        # For large datasets, partition using KMeans first
        coords = diaObjects[['ra', 'dec']].values

        # Determine number of partitions
        n_partitions = max(2, int(np.ceil(n_objects / self.config.maxSubsetSize)))
        self.log.info(f"Partitioning {n_objects} objects into {n_partitions} subsets "
                      "for memory-efficient clustering.")

        # Use MiniBatchKMeans for memory efficiency
        kmeans = MiniBatchKMeans(n_clusters=n_partitions, random_state=42,
                                 batch_size=min(10000, n_objects),
                                 n_init=3)
        partition_labels = kmeans.fit_predict(coords)

        # Initialize cluster labels with unique values per object
        # We'll use a union-find approach to merge clusters across partitions
        cluster_labels = np.arange(n_objects)

        # Build a BallTree for efficient neighbor queries
        # Convert threshold to radians for haversine metric
        threshold_rad = np.radians(self.config.maxClusteringDistance / 3600.0)

        # Process each partition with overlap
        for partition_id in range(n_partitions):
            # Get objects in this partition
            partition_mask = partition_labels == partition_id
            partition_indices = np.where(partition_mask)[0]

            if len(partition_indices) == 0:
                continue

            # Find objects near partition boundaries by checking distance
            # to objects in other partitions
            partition_coords = coords[partition_indices]

            # Find nearby objects from other partitions that could be
            # duplicates of objects in this partition
            other_mask = ~partition_mask
            other_indices = np.where(other_mask)[0]

#            if len(other_indices) > 0:
            if False:
                other_coords = coords[other_indices]

                # Use BallTree to find objects within threshold distance
                tree = BallTree(np.radians(partition_coords), metric='haversine')
                # Query which other objects are within threshold of
                # partition objects
                nearby_indices_list = tree.query_radius(np.radians(other_coords),
                                                        r=threshold_rad)

                # Find which other objects have neighbors in this partition
                nearby_other_mask = np.array([len(x) > 0 for x in nearby_indices_list])
                nearby_other_indices = other_indices[nearby_other_mask]

                # Combine partition objects with nearby boundary objects
                combined_indices = np.concatenate([partition_indices, nearby_other_indices])
            else:
                combined_indices = partition_indices

            if len(combined_indices) == 0:
                continue

            # Create subset DataFrame for clustering
            subset_df = diaObjects.iloc[combined_indices]

            # Cluster this subset
            subset_labels = self._cluster_subset(subset_df)

            # Map subset cluster labels back to global indices.
            # Objects with the same subset label should have the same
            # global label.
            subset_indices = combined_indices
            label_to_global = {}

            for i, (idx, label) in enumerate(zip(subset_indices, subset_labels.values)):
                if label not in label_to_global:
                    # Use the first index we see for this label as the
                    # canonical one
                    label_to_global[label] = cluster_labels[idx]
                else:
                    # Merge: set this object's label to match the canonical one
                    old_label = cluster_labels[idx]
                    new_label = label_to_global[label]
                    if old_label != new_label:
                        # Update all objects with old_label to new_label
                        cluster_labels[cluster_labels == old_label] = new_label

        # Renumber labels to be contiguous
        unique_labels = np.unique(cluster_labels)
        label_map = {old: new for new, old in enumerate(unique_labels)}
        final_labels = np.array([label_map[lbl] for lbl in cluster_labels])

        self.log.info(f"Clustered {n_objects} diaObjects into "
                      f"{len(unique_labels)} clusters.")

        return pd.Series(final_labels, index=diaObjects.index,
                         name="cluster_label")

    def _cluster_subset(self, diaObjects):
        """Apply agglomerative clustering to a subset of DiaObjects.

        This is the core clustering routine applied to subsets that fit
        in memory.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            DataFrame of DiaObjects to cluster. Must include columns:

            - 'ra' : Right ascension in degrees
            - 'dec' : Declination in degrees

        Returns
        -------
        cluster_labels : `pandas.Series`
            Series of cluster labels (integers) indexed to match the
            input DataFrame.
        """
        coords = diaObjects[['ra', 'dec']].values

        # Adjust n_neighbors if subset is smaller than configured value
        n_neighbors = min(self.config.nNeighborsConnectivity, len(diaObjects) - 1)
        if n_neighbors < 1:
            # Single object, just return label 0
            return pd.Series([0], index=diaObjects.index, name="cluster_label")

        connectivity = kneighbors_graph(coords, n_neighbors=n_neighbors,
                                        mode='connectivity', include_self='auto')

        clustering = AgglomerativeClustering(distance_threshold=self.config.maxClusteringDistance/3600.,
                                             connectivity=connectivity,
                                             n_clusters=None).fit(coords)

        return pd.Series(clustering.labels_, index=diaObjects.index,
                         name="cluster_label")

    def count_duplicates(self, catalog, radius_arcsec=None):
        """Count catalog objects with self-matches within given radius.

        This method counts how many objects in the catalog have at least
        one other object within the search radius, providing a measure
        of the duplicate population.

        Parameters
        ----------
        catalog : `pandas.DataFrame`
            Catalog of objects to check for duplicates.
        radius_arcsec : `float`, optional
            Search radius in arcseconds within which to consider
            objects as duplicates. If None, defaults to half of
            config.maxClusteringDistance.

        Returns
        -------
        duplicate_count : `int`
            Number of objects that have at least one neighbor within
            the specified radius (excluding self-matches).

        Notes
        -----
        The method uses Astropy's catalog matching with nthneighbor=2
        to find the nearest neighbor excluding the object itself.
        Objects whose nearest neighbor is within radius_arcsec are
        counted as having duplicates.

        This provides a quick diagnostic count and may not exactly
        match the number of duplicates identified by the clustering
        algorithm.
        """

        if radius_arcsec is None:
            radius_arcsec = self.config.maxClusteringDistance/2

        sc = SkyCoord(catalog['ra'], catalog['dec'], unit='deg')
        idx, d2d, _ = sc.match_to_catalog_sky(sc, nthneighbor=2)
        wmatch = d2d < radius_arcsec*u.arcsecond
        return np.sum(wmatch)
