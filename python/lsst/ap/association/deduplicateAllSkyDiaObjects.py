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
from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import kneighbors_graph

import lsst.dax.apdb as daxApdb
import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as connTypes
from lsst.skymap import BaseSkyMap
import lsst.sphgeom

from lsst.ap.association.utils import paddedRegion

__all__ = ("DeduplicateAllSkyDiaObjectsTask", "DeduplicateAllSkyDiaObjectsConfig")


class DeduplicateAllSkyDiaObjectsConnections(pipeBase.PipelineTaskConnections,
                                       dimensions=("tract", "patch", "skymap")):
    skyMap = pipeBase.connectionTypes.Input(
        doc="Geometry of the tracts and patches that the coadds are defined on.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        dimensions=("skymap",),
        storageClass="SkyMap",
    )
    diaObjectDeduplicationMap = connTypes.Output(
        doc="DiaSources preloaded from the APDB.",
        name="apdb_diaObject_deduplication_map",
        storageClass="ArrowAstropy",
        dimensions=("tract", "patch"),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)


class DeduplicateAllSkyDiaObjectsConfig(pexConfig.Config):
    """Configuration for DeduplicateAllSkyDiaObjectsTask.
    
    Parameters
    ----------
    maxClusteringDistance : `float`
        Maximum distance to merge clusters of duplicate DiaObjects, in arcseconds.
        Default is 1.5 arcseconds.
    earliestMidpointMjdTai : `float`
        MidpointMjdTai (Modified Julian Date in TAI scale) of the earliest
        DiaSource to be reassigned during deduplication. Default corresponds
        to 2025-09-01.
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


class DeduplicateAllSkyDiaObjectsTask(pipeBase.Task):
    """Deduplicate DiaObjects by clustering them spatially and update the APDB.
    
    This task identifies and merges duplicate DiaObjects in the Alert Production
    Database (APDB) that are spatially coincident. It uses agglomerative clustering
    to group nearby objects, then reassigns DiaSources from duplicate objects to
    a single kept DiaObject (the oldest one in each cluster).
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
                Table mapping duplicate DiaObject IDs to the IDs they should be
                merged into. Contains columns 'removedDiaObjectId' and
                'keptDiaObjectId'.
        
        Raises
        ------
        lsst.pipe.base.NoWorkFound
            If no duplicate DiaObjects are found in the APDB.
        """

        # TODO: consider adding a config for the "since" argument
        all_diaObjects = self.apdb.getDiaObjectsForDedup()
        # fix getDiaObjectsForDedup returning a non-unique pandas index
        all_diaObjects = all_diaObjects.reset_index().drop('index', axis=1)
        self.log.info(f"Loaded {len(all_diaObjects)} DiaObjects.")

        # only keep the latest versions of the DIAObjects
        all_diaObjects.sort_values(['diaObjectId', 'validityStartMjdTai'], inplace=True)
        wdup = all_diaObjects.duplicated(subset='diaObjectId', keep='last')
        # copy to avoid setting on copy errors
        diaObjects = all_diaObjects.loc[~wdup, :].copy()
        self.log.info(f"Filtered to {len(all_diaObjects)} latest DiaObjects.")

        duplicate_count = self.count_duplicates(diaObjects)

        if duplicate_count == 0:
            raise pipeBase.NoWorkFound("No duplicate DiaObjects found.")
        else:
            self.log.info(f"Found {duplicate_count} duplicates in {len(diaObjects)} DiaObjects.")

        cluster_labels = self.cluster(diaObjects)

        diaObjects.loc[:,'cluster_label'] = cluster_labels

        deduplication_map = self.remap_clusters(diaObjects)

        self.log.info(f"Reassigned {len(deduplication_map)} duplicates to "
                       "{len(np.unique(deduplication_map['keptDiaObjectId']))} DiaObjects.")

        removed_ids = deduplication_map['removedDiaObjectId'].to_list()
        wremoved = diaObjects['diaObjectId'].apply(lambda x: x in removed_ids)
        kept_diaObjects = diaObjects.loc[~wremoved,:]
        updated_duplicate_count = self.count_duplicates(kept_diaObjects)
        self.log.info(f"After deduplication, {updated_duplicate_count} duplicates "
                       "remain of {len(kept_diaObjects)} DiaObjects.")

        # TODO: consider a "dry run" configuration
        # self.remove_apdb_duplicates(clustered_diaObjects, deduplication_map)

        return pipeBase.Struct(
            diaObjectDeduplicationMap=Table.from_df(deduplication_map))

    def remove_apdb_duplicates(self, diaObjects, deduplication_map):
        """Reassign DiaSources and remove DiaObjects according to provided map.
        
        This method performs the actual database modifications to deduplicate
        DiaObjects. It reassigns DiaSources from duplicate DiaObjects to their
        corresponding kept objects, then marks the duplicate DiaObjects as
        invalid by setting their validity end time.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            DataFrame containing DiaObjects, including both kept and removed
            objects. Must include columns: 'diaObjectId', 'validityStartMjdTai',
            'ra', 'dec'.
        deduplication_map : `pandas.DataFrame`
            Table specifying how to reassociate duplicates. Must contain columns:
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
        for old_id, new_id in deduplication_map.iterrows():
            w = diaSourcesToReassign['diaObjectId'] == old_id
            reassign_count = np.sum(w)

            if reassign_count:
                for row in diaSourcesToReassign.loc[w].itertuples():
                    ds = daxApdb.recordIds.DiaSourceId.from_named_tuple(row)
                    id_map[ds] = new_id
                self.log.verbose('Reassigned {reassign_count} diaSources from diaObject {old_id} to {new_id}')

            else:
                self.log.verbose('No diaSources found for diaObject {old_id}')

        self.apdb.reassignDiaSourcesToDiaObjects(id_map)

        self.apdb.setValidityEnd(DiaObjectsToRemove, Time.now())

        # TODO: turn on
        # self.apdb.resetDedup()

    def remap_clusters(self, diaObjects):
        """Use cluster labels to determine which DiaObjects to remap.
        
        For each cluster containing multiple DiaObjects, this method identifies
        the oldest object (by validityStartMjdTai) to keep and marks all others
        for removal.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            DataFrame of DiaObjects with cluster assignments. Must include columns:
            
            - 'cluster_label' : Cluster ID assigned by the clustering algorithm
            - 'diaObjectId' : Unique identifier for each DiaObject
            - 'validityStartMjdTai' : Start time of object validity 

        Returns
        -------
        deduplication_map : `pandas.DataFrame`
            Mapping of removed DiaObject IDs to kept IDs. Contains columns:
            
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
                    if idi == id_keep:
                        import ipdb; ipdb.set_trace()
                    dedup_map.append((idi, id_keep))

        return pd.DataFrame(dedup_map, columns=['removedDiaObjectId', 'keptDiaObjectId'])

    def cluster(self, diaObjects):
        """Use agglomerative clustering to identify groups of duplicates.
        
        This method applies hierarchical agglomerative clustering to DiaObjects
        based on their sky positions (RA, Dec). Objects within the configured
        maximum distance are grouped into clusters representing duplicate sets.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            DataFrame of DiaObjects to cluster. Must include columns:
            
            - 'ra' : Right ascension in degrees
            - 'dec' : Declination in degrees

        Returns
        -------
        cluster_labels : `pandas.Series`
            Series of cluster labels (integers) indexed to match the input
            DataFrame. Objects with the same label belong to the same cluster
            and are considered duplicates. 
        
        Notes
        -----
        A k-nearest neighbors connectivity graph is used to improve performance.
        However, it prevents objects from being grouped together if they are not among
        each other's k nearest neighbors, even if they fall within the desired distance
        threshold.  `radius_neighbors_graph with a distance threshold would 
        avoid this issue but has proven computationally infeasible.
        """

        X = diaObjects[['ra', 'dec']].values
        connectivity = kneighbors_graph(X, n_neighbors=30, 
                                        mode='connectivity', include_self='auto')

        clustering = AgglomerativeClustering(distance_threshold=self.config.maxClusteringDistance/3600.,
                                             connectivity=connectivity,
                                             n_clusters=None).fit(X)

        self.log.info(f"Clustered {len(clustering.labels_)} diaObjects into "
                      f"{len(np.unique(clustering.labels_))} clusters.")

        return pd.Series(clustering.labels_, index=diaObjects.index,
                         name="cluster_label")

    def count_duplicates(self, catalog, radius_arcsec=None):
        """Count catalog objects with self-matches within specified radius.
        
        This method counts how many objects in the catalog have at least one
        other object within the search radius, providing a measure of the
        duplicate population.

        Parameters
        ----------
        catalog : `pandas.DataFrame`
            Catalog of objects to check for duplicates. 
        radius_arcsec : `float`, optional
            Search radius in arcseconds within which to consider objects as
            duplicates. If None, defaults to half of config.maxClusteringDistance.

        Returns
        -------
        duplicate_count : `int`
            Number of objects that have at least one neighbor within the
            specified radius (excluding self-matches).
        
        Notes
        -----
        The method uses Astropy's catalog matching with nthneighbor=2 to find
        the nearest neighbor excluding the object itself. Objects whose nearest
        neighbor is within radius_arcsec are counted as having duplicates.
        
        This provides a quick diagnostic count and may not exactly match the
        number of duplicates identified by the clustering algorithm.
        """

        if radius_arcsec is None:
            radius_arcsec = self.config.maxClusteringDistance/2

        sc = SkyCoord(catalog['ra'], catalog['dec'], unit='deg')
        idx, d2d, _ = sc.match_to_catalog_sky(sc, nthneighbor=2)
        wmatch = d2d < radius_arcsec*u.arcsecond
        return np.sum(wmatch)
