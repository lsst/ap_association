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

"""Deduplicate DiaObjects exported by exportDiaCatalogsTask
"""


from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering

import lsst.dax.apdb as daxApdb
import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as connTypes
from lsst.skymap import BaseSkyMap
import lsst.sphgeom

from lsst.ap.association.utils import paddedRegion

__all__ = ("DeduplicateDiaObjectsTask", "DeduplicateDiaObjectsConfig")


class DeduplicateDiaObjectsConnections(pipeBase.PipelineTaskConnections,
                                       dimensions=("tract", "patch", "skymap")):

    skyMap = pipeBase.connectionTypes.Input(
        doc="Geometry of the tracts and patches that the coadds are defined on.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        dimensions=("skymap",),
        storageClass="SkyMap",
    )
    diaObjects = connTypes.Input(
        doc="DiaObjects preloaded from the APDB.",
        name="apdb_export_diaObjects",
        storageClass="ArrowAstropy",
        dimensions=("tract", "patch"),
    )
    diaObjectDeduplicationMap = connTypes.Output(
        doc="DiaSources preloaded from the APDB.",
        name="apdb_diaObject_deduplication_map",
        storageClass="ArrowAstropy",
        dimensions=("tract", "patch"),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)


class DeduplicateDiaObjectsConfig(pipeBase.PipelineTaskConfig,
                                  pipelineConnections=DeduplicateDiaObjectsConnections):
    """Config class for DeduplicateDiaObjectsConfig.
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
    angleMargin = pexConfig.RangeField(
        doc="Padding to add to the inner bbox (arcseconds)",
        dtype=float,
        default=2,
        min=0,
    )
    earliestMidpointMjdTai = pexConfig.Field(
        dtype=float,
        default=Time('2025-09-01', scale='tai').mjd,
        doc="MidpointMjdTai of earliest DIASource to be reassigned during "
            "deduplication.",
    )


class DeduplicateDiaObjectsTask(pipeBase.PipelineTask):
    """Deduplicate DiaObjects by clustering them spatially.
    """
    ConfigClass = DeduplicateDiaObjectsConfig
    _DefaultName = "deduplicateDiaObjects"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.apdb = daxApdb.Apdb.from_uri(self.config.apdb_config_url)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        diaObjects = inputs.pop("diaObjects")
        skymap = inputs.pop("skyMap")
        # TODO: this cannot be the best way to do this
        for _, ref in inputRefs:
            dataId = ref.dataId
            if 'tract' in dataId:
                break
            else:
                continue

        outputs = self.run(diaObjects, skymap, dataId["tract"], dataId["patch"])
        butlerQC.put(outputs, outputRefs)

    def run(self, diaObjects, skymap, tract, patch):
        """Take a set of DiaObjects and cluster them to remove duplicates.

        Parameters
        ----------

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components.
            - ``diaObjectDeduplicationMap`` : map of diaObjectIds to their
              replacements (`dict`)
        """

        trimmed_diaObjects = self._trimToPaddedBBox(diaObjects, skymap,
                                                    tract, patch)

        duplicate_count = self.count_duplicates(trimmed_diaObjects)

        if duplicate_count == 0:
            raise pipeBase.NoWorkFound("No duplicate DiaObjects found.")
        else:
            self.log.info(f"Found {duplicate_count} duplicates in {len(trimmed_diaObjects)} DiaObjects.")

        cluster_labels = self.cluster(trimmed_diaObjects)

        clustered_diaObjects = trimmed_diaObjects.join(cluster_labels)

        deduplication_map = self.remap_clusters(clustered_diaObjects)

        self.remove_apdb_duplicates(clustered_diaObjects, deduplication_map)

        return pipeBase.Struct(
            diaObjectDeduplicationMap=Table.from_df(deduplication_map))

    def remove_apdb_duplicates(self, diaObjects, deduplication_map):
        """Reassign diaSources and remove diaObjects according to provided map.

        Parameters
        ----------
        diaObjects : `pd.DataFrame`
        deduplication_map : `pd.DataFrame`
            Table with columns specifying how to reassociate duplicates.
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

        # for use with nightly deduplication
        # self.apdb.resetDedup()

    def remap_clusters(self, diaObjects):
        """Use cluster labels to determine which diaObjects to remap.
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
        """Use AgglomerativeClustering to identify groups of duplicates.
        """

        X = diaObjects[['ra', 'dec']].values
        clustering = AgglomerativeClustering(distance_threshold=self.config.maxClusteringDistance/3600.,
                                             n_clusters=None).fit(X)

        self.log.info(f"Clustered {len(clustering.labels_)} diaObjects into "
                      f"{len(np.unique(clustering.labels_))} clusters.")

        return pd.Series(clustering.labels_, index=diaObjects.index,
                         name="cluster_label")

    def count_duplicates(self, catalog, radius_arcsec=None):
        """Count catalog objects with self-matches within radius_arcsec.
        """

        if radius_arcsec is None:
            radius_arcsec = self.config.maxClusteringDistance/2

        sc = SkyCoord(catalog['ra'], catalog['dec'], unit='deg')
        idx, d2d, _ = sc.match_to_catalog_sky(sc, nthneighbor=2)
        wmatch = d2d < radius_arcsec*u.arcsecond
        return np.sum(wmatch)

    def _trimToPaddedBBox(self, diaObjects, skymap, tract, patch):
        """Trim catalog to the inner BBox padded by angleMargin.

        The APDB load in ExportDiaCatalogsTask returns a region larger than
        the inner patch BBox + angleMargin, so re-trim the data.

        Parameters
        ----------

        Returns
        -------
        trimmed_diaObjects: `pandas.DataFrame`
        """
        wcs = skymap[tract].getWcs()
        bbox = skymap[tract][patch].getInnerBBox()
        patchCorners = wcs.pixelToSky(lsst.geom.Box2D(bbox).getCorners())
        region = lsst.sphgeom.ConvexPolygon([pp.getVector() for pp in patchCorners])
        pad_region = paddedRegion(region,
                                  lsst.sphgeom.Angle.fromDegrees(self.config.angleMargin/3600.))

        diaObjectsDf = diaObjects.to_pandas()

        bbox_mask = []
        for ra, dec in zip(diaObjectsDf['ra'], diaObjectsDf['dec']):
            p = lsst.sphgeom.UnitVector3d(lsst.sphgeom.LonLat.fromDegrees(ra, dec))
            bbox_mask.append(pad_region.contains(p))

        self.log.info(f"Trimmed {len(bbox_mask)} diaObjects to {np.sum(bbox_mask)}.")

        return diaObjectsDf.iloc[bbox_mask]
