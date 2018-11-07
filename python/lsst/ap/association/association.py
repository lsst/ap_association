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

"""A simple implementation of source association task for ap_verify.
"""

__all__ = ["AssociationConfig", "AssociationTask"]

import numpy as np
from scipy.spatial import cKDTree

from lsst.meas.algorithms.indexerRegistry import IndexerRegistry
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
from lsst.daf.base import DateTime

from .afwUtils import make_dia_object_schema


def _set_mean_position(dia_object_record, dia_sources):
    """Compute and set the mean position of the input dia_object_record using
    the positions of the input catalog of DIASources.

    Parameters
    ----------
    dia_object_record : `lsst.afw.table.SourceRecord`
        SourceRecord of the DIAObject to edit.
    dia_sources : `lsst.afw.table.SourceCatalog`
        Catalog of DIASources to compute a mean position from.

    Returns
    -------
    ave_coord : `lsst.afw.geom.SpherePoint`
        Average position of the dia_sources.
    """
    coord_list = [src.getCoord() for src in dia_sources]
    ave_coord = afwGeom.averageSpherePoint(coord_list)
    dia_object_record.setCoord(ave_coord)

    return ave_coord


def _set_flux_stats(dia_object_record, dia_sources, filter_name, filter_id):
    """Compute the mean, standard error, and variance of a DIAObject for
    a given band.

    Parameters
    ----------
    dia_object_record : `lsst.afw.table.SourceRecord`
        SourceRecord of the DIAObject to edit.
    dia_sources : `lsst.afw.table.SourceCatalog`
        Catalog of DIASources to compute a mean position from.
    filter_name : `str`
        Name of the band pass filter to update.
    filter_id : `int`
        id of the filter in the AssociationDB.
    """
    n_sources = len(dia_sources)
    if n_sources == 1 and dia_sources[0]["filterId"] == filter_id:
        dia_object_record['%sPSFluxMean' % filter_name] = dia_sources[0]['psFlux']
        dia_object_record['%sPSFluxSigma' % filter_name] = np.nan
        dia_object_record['%sPSFluxMeanErr' % filter_name] = np.nan
        dia_object_record["%sPSFluxChi2" % filter_name] = np.nan
        dia_object_record['%sPSFluxNdata' % filter_name] = n_sources
    else:
        fluxes = dia_sources.get("psFlux")[
            dia_sources.get("filterId") == filter_id]
        fluxErrors = dia_sources.get("psFluxErr")[
            dia_sources.get("filterId") == filter_id]
        fluxMean = np.mean(fluxes)
        dia_object_record['%sPSFluxMean' % filter_name] = np.mean(fluxes)
        dia_object_record['%sPSFluxSigma' % filter_name] = np.std(fluxes, ddof=1)
        dia_object_record['%sPSFluxMeanErr' % filter_name] = \
            dia_object_record['%sPSFluxSigma' % filter_name] / len(fluxes)
        dia_object_record["%sPSFluxChi2" % filter_name] = np.sum(
            ((fluxMean - fluxes) / fluxErrors) ** 2)
        dia_object_record['%sPSFluxNdata' % filter_name] = len(fluxes)


class AssociationConfig(pexConfig.Config):
    """Config class for AssociationTask.
    """
    maxDistArcSeconds = pexConfig.Field(
        dtype=float,
        doc='Maximum distance in arcseconds to test for a DIASource to be a '
        'match to a DIAObject.',
        default=1.0,
    )
    indexer = IndexerRegistry.makeField(
        doc='Select the spatial indexer to use within the database.',
        default='HTM'
    )


class AssociationTask(pipeBase.Task):
    """Associate DIAOSources into existing DIAObjects.

    This task performs the association of detected DIASources in a visit
    with the previous DIAObjects detected over time. It also creates new
    DIAObjects out of DIASources that cannot be associated with previously
    detected DIAObjects.
    """

    ConfigClass = AssociationConfig
    _DefaultName = "association"

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.indexer = IndexerRegistry[self.config.indexer.name](
            self.config.indexer.active)
        self.dia_object_schema = make_dia_object_schema()

    @pipeBase.timeMethod
    def run(self, dia_sources, exposure, ppdb):
        """Load DIAObjects from the database, associate the sources, and
        persist the results into the L1 database.

        Parameters
        ----------
        dia_sources : `lsst.afw.table.SourceCatalog`
            DIASources to be associated with existing DIAObjects.
        exposure : `lsst.afw.image`
            Input exposure representing the region of the sky the dia_sources
            were detected on. Should contain both the solved WCS and a bounding
            box of the ccd.
        ppdb : `lsst.dax.ppdb.Ppdb`
            Ppdb connection object to retrieve DIASources/Objects from and
            write to.
        """
        # Assure we have a Box2D and can use the getCenter method.

        dia_objects = self.retrieve_dia_objects(exposure, ppdb)

        updated_obj_ids = self.associate_sources(dia_objects, dia_sources)

        # Store newly associated DIASources.
        visitInfo = exposure.getInfo().getVisitInfo()
        ppdb.saveVisit(visitInfo.getExposureId(),
                       visitInfo.getDate().toPython())
        ppdb.storeDiaSources(dia_sources)
        # Update previously existing DIAObjects with the information from their
        # newly association DIASources and create new DIAObjects from
        # unassociated sources.
        self.update_dia_objects(dia_objects,
                                updated_obj_ids,
                                exposure,
                                ppdb)

    @pipeBase.timeMethod
    def retrieve_dia_objects(self, exposure, ppdb):
        """Convert the exposure object into HTM pixels and retrieve DIAObjects
        contained within the exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            An exposure specifying a bounding region with a WCS to load
            DIAOjbects within.
        ppdb : `lsst.dax.ppdb.Ppdb`
            Ppdb connection object to retrieve DIAObjects from.

        Returns
        -------
        diaObjects : `lsst.afw.table.SourceCatalog`
        """
        bbox = afwGeom.Box2D(exposure.getBBox())
        wcs = exposure.getWcs()

        ctr_coord = wcs.pixelToSky(bbox.getCenter())
        max_radius = max(
            ctr_coord.separation(wcs.pixelToSky(pp))
            for pp in bbox.getCorners())

        indexer_indices, on_boundry = self.indexer.getShardIds(
            ctr_coord, max_radius)
        # Index types must be cast to int to work with dax_ppdb.
        index_ranges = [[int(indexer_idx), int(indexer_idx) + 1]
                        for indexer_idx in indexer_indices]
        covering_dia_objects = ppdb.getDiaObjects(index_ranges)

        output_dia_objects = afwTable.SourceCatalog(
            covering_dia_objects.getSchema())
        for cov_dia_object in covering_dia_objects:
            if self._check_dia_object_position(cov_dia_object, bbox, wcs):
                output_dia_objects.append(cov_dia_object)

        # Return deep copy to enforce contiguity.
        return output_dia_objects.copy(deep=True)

    def _check_dia_object_position(self, dia_object_record, bbox, wcs):
        """Check the RA, DEC position of the current dia_object_record against
        the bounding box of the exposure.

        Parameters
        ----------
        dia_object_record : `lsst.afw.table.SourceRecord`
            A SourceRecord object containing the DIAObject we would like to
            test against our bounding box.
        bbox : `lsst.geom.Box2D`
            Bounding box of exposure.
        wcs : `lsst.geom.SkyWcs`
            WCS of exposure.

        Return
        ------
        is_contained : `bool`
            Object position is contained within the bounding box of expMd.
        """
        point = wcs.skyToPixel(dia_object_record.getCoord())
        return bbox.contains(point)

    @pipeBase.timeMethod
    def associate_sources(self, dia_objects, dia_sources):
        """Associate the input DIASources with the catalog of DIAObjects.

        Parameters
        ----------
        dia_objects : `lsst.afw.table.SourceCatalog`
            Catalog of DIAObjects to attempt to associate the input
            DIASources into.
        dia_sources : `lsst.afw.table.SourceCatalog`
            DIASources to associate into the DIAObjectCollection.

        Returns
        -------
        updated_ids : array-like of `int`
            Ids of the DIAObjects that the DIASources associated to including
            the ids of newly created DIAObjects.
        """

        scores = self.score(
            dia_objects, dia_sources,
            self.config.maxDistArcSeconds * afwGeom.arcseconds)
        match_result = self.match(dia_objects, dia_sources, scores)

        self._add_association_meta_data(match_result)

        return match_result.associated_dia_object_ids

    @pipeBase.timeMethod
    def update_dia_objects(self, dia_objects, updated_obj_ids, exposure, ppdb):
        """Update select dia_objects currently stored within the database or
        create new ones.

        Parameters
        ----------
        dia_objects : `lsst.afw.table.SourceCatalog`
            Pre-existing/loaded DIAObjects to copy values that are not updated
            from.
        updated_obj_ids : array-like of `int`
            Ids of the dia_objects that should be updated.
        exposure : `lsst.afw.image.Exposure`
            Input exposure representing the region of the sky the dia_sources
            were detected on. Should contain both the solved WCS and a bounding
            box of the ccd.
        ppdb : `lsst.dax.ppdb.Ppdb`
            Ppdb connection object to retrieve DIASources from and
            write DIAObjects to.
        """
        updated_dia_objects = afwTable.SourceCatalog(
            self.dia_object_schema)
        mapper = afwTable.SchemaMapper(dia_objects.schema, updated_dia_objects.schema)
        mapper.addMapping(dia_objects.schema.find("id")[0], "id", True)

        filter_name = exposure.getFilter().getName()
        filter_id = exposure.getFilter().getId()

        dateTime = exposure.getInfo().getVisitInfo().getDate()

        dia_sources = ppdb.getDiaSources(updated_obj_ids, dateTime.toPython())
        diaObjectId_key = dia_sources.schema['diaObjectId'].asKey()
        dia_sources.sort(diaObjectId_key)

        for obj_id in updated_obj_ids:
            dia_object = dia_objects.find(obj_id)
            if dia_object is None:
                dia_object = updated_dia_objects.addNew()
                dia_object.set('id', obj_id)
            else:
                updated_dia_objects.append(dia_object)

            # Select the dia_sources associated with this DIAObject id and
            # copy the copy the subcatalog for fast slicing.
            start_idx = dia_sources.lower_bound(obj_id, diaObjectId_key)
            end_idx = dia_sources.upper_bound(obj_id, diaObjectId_key)
            obj_dia_sources = dia_sources[start_idx:end_idx].copy(deep=True)

            ave_coord = _set_mean_position(dia_object, obj_dia_sources)
            indexer_id = self.indexer.indexPoints(
                [ave_coord.getRa().asDegrees()],
                [ave_coord.getDec().asDegrees()])[0]
            dia_object.set('pixelId', indexer_id)
            dia_object.set("radecTai", dateTime.get(system=DateTime.MJD))
            dia_object.set("nDiaSources", len(obj_dia_sources))
            _set_flux_stats(dia_object,
                            obj_dia_sources,
                            filter_name,
                            filter_id)

        ppdb.storeDiaObjects(updated_dia_objects, dateTime.toPython())

    @pipeBase.timeMethod
    def score(self, dia_objects, dia_sources, max_dist):
        """Compute a quality score for each dia_source/dia_object pair
        between this catalog of DIAObjects and the input DIASource catalog.

        ``max_dist`` sets maximum separation in arcseconds to consider a
        dia_source a possible match to a dia_object. If the pair is
        beyond this distance no score is computed.

        Parameters
        ----------
        dia_objects : `lsst.afw.table.SourceCatalog`
            A contiguous catalog of DIAObjects to score against dia_sources.
        dia_sources : `lsst.afw.table.SourceCatalog`
            A contiguous catalog of dia_sources to "score" based on distance
            and (in the future) other metrics.
        max_dist : `lsst.afw.geom.Angle`
            Maximum allowed distance to compute a score for a given DIAObject
            DIASource pair.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components:

            - ``scores``: array of floats of match quality updated DIAObjects
                (array-like of `float`).
            - ``obj_idxs``: indexes of the matched DIAObjects in the catalog.
                (array-like of `int`)
            - ``obj_ids``: array of floats of match quality updated DIAObjects
                (array-like of `int`).

            Default values for these arrays are
            INF, -1, and -1 respectively for unassociated sources.
        """
        scores = np.full(len(dia_sources), np.inf, dtype=np.float64)
        obj_idxs = np.full(len(dia_sources), -1, dtype=np.int)
        obj_ids = np.full(len(dia_sources), -1, dtype=np.int)

        if len(dia_objects) == 0:
            return pipeBase.Struct(
                scores=scores,
                obj_idxs=obj_idxs,
                obj_ids=obj_ids)

        spatial_tree = self._make_spatial_tree(dia_objects)

        max_dist_rad = max_dist.asRadians()

        for src_idx, dia_source in enumerate(dia_sources):

            src_point = dia_source.getCoord().getVector()
            dist, obj_idx = spatial_tree.query(src_point)
            if dist < max_dist_rad:
                scores[src_idx] = dist
                obj_ids[src_idx] = dia_objects[obj_idx].getId()
                obj_idxs[src_idx] = obj_idx

        return pipeBase.Struct(
            scores=scores,
            obj_idxs=obj_idxs,
            obj_ids=obj_ids)

    def _make_spatial_tree(self, dia_objects):
        """Create a searchable kd-tree the input dia_object positions.

        Parameters
        ----------
        dia_objects : `lsst.afw.table.SourceCatalog`
            A catalog of DIAObjects to create the tree from.

        Returns
        -------
        kd_tree : `scipy.spatical.cKDTree`
            Searchable kd-tree created from the positions of the DIAObjects.
        """

        coord_key = dia_objects.getCoordKey()
        coord_vects = np.empty((len(dia_objects), 3))

        for obj_idx, dia_object in enumerate(dia_objects):
            coord_vects[obj_idx] = dia_object[coord_key].getVector()

        return cKDTree(coord_vects)

    @pipeBase.timeMethod
    def match(self, dia_objects, dia_sources, score_struct):
        """Match DIAsources to DIAObjects given a score and create new
        DIAObject Ids for new unassociated DIASources.

        Parameters
        ----------
        dia_objects : `lsst.afw.table.SourceCatalog`
            A SourceCatalog of DIAObjects to associate to DIASources.
        dia_sources : `lsst.afw.table.SourceCatalog`
            A contiguous catalog of dia_sources for which the set of scores
            has been computed on with DIAObjectCollection.score.
        score_struct : `lsst.pipe.base.Struct`
            Results struct with components:

            - ``scores``: array of floats of match quality
                updated DIAObjects (array-like of `float`).
            - ``obj_ids``: array of floats of match quality
                updated DIAObjects (array-like of `int`).
            - ``obj_idxs``: indexes of the matched DIAObjects in the catalog.
                (array-like of `int`)

            Default values for these arrays are
            INF, -1 and -1 respectively for unassociated sources.

        Returns
        -------
        result : `lsst.pipeBase.Struct`
            Results struct with components:

            - ``updated_and_new_dia_object_ids`` : ids of new and updated
              dia_objects as the result of association. (`list` of `int`).
            - ``n_updated_dia_objects`` : Number of previously know dia_objects
              with newly associated DIASources. (`int`).
            - ``n_new_dia_objects`` : Number of newly created DIAObjects from
              unassociated DIASources (`int`).
            - ``n_unupdated_dia_objects`` : Number of previous DIAObjects that
              were not associated to a new DIASource (`int`).
        """

        n_previous_dia_objects = len(dia_objects)
        used_dia_object = np.zeros(n_previous_dia_objects, dtype=np.bool)
        used_dia_source = np.zeros(len(dia_sources), dtype=np.bool)
        associated_dia_object_ids = np.zeros(len(dia_sources),
                                             dtype=np.uint64)

        n_updated_dia_objects = 0
        n_new_dia_objects = 0

        # We sort from best match to worst to effectively perform a
        # "handshake" match where both the DIASources and DIAObjects agree
        # their the best match. By sorting this way, scores with NaN (those
        # sources that have no match and will create new DIAObjects) will be
        # placed at the end of the array.
        score_args = score_struct.scores.argsort(axis=None)
        for score_idx in score_args:
            if not np.isfinite(score_struct.scores[score_idx]):
                # Thanks to the sorting the rest of the sources will be
                # NaN for their score. We therefore exit the loop to append
                # sources to a existing DIAObject, leaving these for
                # the loop creating new objects.
                break
            dia_obj_idx = score_struct.obj_idxs[score_idx]
            if used_dia_object[dia_obj_idx]:
                continue
            used_dia_object[dia_obj_idx] = True
            used_dia_source[score_idx] = True
            obj_id = score_struct.obj_ids[score_idx]
            associated_dia_object_ids[score_idx] = obj_id
            n_updated_dia_objects += 1
            dia_sources[int(score_idx)].set("diaObjectId", int(obj_id))

        # Argwhere returns a array shape (N, 1) so we access the index
        # thusly to retrieve the value rather than the tuple
        for (src_idx,) in np.argwhere(np.logical_not(used_dia_source)):
            src_id = dia_sources[int(src_idx)].getId()
            associated_dia_object_ids[src_idx] = src_id
            dia_sources[int(src_idx)].set("diaObjectId", int(src_id))
            n_new_dia_objects += 1

        # Return the ids of the DIAObjects in this DIAObjectCollection that
        # were updated or newly created.
        n_unassociated_dia_objects = \
            n_previous_dia_objects - n_updated_dia_objects
        return pipeBase.Struct(
            associated_dia_object_ids=associated_dia_object_ids,
            n_updated_dia_objects=n_updated_dia_objects,
            n_new_dia_objects=n_new_dia_objects,
            n_unassociated_dia_objects=n_unassociated_dia_objects,)

    def _add_association_meta_data(self, match_result):
        """Store summaries of the association step in the task metadata.

        Parameters
        ----------
        match_result : `lsst.pipeBase.Struct`
            Results struct with components:

            - ``updated_and_new_dia_object_ids`` : ids new and updated
              dia_objects in the collection (`list` of `int`).
            - ``n_updated_dia_objects`` : Number of previously know dia_objects
              with newly associated DIASources. (`int`).
            - ``n_new_dia_objects`` : Number of newly created DIAObjects from
              unassociated DIASources (`int`).
            - ``n_unupdated_dia_objects`` : Number of previous DIAObjects that
              were not associated to a new DIASource (`int`).
        """
        self.metadata.add('numUpdatedDiaObjects',
                          match_result.n_updated_dia_objects)
        self.metadata.add('numNewDiaObjects',
                          match_result.n_new_dia_objects)
        self.metadata.add('numUnassociatedDiaObjects',
                          match_result.n_unassociated_dia_objects)
