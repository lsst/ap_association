#
# LSST Data Management System
# Copyright 2017 LSST/AURA.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

"""A simple implementation of source association task for ap_verify.
"""

from __future__ import absolute_import, division, print_function

import numpy as np
from scipy.spatial import cKDTree

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
from .assoc_db_sqlite import AssociationDBSqliteTask

__all__ = ["AssociationConfig", "AssociationTask"]


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
    if len(dia_sources) == 1:
        dia_object_record['psFluxMean_%s' % filter_name] = dia_sources[0]['psFlux']
        dia_object_record['psFluxSigma_%s' % filter_name] = dia_sources[0]['psFluxErr']
        dia_object_record['psFluxMeanErr_%s' % filter_name] = dia_sources[0]['psFluxErr']
    else:
        fluxes = dia_sources.get("psFlux")[
            dia_sources.get('filterId') == filter_id]
        dia_object_record['psFluxMean_%s' % filter_name] = np.mean(fluxes)
        dia_object_record['psFluxSigma_%s' % filter_name] = np.std(fluxes, 1)
        dia_object_record['psFluxMeanErr_%s' % filter_name] = \
            dia_object_record['psFluxSigma_%s' % filter_name] / len(dia_sources)


class AssociationConfig(pexConfig.Config):
    """Config class for AssociationTask.
    """
    level1_db = pexConfig.ConfigurableField(
        target=AssociationDBSqliteTask,
        doc='Specify where and how to load and store DIAObjects and '
        'DIASources.',
    )
    maxDistArcSeconds = pexConfig.Field(
        dtype=float,
        doc='Maximum distance in arcseconds to test for a DIASource to be a '
        'match to a DIAObject.',
        default=1.0,
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
        self.makeSubtask('level1_db')

    @pipeBase.timeMethod
    def run(self, dia_sources, exposure):
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
        """
        # Assure we have a Box2D and can use the getCenter method.

        dia_objects = self.level1_db.load_dia_objects(exposure)

        updated_obj_ids = self.associate_sources(dia_objects, dia_sources)

        self.level1_db.store_ccd_visit_info(exposure)

        # Create aliases to appropriate flux fields
        dia_sources.getSchema.getAliasMap().set('psFlux', 'base_PsfFlux_flux')
        dia_sources.getSchema.getAliasMap().set('psFluxErr', 'base_PsfFlux_fluxSigma')

        # Store newly associated DIASources.
        self.level1_db.store_dia_sources(
            dia_sources, updated_obj_ids, exposure)
        # Update previously existing DIAObjects with the information from their
        # newly association DIASources and create new DIAObjects from
        # unassociated sources.
        self.update_dia_objects(updated_obj_ids)

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
        updated_ids: array-like of `int`s
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
    def update_dia_objects(self, updated_obj_ids, exposure):
        """Update select dia_objects currently stored within the database or
        create new ones.

        Parameters
        ----------
        updated_obj_ids : array-like of `int`s
            Ids of the dia_objects that should be updated.
        exposure : `lsst.afw.image.Exposure`
            Exposure the newly associated DIASources were detected in.
        """
        updated_dia_objects = afwTable.SourceCatalog(
            self.level1_db.get_dia_object_schema())
        updated_dia_objects.reserve(len(updated_obj_ids))

        filter_name = exposure.getFilter().getName()
        filter_id = self.level1_db.get_db_filter_id_from_name(filter_name)
        for obj_id in updated_obj_ids:
            dia_object = updated_dia_objects.addNew()
            dia_object.set('id', obj_id)

            dia_sources = self.level1_db.load_dia_sources([obj_id])

            dia_object.set('n_dia_sources', len(dia_sources))

            ave_coord = _set_mean_position(dia_object, dia_sources)
            indexer_id = self.level1_db.compute_indexer_id(ave_coord)
            dia_object.set('indexer_id', indexer_id)
            _set_flux_stats(dia_object,
                            dia_sources,
                            filter_name,
                            filter_id)

        self.level1_db.store_dia_objects(updated_dia_objects, False)

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
                (array-like of `float`s).
            - ``obj_idxs``: indexes of the matched DIAObjects in the catalog.
                (array-like of `int`s)
            - ``obj_ids``: array of floats of match quality updated DIAObjects
                (array-like of `int`s).
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
                updated DIAObjects (array-like of `float`s).
            - ``obj_ids``: array of floats of match quality
                updated DIAObjects (array-like of `ints`s).
            - ``obj_idxs``: indexes of the matched DIAObjects in the catalog.
                (array-like of `int`s)
            Default values for these arrays are
            INF, -1 and -1 respectively for unassociated sources.

        Returns
        -------
        result : `lsst.pipeBase.Struct`
            Results struct with components:

            - ``updated_and_new_dia_object_ids`` : ids of new and updated
              dia_objects as the result of association. (`list` of `int`s).
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
            associated_dia_object_ids[score_idx] = score_struct.obj_ids[score_idx]
            n_updated_dia_objects += 1

        # Argwhere returns a array shape (N, 1) so we access the index
        # thusly to retrieve the value rather than the tuple.
        for (src_idx,) in np.argwhere(np.logical_not(used_dia_source)):
            associated_dia_object_ids[src_idx] = dia_sources[int(src_idx)].getId()
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
              dia_objects in the collection (`list` of `int`s).
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
