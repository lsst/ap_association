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

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
from .assoc_db_sqlite import AssociationDBSqliteTask
from .dia_object import DIAObject

__all__ = ["AssociationConfig", "AssociationTask"]


class AssociationConfig(pexConfig.Config):
<<<<<<< HEAD
    """Config class for AssociationTask.
=======
    """
    Config class for AssociationTask.
>>>>>>> ae74366... Debug Sphnix docs.
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
<<<<<<< HEAD
    """Associate DIAOSources into existing DIAObjects.
=======
    """
    Associate DIAOSources into existing DIAObjects.
>>>>>>> ae74366... Debug Sphnix docs.

    This task performs the association of detected DIASources in a visit
    with the previous DIAObjects detected over time. It also creates new
    DIAObjects out of DIASources that cannot be associated with previously
    detected DIAObjects.
    """

    ConfigClass = AssociationConfig
    _DefaultName = "association"

    def __init__(self, **kwargs):
        """ Initialize the the association task and create the database link.
        """
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
        bbox = afwGeom.Box2D(exposure.getBBox())
        wcs = exposure.getWcs()
        expMd = pipeBase.Struct(
            bbox=bbox,
            wcs=wcs,)

        dia_collection = self.level1_db.load(expMd)

        association_result = self.associate_sources(dia_collection,
                                                    dia_sources)

        dia_collection = association_result.dia_collection
        updated_obj_ids = association_result.updated_ids

        dia_collection.update_dia_objects()

        self.level1_db.store_updated(dia_collection, updated_obj_ids)

    @pipeBase.timeMethod
    def associate_sources(self, dia_collection, dia_sources):
        """Associate the input DIASources in to the collection of DIAObjects.

        Parameters
        ----------
        dia_collection : `lsst.ap.association.DIAObjectCollection`
            Collection of DIAObjects to attempt to associate the input
            DIASources into.
        dia_sources : `lsst.afw.table.SourceCatalog`
            DIASources to associate into the DIAObjectCollection.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components:

            - ``dia_collection``: A DIAObjectCollectoin containing the new and
              updated DIAObjects (`lsst.ap.association.DIACollection`).
            - ``updated_ids``: id of the DIAObject in this DIAObjectCollection
              that the given source matched. (`list` of `int`s).
        """

        scores = self.score(
            dia_collection, dia_sources,
            self.config.maxDistArcSeconds * afwGeom.arcseconds)
        match_result = self.match(dia_collection, dia_sources, scores)

        self._add_association_meta_data(match_result)

        return pipeBase.Struct(
            dia_collection=dia_collection,
            updated_ids=match_result.updated_and_new_dia_object_ids,
        )

    def score(self, dia_collection, dia_source_catalog, max_dist):
        """Compute a quality score for each dia_source/dia_object pair
        between this collection and an input diat_source catalog.

        max_dist sets maximum separation in arcseconds to consider a
        dia_source a possible match to a dia_object. If the pair is
        beyond this distance no score is computed.

        Parameters
        ----------
<<<<<<< HEAD
        dia_object_collection : `lsst.ap.association.DIAObjectCollection`
            A DIAObjectCollection to score against dia_sources.
        dia_source_catalog : `lsst.afw.table.SourceCatalog`
=======
        dia_object_collection : lsst.ap.association.DIAObjectCollection
            A DIAObjectCollection to score against dia_sources.
        dia_source_catalog : lsst.afw.SourceCatalog
>>>>>>> ae74366... Debug Sphnix docs.
            A contiguous catalog of dia_sources to "score" based on distance
            and (in the future) other metrics.
        max_dist : `lsst.afw.geom.Angle`
            Maximum allowed distance to compute a score for a given DIAObject
            DIASource pair.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components:

            - ``scores``: array of floats of match quality
                updated DIAObjects (`ndarray` of `float`s).
            - ``obj_ids``: array of floats of match quality
                updated DIAObjects (`ndarray` of `ints`s).
            Default values for these arrays are
            INF and -1 respectively for unassociated sources.
        """
        if not dia_collection._is_valid_tree:
            dia_collection.update_spatial_tree()

        scores = np.ones(len(dia_source_catalog)) * np.inf
        obj_ids = -1 * np.ones(len(dia_source_catalog), dtype=np.int)

        if len(dia_collection.dia_objects) == 0:
            return pipeBase.Struct(
                scores=scores,
                obj_ids=obj_ids)

        for src_idx, dia_source in enumerate(dia_source_catalog):

            src_point = dia_source.getCoord().getVector()
            dist, obj_idx = dia_collection._spatial_tree.query(src_point)
            if dist < max_dist.asRadians():
                scores[src_idx] = dist
                obj_ids[src_idx] = dia_collection.dia_objects[obj_idx].id

        return pipeBase.Struct(
            scores=scores,
            obj_ids=obj_ids)

    def match(self, dia_collection, dia_source_catalog, score_struct):
        """Append DIAsources to DIAObjects given a score and create new
        DIAObjects in this collection from DIASources with poor scores.

        Parameters
        ----------
<<<<<<< HEAD
        dia_object_collection : `lsst.ap.association.DIAObjectCollection`
            A DIAObjectCollection to associate to dia_sources.
        dia_source_catalog : `lsst.afw.table.SourceCatalog`
=======
        dia_object_collection : lsst.ap.association.DIAObjectCollection
            A DIAObjectCollection to associate to dia_sources.
        dia_source_catalog : lsst.afw.SourceCatalog
>>>>>>> ae74366... Debug Sphnix docs.
            A contiguous catalog of dia_sources for which the set of scores
            has been computed on with DIAObjectCollection.score.
        score_struct : `lsst.pipe.base.Struct`
            Results struct with components:

            - ``scores``: array of floats of match quality
                updated DIAObjects (`ndarray` of `float`s).
            - ``obj_ids``: array of floats of match quality
                updated DIAObjects (`ndarray` of `ints`s).
            Default values for these arrays are
            INF and -1 respectively for unassociated sources.

        Returns
        -------
        result : `lsst.pipeBase.Struct`
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

        n_previous_dia_objects = len(dia_collection.dia_objects)
        used_dia_object = np.zeros(n_previous_dia_objects, dtype=np.bool)
        used_dia_source = np.zeros(len(dia_source_catalog), dtype=np.bool)

        updated_dia_objects = []
        new_dia_objects = []

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
            dia_obj_idx = dia_collection._id_to_index[
                score_struct.obj_ids[score_idx]]
            if used_dia_object[dia_obj_idx]:
                continue
            used_dia_object[dia_obj_idx] = True
            used_dia_source[score_idx] = True
            updated_obj_id = score_struct.obj_ids[score_idx]
            updated_dia_objects.append(updated_obj_id)

            dia_collection.dia_objects[dia_obj_idx].append_dia_source(
                dia_source_catalog[int(score_idx)])

        # Argwhere returns a array shape (N, 1) so we access the index
        # thusly to retrieve the value rather than the tuple.
        for (src_idx,) in np.argwhere(np.logical_not(used_dia_source)):
            tmp_src_cat = afwTable.SourceCatalog(dia_source_catalog.schema)
            tmp_src_cat.append(dia_source_catalog[int(src_idx)])
            dia_collection.append(DIAObject(tmp_src_cat))
            new_dia_objects.append(
                dia_collection.dia_objects[-1].id)

        # Return the ids of the DIAObjects in this DIAObjectCollection that
        # were updated or newly created.
        n_updated_dia_objects = len(updated_dia_objects)
        n_unassociated_dia_objects = \
            n_previous_dia_objects - n_updated_dia_objects
        updated_dia_objects.extend(new_dia_objects)
        return pipeBase.Struct(
            updated_and_new_dia_object_ids=updated_dia_objects,
            n_updated_dia_objects=n_updated_dia_objects,
            n_new_dia_objects=len(new_dia_objects),
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
