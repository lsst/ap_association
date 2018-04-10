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

"""Define collections of DIAObjects and how to associate them with
DIASources.
"""

from __future__ import absolute_import, division, print_function

import numpy as np
from scipy.spatial import cKDTree

import lsst.afw.table as afwTable
import lsst.pipe.base as pipeBase

from .dia_object import DIAObject

__all__ = ["DIAObjectCollection"]


class DIAObjectCollection(object):
    """A collection of DIAObjects with convenience functions for scoring and
    matching DIASources into the collection of DIAObjects.

    Parameters
    ----------
    dia_objects : `list` of `lsst.ap.association.DIAObjects`
        List of DIAObjects representing this collection the current (e.g.)
        visit.
    """

    def __init__(self, dia_objects):
        self.dia_objects = dia_objects
        self._id_to_index = {}
        for idx, dia_object in enumerate(self.dia_objects):
            self._id_to_index[dia_object.id] = idx
        self._is_updated = False
        self._is_valid_tree = False
        self.update_dia_objects()
        self.update_spatial_tree()

        # Run internal method to create a spatial index on the dia_objects
        # in this collection for fast pair searching later.

    def get_dia_object(self, id):
        """Retrieve an individual DIAObject from this collection using its
        catalog id.

        Parameters
        ----------
        id : `int`
            id of the DIAObject to retrieve

        Returns
        -------
        dia_object : `lsst.ap.association.DIAObject`
            DIAObject with the ``id`` specified
        """
        return self.dia_objects[self._id_to_index[id]]

    def get_dia_object_ids(self):
        """Retrieve the ids of the DIAObjects stored in this collection.

        Parameters
        ----------
        id : `int`
            id of the DIAObject to retrieve

        Returns
        -------
        dia_object_ids : `list` of `int`s
            List of the ``ids`` of all DIAObjects contained in this collection.
        """
        return list(self._id_to_index.keys())

    def update_dia_objects(self, force=False):
        """Update the summary statistics of all DIAObjects in this
        collection.

        Loop through the DIAObjects that make up this DIAObjectCollection and
        update them as needed. Optional `force` variable forces the DIAObjects
        within the collection to be updated regardless of their `is_updated
        state. We set the variable _is_updated to True after this is run
        to assert that this method has been run and all summary statistics
        in the DIAObejcts are valid for their current associated DIASources.

        Parameters
        ----------
        force : `bool` (optional)
            Force the DIAObjects to update regardless of their internal
            ``is_updated`` status.

        Returns
        -------
        is_updated : `bool`
            Successfully updated
        """
        self._is_updated = False

        for dia_object in self.dia_objects:
            if not dia_object.is_updated or force:
                dia_object.update()
                self._is_valid_tree = False

        self._is_updated = True

        return self._is_updated

    def update_spatial_tree(self):
        """Update the internal search able spatial tree on the DIAObjects.

        Returns
        -------
        is_updated : `bool`
            Successfully updated
        """
        self._is_valid_tree = False
        if not self._is_updated:
            return self._is_valid_tree
        if len(self.dia_objects) == 0:
            self._is_valid_tree = True
            self._spatial_tree = None
            return self._is_valid_tree

        xyzs = np.empty((len(self.dia_objects), 3))
        for obj_idx in range(len(self.dia_objects)):
            tmp_coord = self.dia_objects[obj_idx].dia_object_record.getCoord()
            tmp_vect = tmp_coord.getVector()
            xyzs[obj_idx, 0] = tmp_vect[0]
            xyzs[obj_idx, 1] = tmp_vect[1]
            xyzs[obj_idx, 2] = tmp_vect[2]
        self._spatial_tree = cKDTree(xyzs)

        self._is_valid_tree = True

        return self._is_valid_tree

    def append(self, dia_object):
        """Add a new DIAObject to this collection.

        Parameters
        ----------
        dia_object : `lsst.ap.association.DIAObject`
            Input dia_object to append to this collection.
        """

        self._is_updated = False
        self._is_valid_tree = False

        self._id_to_index[dia_object.id] = len(self.dia_objects)
        self.dia_objects.append(dia_object)

        return None

    def score(self, dia_source_catalog, max_dist):
        """Compute a quality score for each dia_source/dia_object pair
        between this collection and an input diat_source catalog.

        max_dist sets maximum separation in arcseconds to consider a
        dia_source a possible match to a dia_object. If the pair is
        beyond this distance no score is computed.

        Parameters
        ----------
        dia_source_catalog : `lsst.afw.table.SourceCatalog`
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
        if not self._is_valid_tree:
            self.update_spatial_tree()

        scores = np.ones(len(dia_source_catalog)) * np.inf
        obj_ids = -1 * np.ones(len(dia_source_catalog), dtype=np.int)

        if len(self.dia_objects) == 0:
            return pipeBase.Struct(
                scores=scores,
                indices=obj_ids)

        for src_idx, dia_source in enumerate(dia_source_catalog):

            src_point = dia_source.getCoord().getVector()
            dist, obj_idx = self._spatial_tree.query(src_point)
            if dist < max_dist.asRadians():
                scores[src_idx] = dist
                obj_ids[src_idx] = self.dia_objects[obj_idx].id

        return pipeBase.Struct(
            scores=scores,
            obj_ids=obj_ids)

    def match(self, dia_source_catalog, score_struct):
        """Append DIAsources to DIAObjects given a score and create new
        DIAObjects in this collection from DIASources with poor scores.

        Parameters
        ----------
        dia_source_catalog : `lsst.afw.table.SourceCatalog`
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

        n_previous_dia_objects = len(self.dia_objects)
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
            dia_obj_idx = self._id_to_index[score_struct.obj_ids[score_idx]]
            if used_dia_object[dia_obj_idx]:
                continue
            used_dia_object[dia_obj_idx] = True
            used_dia_source[score_idx] = True
            updated_obj_id = score_struct.obj_ids[score_idx]
            updated_dia_objects.append(updated_obj_id)

            self.dia_objects[dia_obj_idx].append_dia_source(
                dia_source_catalog[int(score_idx)])

        # Argwhere returns a array shape (N, 1) so we access the index
        # thusly to retrieve the value rather than the tuple.
        for (src_idx,) in np.argwhere(np.logical_not(used_dia_source)):
            tmp_src_cat = afwTable.SourceCatalog(dia_source_catalog.schema)
            tmp_src_cat.append(dia_source_catalog[int(src_idx)])
            self.append(DIAObject(tmp_src_cat))
            new_dia_objects.append(
                self.dia_objects[-1].id)

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

    @property
    def is_updated(self):
        """Return the status of the internal DIAObjects and if their summary
        statistics have been properly updated.
        """

        return self._is_updated

    @property
    def is_valid_tree(self):
        """Return the status of the internal spatial search tree.

        If the tree has not been updated with the current positions of
        all DIAObjects internal to this collection we return false.
        """

        return self._is_valid_tree
