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

from __future__ import absolute_import, division, print_function

import numpy as np
from scipy.spatial import cKDTree

import lsst.afw.table as afwTable


class DIAObjectCollection(object):
    """ A collection of DIAObjects with convenience functions for scoring and
    matching DIASources into the collection of DIAObjects.

    Attributes
    ----------
    dia_objects : a list of DIAObjects
        List of DIAObjects representing this collection the current (e.g.)
        visit.
    is_updated : bool
        Bool representing that the internal summary statistics of each
        DIAObject has been computed for the current set of DIASources it
        contains.
    is_valid_tree : bool
        Bool represetnting that the internal spatial tree structure is valid
        for the current set of DIAObjects.
    """

    def __init__(self, dia_objects):
        """ Initialize a collection of dia_objects.

        Store and update a list of dia_objects in the collection.

        Parameters
        ----------
        dia_objects : a list of DIAObjects
            List of DIAObjects that represent the collection of variable
            objects for this visit. Each DIAObject will be individually
            updated if not already currently updated with the latest
            DIASources.

        Returns
        -------
        A DIAObjectCollection instance
        """
        self.dia_objects = dia_objects
        self._id_to_index = {}
        for idx, dia_object in enumerate(self.dia_objects):
            self._id_to_index[dia_object.get('id')] = idx
        self._is_updated = False
        self._is_valid_tree = False
        self.update_dia_objects()
        self.update_spatial_tree()

        # Run internal method to create a spatial index on the dia_objects
        # in this collection for fast pair searching later.

    def get_dia_object(self, id):
        """ Retrive an individual DIAObject from this collection using its
        catalog id.

        Parameters
        ----------
        id : int
            id of the DIAObject to retrive
        
        Return
        ------
        A DIAObject
        """
        return self.dia_objects[self._id_to_index[id]]

    def get_dia_object_ids(self):
        """ Retrive the ids of the DIAObjects stored in this collection.

        Parameters
        ----------
        id : int
            id of the DIAObject to retrive
        
        Return
        ------
        A DIAObject
        """
        return self._id_to_index.keys()

    def update_dia_objects(self, force=False):
        """ Update the summary statistics of all DIAObjects in this
        collection.

        Loop through the DIAObjects that make up this DIAObjectCollection and
        update them as needed. Optional `force` variable forces the DIAObjects
        within the collelection to be updated regardless of their `is_updated
        state.

        Parameters
        ----------
        force : bool (optional)
            Force the DIAObjects to update regardless of their internal
            `is_updated` status.

        Returns
        -------
        bool successfully updated
        """
        self._is_updated = False

        for dia_object in self.dia_objects:
            if not dia_object.is_updated or force:
                dia_object.update()
                self._is_valid_tree = False

        self._is_updated = True

        return self._is_updated

    def update_spatial_tree(self):
        """ Update the internal searchable spatial tree on the DIAObjects.

        Returns
        -------
        bool successfully updated
        """
        self._is_valid_tree = False
        if not self._is_updated:
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
        """ Add a new DIAObject to this collection.

        Parameters
        ----------
        dia_object : A DIAObject class instance
            Input dia_object to append to this collection.

        Returns
        -------
        None
        """

        self._is_updated = False
        self._is_valid_tree = False

        self._id_to_index[dia_object.get('id')] = len(self.dia_objects)
        self.dia_objects.append(dia_object)

        return None

    def score(self, dia_source_catalog, max_dist):
        """ Compute a quality score for each dia_source/dia_object pair
        between this collection and an input diat_source catalog.

        max_dist sets maximum seperation in arcseconds to consider a
        dia_source a possible match to a dia_object. If the pair is
        beyond this distance no score is computed.

        Parameters
        ----------
        dia_source_catalog : an lsst.afw.SourceCatalog
            A contiguous catalog of dia_sources to "score" based on distance
            and (in the future) other metrics.
        max_dist : lsst.afw.geom.Angle
            Maximum allowed distance to compute a score for a given DIAObject
            DIASource pair.

        Returns
        -------
        (N DIAObject by M DIASource) float array of between DIAObjects and
        DIASources.
        """
        if not self._is_valid_tree:
            self.update_spatial_tree()
        
        scores = np.ones((len(self.dia_objects),
                          len(dia_source_catalog))) * np.nan
        for src_idx, dia_source in enumerate(dia_source_catalog):

            src_point = dia_source.getCoord().getVector()
            dists, indices = self._spatial_tree.query(
                src_point, k=3)
            scores[indices, src_idx] = np.where(dists < max_dist.asRadians(), 
                                                dists, np.nan)

        return scores

    def match(self, dia_source_catalog, scores):
        """ Append DIAsources to DIAObjects given a score and create new
        DIAObjects in this collection from DIASources with poor scores.

        Parameters
        ----------
        dia_source_catalog : an lsst.afw.SourceCatalog
            A contiguous catalog of dia_sources for which the set of scores
            has been computed on with DIAObjectCollection.score.
        scores : (N DIAObject by M DIASource) float array
            A ndarray containing the scores between each DIAObject and
            each DIASource.

        Returns
        -------
        Indices of newly updated and created DIAObjects
        """

        used_dia_object = np.zeros(len(self.dia_objects), dtype=np.bool)
        used_dia_source = np.zeros(len(self.dia_objects), dtype=np.bool)

        updated_and_new_dia_objects = []

        score_args = scores.argsort(axis=None)
        score_indices = np.array(
            [score_args % len(dia_source_catalog),
             score_args % len(self.dia_objects)]).transpose()
        for score_idx in score_indices:
            print(score_idx)
            print(scores[score_idx])
            if not np.isfinite(scores[score_idx[0], score_idx[1]]):
                break
            if used_dia_object[score_idx[0]] or \
               used_dia_source[score_idx[1]]:
                continue
            used_dia_object[score_idx[0]] = True
            used_dia_source[score_idx[1]] = True
            updated_and_new_dia_objects.append(score_idx[0])

            self.dia_objects[score_idx[0]].append_dia_source(
                dia_source_catalog[score_idx[1]])

        print(used_dia_source)
        print(used_dia_object)

        n_new_objects = 0
        print(np.argwhere(used_dia_source == False))
        for src_idx in np.argwhere(used_dia_source == False):
            tmp_src_cat = afwTable.SourceCatalog(
                dia_source_catalog.schema)
            print(type(dia_source_catalog[src_idx]))
            tmp_src_cat.append(dia_source_catalog[src_idx])
            self.append(DIAObject(tmp_src_cat))
            used_dia_source[src_idx] = True
            n_new_objects += 1

        return np.concatenate(
            [updated_and_new_dia_objects,
             np.arange(len(self.dia_objects),
                       len(self.dia_objects) + n_new_objects,
                       dtype=np.int)])

    @property
    def is_updated(self):
        """ Return the status of the internal DIAObjects and if their summary
        statistics have been properly updated.

        Return
        ------
        bool
        """

        return self._is_updated

    @property
    def is_valid_tree(self):
        """ Return the status of the internal spatial search tree.

        If the tree has not been updated with the current positions of
        all DIAObjects internal to this collection we return false.

        Return
        ------
        bool
        """

        return self._is_valid_tree
