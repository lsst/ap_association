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
import pandas as pd
from scipy.spatial import cKDTree

import lsst.geom as geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.utils.timer import timeMethod

# Enforce an error for unsafe column/array value setting in pandas.
pd.options.mode.chained_assignment = 'raise'


class AssociationConfig(pexConfig.Config):
    """Config class for AssociationTask.
    """
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

    @timeMethod
    def run(self,
            diaSources,
            diaObjects):
        """Associate the new DiaSources with existing DiaObjects.

        Parameters
        ----------
        diaSources : `pandas.DataFrame`
            New DIASources to be associated with existing DIAObjects.
        diaObjects : `pandas.DataFrame`
            Existing diaObjects from the Apdb.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components.

            - ``"matchedDiaSources"`` : DiaSources that were matched. Matched
              Sources have their diaObjectId updated and set to the id of the
              diaObject they were matched to. (`pandas.DataFrame`)
            - ``"unAssocDiaSources"`` : DiaSources that were not matched.
              Unassociated sources have their diaObject set to 0 as they
              were not associated with any existing DiaObjects.
              (`pandas.DataFrame`)
            - ``"nUpdatedDiaObjects"`` : Number of DiaObjects that were
              matched to new DiaSources. (`int`)
            - ``"nUnassociatedDiaObjects"`` : Number of DiaObjects that were
              not matched a new DiaSource. (`int`)
        """
        diaSources = self.check_dia_source_radec(diaSources)
        if len(diaObjects) == 0:
            return pipeBase.Struct(
                matchedDiaSources=pd.DataFrame(columns=diaSources.columns),
                unAssocDiaSources=diaSources,
                nUpdatedDiaObjects=0,
                nUnassociatedDiaObjects=0)

        matchResult = self.associate_sources(diaObjects, diaSources)

        mask = matchResult.diaSources["diaObjectId"] != 0

        return pipeBase.Struct(
            matchedDiaSources=matchResult.diaSources[mask].reset_index(drop=True),
            unAssocDiaSources=matchResult.diaSources[~mask].reset_index(drop=True),
            nUpdatedDiaObjects=matchResult.nUpdatedDiaObjects,
            nUnassociatedDiaObjects=matchResult.nUnassociatedDiaObjects)

    def check_dia_source_radec(self, dia_sources):
        """Check that all DiaSources have non-NaN values for RA/DEC.

        If one or more DiaSources are found to have NaN values, throw a
        warning to the log with the ids of the offending sources. Drop them
        from the table.

        Parameters
        ----------
        dia_sources : `pandas.DataFrame`
            Input DiaSources to check for NaN values.

        Returns
        -------
        trimmed_sources : `pandas.DataFrame`
            DataFrame of DiaSources trimmed of all entries with NaN values for
            RA/DEC.
        """
        nan_mask = (dia_sources.loc[:, "ra"].isnull()
                    | dia_sources.loc[:, "dec"].isnull())
        if np.any(nan_mask):
            nan_idxs = np.argwhere(nan_mask.to_numpy()).flatten()
            for nan_idx in nan_idxs:
                self.log.warning(
                    "DiaSource %i has NaN value for RA/DEC, "
                    "dropping from association." %
                    dia_sources.loc[nan_idx, "diaSourceId"])
            dia_sources = dia_sources[~nan_mask]
        return dia_sources

    @timeMethod
    def associate_sources(self, dia_objects, dia_sources):
        """Associate the input DIASources with the catalog of DIAObjects.

        DiaObject DataFrame must be indexed on ``diaObjectId``.

        Parameters
        ----------
        dia_objects : `pandas.DataFrame`
            Catalog of DIAObjects to attempt to associate the input
            DIASources into.
        dia_sources : `pandas.DataFrame`
            DIASources to associate into the DIAObjectCollection.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components.

            - ``"diaSources"`` : Full set of diaSources both matched and not.
              (`pandas.DataFrame`)
            - ``"nUpdatedDiaObjects"`` : Number of DiaObjects that were
              associated. (`int`)
            - ``"nUnassociatedDiaObjects"`` : Number of DiaObjects that were
              not matched a new DiaSource. (`int`)
        """
        scores = self.score(
            dia_objects, dia_sources,
            self.config.maxDistArcSeconds * geom.arcseconds)
        match_result = self.match(dia_objects, dia_sources, scores)

        return match_result

    @timeMethod
    def score(self, dia_objects, dia_sources, max_dist):
        """Compute a quality score for each dia_source/dia_object pair
        between this catalog of DIAObjects and the input DIASource catalog.

        ``max_dist`` sets maximum separation in arcseconds to consider a
        dia_source a possible match to a dia_object. If the pair is
        beyond this distance no score is computed.

        Parameters
        ----------
        dia_objects : `pandas.DataFrame`
            A contiguous catalog of DIAObjects to score against dia_sources.
        dia_sources : `pandas.DataFrame`
            A contiguous catalog of dia_sources to "score" based on distance
            and (in the future) other metrics.
        max_dist : `lsst.geom.Angle`
            Maximum allowed distance to compute a score for a given DIAObject
            DIASource pair.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components:

            - ``"scores"``: array of floats of match quality updated DIAObjects
                (array-like of `float`).
            - ``"obj_idxs"``: indexes of the matched DIAObjects in the catalog.
                (array-like of `int`)
            - ``"obj_ids"``: array of floats of match quality updated DIAObjects
                (array-like of `int`).

            Default values for these arrays are
            INF, -1, and -1 respectively for unassociated sources.
        """
        scores = np.full(len(dia_sources), np.inf, dtype=np.float64)
        obj_idxs = np.full(len(dia_sources), -1, dtype=np.int64)
        obj_ids = np.full(len(dia_sources), 0, dtype=np.int64)

        if len(dia_objects) == 0:
            return pipeBase.Struct(
                scores=scores,
                obj_idxs=obj_idxs,
                obj_ids=obj_ids)

        spatial_tree = self._make_spatial_tree(dia_objects)

        max_dist_rad = max_dist.asRadians()

        vectors = self._radec_to_xyz(dia_sources)

        scores, obj_idxs = spatial_tree.query(
            vectors,
            distance_upper_bound=max_dist_rad)
        matched_src_idxs = np.argwhere(np.isfinite(scores))
        obj_ids[matched_src_idxs] = dia_objects.index.to_numpy()[
            obj_idxs[matched_src_idxs]]

        return pipeBase.Struct(
            scores=scores,
            obj_idxs=obj_idxs,
            obj_ids=obj_ids)

    def _make_spatial_tree(self, dia_objects):
        """Create a searchable kd-tree the input dia_object positions.

        Parameters
        ----------
        dia_objects : `pandas.DataFrame`
            A catalog of DIAObjects to create the tree from.

        Returns
        -------
        kd_tree : `scipy.spatical.cKDTree`
            Searchable kd-tree created from the positions of the DIAObjects.
        """
        vectors = self._radec_to_xyz(dia_objects)
        return cKDTree(vectors)

    def _radec_to_xyz(self, catalog):
        """Convert input ra/dec coordinates to spherical unit-vectors.

        Parameters
        ----------
        catalog : `pandas.DataFrame`
            Catalog to produce spherical unit-vector from.

        Returns
        -------
        vectors : `numpy.ndarray`, (N, 3)
            Output unit-vectors
        """
        ras = np.radians(catalog["ra"])
        decs = np.radians(catalog["dec"])
        vectors = np.empty((len(ras), 3))

        sin_dec = np.sin(np.pi / 2 - decs)
        vectors[:, 0] = sin_dec * np.cos(ras)
        vectors[:, 1] = sin_dec * np.sin(ras)
        vectors[:, 2] = np.cos(np.pi / 2 - decs)

        return vectors

    @timeMethod
    def match(self, dia_objects, dia_sources, score_struct):
        """Match DIAsources to DiaObjects given a score.

        Parameters
        ----------
        dia_objects : `pandas.DataFrame`
            A SourceCatalog of DIAObjects to associate to DIASources.
        dia_sources : `pandas.DataFrame`
            A contiguous catalog of dia_sources for which the set of scores
            has been computed on with DIAObjectCollection.score.
        score_struct : `lsst.pipe.base.Struct`
            Results struct with components:

            - ``"scores"``: array of floats of match quality
                updated DIAObjects (array-like of `float`).
            - ``"obj_ids"``: array of floats of match quality
                updated DIAObjects (array-like of `int`).
            - ``"obj_idxs"``: indexes of the matched DIAObjects in the catalog.
                (array-like of `int`)

            Default values for these arrays are
            INF, -1 and -1 respectively for unassociated sources.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components.

            - ``"diaSources"`` : Full set of diaSources both matched and not.
              (`pandas.DataFrame`)
            - ``"nUpdatedDiaObjects"`` : Number of DiaObjects that were
              associated. (`int`)
            - ``"nUnassociatedDiaObjects"`` : Number of DiaObjects that were
              not matched a new DiaSource. (`int`)
        """
        n_previous_dia_objects = len(dia_objects)
        used_dia_object = np.zeros(n_previous_dia_objects, dtype=bool)
        used_dia_source = np.zeros(len(dia_sources), dtype=bool)
        associated_dia_object_ids = np.zeros(len(dia_sources),
                                             dtype=np.uint64)
        n_updated_dia_objects = 0

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
            dia_sources.loc[score_idx, "diaObjectId"] = obj_id
            n_updated_dia_objects += 1

        return pipeBase.Struct(
            diaSources=dia_sources,
            nUpdatedDiaObjects=n_updated_dia_objects,
            nUnassociatedDiaObjects=(n_previous_dia_objects
                                     - n_updated_dia_objects))
