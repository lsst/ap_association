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

import itertools

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
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
        doc="Maximum distance in arcseconds for a DIASource to be matched "
        "to a DIAObject. This is the sole association radius: every "
        "DIAObject within this distance is a match candidate, ranked by "
        "position chi^2 when uncertainties are available and by angular "
        "distance otherwise.",
        default=1.0,
    )
    sigmaFloorArcSeconds = pexConfig.Field(
        dtype=float,
        doc="Floor on the per-axis position uncertainty (arcsec) used when "
        "computing chi^2.",
        default=0.05,
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

            - ``matchedDiaSources`` : DiaSources that were matched. Matched
              Sources have their diaObjectId updated and set to the id of the
              diaObject they were matched to. (`pandas.DataFrame`)
            - ``unAssocDiaSources`` : DiaSources that were not matched.
              Unassociated sources have their diaObject set to 0 as they
              were not associated with any existing DiaObjects.
              (`pandas.DataFrame`)
            - ``nUpdatedDiaObjects`` : Number of DiaObjects that were
              matched to new DiaSources. (`int`)
            - ``nUnassociatedDiaObjects`` : Number of DiaObjects that were
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
        nan_mask = dia_sources["ra"].isnull() | dia_sources["dec"].isnull()
        if nan_mask.any():
            nan_ids = dia_sources.loc[nan_mask, "diaSourceId"]
            for nan_id in nan_ids:
                self.log.warning(
                    "DiaSource %i has NaN value for RA/DEC, "
                    "dropping from association.", nan_id)
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

            - ``diaSources`` : Full set of diaSources both matched and not.
              (`pandas.DataFrame`)
            - ``nUpdatedDiaObjects`` : Number of DiaObjects that were
              associated. (`int`)
            - ``nUnassociatedDiaObjects`` : Number of DiaObjects that were
              not matched a new DiaSource. (`int`)
        """
        scores = self.score(
            dia_objects, dia_sources,
            self.config.maxDistArcSeconds * geom.arcseconds)
        match_result = self.match(dia_objects, dia_sources, scores)

        return match_result

    @timeMethod
    def score(self, dia_objects, dia_sources, max_dist):
        """Build the candidate (DIASource, DIAObject) match table and
        score every pair.

        For each DIASource, all DIAObjects within ``max_dist`` are retrieved
        from a kd-tree on unit vectors. Each candidate pair is then scored:

        - If both inputs carry usable ``raErr``/``decErr`` columns, the
          score is the 2-DOF position chi^2, so the match prefers the
          best-fitting object rather than the merely-nearest one.
        - Otherwise, the distance (in radians) is used as the score.

        No candidates are dropped by the score itself: every pair within
        ``max_dist`` is retained, and the score is used only to rank them.

        ``raErr`` and ``decErr`` are taken to follow the LSST DPDD
        convention: each is the marginal uncertainty of the catalog
        coordinate itself in degrees (no cos(dec) factor folded into
        ``raErr``). Under that convention the cos(dec) factor cancels
        between residual and uncertainty, and chi^2 reduces to
        ``dRA^2 / sum(raErr^2) + dDec^2 / sum(decErr^2)``.

        ``max_dist`` is both the candidate pre-filter and the
        association radius: every pair within it is retained, and the
        score is used only to rank candidates in the downstream match.

        Parameters
        ----------
        dia_objects, dia_sources : `pandas.DataFrame`
            Must contain ``ra`` and ``dec``; ``raErr`` and ``decErr`` are
            used when present.
        max_dist : `lsst.geom.Angle`
            Hard angular upper bound on candidate pairs.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Flat candidate-pair table:

            - ``src_idx`` : `numpy.ndarray` of `int`
                Positional source index for each surviving pair.
            - ``obj_idx`` : `numpy.ndarray` of `int`
                Positional object index for each surviving pair.
            - ``scores`` : `numpy.ndarray` of `float`
                Cost of each pair (chi^2 if uncertainty-based, chord
                distance in radians otherwise). Lower is better.
            - ``unmatched_cost`` : `float`
                Cost to assign to the synthetic 'no-match' alternative
                in the linear-assignment match — set so that any
                surviving real candidate is preferred.
        """
        n_src = len(dia_sources)
        n_obj = len(dia_objects)
        empty_int = np.empty(0, dtype=np.int64)
        empty_float = np.empty(0, dtype=np.float64)
        max_dist_rad = max_dist.asRadians()
        # Used as the no-match cost in distance mode; always strictly
        # above any real candidate (which the kd-tree caps at max_dist_rad).
        chord_unmatched_cost = max_dist_rad * 1.01 + 1e-300

        if n_obj == 0 or n_src == 0:
            return pipeBase.Struct(
                src_idx=empty_int,
                obj_idx=empty_int,
                scores=empty_float,
                unmatched_cost=chord_unmatched_cost)

        spatial_tree = self._make_spatial_tree(dia_objects)
        src_vectors = self._radec_to_xyz(dia_sources)

        candidate_lists = spatial_tree.query_ball_point(
            src_vectors, r=max_dist_rad)
        counts = np.fromiter(
            (len(c) for c in candidate_lists), dtype=np.int64, count=n_src)
        n_pairs = int(counts.sum())
        if n_pairs == 0:
            return pipeBase.Struct(
                src_idx=empty_int,
                obj_idx=empty_int,
                scores=empty_float,
                unmatched_cost=chord_unmatched_cost)

        src_idx = np.repeat(np.arange(n_src, dtype=np.int64), counts)
        obj_idx = np.fromiter(
            itertools.chain.from_iterable(candidate_lists),
            dtype=np.int64, count=n_pairs)

        if (self._has_position_errors(dia_sources)
                and self._has_position_errors(dia_objects)):
            scores = self._chi2_position(
                dia_sources, dia_objects, src_idx, obj_idx)
            # ``max_dist`` is the sole association gate: every candidate
            # pair is kept and the chi^2 only ranks them. Price the
            # 'no-match' alternative just above the worst candidate so a
            # real match is always preferred. A diaSource will only not be
            # associated if there are more diaSources than diaObjects.
            unmatchedCostDelta = 1.0
            unmatched_cost = (float(scores.max()) + unmatchedCostDelta)
        else:
            obj_vectors = self._radec_to_xyz(dia_objects)
            diffs = src_vectors[src_idx] - obj_vectors[obj_idx]
            scores = np.linalg.norm(diffs, axis=1)
            unmatched_cost = chord_unmatched_cost

        return pipeBase.Struct(
            src_idx=src_idx,
            obj_idx=obj_idx,
            scores=scores,
            unmatched_cost=unmatched_cost)

    @staticmethod
    def _has_position_errors(catalog):
        """Return True iff ``catalog`` carries ``raErr`` and ``decErr``
        columns with at least one finite, positive value in each.
        """
        if "raErr" not in catalog.columns or "decErr" not in catalog.columns:
            return False
        raErr = catalog["raErr"].to_numpy()
        decErr = catalog["decErr"].to_numpy()
        return (bool(np.any(np.isfinite(raErr) & (raErr > 0.0)))
                and bool(np.any(np.isfinite(decErr) & (decErr > 0.0))))

    def _chi2_position(self, dia_sources, dia_objects, src_idx, obj_idx):
        """Return the 2-DOF position chi^2 for paired DIASources/DIAObjects.

        Non-finite or non-positive per-row uncertainties are replaced
        with ``self.config.sigmaFloorArcSeconds``; the combined per-axis
        variance is also floored at that same value to guard against
        pathologically small reported errors.

        Parameters
        ----------
        dia_sources, dia_objects : `pandas.DataFrame`
            Catalogs containing ``ra``, ``dec``, ``raErr``, ``decErr``
            (all in degrees).
        src_idx, obj_idx : `numpy.ndarray` of `int`
            Paired positional indices; ``src_idx[k]`` is matched against
            ``obj_idx[k]``.

        Returns
        -------
        chi2 : `numpy.ndarray` of `float`
             2 degrees of freedom position chi^2, one value per pair.
        """
        sigma_floor_sq_deg = (self.config.sigmaFloorArcSeconds / 3600.0) ** 2

        def err_sq(catalog, col, idx):
            arr = catalog[col].to_numpy()[idx]
            sq = np.square(arr)
            return np.where(np.isfinite(sq) & (arr > 0.0),
                            sq, sigma_floor_sq_deg)

        src_ra = dia_sources["ra"].to_numpy()[src_idx]
        src_dec = dia_sources["dec"].to_numpy()[src_idx]
        obj_ra = dia_objects["ra"].to_numpy()[obj_idx]
        obj_dec = dia_objects["dec"].to_numpy()[obj_idx]

        # Make sure that the RA difference is never greater than +/-180 in
        # either direction.
        dra = ((src_ra - obj_ra) + 180.0) % 360.0 - 180.0
        ddec = src_dec - obj_dec

        var_ra = (err_sq(dia_sources, "raErr", src_idx)
                  + err_sq(dia_objects, "raErr", obj_idx))
        var_dec = (err_sq(dia_sources, "decErr", src_idx)
                   + err_sq(dia_objects, "decErr", obj_idx))
        var_ra = np.maximum(var_ra, sigma_floor_sq_deg)
        var_dec = np.maximum(var_dec, sigma_floor_sq_deg)

        return dra * dra / var_ra + ddec * ddec / var_dec

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
        """Solve a min-cost bipartite matching between sources and objects.

        Each DIASource is given a synthetic 'no-match' alternative (a
        per-source 'ghost' column) carrying ``score_struct.unmatched_cost``.
        A min-weight full bipartite matching is then solved on the sparse
        cost matrix made from real candidate pairs and ghost edges.
        Sources matched to a ghost are reported as unassociated; sources
        matched to a real object inherit that object's ``diaObjectId``.

        When two sources compete for the same object, the source with a
        strictly worse alternative gets that object, and the other source falls
        back to its second-best candidate rather than creating a new DIAObject.

        Parameters
        ----------
        dia_objects, dia_sources : `pandas.DataFrame`
        score_struct : `lsst.pipe.base.Struct`
            Output of `score`: ``src_idx``, ``obj_idx``, ``scores``,
            ``unmatched_cost``.

        Returns
        -------
        result : `lsst.pipe.base.Struct`

            - ``diaSources`` : input source table with ``diaObjectId``
              populated (0 for unmatched). (`pandas.DataFrame`)
            - ``nUpdatedDiaObjects`` : number of DIAObjects matched to a
              new DIASource. (`int`)
            - ``nUnassociatedDiaObjects`` : number of preloaded DIAObjects
              with no matching DIASource. (`int`)
        """
        n_src = len(dia_sources)
        n_obj = len(dia_objects)
        if not pd.api.types.is_integer_dtype(dia_sources["diaObjectId"]):
            raise ValueError(f"diaSource column diaObjectId must be an integer, "
                             f"got {dia_sources['diaObjectId'].dtype} instead")
        # Preserve the dtype of the incoming diaObjectId column so the
        # output is concat-stable downstream (mixing uint64 and int64
        # forces a silent promotion to float64 in pd.concat).
        associated_dia_object_ids = np.zeros(n_src, dtype=dia_sources["diaObjectId"].dtype)
        n_matched = 0

        if n_src > 0:
            # Per-source ghost edges to columns [n_obj, n_obj + n_src).
            ghost_src = np.arange(n_src, dtype=np.int64)
            ghost_obj = n_obj + ghost_src
            ghost_cost = np.full(
                n_src, score_struct.unmatched_cost, dtype=np.float64)

            # scipy's min_weight_full_bipartite_matching removes
            # explicit-zero edges before solving. Nudge them off zero
            # so that perfect-match candidates (chi^2 = 0) are still
            # seen by the LAP.
            real_weights = np.maximum(
                score_struct.scores.astype(np.float64, copy=False), 1e-300)
            rows = np.concatenate(
                [score_struct.src_idx.astype(np.int64, copy=False),
                 ghost_src])
            cols = np.concatenate(
                [score_struct.obj_idx.astype(np.int64, copy=False),
                 ghost_obj])
            weights = np.concatenate([real_weights, ghost_cost])

            biadj = csr_matrix(
                (weights, (rows, cols)),
                shape=(n_src, n_obj + n_src))
            match_src, match_dst = min_weight_full_bipartite_matching(biadj)

            is_real = match_dst < n_obj
            matched_src = match_src[is_real]
            matched_obj = match_dst[is_real]
            n_matched = int(matched_src.size)
            if n_matched > 0:
                associated_dia_object_ids[matched_src] = (
                    dia_objects.index.to_numpy()[matched_obj])

        dia_sources = dia_sources.copy()
        dia_sources["diaObjectId"] = associated_dia_object_ids

        return pipeBase.Struct(
            diaSources=dia_sources,
            nUpdatedDiaObjects=n_matched,
            nUnassociatedDiaObjects=int(n_obj - n_matched))
