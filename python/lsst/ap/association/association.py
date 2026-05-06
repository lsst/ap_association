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
        """Compute a quality score for each dia_source/dia_object pair.

        The single nearest DIAObject within ``max_dist`` is selected per
        DIASource via a 3D kd-tree on unit vectors. When both inputs
        carry usable ``raErr``/``decErr`` columns, the score is the 2-DOF
        position chi^2 (with non-finite or non-positive uncertainties
        replaced by ``sigmaFloorArcSeconds`` and the combined per-axis
        variance also floored at the same value), which is used only to
        rank the match. ``max_dist`` is the sole association gate. When
        the uncertainty columns are absent or empty, the kd-tree chord
        distance (in radians) is used as the score instead — preserving
        the pre-existing behaviour.

        ``raErr`` and ``decErr`` are taken to follow the LSST DPDD
        convention: each is the marginal uncertainty of the catalog
        coordinate itself in degrees (no cos(dec) factor folded into
        ``raErr``). Under that convention the cos(dec) factor cancels
        between residual and uncertainty, and chi^2 reduces to
        ``dRA^2 / sum(raErr^2) + dDec^2 / sum(decErr^2)``.

        Parameters
        ----------
        dia_objects : `pandas.DataFrame`
            DIAObjects to score against ``dia_sources``. Must contain
            ``ra`` and ``dec``; ``raErr`` and ``decErr`` are used when
            present.
        dia_sources : `pandas.DataFrame`
            DIASources to score. Must contain ``ra`` and ``dec``;
            ``raErr`` and ``decErr`` are used when present.
        max_dist : `lsst.geom.Angle`
            Maximum allowed angular separation; pairs farther than this
            are not considered regardless of their uncertainties.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components:

            - ``scores`` : `numpy.ndarray` of `float`
                Match quality (chi^2 if uncertainty-based, chord
                distance otherwise). INF for unmatched sources.
            - ``obj_idxs`` : `numpy.ndarray` of `int`
                Positional indices of matched DIAObjects; -1 if
                unmatched.
            - ``obj_ids`` : `numpy.ndarray` of `int`
                ``diaObjectId`` of matched DIAObjects; 0 if unmatched.
        """
        n_src = len(dia_sources)
        scores = np.full(n_src, np.inf, dtype=np.float64)
        obj_idxs = np.full(n_src, -1, dtype=np.int64)
        obj_ids = np.full(n_src, 0, dtype=np.int64)

        if len(dia_objects) == 0:
            return pipeBase.Struct(
                scores=scores,
                obj_idxs=obj_idxs,
                obj_ids=obj_ids)

        spatial_tree = self._make_spatial_tree(dia_objects)
        max_dist_rad = max_dist.asRadians()
        vectors = self._radec_to_xyz(dia_sources)

        chord_dists, candidate_obj_idxs = spatial_tree.query(
            vectors,
            distance_upper_bound=max_dist_rad)
        matched = np.isfinite(chord_dists)
        if not np.any(matched):
            return pipeBase.Struct(
                scores=scores,
                obj_idxs=obj_idxs,
                obj_ids=obj_ids)

        use_chi2 = (self._has_position_errors(dia_sources)
                    and self._has_position_errors(dia_objects))
        if use_chi2:
            src_idx = np.flatnonzero(matched)
            obj_idx = candidate_obj_idxs[src_idx]
            chi2 = self._chi2_position(
                dia_sources, dia_objects, src_idx, obj_idx)
            scores[src_idx] = chi2
            obj_idxs[src_idx] = obj_idx
            obj_ids[src_idx] = dia_objects.index.to_numpy()[obj_idx]
        else:
            scores[matched] = chord_dists[matched]
            obj_idxs[matched] = candidate_obj_idxs[matched]
            obj_ids[matched] = dia_objects.index.to_numpy()[
                candidate_obj_idxs[matched]]

        return pipeBase.Struct(
            scores=scores,
            obj_idxs=obj_idxs,
            obj_ids=obj_ids)

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
        # their the best match. Sources with non-finite scores have no
        # match and are skipped — they are left for the new-DIAObject loop.
        finite_idx = np.flatnonzero(np.isfinite(score_struct.scores))
        score_args = finite_idx[np.argsort(score_struct.scores[finite_idx])]
        for score_idx in score_args:
            dia_obj_idx = score_struct.obj_idxs[score_idx]
            if used_dia_object[dia_obj_idx]:
                continue
            used_dia_object[dia_obj_idx] = True
            used_dia_source[score_idx] = True
            obj_id = score_struct.obj_ids[score_idx]
            associated_dia_object_ids[score_idx] = obj_id
            n_updated_dia_objects += 1

        # Assign positionally rather than by label — score_idx values from
        # argsort are 0..N-1, but dia_sources may have a non-contiguous
        # label index after upstream NaN filtering.
        dia_sources = dia_sources.copy()
        dia_sources["diaObjectId"] = associated_dia_object_ids

        return pipeBase.Struct(
            diaSources=dia_sources,
            nUpdatedDiaObjects=n_updated_dia_objects,
            nUnassociatedDiaObjects=(n_previous_dia_objects
                                     - n_updated_dia_objects))
