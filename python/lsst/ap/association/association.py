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

from .diaCalculation import DiaObjectCalculationTask

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
    diaCalculation = pexConfig.ConfigurableField(
        target=DiaObjectCalculationTask,
        doc="Task to compute summary statistics for DiaObjects.",
    )

    def setDefaults(self):
        self.diaCalculation.plugins = ["ap_meanPosition",
                                       "ap_HTMIndex",
                                       "ap_nDiaSources",
                                       "ap_diaObjectFlag",
                                       "ap_meanFlux",
                                       "ap_percentileFlux",
                                       "ap_sigmaFlux",
                                       "ap_chi2Flux",
                                       "ap_madFlux",
                                       "ap_skewFlux",
                                       "ap_minMaxFlux",
                                       "ap_maxSlopeFlux",
                                       "ap_meanErrFlux",
                                       "ap_linearFit",
                                       "ap_stetsonJ",
                                       "ap_meanTotFlux",
                                       "ap_sigmaTotFlux"]

    def validate(self):
        if "ap_HTMIndex" not in self.diaCalculation.plugins:
            raise ValueError("AssociationTask requires the ap_HTMIndex plugin "
                             "be enabled for proper insertion into the Apdb.")


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
        self.makeSubtask("diaCalculation")

    @pipeBase.timeMethod
    def run(self,
            diaSources,
            diaObjects,
            diaSourceHistory):
        """Associate the new DiaSources with existing or new DiaObjects,
        updating the DiaObjects.

        Parameters
        ----------
        diaSources : `pandas.DataFrame`
            New DIASources to be associated with existing DIAObjects.
        diaObjects : `pandas.DataFrame`
            Existing diaObjects from the Apdb.
        diaSourceHistory : `pandas.DataFrame`
            12 month DiaSource history of the loaded ``diaObjects``.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components.

            - ``diaObjects`` : Complete set of dia_objects covering the input
              exposure. Catalog contains newly created, updated, and untouched
              diaObjects. (`pandas.DataFrame`)
            - ``updatedDiaObjects`` : Subset of DiaObjects that were updated
              or created during processing. (`pandas.DataFrame`)
            - ``diaSources`` : DiaSources detected in this ccdVisit with
              associated diaObjectIds. (`pandas.DataFrame`)
        """
        diaSources = self.check_dia_source_radec(diaSources)

        matchResult = self.associate_sources(diaObjects, diaSources)

        diaObjects = diaObjects.append(matchResult.new_dia_objects,
                                       sort=True)

        if diaObjects.index.has_duplicates:
            self.log.warn(
                "Duplicate DiaObjects created after association. This may "
                "cause downstream pipeline issues. Dropping duplicated rows.")
            # Drop duplicates via index and keep the first appearance.
            diaObjects = diaObjects.groupby(diaObjects.index).first()

        # Now that we know the DiaObjects our new DiaSources are associated
        # with, we index the new DiaSources the same way as the full history
        # and merge the tables.
        diaSources.set_index(["diaObjectId", "filterName", "diaSourceId"],
                             drop=False,
                             inplace=True)
        mergedDiaSourceHistory = diaSourceHistory.append(diaSources, sort=True)
        if mergedDiaSourceHistory.index.has_duplicates:
            self.log.warn(
                "Duplicate DiaSources created after association and merging "
                "with history. This may cause downstream pipeline issues. "
                "Dropping duplicated rows.")
            # Drop duplicates via index and keep the first appearance. Reset
            # due to the index shape being slight different thatn expected.
            mergedDiaSourceHistory = mergedDiaSourceHistory.groupby(
                mergedDiaSourceHistory.index).first().reset_index(drop=True)
            mergedDiaSourceHistory.set_index(
                ["diaObjectId", "filterName", "diaSourceId"],
                drop=False,
                inplace=True)

        # Get the current filter being processed.
        filterName = diaSources["filterName"].iat[0]

        # Update previously existing DIAObjects with the information from their
        # newly association DIASources and create new DIAObjects from
        # unassociated sources.
        updatedResults = self.diaCalculation.run(
            diaObjects,
            mergedDiaSourceHistory,
            matchResult.associated_dia_object_ids,
            filterName)

        allDiaObjects = updatedResults.diaObjectCat
        updatedDiaObjects = updatedResults.updatedDiaObjects
        if allDiaObjects.index.has_duplicates:
            self.log.warn(
                "Duplicate DiaObjects (loaded + updated) created after "
                "DiaCalculation. This may cause downstream pipeline issues. "
                "Dropping duplicated rows.")
            # Drop duplicates via index and keep the first appearance.
            allDiaObjects = allDiaObjects.groupby(allDiaObjects.index).first()
        if updatedDiaObjects.index.has_duplicates:
            self.log.warn(
                "Duplicate DiaObjects (updated) created after "
                "DiaCalculation. This may cause downstream pipeline issues. "
                "Dropping duplicated rows.")
            # Drop duplicates via index and keep the first appearance.
            updatedDiaObjects = updatedDiaObjects.groupby(
                updatedDiaObjects.index).first()

        return pipeBase.Struct(
            diaObjects=allDiaObjects,
            updatedDiaObjects=updatedDiaObjects,
            diaSources=diaSources,
        )

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
                    | dia_sources.loc[:, "decl"].isnull())
        if np.any(nan_mask):
            nan_idxs = np.argwhere(nan_mask.to_numpy()).flatten()
            for nan_idx in nan_idxs:
                self.log.warning(
                    "DiaSource %i has NaN value for RA/DEC, "
                    "dropping from association." %
                    dia_sources.loc[nan_idx, "diaSourceId"])
            dia_sources = dia_sources[~nan_mask]
        return dia_sources

    @pipeBase.timeMethod
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
        result : `lsst.pipeBase.Struct`
            Results struct with components:

            - ``updated_and_new_dia_object_ids`` : ids of new and updated
              dia_objects as the result of association. (`list` of `int`).
            - ``new_dia_objects`` : Newly created DiaObjects from
              unassociated diaSources. (`pandas.DataFrame`)
            - ``n_updated_dia_objects`` : Number of previously known
              dia_objects with newly associated DIASources. (`int`).
            - ``n_new_dia_objects`` : Number of newly created DIAObjects from
              unassociated DIASources (`int`).
            - ``n_unupdated_dia_objects`` : Number of previous DIAObjects that
              were not associated to a new DIASource (`int`).
        """

        scores = self.score(
            dia_objects, dia_sources,
            self.config.maxDistArcSeconds * geom.arcseconds)
        match_result = self.match(dia_objects, dia_sources, scores)

        self._add_association_meta_data(match_result)

        return match_result

    @pipeBase.timeMethod
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
        decs = np.radians(catalog["decl"])
        vectors = np.empty((len(ras), 3))

        sin_dec = np.sin(np.pi / 2 - decs)
        vectors[:, 0] = sin_dec * np.cos(ras)
        vectors[:, 1] = sin_dec * np.sin(ras)
        vectors[:, 2] = np.cos(np.pi / 2 - decs)

        return vectors

    @pipeBase.timeMethod
    def match(self, dia_objects, dia_sources, score_struct):
        """Match DIAsources to DIAObjects given a score and create new
        DIAObject Ids for new unassociated DIASources.

        Parameters
        ----------
        dia_objects : `pandas.DataFrame`
            A SourceCatalog of DIAObjects to associate to DIASources.
        dia_sources : `pandas.DataFrame`
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
            - ``new_dia_objects`` : Newly created DiaObjects from unassociated
              diaSources. (`pandas.DataFrame`)
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
        new_dia_objects = []

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
            dia_sources.loc[score_idx, "diaObjectId"] = obj_id

        # Argwhere returns a array shape (N, 1) so we access the index
        # thusly to retrieve the value rather than the tuple
        for (src_idx,) in np.argwhere(np.logical_not(used_dia_source)):
            src_id = dia_sources.loc[src_idx, "diaSourceId"]
            new_dia_objects.append(self._initialize_dia_object(src_id))
            associated_dia_object_ids[src_idx] = src_id
            dia_sources.loc[src_idx, "diaObjectId"] = src_id
            n_new_dia_objects += 1

        if len(new_dia_objects) > 0:
            new_dia_objects = pd.DataFrame(data=new_dia_objects)
        else:
            # Create a junk DiaObject to get the columns.
            tmpObj = self._initialize_dia_object(0)
            new_dia_objects = pd.DataFrame(data=new_dia_objects,
                                           columns=tmpObj.keys())
        new_dia_objects.set_index("diaObjectId", inplace=True, drop=False)

        # Return the ids of the DIAObjects in this DIAObjectCollection that
        # were updated or newly created.
        n_unassociated_dia_objects = \
            n_previous_dia_objects - n_updated_dia_objects
        return pipeBase.Struct(
            associated_dia_object_ids=associated_dia_object_ids,
            new_dia_objects=new_dia_objects,
            n_updated_dia_objects=n_updated_dia_objects,
            n_new_dia_objects=n_new_dia_objects,
            n_unassociated_dia_objects=n_unassociated_dia_objects,)

    def _initialize_dia_object(self, objId):
        """Create a new DiaObject with values required to be initialized by the
        Ppdb.

        Parameters
        ----------
        objid : `int`
            ``diaObjectId`` value for the of the new DiaObject.

        Returns
        -------
        diaObject : `dict`
            Newly created DiaObject with keys:

            ``diaObjectId``
                Unique DiaObjectId (`int`).
            ``pmParallaxNdata``
                Number of data points used for parallax calculation (`int`).
            ``nearbyObj1``
                Id of the a nearbyObject in the Object table (`int`).
            ``nearbyObj2``
                Id of the a nearbyObject in the Object table (`int`).
            ``nearbyObj3``
                Id of the a nearbyObject in the Object table (`int`).
            ``?PSFluxData``
                Number of data points used to calculate point source flux
                summary statistics in each bandpass (`int`).
        """
        new_dia_object = {"diaObjectId": objId,
                          "pmParallaxNdata": 0,
                          "nearbyObj1": 0,
                          "nearbyObj2": 0,
                          "nearbyObj3": 0,
                          "flags": 0}
        for f in ["u", "g", "r", "i", "z", "y"]:
            new_dia_object["%sPSFluxNdata" % f] = 0
        return new_dia_object

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
