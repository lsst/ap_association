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

import numpy as np
import pandas as pd
import unittest

import lsst.geom as geom
import lsst.utils.tests
from lsst.ap.association import AssociationTask
from lsst.ap.association.association import AssociationConfig


class TestAssociationTask(unittest.TestCase):

    def setUp(self):
        """Create sets of diaSources and diaObjects.
        """
        rng = np.random.default_rng(1234)
        self.nObjects = 5
        scatter = 0.1/3600
        self.diaObjects = pd.DataFrame(data=[
            {"ra": 0.04*(idx + 1), "dec": 0.04*(idx + 1),
             "diaObjectId": idx + 1}
            for idx in range(self.nObjects)])
        self.diaObjects.set_index("diaObjectId", drop=False, inplace=True)
        self.nSources = 5
        self.diaSources = pd.DataFrame(data=[
            {"ra": 0.04*idx + scatter*rng.uniform(-1, 1),
             "dec": 0.04*idx + scatter*rng.uniform(-1, 1),
             "diaSourceId": idx + 1 + self.nObjects, "diaObjectId": 0, "trailLength": 5.5*idx,
             "flags": 0}
            for idx in range(self.nSources)])
        self.diaSourceZeroScatter = pd.DataFrame(data=[
            {"ra": 0.04*idx,
             "dec": 0.04*idx,
             "diaSourceId": idx + 1 + self.nObjects, "diaObjectId": 0, "trailLength": 5.5*idx,
             "flags": 0}
            for idx in range(self.nSources)])

    def test_run(self):
        """Test the full task by associating a set of diaSources to
        existing diaObjects.
        """
        config = AssociationTask.ConfigClass()
        assocTask = AssociationTask(config=config)
        results = assocTask.run(self.diaSources, self.diaObjects)

        self.assertEqual(results.nUpdatedDiaObjects, len(self.diaObjects) - 1)
        self.assertEqual(results.nUnassociatedDiaObjects, 1)
        self.assertEqual(len(results.matchedDiaSources),
                         len(self.diaObjects) - 1)
        self.assertEqual(len(results.unAssocDiaSources), 1)
        np.testing.assert_array_equal(results.matchedDiaSources["diaObjectId"].values, [1, 2, 3, 4])
        np.testing.assert_array_equal(results.unAssocDiaSources["diaObjectId"].values, [0])

    def test_run_no_existing_objects(self):
        """Test the run method with a completely empty database.
        """
        assocTask = AssociationTask()
        results = assocTask.run(
            self.diaSources,
            pd.DataFrame(columns=["ra", "dec", "diaObjectId", "trailLength"]))
        self.assertEqual(results.nUpdatedDiaObjects, 0)
        self.assertEqual(results.nUnassociatedDiaObjects, 0)
        self.assertEqual(len(results.matchedDiaSources), 0)
        self.assertTrue(np.all(results.unAssocDiaSources["diaObjectId"] == 0))

    def test_associate_sources(self):
        """Test the performance of the associate_sources method in
        AssociationTask.
        """
        assoc_task = AssociationTask()
        assoc_result = assoc_task.associate_sources(
            self.diaObjects, self.diaSources)

        for test_obj_id, expected_obj_id in zip(
                assoc_result.diaSources["diaObjectId"].to_numpy(),
                [0, 1, 2, 3, 4]):
            self.assertEqual(test_obj_id, expected_obj_id)
        np.testing.assert_array_equal(assoc_result.diaSources["diaObjectId"].values, [0, 1, 2, 3, 4])

    def test_score_and_match(self):
        """Test association between a set of sources and an existing
        DIAObjectCollection.
        """

        assoc_task = AssociationTask()
        score_struct = assoc_task.score(self.diaObjects,
                                        self.diaSourceZeroScatter,
                                        1.0 * geom.arcseconds)
        # Source 0 has no nearby DIAObject (closest is ~144" away).
        self.assertNotIn(0, score_struct.src_idx)
        # Sources 1..4 each have at least one candidate at ~zero score.
        for src_idx in range(1, len(self.diaSources)):
            mask = score_struct.src_idx == src_idx
            self.assertTrue(np.any(mask))
            self.assertAlmostEqual(score_struct.scores[mask].min(), 0.0,
                                   places=10)

        # Linear-assignment match: 4 sources match 4 distinct objects;
        # the fifth object is left unassociated.
        match_result = assoc_task.match(
            self.diaObjects, self.diaSources, score_struct)
        self.assertEqual(match_result.nUpdatedDiaObjects, 4)
        self.assertEqual(match_result.nUnassociatedDiaObjects, 1)

    def test_remove_nan_dia_sources(self):
        """Test removing DiaSources with NaN locations.
        """
        self.diaSources.loc[2, "ra"] = np.nan
        self.diaSources.loc[3, "dec"] = np.nan
        self.diaSources.loc[4, "ra"] = np.nan
        self.diaSources.loc[4, "dec"] = np.nan
        assoc_task = AssociationTask()
        out_dia_sources = assoc_task.check_dia_source_radec(self.diaSources)
        self.assertEqual(len(out_dia_sources), len(self.diaSources) - 3)

    def test_score_falls_back_without_error_columns(self):
        """Score uses chord distance when raErr/decErr columns are absent.
        """
        task = AssociationTask()
        result = task.score(self.diaObjects, self.diaSources,
                            1.0 * geom.arcseconds)
        # 4 of 5 sources have at least one candidate.
        self.assertEqual(len(np.unique(result.src_idx)), 4)
        # Chord distance on the unit sphere — values are radians.
        self.assertTrue(np.all(result.scores < 1e-3))

    def test_score_falls_back_when_errors_all_nan(self):
        """Empty (all-NaN) raErr/decErr columns trigger the chord fallback.
        """
        objects = self.diaObjects.copy()
        objects["raErr"] = np.nan
        objects["decErr"] = np.nan
        sources = self.diaSources.copy()
        sources["raErr"] = np.nan
        sources["decErr"] = np.nan
        task = AssociationTask()
        result = task.score(objects, sources, 1.0 * geom.arcseconds)
        self.assertEqual(len(np.unique(result.src_idx)), 4)
        self.assertTrue(np.all(result.scores < 1e-3))

    def test_chi2_accepts_within_uncertainty(self):
        """A 0.05" separation with 0.05" per-axis sigma yields chi^2 ~ 0.5.
        """
        sig_deg = 0.05 / 3600.0
        objects = pd.DataFrame([
            {"ra": 1.0, "dec": 1.0, "raErr": sig_deg, "decErr": sig_deg,
             "diaObjectId": 1}
        ]).set_index("diaObjectId", drop=False)
        sources = pd.DataFrame([
            {"ra": 1.0, "dec": 1.0 + 0.05 / 3600.0,
             "raErr": sig_deg, "decErr": sig_deg,
             "diaSourceId": 100, "diaObjectId": 0}
        ])
        task = AssociationTask()
        result = task.score(objects, sources, 1.0 * geom.arcseconds)
        self.assertEqual(len(result.src_idx), 1)
        self.assertEqual(result.src_idx[0], 0)
        self.assertEqual(result.obj_idx[0], 0)
        # var per axis = 2 * sig^2 (above the floor); dra=0, ddec=0.05".
        # chi^2 = (0.05)^2 / (2 * 0.05^2) = 0.5, and
        # NLL = 0.5*chi^2 + 0.5*ln(var_ra * var_dec).
        var = 2 * sig_deg ** 2
        expected_nll = 0.5 * 0.5 + 0.5 * np.log(var * var)
        self.assertAlmostEqual(result.scores[0], expected_nll, places=6)

    def test_within_maxdist_matches_despite_high_chi2(self):
        """A 0.5" separation with 0.05" per-axis sigma (chi^2 ~ 50, i.e.
        ~7 sigma) is still matched: the chi^2 only ranks candidates, and
        ``maxDistArcSeconds`` is the sole association gate.
        """
        sig_deg = 0.05 / 3600.0
        objects = pd.DataFrame([
            {"ra": 1.0, "dec": 1.0, "raErr": sig_deg, "decErr": sig_deg,
             "diaObjectId": 1}
        ]).set_index("diaObjectId", drop=False)
        sources = pd.DataFrame([
            {"ra": 1.0, "dec": 1.0 + sig_deg*10,  # 10 sigma in dec
             "raErr": sig_deg, "decErr": sig_deg,
             "diaSourceId": 100, "diaObjectId": 0}
        ])
        task = AssociationTask()
        result = task.score(objects, sources, 1.0 * geom.arcseconds)
        # The pair is kept as a candidate despite the large offset
        # (chi^2 ~ 50) since there is no significance cut.
        self.assertEqual(len(result.src_idx), 1)
        self.assertTrue(np.isfinite(result.scores[0]))
        # End-to-end the source is associated to the object.
        run_result = task.run(sources, objects)
        self.assertEqual(run_result.nUpdatedDiaObjects, 1)
        matched = run_result.matchedDiaSources.set_index("diaSourceId")
        self.assertEqual(int(matched.loc[100, "diaObjectId"]), 1)

    def test_nll_prefers_better_localized_object(self):
        """With two candidate objects, the NLL cost associates to the
        better-localized (tighter) object even though a bare chi^2 would
        pick the looser, farther one, whose larger variance deflates its
        chi^2.
        """
        src_sig = 0.05 / 3600.0
        objects = pd.DataFrame([
            {"ra": 1.0, "dec": 1.0 + 0.15 / 3600.0,  # tight and close
             "raErr": 0.05 / 3600.0, "decErr": 0.05 / 3600.0,
             "diaObjectId": 1},
            {"ra": 1.0, "dec": 1.0 + 0.30 / 3600.0,  # loose and far
             "raErr": 0.5 / 3600.0, "decErr": 0.5 / 3600.0,
             "diaObjectId": 2},
        ]).set_index("diaObjectId", drop=False)
        sources = pd.DataFrame([
            {"ra": 1.0, "dec": 1.0,
             "raErr": src_sig, "decErr": src_sig,
             "diaSourceId": 100, "diaObjectId": 0}
        ])
        # A bare chi^2 would prefer the loose object 2 (smaller chi^2)...
        chi2_tight = (0.15 / 3600.0) ** 2 / (2 * (0.05 / 3600.0) ** 2)
        chi2_loose = (0.30 / 3600.0) ** 2 / ((0.05 / 3600.0) ** 2
                                             + (0.5 / 3600.0) ** 2)
        self.assertLess(chi2_loose, chi2_tight)

        task = AssociationTask()
        score_struct = task.score(objects, sources, 1.0 * geom.arcseconds)
        # ...but the NLL cost prefers the tight object 1. obj_idx are
        # positional: 0 -> diaObjectId 1 (tight), 1 -> diaObjectId 2.
        nll = {int(obj): float(score) for obj, score
               in zip(score_struct.obj_idx, score_struct.scores)}
        self.assertLess(nll[0], nll[1])
        # End-to-end the source associates to the tighter object 1.
        result = task.run(sources, objects)
        matched = result.matchedDiaSources.set_index("diaSourceId")
        self.assertEqual(int(matched.loc[100, "diaObjectId"]), 1)

    def test_chi2_ra_wraparound(self):
        """Pairs straddling RA=0/360 are scored correctly.
        """
        sig_deg = 0.05 / 3600.0
        objects = pd.DataFrame([
            {"ra": 359.99999, "dec": 0.0, "raErr": sig_deg, "decErr": sig_deg,
             "diaObjectId": 1}
        ]).set_index("diaObjectId", drop=False)
        sources = pd.DataFrame([
            {"ra": 0.00001, "dec": 0.0,
             "raErr": sig_deg, "decErr": sig_deg,
             "diaSourceId": 100, "diaObjectId": 0}
        ])
        task = AssociationTask()
        result = task.score(objects, sources, 1.0 * geom.arcseconds)
        # Separation is 0.02" of RA across the wrap — within both cuts.
        self.assertEqual(len(result.src_idx), 1)
        self.assertEqual(result.src_idx[0], 0)
        self.assertEqual(result.obj_idx[0], 0)

    def test_lap_recovers_match_lost_to_greedy(self):
        """LAP keeps a source matched when its nearest object is taken
        by a better-fitting competitor — the case that single-nearest
        greedy used to drop into a new DIAObject.
        """
        sig_deg = 0.5 / 3600.0  # 0.5" per axis
        objects = pd.DataFrame([
            {"ra": 1.0, "dec": 1.0,
             "raErr": sig_deg, "decErr": sig_deg, "diaObjectId": 1},
            {"ra": 1.0, "dec": 1.0 + 0.99 / 3600.0,
             "raErr": sig_deg, "decErr": sig_deg, "diaObjectId": 2},
        ]).set_index("diaObjectId", drop=False)
        # Source 100: closer to obj 1, but obj 2 is also a viable candidate.
        # Source 101: sits exactly on obj 1.
        sources = pd.DataFrame([
            {"ra": 1.0, "dec": 1.0 + 0.4 / 3600.0,
             "raErr": sig_deg, "decErr": sig_deg, "diaSourceId": 100,
             "diaObjectId": 0},
            {"ra": 1.0, "dec": 1.0,
             "raErr": sig_deg, "decErr": sig_deg, "diaSourceId": 101,
             "diaObjectId": 0},
        ])
        task = AssociationTask()
        result = task.run(sources, objects)
        # Both sources matched; LAP picks src 100 → obj 2, src 101 → obj 1.
        self.assertEqual(result.nUpdatedDiaObjects, 2)
        self.assertEqual(result.nUnassociatedDiaObjects, 0)
        matched = result.matchedDiaSources.set_index("diaSourceId")
        self.assertEqual(int(matched.loc[100, "diaObjectId"]), 2)
        self.assertEqual(int(matched.loc[101, "diaObjectId"]), 1)

    def test_match_preserves_diaObjectId_dtype(self):
        """match() preserves the incoming diaObjectId dtype rather than
        overwriting it with uint64. This keeps pd.concat dtype-stable
        downstream — mixing uint64 and int64 silently promotes to
        float64.
        """
        self.assertEqual(self.diaSources["diaObjectId"].dtype, np.int64)
        task = AssociationTask()
        result = task.run(self.diaSources, self.diaObjects)
        self.assertEqual(
            result.matchedDiaSources["diaObjectId"].dtype, np.int64)
        self.assertEqual(
            result.unAssocDiaSources["diaObjectId"].dtype, np.int64)

        sources_u64 = self.diaSources.copy()
        sources_u64["diaObjectId"] = sources_u64["diaObjectId"].astype(
            np.uint64)
        result = task.run(sources_u64, self.diaObjects)
        self.assertEqual(
            result.matchedDiaSources["diaObjectId"].dtype, np.uint64)
        self.assertEqual(
            result.unAssocDiaSources["diaObjectId"].dtype, np.uint64)
        with self.assertRaises(ValueError):
            sources_float = self.diaSources.copy()
            sources_float["diaObjectId"] = sources_float["diaObjectId"].astype(float)
            task.run(sources_float, self.diaObjects)

    def test_fallback_sigma_scores_missing_error_pair(self):
        """A pair whose errors are missing is scored using
        ``fallbackSigmaArcSeconds``: the substituted per-axis uncertainty
        sets the pair's score. The fallback affects only the score
        (ranking), not whether the pair is a candidate.
        """
        sig_deg = 0.05 / 3600.0
        # One row in each catalog carries a real error so chi^2 mode
        # triggers; the rows under test (id 100 / obj 1) have NaN errors.
        objects = pd.DataFrame([
            {"ra": 1.0, "dec": 1.0,
             "raErr": np.nan, "decErr": np.nan, "diaObjectId": 1},
            {"ra": 2.0, "dec": 2.0,
             "raErr": sig_deg, "decErr": sig_deg, "diaObjectId": 2},
        ]).set_index("diaObjectId", drop=False)
        sources = pd.DataFrame([
            {"ra": 1.0, "dec": 1.0 + 0.5 / 3600.0,
             "raErr": np.nan, "decErr": np.nan,
             "diaSourceId": 100, "diaObjectId": 0},
            {"ra": 2.0, "dec": 2.0,
             "raErr": sig_deg, "decErr": sig_deg,
             "diaSourceId": 200, "diaObjectId": 0},
        ])

        def missing_pair_score(fallback):
            cfg = AssociationConfig()
            cfg.fallbackSigmaArcSeconds = fallback
            result = AssociationTask(config=cfg).score(
                objects, sources, 1.0 * geom.arcseconds)
            # The missing-error pair is source row 0 -> object row 0.
            mask = (result.src_idx == 0) & (result.obj_idx == 0)
            self.assertEqual(int(mask.sum()), 1)
            return float(result.scores[mask][0])

        def expected_nll(fallback, chi2):
            # Both rows have missing errors, so the combined per-axis
            # variance is 2 * fallback^2 (deg^2), above the 0.05" floor.
            var = 2 * (fallback / 3600.0) ** 2
            return 0.5 * chi2 + 0.5 * np.log(var * var)

        # ddec = 0.5"; fallback 0.05" -> chi^2 = 50, fallback 0.2" -> 3.125.
        self.assertAlmostEqual(
            missing_pair_score(0.05), expected_nll(0.05, 50.0), places=6)
        self.assertAlmostEqual(
            missing_pair_score(0.2), expected_nll(0.2, 3.125), places=6)

        # Even with the small fallback (large chi^2) the pair is within
        # maxDist, so it is still matched to its object.
        cfg = AssociationConfig()
        cfg.fallbackSigmaArcSeconds = 0.05
        result = AssociationTask(config=cfg).run(sources, objects)
        matched = result.matchedDiaSources.set_index("diaSourceId")
        self.assertEqual(int(matched.loc[100, "diaObjectId"]), 1)
        self.assertEqual(int(matched.loc[200, "diaObjectId"]), 2)

    def test_contention_leaves_extra_source_unassociated(self):
        """When two sources fall within ``maxDistArcSeconds`` of a single
        object, the better-fitting one is matched and the other is left
        unassociated rather than forced onto the same object.
        """
        sig_deg = 0.1 / 3600.0
        objects = pd.DataFrame([
            {"ra": 1.0, "dec": 1.0,
             "raErr": sig_deg, "decErr": sig_deg, "diaObjectId": 1},
        ]).set_index("diaObjectId", drop=False)
        sources = pd.DataFrame([
            {"ra": 1.0, "dec": 1.0 + 0.2 / 3600.0,  # closer to the object
             "raErr": sig_deg, "decErr": sig_deg, "diaSourceId": 100,
             "diaObjectId": 0},
            {"ra": 1.0, "dec": 1.0 + 0.5 / 3600.0,  # farther
             "raErr": sig_deg, "decErr": sig_deg, "diaSourceId": 101,
             "diaObjectId": 0},
        ])
        task = AssociationTask()
        result = task.run(sources, objects)
        # One object, so at most one source can match it.
        self.assertEqual(result.nUpdatedDiaObjects, 1)
        self.assertEqual(result.nUnassociatedDiaObjects, 0)
        self.assertEqual(len(result.matchedDiaSources), 1)
        self.assertEqual(len(result.unAssocDiaSources), 1)
        # The closer source wins the object; the other is unassociated.
        matched = result.matchedDiaSources.set_index("diaSourceId")
        self.assertEqual(int(matched.loc[100, "diaObjectId"]), 1)
        self.assertEqual(
            int(result.unAssocDiaSources["diaSourceId"].iloc[0]), 101)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
