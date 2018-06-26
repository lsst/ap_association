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
import os
import tempfile
import unittest

from lsst.afw.cameraGeom.testUtils import DetectorWrapper
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.image.utils as afwImageUtils
import lsst.daf.base as dafBase
import lsst.pipe.base as pipeBase
from lsst.ap.association import \
    AssociationDBSqliteTask, \
    AssociationDBSqliteConfig, \
    AssociationTask, \
    AssociationConfig, \
    make_minimal_dia_object_schema
import lsst.utils.tests
from test_association_db_sqlite import \
    create_test_points


class TestAssociationTask(unittest.TestCase):

    def setUp(self):
        """Create a sqlite3 database with default tables and schemas.
        """
        # CFHT Filters from the camera mapper.
        afwImageUtils.resetFilters()
        afwImageUtils.defineFilter('u', lambdaEff=374, alias="u.MP9301")
        afwImageUtils.defineFilter('g', lambdaEff=487, alias="g.MP9401")
        afwImageUtils.defineFilter('r', lambdaEff=628, alias="r.MP9601")
        afwImageUtils.defineFilter('i', lambdaEff=778, alias="i.MP9701")
        afwImageUtils.defineFilter('z', lambdaEff=1170, alias="z.MP9801")
        self.tmp_file, self.db_file = tempfile.mkstemp(dir=os.path.dirname(__file__))
        self.filter_names = ['g', 'r']
        self.dia_object_schema = make_minimal_dia_object_schema(
            self.filter_names)
        assoc_db_config = AssociationDBSqliteConfig()
        assoc_db_config.db_name = self.db_file
        assoc_db_config.filter_names = self.filter_names

        assoc_db = AssociationDBSqliteTask(config=assoc_db_config)
        assoc_db.create_tables()
        assoc_db.close()

        # metadata taken from CFHT data
        # v695856-e0/v695856-e0-c000-a00.sci_img.fits

        self.metadata = dafBase.PropertySet()

        self.metadata.set("SIMPLE", "T")
        self.metadata.set("BITPIX", -32)
        self.metadata.set("NAXIS", 2)
        self.metadata.set("NAXIS1", 1024)
        self.metadata.set("NAXIS2", 1153)
        self.metadata.set("RADECSYS", 'FK5')
        self.metadata.set("EQUINOX", 2000.)

        self.metadata.setDouble("CRVAL1", 215.604025685476)
        self.metadata.setDouble("CRVAL2", 53.1595451514076)
        self.metadata.setDouble("CRPIX1", 1109.99981456774)
        self.metadata.setDouble("CRPIX2", 560.018167811613)
        self.metadata.set("CTYPE1", 'RA---SIN')
        self.metadata.set("CTYPE2", 'DEC--SIN')

        self.metadata.setDouble("CD1_1", 5.10808596133527E-05)
        self.metadata.setDouble("CD1_2", 1.85579539217196E-07)
        self.metadata.setDouble("CD2_2", -5.10281493481982E-05)
        self.metadata.setDouble("CD2_1", -8.27440751733828E-07)

        self.wcs = afwGeom.makeSkyWcs(self.metadata)
        self.exposure = afwImage.makeExposure(
            afwImage.makeMaskedImageFromArrays(np.ones((1024, 1153))),
            self.wcs)
        detector = DetectorWrapper(id=23, bbox=self.exposure.getBBox()).detector
        visit = afwImage.VisitInfo(
            exposureId=1234,
            exposureTime=200.,
            date=dafBase.DateTime(nsecs=1400000000 * 10**9))
        self.exposure.setDetector(detector)
        self.exposure.getInfo().setVisitInfo(visit)
        self.exposure.setFilter(afwImage.Filter('g'))
        self.flux0 = 10000
        self.flux0_err = 100
        self.exposure.getCalib().setFluxMag0((self.flux0, self.flux0_err))

        bbox = afwGeom.Box2D(self.exposure.getBBox())
        wcs = self.exposure.getWcs()
        self.expMd = pipeBase.Struct(
            bbox=bbox,
            wcs=wcs,)

    def tearDown(self):
        """Delete the database after we are done with it.
        """
        del self.tmp_file
        os.remove(self.db_file)
        del self.db_file
        del self.metadata
        del self.wcs
        del self.exposure

    def test_run(self):
        """Test the run method with a database that already exists and
        contains DIAObjects and Sources.
        """
        dia_objects = self._run_association_and_retrieve_objects(True)
        not_updated_idx = 0
        updated_idx_start = 1
        new_idx_start = 5
        total_expected_dia_objects = 10
        self.assertEqual(len(dia_objects), total_expected_dia_objects)

        # Test to make sure the number of DIAObjects have been properly
        # associated within the db.
        for obj_idx, dia_object in enumerate(dia_objects):
            if obj_idx == not_updated_idx:
                # Test the DIAObject we expect to not be associated with any
                # new DIASources.
                self.assertEqual(dia_object['nDiaSources'], 2)
                self.assertEqual(dia_object.getId(), obj_idx)
            elif updated_idx_start <= obj_idx < new_idx_start:
                # Test that associating to the existing DIAObjects went
                # as planned and test that the IDs of the newly associated
                # DIASources is correct.
                self.assertEqual(dia_object['nDiaSources'], 3)
                self.assertEqual(dia_object.getId(), obj_idx)
            else:
                self.assertEqual(dia_object['nDiaSources'], 1)
                self.assertEqual(dia_object.getId(), obj_idx + 4 + 5)

    def test_run_no_existing_objects(self):
        """Test the run method with a completely empty database.
        """
        dia_objects = self._run_association_and_retrieve_objects(False)
        total_expected_dia_objects = 9
        self.assertEqual(len(dia_objects),
                         total_expected_dia_objects)
        for obj_idx, output_dia_object in enumerate(dia_objects):
            self.assertEqual(output_dia_object['nDiaSources'], 1)
            self.assertEqual(output_dia_object.getId(),
                             obj_idx + 10)
            self.assertEqual(output_dia_object.getId(), obj_idx + 10)

    def _run_association_and_retrieve_objects(self, create_objects=False):
        """Convenience method for testing the Association run method.

        Parameters
        ----------
        create_objects : `bool`
            Boolean specifying if seed DIAObjects and DIASources should be
            inserted into the database before association.

        Return
        ------
        dia_objects : `lsst.afw.table.SourceCatalog`
            Final set of DIAObjects to be tested.
        """
        if create_objects:
            self._store_dia_objects_and_sources()

        source_centers = [
            [self.wcs.pixelToSky(idx, idx).getRa().asDegrees(),
             self.wcs.pixelToSky(idx, idx).getDec().asDegrees()]
            for idx in np.linspace(1, 1000, 10)[1:]]
        dia_sources = create_test_points(
            point_locs_deg=source_centers,
            start_id=10,
            scatter_arcsec=-1)
        for dia_source in dia_sources:
            dia_source["psFlux"] = 10000.
            dia_source["psFluxErr"] = 100.

        assoc_config = AssociationConfig()
        assoc_config.level1_db.value.db_name = self.db_file
        assoc_config.level1_db.value.filter_names = self.filter_names
        assoc_task = AssociationTask(config=assoc_config)

        assoc_task.run(dia_sources, self.exposure)
        assoc_task.level1_db.close()

        assoc_db_config = AssociationDBSqliteConfig()
        assoc_db_config.db_name = self.db_file
        assoc_db_config.filter_names = self.filter_names
        assoc_db = AssociationDBSqliteTask(config=assoc_db_config)

        dia_objects = assoc_db.load_dia_objects(self.exposure)
        assoc_db.close()
        return dia_objects

    def _store_dia_objects_and_sources(self):
        """Method for storing a set of test DIAObjects and sources into
        the L1 database.
        """

        # This should create a DB of 5 DIAObjects with 2 DIASources associated
        # to them. The DIASources are "observed" in g and r.

        # Create an empty database
        assoc_db_config = AssociationDBSqliteConfig()
        assoc_db_config.db_name = self.db_file
        assoc_db_config.filter_names = self.filter_names
        assoc_db = AssociationDBSqliteTask(config=assoc_db_config)
        assoc_db.create_tables()

        # Create DIObjects, give them fluxes, and store them
        n_objects = 5
        object_centers = np.array([
            [self.wcs.pixelToSky(idx, idx).getRa().asDegrees(),
             self.wcs.pixelToSky(idx, idx).getDec().asDegrees()]
            for idx in np.linspace(1, 1000, 10)])
        dia_objects = create_test_points(
            point_locs_deg=object_centers[:n_objects],
            start_id=0,
            schema=self.dia_object_schema,
            scatter_arcsec=-1,)
        # Set the DIAObject fluxes and number of associated sources.
        for dia_object in dia_objects:
            dia_object['nDiaSources'] = 2
            for filter_name in self.filter_names:
                dia_object['psFluxMean_%s' % filter_name] = 1
                dia_object['psFluxMeanErr_%s' % filter_name] = 1
                dia_object['psFluxSigma_%s' % filter_name] = 1
        assoc_db.store_dia_objects(dia_objects, True)

        # Create DIASources, update their ccdVisitId and fluxes, and store
        # them.
        dia_sources = create_test_points(
            point_locs_deg=np.concatenate(
                [object_centers[:n_objects], object_centers[:n_objects]]),
            start_id=0,
            scatter_arcsec=-1,
            associated_ids=[0, 1, 2, 3, 4,
                            0, 1, 2, 3, 4])
        for src_idx, dia_source in enumerate(dia_sources):
            dia_source['ccdVisitId'] = 1232
            if src_idx >= n_objects:
                dia_source['ccdVisitId'] += 1
            dia_source["psFlux"] = 10000 / self.flux0
            dia_source["psFluxErr"] = np.sqrt(
                (100 / self.flux0) ** 2 +
                (10000 * self.flux0_err / self.flux0 ** 2) ** 2)
            if src_idx < n_objects:
                dia_source["filterName"] = 'g'
                dia_source["filterId"] = 1
            else:
                dia_source["filterName"] = 'r'
                dia_source["filterId"] = 2
        assoc_db.store_dia_sources(dia_sources)
        assoc_db.close()

    def test_update_dia_objects(self):
        """Test the update_dia_objects method.
        """
        self._store_dia_objects_and_sources()
        # Create new DIAObjects
        n_sources = 5
        object_centers = np.array([
            [self.wcs.pixelToSky(idx, idx).getRa().asDegrees(),
             self.wcs.pixelToSky(idx, idx).getDec().asDegrees()]
            for idx in np.linspace(1, 1000, 10)])
        dia_sources = create_test_points(
            point_locs_deg=object_centers[:n_sources],
            start_id=10,
            scatter_arcsec=-1)
        # Store raw uncalibrated fluxes for these sources.
        for dia_source in dia_sources:
            dia_source["psFlux"] = 20000
            dia_source["psFluxErr"] = 100
            dia_source["filterName"] = "g"
            dia_source["filterId"] = 1

        # Store them in the DB containing pre-existing objects.
        assoc_db_config = AssociationDBSqliteConfig()
        assoc_db_config.db_name = self.db_file
        assoc_db_config.filter_names = self.filter_names
        assoc_db = AssociationDBSqliteTask(config=assoc_db_config)
        assoc_db.create_tables()
        # Load the Existing DIAObjects for use in the update method.
        loaded_dia_objects = assoc_db.load_dia_objects(self.exposure)

        # Store the new DIASources with associations and an exposure.
        assoc_db.store_dia_sources(dia_sources,
                                   [1, 2, 3, 4, 14],
                                   self.exposure)
        assoc_db.close()
        del assoc_db

        # Create our task and update the stored DIAObjects.
        assoc_config = AssociationConfig()
        assoc_config.level1_db.value.db_name = self.db_file
        assoc_config.level1_db.filter_names = self.filter_names
        assoc_task = AssociationTask(config=assoc_config)
        assoc_task.update_dia_objects(loaded_dia_objects,
                                      [1, 2, 3, 4, 14],
                                      self.exposure)

        # Retrieve the DIAObjects from the DB.
        output_dia_objects = assoc_task.level1_db.load_dia_objects(self.exposure)

        # Data and column names to test.
        test_column_names = [
            'id', 'nDiaSources',
            'psFluxMean_g', 'psFluxMeanErr_g', 'psFluxSigma_g',
            'psFluxMean_r', 'psFluxMeanErr_r', 'psFluxSigma_r']
        test_dia_object_values = [
            {'id': 0, 'nDiaSources': 2,
             'psFluxMean_g': 1., 'psFluxMeanErr_g': 1., 'psFluxSigma_g': 1.,
             'psFluxMean_r': 1., 'psFluxMeanErr_r': 1., 'psFluxSigma_r': 1.},
            {'id': 1, 'nDiaSources': 3,
             'psFluxMean_g': 1.5, 'psFluxMeanErr_g': 0.23570226, 'psFluxSigma_g': 0.70710678,
             'psFluxMean_r': 1., 'psFluxMeanErr_r': 1., 'psFluxSigma_r': 1.},
            {'id': 2, 'nDiaSources': 3,
             'psFluxMean_g': 1.5, 'psFluxMeanErr_g': 0.23570226, 'psFluxSigma_g': 0.70710678,
             'psFluxMean_r': 1., 'psFluxMeanErr_r': 1., 'psFluxSigma_r': 1.},
            {'id': 3, 'nDiaSources': 3,
             'psFluxMean_g': 1.5, 'psFluxMeanErr_g': 0.23570226, 'psFluxSigma_g': 0.70710678,
             'psFluxMean_r': 1., 'psFluxMeanErr_r': 1., 'psFluxSigma_r': 1.},
            {'id': 4, 'nDiaSources': 3,
             'psFluxMean_g': 1.5, 'psFluxMeanErr_g': 0.23570226, 'psFluxSigma_g': 0.70710678,
             'psFluxMean_r': 1., 'psFluxMeanErr_r': 1., 'psFluxSigma_r': 1.},
            {'id': 14, 'nDiaSources': 1,
             'psFluxMean_g': 2., 'psFluxMeanErr_g': np.nan, 'psFluxSigma_g': np.nan,
             'psFluxMean_r': np.nan, 'psFluxMeanErr_r': np.nan, 'psFluxSigma_r': np.nan}
        ]

        # Test that the stored values are as expected.
        self.assertEqual(len(output_dia_objects), 6)
        for dia_object, values in zip(output_dia_objects,
                                      test_dia_object_values):
            for test_name in test_column_names:
                if np.isnan(values[test_name]):
                    self.assertTrue(np.isnan(dia_object[test_name]))
                elif test_name == 'id' or test_name == 'nDiaSources':
                    self.assertEqual(dia_object[test_name],
                                     values[test_name])
                else:
                    self.assertAlmostEqual(dia_object[test_name],
                                           values[test_name])

    def test_associate_sources(self):
        """Test the performance of the associate_sources method in
        AssociationTask.
        """
        n_objects = 5
        dia_objects = create_test_points(
            point_locs_deg=[[0.04 * obj_idx, 0.04 * obj_idx]
                            for obj_idx in range(n_objects)],
            start_id=0,
            schema=self.dia_object_schema,
            scatter_arcsec=-1,)

        n_sources = 5
        dia_sources = create_test_points(
            point_locs_deg=[
                [0.04 * (src_idx + 1),
                 0.04 * (src_idx + 1)]
                for src_idx in range(n_sources)],
            start_id=n_objects,
            scatter_arcsec=0.1)

        assoc_task = AssociationTask()
        assoc_result = assoc_task.associate_sources(
            dia_objects, dia_sources)

        for test_obj_id, expected_obj_id in zip(assoc_result, [1, 2, 3, 4, 9]):
            self.assertEqual(test_obj_id, expected_obj_id)

    def test_score_and_match(self):
        """Test association between a set of sources and an existing
        DIAObjectCollection.

        This also tests that a DIASource that can't be associated within
        tolerance is appended to the DIAObjectCollection as a new
        DIAObject.
        """

        assoc_task = AssociationTask()
        # Create a set of DIAObjects that contain only one DIASource
        n_objects = 5
        dia_objects = create_test_points(
            point_locs_deg=[[0.04 * obj_idx, 0.04 * obj_idx]
                            for obj_idx in range(n_objects)],
            start_id=0,
            schema=self.dia_object_schema,
            scatter_arcsec=-1,)

        n_sources = 5
        dia_sources = create_test_points(
            point_locs_deg=[
                [0.04 * (src_idx + 1),
                 0.04 * (src_idx + 1)]
                for src_idx in range(n_sources)],
            start_id=n_objects,
            scatter_arcsec=-1)

        score_struct = assoc_task.score(dia_objects,
                                        dia_sources,
                                        1.0 * afwGeom.arcseconds)
        self.assertFalse(np.isfinite(score_struct.scores[-1]))
        for src_idx in range(4):
            # Our scores should be extremely close to 0 but not exactly so due
            # to machine noise.
            self.assertAlmostEqual(score_struct.scores[src_idx], 0.0,
                                   places=16)

        # After matching each DIAObject should now contain 2 DIASources
        # except the last DIAObject in this collection which should be
        # newly created during the matching step and contain only one
        # DIASource.
        match_result = assoc_task.match(dia_objects, dia_sources, score_struct)
        updated_ids = match_result.associated_dia_object_ids
        self.assertEqual(len(updated_ids), 5)
        self.assertEqual(match_result.n_updated_dia_objects, 4)
        self.assertEqual(match_result.n_new_dia_objects, 1)
        self.assertEqual(match_result.n_unassociated_dia_objects, 1)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
