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
import unittest

from lsst.ap.association import \
    AssociationDBSqliteTask, \
    AssociationDBSqliteConfig, \
    make_minimal_dia_source_schema
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
import lsst.afw.image as afwImage
import lsst.afw.image.utils as afwImageUtils
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
from lsst.pex.exceptions import InvalidParameterError
import lsst.pipe.base as pipeBase
import lsst.utils.tests


def create_test_points(point_locs_deg,
                       start_id=0,
                       schema=None,
                       scatter_arcsec=1.0,
                       indexer_ids=None,
                       associated_ids=None):
    """Create dummy DIASources or DIAObjects for use in our tests.

    Parameters
    ----------
    point_locs_deg : array-like (N, 2) of `float`s
        Positions of the test points to create in RA, DEC.
    start_id : `int`
        Unique id of the first object to create. The remaining sources are
        incremented by one from the first id.
    schema : `lsst.afw.table.Schema`
        Schema of the objects to create. Defaults to the DIASource schema.
    scatter_arcsec : `float`
        Scatter to add to the position of each DIASource.
    indexer_ids : `list` of `ints`s
        Id numbers of pixelization indexer to store. Must be the same length
        as the first dimension of point_locs_deg.
    associated_ids : `list` of `ints`s
        Id numbers of associated DIAObjects to store. Must be the same length
        as the first dimension of point_locs_deg.

    Returns
    -------
    test_points : `lsst.afw.table.SourceCatalog`
        Catalog of points to test.
    """
    if schema is None:
        schema = make_minimal_dia_source_schema()
    sources = afwTable.SourceCatalog(schema)

    for src_idx, (ra, dec,) in enumerate(point_locs_deg):
        src = sources.addNew()
        src['id'] = src_idx + start_id
        coord = afwGeom.SpherePoint(ra, dec, afwGeom.degrees)
        if scatter_arcsec > 0.0:
            coord = coord.offset(
                np.random.rand() * 360 * afwGeom.degrees,
                np.random.rand() * scatter_arcsec * afwGeom.arcseconds)
        if indexer_ids is not None:
            src['indexer_id'] = indexer_ids[src_idx]
        if associated_ids is not None:
            src['diaObjectId'] = associated_ids[src_idx]
        src.setCoord(coord)

    return sources


class TestAssociationDBSqlite(unittest.TestCase):

    def setUp(self):
        """Initialize an empty database.
        """

        # CFHT Filters from the camera mapper.
        afwImageUtils.resetFilters()
        afwImageUtils.defineFilter('u', lambdaEff=374, alias="u.MP9301")
        afwImageUtils.defineFilter('g', lambdaEff=487, alias="g.MP9401")
        afwImageUtils.defineFilter('r', lambdaEff=628, alias="r.MP9601")
        afwImageUtils.defineFilter('i', lambdaEff=778, alias="i.MP9701")
        afwImageUtils.defineFilter('z', lambdaEff=1170, alias="z.MP9801")

        assoc_db_config = AssociationDBSqliteConfig()
        assoc_db_config.filter_names = ['u', 'g', 'r', 'i', 'z']
        self.assoc_db = AssociationDBSqliteTask(config=assoc_db_config)
        self.assoc_db.create_tables()
        self.assoc_db._commit()

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
            exposureId=4321,
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
        """Close the database connection and delete the object.
        """
        self.assoc_db.close()
        del self.assoc_db

    def _compare_source_records(self, record_a, record_b):
        """Compare the values stored in two source records.

        This comparison assumes that the schema for record_a is a
        subset of or equal to the schema of record_b.

        Parameters
        ----------
        record_a : `lsst.afw.table.SourceRecord`
        record_b : `lsst.afw.table.SourceRecord`
        """
        for sub_schema in record_a.schema:
            value_a = record_a[sub_schema.getKey()]
            value_b = record_a[sub_schema.getKey()]
            if sub_schema.getField().getTypeString() == 'Angle':
                value_a = value_a.asDegrees()
                value_b = value_b.asDegrees()

            if sub_schema.getField().getTypeString()[0] == 'S':
                self.assertEqual(value_a, value_b)
            elif np.isfinite(value_a) and np.isfinite(value_b):
                if sub_schema.getField().getTypeString() == 'L':
                    self.assertEqual(value_a, value_b)
                else:
                    self.assertAlmostEqual(value_a, value_b)
            else:
                self.assertFalse(np.isfinite(value_a))
                self.assertFalse(np.isfinite(value_b))

    def test_load_dia_objects(self):
        """Test the retrieval of DIAObjects from the database.
        """
        # Create DIAObjects with real positions on the sky with the first
        # point out of the CCD bounding box.
        n_objects = 10
        n_missing_objects = 1
        # Loop backward so the missing point is last.
        object_centers = [
            [self.wcs.pixelToSky(idx, idx).getRa().asDegrees(),
             self.wcs.pixelToSky(idx, idx).getDec().asDegrees()]
            for idx in reversed(np.linspace(-10, 1000, n_objects))]
        dia_objects = create_test_points(
            point_locs_deg=object_centers,
            start_id=0,
            schema=self.assoc_db.get_dia_object_schema(),
            scatter_arcsec=-1)
        for src_idx, dia_object in enumerate(dia_objects):
            dia_object['psFluxMean_g'] = 10000. + np.random.randn() * 100.
            dia_object['psFluxMeanErr_g'] = 100. + np.random.randn() * 10.
            dia_object['psFluxSigma_g'] = 100. + np.random.randn() * 10.

        # Store the DIAObjects.
        self.assoc_db.store_dia_objects(dia_objects, True)

        # Load the DIAObjects using the bounding box and WCS associated with
        # them.
        output_dia_objects = self.assoc_db.load_dia_objects(self.exposure)
        # One of the objects should be outside of the bounding box and will
        # therefore not be loaded.
        self.assertEqual(len(output_dia_objects),
                         n_objects - n_missing_objects)

        # Loop over the 9 output_dia_objects
        for dia_object, created_object in zip(output_dia_objects, dia_objects):
            # HTM trixel for this CCD at level 7.
            created_object["indexer_id"] = 225823
            self._compare_source_records(dia_object, created_object)

    def test_store_dia_objects_no_indexer_id_update(self):
        """Test the storage and retrieval of DIAObjects from the database
        without updating their HTM index.
        """
        # Create DIAObjects with real positions on the sky.
        n_objects = 5
        object_centers = [
            [self.wcs.pixelToSky(idx, idx).getRa().asDegrees(),
             self.wcs.pixelToSky(idx, idx).getDec().asDegrees()]
            for idx in np.linspace(1, 1000, 10)[:n_objects]]
        dia_objects = create_test_points(
            point_locs_deg=object_centers,
            start_id=0,
            schema=self.assoc_db.get_dia_object_schema(),
            scatter_arcsec=1.0)
        for src_idx, dia_object in enumerate(dia_objects):
            dia_object['psFluxMean_g'] = 10000. + np.random.randn() * 100.
            dia_object['psFluxMeanErr_g'] = 100. + np.random.randn() * 10.
            dia_object['psFluxSigma_g'] = 100. + np.random.randn() * 10.

        # Store their values and test if they are preserved after round tripping
        # to the DB.
        self.assoc_db.store_dia_objects(dia_objects, False)
        output_dia_objects = self._retrieve_source_catalog(
            self.assoc_db._dia_object_converter)
        self.assertEqual(len(output_dia_objects), len(dia_objects))
        for dia_object, created_object in zip(output_dia_objects, dia_objects):
            self._compare_source_records(dia_object, created_object)

    def test_store_dia_objects_indexer_id_update(self):
        """Test the storage and retrieval of DIAObjects from the database
        while updating their HTM index.
        """

        # Create DIAObjects with real positions on the sky.
        n_objects = 5
        object_centers = [
            [self.wcs.pixelToSky(idx, idx).getRa().asDegrees(),
             self.wcs.pixelToSky(idx, idx).getDec().asDegrees()]
            for idx in np.linspace(1, 1000, 10)[:n_objects]]
        dia_objects = create_test_points(
            point_locs_deg=object_centers,
            start_id=0,
            schema=self.assoc_db.get_dia_object_schema(),
            scatter_arcsec=1.0)
        # Store and overwrite the same sources this time updating their HTM
        # index.
        for src_idx, dia_object in enumerate(dia_objects):
            dia_object['psFluxMean_g'] = 10000. + np.random.randn() * 100.
            dia_object['psFluxMeanErr_g'] = 100. + np.random.randn() * 10.
            dia_object['psFluxSigma_g'] = 100. + np.random.randn() * 10.
        self.assoc_db.store_dia_objects(dia_objects, True)

        # Retrieve the DIAObjects again and test that their HTM index has
        # been updated properly.
        output_dia_objects = self._retrieve_source_catalog(
            self.assoc_db._dia_object_converter)
        self.assertEqual(len(output_dia_objects), len(dia_objects))
        for dia_object, created_object in zip(output_dia_objects, dia_objects):
            # HTM trixel for this CCD at level 7.
            created_object["indexer_id"] = 225823
            self._compare_source_records(dia_object, created_object)

    def test_indexer_ids(self):
        """Test that the returned HTM pixel indices are returned as expected.
        """
        n_objects = 5
        object_centers = [[0.1 * idx, 0.1 * idx] for idx in range(n_objects)]
        dia_objects = create_test_points(
            point_locs_deg=object_centers,
            start_id=0,
            schema=self.assoc_db.get_dia_object_schema(),
            scatter_arcsec=-1)
        expected_ids = [131072, 253952, 253952, 253952, 253955]
        for obj, indexer_id in zip(dia_objects, expected_ids):
            self.assertEqual(self.assoc_db.compute_indexer_id(obj.getCoord()),
                             indexer_id)

    def test_get_db_filter_id_and_name(self):
        """Test that the filter name mapper to id is working.
        """
        for filter_name, filter_id in zip(self.assoc_db.config.filter_names,
                                          range(len(self.assoc_db.config.filter_names))):
            self.assertEqual(
                self.assoc_db.get_db_filter_id_from_name(filter_name),
                filter_id)
            self.assertEqual(
                self.assoc_db.get_db_filter_name_from_id(filter_id),
                filter_name)
        with self.assertRaises(InvalidParameterError):
            self.assoc_db.get_db_filter_id_from_name("J")
            self.assoc_db.get_db_filter_name_from_id(100)

    def test_store_ccd_visit_info(self):
        """Test storing and retrieving CcdVisit info.
        """
        self.assoc_db.store_ccd_visit_info(self.exposure)
        self.assoc_db._db_cursor.execute("SELECT * FROM CcdVisit")
        stored_values = self.assoc_db._get_ccd_visit_info_from_exposure(
            self.exposure)
        rows = self.assoc_db._db_cursor.fetchall()
        for row in rows:
            for db_value, value in zip(row, stored_values.values()):
                self.assertEqual(db_value, value)

    def test_load_dia_sources(self):
        """Test the retrieval of DIASources from the database.
        """
        n_sources = 5
        dia_sources = create_test_points(
            point_locs_deg=[[0.1, 0.1] for idx in range(n_sources)],
            start_id=0,
            schema=self.assoc_db.get_dia_source_schema(),
            scatter_arcsec=1.0,
            associated_ids=range(n_sources))

        for dia_source in dia_sources:
            dia_source['psFlux'] = 10000. + np.random.randn() * 100.
            dia_source['psFluxErr'] = 100. + np.random.randn() * 10.
            dia_source['filterName'] = self.exposure.getFilter().getName()

        # Store the first set of DIASources and retrieve them using their
        # associated DIAObject id.
        self.assoc_db.store_dia_sources(dia_sources,
                                        range(n_sources),
                                        self.exposure)

        for dia_source in dia_sources:
            tmp_flux = dia_source['psFlux']
            tmp_flux_err = dia_source['psFluxErr']
            dia_source['psFlux'] = tmp_flux / self.flux0
            dia_source['psFluxErr'] = np.sqrt(
                (tmp_flux_err / self.flux0) ** 2 +
                (tmp_flux * self.flux0_err / self.flux0 ** 2) ** 2)
            dia_source['ccdVisitId'] = \
                self.exposure.getInfo().getVisitInfo().getExposureId()

        for dia_object_id, dia_source in zip(range(n_sources), dia_sources):
            stored_dia_sources = self.assoc_db.load_dia_sources([dia_object_id])
            # Should load only one object.
            self.assertEqual(len(stored_dia_sources), 1)
            self._compare_source_records(stored_dia_sources[0], dia_source)

        # Load all stored DIASources at once.
        stored_dia_sources = self.assoc_db.load_dia_sources(range(n_sources))
        self.assertEqual(len(stored_dia_sources), n_sources)
        for dia_source, created_source in zip(stored_dia_sources, dia_sources):
            self._compare_source_records(dia_source, created_source)

        # Test that asking for an id that has no associated sources returns
        # and empty catalog.
        empty_dia_sources = self.assoc_db.load_dia_sources([6])
        self.assertEqual(len(empty_dia_sources), 0)

    def test_store_dia_sources_different_schema(self):
        """Test the storage of DIASources in the database.
        """
        # Create test associated DIAObjects and DIASources.
        schema = afwTable.SourceTable.makeMinimalSchema()
        schema.addField('base_PsfFlux_flux', type='D')
        schema.addField('base_PsfFlux_fluxSigma', type='D')
        schema.addField('junk1', type='L')
        schema.addField('junk2', type='D')
        schema.addField('junk3', type='L')
        schema.addField('junk4', type='D')

        n_sources = 5
        source_centers = [[1. * idx, 1. * idx] for idx in range(n_sources)]
        obj_ids = [idx for idx in range(n_sources)]
        dia_sources = create_test_points(
            point_locs_deg=source_centers,
            start_id=0,
            schema=schema,
            scatter_arcsec=-1)
        for src_idx, dia_source in enumerate(dia_sources):
            dia_source['base_PsfFlux_flux'] = 10000.
            dia_source['base_PsfFlux_fluxSigma'] = 100.

        # Check the DIASources round trip properly. We don't need to be
        # as complex here as the call signature has been almost fully tested
        # here by the ``test_store_catalog_dia_sources`` tests.
        self.assoc_db.store_dia_sources(dia_sources,
                                        range(5),
                                        self.exposure)
        round_trip_dia_source_catalog = self._retrieve_source_catalog(
            self.assoc_db._dia_source_converter)

        # Remake the DIASources with the correct values and columns for
        # comparison.
        dia_sources = create_test_points(
            point_locs_deg=source_centers,
            start_id=0,
            schema=make_minimal_dia_source_schema(),
            scatter_arcsec=-1,
            associated_ids=range(5))
        for src_idx, dia_source in enumerate(dia_sources):
            dia_source['filterName'] = self.exposure.getFilter().getName()
            dia_source['ccdVisitId'] = \
                self.exposure.getInfo().getVisitInfo().getExposureId()
            dia_source['psFlux'] = 10000. / self.flux0
            dia_source['psFluxErr'] = np.sqrt(
                (100. / self.flux0) ** 2 +
                (10000. * self.flux0_err / self.flux0 ** 2) ** 2)

        for stored_dia_source, dia_source in zip(round_trip_dia_source_catalog,
                                                 dia_sources,):

            self._compare_source_records(stored_dia_source, dia_source)

    def test_store_dia_sources(self):
        """Test the storage of DIASources in the database.
        """
        # Create test associated DIAObjects and DIASources.
        n_sources = 5
        source_centers = [[1. * idx, 1. * idx] for idx in range(n_sources)]
        obj_ids = [idx for idx in range(n_sources)]
        dia_sources = create_test_points(
            point_locs_deg=source_centers,
            start_id=0,
            schema=self.assoc_db.get_dia_source_schema(),
            scatter_arcsec=1.0,
            associated_ids=range(5))
        for src_idx, dia_source in enumerate(dia_sources):
            dia_source['psFlux'] = 10000. + np.random.randn() * 100.
            dia_source['psFluxErr'] = 100. + np.random.randn() * 10.

        # Check the DIASources round trip properly. We don't need to be
        # as complex here as the call signature has been almost fully tested
        # here by the ``test_store_catalog_dia_sources`` tests.
        self.assoc_db.store_dia_sources(dia_sources,
                                        range(5),
                                        self.exposure)
        round_trip_dia_source_catalog = self._retrieve_source_catalog(
            self.assoc_db._dia_source_converter)

        for stored_dia_source, dia_source, obj_id, filter_name in zip(
                round_trip_dia_source_catalog,
                dia_sources,
                obj_ids,
                self.assoc_db.config.filter_names):
            dia_source['diaObjectId'] = obj_id
            tmp_flux = dia_source['psFlux']
            tmp_flux_err = dia_source['psFluxErr']
            dia_source['psFlux'] = tmp_flux / self.flux0
            dia_source['psFluxErr'] = np.sqrt(
                (tmp_flux_err / self.flux0) ** 2 +
                (tmp_flux * self.flux0_err / self.flux0 ** 2) ** 2)
            dia_source['filterName'] = self.exposure.getFilter().getName()
            dia_source['ccdVisitId'] = \
                self.exposure.getInfo().getVisitInfo().getExposureId()

            self._compare_source_records(stored_dia_source, dia_source)

    def test_store_catalog_dia_objects(self):
        """Test storing a SourceRecord object in either the dia_objects and
        dia_sources table.
        """

        # Create test associated DIAObjects and DIASources.
        n_objects = 5
        object_centers = [[1. * idx, 1. * idx] for idx in range(n_objects)]
        dia_objects = create_test_points(
            point_locs_deg=object_centers,
            start_id=0,
            schema=self.assoc_db.get_dia_object_schema(),
            scatter_arcsec=1.0)

        for src_idx, dia_object in enumerate(dia_objects):
            dia_object['psFluxMean_g'] = 10000. + np.random.randn() * 100.
            dia_object['psFluxMeanErr_g'] = 100. + np.random.randn() * 10.
            dia_object['psFluxSigma_g'] = 100. + np.random.randn() * 10.

        # Check the DIAObjects round trip properly.
        self.assoc_db._store_catalog(
            dia_objects, self.assoc_db._dia_object_converter)
        round_trip_dia_object_catalog = self._retrieve_source_catalog(
            self.assoc_db._dia_object_converter)
        for stored_dia_object, dia_object in zip(round_trip_dia_object_catalog,
                                                 dia_objects):
            self._compare_source_records(stored_dia_object, dia_object)

    def test_store_catalog_dia_sources(self):
        """Test storing a DIASources with the full functionality of store
        catalogs.
        """
        self._store_catalog_dia_sources(True, True)

    def test_store_catalog_dia_sources_no_id(self):
        """Test storing a DIASources with without updating the associated
        DIAObject ids.
        """
        self._store_catalog_dia_sources(False, True)

    def test_store_catalog_dia_sources_no_exposure(self):
        """Test storing a DIASources with without updating the exposure
        properties
        """
        self._store_catalog_dia_sources(True, False)

    def test_store_catalog_dia_sources_no_id_no_exposure(self):
        """Test storing a DIASources with without updating the associated
        DIAObject ids or the exposure properties.
        """
        self._store_catalog_dia_sources(False, False)

    def _store_catalog_dia_sources(self, use_ids=False, use_exposure=False):
        """Test storing a SourceRecord object in either the dia_objects and
        dia_sources table.
        """

        # Create test associated DIAObjects and DIASources.
        n_sources = 5
        source_centers = [[1. * idx, 1. * idx] for idx in range(n_sources)]
        obj_ids = [idx for idx in range(n_sources)]
        dia_sources = create_test_points(
            point_locs_deg=source_centers,
            start_id=0,
            schema=self.assoc_db.get_dia_source_schema(),
            scatter_arcsec=1.0)
        for src_idx, dia_source in enumerate(dia_sources):
            if not use_ids:
                dia_source['diaObjectId'] = src_idx
            dia_source['psFlux'] = 10000. + np.random.randn() * 100.
            dia_source['psFluxErr'] = 100. + np.random.randn() * 10.
            if not use_exposure:
                dia_source['filterName'] = self.assoc_db.config.filter_names[src_idx]
                dia_source['ccdVisitId'] = 1234 + src_idx

        # Check the DIASources round trip properly.
        if use_ids and not use_exposure:
            self.assoc_db._store_catalog(dia_sources,
                                         self.assoc_db._dia_source_converter,
                                         obj_ids=obj_ids)
        elif not use_ids and use_exposure:
            self.assoc_db._store_catalog(dia_sources,
                                         self.assoc_db._dia_source_converter,
                                         exposure=self.exposure)
        elif use_ids and use_exposure:
            self.assoc_db._store_catalog(dia_sources,
                                         self.assoc_db._dia_source_converter,
                                         obj_ids=obj_ids,
                                         exposure=self.exposure)
        else:
            self.assoc_db._store_catalog(dia_sources,
                                         self.assoc_db._dia_source_converter)
        self.assoc_db._commit()
        round_trip_dia_source_catalog = self._retrieve_source_catalog(
            self.assoc_db._dia_source_converter)

        for stored_dia_source, dia_source, obj_id, filter_name in zip(
                round_trip_dia_source_catalog,
                dia_sources,
                obj_ids,
                self.assoc_db.config.filter_names):
            if use_ids:
                dia_source['diaObjectId'] = obj_id
            if use_exposure:
                tmp_flux = dia_source['psFlux']
                tmp_flux_err = dia_source['psFluxErr']
                dia_source['psFlux'] = tmp_flux / self.flux0
                dia_source['psFluxErr'] = np.sqrt(
                    (tmp_flux_err / self.flux0) ** 2 +
                    (tmp_flux * self.flux0_err / self.flux0 ** 2) ** 2)
                dia_source['filterName'] = self.exposure.getFilter().getName()
                dia_source['ccdVisitId'] = \
                    self.exposure.getInfo().getVisitInfo().getExposureId()

            self._compare_source_records(stored_dia_source, dia_source)

    def _retrieve_source_catalog(self, converter):
        """Convenience method for retrieving a source catalog object from the
        DB.

        Parameters
        ----------
        converter : `lsst.ap.association.SqliteDBConverter`
            converter defining the table and schema to store.

        Returns
        -------
        source_catalog : `lsst.afw.table.SourceCatalog`
            SourceCatalog of the requested objects.
        """
        self.assoc_db._db_cursor.execute(
            "SELECT * FROM %s" % converter.table_name)

        rows = self.assoc_db._db_cursor.fetchall()
        output_source_catalog = afwTable.SourceCatalog(converter.schema)
        output_source_catalog.reserve(len(rows))
        for row in rows:
            output_source_catalog.append(converter.source_record_from_db_row(row))

        return output_source_catalog


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
