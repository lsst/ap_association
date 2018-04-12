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
    make_minimal_dia_object_schema, \
    make_minimal_dia_source_schema
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
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
        self.assoc_db = AssociationDBSqliteTask()
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
            if sub_schema.getField().getTypeString() == 'L':
                self.assertEqual(record_a[sub_schema.getKey()],
                                 record_b[sub_schema.getKey()])
            elif sub_schema.getField().getTypeString() == 'Angle':
                self.assertAlmostEqual(
                    record_a[sub_schema.getKey()].asDegrees(),
                    record_b[sub_schema.getKey()].asDegrees())
            else:
                self.assertAlmostEqual(record_a[sub_schema.getKey()],
                                       record_b[sub_schema.getKey()])

    def test_load_dia_objects(self):
        """Test the retrieval of DIAObjects from the database.
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
            schema=make_minimal_dia_object_schema(),
            scatter_arcsec=-1)

        # Store the DIAObjects.
        self.assoc_db.store_dia_objects(dia_objects, True)

        # Load the DIAObjects using the bounding box and WCS associated with
        # them.
        output_dia_objects = self.assoc_db.load_dia_objects(self.expMd)
        self.assertEqual(len(output_dia_objects), len(dia_objects))
        for dia_object, created_object in zip(output_dia_objects, dia_objects):
            # HTM trixel for this CCD at level 7.
            created_object["indexer_id"] = 225823
            self._compare_source_records(dia_object, created_object)

        # Add new DIAObjects outside of the CCD area and test that they are
        # not loaded.
        n_objects = 5
        object_centers = [
            [self.wcs.pixelToSky(idx, idx).getRa().asDegrees(),
             self.wcs.pixelToSky(idx, idx).getDec().asDegrees()]
            for idx in np.linspace(-100, 0, 10)[:n_objects]]
        outside_dia_objects = create_test_points(
            point_locs_deg=object_centers,
            start_id=0,
            schema=make_minimal_dia_object_schema(),
            scatter_arcsec=-1)

        # Store the DIAObjects.
        self.assoc_db.store_dia_objects(outside_dia_objects, True)
        output_dia_objects = self.assoc_db.load_dia_objects(self.expMd)

        # The DIAObjects loaded should be the same as previously stored
        # objects.
        self.assertEqual(len(output_dia_objects), len(dia_objects))
        for dia_object, created_object in zip(output_dia_objects, dia_objects):
            # HTM trixel for this CCD at level 7.
            created_object["indexer_id"] = 225823
            self._compare_source_records(dia_object, created_object)

    def test_store_dia_objects(self):
        """Test the storage and retrieval of DIAObjects from the database.
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
            schema=make_minimal_dia_object_schema(),
            scatter_arcsec=1.0)

        # Store their values and test if they are preserved after round tripping
        # to the DB.
        self.assoc_db.store_dia_objects(dia_objects, False)
        output_dia_objects = self._retrieve_source_catalog(
            self.assoc_db._dia_object_converter)
        self.assertEqual(len(output_dia_objects), len(dia_objects))
        for dia_object, created_object in zip(output_dia_objects, dia_objects):
            self._compare_source_records(dia_object, created_object)

        # Store and overwrite the same sources this time updating their HTM
        # index.
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
            schema=make_minimal_dia_object_schema(),
            scatter_arcsec=-1)
        for obj in dia_objects:
            print(self.assoc_db.compute_indexer_id(obj.getCoord()))

    def test_load_dia_sources(self):
        """Test the retrieval of DIASources from the database.
        """
        n_sources = 5
        dia_sources = create_test_points(
            point_locs_deg=[[0.1, 0.1] for idx in range(n_sources)],
            start_id=0,
            schema=make_minimal_dia_source_schema(),
            scatter_arcsec=1.0)

        # Store the first set of DIASources and retrieve them using their
        # associated DIAObject id.
        self.assoc_db.store_dia_sources(dia_sources, [0] * n_sources)
        stored_dia_sources = self.assoc_db.load_dia_sources([0])
        self.assertEqual(len(stored_dia_sources), n_sources)
        for dia_source, created_source in zip(stored_dia_sources, dia_sources):
            self._compare_source_records(dia_source, created_source)

        # Test that asking for an id that has no associated sources returns
        # and empty catalog.
        empty_dia_sources = self.assoc_db.load_dia_sources([1])
        self.assertEqual(len(empty_dia_sources), 0)

        # Create new DIASoures associated with a different DIAObject id.
        dia_sources = create_test_points(
            point_locs_deg=[[0.1, 0.1] for idx in range(n_sources)],
            start_id=n_sources,
            schema=make_minimal_dia_source_schema(),
            scatter_arcsec=1.0)
        self.assoc_db.store_dia_sources(dia_sources, [1] * n_sources)

        # Load all the associated DIASources and test that the returned catalog
        # is the correct length.
        stored_dia_sources = self.assoc_db.load_dia_sources([0, 1])
        self.assertEqual(len(stored_dia_sources), 2 * n_sources)

        # Load the DIASources associated with the second DIAObject and test
        # their values.
        stored_dia_sources = self.assoc_db.load_dia_sources([1])
        self.assertEqual(len(stored_dia_sources), n_sources)
        for dia_source, created_source in zip(stored_dia_sources, dia_sources):
            self._compare_source_records(dia_source, created_source)

    def test_store_dia_sources(self):
        """Test the storage of DIASources in the database.
        """
        # Create DIASources
        n_sources = 5
        dia_sources = create_test_points(
            point_locs_deg=[[0.1, 0.1] for idx in range(n_sources)],
            start_id=0,
            schema=make_minimal_dia_source_schema(),
            scatter_arcsec=1.0,
            associated_ids=[0 for idx in range(n_sources)])

        # Store the DIASources
        self.assoc_db.store_dia_sources(dia_sources, range(n_sources))

        # Retrieve and test DIASources.
        stored_dia_sources = self._retrieve_source_catalog(
            self.assoc_db._dia_object_converter)
        self.assertEqual(len(stored_dia_sources), n_sources)
        for dia_source, created_source in zip(stored_dia_sources, dia_sources):
            self._compare_source_records(dia_source, created_source)

    def test_store_catalog_objects(self):
        """Test storing a SourceRecord object in either the dia_objects and
        dia_sources table.
        """

        # Create test associated DIAObjects and DIASources.
        dia_objects = create_test_points(
            point_locs_deg=[[0.0, 0.0],
                            [1.0, 1.0]],
            start_id=0,
            schema=make_minimal_dia_object_schema(),
            scatter_arcsec=1.0,
            associated_ids=None)
        dia_sources = create_test_points(
            point_locs_deg=[[0.0, 0.0],
                            [1.0, 1.0]],
            start_id=0,
            schema=make_minimal_dia_source_schema(),
            scatter_arcsec=1.0,
            associated_ids=[0, 1])

        # Check the DIAObjects round trip properly.
        self._store_source_catalog(dia_objects,
                                   self.assoc_db._dia_object_converter)
        round_trip_dia_object_catalog = self._retrieve_source_catalog(
            self.assoc_db._dia_object_converter)
        for stored_dia_object, dia_object in zip(round_trip_dia_object_catalog,
                                                 dia_objects):
            self._compare_source_records(stored_dia_object, dia_object)

        # Check the DIASources round trip properly.
        self._store_source_catalog(dia_sources,
                                   self.assoc_db._dia_source_converter)
        round_trip_dia_source_catalog = self._retrieve_source_catalog(
            self.assoc_db._dia_source_converter)
        for stored_dia_source, dia_source in zip(round_trip_dia_source_catalog,
                                                 dia_sources):
            self._compare_source_records(stored_dia_source, dia_source)

    def _store_source_catalog(self, source_catalog, converter):
        """Convenience method for storing a source catalog object in the DB.

        Parameters
        ----------
        source_catalog : `lsst.afw.table.SourceCatalog`
            SourceCatalog to store.
        converter : `lsst.ap.association.SqliteDBConverter`
            converter defining the table and schema to store.
        """
        self.assoc_db._store_catalog(
            source_catalog, converter)
        self.assoc_db._commit()

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
