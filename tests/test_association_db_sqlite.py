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


def create_test_points(n_points=5,
                       start_id=0,
                       schema=None,
                       point_locs_deg=[[0.0, 0.0]],
                       scatter_arcsec=1.0,
                       indexer_ids=None,
                       associated_ids=None):
    """ Create dummy DIASources or DIAObjects for use in our tests.

    Parameters
    ----------
    n_points : `int`
        Number of fake sources to create for testing.
    start_id : `int`
        Unique id of the first object to create. The remaining sources are
        incremented by one from the first id.
    schema : `lsst.afw.table.Schema`
        Schema of the objects to create. Defaults to the DIASource schema.
    point_locs_deg : `list` of `float`s
        Positions of the test points to create.
    scatter_arcsec : `float`
        Scatter to add to the position of each DIASource.
    indexer_ids : `list` of `ints`s
        Id numbers of pixelization indexer to store.
    associated_ids : `list` of `ints`s
        Id numbers of associated DIAObjects to store..

    Returns
    -------
    `lsst.afw.table.SourceCatalog`
    """
    if schema is None:
        schema = make_minimal_dia_source_schema()
    sources = afwTable.SourceCatalog(schema)

    for src_idx in range(n_points):
        src = sources.addNew()
        src['id'] = src_idx + start_id
        coord = afwGeom.SpherePoint(point_locs_deg[src_idx][0],
                                    point_locs_deg[src_idx][1],
                                    afwGeom.degrees)
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

    def test_store_and_load_dia_objects(self):
        """Test the storage and retrieval of DIAObjects from the database.
        """
        n_objects = 5
        object_centers = [
            [self.wcs.pixelToSky(idx, idx).getRa().asDegrees(),
             self.wcs.pixelToSky(idx, idx).getDec().asDegrees()]
            for idx in np.linspace(1, 1000, 10)[:n_objects]]
        dia_objects = create_test_points(
            n_points=n_objects,
            start_id=0,
            schema=make_minimal_dia_object_schema(),
            point_locs_deg=object_centers,
            scatter_arcsec=0.0)

        self.assoc_db.store_dia_objects(dia_objects, True)

        # Get exposure bounding box and load data.
        bbox = afwGeom.Box2D(self.exposure.getBBox())
        wcs = self.exposure.getWcs()
        expMd = pipeBase.Struct(
            bbox=bbox,
            wcs=wcs,)
        output_dia_objects = self.assoc_db.load_dia_objects(expMd)

        self.assertEqual(len(output_dia_objects), len(dia_objects))
        for dia_object, created_object in zip(output_dia_objects, dia_objects):
            # HTM trixel for this CCD.
            created_object["indexer_id"] = 225823
            self._compare_source_records(dia_object, created_object)

        # Test overwrite and storage with empty indexer_ids
        dia_objects = create_test_points(
            n_points=n_objects,
            start_id=0,
            schema=make_minimal_dia_object_schema(),
            point_locs_deg=object_centers,
            scatter_arcsec=1.0,
            indexer_ids=[0 for idx in range(n_objects)],)
        self.assoc_db.store_dia_objects(dia_objects, False)

        # This should load an empty catalog.
        output_dia_objects = self.assoc_db.load_dia_objects(expMd)
        self.assertEqual(0, len(output_dia_objects))
        # This should load all objects.
        output_dia_objects = self.assoc_db._get_dia_object_catalog([0])
        self.assertEqual(n_objects, len(output_dia_objects))

    def test_store_and_get_dia_object_catalog(self):
        """Test the storage and retrieval of DIAObjects from the database.
        """
        n_objects = 5
        dia_objects = create_test_points(
            n_points=n_objects,
            start_id=0,
            schema=make_minimal_dia_object_schema(),
            point_locs_deg=[[0.0, 0.0] for idx in range(n_objects)],
            scatter_arcsec=1.0,
            indexer_ids=[0 for idx in range(n_objects)],)

        self.assoc_db.store_dia_objects(dia_objects, False)

        output_dia_objects = self.assoc_db._get_dia_object_catalog([0])
        self.assertEqual(len(output_dia_objects), n_objects)
        for dia_object, created_object in zip(output_dia_objects, dia_objects):
            self._compare_source_records(dia_object, created_object)

    def test_store_and_load_dia_sources(self):
        """Test the retrieval of DIASources from the database.
        """
        n_sources = 5
        dia_sources = create_test_points(
            n_points=n_sources,
            start_id=0,
            schema=make_minimal_dia_source_schema(),
            point_locs_deg=[[0.1, 0.1] for idx in range(n_sources)],
            scatter_arcsec=1.0,
            associated_ids=[0 for idx in range(n_sources)])

        self.assoc_db.store_dia_sources(dia_sources, [0 for idx in range(n_sources)])

        src_cat = self.assoc_db.load_dia_sources([0])
        self.assertEqual(len(src_cat), n_sources)
        for dia_source, created_source in zip(src_cat, dia_sources):
            self._compare_source_records(dia_source, created_source)

        dia_sources = create_test_points(
            n_points=n_sources,
            start_id=5,
            schema=make_minimal_dia_source_schema(),
            point_locs_deg=[[0.0, 0.0] for idx in range(n_sources)],
            scatter_arcsec=1.0,
            associated_ids=[1 for idx in range(n_sources)])

        self.assoc_db.store_dia_sources(dia_sources)

        src_cat = self.assoc_db.load_dia_sources([1])
        self.assertEqual(len(src_cat), n_sources)
        for dia_source, created_source in zip(src_cat, dia_sources):
            self._compare_source_records(dia_source, created_source)

    def test_store_record_objects(self):
        """Test storing a SourceRecord object in either the dia_objects and
        dia_sources table.
        """
        dia_objects = create_test_points(
            n_points=1,
            start_id=0,
            schema=make_minimal_dia_object_schema(),
            point_locs_deg=[[0.0, 0.0]],
            scatter_arcsec=1.0,
            associated_ids=None)
        dia_sources = create_test_points(
            n_points=1,
            start_id=0,
            schema=make_minimal_dia_source_schema(),
            point_locs_deg=[[0.0, 0.0]],
            scatter_arcsec=1.0,
            associated_ids=[0])

        dia_object_record = self._store_and_retrieve_source_record(
            dia_objects[0],
            self.assoc_db._dia_object_converter)

        dia_source_record = self._store_and_retrieve_source_record(
            dia_sources[0],
            self.assoc_db._dia_source_converter)

        self._compare_source_records(dia_object_record, dia_objects[0])
        self._compare_source_records(dia_source_record, dia_sources[0])

    def _store_and_retrieve_source_record(self,
                                          source_record,
                                          converter):
        """Convenience method for round tripping a source record object.

        Parameters
        ----------
        source_record : `lsst.afw.table.SourceRecord`
            SourceRecord to store.
        converter : `lsst.ap.association.SqliteDBConverter`
            converter defining the table and schema to store.

        Returns
        -------
        source_record : `lsst.afw.table.SourceRecord`
            SourceRecord of the requested object.
        """
        self.assoc_db._store_record(
            source_record, converter)
        self.assoc_db._commit()

        self.assoc_db._db_cursor.execute(
            "SELECT * FROM %s" % converter.table_name)
        return converter.source_record_from_db_row(
            self.assoc_db._db_cursor.fetchone())


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
