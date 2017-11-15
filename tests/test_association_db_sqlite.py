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
    DIAObjectCollection
from lsst.afw.coord import Coord
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.daf.base as dafBase
import lsst.pipe.base as pipeBase
import lsst.utils.tests
from test_dia_collection import create_test_dia_objects


class TestAssociationDBSqlite(unittest.TestCase):

    def setUp(self):
        """ Initialize an empty database.
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

        self.wcs = afwImage.makeWcs(self.metadata)
        self.exposure = afwImage.makeExposure(
            afwImage.makeMaskedImageFromArrays(np.ones((1024, 1153))),
            self.wcs)

    def tearDown(self):
        """ Close the database connection and delete the object.
        """
        self.assoc_db.close()
        del self.assoc_db

    def test_load(self):
        """ Test loading of DIAObjects.
        """

        n_objects = 5
        n_sources_per_object = 2
        object_centers = [
            [self.wcs.pixelToSky(idx, idx).getRa().asDegrees(),
             self.wcs.pixelToSky(idx, idx).getDec().asDegrees()]
            for idx in np.linspace(1, 1000, 10)[:n_objects]]
        dia_objects = create_test_dia_objects(
            n_objects=n_objects,
            n_sources=n_sources_per_object,
            start_id=0,
            object_centers_degrees=object_centers,
            scatter_arcsec=-1.)
        dia_collection = DIAObjectCollection(dia_objects)

        self.assoc_db.store(dia_collection, True)

        bbox = afwGeom.Box2D(self.exposure.getBBox())
        wcs = self.exposure.getWcs()
        expMd = pipeBase.Struct(
            bbox=bbox,
            wcs=wcs,)
        output_dia_collection = self.assoc_db.load(expMd)

        for obj_idx in range(5):
            self.assertEqual(
                output_dia_collection.dia_objects[obj_idx].n_dia_sources,
                n_sources_per_object)
            input_dia_object = output_dia_collection.get_dia_object(
                output_dia_collection.dia_objects[obj_idx].id)
            self._compare_source_records(
                output_dia_collection.dia_objects[obj_idx].dia_object_record,
                input_dia_object.dia_object_record)

            output_src_cat = output_dia_collection.dia_objects[
                obj_idx].dia_source_catalog
            output_src_cat.sort(
                output_src_cat.getSchema().find('id').key)

            for record_a, record_b in zip(
                    output_src_cat,
                    input_dia_object.dia_source_catalog):
                self._compare_source_records(
                    record_a, record_b)

    def _compare_source_records(self, record_a, record_b):
        """ Compare the values stored in two source records.

        This comparison assumes that the schema for record_a is a
        subset of or equal to the schema of record_b.

        Parameters
        ----------
        record_a : lsst.afw.table.SourceRecord
        record_b : lsst.afw.table.SourceRecord
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

    def test_store(self):
        """ Test storing of a DIACollection.
        """
        self._test_store_index_option(True)

    def test_store_no_index_update(self):
        """ Test storing of a DIACollection without updating the spatial index
        of the stored DIAObjects.
        """
        self._test_store_index_option(False)

    def _test_store_index_option(self, update_spatial_index):
        """ Convenience function for testing the store method.

        Parameters
        ---------
        update_spatial_index : bool
            Specify whether to update the spatial index of the DIAObject
            before storage.
        """
        dia_objects = create_test_dia_objects(
            n_objects=1,
            n_sources=1,
            start_angle_degrees=0.1,
            scatter_arcsec=0.0)
        if not update_spatial_index:
            # Set the spatial index to some arbitrary value
            dia_objects[0].dia_object_record.set('indexer_id', 10)
        dia_collection = DIAObjectCollection(dia_objects)

        self.assoc_db.store(dia_collection, update_spatial_index)

        self.assoc_db._db_cursor.execute(
            "SELECT * FROM dia_objects")
        for row in self.assoc_db._db_cursor.fetchall():
            if update_spatial_index:
                # Value is HTM cell number at level=7, RA,DEC=0.1
                dia_objects[0].dia_object_record.set('indexer_id', 253952)
            round_trip_object = \
                self.assoc_db._dia_object_converter.source_record_from_db_row(
                    row)
            self._compare_source_records(
                round_trip_object,
                dia_objects[0].dia_object_record)

    def test_store_updated(self):
        """ Test the storage of newly associated DIAObjects and DIASources.
        """
        dia_objects = create_test_dia_objects(
            n_objects=1,
            n_sources=1,
            start_id=0,
            start_angle_degrees=0.0)
        dia_collection = DIAObjectCollection(dia_objects)

        self.assoc_db.store(dia_collection, True)
        self.assoc_db._commit()

        new_dia_object = create_test_dia_objects(
            n_objects=1,
            n_sources=1,
            start_id=1,
            start_angle_degrees=0.1)
        dia_collection.append(new_dia_object[0])

        # We grab a new source to append to our first source.
        tmp_dia_objects = create_test_dia_objects(
            n_objects=1,
            n_sources=1,
            start_id=2,
            start_angle_degrees=0.0)
        new_src = tmp_dia_objects[0].dia_source_catalog[0]
        dia_collection.dia_objects[0].append_dia_source(new_src)

        dia_collection.update_dia_objects()
        dia_collection.update_spatial_tree()

        self.assoc_db.store_updated(dia_collection, [0, 1])

        self.assoc_db._db_cursor.execute(
            "SELECT indexer_id FROM dia_objects")
        indexer_ids = np.array(
            self.assoc_db._db_cursor.fetchall(), np.int).flatten()

        output_dia_objects = self.assoc_db._get_dia_objects(indexer_ids)

        for obj_idx in range(2):
            if obj_idx == 0:
                self.assertEqual(
                    output_dia_objects[obj_idx].n_dia_sources, 2)
            else:
                self.assertEqual(
                    output_dia_objects[obj_idx].n_dia_sources, 1)
            input_dia_object = dia_collection.get_dia_object(
                dia_collection.dia_objects[obj_idx].id)
            self._compare_source_records(
                output_dia_objects[obj_idx].dia_object_record,
                input_dia_object.dia_object_record)

            output_src_cat = output_dia_objects[
                obj_idx].dia_source_catalog
            output_src_cat.sort(
                output_src_cat.getSchema().find('id').key)

            for record_a, record_b in zip(
                    output_src_cat,
                    input_dia_object.dia_source_catalog):
                self._compare_source_records(
                    record_a, record_b)

    def test_get_dia_objects(self):
        """ Test the retrieval of DIAObjects from the database.
        """
        dia_objects = create_test_dia_objects(
            n_objects=2, n_sources=2, scatter_arcsec=0.0)
        dia_collection = DIAObjectCollection(dia_objects)
        self.assoc_db.store(dia_collection, True)

        self.assoc_db._commit()

        self.assoc_db._db_cursor.execute(
            "SELECT indexer_id FROM dia_objects")
        indexer_ids = np.array(
            self.assoc_db._db_cursor.fetchall(), np.int).flatten()

        output_dia_objects = self.assoc_db._get_dia_objects(indexer_ids)

        for obj_idx in range(2):
            self.assertEqual(
                output_dia_objects[obj_idx].n_dia_sources, 2)
            input_dia_object = dia_collection.get_dia_object(
                dia_collection.dia_objects[obj_idx].id)
            self._compare_source_records(
                output_dia_objects[obj_idx].dia_object_record,
                input_dia_object.dia_object_record)

            output_src_cat = output_dia_objects[
                obj_idx].dia_source_catalog
            output_src_cat.sort(
                output_src_cat.getSchema().find('id').key)

            for record_a, record_b in zip(
                    output_src_cat,
                    input_dia_object.dia_source_catalog):
                self._compare_source_records(
                    record_a, record_b)

    def test_get_dia_object_records(self):
        """ Test the retrieval of SourceRecord objects representing the
        summarized DIAObjects from the database.
        """
        dia_objects = create_test_dia_objects(
            n_objects=5, n_sources=1, scatter_arcsec=0.0)
        dia_collection = DIAObjectCollection(dia_objects)
        self.assoc_db.store(dia_collection, True)

        self.assoc_db._commit()

        self.assoc_db._db_cursor.execute(
            "SELECT indexer_id FROM dia_objects")
        indexer_ids = np.array(
            self.assoc_db._db_cursor.fetchall(), np.int).flatten()

        dia_object_catalog = self.assoc_db._get_dia_object_records(indexer_ids)

        for obj_idx, dia_object_record in enumerate(dia_object_catalog):
            self._compare_source_records(
                dia_object_record,
                dia_collection.dia_objects[obj_idx].dia_object_record)

    def test_get_dia_sources(self):
        """ Test the retrieval of DIASources from the database.
        """
        dia_objects = create_test_dia_objects(
            n_objects=1, n_sources=5)
        dia_collection = DIAObjectCollection(dia_objects)

        self.assoc_db.store(dia_collection, True)

        self.assoc_db._commit()

        src_cat = self.assoc_db._get_dia_sources(0)
        for dia_source, created_source in zip(
                src_cat, dia_collection.dia_objects[0].dia_source_catalog):
            self._compare_source_records(dia_source, created_source)

    def test_store_dia_object_dia_source_pair(self):
        """ Test storing the ids of associated DIAObjects and DIASources.
        """
        for obj_id in range(2):
            for src_id in range(5):
                self.assoc_db._store_dia_object_source_pair(
                    obj_id, src_id + (obj_id * 5))
        self.assoc_db._commit()

        self.assoc_db._db_cursor.execute(
            "SELECT * FROM dia_objects_to_dia_sources")
        obj_id = -1
        for row_idx, row in enumerate(self.assoc_db._db_cursor.fetchall()):
            if row_idx % 5 == 0:
                obj_id += 1
            self.assertEqual(row, (row_idx, obj_id))

    def test_store_record_objects(self):
        """ Test storing a SourceRecord object in either the dia_objects and
        dia_sources table.
        """
        dia_objects = create_test_dia_objects(
            n_objects=1, n_sources=1, scatter_arcsec=0.0)
        dia_object_record = self._store_and_retrieve_source_record(
            dia_objects[0].dia_object_record,
            self.assoc_db._dia_object_converter)
        dia_source_record = self._store_and_retrieve_source_record(
            dia_objects[0].dia_source_catalog[0],
            self.assoc_db._dia_object_converter)

        self._compare_source_records(dia_object_record,
                                     dia_objects[0].dia_object_record)
        self._compare_source_records(
            dia_source_record, dia_objects[0].dia_source_catalog[0])

    def _store_and_retrieve_source_record(self,
                                          source_record,
                                          converter):
        """ Convenience method for round tripping a source record object.

        Parameters
        ----------
        source_record : lsst.afw.table.SourceRecord
            SourceRecord to store.
        converter : lsst.ap.association.SqliteDBConverter
            converter defining the table and schema to store.

        Return
        ------
        lsst.afw.table.SourceRecord
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
