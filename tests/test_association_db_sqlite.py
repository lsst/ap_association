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

from lsst.ap.association import *
from lsst.afw.coord import Coord
import lsst.afw.geom as afwGeom
import lsst.utils.tests
from test_dia_collection import create_test_dia_objects
from test_dia_object import create_test_dia_sources

class TestAssociationDBSqlite(unittest.TestCase):

    def setup(self):
        pass

    def tearDown(self):
        pass

    def test_init(self):
        """ Test instatiation of the task.
        """
        assoc_db = AssociationDBSqliteTask()
        assoc_db.close()

    def test_create_tables(self):
        """ Test the creation of an empty table.
        """
        assoc_db = AssociationDBSqliteTask()
        assoc_db.create_tables()
        assoc_db.commit()
        assoc_db.close()

    def test_load(self):
        """ Test loading of DIAObjects
        """
        dia_objects = create_test_dia_objects(
            n_objects=5, n_src=2, increment_degrees=0.04)
        dia_collection = DIAObjectCollection(dia_objects)
        assoc_db = AssociationDBSqliteTask()
        assoc_db.create_tables()

        assoc_db.store(dia_collection, True)

        ctr_point = Coord(
            afwGeom.Angle(0.08, units=afwGeom.degrees), 
            afwGeom.Angle(0.08, units=afwGeom.degrees))
        output_dia_collection = assoc_db.load(
            ctr_point, afwGeom.Angle(0.2))

        for obj_idx in xrange(5):
            self.assertEqual(
                output_dia_collection.dia_objects[obj_idx].n_dia_sources, 2)
            obj_id = output_dia_collection.dia_objects[obj_idx].get('id')
            self.assertAlmostEqual(
                output_dia_collection.dia_objects[obj_idx].ra.asDegrees(),
                obj_id / 2 * 0.04)
            self.assertAlmostEqual(
                output_dia_collection.dia_objects[obj_idx].dec.asDegrees(),
                obj_id / 2 * 0.04)
            for src_idx, src_record in enumerate(
                    output_dia_collection.dia_objects[
                        obj_idx].dia_source_catalog):
                self.assertEqual(src_record.getId(), src_idx +  obj_id)
                self.assertAlmostEqual(
                    src_record.getRa().asDegrees(), obj_id / 2 * 0.04)
                self.assertAlmostEqual(
                    src_record.getDec().asDegrees(), obj_id / 2 * 0.04)

        assoc_db.close()

    def test_store(self):
        """ Test storing of a DIACollection.
        """
        dia_objects = create_test_dia_objects(n_objects=1, n_src=1)
        dia_collection = DIAObjectCollection(dia_objects)
        assoc_db = AssociationDBSqliteTask()
        assoc_db.create_tables()

        assoc_db.store(dia_collection, True)

        assoc_db.commit()
        assoc_db.close()

    def test_store_update(self):
        """ Test the retrieval of DIASources from the database.
        """
        dia_objects = create_test_dia_objects(n_objects=1, n_src=1)
        dia_collection = DIAObjectCollection(dia_objects)
        assoc_db = AssociationDBSqliteTask()
        assoc_db.create_tables()

        assoc_db.store(dia_collection, True)
        assoc_db.commit()

        new_dia_object = create_test_dia_objects(
            n_objects=1, n_src=1, start_id=1)
        dia_collection.append(new_dia_object[0])
        dia_collection.update_dia_objects()
        dia_collection.update_spatial_tree()

        assoc_db.store_updated(dia_collection, [1])

        assoc_db.close()

    def test_get_dia_objects(self):
        """ Test the retrieval of DIAObjects from the database.
        """
        dia_objects = create_test_dia_objects(
            n_objects=2, n_src=2, scatter_arcsec=0.0)
        dia_collection = DIAObjectCollection(dia_objects)
        assoc_db = AssociationDBSqliteTask()
        assoc_db.create_tables()
        assoc_db.store(dia_collection, True)

        assoc_db.commit()

        assoc_db._db_cursor.execute(
            "SELECT indexer_id FROM dia_objects")
        indexer_ids = np.array(
            assoc_db._db_cursor.fetchall(), np.int).flatten()

        output_dia_collection = assoc_db.get_dia_objects(indexer_ids)

        for obj_idx in xrange(2):
            self.assertEqual(
                output_dia_collection.dia_objects[obj_idx].n_dia_sources, 2)
            obj_id = output_dia_collection.dia_objects[obj_idx].get('id')
            self.assertAlmostEqual(
                output_dia_collection.dia_objects[obj_idx].ra.asDegrees(),
                obj_id / 2 * 0.1)
            self.assertAlmostEqual(
                output_dia_collection.dia_objects[obj_idx].dec.asDegrees(),
                obj_id / 2 * 0.1)
            for src_idx, src_record in enumerate(
                    output_dia_collection.dia_objects[
                        obj_idx].dia_source_catalog):
                self.assertEqual(src_record.getId(), src_idx +  obj_id)
                self.assertAlmostEqual(
                    src_record.getRa().asDegrees(), obj_id / 2 * 0.1)
                self.assertAlmostEqual(
                    src_record.getDec().asDegrees(), obj_id / 2 * 0.1)

        assoc_db.close()

    def test_get_dia_object_records(self):
        """ Test the retrieval of SourceRecord objects representing the
        summarized DIAObjects from the database.
        """
        dia_objects = create_test_dia_objects(
            n_objects=5, n_src=1, scatter_arcsec=0.0)
        dia_collection = DIAObjectCollection(dia_objects)
        assoc_db = AssociationDBSqliteTask()
        assoc_db.create_tables()
        assoc_db.store(dia_collection, True)

        assoc_db.commit()

        assoc_db._db_cursor.execute(
            "SELECT indexer_id FROM dia_objects")
        indexer_ids = np.array(
            assoc_db._db_cursor.fetchall(), np.int).flatten()

        assoc_db.get_dia_object_records(indexer_ids)

        for obj_idx in xrange(2):
            self.assertEqual(
                dia_collection.dia_objects[obj_idx].n_dia_sources, 2)
            for src_id, src_record in enumerate(
                    dia_collection.dia_objects[obj_idx].dia_source_catalog):
                self.assertEqual(src_record.getId(), src_id + obj_idx * 2)

        assoc_db.close()

    def test_get_dia_sources(self):
        """ Test the retrieval of DIASources from the database.
        """
        dia_objects = create_test_dia_objects(n_objects=1, n_src=5)
        dia_collection = DIAObjectCollection(dia_objects)

        assoc_db = AssociationDBSqliteTask()
        assoc_db.create_tables()
        assoc_db.store(dia_collection, True)

        assoc_db.commit()
        src_cat = assoc_db.get_dia_sources([0, 1, 2, 3, 4])
        for src_idx, src in enumerate(src_cat):
            self.assertEqual(src['id'], src_idx)

        assoc_db.close()

    def test_store_dia_object_dia_source_pair(self):
        """ Test storing a DIAObject.
        """
        dia_objects = create_test_dia_objects(n_objects=1, n_src=1)
        assoc_db = AssociationDBSqliteTask()
        assoc_db.create_tables()

        for obj_id in range(2):
            for src_id in range(5):
                assoc_db.store_dia_object_source_pair(
                    obj_id, src_id)

        assoc_db.commit()
        assoc_db.close()

    def test_store_dia_object(self):
        """ Test storing a DIAObject.
        """
        dia_objects = create_test_dia_objects(n_objects=1, n_src=1)
        assoc_db = AssociationDBSqliteTask()
        assoc_db.create_tables()

        assoc_db.store_dia_object(dia_objects[0])

        assoc_db.commit()
        assoc_db.close()

    def test_store_dia_sources(self):
        """ Test storing a Source
        """
        dia_sources = create_test_dia_sources(1)
        assoc_db = AssociationDBSqliteTask()
        assoc_db.create_tables()

        assoc_db.store_dia_source(dia_sources[0])

        assoc_db.commit()
        assoc_db.close()


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":

    lsst.utils.tests.init()
    unittest.main()
