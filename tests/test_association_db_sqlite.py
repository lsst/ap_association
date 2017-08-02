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
        pass

    def test_store(self):
        """ Test storing of a DIACollection.
        """
        dia_objects = create_test_dia_objects(n_objects=1, n_src=1)
        dia_collection = DIAObjectCollection(dia_objects)
        assoc_db = AssociationDBSqliteTask()
        assoc_db.create_tables()

        assoc_db.store(dia_collection, [0])

        assoc_db.commit()
        assoc_db.close()

    def test_get_dia_objects(self):
        pass

    def test_get_dia_sources(self):
        pass

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
