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
import unittest

import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
from lsst.afw.coord import Coord
import lsst.daf.base as dafBase
from lsst.ap.association import \
    AssociationDBSqliteTask, \
    AssociationDBSqliteConfig, \
    AssociationTask, \
    AssociationConfig, \
    DIAObjectCollection, \
    make_minimal_dia_source_schema
import lsst.utils.tests
from test_dia_collection import create_test_dia_objects


def create_test_dia_sources(n_sources=5, start_id=0, source_locs_deg=[0.0],
                            scatter_arcsec=1.0):
    """ Create dummy DIASources for use in our tests.

    Parameters
    ----------
    n_sources : int (optional)
        Number of fake sources to create for testing.

    Returns
    -------
    A lsst.afw.SourceCatalog
    """
    sources = afwTable.SourceCatalog(make_minimal_dia_source_schema())

    for src_idx in range(n_sources):
        src = sources.addNew()
        src['id'] = src_idx + start_id
        coord = Coord(source_locs_deg[src_idx] * afwGeom.degrees,
                      source_locs_deg[src_idx] * afwGeom.degrees)
        if scatter_arcsec > 0.0:
            coord.offset(
                np.random.rand() * 360 * afwGeom.degrees,
                np.random.rand() * scatter_arcsec * afwGeom.arcseconds)
        src.setCoord(coord)

    return sources


class TestAssociationTask(unittest.TestCase):

    def setUp(self):
        """ Create a sqlite3 database with default tables and schemas.
        """
        self.db_file = \
            os.path.join(os.path.dirname(__file__), '/tmp_db.sqlite3')
        assoc_db_config = AssociationDBSqliteConfig()
        assoc_db_config.db_name = self.db_file
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

        self.assoc_task_config = AssociationConfig()
        self.assoc_task = AssociationTask()
    
    def tearDown(self):
        """ Delete the database after we are done with it.
        """
        pass
        os.remove(self.db_file)
        del self.db_file
        del self.metadata

    def test_run(self):
        """ Test the run method with a database that already exists and
        contains DIAObjects and Sources.
        """
        pass

    def test_run_no_existing_objects(self):
        """ Test the run method with a completely empty database.
        """
        pass

    def _association_tester(self, create_objects=False):
        """ Convienience method for testing the Association run method.
        """
        pass

    def test_associate_sources(self):
        """ Test if performance of associate_sources method in
        AssociationTask.
        """
        n_objects = 5
        n_sources_per_object = 2
        dia_objects = create_test_dia_objects(
            n_objects=n_objects,
            n_sources=n_sources_per_object,
            start_angle_degrees=0.0,
            increment_degrees=0.04,
            scatter_arcsec=0.0)
        dia_collection = DIAObjectCollection(dia_objects)

        dia_sources = create_test_dia_sources(
            n_sources=9,
            start_id=10,
            source_locs_deg=[0.04 * (src_idx + 1) for src_idx in range(9)],
            scatter_arcsec=1.0)

        import pdb; pdb.set_trace()

        assoc_task = AssociationTask()
        assoc_result = assoc_task.associate_sources(dia_collection,
                                                    dia_sources)

        not_updated_idx = 0
        updated_idx_start = 1
        new_idx_start = 5
        self.assertEqual(len(assoc_result.dia_collection.dia_objects),
                         10)
        for obj_idx, output_dia_object in \
                enumerate(assoc_result.dia_collection.dia_objects):
            if obj_idx == not_updated_idx:
                self.assertEqual(output_dia_object.n_dia_sources, 2)
            elif updated_idx_start <= obj_idx < new_idx_start:
                self.assertEqual(output_dia_object.n_dia_sources, 3)
                self.assertEqual(output_dia_object.dia_source_catalog[-1].getId(),
                                 obj_idx - 1 + 10)
            else:
                self.assertEqual(output_dia_object.n_dia_sources, 1)
                self.assertEqual(output_dia_object.dia_source_catalog[-1].getId(),
                                 obj_idx - 1 + 10)
                self.assertEqual(output_dia_object.id, obj_idx + 5 + 4)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
