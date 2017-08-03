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
from test_dia_object import create_test_dia_sources


def create_test_dia_objects(n_objects=5, n_src=5, scatter_arcsec=1.0,
                            start_id=0, start_angle_degrees=0.0):
    """ Create DIAObjects with a specified number of DIASources attached.

    Parameters
    ----------
    n_objects : int
        Number of DIAObjects to generate.
    n_src : int
        Number of DIASources to generate for each DIAObject.
    scatter_arcsec : float
        Scatter to add to the position of each DIASource.

    Returns
    -------
    A list of DIAObjects
    """
    output_dia_objects = []
    for obj_idx in range(n_objects):
        src_cat = create_test_dia_sources(n_src)
        for src_idx in range(n_src):
            edit_source_record(
                src_cat[src_idx], start_id + n_src * obj_idx + src_idx,
                start_angle_degrees + 0.1 * obj_idx,
                start_angle_degrees + 0.1 * obj_idx,
                scatter_arcsec)
        output_dia_objects.append(DIAObject(src_cat))
    return output_dia_objects


def edit_source_record(src, src_id, ra_degrees, dec_degrees, scatter_arcsec):
    """ Edit the center coordinate and id of a source record in place.

    Parameters
    ----------
    src : lsst.afw.table.SourceRecord
        Input SourceRecord to edit.
    id : int
        Integer value to set the record id to.
    ra_degrees : float
        RA value to set the record coordinate to.
    dec_degrees : float
        DEC value to set the record coordinate to.
    scatter_arcsec : float
        Arcsecond scatter to add to the position of the source record coord.
    """
    coord = Coord(afwGeom.Angle(ra_degrees, units=afwGeom.degrees),
                  afwGeom.Angle(dec_degrees, units=afwGeom.degrees))
    if scatter_arcsec > 0.0:
        coord.offset(
            afwGeom.Angle(np.random.rand() * 360, units=afwGeom.degrees),
            afwGeom.Angle(np.random.rand() * scatter_arcsec,
                          units=afwGeom.arcseconds))
    src.setCoord(coord)
    src['id'] = src_id


class TestDIAObjectCollection(unittest.TestCase):

    def test_init(self):
        """ Test that we can properly create a DIAObjectCollection from a list
        of DIAObjects.
        """
        obj_list = create_test_dia_objects(1, 1, 0.1)
        obj_collection = DIAObjectCollection(obj_list)

        self.assertTrue(obj_collection.is_updated)
        self.assertTrue(obj_collection.is_valid_tree)

        self.assertEqual(obj_collection.get_dia_object_ids()[0], 0)
        self.assertEqual(obj_collection.get_dia_object(0).get('id'), 0)

    def test_append_and_update(self):
        """ Test that we can add a new DIAObject to an existing
        DIAObjectCollection.
        """
        obj_list = create_test_dia_objects(1, 1, 0.1)
        obj_collection = DIAObjectCollection(obj_list)

        new_dia_obj = create_test_dia_objects(1, 1, 0.1, 1, 0.1)[0]
        obj_collection.append(new_dia_obj)
        self.assertFalse(obj_collection.is_updated)
        self.assertFalse(obj_collection.is_valid_tree)

        obj_collection.update_dia_objects()
        self.assertTrue(obj_collection.is_updated)
        obj_collection.update_spatial_tree()
        self.assertTrue(obj_collection.is_valid_tree)

        self.assertEqual(len(obj_collection.get_dia_object_ids()), 2)
        self.assertEqual(obj_collection.get_dia_object(1).get('id'), 1)

    def test_score_and_match(self):
        """ Test association between a set of sources and an existing
        DIAObjectCollection.

        This also tests that a DIASource that can't be associated within
        tolerance is addpended to the DIAObjectCollection as a new
        DIAObject.
        """
        obj_collection = DIAObjectCollection(
            create_test_dia_objects(4, 1, -1))
        src_cat = create_test_dia_sources(5)
        for src_idx, src in enumerate(src_cat):
            edit_source_record(
                src, src_idx + 4, 0.1 * src_idx, 0.1 * src_idx, -1)
        score_struct = obj_collection.score(
            src_cat,  afwGeom.Angle(1.0, units=afwGeom.arcseconds))

        self.assertFalse(np.isfinite(score_struct.scores[-1]))
        for src_idx in range(4):
            self.assertAlmostEqual(score_struct.scores[src_idx], 0.0,
                                   places=9)

        updated_indices = obj_collection.match(src_cat, score_struct)
        self.assertEqual(len(obj_collection.dia_objects), 5)

        for idx, obj_id in enumerate(obj_collection.get_dia_object_ids()):
            self.assertEqual(idx, updated_indices[idx])
            if idx != 4:
                self.assertEqual(
                    obj_collection.get_dia_object(obj_id).n_dia_sources, 2)
            else:
                self.assertEqual(
                    obj_collection.get_dia_object(obj_id).n_dia_sources, 1)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":

    lsst.utils.tests.init()
    unittest.main()
