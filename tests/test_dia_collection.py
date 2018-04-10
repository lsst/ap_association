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

from lsst.ap.association import DIAObject, DIAObjectCollection
import lsst.afw.geom as afwGeom
import lsst.utils.tests
from test_dia_object import create_test_dia_sources


def create_test_dia_objects(n_objects=5, n_sources=5, start_id=0,
                            start_angle_degrees=0.0, increment_degrees=0.1,
                            scatter_arcsec=1.0):
    """ Create DIAObjects with a specified number of DIASources attached.

    Parameters
    ----------
    n_objects : int
        Number of DIAObjects to generate.
    n_src : int
        Number of DIASources to generate for each DIAObject.
    start_id : int
        Starting index to increment the created DIAObjects from.
    start_angle_degrees : float
        Starting position of the objects. Additional objects are
        incremented from this position by 0.1 degrees. The position
        of the first object will be RA=start_angle_degrees,
        DEC=start_angle_degreesi
    increment_degrees : float
        Amount to increment RA and DEC by for each new DIAObject
    scatter_arcsec : float
        Scatter to add to the position of each DIASource.

    Returns
    -------
    `list` of `lsst.ap.association.DIAObjects`
    """
    output_dia_objects = []
    for obj_idx in range(n_objects):
        src_cat = create_test_dia_sources(n_sources)
        for src_idx in range(n_sources):
            edit_and_offset_source_record(
                src_cat[src_idx],
                start_id + n_sources * obj_idx + src_idx,
                start_angle_degrees + increment_degrees * obj_idx,
                start_angle_degrees + increment_degrees * obj_idx,
                scatter_arcsec)
        output_dia_objects.append(DIAObject(src_cat))
    return output_dia_objects


def edit_and_offset_source_record(src, src_id, ra_degrees, dec_degrees,
                                  scatter_arcsec):
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
    coord = lsst.afw.geom.SpherePoint(ra_degrees, dec_degrees, afwGeom.degrees)
    if scatter_arcsec > 0.0:
        coord = coord.offset(
            np.random.rand() * 360 * afwGeom.degrees,
            np.random.rand() * scatter_arcsec * afwGeom.arcseconds)
    src.setCoord(coord)
    src['id'] = src_id


class TestDIAObjectCollection(unittest.TestCase):

    def test_init(self):
        """ Test that we can properly create a DIAObjectCollection from a list
        of DIAObjects.
        """
        obj_list = create_test_dia_objects(n_objects=1, n_sources=1,
                                           start_id=0)
        obj_collection = DIAObjectCollection(obj_list)

        self.assertTrue(obj_collection.is_updated)
        self.assertTrue(obj_collection.is_valid_tree)

        self.assertEqual(obj_collection.get_dia_object_ids(), [0, ])
        self.assertEqual(obj_collection.get_dia_object(0).get('id'), 0)

    def test_append_and_update(self):
        """ Test that we can add a new DIAObject to an existing
        DIAObjectCollection.
        """
        obj_list = create_test_dia_objects(n_objects=1, n_sources=1, start_id=0)
        obj_collection = DIAObjectCollection(obj_list)

        new_dia_obj = create_test_dia_objects(n_objects=1,
                                              n_sources=1,
                                              start_id=1,
                                              start_angle_degrees=0.1)[0]
        obj_collection.append(new_dia_obj)
        self.assertFalse(obj_collection.is_updated)
        self.assertFalse(obj_collection.is_valid_tree)

        obj_collection.update_dia_objects()
        self.assertTrue(obj_collection.is_updated)
        self.assertFalse(obj_collection.is_valid_tree)

        obj_collection.update_spatial_tree()
        self.assertTrue(obj_collection.is_updated)
        self.assertTrue(obj_collection.is_valid_tree)

        self.assertEqual(obj_collection.get_dia_object_ids(), [0, 1])

    def test_empty_dia_collection(self):
        """ Test the creation and appending to a empty DIAObjectCollection.
        """
        dia_collection = DIAObjectCollection([])
        src_cat = create_test_dia_sources(5)
        for src_idx, src in enumerate(src_cat):
            edit_and_offset_source_record(
                src,
                src_idx,
                0.1 * src_idx,
                0.1 * src_idx,
                -1)

        score_struct = dia_collection.score(
            src_cat, afwGeom.Angle(1.0, units=afwGeom.arcseconds))
        for src_idx in range(5):
            # Our scores should be extremely close to 0 but not exactly so due
            # to machine noise.
            self.assertTrue(np.isinf(score_struct.scores[src_idx]))

        match_result = dia_collection.match(src_cat, score_struct)
        updated_ids = match_result.updated_and_new_dia_object_ids
        self.assertEqual(len(dia_collection.dia_objects), 5)
        self.assertEqual(match_result.n_updated_dia_objects, 0)
        self.assertEqual(match_result.n_new_dia_objects, 5)
        self.assertEqual(match_result.n_unassociated_dia_objects, 0)

        for idx, obj_id in enumerate(dia_collection.get_dia_object_ids()):
            self.assertEqual(obj_id, updated_ids[idx])
            # We created a new DIAObject in the collection hence the last
            # DIAObject in this collection is new and contains only one
            # DIASource.
            tmp_dia_obj = dia_collection.get_dia_object(obj_id)
            self.assertEqual(tmp_dia_obj.n_dia_sources, 1)
            self.assertEqual(tmp_dia_obj.id,
                             tmp_dia_obj.dia_source_catalog[-1].getId())


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
