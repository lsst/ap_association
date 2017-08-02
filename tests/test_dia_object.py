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

import unittest

from lsst.ap.association import DIAObject
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.utils.tests


def create_test_dia_sources(n_sources=5):
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
        src['id'] = src_idx
        src['coord_ra'] = afwGeom.Angle(0.0 + 1. * src_idx,
                                        units=afwGeom.degrees)
        src['coord_dec'] = afwGeom.Angle(0.0 + 1. * src_idx,
                                         units=afwGeom.degrees)
        # Add a flux at some point

    return sources


class TestDIAObject(unittest.TestCase):

    def test_init(self):
        """ Test DIAObject creation and if we can properly instantiate a
        DIAObject from an already create dia_object_record.
        """
        single_source = create_test_dia_sources(1)

        dia_obj = DIAObject(single_source, None)
        dia_obj_dup = DIAObject(single_source, dia_obj.dia_object_record)

        compare_values = {
            "coord_ra": 0.0,
            "coord_dec": 0.0
        }

        self._compare_dia_object_values(dia_obj, 0, compare_values)
        self._compare_dia_object_values(dia_obj_dup, dia_obj.id,
                                        compare_values)

    def _compare_dia_object_values(self, dia_object, expected_id,
                                   expected_dict):
        """ Compare values computed in the compute_summary_statistics
        DIAObject class method with those expected.

        Parameters
        ----------
        dia_object : lsst.ap.association.DIAObject
            Input DIAObect to test
        expected_id : int
            Expected id of the DIAObject
        expected_dict : dict
            Dictionary of field name and expected value
        """

        self.assertEqual(dia_object.get('id'), expected_id)
        for key in expected_dict.keys():
            if isinstance(dia_object.get(key), afwGeom.Angle):
                self.assertAlmostEqual(dia_object.get(key).asDegrees(),
                                       expected_dict[key])
            else:
                self.assertAlmostEqual(dia_object.get(key),
                                       expected_dict[key])

    def test_update(self):
        """ If we instantiate the class with a set of DIASources we test to
        make sure that the DIAObject is intantiated correctly and that we
        compute the summary statistics in the expected way.
        """
        sources = create_test_dia_sources(5)
        dia_obj = DIAObject(sources, None)

        compare_values = {
            "coord_ra": 1.9987807133764057,
            "coord_dec": 2.000608742419802
        }

        self._compare_dia_object_values(dia_obj, 0, compare_values)

    def test_dia_source_append_and_update(self):
        """ Test the appending of a DIASource to a DIAObject. We also
        test that the update function works properly and that the summary
        statistics are computed as expected.
        """
        single_source = create_test_dia_sources(1)
        dia_obj = DIAObject(single_source, None)
        self.assertEqual(dia_obj.id, single_source[0].getId())
        self.assertTrue(dia_obj.is_updated)

        sources = create_test_dia_sources(2)
        dia_obj.append_dia_source(sources[1])
        self.assertFalse(dia_obj.is_updated)

        associated_sources = dia_obj.dia_source_catalog
        self.assertEqual(len(associated_sources), 2)
        self.assertEqual(associated_sources[-1].getId(),
                         sources[-1].getId())
        self.assertEqual(associated_sources[-1].getCoord(),
                         sources[-1].getCoord())

        compare_values = {
            "coord_ra": 0.4999619199226218,
            "coord_dec": 0.5000190382261059
        }

        dia_obj.update()
        self._compare_dia_object_values(dia_obj, 0, compare_values)

    def test_compute_light_curve(self):
        """ Not implemented yet.
        """
        pass


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":

    lsst.utils.tests.init()
    unittest.main()
