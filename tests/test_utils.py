# This file is part of ap_association.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
import unittest

import pandas as pd

import lsst.daf.butler as dafButler

from lsst.ap.association.utils import readSdmSchemaFile, make_empty_catalog, convertTableToSdmSchema, \
    getMidpointFromTimespan, objID_to_ssObjectID, ssObjectID_to_objID
from utils_tests import makeExposure, makeRegionTime


class TestUtils(unittest.TestCase):

    def test_make_empty_catalog(self):
        """Check that an empty catalog has the correct format.
        """
        schemaFile = os.path.join("${SDM_SCHEMAS_DIR}", "yml", "apdb.yaml")
        schemaName = "ApdbSchema"
        schema = readSdmSchemaFile(schemaFile, schemaName)

        tableNames = ["DiaObject", "DiaSource", "DiaForcedSource"]
        for tableName in tableNames:
            emptyDiaObjects = make_empty_catalog(schema, tableName=tableName)
            self.assertTrue(emptyDiaObjects.empty)

            emptyColumns = set(emptyDiaObjects.columns)
            self.assertIn("ra", emptyColumns)
            self.assertIn("dec", emptyColumns)
            self.assertIn("diaObjectId", emptyColumns)

            emptyDf = pd.DataFrame(columns=["diaObjectId",])
            emptyDf.set_index("diaObjectId")
            convertedEmptyDiaObjects = convertTableToSdmSchema(schema, emptyDf, tableName=tableName)
            # TODO: we have no tests of convertTableToSdmSchema, so it's dangerous to use it as an oracle.
            emptyTypes = dict(zip(emptyDiaObjects.columns, emptyDiaObjects.dtypes))
            convertedEmptyTypes = dict(zip(convertedEmptyDiaObjects.columns, convertedEmptyDiaObjects.dtypes))
            self.assertEqual(emptyTypes, convertedEmptyTypes)

    def test_regionTime_timespan(self):
        """Check that the midpoint from a RegionTimeInfo matches the time from
        the visitInfo of the exposure.
        """
        exposure = makeExposure()
        regionTime = makeRegionTime(exposure)
        visitTime = exposure.visitInfo.date.toAstropy()
        midpoint = getMidpointFromTimespan(regionTime.timespan)
        self.assertEqual(visitTime.jd, midpoint.jd)

    def test_invalidTimespan(self):
        """Test the error handling of getMidpointFromTimespan.
        """
        exposure = makeExposure()
        visitTime = exposure.visitInfo.date.toAstropy()
        timespan_no_end = dafButler.Timespan(begin=visitTime, end=dafButler.Timespan.EMPTY)
        timespan_no_begin = dafButler.Timespan(begin=dafButler.Timespan.EMPTY, end=visitTime)
        with self.assertRaises(ValueError):
            getMidpointFromTimespan(timespan_no_end)
        with self.assertRaises(ValueError):
            getMidpointFromTimespan(timespan_no_begin)

        timespan_none_end = dafButler.Timespan(begin=visitTime, end=None)
        timespan_none_begin = dafButler.Timespan(begin=None, end=visitTime)
        with self.assertRaises(ValueError):
            getMidpointFromTimespan(timespan_none_end, allowUnbounded=False)
        with self.assertRaises(ValueError):
            getMidpointFromTimespan(timespan_none_begin, allowUnbounded=False)
        self.assertEqual(getMidpointFromTimespan(timespan_none_begin).jd, visitTime.jd)
        self.assertEqual(getMidpointFromTimespan(timespan_none_end).jd, visitTime.jd)

        timespan_none_both = dafButler.Timespan(begin=None, end=None)
        with self.assertRaises(ValueError):
            getMidpointFromTimespan(timespan_none_both, allowUnbounded=True)
        with self.assertRaises(ValueError):
            getMidpointFromTimespan(timespan_none_both, allowUnbounded=False)

    def test_ssObjectID_to_objID_and_objID_to_ssObjectID(self):
        """Convert between ssObjectIDs and MPC packed designations
        """
        allowed_strings = ['J95X00A', 'J95X01L', 'J95F13B', 'J98SA8Q', 'J98SC7V', 'J98SG2S'] \
            + ['K99AJ3Z', 'K08Aa0A', 'K07Tf8A', 'PLS2040', 'T1S3138', 'T2S1010', 'T3S4101'] \
            + ['       ', '\x00\x00\x00\x00\x00\x00\x00']
        allowed_flags = [i for i in range(0, 256)]
        allowed_ssObjectIDs = [0, 1 << 64 - 1] + [1 << n for n in range(64)]
        for allowed_string in allowed_strings:
            for allowed_flag in allowed_flags:
                returned_string, returned_flag = ssObjectID_to_objID(
                    objID_to_ssObjectID(allowed_string, allowed_flag))
                self.assertEqual((allowed_string, allowed_flag), (returned_string, returned_flag))
        for allowed_ssObjectID in allowed_ssObjectIDs:
            returned_ssObjectID = objID_to_ssObjectID(*ssObjectID_to_objID(allowed_ssObjectID))
            self.assertEqual(allowed_ssObjectID, returned_ssObjectID)

    def test_invalid_ssObjectID_to_objID_and_objID_to_ssObjectID(self):
        """Convert between ssObjectIDs and MPC packed designations
        """
        allowed_strings = ['J95X00A', 'J95X01L', 'J95F13B', 'J98SA8Q', 'J98SC7V', 'J98SG2S'] \
            + ['K99AJ3Z', 'K08Aa0A', 'K07Tf8A', 'PLS2040', 'T1S3138', 'T2S1010', 'T3S4101']
        allowed_flags = [i for i in range(0, 256)]
        disallowed_flags = [-999999999, -512, -256, -255, -1, 256, 512, 99999999]
        disallowed_strings = [''] + [ch for ch in 'ABCDEFGHIJKMNOPQRSTUVWXYZ0123456789 -'] \
            + ['A' * i for i in range(2, 7)] + ['Z' * i for i in range(2, 7)] \
            + ['Ā', '🔭', 'A' * 8, ' ' * 8, 'A' * 128]
        disallowed_ssObjectIDs = [-1, 1 << 64, 1 << 64 + 1, 2 << 65]
        for allowed_string in allowed_strings:
            for disallowed_flag in disallowed_flags:
                with self.assertRaises(ValueError):
                    objID_to_ssObjectID(allowed_string, disallowed_flag)
        for disallowed_string in disallowed_strings:
            for allowed_flag in allowed_flags:
                with self.assertRaises(ValueError):
                    objID_to_ssObjectID(disallowed_string, allowed_flag)
        for disallowed_string in disallowed_strings:
            for disallowed_flag in disallowed_flags:
                with self.assertRaises(ValueError):
                    objID_to_ssObjectID(disallowed_string, disallowed_flag)
        for disallowed_ssObjectID in disallowed_ssObjectIDs:
            with self.assertRaises(ValueError):
                ssObjectID_to_objID(disallowed_ssObjectID)
