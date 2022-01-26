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

from unittest.mock import patch
import numpy as np
import os
import pandas as pd
import unittest

import lsst.ap.association.skyBotEphemerisQuery as ephQ
import lsst.geom as geom
import lsst.pipe.base as pipeBase
from lsst.utils import getPackageDir
import lsst.utils.tests


class TestSkyBotEphemerisQuery(unittest.TestCase):

    def test_skyBotConeSearch(self):
        """Test that our parsing of SkyBot return data succeeds and produces
        consistent hashed dataIds.
        """
        def requestReplace(input1, input2):
            """Junk wrapper for replacing the external internal call with an
            internel data load.
            """
            with open(os.path.join(getPackageDir("ap_association"),
                                   "tests",
                                   "data",
                                   "testSSObjects.txt"),
                      "r") as f:
                outputText = f.read()
            return pipeBase.Struct(text=outputText)
        with patch('lsst.ap.association.skyBotEphemerisQuery.requests.request',
                   new=requestReplace):
            ephTask = ephQ.SkyBotEphemerisQueryTask()
            testOut = ephTask._skybotConeSearch(geom.SpherePoint(0, 0, geom.degrees), 57071, 1.7)
        testData = pd.read_parquet(
            os.path.join(getPackageDir("ap_association"),
                         "tests",
                         "data",
                         "testSSObjects.parq")
        )
        self.assertEqual(len(testData), len(testOut))
        self.assertTrue(np.all(np.equal(testData["ssObjectId"], testData["ssObjectId"])))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
