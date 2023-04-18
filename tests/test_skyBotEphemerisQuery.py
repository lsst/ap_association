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
import lsst.daf.base as dafBase
import lsst.geom as geom
import lsst.afw.image as afwImage
import lsst.pipe.base as pipeBase
from lsst.utils import getPackageDir
import lsst.utils.tests


class TestSkyBotEphemerisQuery(unittest.TestCase):

    def setUp(self):
        super().setUp()

        # Explicit date calculation to avoid errors from misuse of time libraries.
        mjd = 57071.0
        self.utc_jd = mjd + 2_400_000.5 - 35.0 / (24.0 * 60.0 * 60.0)

        self.visitId = 42
        self.visitInfo = afwImage.VisitInfo(
            # Incomplete VisitInfo; Python constructor allows any value to
            # be defaulted.
            exposureTime=30.0,
            darkTime=3.0,
            date=dafBase.DateTime(mjd, system=dafBase.DateTime.MJD),
            boresightRaDec=geom.SpherePoint(0.0, 0.0, geom.degrees),
        )

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
            testOut = ephTask._skybotConeSearch(self.visitInfo.boresightRaDec,
                                                self.visitInfo.date.get(),
                                                1.7)
        testData = pd.read_parquet(
            os.path.join(getPackageDir("ap_association"),
                         "tests",
                         "data",
                         "testSSObjects.parq")
        )
        self.assertEqual(len(testData), len(testOut))
        self.assertTrue(np.all(np.equal(testOut["ssObjectId"], testData["ssObjectId"])))

    def test_skybotRun(self):
        """Test that the correct upload is requested.
        """
        task = ephQ.SkyBotEphemerisQueryTask()
        with patch.object(task, '_skybotConeSearch') as mockSearch:
            task.run([pipeBase.InMemoryDatasetHandle(self.visitInfo)], self.visitId)
            mockSearch.assert_called_once()
            self.assertEqual(len(mockSearch.call_args.args), 3)
            self.assertEqual(mockSearch.call_args.args[0], self.visitInfo.boresightRaDec)
            self.assertAlmostEqual(mockSearch.call_args.args[1], self.utc_jd)
            self.assertAlmostEqual(mockSearch.call_args.args[2], task.config.queryRadiusDegrees)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
