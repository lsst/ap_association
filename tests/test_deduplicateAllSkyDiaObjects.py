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
import astropy.units
import numpy as np
import pandas as pd
import tempfile
import unittest
import yaml

from lsst.ap.association import DeduplicateAllSkyDiaObjectsTask
from lsst.ap.association.utils import getMidpointFromTimespan, readSchemaFromApdb
from lsst.dax.apdb import Apdb, ApdbSql, ApdbTables
from lsst.utils import getPackageDir
import lsst.utils.tests
from utils_tests import makeExposure, makeDiaObjects, makeDiaSources, makeDiaForcedSources, makeRegionTime, \
    getRegion


def _data_file_name(basename, module_name):
    """Return path name of a data file.

    Parameters
    ----------
    basename : `str`
        Name of the file to add to the path string.
    module_name : `str`
        Name of lsst stack package environment variable.

    Returns
    -------
    data_file_path : `str`
       Full path of the file to load from the "data" directory in a given
       repository.
    """
    return os.path.join(getPackageDir(module_name), "data", basename)


class TestDeduplicateAllSkyDiaObjects(unittest.TestCase):

    def setUp(self):
        # Create an instance of random generator with fixed seed.
        rng = np.random.default_rng(1234)

        self.db_file_fd, self.db_file = tempfile.mkstemp(
            dir=os.path.dirname(__file__))
        self.addCleanup(os.remove, self.db_file)
        self.addCleanup(os.close, self.db_file_fd)

        self.apdbConfig = ApdbSql.init_database(db_url="sqlite:///" + self.db_file)
        self.config_file = tempfile.NamedTemporaryFile()
        self.addCleanup(self.config_file.close)
        self.apdbConfig.save(self.config_file.name)
        self.apdb = Apdb.from_config(self.apdbConfig)
        self.schema = readSchemaFromApdb(self.apdb)

        self.exposure = makeExposure(False, False)
        self.regionTime = makeRegionTime(exposure=self.exposure)
        self.dateTime = getMidpointFromTimespan(self.regionTime.timespan)

        original_diaObjects = makeDiaObjects(50, self.exposure, rng)

        # create duplicates
        dup_0_0 = {"ra": original_diaObjects.iloc[0]['ra'],
                   "dec": original_diaObjects.iloc[0]['dec'] + 0.2/3600,
                   "diaObjectId": 1000,
                    "nDiaSources": 1}
        for f in ["u", "g", "r", "i", "z", "y"]:
            dup_0_0["%s_psfFluxNdata" % f] = 0
        dup_0_1 = {"ra": original_diaObjects.iloc[0]['ra'] + 0.1/3600,
                   "dec": original_diaObjects.iloc[0]['dec'] + 0.2/3600,
                   "diaObjectId": 1001,
                    "nDiaSources": 1}
        for f in ["u", "g", "r", "i", "z", "y"]:
            dup_0_1["%s_psfFluxNdata" % f] = 0
        dup_0_2 = {"ra": original_diaObjects.iloc[0]['ra'],
                   "dec": original_diaObjects.iloc[0]['dec'] - 0.1/3600,
                   "diaObjectId": 1002,
                    "nDiaSources": 1}
        for f in ["u", "g", "r", "i", "z", "y"]:
            dup_0_2["%s_psfFluxNdata" % f] = 0
        dup_0_3 = {"ra": original_diaObjects.iloc[0]['ra'] + 0.2/3600,
                   "dec": original_diaObjects.iloc[0]['dec'] - 0.1/3600,
                   "diaObjectId": 1003,
                    "nDiaSources": 1}
        for f in ["u", "g", "r", "i", "z", "y"]:
            dup_0_3["%s_psfFluxNdata" % f] = 0

            
        dup_1_0 = {"ra": original_diaObjects.iloc[1]['ra'],
                   "dec": original_diaObjects.iloc[1]['dec'] + 0.2/3600,
                   "diaObjectId": 2000,
                    "nDiaSources": 1}
        for f in ["u", "g", "r", "i", "z", "y"]:
            dup_1_0["%s_psfFluxNdata" % f] = 0
        dup_1_1 = {"ra": original_diaObjects.iloc[1]['ra'] + 0.1/3600,
                   "dec": original_diaObjects.iloc[1]['dec'] + 0.2/3600,
                   "diaObjectId": 2001,
                    "nDiaSources": 1}
        for f in ["u", "g", "r", "i", "z", "y"]:
            dup_1_1["%s_psfFluxNdata" % f] = 0
        dup_1_2 = {"ra": original_diaObjects.iloc[1]['ra'],
                   "dec": original_diaObjects.iloc[1]['dec'] - 0.1/3600,
                   "diaObjectId": 2002,
                    "nDiaSources": 1}
        for f in ["u", "g", "r", "i", "z", "y"]:
            dup_1_2["%s_psfFluxNdata" % f] = 0    


        dup_2_0 = {"ra": original_diaObjects.iloc[2]['ra'],
                   "dec": original_diaObjects.iloc[2]['dec'] + 0.2/3600,
                   "diaObjectId": 3000,
                    "nDiaSources": 1}
        for f in ["u", "g", "r", "i", "z", "y"]:
            dup_2_0["%s_psfFluxNdata" % f] = 0
        
        data = [dup_0_0, dup_0_1, dup_0_2, dup_0_3,
                dup_1_0, dup_1_1, dup_1_2,
                dup_2_0]
        duplicate_diaObjects = pd.DataFrame(data=data).set_index("diaObjectId", drop=False)

        self.diaObjects = pd.concat([original_diaObjects, duplicate_diaObjects])

        self.diaSources = makeDiaSources(
            100, self.diaObjects["diaObjectId"].to_numpy(), self.exposure, rng)
        self.diaForcedSources = makeDiaForcedSources(
            200, self.diaObjects["diaObjectId"].to_numpy(), self.exposure, rng)

        # Store the test diaSources as though they were observed a month before
        # the current exposure.
        dateTime = self.regionTime.timespan.begin.tai - 30 * astropy.units.day
        self.apdb.store(dateTime,
                        self.diaObjects,
                        self.diaSources,
                        self.diaForcedSources)

        # These columns are not in the DPDD, yet do appear in DiaSource.yaml.
        # We don't need to check them against the default APDB schema.
        self.ignoreColumns = ["band", "bboxSize", "isDipole", "flags"]

    def _makeConfig(self, **kwargs):
        config = DeduplicateAllSkyDiaObjectsTask.ConfigClass()
        config.apdb_config_url = self.config_file.name
        config.update(**kwargs)
        return config

    def testRun(self):
        """Test the full run method for the loader.
        """
        diaConfig = self._makeConfig()
        task = DeduplicateAllSkyDiaObjectsTask(config=diaConfig)
        result = task.run()

#        self.assertEqual(len(result.diaObjects), len(self.diaObjects))
#        self.assertEqual(len(result.diaSources), len(self.diaSources))
#        self.assertEqual(len(result.diaForcedSources),
#                         len(self.diaForcedSources))

    def test_remap_clusters(self):
        pass

    def test_cluster(self):
        pass    

    def test_count_duplicates(self):
        pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
