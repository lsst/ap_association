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
import tempfile
import unittest
import yaml

from lsst.ap.association import LoadDiaCatalogsTask
from lsst.ap.association.utils import getMidpointFromTimespan, readSchemaFromApdb
from lsst.dax.apdb import Apdb, ApdbSql, ApdbTables
from lsst.resources import ResourcePath
import lsst.utils.tests
from utils_tests import makeExposure, makeDiaObjects, makeDiaSources, makeDiaForcedSources, makeRegionTime, \
    getRegion


class TestLoadDiaCatalogs(unittest.TestCase):

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

        self.diaObjects = makeDiaObjects(20, self.exposure, rng)
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
        config = LoadDiaCatalogsTask.ConfigClass()
        config.apdb_config_url = self.config_file.name
        config.update(**kwargs)
        return config

    def testRun(self):
        """Test the full run method for the loader.
        """
        diaConfig = self._makeConfig()
        diaLoader = LoadDiaCatalogsTask(config=diaConfig)
        result = diaLoader.run(self.regionTime)

        self.assertEqual(len(result.diaObjects), len(self.diaObjects))
        self.assertEqual(len(result.diaSources), len(self.diaSources))
        self.assertEqual(len(result.diaForcedSources),
                         len(self.diaForcedSources))

    def testLoadDiaObjects(self):
        """Test that the correct number of diaObjects are loaded.
        """
        diaConfig = self._makeConfig()
        diaLoader = LoadDiaCatalogsTask(config=diaConfig)
        region = getRegion(self.exposure)
        diaObjects = diaLoader.loadDiaObjects(region,
                                              self.schema)
        self.assertEqual(len(diaObjects), len(self.diaObjects))

    def testLoadDiaForcedSources(self):
        """Test that the correct number of diaForcedSources are loaded.
        """
        diaConfig = self._makeConfig()
        diaLoader = LoadDiaCatalogsTask(config=diaConfig)
        region = getRegion(self.exposure)
        diaForcedSources = diaLoader.loadDiaForcedSources(
            self.diaObjects,
            region,
            self.dateTime,
            self.schema)
        self.assertEqual(len(diaForcedSources), len(self.diaForcedSources))

    def testLoadDiaSources(self):
        """Test that the correct number of diaSources are loaded.

        Also check that they can be properly loaded both by location and
        ``diaObjectId``.
        """
        diaConfig = self._makeConfig()
        diaLoader = LoadDiaCatalogsTask(config=diaConfig)

        region = getRegion(self.exposure)
        diaSources = diaLoader.loadDiaSources(self.diaObjects,
                                              region,
                                              self.dateTime,
                                              self.schema)
        self.assertEqual(len(diaSources), len(self.diaSources))

    def test_apdbSchema(self):
        """Test that the default DiaSource schema from dax_apdb agrees with the
        column names defined here in ap_association/data/DiaSource.yaml.
        """
        tableDef = self.apdb.tableDef(ApdbTables.DiaSource)
        apdbSchemaColumns = [column.name for column in tableDef.columns]

        functorFile = ResourcePath("resource://lsst.ap.association/resources/data/DiaSource.yaml")
        with functorFile.open("r") as yaml_stream:
            diaSourceFunctor = yaml.safe_load_all(yaml_stream)
            for functor in diaSourceFunctor:
                diaSourceColumns = [column for column in list(functor['funcs'].keys())
                                    if column not in self.ignoreColumns]
            self.assertLess(set(diaSourceColumns), set(apdbSchemaColumns))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
