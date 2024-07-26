import random
import unittest
import numpy as np
import tempfile
import os

import lsst.utils.tests
from utils_tests import makeExposure, makeDiaObjects, makeDiaSources
from lsst.ap.association.satFilter import (SatelliteFilterTask,
                                           SatelliteFilterConfig)
from lsst.ap.association.utils import readSchemaFromApdb
from lsst.dax.apdb import Apdb, ApdbSql, ApdbTables


class TestSatelliteFilterTask(unittest.TestCase):

    def setUp(self):
        # Generate 200 random pairs within a fake image. Only 100 of the pairs
        # will be used as the satellite filte, the rest will be the normal
        # 'sources'. This is for testing only right now, as we need this not to
        # be random during tests.
        self.psf = 0.4
        rng = np.random.default_rng(1234)
        self.db_file_fd, self.db_file = tempfile.mkstemp(
            dir=os.path.dirname(__file__))
        self.addCleanup(os.remove, self.db_file)
        self.addCleanup(os.close, self.db_file_fd)

        self.apdbConfig = ApdbSql.init_database(
            db_url="sqlite:///" + self.db_file)
        self.apdb = Apdb.from_config(self.apdbConfig)
        self.schema = readSchemaFromApdb(self.apdb)

        self.exposure = makeExposure(False, False)

        self.diaObjects = makeDiaObjects(20, self.exposure, rng)
        self.diaSources = makeDiaSources(
            100, self.diaObjects["diaObjectId"].to_numpy(), self.exposure, rng)

        self.random_pairs = self._generate_random_pairs(5,0,100,0,90,20)

        self.sat_coords = []
        for i in range(5):
            self.sat_coords.append((np.array([self.diaSources['ra'][i]+1, self.diaSources['dec'][i]+1]),
                                    np.array([self.diaSources['ra'][i]-1, self.diaSources['dec'][i]-1])))


    def test_filter(self):
        satelliteFilterTask = SatelliteFilterTask()

        result = satelliteFilterTask.run( self.diaSources,
                                         self.psf, np.array(self.sat_coords))

    def _generate_random_pairs(self, num_pairs, xmin, xmax, ymin, ymax, max_distance):
        pairs = []
        for _ in range(num_pairs):
            x1 = random.randint(xmin, xmax)
            y1 = random.randint(ymin, ymax)

            x2 = random.randint(max(xmin, x1 - max_distance),
                                min(xmax, x1 + max_distance))
            y2 = random.randint(max(ymin, y1 - max_distance),
                                min(ymax, y1 + max_distance))

            pairs.append((np.array([x1, y1]), np.array([x2, y2])))

        return np.array(pairs)


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()