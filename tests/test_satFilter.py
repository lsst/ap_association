import random
import unittest

import lsst.utils.tests
import utils_tests
import lsst.meas.base.tests as measTests
import lsst.geom as geom
from lsst.ap.association.satFilter import (SatelliteFilterTask,
                                                        SatelliteFilterConfig)

class TestSatelliteFilterTask(unittest.TestCase):

    def setUp(self):
        # Generate 200 random pairs within a fake image. Only 100 of the pairs
        # will be used as the satellite filte, the rest will be the normal
        # 'sources'. This is for testing only right now, as we need this not to
        # be random during tests.

        # This needs to be in ra and dec

        import pydevd_pycharm
        pydevd_pycharm.settrace('localhost', port=8888, stdoutToServer=True,
                                stderrToServer=True)
        num_pairs = 200
        xmin, xmax = 0, 4176
        ymin, ymax = 0, 4176
        max_distance = 500
        psf = 5.0

        random_pairs = self._generate_random_pairs(num_pairs, xmin, xmax, ymin,
                                                   ymax, max_distance)

        self.nSources = 200
        self.yLoc = 100
        self.expId = 4321
        self.bbox = geom.Box2I(geom.Point2I(0, 0),
                               geom.Extent2I(4176, 4176))
        dataset = measTests.TestDataset(self.bbox)
        for srcIdx in range(self.nSources):
            dataset.addSource(10000.0, geom.Point2D(srcIdx, self.yLoc))
        schema = dataset.makeMinimalSchema()
        schema.addField("sky_source", type="Flag", doc="Sky objects.")
        schema.addField('base_SdssCentroid_x', type="F", doc="Trail extends off image")
        schema.addField('base_SdssCentroid_y', type="F", doc="Trail extends off image")
        _, self.diaSourceCat = dataset.realize(10.0, schema, randomSeed=1234)

        # I have more pairs than I need, I will only make "sources" out of the
        # first in each pair for now
        for i, pair in enumerate(random_pairs, 0):
            self.diaSourceCat[i]['base_SdssCentroid_x'] = pair[0][0]
            self.diaSourceCat[i]['base_SdssCentroid_y'] = pair[0][1]

        mjd = 57071.0
        self.utc_jd = mjd + 2_400_000.5 - 35.0 / (24.0 * 60.0 * 60.0)

        # self.visitInfo = afwImage.VisitInfo(
            # This incomplete visitInfo is sufficient for testing because the
            # Python constructor sets all other required values to some
            # default.
       #     exposureTime=30.0,
       #     darkTime=3.0,
       #     date=dafBase.DateTime(mjd, system=dafBase.DateTime.MJD),
       #     boresightRaDec=geom.SpherePoint(0.0, 0.0, geom.degrees),
       # )
    def _generate_random_pairs(self, num_pairs, xmin, xmax, ymin, ymax, max_distance):
        pairs = []
        for _ in range(num_pairs):
            x1 = random.randint(xmin, xmax)
            y1 = random.randint(ymin, ymax)

            # Generate x2 and y2 such that they are at most max_distance away from x1, y1
            x2 = random.randint(max(xmin, x1 - max_distance),
                                min(xmax, x1 + max_distance))
            y2 = random.randint(max(ymin, y1 - max_distance),
                                min(ymax, y1 + max_distance))

            pairs.append(((x1, y1), (x2, y2)))

        return pairs

    def test_filter(self):

        import pydevd_pycharm
        pydevd_pycharm.settrace('localhost', port=8888, stdoutToServer=True,
                                stderrToServer=True)

        SatelliteFilterTask
        satelliteFilterTask = SatelliteFilterTask()

        result = satelliteFilterTask.run(self.diaSourceCatdia_sources,
                                                self.psf, self.random_pairs[0:100])

def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()