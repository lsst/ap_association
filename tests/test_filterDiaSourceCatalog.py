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

import unittest

from lsst.ap.association.filterDiaSourceCatalog import (FilterDiaSourceCatalogConfig,
                                                        FilterDiaSourceCatalogTask)
import lsst.geom as geom
import lsst.meas.base.tests as measTests
import lsst.utils.tests


class TestFilterDiaSourceCatalogTask(unittest.TestCase):

    def setUp(self):
        self.nSources = 10
        self.nSkySources = 5
        self.yLoc = 100
        self.expId = 4321
        self.bbox = geom.Box2I(geom.Point2I(0, 0),
                               geom.Extent2I(1024, 1153))
        dataset = measTests.TestDataset(self.bbox)
        for srcIdx in range(self.nSources):
            dataset.addSource(10000.0, geom.Point2D(srcIdx, self.yLoc))
        schema = dataset.makeMinimalSchema()
        schema.addField("sky_source", type="Flag", doc="Sky objects.")
        _, self.diaSourceCat = dataset.realize(10.0, schema, randomSeed=1234)
        self.diaSourceCat[0:self.nSkySources]["sky_source"] = True
        self.config = FilterDiaSourceCatalogConfig()

    def test_run_without_filter(self):
        self.config.doRemoveSkySources = False
        self.config.doWriteRejectedSources = True
        filterDiaSourceCatalogTask = FilterDiaSourceCatalogTask(config=self.config)
        result = filterDiaSourceCatalogTask.run(self.diaSourceCat)
        self.assertEqual(len(result.filteredDiaSourceCat), len(self.diaSourceCat))
        self.assertEqual(len(result.rejectedDiaSources), 0)
        self.assertEqual(len(self.diaSourceCat), self.nSources)

    def test_run_with_filter(self):
        self.config.doRemoveSkySources = True
        self.config.doWriteRejectedSources = True
        filterDiaSourceCatalogTask = FilterDiaSourceCatalogTask(config=self.config)
        result = filterDiaSourceCatalogTask.run(self.diaSourceCat)
        nExpectedFilteredSources = self.nSources - self.nSkySources
        self.assertEqual(len(result.filteredDiaSourceCat),
                         len(self.diaSourceCat[~self.diaSourceCat['sky_source']]))
        self.assertEqual(len(result.filteredDiaSourceCat), nExpectedFilteredSources)
        self.assertEqual(len(result.rejectedDiaSources), self.nSkySources)
        self.assertEqual(len(self.diaSourceCat), self.nSources)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
