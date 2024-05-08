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
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase


class TestFilterDiaSourceCatalogTask(unittest.TestCase):

    def setUp(self):
        self.nSources = 15
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
        schema.addField('ext_trailedSources_Naive_flag_off_image', type="Flag",
                        doc="Trail extends off image")
        schema.addField('ext_trailedSources_Naive_flag_suspect_long_trail',
                        type="Flag", doc="Trail length is greater than three times the psf radius")
        schema.addField('ext_trailedSources_Naive_flag_edge', type="Flag",
                        doc="Trail contains edge pixels")
        schema.addField('ext_trailedSources_Naive_flag_nan', type="Flag",
                        doc="One or more trail coordinates are missing")
        schema.addField('ext_trailedSources_Naive_length', type="F",
                        doc="Length of the source trail")
        _, self.diaSourceCat = dataset.realize(10.0, schema, randomSeed=1234)
        self.diaSourceCat[0:self.nSkySources]["sky_source"] = True
        # The last 10 sources will all contained trail length measurements,
        # increasing in size by 1.5 arcseconds. Only the last three will have
        # lengths which are too long and will be filtered out.
        self.nFilteredTrailedSources = 0
        for srcIdx in range(5, 15):
            self.diaSourceCat[srcIdx]["ext_trailedSources_Naive_length"] = 1.5*(srcIdx-4)
            if 1.5*(srcIdx-4) > 36000/3600.0/24.0 * 30.0:
                self.nFilteredTrailedSources += 1
        # Setting a combination of flags for filtering in tests
        self.diaSourceCat[5]["ext_trailedSources_Naive_flag_off_image"] = True
        self.diaSourceCat[6]["ext_trailedSources_Naive_flag_suspect_long_trail"] = True
        self.diaSourceCat[6]["ext_trailedSources_Naive_flag_edge"] = True
        # As only two of these flags are set, the total number of filtered
        # sources will be self.nFilteredTrailedSources + 2
        self.nFilteredTrailedSources += 2
        self.config = FilterDiaSourceCatalogConfig()
        mjd = 57071.0
        self.utc_jd = mjd + 2_400_000.5 - 35.0 / (24.0 * 60.0 * 60.0)

        self.visitInfo = afwImage.VisitInfo(
            # This incomplete visitInfo is sufficient for testing because the
            # Python constructor sets all other required values to some
            # default.
            exposureTime=30.0,
            darkTime=3.0,
            date=dafBase.DateTime(mjd, system=dafBase.DateTime.MJD),
            boresightRaDec=geom.SpherePoint(0.0, 0.0, geom.degrees),
        )

    def test_run_without_filter(self):
        """Test that when all filters are turned off all sources in the catalog
        are returned.
        """
        self.config.doRemoveSkySources = False
        self.config.doWriteRejectedSkySources = False
        self.config.doTrailedSourceFilter = False
        filterDiaSourceCatalogTask = FilterDiaSourceCatalogTask(config=self.config)
        result = filterDiaSourceCatalogTask.run(self.diaSourceCat, self.visitInfo)
        self.assertEqual(len(result.filteredDiaSourceCat), len(self.diaSourceCat))
        self.assertEqual(len(result.rejectedDiaSources), 0)
        self.assertEqual(len(self.diaSourceCat), self.nSources)

    def test_run_with_filter_sky_only(self):
        """Test that when only the sky filter is turned on the first five
        sources which are flagged as sky objects are filtered out of the
        catalog and the rest are returned.
        """
        self.config.doRemoveSkySources = True
        self.config.doWriteRejectedSkySources = True
        self.config.doTrailedSourceFilter = False
        filterDiaSourceCatalogTask = FilterDiaSourceCatalogTask(config=self.config)
        result = filterDiaSourceCatalogTask.run(self.diaSourceCat, self.visitInfo)
        nExpectedFilteredSources = self.nSources - self.nSkySources
        self.assertEqual(len(result.filteredDiaSourceCat),
                         len(self.diaSourceCat[~self.diaSourceCat['sky_source']]))
        self.assertEqual(len(result.filteredDiaSourceCat), nExpectedFilteredSources)
        self.assertEqual(len(result.rejectedDiaSources), self.nSkySources)
        self.assertEqual(len(self.diaSourceCat), self.nSources)

    def test_run_with_filter_trailed_sources_only(self):
        """Test that when only the trail filter is turned on the correct number
        of sources are filtered out. The filtered sources should be the last
        three sources which have long trails, one source where both the suspect
        trail and edge trail flag are set, and one source where off_image is
        set. All sky objects should remain in the catalog.
        """
        self.config.doRemoveSkySources = False
        self.config.doWriteRejectedSkySources = False
        self.config.doTrailedSourceFilter = True
        filterDiaSourceCatalogTask = FilterDiaSourceCatalogTask(config=self.config)
        result = filterDiaSourceCatalogTask.run(self.diaSourceCat, self.visitInfo)
        nExpectedFilteredSources = self.nSources - self.nFilteredTrailedSources
        self.assertEqual(len(result.filteredDiaSourceCat), nExpectedFilteredSources)
        self.assertEqual(len(self.diaSourceCat), self.nSources)

    def test_run_with_all_filters(self):
        """Test that all sources are filtered out correctly. Only six sources
        should remain in the catalog after filtering.
        """
        self.config.doRemoveSkySources = True
        self.config.doWriteRejectedSkySources = True
        self.config.doTrailedSourceFilter = True
        filterDiaSourceCatalogTask = FilterDiaSourceCatalogTask(config=self.config)
        result = filterDiaSourceCatalogTask.run(self.diaSourceCat, self.visitInfo)
        nExpectedFilteredSources = self.nSources - self.nSkySources - self.nFilteredTrailedSources
        # 5 filtered out sky sources
        # 4 filtered out trailed sources, 2 with long trails 2 with flags
        # 6 sources left
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
