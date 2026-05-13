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
import numpy as np

from lsst.ap.association.filterDiaSourceCatalog import (FilterDiaSourceCatalogConfig,
                                                        FilterDiaSourceCatalogTask,
                                                        FilterDiaSourceReliabilityConfig,
                                                        FilterDiaSourceReliabilityTask)
import lsst.geom as geom
import lsst.afw.detection as afwDetect
import lsst.afw.geom as afwGeom
import lsst.meas.base.tests as measTests
import lsst.utils.tests
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.daf.base as dafBase


class TestFilterDiaSourceCatalogTask(unittest.TestCase):

    def setUp(self):
        self.config = FilterDiaSourceCatalogConfig()

        self.nSkySources = 5
        self.nCrCenterSources = 6
        self.nFakeFlagSources = 4
        self.nNegativeSources = 7
        self.nTrailedSources = 10
        self.pixelScale = 0.2  # arcseconds/pixel
        self.nSources = (self.nSkySources + self.nCrCenterSources + self.nFakeFlagSources
                         + self.nNegativeSources + self.nTrailedSources)
        self.yLoc = 100
        self.expId = 4321
        self.bbox = geom.Box2I(geom.Point2I(0, 0),
                               geom.Extent2I(1024, 1153))
        dataset = measTests.TestDataset(self.bbox)
        for srcIdx in range(self.nSources):
            dataset.addSource(10000.0, geom.Point2D(srcIdx, self.yLoc))
        schema = dataset.makeMinimalSchema()
        schema.addField("sky_source", type="Flag", doc="Sky objects.")
        schema.addField("slot_Centroid_flag", type="Flag", doc="The centroid calculation "
                        "failed. The source position is incorrect.")
        schema.addField("base_PixelFlags_flag_crCenter", type="Flag", doc="A cosmic ray was detected "
                        "and interpolated in this object's center.")
        schema.addField("base_PixelFlags_flag_high_varianceCenterAll", type="Flag", doc="The object was "
                        "detected in a region with exceptionally high template variance.")
        schema.addField("fakeBadFlag", type="Flag", doc="A fake flag to test a badFlagList longer "
                        "than one item.")
        schema.addField("ip_diffim_forced_PsfFlux_instFlux", type="F",
                        doc="Forced photometry flux for a point source model measured on the visit image "
                        "centered at DiaSource position.")
        schema.addField("ip_diffim_forced_PsfFlux_instFluxErr", type="F",
                        doc="Estimated uncertainty of ip_diffim_forced_PsfFlux_instFlux.")
        schema.addField('ext_trailedSources_Naive_flag', type="Flag",
                        doc="General trail measurement failure flag")
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
        schema.addField("reliability", type="F",
                        doc="Reliability of the source")
        schema.addField("reliabilityVersion", type=str, size=7,
                        doc="Version of the reliability model")
        _, self.diaSourceCat = dataset.realize(10.0, schema, randomSeed=1234)

        # set the sky_source flag for the first set
        self.diaSourceCat[0:self.nSkySources]["sky_source"] = True

        # set the pixelFlags_crCenter flag
        crCenter_offset = self.nSkySources
        self.diaSourceCat[crCenter_offset:crCenter_offset+self.nCrCenterSources][
            "base_PixelFlags_flag_crCenter"
        ] = True

        # set the fakeBadFlag flag
        fakeFlag_offset = crCenter_offset + self.nCrCenterSources
        self.diaSourceCat[fakeFlag_offset:fakeFlag_offset+self.nFakeFlagSources][
            "fakeBadFlag"
        ] = True

        # create increasingly negative ip_diffim_forced_PsfFlux_instFlux/ip_diffim_forced_PsfFlux_instFluxErr
        self.nRemovedNegativeSources = 0
        negativeSources_offset = fakeFlag_offset + self.nFakeFlagSources
        for i, srcIdx in enumerate(range(negativeSources_offset,
                                         negativeSources_offset+self.nNegativeSources)):
            self.diaSourceCat[srcIdx]["ip_diffim_forced_PsfFlux_instFlux"] = -0.5 * i
            self.diaSourceCat[srcIdx]["ip_diffim_forced_PsfFlux_instFluxErr"] = 1.01
            if (-0.5 * i)/1.01 < self.config.minAllowedDirectSnr:
                self.nRemovedNegativeSources += 1

        # The last 10 sources will all contained trail length measurements,
        # increasing in size by 1.5 arcseconds. Only the last three will have
        # lengths which are too long and will be filtered out.
        self.nFilteredTrailedSources = 0
        trail_offset = negativeSources_offset + self.nNegativeSources
        for i, srcIdx in enumerate(range(trail_offset, trail_offset+self.nTrailedSources)):
            self.diaSourceCat[srcIdx]["ext_trailedSources_Naive_length"] = 1.5*(i+1)/self.pixelScale
            if 1.5*(i+1) > 36000/3600.0/24.0 * 30.0:
                self.nFilteredTrailedSources += 1
        # Setting a combination of flags for filtering in tests
        self.diaSourceCat[trail_offset+1]["ext_trailedSources_Naive_flag_off_image"] = True
        self.diaSourceCat[trail_offset+2]["ext_trailedSources_Naive_flag_suspect_long_trail"] = True
        self.diaSourceCat[trail_offset+2]["ext_trailedSources_Naive_flag_edge"] = True
        # As only two of these flags are set, the total number of filtered
        # sources will be self.nFilteredTrailedSources + 2
        self.nFilteredTrailedSources += 2

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
        self.config.badFlagList = []
        self.config.doRemoveNegativeDirectImageSources = False
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
        self.config.badFlagList = []
        self.config.doRemoveNegativeDirectImageSources = False
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

    def test_run_with_filter_defaultBadFlagList_only(self):
        """Test that when only the CR center filter is turned on the six sources which are flagged
        as base_PixelFlags_flag_crCenter are filtered out of the catalog and the rest are returned.
        """
        self.config.doRemoveSkySources = False
        self.config.doRemoveNegativeDirectImageSources = False
        self.config.doWriteRejectedSkySources = False
        self.config.doTrailedSourceFilter = False
        filterDiaSourceCatalogTask = FilterDiaSourceCatalogTask(config=self.config)
        result = filterDiaSourceCatalogTask.run(self.diaSourceCat, self.visitInfo)
        nExpectedFilteredSources = self.nSources - self.nCrCenterSources
        self.assertEqual(len(result.filteredDiaSourceCat),
                         len(self.diaSourceCat[~self.diaSourceCat['base_PixelFlags_flag_crCenter']]))
        self.assertEqual(len(result.filteredDiaSourceCat), nExpectedFilteredSources)
        self.assertEqual(len(result.rejectedDiaSources), self.nCrCenterSources)
        self.assertEqual(len(self.diaSourceCat), self.nSources)

    def test_run_with_filter_nonDefaultBadFlagList_only(self):
        """Test that the badFlagList filters appropriately when it is not when the default configuration.
        The six sources flagged base_PixelFlags_flag_crCenter and the four sources flagged fakeBadFlag
        should be filtered out of the catalog and the rest are returned.
        """
        self.config.doRemoveSkySources = False
        self.config.badFlagList = ["base_PixelFlags_flag_crCenter", "fakeBadFlag"]
        self.config.doRemoveNegativeDirectImageSources = False
        self.config.doWriteRejectedSkySources = False
        self.config.doTrailedSourceFilter = False
        filterDiaSourceCatalogTask = FilterDiaSourceCatalogTask(config=self.config)
        result = filterDiaSourceCatalogTask.run(self.diaSourceCat, self.visitInfo)
        nExpectedFilteredSources = self.nSources - self.nCrCenterSources - self.nFakeFlagSources
        self.assertEqual(len(result.filteredDiaSourceCat),
                         len(self.diaSourceCat[~self.diaSourceCat['base_PixelFlags_flag_crCenter']
                             & ~self.diaSourceCat['fakeBadFlag']]))
        self.assertEqual(len(result.filteredDiaSourceCat), nExpectedFilteredSources)
        self.assertEqual(len(result.rejectedDiaSources), self.nCrCenterSources + self.nFakeFlagSources)
        self.assertEqual(len(self.diaSourceCat), self.nSources)

    def test_run_with_filter_negative_only(self):
        """Test that when only the negative filter is turned on then
        sources which below the negtive snr cut are filtered out of the
        catalog and the rest are returned.
        """
        self.config.doRemoveSkySources = False
        self.config.badFlagList = []
        self.config.doRemoveNegativeDirectImageSources = True
        self.config.doWriteRejectedSkySources = True
        self.config.doTrailedSourceFilter = False
        filterDiaSourceCatalogTask = FilterDiaSourceCatalogTask(config=self.config)
        result = filterDiaSourceCatalogTask.run(self.diaSourceCat, self.visitInfo)
        nExpectedFilteredSources = self.nSources - self.nRemovedNegativeSources
        self.assertEqual(len(result.filteredDiaSourceCat), nExpectedFilteredSources)
        self.assertEqual(len(result.rejectedDiaSources), self.nRemovedNegativeSources)
        self.assertEqual(len(self.diaSourceCat), self.nSources)
        self.assertEqual(np.sum(result.filteredDiaSourceCat['ip_diffim_forced_PsfFlux_instFlux']
                                / result.filteredDiaSourceCat['ip_diffim_forced_PsfFlux_instFluxErr']
                                < self.config.minAllowedDirectSnr), 0)
        self.assertEqual(np.sum(result.rejectedDiaSources['ip_diffim_forced_PsfFlux_instFlux']
                                / result.rejectedDiaSources['ip_diffim_forced_PsfFlux_instFluxErr']
                                < self.config.minAllowedDirectSnr),
                         self.nRemovedNegativeSources)

    def test_run_with_filter_negative_and_sky(self):
        """Test concatenating rejects when both sky and negative filtering
        are on.
        """
        self.config.doRemoveSkySources = True
        self.config.badFlagList = []
        self.config.doRemoveNegativeDirectImageSources = True
        self.config.doWriteRejectedSkySources = True
        self.config.doTrailedSourceFilter = False
        filterDiaSourceCatalogTask = FilterDiaSourceCatalogTask(config=self.config)
        result = filterDiaSourceCatalogTask.run(self.diaSourceCat, self.visitInfo)
        nExpectedFilteredSources = self.nSources - self.nSkySources - self.nRemovedNegativeSources
        self.assertEqual(len(result.filteredDiaSourceCat), nExpectedFilteredSources)
        self.assertEqual(len(result.rejectedDiaSources), self.nSkySources + self.nRemovedNegativeSources)
        self.assertEqual(len(self.diaSourceCat), self.nSources)

    def test_run_with_filter_reliability_only(self):
        """Test that when only the reliability filter is turned on,
        sources below the reliability threshold are filtered out."""

        reliability_threshold = 0.7
        config = FilterDiaSourceReliabilityConfig()
        config.minReliability = reliability_threshold

        schema = afwTable.SourceTable.makeMinimalSchema()
        schema.addField("score", type="F",
                        doc="Reliability of the source")
        schema.addField("version", type=str, size=7,
                        doc="Version of the reliability model.")
        reliabilityCat = afwTable.SourceCatalog(schema)
        reliabilityCat.reserve(self.nSources)

        # Set reliability: first half below threshold, second half above
        for srcIdx in range(self.nSources):
            reliabilityCat.addNew()
            if srcIdx < self.nSources // 2:
                reliabilityCat[srcIdx]["score"] = 0.25
            else:
                reliabilityCat[srcIdx]["score"] = 0.95
        reliabilityCat['id'] = self.diaSourceCat['id']
        nLowReliability = np.sum(reliabilityCat["score"] < 0.5)
        filterDiaSourceCatalogTask = FilterDiaSourceReliabilityTask(config=config)
        result = filterDiaSourceCatalogTask.run(self.diaSourceCat, reliabilityCat)
        self.assertEqual(len(result.filteredDiaSources), self.nSources - nLowReliability)
        self.assertEqual(len(result.rejectedDiaSources), nLowReliability)
        self.assertTrue(np.all(result.filteredDiaSources["reliability"] >= reliability_threshold))
        self.assertTrue(np.all(result.rejectedDiaSources["reliability"] < reliability_threshold))

    def test_run_with_filter_trailed_sources_only(self):
        """Test that when only the trail filter is turned on the correct number
        of sources are filtered out. The filtered sources should be the last
        three sources which have long trails, one source where both the suspect
        trail and edge trail flag are set, and one source where off_image is
        set. All sky objects should remain in the catalog.
        """
        self.config.doRemoveSkySources = False
        self.config.badFlagList = []
        self.config.doRemoveNegativeDirectImageSources = False
        self.config.doWriteRejectedSkySources = False
        self.config.doTrailedSourceFilter = True
        filterDiaSourceCatalogTask = FilterDiaSourceCatalogTask(config=self.config)
        result = filterDiaSourceCatalogTask.run(self.diaSourceCat, self.visitInfo)
        nExpectedFilteredSources = self.nSources - self.nFilteredTrailedSources
        self.assertEqual(len(result.filteredDiaSourceCat), nExpectedFilteredSources)
        self.assertEqual(len(self.diaSourceCat), self.nSources)

    def test_run_with_all_filters(self):
        """Test that all sources are filtered out correctly. Only 15 sources
        should remain in the catalog after filtering.
        """
        self.config.doRemoveSkySources = True
        self.config.doRemoveNegativeDirectImageSources = True
        self.config.doWriteRejectedSkySources = True
        self.config.doTrailedSourceFilter = True
        filterDiaSourceCatalogTask = FilterDiaSourceCatalogTask(config=self.config)
        result = filterDiaSourceCatalogTask.run(self.diaSourceCat, self.visitInfo)
        nExpectedFilteredSources = (self.nSources - self.nSkySources
                                    - self.nCrCenterSources
                                    - self.nFilteredTrailedSources
                                    - self.nRemovedNegativeSources)
        nExpectedRejectedSources = (self.nSkySources
                                    + self.nCrCenterSources
                                    + self.nRemovedNegativeSources)
        # 32 total sources
        # 5 filtered out sky sources
        # 6 filtered out sources with cosmic ray detections
        # 2 filtered out negative sources
        # 4 filtered out trailed sources, 2 with long trails 2 with flags
        # 15 sources left
        self.assertEqual(len(result.filteredDiaSourceCat), nExpectedFilteredSources)
        # 17 sources rejected, 4 trailed sources not included, 13 rejected sources in catalog.
        self.assertEqual(len(result.rejectedDiaSources), nExpectedRejectedSources)
        self.assertEqual(len(self.diaSourceCat), self.nSources)

    def test_pixelScale_calculation(self):
        """Check the calculation of the pixel scale from the input catalog.
        """
        self.config.doTrailedSourceFilter = True
        filterDiaSourceCatalogTask = FilterDiaSourceCatalogTask(config=self.config)
        scale = filterDiaSourceCatalogTask._estimate_pixel_scale(self.diaSourceCat)
        # Should be almost but not actually equal
        self.assertNotEqual(self.config.estimatedPixelScale, scale)
        self.assertAlmostEqual(self.config.estimatedPixelScale, scale, places=6)

        # If the estimatedPixelScale is very different, that value should be
        #  used exactly and it should not raise an error.
        self.config.estimatedPixelScale = 1.2
        filterDiaSourceCatalogTask = FilterDiaSourceCatalogTask(config=self.config)
        scale = filterDiaSourceCatalogTask._estimate_pixel_scale(self.diaSourceCat)
        self.assertEqual(self.config.estimatedPixelScale, scale)

    def test_run_with_filter_centroid_flag(self):
        """Test that sources with slot_Centroid_flag set are filtered with the
        default badFlagList configuration.
        """
        self.config.doRemoveSkySources = False
        self.config.doRemoveNegativeDirectImageSources = False
        self.config.doWriteRejectedSkySources = False
        self.config.doTrailedSourceFilter = False
        # Default badFlagList includes slot_Centroid_flag
        self.assertIn("slot_Centroid_flag", self.config.badFlagList)

        # Set slot_Centroid_flag on a few sources that aren't already flagged
        nCentroidFlagged = 3
        # Pick sources at the end that have no other badFlagList flags
        trail_offset = (self.nSkySources + self.nCrCenterSources
                        + self.nFakeFlagSources + self.nNegativeSources)
        for i in range(nCentroidFlagged):
            self.diaSourceCat[trail_offset + i]["slot_Centroid_flag"] = True

        filterDiaSourceCatalogTask = FilterDiaSourceCatalogTask(config=self.config)
        result = filterDiaSourceCatalogTask.run(self.diaSourceCat, self.visitInfo)
        # Default badFlagList filters crCenter + slot_Centroid_flag + high_varianceCenterAll.
        # crCenter flags were already set for nCrCenterSources.
        # Our added centroid flags are on previously unflagged sources.
        nExpectedFiltered = self.nSources - self.nCrCenterSources - nCentroidFlagged
        self.assertEqual(len(result.filteredDiaSourceCat), nExpectedFiltered)

    def test_check_dia_source_trail_bbox_no_flag(self):
        """Test that sources without ext_trailedSources_Naive_flag are not
        flagged by the bbox check, regardless of bbox size.
        """
        self.config.doTrailedSourceFilter = True
        task = FilterDiaSourceCatalogTask(config=self.config)
        exposure_time = self.visitInfo.getExposureTime()
        # None of the sources have the trail flag set by default
        bbox_mask = task._check_dia_source_trail_bbox(self.diaSourceCat, exposure_time)
        self.assertEqual(np.sum(bbox_mask), 0)

    def test_check_dia_source_trail_bbox_no_flag_large_bbox(self):
        """Test that sources without ext_trailedSources_Naive_flag are not
        flagged by the bbox check even when the bounding box is large enough
        that it would otherwise trigger the bbox filter.
        """
        self.config.doTrailedSourceFilter = True
        task = FilterDiaSourceCatalogTask(config=self.config)
        exposure_time = self.visitInfo.getExposureTime()
        pixelScale = task._estimate_pixel_scale(self.diaSourceCat)
        max_length_pixels = self.config.max_trail_length * exposure_time / pixelScale

        # Give source 0 a large footprint exceeding the threshold, but do NOT
        # set ext_trailedSources_Naive_flag — the bbox check should be skipped.
        srcIdx = 0
        largeSpan = afwGeom.SpanSet(
            geom.Box2I(geom.Point2I(0, 0),
                       geom.Extent2I(int(max_length_pixels) + 10, 5))
        )
        footprint = afwDetect.Footprint(largeSpan)
        self.diaSourceCat[srcIdx].setFootprint(footprint)
        self.assertFalse(self.diaSourceCat[srcIdx]["ext_trailedSources_Naive_flag"])

        bbox_mask = task._check_dia_source_trail_bbox(self.diaSourceCat, exposure_time)
        self.assertFalse(bbox_mask[srcIdx])
        self.assertEqual(np.sum(bbox_mask), 0)

    def test_check_dia_source_trail_bbox_flag_small_bbox(self):
        """Test that sources with ext_trailedSources_Naive_flag set but a small bounding
        box are not flagged.
        """
        self.config.doTrailedSourceFilter = True
        task = FilterDiaSourceCatalogTask(config=self.config)
        exposure_time = self.visitInfo.getExposureTime()
        # Set the trail flag on a source with a small footprint (PSF-sized)
        self.diaSourceCat[0]["ext_trailedSources_Naive_flag"] = True
        bbox_mask = task._check_dia_source_trail_bbox(self.diaSourceCat, exposure_time)
        # The default PSF footprint is small (~7 pixels), well below the
        # threshold for max_trail_length * exposure_time / pixelScale
        self.assertFalse(bbox_mask[0])

    def test_check_dia_source_trail_bbox_flag_large_bbox(self):
        """Test that sources with ext_trailedSources_Naive_flag set and a
        large bounding box are flagged.
        """
        self.config.doTrailedSourceFilter = True
        task = FilterDiaSourceCatalogTask(config=self.config)
        exposure_time = self.visitInfo.getExposureTime()
        pixelScale = task._estimate_pixel_scale(self.diaSourceCat)
        max_length_pixels = self.config.max_trail_length * exposure_time / pixelScale

        # Set the trail flag on a source and give it a large footprint
        srcIdx = 0
        self.diaSourceCat[srcIdx]["ext_trailedSources_Naive_flag"] = True
        largeSpan = afwGeom.SpanSet(
            geom.Box2I(geom.Point2I(0, 0),
                       geom.Extent2I(int(max_length_pixels) + 10, 5))
        )
        footprint = afwDetect.Footprint(largeSpan)
        self.diaSourceCat[srcIdx].setFootprint(footprint)

        bbox_mask = task._check_dia_source_trail_bbox(self.diaSourceCat, exposure_time)
        self.assertTrue(bbox_mask[srcIdx])

    def test_check_dia_source_trail_bbox_diagonal(self):
        """Test a source with a large diagonal but small width/height should still
        be flagged.
        """
        self.config.doTrailedSourceFilter = True
        task = FilterDiaSourceCatalogTask(config=self.config)
        exposure_time = self.visitInfo.getExposureTime()
        pixelScale = task._estimate_pixel_scale(self.diaSourceCat)
        max_length_pixels = self.config.max_trail_length * exposure_time / pixelScale

        srcIdx = 1
        self.diaSourceCat[srcIdx]["ext_trailedSources_Naive_flag"] = True
        # Create a bbox whose individual sides are below max_length_pixels
        # but whose diagonal exceeds it.
        side = int(max_length_pixels * 0.8)
        largeSpan = afwGeom.SpanSet(
            geom.Box2I(geom.Point2I(0, 0), geom.Extent2I(side, side))
        )
        footprint = afwDetect.Footprint(largeSpan)
        self.diaSourceCat[srcIdx].setFootprint(footprint)

        bbox_mask = task._check_dia_source_trail_bbox(self.diaSourceCat, exposure_time)
        # diagonal = side * sqrt(2) ~ max_length_pixels * 1.13 > max_length_pixels
        self.assertTrue(bbox_mask[srcIdx])

    def test_run_with_trail_and_bbox_filter(self):
        """Test doTrailedSourceFilter=True filters sources
        via both trail length and bbox fallback.
        """
        self.config.doRemoveSkySources = False
        self.config.badFlagList = []
        self.config.doRemoveNegativeDirectImageSources = False
        self.config.doWriteRejectedSkySources = False
        self.config.doTrailedSourceFilter = True
        task = FilterDiaSourceCatalogTask(config=self.config)
        exposure_time = self.visitInfo.getExposureTime()
        pixelScale = task._estimate_pixel_scale(self.diaSourceCat)
        max_length_pixels = self.config.max_trail_length * exposure_time / pixelScale

        # Set the trail flag and give a large footprint to a source that
        # wouldn't be caught by the trail length check (first source).
        srcIdx = 0
        self.diaSourceCat[srcIdx]["ext_trailedSources_Naive_flag"] = True
        largeSpan = afwGeom.SpanSet(
            geom.Box2I(geom.Point2I(0, 0),
                       geom.Extent2I(int(max_length_pixels) + 10, 5))
        )
        footprint = afwDetect.Footprint(largeSpan)
        self.diaSourceCat[srcIdx].setFootprint(footprint)

        result = task.run(self.diaSourceCat, self.visitInfo)
        # The bbox-flagged source should have been removed
        nExpectedFiltered = self.nSources - self.nFilteredTrailedSources - 1
        self.assertEqual(len(result.filteredDiaSourceCat), nExpectedFiltered)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
