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

import numpy as np
import os
import unittest

from lsst.ap.association import (
    MapApDataConfig,
    MapApDataTask,
    MapDiaSourceConfig,
    MapDiaSourceTask,
    UnpackApdbFlags)
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.image.utils as afwImageUtils
import lsst.geom as geom
from lsst.utils import getPackageDir
import lsst.utils.tests


def make_input_source_catalog(n_objects, add_flags=False):
    """Create tests objects to map into apData products.

    Parameters
    ----------
    n_objects: `int`
        Number of objects to create.
    """
    schema = afwTable.SourceTable.makeMinimalSchema()
    schema.addField("base_NaiveCentroid_x", type="D")
    schema.addField("base_NaiveCentroid_y", type="D")
    schema.addField("base_PsfFlux_instFlux", type="D")
    schema.addField("base_PsfFlux_instFluxErr", type="D")
    schema.addField("ip_diffim_DipoleFit_separation", type="D")
    schema.addField("ip_diffim_DipoleFit_orientation", type="D")
    schema.addField("ip_diffim_DipoleFit_neg_instFlux", type="D")
    schema.addField("ip_diffim_DipoleFit_neg_instFluxErr", type="D")
    schema.addField("ip_diffim_DipoleFit_pos_instFlux", type="D")
    schema.addField("ip_diffim_DipoleFit_pos_instFluxErr", type="D")
    schema.addField("ip_diffim_forced_PsfFlux_instFlux", type="D")
    schema.addField("ip_diffim_forced_PsfFlux_instFluxErr", type="D")
    if add_flags:
        schema.addField("base_PixelFlags_flag", type="Flag")
        schema.addField("base_PixelFlags_flag_offimage", type="Flag")

    objects = afwTable.SourceCatalog(schema)
    objects.preallocate(n_objects)
    objects.definePsfFlux("base_PsfFlux")
    objects.defineCentroid("base_NaiveCentroid")

    for obj_idx in range(n_objects):
        obj = objects.addNew()
        for subSchema in schema:
            if isinstance(obj.get(subSchema.getKey()), geom.Angle):
                obj.set(subSchema.getKey(), 1. * geom.degrees)
            elif subSchema.getField().getName() == "ip_diffim_DipoleFit_neg_instFlux":
                obj.set(subSchema.getKey(), -1)
            else:
                obj.set(subSchema.getKey(), 1)
    return objects


class TestAPDataMapperTask(unittest.TestCase):

    def setUp(self):
        # CFHT Filters from the camera mapper.
        afwImageUtils.resetFilters()
        afwImageUtils.defineFilter('u', lambdaEff=374, alias="u.MP9301")
        afwImageUtils.defineFilter('g', lambdaEff=487, alias="g.MP9401")
        afwImageUtils.defineFilter('r', lambdaEff=628, alias="r.MP9601")
        afwImageUtils.defineFilter('i', lambdaEff=778, alias="i.MP9701")
        afwImageUtils.defineFilter('z', lambdaEff=1170, alias="z.MP9801")

        self.metadata = dafBase.PropertySet()

        self.metadata.set("SIMPLE", "T")
        self.metadata.set("BITPIX", -32)
        self.metadata.set("NAXIS", 2)
        self.metadata.set("NAXIS1", 1024)
        self.metadata.set("NAXIS2", 1153)
        self.metadata.set("RADECSYS", 'FK5')
        self.metadata.set("EQUINOX", 2000.)

        self.metadata.setDouble("CRVAL1", 215.604025685476)
        self.metadata.setDouble("CRVAL2", 53.1595451514076)
        self.metadata.setDouble("CRPIX1", 1109.99981456774)
        self.metadata.setDouble("CRPIX2", 560.018167811613)
        self.metadata.set("CTYPE1", 'RA---SIN')
        self.metadata.set("CTYPE2", 'DEC--SIN')

        self.metadata.setDouble("CD1_1", 5.10808596133527E-05)
        self.metadata.setDouble("CD1_2", 1.85579539217196E-07)
        self.metadata.setDouble("CD2_2", -5.10281493481982E-05)
        self.metadata.setDouble("CD2_1", -8.27440751733828E-07)

        self.wcs = afwGeom.makeSkyWcs(self.metadata)
        self.exposure = afwImage.makeExposure(
            afwImage.makeMaskedImageFromArrays(np.ones((1024, 1153))),
            self.wcs)
        detector = DetectorWrapper(id=23, bbox=self.exposure.getBBox()).detector
        visit = afwImage.VisitInfo(
            exposureId=4321,
            exposureTime=200.,
            date=dafBase.DateTime(nsecs=1400000000 * 10**9))
        self.exposure.setDetector(detector)
        self.exposure.getInfo().setVisitInfo(visit)
        self.exposure.setFilter(afwImage.Filter('g.MP9401'))
        scale = 2
        scaleErr = 1
        self.photoCalib = afwImage.PhotoCalib(scale, scaleErr)
        self.exposure.setPhotoCalib(self.photoCalib)

        self.inputCatalogNoFlags = make_input_source_catalog(10, False)
        self.inputCatalog = make_input_source_catalog(10, True)

    def test_run(self):
        """Test the generic data product mapper.
        """
        outSchema = afwTable.SourceTable.makeMinimalSchema()
        outSchema.addField("psFlux", type="D")
        outSchema.addField("psFluxErr", type="D")

        mapApDConfig = MapApDataConfig()
        mapApDConfig.copyColumns = {
            "id": "id",
            "parent": "parent",
            "coord_ra": "coord_ra",
            "coord_dec": "coord_dec",
            "slot_PsfFlux_instFlux": "psFlux",
            "slot_PsfFlux_instFluxErr": "psFluxErr"
        }

        mapApD = MapApDataTask(inputSchema=self.inputCatalog.schema,
                               outputSchema=outSchema,
                               config=mapApDConfig)
        outputCatalog = mapApD.run(self.inputCatalog)

        for inObj, outObj in zip(self.inputCatalog, outputCatalog):
            for inputName, outputName in mapApDConfig.copyColumns.items():
                self.assertEqual(inObj[inputName], outObj[outputName])

    def test_run_dia_source(self):
        """Test the DiaSource specific data product mapper/calibrator.
        """
        mapApDConfig = self._create_map_dia_source_config()
        mapApD = MapDiaSourceTask(inputSchema=self.inputCatalog.schema,
                                  config=mapApDConfig)
        outputCatalog = mapApD.run(self.inputCatalog, self.exposure)

        expectedMeanDip = 2.
        expectedDiffFlux = 0.
        expectedLength = 0.18379083

        for inObj, outObj in zip(self.inputCatalog, outputCatalog):
            self.assertEqual(
                outObj["ccdVisitId"],
                self.exposure.getInfo().getVisitInfo().getExposureId())
            self.assertEqual(
                outObj["midPointTai"],
                self.exposure.getInfo().getVisitInfo().getDate().get(
                    system=dafBase.DateTime.MJD))
            self.assertEqual(
                outObj["flags"],
                1 * 2 ** 0 + 1 * 2 ** 1)
            for inputName, outputName in mapApDConfig.copyColumns.items():
                if inputName.startswith("slot_PsfFlux"):
                    self._test_calibrated_flux(inObj, outObj)
                else:
                    self.assertEqual(inObj[inputName], outObj[outputName])
                self.assertEqual(outObj["dipMeanFlux"], expectedMeanDip)
                self.assertEqual(outObj["dipFluxDiff"], expectedDiffFlux)
                self.assertAlmostEqual(outObj["dipLength"], expectedLength)
        # Mapper should always emit standardized filters
        self.assertEqual(outObj["filterName"], 'g')

    def _create_map_dia_source_config(self):
        """Create a test config for use in MapDiaSourceTask.

        Returns
        -------
        configurable : `lsst.pex.config.Config`
            Configurable for use in MapDiaSourceTask.
        """
        configurable = MapDiaSourceConfig()
        configurable.copyColumns = {
            "id": "id",
            "parent": "parent",
            "coord_ra": "coord_ra",
            "coord_dec": "coord_dec",
            "slot_PsfFlux_instFlux": "psFlux",
            "slot_PsfFlux_instFluxErr": "psFluxErr"
        }
        configurable.calibrateColumns = ["slot_PsfFlux"]
        configurable.flagMap = os.path.join(
            getPackageDir("ap_association"),
            "tests",
            "test-flag-map.yaml")

        return configurable

    def test_run_dia_source_wrong_flags(self):
        """Test that the proper errors are thrown when requesting flag columns
        that are not in the input schema.
        """
        mapApDConfig = self._create_map_dia_source_config()
        with self.assertRaises(KeyError):
            MapDiaSourceTask(inputSchema=self.inputCatalogNoFlags.schema,
                             config=mapApDConfig)

    def test_calibrateFluxes(self):
        """Test that flux calibration works as expected.
        """
        outSchema = afwTable.SourceTable.makeMinimalSchema()
        outSchema.addField("psFlux", type="D")
        outSchema.addField("psFluxErr", type="D")

        outputCatalog = afwTable.SourceCatalog(outSchema)
        outRecord = outputCatalog.addNew()

        mapApDConfig = self._create_map_dia_source_config()
        mapApD = MapDiaSourceTask(inputSchema=self.inputCatalog.schema,
                                  config=mapApDConfig)

        mapApD.calibrateFluxes(self.inputCatalog[0],
                               outRecord,
                               self.photoCalib)
        self._test_calibrated_flux(self.inputCatalog[0], outRecord)

    def _test_calibrated_flux(self, inputRecord, outputRecord):
        """Compare calibrated fluxes to expectation from zero point.

        Parameters
        ----------
        inputRecord: `lsst.afw.table.SourceRecord`
            Input source record with uncalibrated flux values.
        outputRecord: `lsst.afw.table.SourceRecord`
            Source record with calibrated fluxes.
        """
        expected = self.photoCalib.instFluxToNanojansky(inputRecord["slot_PsfFlux_instFlux"],
                                                        inputRecord["slot_PsfFlux_instFluxErr"])
        self.assertAlmostEqual(outputRecord["psFlux"], expected.value)
        self.assertAlmostEqual(outputRecord["psFluxErr"], expected.error)

    def test_bit_unpacker(self):
        """Test that the integer bit packer is functioning correctly.
        """
        mapApDConfig = self._create_map_dia_source_config()
        mapApD = MapDiaSourceTask(inputSchema=self.inputCatalog.schema,
                                  config=mapApDConfig)
        for idx, obj in enumerate(self.inputCatalog):
            if idx in [1, 3, 5]:
                obj.set("base_PixelFlags_flag", 0)
            if idx in [1, 4, 6]:
                obj.set("base_PixelFlags_flag_offimage", 0)
        outputCatalog = mapApD.run(self.inputCatalog, self.exposure)

        unpacker = UnpackApdbFlags(mapApDConfig.flagMap, "DiaSource")
        flag_values = unpacker.unpack(outputCatalog.get("flags"), "flags")

        for idx, flag in enumerate(flag_values):
            if idx in [1, 3, 5]:
                self.assertFalse(flag['base_PixelFlags_flag'])
            else:
                self.assertTrue(flag['base_PixelFlags_flag'])

            if idx in [1, 4, 6]:
                self.assertFalse(flag['base_PixelFlags_flag_offimage'])
            else:
                self.assertTrue(flag['base_PixelFlags_flag_offimage'])


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
