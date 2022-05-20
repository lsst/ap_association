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
import pandas as pd
import unittest

from lsst.ap.association.ssoAssociation import SolarSystemAssociationTask
import astshim as ast
import lsst.geom as geom
import lsst.afw.geom as afwGeom
import lsst.meas.base.tests as measTests
import lsst.utils.tests


class TestSolarSystemAssociation(unittest.TestCase):

    def setUp(self):
        # Make fake sources
        self.nSources = 10
        self.bbox = geom.Box2I(geom.Point2I(0, 0),
                               geom.Extent2I(1024, 1153))
        self.xyLoc = 100
        dataset = measTests.TestDataset(self.bbox)
        for srcIdx in range(self.nSources):
            dataset.addSource(100000.0,
                              geom.Point2D(srcIdx*self.xyLoc,
                                           srcIdx*self.xyLoc))
        schema = dataset.makeMinimalSchema()
        schema.addField("base_PixelFlags_flag", type="Flag")
        schema.addField("base_PixelFlags_flag_offimage", type="Flag")
        self.exposure, catalog = dataset.realize(
            10.0, schema, randomSeed=1234)
        for src in catalog:
            src.setCoord(self.exposure.getWcs().pixelToSky(src.getCentroid()))
        # Non-invertible WCS to test robustness to distortions.
        # Coefficients transform (x - 1e-7*x^3 -> x, y -> y); see docs for PolyMap.
        pixelCoeffs = np.array([[-1.0e-7, 1, 3, 0],
                                [1.0, 1, 1, 0],
                                [1.0, 2, 0, 1],
                                ])
        self.exposure.setWcs(afwGeom.makeModifiedWcs(
            afwGeom.TransformPoint2ToPoint2(ast.PolyMap(pixelCoeffs, 2, options="IterInverse=1")),
            self.exposure.wcs,
            modifyActualPixels=False
        ))

        # Convert to task required format
        self.testDiaSources = catalog.asAstropy().to_pandas()
        self.testDiaSources.rename(columns={"coord_ra": "ra",
                                            "coord_dec": "decl"},
                                   inplace=True,)
        self.testDiaSources.loc[:, "ra"] = np.rad2deg(self.testDiaSources["ra"])
        self.testDiaSources.loc[:, "decl"] = np.rad2deg(self.testDiaSources["decl"])
        self.testDiaSources["ssObjectId"] = 0

        # Grab a subset to treat as solar system objects
        self.testSsObjects = self.testDiaSources[2:8].reset_index()
        # Assign them ids starting from 1.
        self.testSsObjects.loc[:, "ssObjectId"] = np.arange(
            1, len(self.testSsObjects) + 1, dtype=int,)
        self.testSsObjects["Err(arcsec)"] = np.ones(len(self.testSsObjects))

    def test_run(self):
        """Test that association and id assignment work as expected.
        """
        ssAssocTask = SolarSystemAssociationTask()
        results = ssAssocTask.run(self.testDiaSources,
                                  self.testSsObjects,
                                  self.exposure)
        self.assertEqual(self.nSources,
                         len(results.ssoAssocDiaSources)
                         + len(results.unAssocDiaSources))
        self.assertEqual(len(self.testSsObjects),
                         len(results.ssoAssocDiaSources))
        self.assertEqual(self.nSources - len(self.testSsObjects),
                         len(results.unAssocDiaSources))
        for idx, ssObject in self.testSsObjects.iterrows():
            self.assertEqual(
                ssObject["ssObjectId"],
                results.ssoAssocDiaSources.iloc[idx]["ssObjectId"])

    def test_mask(self):
        """Test that masking against the CCD bounding box works as expected.
        """
        ssAssocTask = SolarSystemAssociationTask()
        # Test will all inside ccd
        maskedObjects = ssAssocTask._maskToCcdRegion(self.testSsObjects,
                                                     self.exposure,
                                                     1.0)
        self.assertEqual(len(maskedObjects), len(self.testSsObjects))

        # Add a new SolarSystemObjects outside of the bbox and test that it
        # is excluded.
        testObjects = self.testSsObjects.loc[:2].reset_index(drop=True)
        testObjects.loc[0, "ra"] = 150
        testObjects.loc[0, "decl"] = 80
        testObjects.loc[1, "ra"] = 150
        testObjects.loc[1, "decl"] = -80
        # Coordinates are chosen so that the inverse WCS erroneously maps them
        # to inside the box (to (74.5, 600.6) instead of around (1745, 600.6)).
        testObjects.loc[2, "ra"] = 44.91215199831453
        testObjects.loc[2, "decl"] = 45.001331943391406
        maskedObjects = ssAssocTask._maskToCcdRegion(
            pd.concat([self.testSsObjects, testObjects]),
            self.exposure,
            1.0)
        self.assertEqual(len(maskedObjects), len(self.testSsObjects))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
