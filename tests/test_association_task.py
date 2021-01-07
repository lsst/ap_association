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

from lsst.afw.cameraGeom.testUtils import DetectorWrapper
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.image.utils as afwImageUtils
import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
import lsst.geom as geom
import lsst.sphgeom as sphgeom
import lsst.utils.tests

from lsst.ap.association import \
    AssociationTask, \
    make_dia_source_schema, \
    make_dia_object_schema


def create_test_points(point_locs_deg,
                       wcs=None,
                       start_id=0,
                       schema=None,
                       scatter_arcsec=1.0,
                       indexer_ids=None,
                       associated_ids=None):
    """Create dummy DIASources or DIAObjects for use in our tests.
    Parameters
    ----------
    point_locs_deg : array-like (N, 2) of `float`s
        Positions of the test points to create in RA, DEC.
    wcs : `lsst.afw.geom.SkyWcs`
        Wcs to convert RA/Dec to x/y if provided.
    start_id : `int`
        Unique id of the first object to create. The remaining sources are
        incremented by one from the first id.
    schema : `lsst.afw.table.Schema`
        Schema of the objects to create. Defaults to the DIASource schema.
    scatter_arcsec : `float`
        Scatter to add to the position of each DIASource.
    indexer_ids : `list` of `ints`s
        Id numbers of pixelization indexer to store. Must be the same length
        as the first dimension of point_locs_deg.
    associated_ids : `list` of `ints`s
        Id numbers of associated DIAObjects to store. Must be the same length
        as the first dimension of point_locs_deg.
    Returns
    -------
    test_points : `lsst.afw.table.SourceCatalog`
        Catalog of points to test.
    """
    if schema is None:
        schema = make_dia_source_schema()
    sources = afwTable.SourceCatalog(schema)

    for src_idx, (ra, dec,) in enumerate(point_locs_deg):
        src = sources.addNew()
        src['id'] = src_idx + start_id
        coord = geom.SpherePoint(ra, dec, geom.degrees)
        if scatter_arcsec > 0.0:
            coord = coord.offset(
                np.random.rand() * 360 * geom.degrees,
                np.random.rand() * scatter_arcsec * geom.arcseconds)
        if indexer_ids is not None:
            src['pixelId'] = indexer_ids[src_idx]
        if associated_ids is not None:
            src['diaObjectId'] = associated_ids[src_idx]
        src.setCoord(coord)

        if wcs is not None:
            xyCentroid = wcs.skykToPixel(coord)
            src.set("x", xyCentroid.getX())
            src.set("y", xyCentroid.getY())

    return sources


def create_test_points_pandas(point_locs_deg,
                              wcs=None,
                              start_id=0,
                              schema=None,
                              scatter_arcsec=1.0,
                              indexer_ids=None,
                              associated_ids=None):
    """Create dummy DIASources or DIAObjects for use in our tests.
    Parameters
    ----------
    point_locs_deg : array-like (N, 2) of `float`s
        Positions of the test points to create in RA, DEC.
    wcs : `lsst.afw.geom.SkyWcs`
        Wcs to convert RA/Dec to x/y if provided.
    start_id : `int`
        Unique id of the first object to create. The remaining sources are
        incremented by one from the first id.
    schema : `lsst.afw.table.Schema`
        Schema of the objects to create. Defaults to the DIASource schema.
    scatter_arcsec : `float`
        Scatter to add to the position of each DIASource.
    indexer_ids : `list` of `ints`s
        Id numbers of pixelization indexer to store. Must be the same length
        as the first dimension of point_locs_deg.
    associated_ids : `list` of `ints`s
        Id numbers of associated DIAObjects to store. Must be the same length
        as the first dimension of point_locs_deg.
    Returns
    -------
    test_points : `pandas.DataFrame`
        Catalog of points to test.
    """
    if schema is None:
        schema = make_dia_source_schema()
    sources = afwTable.SourceCatalog(schema)

    for src_idx, (ra, dec,) in enumerate(point_locs_deg):
        src = sources.addNew()
        src['id'] = src_idx + start_id
        coord = geom.SpherePoint(ra, dec, geom.degrees)
        if scatter_arcsec > 0.0:
            coord = coord.offset(
                np.random.rand() * 360 * geom.degrees,
                np.random.rand() * scatter_arcsec * geom.arcseconds)
        if indexer_ids is not None:
            src['pixelId'] = indexer_ids[src_idx]
        if associated_ids is not None:
            src['diaObjectId'] = associated_ids[src_idx]
        src.setCoord(coord)

        if wcs is not None:
            xyCentroid = wcs.skykToPixel(coord)
            src.set("x", xyCentroid.getX())
            src.set("y", xyCentroid.getY())

    sources = sources.asAstropy().to_pandas()

    return sources


class TestAssociationTask(unittest.TestCase):

    def setUp(self):
        """Create a sqlite3 database with default tables and schemas.
        """
        # CFHT Filters from the camera mapper.
        self.filter_names = ["u", "g", "r", "i", "z"]
        afwImageUtils.resetFilters()
        afwImageUtils.defineFilter('u', lambdaEff=374, alias="u.MP9301")
        afwImageUtils.defineFilter('g', lambdaEff=487, alias="g.MP9401")
        afwImageUtils.defineFilter('r', lambdaEff=628, alias="r.MP9601")
        afwImageUtils.defineFilter('i', lambdaEff=778, alias="i.MP9701")
        afwImageUtils.defineFilter('z', lambdaEff=1170, alias="z.MP9801")

        self.dia_object_schema = make_dia_object_schema()

        # metadata taken from CFHT data
        # v695856-e0/v695856-e0-c000-a00.sci_img.fits

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
            exposureId=1234,
            exposureTime=200.,
            date=dafBase.DateTime("2014-05-13T17:00:00.000000000",
                                  dafBase.DateTime.Timescale.TAI))
        self.exposure.setDetector(detector)
        self.exposure.getInfo().setVisitInfo(visit)
        self.exposure.setFilter(afwImage.Filter('g'))
        self.flux0 = 10000
        self.flux0_err = 100
        self.exposure.setPhotoCalib(
            afwImage.PhotoCalib(self.flux0, self.flux0_err))

        bbox = geom.Box2D(self.exposure.getBBox())
        wcs = self.exposure.getWcs()

        self.pixelator = sphgeom.HtmPixelization(20)
        region = sphgeom.ConvexPolygon([wcs.pixelToSky(pp).getVector()
                                        for pp in bbox.getCorners()])

        indices = self.pixelator.envelope(region, 64)
        # Index types must be cast to int to work with dax_apdb.
        self.index_ranges = indices.ranges()

    def tearDown(self):
        """Delete the database after we are done with it.
        """
        del self.metadata
        del self.wcs
        del self.exposure

    def test_run(self):
        """Test the run method with a database that already exists and
        contains DIAObjects and Sources.
        """
        dia_objects = self._run_association_and_retrieve_objects(True)
        not_updated_idx = 0
        updated_idx_start = 1
        new_idx_start = 5
        total_expected_dia_objects = 10
        self.assertEqual(len(dia_objects), total_expected_dia_objects)

        # Test to make sure the number of DIAObjects have been properly
        # associated within the db.
        for obj_idx, (df_idx, dia_object) in enumerate(dia_objects.iterrows()):
            if df_idx == not_updated_idx:
                # Test the DIAObject we expect to not be associated with any
                # new DIASources.
                self.assertEqual(dia_object['gPSFluxNdata'], 1)
                self.assertEqual(dia_object['rPSFluxNdata'], 1)
                self.assertEqual(dia_object['nDiaSources'], 2)
                self.assertEqual(df_idx, obj_idx)
            elif updated_idx_start <= df_idx < new_idx_start:
                # Test that associating to the existing DIAObjects went
                # as planned and test that the IDs of the newly associated
                # DIASources is correct.
                self.assertEqual(dia_object['gPSFluxNdata'], 2)
                self.assertEqual(dia_object['rPSFluxNdata'], 1)
                self.assertEqual(dia_object['nDiaSources'], 3)
                self.assertEqual(df_idx, obj_idx)
            else:
                self.assertEqual(dia_object['gPSFluxNdata'], 1)
                self.assertEqual(dia_object['nDiaSources'], 1)
                self.assertEqual(df_idx, obj_idx + 4 + 5)

    def test_run_no_existing_objects(self):
        """Test the run method with a completely empty database.
        """
        dia_objects = self._run_association_and_retrieve_objects(False)
        total_expected_dia_objects = 9
        self.assertEqual(len(dia_objects),
                         total_expected_dia_objects)
        for obj_idx, (df_idx, output_dia_object) in enumerate(dia_objects.iterrows()):
            self.assertEqual(output_dia_object['gPSFluxNdata'], 1)
            self.assertEqual(df_idx, obj_idx + 10)

    def _run_association_and_retrieve_objects(self, create_objects=False):
        """Convenience method for testing the Association run method.

        Parameters
        ----------
        create_objects : `bool`
            Boolean specifying if seed DIAObjects and DIASources should be
            inserted into the database before association.

        Return
        ------
        dia_objects : `lsst.afw.table.SourceCatalog`
            Final set of DIAObjects to be tested.
        """
        if create_objects:
            diaObjects, diaSourceHistory = \
                self._create_dia_objects_and_sources()
        else:
            diaObjects = pd.DataFrame(columns=["diaObjectId"])
            diaSourceHistory = pd.DataFrame(columns=["diaObjectId",
                                                     "filterName",
                                                     "diaSourceId"])
        diaObjects.set_index("diaObjectId",
                             inplace=True,
                             drop=False)
        diaSourceHistory.set_index(["diaObjectId",
                                    "filterName",
                                    "diaSourceId"],
                                   inplace=True,
                                   drop=False)

        source_centers = [
            [self.wcs.pixelToSky(idx, idx).getRa().asDegrees(),
             self.wcs.pixelToSky(idx, idx).getDec().asDegrees()]
            for idx in np.linspace(1, 1000, 10)[1:]]
        dia_sources = create_test_points(
            point_locs_deg=source_centers,
            start_id=10,
            scatter_arcsec=-1)
        for dia_source in dia_sources:
            self._set_source_values(
                dia_source=dia_source,
                flux=10000,
                fluxErr=100,
                # TODO DM-27170: fix this [0] workaround which gets a
                # single character representation of the band.
                filterName=self.exposure.getFilter().getCanonicalName()[0],
                ccdVisitId=self.exposure.getInfo().getVisitInfo().getExposureId(),
                midPointTai=self.exposure.getInfo().getVisitInfo().getDate().get(system=dafBase.DateTime.MJD))

        assoc_task = AssociationTask()

        diaSources = dia_sources.asAstropy().to_pandas()
        diaSources.rename(columns={"coord_ra": "ra",
                                   "coord_dec": "decl",
                                   "id": "diaSourceId",
                                   "parent": "parentDiaSourceId"},
                          inplace=True)
        diaSources["ra"] = np.degrees(diaSources["ra"])
        diaSources["decl"] = np.degrees(diaSources["decl"])

        if len(diaObjects) == 0:
            diaSourceHistory = pd.DataFrame(columns=["diaObjectId",
                                                     "filterName",
                                                     "diaSourceId"])
        diaSourceHistory.set_index(
            ["diaObjectId", "filterName", "diaSourceId"],
            drop=False,
            inplace=True)

        results = assoc_task.run(diaSources,
                                 diaObjects,
                                 diaSourceHistory)
        return results.diaObjects

    def _set_source_values(self, dia_source, flux, fluxErr, filterName,
                           ccdVisitId, midPointTai):
        """Set fluxes and visit info for DiaSources.

        Parameters
        ----------
        dia_source : `lsst.afw.table.SourceRecord`
            SourceRecord object to edit.
        flux : `double`
            Flux of DiaSource
        fluxErr : `double`
            Flux error of DiaSource
        filterName : `string`
            Name of filter for flux.
        ccdVisitId : `int`
            Integer id of this ccd/visit.
        midPointTai : `double`
            Time of observation
        """
        dia_source['ccdVisitId'] = ccdVisitId
        dia_source["midPointTai"] = midPointTai
        dia_source["psFlux"] = flux / self.flux0
        dia_source["psFluxErr"] = np.sqrt(
            (fluxErr / self.flux0) ** 2
            + (flux * self.flux0_err / self.flux0 ** 2) ** 2)
        dia_source["apFlux"] = flux / self.flux0
        dia_source["apFluxErr"] = np.sqrt(
            (fluxErr / self.flux0) ** 2
            + (flux * self.flux0_err / self.flux0 ** 2) ** 2)
        dia_source["totFlux"] = flux / self.flux0
        dia_source["totFluxErr"] = np.sqrt(
            (fluxErr / self.flux0) ** 2
            + (flux * self.flux0_err / self.flux0 ** 2) ** 2)
        dia_source["filterName"] = filterName
        dia_source["x"] = 0.
        dia_source["y"] = 0.

    def _create_dia_objects_and_sources(self):
        """Method for storing a set of test DIAObjects and sources into
        the L1 database.
        """

        # This should create a DB of 5 DIAObjects with 2 DIASources associated
        # to them. The DIASources are "observed" in g and r.

        # Create DIObjects, give them fluxes, and store them
        n_objects = 5
        object_centers = np.array([
            [self.wcs.pixelToSky(idx, idx).getRa().asDegrees(),
             self.wcs.pixelToSky(idx, idx).getDec().asDegrees()]
            for idx in np.linspace(1, 1000, 10)])
        dia_objects = create_test_points(
            point_locs_deg=object_centers[:n_objects],
            start_id=0,
            schema=self.dia_object_schema,
            scatter_arcsec=-1,)
        # Set the DIAObject fluxes and number of associated sources.
        for dia_object in dia_objects:
            dia_object["nDiaSources"] = 2
            for filter_name in self.filter_names:
                sphPoint = geom.SpherePoint(dia_object.getCoord())
                htmIndex = self.pixelator.index(sphPoint.getVector())
                dia_object["pixelId"] = htmIndex
                dia_object['%sPSFluxMean' % filter_name] = 1
                dia_object['%sPSFluxMeanErr' % filter_name] = 1
                dia_object['%sPSFluxSigma' % filter_name] = 1
                dia_object['%sPSFluxNdata' % filter_name] = 1
        dia_objects = dia_objects.asAstropy().to_pandas()
        dia_objects.rename(columns={"coord_ra": "ra",
                                    "coord_dec": "decl",
                                    "id": "diaObjectId"},
                           inplace=True)
        dia_objects["ra"] = np.degrees(dia_objects["ra"])
        dia_objects["decl"] = np.degrees(dia_objects["decl"])

        dateTime = dafBase.DateTime("2014-05-13T16:00:00.000000000",
                                    dafBase.DateTime.Timescale.TAI)

        # Create DIASources, update their ccdVisitId and fluxes, and store
        # them.
        dia_sources = create_test_points(
            point_locs_deg=np.concatenate(
                [object_centers[:n_objects], object_centers[:n_objects]]),
            start_id=0,
            scatter_arcsec=-1,
            associated_ids=[0, 1, 2, 3, 4,
                            0, 1, 2, 3, 4])
        for src_idx, dia_source in enumerate(dia_sources):
            if src_idx < n_objects:
                self._set_source_values(
                    dia_source=dia_source,
                    flux=10000,
                    fluxErr=100,
                    filterName='g',
                    ccdVisitId=1232,
                    midPointTai=dateTime.get(system=dafBase.DateTime.MJD))
            else:
                self._set_source_values(
                    dia_source=dia_source,
                    flux=10000,
                    fluxErr=100,
                    filterName='r',
                    ccdVisitId=1233,
                    midPointTai=dateTime.get(system=dafBase.DateTime.MJD))
        dia_sources = dia_sources.asAstropy().to_pandas()
        dia_sources.rename(columns={"coord_ra": "ra",
                                    "coord_dec": "decl",
                                    "id": "diaSourceId",
                                    "parent": "parentDiaSourceId"},
                           inplace=True)
        dia_sources["ra"] = np.degrees(dia_sources["ra"])
        dia_sources["decl"] = np.degrees(dia_sources["decl"])
        return dia_objects, dia_sources

    def test_associate_sources(self):
        """Test the performance of the associate_sources method in
        AssociationTask.
        """
        n_objects = 5
        dia_objects = create_test_points_pandas(
            point_locs_deg=[[0.04 * obj_idx, 0.04 * obj_idx]
                            for obj_idx in range(n_objects)],
            start_id=0,
            schema=self.dia_object_schema,
            scatter_arcsec=-1,)
        dia_objects.rename(columns={"coord_ra": "ra",
                                    "coord_dec": "decl",
                                    "id": "diaObjectId"},
                           inplace=True)

        n_sources = 5
        dia_sources = create_test_points_pandas(
            point_locs_deg=[
                [0.04 * (src_idx + 1),
                 0.04 * (src_idx + 1)]
                for src_idx in range(n_sources)],
            start_id=n_objects,
            scatter_arcsec=0.1)
        dia_sources.rename(columns={"coord_ra": "ra",
                                    "coord_dec": "decl",
                                    "id": "diaSourceId"},
                           inplace=True)

        assoc_task = AssociationTask()
        assoc_result = assoc_task.associate_sources(
            dia_objects, dia_sources)

        for test_obj_id, expected_obj_id in zip(
                assoc_result.associated_dia_object_ids,
                [1, 2, 3, 4, 9]):
            self.assertEqual(test_obj_id, expected_obj_id)

    def test_score_and_match(self):
        """Test association between a set of sources and an existing
        DIAObjectCollection.

        This also tests that a DIASource that can't be associated within
        tolerance is appended to the DIAObjectCollection as a new
        DIAObject.
        """

        assoc_task = AssociationTask()
        # Create a set of DIAObjects that contain only one DIASource
        n_objects = 5
        dia_objects = create_test_points_pandas(
            point_locs_deg=[[0.04 * obj_idx, 0.04 * obj_idx]
                            for obj_idx in range(n_objects)],
            start_id=0,
            schema=self.dia_object_schema,
            scatter_arcsec=-1,)
        dia_objects.rename(columns={"coord_ra": "ra",
                                    "coord_dec": "decl",
                                    "id": "diaObjectId"},
                           inplace=True)

        n_sources = 5
        dia_sources = create_test_points_pandas(
            point_locs_deg=[
                [0.04 * (src_idx + 1),
                 0.04 * (src_idx + 1)]
                for src_idx in range(n_sources)],
            start_id=n_objects,
            scatter_arcsec=-1)
        dia_sources.rename(columns={"coord_ra": "ra",
                                    "coord_dec": "decl",
                                    "id": "diaSourceId"},
                           inplace=True)

        score_struct = assoc_task.score(dia_objects,
                                        dia_sources,
                                        1.0 * geom.arcseconds)
        self.assertFalse(np.isfinite(score_struct.scores[-1]))
        for src_idx in range(4):
            # Our scores should be extremely close to 0 but not exactly so due
            # to machine noise.
            self.assertAlmostEqual(score_struct.scores[src_idx], 0.0,
                                   places=16)

        # After matching each DIAObject should now contain 2 DIASources
        # except the last DIAObject in this collection which should be
        # newly created during the matching step and contain only one
        # DIASource.
        match_result = assoc_task.match(dia_objects, dia_sources, score_struct)
        updated_ids = match_result.associated_dia_object_ids
        self.assertEqual(len(updated_ids), 5)
        self.assertEqual(match_result.n_updated_dia_objects, 4)
        self.assertEqual(match_result.n_new_dia_objects, 1)
        self.assertEqual(match_result.n_unassociated_dia_objects, 1)

        # Test updating all DiaObjects
        n_objects = 4
        dia_objects = create_test_points_pandas(
            point_locs_deg=[[0.04 * obj_idx, 0.04 * obj_idx]
                            for obj_idx in range(n_objects)],
            start_id=0,
            schema=self.dia_object_schema,
            scatter_arcsec=-1,)
        dia_objects.rename(columns={"coord_ra": "ra",
                                    "coord_dec": "decl",
                                    "id": "diaObjectId"},
                           inplace=True)

        n_sources = 4
        dia_sources = create_test_points_pandas(
            point_locs_deg=[
                [0.04 * src_idx,
                 0.04 * src_idx]
                for src_idx in range(n_sources)],
            start_id=n_objects,
            scatter_arcsec=-1)

        dia_sources.rename(columns={"coord_ra": "ra",
                                    "coord_dec": "decl",
                                    "id": "diaSourceId"},
                           inplace=True)
        score_struct = assoc_task.score(dia_objects[1:],
                                        dia_sources[:-1],
                                        1.0 * geom.arcseconds)
        match_result = assoc_task.match(dia_objects, dia_sources, score_struct)
        updated_ids = match_result.associated_dia_object_ids
        self.assertEqual(len(updated_ids), 4)

    def test_remove_nan_dia_sources(self):
        n_sources = 6
        dia_sources = create_test_points_pandas(
            point_locs_deg=[
                [0.04 * (src_idx + 1),
                 0.04 * (src_idx + 1)]
                for src_idx in range(n_sources)],
            start_id=0,
            scatter_arcsec=-1)
        dia_sources.rename(columns={"coord_ra": "ra",
                                    "coord_dec": "decl",
                                    "id": "diaSourceId"},
                           inplace=True)

        dia_sources.loc[2, "ra"] = np.nan
        dia_sources.loc[3, "decl"] = np.nan
        dia_sources.loc[4, "ra"] = np.nan
        dia_sources.loc[4, "decl"] = np.nan
        assoc_task = AssociationTask()
        out_dia_sources = assoc_task.check_dia_source_radec(dia_sources)
        self.assertEqual(len(out_dia_sources), n_sources - 3)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
