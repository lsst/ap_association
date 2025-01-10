#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#

"""Standalone pipelinetask to populate an APDB with simulated data.
"""

__all__ = ("TestApdbConfig",
           "TestApdbTask",
           )


import numpy as np
import pandas as pd

from lsst.daf.base import DateTime
import lsst.dax.apdb as daxApdb
import lsst.geom
from lsst.meas.base import DetectorVisitIdGeneratorConfig, IdGenerator
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.ap.association.loadDiaCatalogs import LoadDiaCatalogsTask
from lsst.ap.association.diaPipe import DiaPipelineTask
from lsst.ap.association.utils import (
    convertTableToSdmSchema,
    readSchemaFromApdb,
    column_dtype,
    make_empty_catalog,
    dropEmptyColumns,
)
import lsst.sphgeom


class TestApdbConnections(
        pipeBase.PipelineTaskConnections,
        dimensions=("instrument", "visit", "detector")):
    """Butler connections for TestApdbTask.
    """

    apdbTestMarker = pipeBase.connectionTypes.Output(
        doc="Marker dataset storing the configuration of the Apdb for each "
            "visit/detector. Used to signal the completion of the pipeline.",
        name="apdbTest_marker",
        storageClass="Config",
        dimensions=("instrument", "visit", "detector"),
    )


class TestApdbConfig(pipeBase.PipelineTaskConfig,
                     pipelineConnections=TestApdbConnections):
    """Config for TestApdbTask.
    """
    apdb_config_url = pexConfig.Field(
        dtype=str,
        default=None,
        optional=False,
        doc="A config file specifying the APDB and its connection parameters, "
            "typically written by the apdb-cli command-line utility. "
            "The database must already be initialized.",
    )
    survey_area = pexConfig.Field(
        dtype=float,
        default=20000,
        doc="Area (in degrees) of the simulated survey",
    )
    fov = pexConfig.Field(
        dtype=float,
        default=9.6,
        doc="Field of view of the camera, in square degrees.",
    )
    stellar_density = pexConfig.Field(
        dtype=float,
        default=1750,
        doc="Average number of real transient and variable objects"
            " detected per square degree. For this simulation, these will"
            " always be detected, and detected in the same location."
            "The default is chosen such that:"
            " Density x Rubin Fov (9.6) x # of visits per night (~600) ~ 10M",
    )
    false_positive_ratio = pexConfig.Field(
        dtype=float,
        default=4,
        doc="Average ratio of false detections to real sources."
            "These will be detected in random locations.",
    )
    false_positive_variability = pexConfig.Field(
        dtype=float,
        default=100,
        doc="Parameter characterizing the variability in the rate of false"
            " positives.",
    )
    sky_seed = pexConfig.Field(
        dtype=int,
        default=37,
        doc="Seed used to simulate the real sources.",
    )
    historyThreshold = pexConfig.Field(
        dtype=int,
        doc="Minimum number of detections of a diaObject required "
            "to run forced photometry. Set to 1 to include all diaObjects.",
        default=2,
    )
    objId_start = pexConfig.Field(
        dtype=int,
        default=1,
        doc="Starting diaObject ID number of real objects.",
    )
    maximum_table_length = pexConfig.Field(
        dtype=int,
        default=65535,
        doc="Maximum length of tables allowed to be written in one operation"
            " to the Cassandra APDB",
    )

    idGenerator = DetectorVisitIdGeneratorConfig.make_field()
    idGeneratorFakes = DetectorVisitIdGeneratorConfig.make_field()
    idGeneratorForced = DetectorVisitIdGeneratorConfig.make_field()

    def setDefaults(self):
        self.idGenerator = DetectorVisitIdGeneratorConfig(release_id=0, n_releases=3)
        self.idGeneratorFakes = DetectorVisitIdGeneratorConfig(release_id=1, n_releases=3)
        self.idGeneratorForced = DetectorVisitIdGeneratorConfig(release_id=2, n_releases=3)


class TestApdbTask(LoadDiaCatalogsTask):
    """Task for loading, associating and storing Difference Image Analysis
    (DIA) Objects and Sources.
    """
    ConfigClass = TestApdbConfig
    _DefaultName = "apdbTest"

    def __init__(self, initInputs=None, **kwargs):
        super().__init__(**kwargs)
        self.apdb = daxApdb.Apdb.from_uri(self.config.apdb_config_url)
        self.schema = readSchemaFromApdb(self.apdb)
        self.prepareSurvey()

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        inputs["visit"] = butlerQC.quantum.dataId["visit"]
        inputs["detector"] = butlerQC.quantum.dataId["detector"]
        inputs["idGenerator"] = self.config.idGenerator.apply(butlerQC.quantum.dataId)
        inputs["idGeneratorFakes"] = self.config.idGeneratorFakes.apply(butlerQC.quantum.dataId)
        inputs["idGeneratorForced"] = self.config.idGeneratorForced.apply(butlerQC.quantum.dataId)
        self.run(**inputs)
        # Note that the commented-out code below is intentionally omitted
        # We don't want this to write anything to the Butler!
        butlerQC.put(pipeBase.Struct(), outputRefs)

    def prepareSurvey(self):
        """Prepare survey parameters common to all observations
        """
        deg2ToSteradian = np.deg2rad(1)**2  # Conversion factor from square degrees to steradians
        sphereDegrees = 4*np.pi/deg2ToSteradian  # area of a unit sphere in square degrees
        sphereFraction = self.config.survey_area/sphereDegrees
        survey_area = self.config.survey_area*deg2ToSteradian
        # Avoid the branch cut of the arcsin function at the equator
        # Calculate the angle from the pole needed to cover `survey_area` square
        #  degrees of a sphere
        if sphereFraction < 0.5:
            declinationMax = np.arcsin(1 - survey_area/(2*np.pi))
        else:
            declinationMax = np.pi - np.arcsin(1 - (4*np.pi - survey_area)/(2*np.pi))
        self.radiusMax = stereographicRaDec2XY(0, declinationMax)[0]
        self.nReal = int(self.config.survey_area*self.config.stellar_density)
        # Radius of the focal plane, in radians
        self.fpRadius = np.radians(np.sqrt(self.config.fov/np.pi))

    def run(self, visit, detector,
            idGenerator=IdGenerator(),
            idGeneratorFakes=IdGenerator(),
            idGeneratorForced=IdGenerator()):
        """Generate a full focal plane simulation of real and fake sources,
        run association and write to the APDB

        Parameters
        ----------
        visit : `int`
            Visit ID.
        detector : `int`
            Detector number, used for the ID generator.
            Note that each instance simulates an entire focal plane.
            The detector dimension is used as a convenience to allow re-using
            code and infrastructure.
        idGenerator : `lsst.meas.base,IdGenerator`, optional
            ID generator used for "real" sources.
        idGeneratorFakes : `lsst.meas.base,IdGenerator`, optional
            ID generator for random "fake" sources.
            Separate ID generators are needed to avoid overflows between visits,
            since an entire focal plane is simulated per instance.
        idGeneratorForced : `lsst.meas.base,IdGenerator`, optional
            ID generator for forced sources at existing diaObject locations.
        """
        idGen = idGenerator.make_table_id_factory()
        idGenFakes = idGeneratorFakes.make_table_id_factory()
        idGenForced = idGeneratorForced.make_table_id_factory()
        seed = int(visit*1000 + detector)
        x, y = randomCircleXY(self.radiusMax, seed)
        ra, dec = stereographicXY2RaDec(x, y)
        if ra < 0:
            ra += 2*np.pi
        delta0, _ = stereographicRaDec2XY(0, dec - self.fpRadius)
        delta1, _ = stereographicRaDec2XY(0, dec + self.fpRadius)
        delta = abs(delta1 - delta0)/2.
        x0 = x - delta
        x1 = x + delta
        y0 = y - delta
        y1 = y + delta
        ra0, dec0 = stereographicXY2RaDec(x0, y0)
        ra1, dec1 = stereographicXY2RaDec(x1, y1)
        region = lsst.sphgeom.ConvexPolygon([lsst.geom.SpherePoint(ra0, dec0, lsst.geom.radians).getVector(),
                                            lsst.geom.SpherePoint(ra0, dec1, lsst.geom.radians).getVector(),
                                            lsst.geom.SpherePoint(ra1, dec1, lsst.geom.radians).getVector(),
                                            lsst.geom.SpherePoint(ra1, dec0, lsst.geom.radians).getVector()]
                                            )
        self.log.info(f"Simulating {self.nReal} sources around RA={np.degrees(ra)}, Dec={np.degrees(dec)}")
        xS, yS = randomCircleXY(self.radiusMax, self.config.sky_seed, n=self.nReal)
        diaObjIds = np.arange(self.config.objId_start, self.nReal + self.config.objId_start)
        inds = (xS > x0) & (xS < x1) & (yS > y0) & (yS < y1)
        nSim = np.sum(inds)
        diaObjIds = diaObjIds[inds]
        self.log.info(f"{nSim} sources made the spatial cut")
        diaSourcesReal = self.createDiaSources(*stereographicXY2RaDec(xS[inds], yS[inds]),
                                               idGenerator=idGen,
                                               diaObjectIds=diaObjIds)

        rng = np.random.RandomState(seed)
        scale = nSim*self.config.false_positive_ratio/self.config.false_positive_variability
        nBogus = int(rng.standard_gamma(scale)*self.config.false_positive_variability)

        diaSourcesBogus = self.createDiaSources(*self.generateFalseDetections(x0, x1, y0, y1, nBogus, seed),
                                                idGenerator=idGenFakes)
        diaSourcesRaw = pd.concat([diaSourcesReal, diaSourcesBogus])

        diaSources = convertTableToSdmSchema(self.schema, diaSourcesRaw, tableName="DiaSource")

        diaObjects = self.loadDiaObjects(region, self.schema)

        # Associate DiaSources with DiaObjects
        associatedDiaSources, newDiaObjects = self.associateDiaSources(diaSources, diaObjects)
        # Merge new and preloaded diaObjects
        mergedDiaObjects = self.mergeAssociatedCatalogs(diaObjects, newDiaObjects)

        nObj = len(mergedDiaObjects)
        nSrc = len(associatedDiaSources)
        dateTime = DateTime.now().toAstropy()
        ind = 0
        # Note that nObj must always be equal to or greater than nSrc
        for start in range(0, nObj, self.config.maximum_table_length):
            end = min(start + self.config.maximum_table_length, nObj)
            diaObjectsChunk = mergedDiaObjects.iloc[start:end]
            self.log.info(f"Writing diaObject chunk {ind} of length {len(diaObjectsChunk)} to the APDB")
            srcEnd = min(start + self.config.maximum_table_length, nSrc)
            if srcEnd <= start:
                finalDiaSources = None
            else:
                diaSourcesChunk = associatedDiaSources.iloc[start:srcEnd]
                finalDiaSources = convertTableToSdmSchema(self.schema, diaSourcesChunk, tableName="DiaSource")
                self.log.info(f"Writing diaOSource chunk {ind} of length {len(diaSourcesChunk)} to the APDB")
            diaForcedSources = self.runForcedMeasurement(diaObjectsChunk, idGenForced, visit, detector)

            finalDiaObjects = convertTableToSdmSchema(self.schema, diaObjectsChunk, tableName="DiaObject")
            finalDiaForcedSources = convertTableToSdmSchema(self.schema, diaForcedSources,
                                                            tableName="DiaForcedSource")
            self.writeToApdb(finalDiaObjects, finalDiaSources, finalDiaForcedSources, dateTime)
            ind += 1
        marker = pexConfig.Config()
        return pipeBase.Struct(apdbTestMarker=marker)

    def createDiaSources(self, raVals, decVals, idGenerator, diaObjectIds=None):
        """Create diaSources with the supplied coordinates.

        Parameters
        ----------
        raVals : `numpy.ndarray`
            Right Ascension (RA) values of the desired sources.
        decVals : `numpy.ndarray`
            Declination (Dec) values of the desired sources.
        idGenerator : `lsst.meas.base,IdGenerator`, optional
            ID generator used for the source IDs
        diaObjectIds : `numpy.ndarray`, optional
            diaObjectIds to assign to the sources. If not supplied, the object
            ID is set to the diaSourceId

        Returns
        -------
        diaSources : `pandas.DataFrame`
            Table of sources with the supplied coordinates and IDs assigned.
        """
        n = len(raVals)
        raVals[raVals < 0] += 2*np.pi
        ra = pd.Series(np.degrees(raVals), name='ra')
        dec = pd.Series(np.degrees(decVals), name='dec')
        diaSourceId = pd.Series([idGenerator() for i in range(n)], name='diaSourceId')
        if diaObjectIds is None:
            diaObjectId = pd.Series(diaSourceId, name='diaObjectId')
        else:
            diaObjectId = pd.Series(diaObjectIds, name='diaObjectId')
        diaSources = pd.concat([diaSourceId, diaObjectId, ra, dec], axis=1)
        # Do *not* set the index to the diaSourceId. It will need to be diaObjectId for matching (later)
        # diaSources.set_index("diaSourceId", inplace=True)
        return diaSources

    def generateFalseDetections(self, x0, x1, y0, y1, nBogus, seed):
        """Generate random coordinates within the ranges provided.

        Parameters
        ----------
        x0 : `int`
            Minimum projected x coordinate
        x1 : `int`
            Maximum projected x coordinate
        y0 : `int`
            Minimum projected y coordinate
        y1 : `int`
            Maximum projected y coordinate
        nBogus : `int`
            Number of "fake" sources within the region.
        rng : `int`
            Seed value for the random number generator to provide repeatable results.

        Returns
        -------
        ra, dec : `numpy.ndarray`
            Coordinates matching the randomly generated locations.
        """
        self.log.info(f"Simulating {nBogus} false detections within region.")
        rng = np.random.RandomState(seed)
        x = rng.random_sample(nBogus)*(x1 - x0) + x0
        y = rng.random_sample(nBogus)*(y1 - y0) + y0
        return stereographicXY2RaDec(x, y)

    def simpleMatch(self, diaSourceTable, diaObjects):
        """Match by pre-defined ID.

        Parameters
        ----------
        diaSourceTable : `pandas.DataFrame`
            New DIASources to be associated with existing DIAObjects.
        diaObjects : `pandas.DataFrame`
            Existing diaObjects from the Apdb.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components.

            - ``matchedDiaSources`` : DiaSources that were matched. Matched
              Sources have their diaObjectId updated and set to the id of the
              diaObject they were matched to. (`pandas.DataFrame`)
            - ``unAssocDiaSources`` : DiaSources that were not matched.
              Unassociated sources have their diaObject set to 0 as they
              were not associated with any existing DiaObjects.
              (`pandas.DataFrame`)
            - ``nUpdatedDiaObjects`` : Number of DiaObjects that were
              matched to new DiaSources. (`int`)
            - ``nUnassociatedDiaObjects`` : Number of DiaObjects that were
              not matched a new DiaSource. (`int`)
        """
        # Only use diaObjectId as the index here, and do not update the index
        # of the original diaSourceTable
        diaSourceTable.set_index("diaObjectId", inplace=True)
        matchedDiaSources = diaSourceTable.loc[diaSourceTable.index.intersection(diaObjects.index)]
        unAssocDiaSources = diaSourceTable.loc[diaSourceTable.index.difference(diaObjects.index)]

        matchedDiaObjectInds = diaObjects.index.intersection(diaSourceTable.index)
        if not diaObjects.empty:
            nDiaSources = diaObjects.loc[matchedDiaObjectInds, "nDiaSources"] + 1
            diaObjects.loc[matchedDiaObjectInds, "nDiaSources"] = nDiaSources

        # Reset the index of the diaSource dataframes, so diaObjectId is still
        #  a valid column
        return pipeBase.Struct(unAssocDiaSources=unAssocDiaSources.reset_index(),
                               matchedDiaSources=matchedDiaSources.reset_index(),
                               nUpdatedDiaObjects=len(matchedDiaSources),
                               nUnassociatedDiaObjects=len(unAssocDiaSources),
                               )

    def associateDiaSources(self, diaSourceTable, diaObjects):
        """Associate DiaSources with DiaObjects.

        Parameters
        ----------
        diaSourceTable : `pandas.DataFrame`
            Newly detected DiaSources.
        diaObjects : `pandas.DataFrame`
            Table of  DiaObjects from preloaded DiaObjects.

        Returns
        -------
        associatedDiaSources : `pandas.DataFrame`
            Associated DiaSources with DiaObjects.
        newDiaObjects : `pandas.DataFrame`
            Table of new DiaObjects after association.
        """
        # Associate new DiaSources with existing DiaObjects.
        assocResults = self.simpleMatch(diaSourceTable, diaObjects)

        toAssociate = []

        # Create new DiaObjects from unassociated diaSources.
        createResults = self.createNewDiaObjects(assocResults.unAssocDiaSources)
        if len(assocResults.matchedDiaSources) > 0:
            toAssociate.append(assocResults.matchedDiaSources)
        toAssociate.append(createResults.diaSources)
        associatedDiaSources = pd.concat(toAssociate)

        self.log.info("%i updated and %i unassociated diaObjects. Creating %i new diaObjects",
                      assocResults.nUpdatedDiaObjects,
                      assocResults.nUnassociatedDiaObjects,
                      createResults.nNewDiaObjects,
                      )

        # Index the DiaSource catalog for this visit after all associations
        # have been made.
        associatedDiaSources.set_index(["diaObjectId",
                                        "band",
                                        "diaSourceId"],
                                       drop=False,
                                       inplace=True)
        return (associatedDiaSources, createResults.newDiaObjects)

    def createNewDiaObjects(self, unAssocDiaSources):
        """Loop through the set of DiaSources and create new DiaObjects
        for unassociated DiaSources.

        Parameters
        ----------
        unAssocDiaSources : `pandas.DataFrame`
            Set of DiaSources to create new DiaObjects from.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            Results struct containing:

            - diaSources : `pandas.DataFrame`
                DiaSource catalog with updated DiaObject ids.
            - newDiaObjects : `pandas.DataFrame`
                Newly created DiaObjects from the unassociated DiaSources.
            - nNewDiaObjects : `int`
                Number of newly created diaObjects.
        """
        if len(unAssocDiaSources) == 0:
            newDiaObjects = make_empty_catalog(self.schema, tableName="DiaObject")
        else:
            # Do *not* set the diaObjectId to the diaSourceId.
            # For this simulation we are using custom diaObjectIds,
            #  and need to preserve them.
            # unAssocDiaSources["diaObjectId"] = unAssocDiaSources["diaSourceId"]
            newDiaObjects = convertTableToSdmSchema(self.schema, unAssocDiaSources,
                                                    tableName="DiaObject")
            newDiaObjects.nDiaSources = 1
        return pipeBase.Struct(diaSources=unAssocDiaSources,
                               newDiaObjects=newDiaObjects,
                               nNewDiaObjects=len(newDiaObjects))

    def mergeAssociatedCatalogs(self, diaObjects, newDiaObjects):
        """Merge the associated diaObjects to their previous history.
        Also update the index of associatedDiaSources in place.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            Table of  DiaObjects from preloaded DiaObjects.
        newDiaObjects : `pandas.DataFrame`
            Table of new DiaObjects after association.

        Returns
        -------
        mergedDiaObjects : `pandas.DataFrame`
            Table of new DiaObjects merged with their history.

        Raises
        ------
        RuntimeError
            Raised if duplicate DiaObjects are found.
        """

        # Append new DiaObjects to their previous history.
        # Do not modify diaSources
        newDiaObjects.set_index("diaObjectId", drop=False, inplace=True)
        if diaObjects.empty:
            mergedDiaObjects = newDiaObjects
        elif not newDiaObjects.empty:
            mergedDiaObjects = pd.concat([diaObjects, newDiaObjects], sort=True)
        else:
            mergedDiaObjects = diaObjects
        if DiaPipelineTask.testDataFrameIndex(mergedDiaObjects):
            raise RuntimeError("Duplicate DiaObjects created after association.")
        return mergedDiaObjects

    def runForcedMeasurement(self, diaObjects, idGenerator, visit, detector):
        """Forced Source Measurement

        Forced photometry on the difference and calibrated exposures using the
        new and updated DiaObject locations.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            Catalog of DiaObjects.
        idGenerator : `lsst.meas.base.IdGenerator`
            Object that generates source IDs and random number generator seeds.
        visit : `int`
            Visit ID.
        detector : `int`
            Detector number, used for the ID generator.
            Note that each instance simulates an entire focal plane.
            The detector dimension is used as a convenience to allow re-using
            code and infrastructure.

        Returns
        -------
        diaForcedSources : `pandas.DataFrame`
            Catalog of calibrated forced photometered fluxes on both the
            difference and direct images at DiaObject locations.
        """
        # Restrict forced source measurement to objects with sufficient history to be reliable.
        objectTable = diaObjects.query(f'nDiaSources >= {self.config.historyThreshold}')
        preserveColumns = ["diaObjectId", "ra", "dec"]
        baseForcedSources = objectTable[preserveColumns].copy()
        baseForcedSources["visit"] = visit
        preserveColumns.append("visit")
        baseForcedSources["detector"] = detector
        preserveColumns.append("detector")
        baseForcedSources["diaForcedSourceId"] = [idGenerator() for i in range(len(baseForcedSources))]
        preserveColumns.append("diaForcedSourceId")
        # Fill the forced sources for each diaObject with random data
        diaForcedSources = fillRandomTable(self.schema, baseForcedSources,
                                           tableName="DiaForcedSource",
                                           preserveColumns=preserveColumns)
        self.log.info(f"Updating {len(diaForcedSources)} diaForcedSources in the APDB")
        diaForcedSources = convertTableToSdmSchema(self.schema, diaForcedSources,
                                                   tableName="DiaForcedSource",
                                                   )
        return diaForcedSources

    def writeToApdb(self, updatedDiaObjects, associatedDiaSources, diaForcedSources, dateTime):
        """Write to the Alert Production Database (Apdb).

        Store DiaSources, updated DiaObjects, and DiaForcedSources in the
        Alert Production Database (Apdb).

        Parameters
        ----------
        updatedDiaObjects : `pandas.DataFrame`
            Catalog of updated DiaObjects.
        associatedDiaSources : `pandas.DataFrame`
            Associated DiaSources with DiaObjects.
        diaForcedSources : `pandas.DataFrame`
            Catalog of calibrated forced photometered fluxes on both the
            difference and direct images at DiaObject locations.
        """
        # Store DiaSources, updated DiaObjects, and DiaForcedSources in the
        # Apdb.
        # Drop empty columns that are nullable in the APDB.
        diaObjectStore = dropEmptyColumns(self.schema, updatedDiaObjects, tableName="DiaObject")
        if associatedDiaSources is None:
            diaSourceStore = None
        else:
            diaSourceStore = dropEmptyColumns(self.schema, associatedDiaSources, tableName="DiaSource")
        diaForcedSourceStore = dropEmptyColumns(self.schema, diaForcedSources, tableName="DiaForcedSource")
        self.apdb.store(
            dateTime,
            diaObjectStore,
            diaSourceStore,
            diaForcedSourceStore)
        self.log.info("APDB updated.")


def randomCircleXY(radius, seed, n=None):
    """Generate random x, y coordinates within a circular region

    Parameters
    ----------
    radius : `float`
        Radius of the circular region.
    seed : `int`
        Seed for the random number generator.
    n : `int`, optional
        Number of points to generate.

    Returns
    -------
    x, y : `float`, or `numpy.ndarray`
        Random coordinates inside the region.
    """
    # Draw the locations from a regular grid in x,y, but use polar
    # coordinates to maintain a circular region in this space.
    rng = np.random.RandomState(seed)
    r2 = rng.random_sample(n)*radius**2
    phi = rng.random_sample(n)*2*np.pi
    x = np.sqrt(r2)*np.cos(phi)
    y = np.sqrt(r2)*np.sin(phi)
    return (x, y)


def stereographicXY2RaDec(x, y):
    """Convert from a grid-like stereographic projection to RA and Dec

    Notes
    -----
    A stereographic projection centered on the pole
    From Eq 56 of Calabretta and Greisen (2002)
    "Representations of celestial coordinates in FITS"

    Parameters
    ----------
    x : `float`
        Column index of the stereographic grid.
    y : `float`
        Row index of the stereographic grid.

    Returns
    -------
    ra, dec : `float`
        Right Ascension and Declination.
    """
    r = np.sqrt(x**2 + y**2)
    dec = np.pi/2 - 2*np.arctan(r/2)
    ra = np.arctan2(y, x)
    return (ra, dec)


def stereographicRaDec2XY(ra, dec):
    """Convert from Ra, Dec to a grid-like stereographic projection

    Notes
    -----
    A stereographic projection centered on the pole
    From Eq 56 of Calabretta and Greisen (2002)
    "Representations of celestial coordinates in FITS"

    Parameters
    ----------
    ra : `float`
        Right Ascension (RA) in radians.
    dec : `float`
        Declination (Dec) in radians.

    Returns
    -------
    x, y : `float`
        Grid coordinates of the stereographic projection
    """
    r = 2*np.cos(dec)/(1 + np.sin(dec))
    x = r*np.cos(ra)
    y = r*np.sin(ra)
    return (x, y)


def fillRandomTable(apdbSchema, sourceTable, tableName, preserveColumns=None):
    """Force a table to conform to the schema defined by the APDB.

    This method uses the table definitions in ``sdm_schemas`` to
    load the schema of the APDB, and does not actually connect to the APDB.

    Parameters
    ----------
    apdbSchema : `dict` [`str`, `lsst.dax.apdb.schema_model.Table`]
        Schema from ``sdm_schemas`` containing the table definition to use.
    sourceTable : `pandas.DataFrame`
        The input table to convert.
    tableName : `str`
        Name of the table in the schema to use.
    preserveColumns : `list` of `str`, optional
        List of columns to copy from the input sourceTable.

    Returns
    -------
    `pandas.DataFrame`
        A table with the correct schema for the APDB and data copied from
        the input ``sourceTable``.
    """
    table = apdbSchema[tableName]

    data = {}
    nSrc = len(sourceTable)
    rng = np.random.default_rng()
    if preserveColumns is None:
        preserveColumns = []

    for columnDef in table.columns:
        if columnDef.name in preserveColumns:
            data[columnDef.name] = sourceTable[columnDef.name]
        dtype = column_dtype(columnDef.datatype)
        if columnDef.name in sourceTable.columns:
            data[columnDef.name] = pd.Series(sourceTable[columnDef.name], dtype=dtype)
        else:
            try:
                dataInit = rng.random(nSrc).astype(dtype)
            except (TypeError, ValueError):
                dataInit = np.zeros(nSrc, dtype=column_dtype(columnDef.datatype))
            data[columnDef.name] = pd.Series(dataInit, index=sourceTable.index)
    return pd.DataFrame(data)
