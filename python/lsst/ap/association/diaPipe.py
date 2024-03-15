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

"""PipelineTask for associating DiaSources with previous DiaObjects.

Additionally performs forced photometry on the calibrated and difference
images at the updated locations of DiaObjects.

Currently loads directly from the Apdb rather than pre-loading.
"""

__all__ = ("DiaPipelineConfig",
           "DiaPipelineTask",
           "DiaPipelineConnections")

import numpy as np
import pandas as pd

from lsst.daf.base import DateTime
import lsst.dax.apdb as daxApdb
from lsst.meas.base import DetectorVisitIdGeneratorConfig, DiaObjectCalculationTask
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as connTypes
from lsst.utils.timer import timeMethod

from lsst.ap.association import (
    AssociationTask,
    DiaForcedSourceTask,
    LoadDiaCatalogsTask,
    PackageAlertsTask)
from lsst.ap.association.ssoAssociation import SolarSystemAssociationTask


class DiaPipelineConnections(
        pipeBase.PipelineTaskConnections,
        dimensions=("instrument", "visit", "detector"),
        defaultTemplates={"coaddName": "deep", "fakesType": ""}):
    """Butler connections for DiaPipelineTask.
    """
    diaSourceTable = connTypes.Input(
        doc="Catalog of calibrated DiaSources.",
        name="{fakesType}{coaddName}Diff_diaSrcTable",
        storageClass="DataFrame",
        dimensions=("instrument", "visit", "detector"),
    )
    solarSystemObjectTable = connTypes.Input(
        doc="Catalog of SolarSolarSystem objects expected to be observable in "
            "this detectorVisit.",
        name="visitSsObjects",
        storageClass="DataFrame",
        dimensions=("instrument", "visit"),
    )
    diffIm = connTypes.Input(
        doc="Difference image on which the DiaSources were detected.",
        name="{fakesType}{coaddName}Diff_differenceExp",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector"),
    )
    exposure = connTypes.Input(
        doc="Calibrated exposure differenced with a template image during "
            "image differencing.",
        name="{fakesType}calexp",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector"),
    )
    template = connTypes.Input(
        doc="Warped template used to create `subtractedExposure`. Not PSF "
            "matched.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="ExposureF",
        name="{fakesType}{coaddName}Diff_templateExp",
    )
    apdbMarker = connTypes.Output(
        doc="Marker dataset storing the configuration of the Apdb for each "
            "visit/detector. Used to signal the completion of the pipeline.",
        name="apdb_marker",
        storageClass="Config",
        dimensions=("instrument", "visit", "detector"),
    )
    associatedDiaSources = connTypes.Output(
        doc="Optional output storing the DiaSource catalog after matching, "
            "calibration, and standardization for insertion into the Apdb.",
        name="{fakesType}{coaddName}Diff_assocDiaSrc",
        storageClass="DataFrame",
        dimensions=("instrument", "visit", "detector"),
    )
    diaForcedSources = connTypes.Output(
        doc="Optional output storing the forced sources computed at the diaObject positions.",
        name="{fakesType}{coaddName}Diff_diaForcedSrc",
        storageClass="DataFrame",
        dimensions=("instrument", "visit", "detector"),
    )
    diaObjects = connTypes.Output(
        doc="Optional output storing the updated diaObjects associated to these sources.",
        name="{fakesType}{coaddName}Diff_diaObject",
        storageClass="DataFrame",
        dimensions=("instrument", "visit", "detector"),
    )
    longTrailedSources = pipeBase.connectionTypes.Output(
        doc="Optional output temporarily storing long trailed diaSources.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="DataFrame",
        name="{fakesType}{coaddName}Diff_longTrailedSrc",
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not config.doWriteAssociatedSources:
            self.outputs.remove("associatedDiaSources")
            self.outputs.remove("diaForcedSources")
            self.outputs.remove("diaObjects")
        if not config.doSolarSystemAssociation:
            self.inputs.remove("solarSystemObjectTable")
        if not config.associator.doTrailedSourceFilter:
            self.outputs.remove("longTrailedSources")

    def adjustQuantum(self, inputs, outputs, label, dataId):
        """Override to make adjustments to `lsst.daf.butler.DatasetRef` objects
        in the `lsst.daf.butler.core.Quantum` during the graph generation stage
        of the activator.

        This implementation checks to make sure that the filters in the dataset
        are compatible with AP processing as set by the Apdb/DPDD schema.

        Parameters
        ----------
        inputs : `dict`
            Dictionary whose keys are an input (regular or prerequisite)
            connection name and whose values are a tuple of the connection
            instance and a collection of associated `DatasetRef` objects.
            The exact type of the nested collections is unspecified; it can be
            assumed to be multi-pass iterable and support `len` and ``in``, but
            it should not be mutated in place.  In contrast, the outer
            dictionaries are guaranteed to be temporary copies that are true
            `dict` instances, and hence may be modified and even returned; this
            is especially useful for delegating to `super` (see notes below).
        outputs : `dict`
            Dict of output datasets, with the same structure as ``inputs``.
        label : `str`
            Label for this task in the pipeline (should be used in all
            diagnostic messages).
        data_id : `lsst.daf.butler.DataCoordinate`
            Data ID for this quantum in the pipeline (should be used in all
            diagnostic messages).

        Returns
        -------
        adjusted_inputs : `dict`
            Dict of the same form as ``inputs`` with updated containers of
            input `DatasetRef` objects.  Connections that are not changed
            should not be returned at all.  Datasets may only be removed, not
            added.  Nested collections may be of any multi-pass iterable type,
            and the order of iteration will set the order of iteration within
            `PipelineTask.runQuantum`.
        adjusted_outputs : `dict`
            Dict of updated output datasets, with the same structure and
            interpretation as ``adjusted_inputs``.

        Raises
        ------
        ScalarError
            Raised if any `Input` or `PrerequisiteInput` connection has
            ``multiple`` set to `False`, but multiple datasets.
        NoWorkFound
            Raised to indicate that this quantum should not be run; not enough
            datasets were found for a regular `Input` connection, and the
            quantum should be pruned or skipped.
        FileNotFoundError
            Raised to cause QuantumGraph generation to fail (with the message
            included in this exception); not enough datasets were found for a
            `PrerequisiteInput` connection.
        """
        _, refs = inputs["diffIm"]
        for ref in refs:
            if ref.dataId["band"] not in self.config.validBands:
                raise ValueError(
                    f"Requested '{ref.dataId['band']}' not in "
                    "DiaPipelineConfig.validBands. To process bands not in "
                    "the standard Rubin set (ugrizy) you must add the band to "
                    "the validBands list in DiaPipelineConfig and add the "
                    "appropriate columns to the Apdb schema.")
        return super().adjustQuantum(inputs, outputs, label, dataId)


class DiaPipelineConfig(pipeBase.PipelineTaskConfig,
                        pipelineConnections=DiaPipelineConnections):
    """Config for DiaPipelineTask.
    """
    coaddName = pexConfig.Field(
        doc="coadd name: typically one of deep, goodSeeing, or dcr",
        dtype=str,
        default="deep",
    )
    apdb = daxApdb.ApdbSql.makeField(
        doc="Database connection for storing associated DiaSources and "
            "DiaObjects. Must already be initialized.",
    )
    validBands = pexConfig.ListField(
        dtype=str,
        default=["u", "g", "r", "i", "z", "y"],
        doc="List of bands that are valid for AP processing. To process a "
            "band not on this list, the appropriate band specific columns "
            "must be added to the Apdb schema in dax_apdb.",
    )
    diaCatalogLoader = pexConfig.ConfigurableField(
        target=LoadDiaCatalogsTask,
        doc="Task to load DiaObjects and DiaSources from the Apdb.",
    )
    associator = pexConfig.ConfigurableField(
        target=AssociationTask,
        doc="Task used to associate DiaSources with DiaObjects.",
    )
    doSolarSystemAssociation = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Process SolarSystem objects through the pipeline.",
    )
    solarSystemAssociator = pexConfig.ConfigurableField(
        target=SolarSystemAssociationTask,
        doc="Task used to associate DiaSources with SolarSystemObjects.",
    )
    diaCalculation = pexConfig.ConfigurableField(
        target=DiaObjectCalculationTask,
        doc="Task to compute summary statistics for DiaObjects.",
    )
    diaForcedSource = pexConfig.ConfigurableField(
        target=DiaForcedSourceTask,
        doc="Task used for force photometer DiaObject locations in direct and "
            "difference images.",
    )
    alertPackager = pexConfig.ConfigurableField(
        target=PackageAlertsTask,
        doc="Subtask for packaging Ap data into alerts.",
    )
    doPackageAlerts = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Package Dia-data into serialized alerts for distribution and "
            "write them to disk.",
    )
    doWriteAssociatedSources = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Write out associated DiaSources, DiaForcedSources, and DiaObjects, "
            "formatted following the Science Data Model.",
    )
    imagePixelMargin = pexConfig.RangeField(
        dtype=int,
        default=10,
        min=0,
        doc="Pad the image by this many pixels before removing off-image "
            "diaObjects for association.",
    )
    idGenerator = DetectorVisitIdGeneratorConfig.make_field()

    def setDefaults(self):
        self.apdb.dia_object_index = "baseline"
        self.apdb.dia_object_columns = []
        self.diaCalculation.plugins = ["ap_meanPosition",
                                       "ap_nDiaSources",
                                       "ap_diaObjectFlag",
                                       "ap_meanFlux",
                                       "ap_percentileFlux",
                                       "ap_sigmaFlux",
                                       "ap_chi2Flux",
                                       "ap_madFlux",
                                       "ap_skewFlux",
                                       "ap_minMaxFlux",
                                       "ap_maxSlopeFlux",
                                       "ap_meanErrFlux",
                                       "ap_linearFit",
                                       "ap_stetsonJ",
                                       "ap_meanTotFlux",
                                       "ap_sigmaTotFlux"]


class DiaPipelineTask(pipeBase.PipelineTask):
    """Task for loading, associating and storing Difference Image Analysis
    (DIA) Objects and Sources.
    """
    ConfigClass = DiaPipelineConfig
    _DefaultName = "diaPipe"

    def __init__(self, initInputs=None, **kwargs):
        super().__init__(**kwargs)
        self.apdb = self.config.apdb.apply()
        self.makeSubtask("diaCatalogLoader")
        self.makeSubtask("associator")
        self.makeSubtask("diaCalculation")
        self.makeSubtask("diaForcedSource")
        if self.config.doPackageAlerts:
            self.makeSubtask("alertPackager")
        if self.config.doSolarSystemAssociation:
            self.makeSubtask("solarSystemAssociator")

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        inputs["idGenerator"] = self.config.idGenerator.apply(butlerQC.quantum.dataId)
        inputs["band"] = butlerQC.quantum.dataId["band"]
        if not self.config.doSolarSystemAssociation:
            inputs["solarSystemObjectTable"] = None

        outputs = self.run(**inputs)

        butlerQC.put(outputs, outputRefs)

    @timeMethod
    def run(self,
            diaSourceTable,
            solarSystemObjectTable,
            diffIm,
            exposure,
            template,
            band,
            idGenerator):
        """Process DiaSources and DiaObjects.

        Load previous DiaObjects and their DiaSource history. Calibrate the
        values in the diaSourceCat. Associate new DiaSources with previous
        DiaObjects. Run forced photometry at the updated DiaObject locations.
        Store the results in the Alert Production Database (Apdb).

        Parameters
        ----------
        diaSourceTable : `pandas.DataFrame`
            Newly detected DiaSources.
        diffIm : `lsst.afw.image.ExposureF`
            Difference image exposure in which the sources in ``diaSourceCat``
            were detected.
        exposure : `lsst.afw.image.ExposureF`
            Calibrated exposure differenced with a template to create
            ``diffIm``.
        template : `lsst.afw.image.ExposureF`
            Template exposure used to create diffIm.
        band : `str`
            The band in which the new DiaSources were detected.
        idGenerator : `lsst.meas.base.IdGenerator`
            Object that generates source IDs and random number generator seeds.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            Results struct with components.

            - ``apdbMaker`` : Marker dataset to store in the Butler indicating
              that this ccdVisit has completed successfully.
              (`lsst.dax.apdb.ApdbConfig`)
            - ``associatedDiaSources`` : Catalog of newly associated
              DiaSources. (`pandas.DataFrame`)
        """
        # Load the DiaObjects and DiaSource history.
        loaderResult = self.diaCatalogLoader.run(diffIm, self.apdb)
        if len(loaderResult.diaObjects) > 0:
            diaObjects = self.purgeDiaObjects(diffIm.getBBox(), diffIm.getWcs(), loaderResult.diaObjects,
                                              buffer=self.config.imagePixelMargin)
        else:
            diaObjects = loaderResult.diaObjects

        # Associate new DiaSources with existing DiaObjects.
        assocResults = self.associator.run(diaSourceTable, diaObjects,
                                           exposure_time=diffIm.visitInfo.exposureTime)

        if self.config.doSolarSystemAssociation:
            ssoAssocResult = self.solarSystemAssociator.run(
                assocResults.unAssocDiaSources,
                solarSystemObjectTable,
                diffIm)
            createResults = self.createNewDiaObjects(
                ssoAssocResult.unAssocDiaSources)
            toAssociate = []
            if len(assocResults.matchedDiaSources) > 0:
                toAssociate.append(assocResults.matchedDiaSources)
            if len(ssoAssocResult.ssoAssocDiaSources) > 0:
                toAssociate.append(ssoAssocResult.ssoAssocDiaSources)
            toAssociate.append(createResults.diaSources)
            associatedDiaSources = pd.concat(toAssociate)
            nTotalSsObjects = ssoAssocResult.nTotalSsObjects
            nAssociatedSsObjects = ssoAssocResult.nAssociatedSsObjects
        else:
            createResults = self.createNewDiaObjects(
                assocResults.unAssocDiaSources)
            toAssociate = []
            if len(assocResults.matchedDiaSources) > 0:
                toAssociate.append(assocResults.matchedDiaSources)
            toAssociate.append(createResults.diaSources)
            associatedDiaSources = pd.concat(toAssociate)
            nTotalSsObjects = 0
            nAssociatedSsObjects = 0

        # Create new DiaObjects from unassociated diaSources.
        self._add_association_meta_data(assocResults.nUpdatedDiaObjects,
                                        assocResults.nUnassociatedDiaObjects,
                                        createResults.nNewDiaObjects,
                                        nTotalSsObjects,
                                        nAssociatedSsObjects)
        # Index the DiaSource catalog for this visit after all associations
        # have been made.
        updatedDiaObjectIds = associatedDiaSources["diaObjectId"][
            associatedDiaSources["diaObjectId"] != 0].to_numpy()
        associatedDiaSources.set_index(["diaObjectId",
                                        "band",
                                        "diaSourceId"],
                                       drop=False,
                                       inplace=True)

        # Append new DiaObjects and DiaSources to their previous history.
        diaObjects = pd.concat(
            [diaObjects,
             createResults.newDiaObjects.set_index("diaObjectId", drop=False)],
            sort=True)
        if self.testDataFrameIndex(diaObjects):
            raise RuntimeError(
                "Duplicate DiaObjects created after association. This is "
                "likely due to re-running data with an already populated "
                "Apdb. If this was not the case then there was an unexpected "
                "failure in Association while matching and creating new "
                "DiaObjects and should be reported. Exiting.")

        if len(loaderResult.diaSources) > 0:
            # We need to coerce the types of loaderResult.diaSources
            # to be the same as associatedDiaSources, thanks to pandas
            # datetime issues (DM-41100). And we may as well coerce
            # all the columns to ensure consistency for future compatibility.
            for name, dtype in associatedDiaSources.dtypes.items():
                if name in loaderResult.diaSources.columns and loaderResult.diaSources[name].dtype != dtype:
                    self.log.debug(
                        "Coercing loaderResult.diaSources column %s from %s to %s",
                        name,
                        str(loaderResult.diaSources[name].dtype),
                        str(dtype),
                    )
                    loaderResult.diaSources[name] = loaderResult.diaSources[name].astype(dtype)

            mergedDiaSourceHistory = pd.concat(
                [loaderResult.diaSources, associatedDiaSources],
                sort=True)
        else:
            mergedDiaSourceHistory = pd.concat([associatedDiaSources], sort=True)

        # Test for DiaSource duplication first. If duplicates are found,
        # this likely means this is duplicate data being processed and sent
        # to the Apdb.
        if self.testDataFrameIndex(mergedDiaSourceHistory):
            raise RuntimeError(
                "Duplicate DiaSources found after association and merging "
                "with history. This is likely due to re-running data with an "
                "already populated Apdb. If this was not the case then there "
                "was an unexpected failure in Association while matching "
                "sources to objects, and should be reported. Exiting.")

        # Compute DiaObject Summary statistics from their full DiaSource
        # history.
        diaCalResult = self.diaCalculation.run(
            diaObjects,
            mergedDiaSourceHistory,
            updatedDiaObjectIds,
            [band])
        # Test for duplication in the updated DiaObjects.
        if self.testDataFrameIndex(diaCalResult.diaObjectCat):
            raise RuntimeError(
                "Duplicate DiaObjects (loaded + updated) created after "
                "DiaCalculation. This is unexpected behavior and should be "
                "reported. Exiting.")
        if self.testDataFrameIndex(diaCalResult.updatedDiaObjects):
            raise RuntimeError(
                "Duplicate DiaObjects (updated) created after "
                "DiaCalculation. This is unexpected behavior and should be "
                "reported. Exiting.")

        # Force photometer on the Difference and Calibrated exposures using
        # the new and updated DiaObject locations.
        diaForcedSources = self.diaForcedSource.run(
            diaCalResult.diaObjectCat,
            diaCalResult.updatedDiaObjects.loc[:, "diaObjectId"].to_numpy(),
            exposure,
            diffIm,
            idGenerator=idGenerator)

        # Store DiaSources, updated DiaObjects, and DiaForcedSources in the
        # Apdb.
        self.apdb.store(
            DateTime.now(),
            diaCalResult.updatedDiaObjects,
            associatedDiaSources,
            diaForcedSources)

        if self.config.doPackageAlerts:
            if len(loaderResult.diaForcedSources) > 1:
                # We need to coerce the types of loaderResult.diaForcedSources
                # to be the same as associatedDiaSources, thanks to pandas
                # datetime issues (DM-41100). And we may as well coerce
                # all the columns to ensure consistency for future compatibility.
                for name, dtype in diaForcedSources.dtypes.items():
                    if name in loaderResult.diaForcedSources.columns and \
                       loaderResult.diaForcedSources[name].dtype != dtype:
                        self.log.debug(
                            "Coercing loaderResult.diaForcedSources column %s from %s to %s",
                            name,
                            str(loaderResult.diaForcedSources[name].dtype),
                            str(dtype),
                        )
                        loaderResult.diaForcedSources[name] = (
                            loaderResult.diaForcedSources[name].astype(dtype)
                        )
                diaForcedSources = pd.concat(
                    [diaForcedSources, loaderResult.diaForcedSources],
                    sort=True)
            if self.testDataFrameIndex(diaForcedSources):
                self.log.warning(
                    "Duplicate DiaForcedSources created after merge with "
                    "history and new sources. This may cause downstream "
                    "problems. Dropping duplicates.")
                # Drop duplicates via index and keep the first appearance.
                # Reset due to the index shape being slight different than
                # expected.
                diaForcedSources = diaForcedSources.groupby(
                    diaForcedSources.index).first()
                diaForcedSources.reset_index(drop=True, inplace=True)
                diaForcedSources.set_index(
                    ["diaObjectId", "diaForcedSourceId"],
                    drop=False,
                    inplace=True)
            self.alertPackager.run(associatedDiaSources,
                                   diaCalResult.diaObjectCat,
                                   loaderResult.diaSources,
                                   diaForcedSources,
                                   diffIm,
                                   template)

        return pipeBase.Struct(apdbMarker=self.config.apdb.value,
                               associatedDiaSources=associatedDiaSources,
                               diaForcedSources=diaForcedSources,
                               diaObjects=diaObjects,
                               longTrailedSources=assocResults.longTrailedSources
                               )

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

            - ``diaSources`` : DiaSource catalog with updated DiaObject ids.
              (`pandas.DataFrame`)
            - ``newDiaObjects`` : Newly created DiaObjects from the
              unassociated DiaSources. (`pandas.DataFrame`)
            - ``nNewDiaObjects`` : Number of newly created diaObjects.(`int`)
        """
        if len(unAssocDiaSources) == 0:
            tmpObj = self._initialize_dia_object(0)
            newDiaObjects = pd.DataFrame(data=[],
                                         columns=tmpObj.keys())
        else:
            newDiaObjects = unAssocDiaSources["diaSourceId"].apply(
                self._initialize_dia_object)
            unAssocDiaSources["diaObjectId"] = unAssocDiaSources["diaSourceId"]
        return pipeBase.Struct(diaSources=unAssocDiaSources,
                               newDiaObjects=newDiaObjects,
                               nNewDiaObjects=len(newDiaObjects))

    def _initialize_dia_object(self, objId):
        """Create a new DiaObject with values required to be initialized by the
        Ppdb.

        Parameters
        ----------
        objid : `int`
            ``diaObjectId`` value for the of the new DiaObject.

        Returns
        -------
        diaObject : `dict`
            Newly created DiaObject with keys:

            ``diaObjectId``
                Unique DiaObjectId (`int`).
            ``pmParallaxNdata``
                Number of data points used for parallax calculation (`int`).
            ``nearbyObj1``
                Id of the a nearbyObject in the Object table (`int`).
            ``nearbyObj2``
                Id of the a nearbyObject in the Object table (`int`).
            ``nearbyObj3``
                Id of the a nearbyObject in the Object table (`int`).
            ``?_psfFluxNdata``
                Number of data points used to calculate point source flux
                summary statistics in each bandpass (`int`).
        """
        new_dia_object = {"diaObjectId": objId,
                          "pmParallaxNdata": 0,
                          "nearbyObj1": 0,
                          "nearbyObj2": 0,
                          "nearbyObj3": 0,
                          "flags": 0}
        for f in ["u", "g", "r", "i", "z", "y"]:
            new_dia_object["%s_psfFluxNdata" % f] = 0
        return pd.Series(data=new_dia_object)

    def testDataFrameIndex(self, df):
        """Test the sorted DataFrame index for duplicates.

        Wrapped as a separate function to allow for mocking of the this task
        in unittesting. Default of a mock return for this test is True.

        Parameters
        ----------
        df : `pandas.DataFrame`
            DataFrame to text.

        Returns
        -------
        `bool`
            True if DataFrame contains duplicate rows.
        """
        return df.index.has_duplicates

    def _add_association_meta_data(self,
                                   nUpdatedDiaObjects,
                                   nUnassociatedDiaObjects,
                                   nNewDiaObjects,
                                   nTotalSsObjects,
                                   nAssociatedSsObjects):
        """Store summaries of the association step in the task metadata.

        Parameters
        ----------
        nUpdatedDiaObjects : `int`
            Number of previous DiaObjects associated and updated in this
            ccdVisit.
        nUnassociatedDiaObjects : `int`
            Number of previous DiaObjects that were not associated or updated
            in this ccdVisit.
        nNewDiaObjects : `int`
            Number of newly created DiaObjects for this ccdVisit.
        nTotalSsObjects : `int`
            Number of SolarSystemObjects within the observable detector
            area.
        nAssociatedSsObjects : `int`
            Number of successfully associated SolarSystemObjects.
        """
        self.metadata.add('numUpdatedDiaObjects', nUpdatedDiaObjects)
        self.metadata.add('numUnassociatedDiaObjects', nUnassociatedDiaObjects)
        self.metadata.add('numNewDiaObjects', nNewDiaObjects)
        self.metadata.add('numTotalSolarSystemObjects', nTotalSsObjects)
        self.metadata.add('numAssociatedSsObjects', nAssociatedSsObjects)

    def purgeDiaObjects(self, bbox, wcs, diaObjCat, buffer=0):
        """Drop diaObjects that are outside the exposure bounding box.

        Parameters
        ----------
        bbox : `lsst.geom.Box2I`
            Bounding box of the exposure.
        wcs : `lsst.afw.geom.SkyWcs`
            Coordinate system definition (wcs) for the exposure.
        diaObjCat : `pandas.DataFrame`
            DiaObjects loaded from the Apdb.
        buffer : `int`, optional
            Width, in pixels, to pad the exposure bounding box.

        Returns
        -------
        diaObjCat : `pandas.DataFrame`
            DiaObjects loaded from the Apdb, restricted to the exposure
            bounding box.
        """
        try:
            bbox.grow(buffer)
            raVals = diaObjCat.ra.to_numpy()
            decVals = diaObjCat.dec.to_numpy()
            xVals, yVals = wcs.skyToPixelArray(raVals, decVals, degrees=True)
            selector = bbox.contains(xVals, yVals)
            nPurged = np.sum(~selector)
            if nPurged > 0:
                diaObjCat = diaObjCat[selector].copy()
                self.log.info(f"Dropped {nPurged} diaObjects that were outside the bbox "
                              f"leaving {len(diaObjCat)} in the catalog")
        except Exception as e:
            self.log.warning("Error attempting to check diaObject history: %s", e)
        return diaObjCat
