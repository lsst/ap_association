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
"""

__all__ = ("DiaPipelineConfig",
           "DiaPipelineTask",
           "DiaPipelineConnections")


import warnings

import lsst.dax.apdb as daxApdb
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as connTypes
import numpy as np
import pandas as pd
from lsst.ap.association import (
    AssociationTask,
    DiaForcedSourceTask,
    PackageAlertsTask)
from lsst.ap.association.ssoAssociation import SolarSystemAssociationTask
from lsst.ap.association.utils import convertTableToSdmSchema, readSchemaFromApdb, dropEmptyColumns, \
    make_empty_catalog, makeEmptyForcedSourceTable
from lsst.daf.base import DateTime
from lsst.meas.base import DetectorVisitIdGeneratorConfig, \
    DiaObjectCalculationTask
from lsst.utils.timer import timeMethod


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
        name="preloaded_SsObjects",
        storageClass="DataFrame",
        dimensions=("instrument", "group", "detector"),
        minimum=0,
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
    preloadedDiaObjects = connTypes.Input(
        doc="DiaObjects preloaded from the APDB.",
        name="preloaded_diaObjects",
        storageClass="DataFrame",
        dimensions=("instrument", "group", "detector"),
    )
    preloadedDiaSources = connTypes.Input(
        doc="DiaSources preloaded from the APDB.",
        name="preloaded_diaSources",
        storageClass="DataFrame",
        dimensions=("instrument", "group", "detector"),
    )
    preloadedDiaForcedSources = connTypes.Input(
        doc="DiaForcedSources preloaded from the APDB.",
        name="preloaded_diaForcedSources",
        storageClass="DataFrame",
        dimensions=("instrument", "group", "detector"),
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

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not config.doWriteAssociatedSources:
            self.outputs.remove("associatedDiaSources")
            self.outputs.remove("diaForcedSources")
            self.outputs.remove("diaObjects")
        elif not config.doRunForcedMeasurement:
            self.outputs.remove("diaForcedSources")
        if not config.doSolarSystemAssociation:
            self.inputs.remove("solarSystemObjectTable")

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
    apdb = pexConfig.ConfigurableField(  # TODO: remove on DM-43419
        target=daxApdb.ApdbSql,
        doc="Database connection for storing associated DiaSources and "
            "DiaObjects. Must already be initialized.",
        deprecated="This field has been replaced by ``apdb_config_url``; set "
                   "``doConfigureApdb=False`` to use it. Will be removed after v28.",
    )
    apdb_config_url = pexConfig.Field(
        dtype=str,
        default=None,
        optional=False,
        doc="A config file specifying the APDB and its connection parameters, "
            "typically written by the apdb-cli command-line utility. "
            "The database must already be initialized.",
    )
    validBands = pexConfig.ListField(
        dtype=str,
        default=["u", "g", "r", "i", "z", "y"],
        doc="List of bands that are valid for AP processing. To process a "
            "band not on this list, the appropriate band specific columns "
            "must be added to the Apdb schema in dax_apdb.",
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
    doRunForcedMeasurement = pexConfig.Field(
        dtype=bool,
        default=True,
        deprecated="Added to allow disabling forced sources for performance "
                   "reasons during the ops rehearsal. "
                   "It is expected to be removed.",
        doc="Run forced measurement on all of the diaObjects? "
            "This should only be turned off for debugging purposes.",
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
    doConfigureApdb = pexConfig.Field(  # TODO: remove on DM-43419
        dtype=bool,
        default=True,
        doc="Use the deprecated ``apdb`` sub-config to set up the APDB, "
            "instead of the new config (``apdb_config_url``). This field is "
            "provided for backward-compatibility ONLY and will be removed "
            "without notice after v28.",
    )

    def setDefaults(self):
        if self.doConfigureApdb:
            self.apdb.dia_object_index = "baseline"
            self.apdb.dia_object_columns = []
        self.diaCalculation.plugins = ["ap_meanPosition",
                                       "ap_nDiaSources",
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

    # TODO: remove on DM-43419
    def validate(self):
        # Sidestep Config.validate to avoid validating uninitialized fields we're not using.
        skip = {"apdb_config_url"} if self.doConfigureApdb else {"apdb"}
        for name, field in self._fields.items():
            if name not in skip:
                field.validate(self)

        # It's possible to use apdb without setting it, bypassing the deprecation warning.
        if self.doConfigureApdb:
            warnings.warn("Config field DiaPipelineConfig.apdb is deprecated: "
                          # Workaround for DM-44051
                          "This field has been replaced by ``apdb_config_url``; set "
                          "``doConfigureApdb=False`` to use it. Will be removed after v28.",
                          FutureWarning)


class DiaPipelineTask(pipeBase.PipelineTask):
    """Task for loading, associating and storing Difference Image Analysis
    (DIA) Objects and Sources.
    """
    ConfigClass = DiaPipelineConfig
    _DefaultName = "diaPipe"

    def __init__(self, initInputs=None, **kwargs):
        super().__init__(**kwargs)
        if self.config.doConfigureApdb:
            self.apdb = self.config.apdb.apply()
        else:
            self.apdb = daxApdb.Apdb.from_uri(self.config.apdb_config_url)
        self.schema = readSchemaFromApdb(self.apdb)
        self.makeSubtask("associator")
        self.makeSubtask("diaCalculation")
        if self.config.doRunForcedMeasurement:
            self.makeSubtask("diaForcedSource")
        if self.config.doPackageAlerts:
            self.makeSubtask("alertPackager")
        if self.config.doSolarSystemAssociation:
            self.makeSubtask("solarSystemAssociator")

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        inputs["idGenerator"] = self.config.idGenerator.apply(butlerQC.quantum.dataId)
        inputs["band"] = butlerQC.quantum.dataId["band"]
        inputs["legacySolarSystemTable"] = None
        if not self.config.doSolarSystemAssociation:
            inputs["solarSystemObjectTable"] = None

        outputs = self.run(**inputs)

        butlerQC.put(outputs, outputRefs)

    @timeMethod
    def run(self,
            diaSourceTable,
            legacySolarSystemTable,
            diffIm,
            exposure,
            template,
            preloadedDiaObjects,
            preloadedDiaSources,
            preloadedDiaForcedSources,
            band,
            idGenerator,
            solarSystemObjectTable=None):
        """Process DiaSources and DiaObjects.

        Load previous DiaObjects and their DiaSource history. Calibrate the
        values in the diaSourceCat. Associate new DiaSources with previous
        DiaObjects. Run forced photometry at the updated DiaObject locations.
        Store the results in the Alert Production Database (Apdb).

        Parameters
        ----------
        diaSourceTable : `pandas.DataFrame`
            Newly detected DiaSources.
        legacySolarSystemTable : `pandas.DataFrame`
            Not used
        diffIm : `lsst.afw.image.ExposureF`
            Difference image exposure in which the sources in ``diaSourceCat``
            were detected.
        exposure : `lsst.afw.image.ExposureF`
            Calibrated exposure differenced with a template to create
            ``diffIm``.
        template : `lsst.afw.image.ExposureF`
            Template exposure used to create diffIm.
        preloadedDiaObjects : `pandas.DataFrame`
            Previously detected DiaObjects, loaded from the APDB.
        preloadedDiaSources : `pandas.DataFrame`
            Previously detected DiaSources, loaded from the APDB.
        preloadedDiaForcedSources : `pandas.DataFrame`
            Catalog of previously detected forced DiaSources, from the APDB
        band : `str`
            The band in which the new DiaSources were detected.
        idGenerator : `lsst.meas.base.IdGenerator`
            Object that generates source IDs and random number generator seeds.
        solarSystemObjectTable : `pandas.DataFrame`
            Preloaded Solar System objects expected to be visible in the image.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            Results struct with components.

            - ``apdbMarker`` : Marker dataset to store in the Butler indicating
              that this ccdVisit has completed successfully.
              (`lsst.dax.apdb.ApdbConfig`)
            - ``associatedDiaSources`` : Catalog of newly associated
              DiaSources. (`pandas.DataFrame`)
            - ``diaForcedSources`` : Catalog of new and previously detected
              forced DiaSources. (`pandas.DataFrame`)
            - ``diaObjects`` : Updated table of DiaObjects. (`pandas.DataFrame`)

        Raises
        ------
        RuntimeError
            Raised if duplicate DiaObjects or duplicate DiaSources are found.
        """
        # Accept either legacySolarSystemTable or optional solarSystemObjectTable.
        if legacySolarSystemTable is not None and solarSystemObjectTable is None:
            solarSystemObjectTable = legacySolarSystemTable

        if not preloadedDiaObjects.empty:
            # Include a small buffer outside the image so that we can associate sources near the edge
            diaObjects, _ = self.purgeDiaObjects(diffIm.getBBox(), diffIm.getWcs(), preloadedDiaObjects,
                                                 buffer=self.config.imagePixelMargin)
        else:
            diaObjects = preloadedDiaObjects

        # Associate DiaSources with DiaObjects
        associatedDiaSources, newDiaObjects = self.associateDiaSources(
            diaSourceTable, solarSystemObjectTable, diffIm, diaObjects
        )

        # Merge associated diaSources
        mergedDiaSourceHistory, mergedDiaObjects, updatedDiaObjectIds = self.mergeAssociatedCatalogs(
            preloadedDiaSources, associatedDiaSources, diaObjects, newDiaObjects, diffIm
        )

        # Compute DiaObject Summary statistics from their full DiaSource
        # history.
        diaCalResult = self.diaCalculation.run(
            mergedDiaObjects,
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

        # Forced source measurement
        if self.config.doRunForcedMeasurement:
            diaForcedSources = self.runForcedMeasurement(
                diaCalResult.diaObjectCat, diaCalResult.updatedDiaObjects, exposure, diffIm, idGenerator
            )

        else:
            # alertPackager needs correct columns
            diaForcedSources = makeEmptyForcedSourceTable(self.schema)

        # Write results to Alert Production Database (APDB)
        self.writeToApdb(diaCalResult.updatedDiaObjects, associatedDiaSources, diaForcedSources)

        # Package alerts
        if self.config.doPackageAlerts:
            # Append new forced sources to the full history
            diaForcedSourcesFull = self.mergeCatalogs(preloadedDiaForcedSources, diaForcedSources,
                                                      "preloadedDiaForcedSources")
            if self.testDataFrameIndex(diaForcedSources):
                self.log.warning(
                    "Duplicate DiaForcedSources created after merge with "
                    "history and new sources. This may cause downstream "
                    "problems. Dropping duplicates.")
                # Drop duplicates via index and keep the first appearance.
                # Reset due to the index shape being slight different than
                # expected.
                diaForcedSourcesFull = diaForcedSourcesFull.groupby(
                    diaForcedSourcesFull.index).first()
                diaForcedSourcesFull.reset_index(drop=True, inplace=True)
                diaForcedSourcesFull.set_index(
                    ["diaObjectId", "diaForcedSourceId"],
                    drop=False,
                    inplace=True)
            self.alertPackager.run(associatedDiaSources,
                                   diaCalResult.diaObjectCat,
                                   preloadedDiaSources,
                                   diaForcedSourcesFull,
                                   diffIm,
                                   exposure,
                                   template,
                                   doRunForcedMeasurement=self.config.doRunForcedMeasurement,
                                   )

        # For historical reasons, apdbMarker is a Config even if it's not meant to be read.
        # A default Config is the cheapest way to satisfy the storage class.
        marker = self.config.apdb.value if self.config.doConfigureApdb else pexConfig.Config()
        return pipeBase.Struct(apdbMarker=marker,
                               associatedDiaSources=associatedDiaSources,
                               diaForcedSources=diaForcedSources,
                               diaObjects=diaCalResult.diaObjectCat,
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
            unAssocDiaSources["diaObjectId"] = unAssocDiaSources["diaSourceId"]
            newDiaObjects = convertTableToSdmSchema(self.schema, unAssocDiaSources,
                                                    tableName="DiaObject")
        return pipeBase.Struct(diaSources=unAssocDiaSources,
                               newDiaObjects=newDiaObjects,
                               nNewDiaObjects=len(newDiaObjects))

    @timeMethod
    def associateDiaSources(self, diaSourceTable, solarSystemObjectTable, diffIm, diaObjects):
        """Associate DiaSources with DiaObjects.

        Associate new DiaSources with existing DiaObjects. Create new
        DiaObjects fron unassociated DiaSources. Index DiaSource catalogue
        after associations. Append new DiaObjects and DiaSources to their
        previous history. Test for DiaSource and DiaObject duplications.
        Compute DiaObject Summary statistics from their full DiaSource
        history. Test for duplication in the updated DiaObjects.

        Parameters
        ----------
        diaSourceTable : `pandas.DataFrame`
            Newly detected DiaSources.
        solarSystemObjectTable : `pandas.DataFrame`
            Preloaded Solar System objects expected to be visible in the image.
        diffIm : `lsst.afw.image.ExposureF`
            Difference image exposure in which the sources in ``diaSourceCat``
            were detected.
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
        assocResults = self.associator.run(diaSourceTable, diaObjects)

        if self.config.doSolarSystemAssociation and solarSystemObjectTable is not None:
            ssoAssocResult = self.solarSystemAssociator.run(
                assocResults.unAssocDiaSources,
                solarSystemObjectTable,
                diffIm)
            # Create new DiaObjects from unassociated diaSources.
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
            # Create new DiaObjects from unassociated diaSources.
            createResults = self.createNewDiaObjects(
                assocResults.unAssocDiaSources)
            toAssociate = []
            if len(assocResults.matchedDiaSources) > 0:
                toAssociate.append(assocResults.matchedDiaSources)
            toAssociate.append(createResults.diaSources)
            associatedDiaSources = pd.concat(toAssociate)
            nTotalSsObjects = 0
            nAssociatedSsObjects = 0

        self._add_association_meta_data(assocResults.nUpdatedDiaObjects,
                                        assocResults.nUnassociatedDiaObjects,
                                        createResults.nNewDiaObjects,
                                        nTotalSsObjects,
                                        nAssociatedSsObjects)
        self.log.info("%i updated and %i unassociated diaObjects. Creating %i new diaObjects",
                      assocResults.nUpdatedDiaObjects,
                      assocResults.nUnassociatedDiaObjects,
                      createResults.nNewDiaObjects,
                      )
        return (associatedDiaSources, createResults.newDiaObjects)

    @timeMethod
    def mergeAssociatedCatalogs(self, preloadedDiaSources, associatedDiaSources, diaObjects, newDiaObjects,
                                diffIm):
        """Merge the associated diaSource and diaObjects to their previous history.

        Parameters
        ----------
        preloadedDiaSources : `pandas.DataFrame`
            Previously detected DiaSources, loaded from the APDB.
        associatedDiaSources : `pandas.DataFrame`
            Associated DiaSources with DiaObjects.
        diaObjects : `pandas.DataFrame`
            Table of  DiaObjects from preloaded DiaObjects.
        newDiaObjects : `pandas.DataFrame`
            Table of new DiaObjects after association.

        Raises
        ------
        RuntimeError
            Raised if duplicate DiaObjects or duplicate DiaSources are found.

        Returns
        -------
        mergedDiaSourceHistory : `pandas.DataFrame`
            The combined catalog, with all of the rows from preloadedDiaSources
            catalog ordered before the rows of associatedDiaSources catalog.
        mergedDiaObjects : `pandas.DataFrame`
            Table of new DiaObjects merged with their history.
        updatedDiaObjectIds : `numpy.Array`
            Object Id's from associated diaSources.
        """
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
        if diaObjects.empty:
            mergedDiaObjects = newDiaObjects.set_index("diaObjectId", drop=False)
        elif not newDiaObjects.empty:
            mergedDiaObjects = pd.concat(
                [diaObjects,
                 newDiaObjects.set_index("diaObjectId", drop=False)],
                sort=True)
        else:
            mergedDiaObjects = diaObjects

        # Exclude any objects that are off the image after association.
        mergedDiaObjects, updatedDiaObjectIds = self.purgeDiaObjects(diffIm.getBBox(), diffIm.getWcs(),
                                                                     mergedDiaObjects,
                                                                     diaObjectIds=updatedDiaObjectIds,
                                                                     buffer=-1)
        if self.testDataFrameIndex(mergedDiaObjects):
            raise RuntimeError(
                "Duplicate DiaObjects created after association. This is "
                "likely due to re-running data with an already populated "
                "Apdb. If this was not the case then there was an unexpected "
                "failure in Association while matching and creating new "
                "DiaObjects and should be reported. Exiting.")

        mergedDiaSourceHistory = self.mergeCatalogs(preloadedDiaSources, associatedDiaSources,
                                                    "preloadedDiaSources")

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
        return (mergedDiaSourceHistory, mergedDiaObjects, updatedDiaObjectIds)

    @timeMethod
    def runForcedMeasurement(self, diaObjects, updatedDiaObjects, exposure, diffIm, idGenerator):
        """Forced Source Measurement

        Forced photometry on the difference and calibrated exposures using the
        new and updated DiaObject locations.

        Parameters
        ----------
        diaObjects : `pandas.DataFrame`
            Catalog of DiaObjects.
        updatedDiaObjects : `pandas.DataFrame`
            Catalog of updated DiaObjects.
        exposure : `lsst.afw.image.ExposureF`
            Calibrated exposure differenced with a template to create
            ``diffIm``.
        diffIm : `lsst.afw.image.ExposureF`
            Difference image exposure in which the sources in ``diaSourceCat``
            were detected.
        idGenerator : `lsst.meas.base.IdGenerator`
            Object that generates source IDs and random number generator seeds.

        Returns
        -------
        diaForcedSources : `pandas.DataFrame`
            Catalog of calibrated forced photometered fluxes on both the
            difference and direct images at DiaObject locations.
        """
        # Force photometer on the Difference and Calibrated exposures using
        # the new and updated DiaObject locations.
        diaForcedSources = self.diaForcedSource.run(
            diaObjects,
            updatedDiaObjects.loc[:, "diaObjectId"].to_numpy(),
            exposure,
            diffIm,
            idGenerator=idGenerator)
        self.log.info(f"Updating {len(diaForcedSources)} diaForcedSources in the APDB")
        diaForcedSources = convertTableToSdmSchema(self.schema, diaForcedSources,
                                                   tableName="DiaForcedSource",
                                                   )
        return diaForcedSources

    @timeMethod
    def writeToApdb(self, updatedDiaObjects, associatedDiaSources, diaForcedSources):
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
        diaSourceStore = dropEmptyColumns(self.schema, associatedDiaSources, tableName="DiaSource")
        diaForcedSourceStore = dropEmptyColumns(self.schema, diaForcedSources, tableName="DiaForcedSource")
        self.apdb.store(
            DateTime.now().toAstropy(),
            diaObjectStore,
            diaSourceStore,
            diaForcedSourceStore)
        self.log.info("APDB updated.")

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

    def purgeDiaObjects(self, bbox, wcs, diaObjCat, diaObjectIds=None, buffer=0):
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
                if diaObjectIds is not None:
                    # We also need to drop any of the associated IDs if this runs after association
                    purgedIds = diaObjCat[~selector].diaObjectId
                    diaObjectIds = diaObjectIds[~np.isin(diaObjectIds, purgedIds)]
                    self.log.info("Dropped %i diaObjects that were outside the bbox "
                                  "after association, leaving %i in the catalog",
                                  nPurged, len(diaObjCat) - nPurged)
                else:
                    self.log.info("Dropped %i diaObjects that were outside the padded bbox "
                                  "before association, leaving %i in the catalog",
                                  nPurged, len(diaObjCat) - nPurged)
                diaObjCat = diaObjCat[selector].copy()
        except Exception as e:
            self.log.warning("Error attempting to check diaObject history: %s", e)
        return diaObjCat, diaObjectIds

    def mergeCatalogs(self, originalCatalog, newCatalog, catalogName):
        """Combine two catalogs, ensuring that the columns of the new catalog
        have the same dtype as the original.

        Parameters
        ----------
        originalCatalog : `pandas.DataFrame`
            The original catalog to be added to.
        newCatalog : `pandas.DataFrame`
            The new catalog to append to `originalCatalog`
        catalogName : `str`, optional
            The name of the catalog to use for logging messages.

        Returns
        -------
        mergedCatalog : `pandas.DataFrame`
            The combined catalog, with all of the rows from ``originalCatalog``
            ordered before the rows of ``newCatalog``
        """
        if len(newCatalog) > 0:
            catalog = newCatalog.copy(deep=True)
            # We need to coerce the types of `newCatalog`
            # to be the same as `originalCatalog`, thanks to pandas
            # datetime issues (DM-41100). And we may as well coerce
            # all the columns to ensure consistency for future compatibility.
            for name, dtype in originalCatalog.dtypes.items():
                if name in newCatalog.columns and newCatalog[name].dtype != dtype:
                    self.log.debug(
                        "Coercing %s column %s from %s to %s",
                        catalogName,
                        name,
                        str(newCatalog[name].dtype),
                        str(dtype),
                    )
                    catalog[name] = newCatalog[name].astype(dtype)

            mergedCatalog = pd.concat([originalCatalog, catalog], sort=True)
        else:
            mergedCatalog = pd.concat([originalCatalog], sort=True)
        return mergedCatalog.loc[:, originalCatalog.columns]
