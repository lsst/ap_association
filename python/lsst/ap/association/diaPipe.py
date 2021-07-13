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

import os

import lsst.dax.apdb as daxApdb
from lsst.meas.base import DiaObjectCalculationTask
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as connTypes

from lsst.ap.association import (
    AssociationTask,
    DiaForcedSourceTask,
    LoadDiaCatalogsTask,
    make_dia_object_schema,
    make_dia_source_schema,
    PackageAlertsTask)

__all__ = ("DiaPipelineConfig",
           "DiaPipelineTask",
           "DiaPipelineConnections")


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
    diffIm = connTypes.Input(
        doc="Difference image on which the DiaSources were detected.",
        name="{fakesType}{coaddName}Diff_differenceExp",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector"),
    )
    exposure = connTypes.Input(
        doc="Calibrated exposure differenced with a template image during "
            "image differencing.",
        name="calexp",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector"),
    )
    warpedExposure = connTypes.Input(
        doc="Warped template used to create `subtractedExposure`. Not PSF "
            "matched.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="ExposureF",
        name="{fakesType}{coaddName}Diff_warpedExp",
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
            "calibration, and standardization for insertation into the Apdb.",
        name="{fakesType}{coaddName}Diff_assocDiaSrc",
        storageClass="DataFrame",
        dimensions=("instrument", "visit", "detector"),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not config.doWriteAssociatedSources:
            self.outputs.remove("associatedDiaSources")

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
    apdb = pexConfig.ConfigurableField(
        target=daxApdb.Apdb,
        ConfigClass=daxApdb.ApdbConfig,
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
        default=False,
        doc="Write out associated and SDMed DiaSources.",
    )

    def setDefaults(self):
        self.apdb.dia_object_index = "baseline"
        self.apdb.dia_object_columns = []
        self.apdb.extra_schema_file = os.path.join(
            "${AP_ASSOCIATION_DIR}",
            "data",
            "apdb-ap-pipe-schema-extra.yaml")
        self.diaCalculation.plugins = ["ap_meanPosition",
                                       "ap_HTMIndex",
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

    def validate(self):
        pexConfig.Config.validate(self)
        if self.diaCatalogLoader.htmLevel != \
           self.diaCalculation.plugins["ap_HTMIndex"].htmLevel:
            raise ValueError("HTM index level in LoadDiaCatalogsTask must be "
                             "equal to HTMIndexDiaCalculationPlugin index "
                             "level.")
        if "ap_HTMIndex" not in self.diaCalculation.plugins:
            raise ValueError("DiaPipe requires the ap_HTMIndex plugin "
                             "be enabled for proper insertion into the Apdb.")


class DiaPipelineTask(pipeBase.PipelineTask):
    """Task for loading, associating and storing Difference Image Analysis
    (DIA) Objects and Sources.
    """
    ConfigClass = DiaPipelineConfig
    _DefaultName = "diaPipe"
    RunnerClass = pipeBase.ButlerInitializedTaskRunner

    def __init__(self, initInputs=None, **kwargs):
        super().__init__(**kwargs)
        self.apdb = self.config.apdb.apply(
            afw_schemas=dict(DiaObject=make_dia_object_schema(),
                             DiaSource=make_dia_source_schema()))
        self.makeSubtask("diaCatalogLoader")
        self.makeSubtask("associator")
        self.makeSubtask("diaCalculation")
        self.makeSubtask("diaForcedSource")
        if self.config.doPackageAlerts:
            self.makeSubtask("alertPackager")

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        expId, expBits = butlerQC.quantum.dataId.pack("visit_detector",
                                                      returnMaxBits=True)
        inputs["ccdExposureIdBits"] = expBits
        inputs["band"] = butlerQC.quantum.dataId["band"]

        outputs = self.run(**inputs)

        butlerQC.put(outputs, outputRefs)

    @pipeBase.timeMethod
    def run(self,
            diaSourceTable,
            diffIm,
            exposure,
            warpedExposure,
            ccdExposureIdBits,
            band):
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
        warpedExposure : `lsst.afw.image.ExposureF`
            Template exposure used to create diffIm.
        ccdExposureIdBits : `int`
            Number of bits used for a unique ``ccdVisitId``.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            Results struct with components.

            - ``apdb_maker`` : Marker dataset to store in the Butler indicating
              that this ccdVisit has completed successfully.
              (`lsst.dax.apdb.ApdbConfig`)
        """
        self.log.info("Running DiaPipeline...")
        # Put the SciencePipelines through a SDMification step and return
        # calibrated columns with the expect output database names.

        # Load the DiaObjects and DiaSource history.
        loaderResult = self.diaCatalogLoader.run(diffIm, self.apdb)

        # Associate new DiaSources with existing DiaObjects and update
        # DiaObject summary statistics using the full DiaSource history.
        assocResults = self.associator.run(diaSourceTable,
                                           loaderResult.diaObjects,
                                           loaderResult.diaSources)

        mergedDiaSourceHistory = loaderResult.diaSources.append(
            assocResults.diaSources,
            sort=True)
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

        diaCalResult = self.diaCalculation.run(
            assocResults.diaObjects,
            mergedDiaSourceHistory,
            assocResults.matchedDiaObjectIds,
            [band])
        if self.testDataFrameIndex(diaCalResult.diaObjectCat):
            raise RuntimeError(
                "Duplicate DiaObjects (loaded + updated) created after "
                "DiaCalculation. This is unexpected behavior and should be "
                "reported. Existing.")
        if self.testDataFrameIndex(diaCalResult.updatedDiaObjects):
            raise RuntimeError(
                "Duplicate DiaObjects (updated) created after "
                "DiaCalculation. This is unexpected behavior and should be "
                "reported. Existing.")

        # Force photometer on the Difference and Calibrated exposures using
        # the new and updated DiaObject locations.
        diaForcedSources = self.diaForcedSource.run(
            diaCalResult.diaObjectCat,
            diaCalResult.updatedDiaObjects.loc[:, "diaObjectId"].to_numpy(),
            ccdExposureIdBits,
            exposure,
            diffIm)

        # Store DiaSources and updated DiaObjects in the Apdb.
        self.apdb.storeDiaSources(assocResults.diaSources)
        self.apdb.storeDiaObjects(
            diaCalResult.updatedDiaObjects,
            exposure.getInfo().getVisitInfo().getDate().toPython())
        self.apdb.storeDiaForcedSources(diaForcedSources)

        if self.config.doPackageAlerts:
            if len(loaderResult.diaForcedSources) > 1:
                diaForcedSources = diaForcedSources.append(
                    loaderResult.diaForcedSources,
                    sort=True)
            if self.testDataFrameIndex(diaForcedSources):
                self.log.warn(
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
            self.alertPackager.run(assocResults.diaSources,
                                   diaCalResult.diaObjectCat,
                                   loaderResult.diaSources,
                                   diaForcedSources,
                                   diffIm,
                                   warpedExposure,
                                   ccdExposureIdBits)

        return pipeBase.Struct(apdbMarker=self.config.apdb.value,
                               associatedDiaSources=assocResults.diaSources)

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
