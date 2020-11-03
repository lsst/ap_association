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
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as connTypes
from lsst.utils import getPackageDir

from lsst.ap.association import (
    AssociationTask,
    DiaForcedSourceTask,
    LoadDiaCatalogsTask,
    MapDiaSourceTask,
    make_dia_object_schema,
    make_dia_source_schema,
    PackageAlertsTask)

__all__ = ("DiaPipelineConfig",
           "DiaPipelineTask",
           "DiaPipelineConnections")


class DiaPipelineConnections(pipeBase.PipelineTaskConnections,
                             dimensions=("instrument", "visit", "detector"),
                             defaultTemplates={"coaddName": "deep", "fakesType": ""}):
    """Butler connections for DiaPipelineTask.
    """
    diaSourceSchema = connTypes.InitInput(
        doc="Schema of the DiaSource catalog produced during image "
            "differencing",
        name="{fakesType}{coaddName}Diff_diaSrc_schema",
        storageClass="SourceCatalog",
        multiple=True
    )
    diaSourceCat = connTypes.Input(
        doc="Catalog of DiaSources produced during image differencing.",
        name="{fakesType}{coaddName}Diff_diaSrc",
        storageClass="SourceCatalog",
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
        doc="Optional output storing the DiaSource catalog after matching and "
            "SDMification.",
        name="{fakesType}{coaddName}Diff_assocDiaSrc",
        storageClass="DataFrame",
        dimensions=("instrument", "visit", "detector"),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not config.doWriteAssociatedSources:
            self.outputs.remove("associatedDiaSources")


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
    diaSourceDpddifier = pexConfig.ConfigurableField(
        target=MapDiaSourceTask,
        doc="Task for assigning columns from the raw output of ip_diffim into "
            "a schema that more closely resembles the DPDD.",
    )
    diaCatalogLoader = pexConfig.ConfigurableField(
        target=LoadDiaCatalogsTask,
        doc="Task to load DiaObjects and DiaSources from the Apdb.",
    )
    associator = pexConfig.ConfigurableField(
        target=AssociationTask,
        doc="Task used to associate DiaSources with DiaObjects.",
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
            getPackageDir("ap_association"),
            "data",
            "apdb-ap-pipe-schema-extra.yaml")

    def validate(self):
        pexConfig.Config.validate(self)
        if self.diaCatalogLoader.htmLevel != \
           self.associator.diaCalculation.plugins["ap_HTMIndex"].htmLevel:
            raise ValueError("HTM index level in LoadDiaCatalogsTask must be "
                             "equal to HTMIndexDiaCalculationPlugin index "
                             "level.")


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
        self.makeSubtask("diaSourceDpddifier",
                         inputSchema=initInputs["diaSourceSchema"].schema)
        self.makeSubtask("diaCatalogLoader")
        self.makeSubtask("associator")
        self.makeSubtask("diaForcedSource")
        if self.config.doPackageAlerts:
            self.makeSubtask("alertPackager")

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        expId, expBits = butlerQC.quantum.dataId.pack("visit_detector",
                                                      returnMaxBits=True)
        inputs["ccdExposureIdBits"] = expBits

        outputs = self.run(**inputs)

        butlerQC.put(outputs, outputRefs)

    @pipeBase.timeMethod
    def run(self, diaSourceCat, diffIm, exposure, warpedExposure, ccdExposureIdBits):
        """Process DiaSources and DiaObjects.

        Load previous DiaObjects and their DiaSource history. Calibrate the
        values in the diaSourceCat. Associate new DiaSources with previous
        DiaObjects. Run forced photometry at the updated DiaObject locations.
        Store the results in the Alert Production Database (Apdb).

        Parameters
        ----------
        diaSourceCat : `lsst.afw.table.SourceCatalog`
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
        diaSources = self.diaSourceDpddifier.run(diaSourceCat,
                                                 diffIm,
                                                 return_pandas=True)

        # Load the DiaObjects and DiaSource history.
        loaderResult = self.diaCatalogLoader.run(diffIm, self.apdb)

        # Associate new DiaSources with existing DiaObjects and update
        # DiaObject summary statistics using the full DiaSource history.
        assocResults = self.associator.run(diaSources,
                                           loaderResult.diaObjects,
                                           loaderResult.diaSources)

        # Force photometer on the Difference and Calibrated exposures using
        # the new and updated DiaObject locations.
        diaForcedSources = self.diaForcedSource.run(
            assocResults.diaObjects,
            assocResults.updatedDiaObjects.loc[:, "diaObjectId"].to_numpy(),
            ccdExposureIdBits,
            exposure,
            diffIm)

        # Store DiaSources and updated DiaObjects in the Apdb.
        self.apdb.storeDiaSources(assocResults.diaSources)
        self.apdb.storeDiaObjects(
            assocResults.updatedDiaObjects,
            exposure.getInfo().getVisitInfo().getDate().toPython())
        self.apdb.storeDiaForcedSources(diaForcedSources)
        if self.config.doPackageAlerts:
            if len(loaderResult.diaForcedSources) > 1:
                diaForcedSources = diaForcedSources.append(
                    loaderResult.diaForcedSources,
                    sort=True)
            self.alertPackager.run(assocResults.diaSources,
                                   assocResults.diaObjects,
                                   loaderResult.diaSources,
                                   diaForcedSources,
                                   diffIm,
                                   warpedExposure,
                                   ccdExposureIdBits)

        return pipeBase.Struct(apdbMarker=self.config.apdb.value,
                               associatedDiaSources=assocResults.diaSources)
