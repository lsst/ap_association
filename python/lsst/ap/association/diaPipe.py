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

import lsst.dax.apdb as daxApdb
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as connTypes
from lsst.ap.association import (
    AssociationTask,
    DiaForcedSourceTask,
    LoadDiaCatalogsTask,
    MapDiaSourceTask,
    make_dia_object_schema,
    make_dia_source_schema)


class DiaPipelineConnections(pipeBase.PipelineTaskConnections,
                             dimensions=("instrument", "visit", "detector"),
                             defaultTemplates={}):
    """
    """
    diaSourceSchema = connTypes.InitInput(
        doc="",
        name="deepDiff_diaSrc_schema",
        storageClass="SourceCatalog",
        multiple=True
    )
    diaSourceCat = connTypes.Input(
        doc="",
        name="deepDiff_diaSrc",
        storageClass="SourceCatalog",
        dimensions=("instrument", "visit", "detector"),
    )
    diffIm = connTypes.Input(
        doc="",
        name="deepDiff_differenceExp",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector"),
    )
    exposure = connTypes.Input(
        doc="",
        name="calexp",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector"),
    )
    apdbMarker = connTypes.Output(
        doc="",
        name="apdb_marker",
        storageClass="",
        dimensions=("instrument", "visit", "detector"),
    )


class DiaPipelineConfig(pipeBase.PipelineTaskConfig,
                        pipelineConnections=DiaPipelineConnections):
    """
    """
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

    def validate(self):
        pass


class DiaPipelineTask(pipeBase.PipelineTask):
    """
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
                         inputSchema=initInputs["diaSourceSchema"])
        self.makeSubtask("diaCatalogLoader")
        self.makeSubtask("associator")
        self.makeSubtask("diaForcedSource")

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        expId, expBits = butlerQC.quantum.dataId.pack("visit_detector",
                                                      returnMaxBits=True)
        inputs["ccdExposureIdBits"] = expBits

        outputs = self.run(**inputs)

        butlerQC.put(outputs, outputRefs)

    @pipeBase.timedMethod
    def run(self, diaSourceCat, diffIm, exposure, ccdExposureIdBits):
        self.log.info("Running Association...")
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
                                           loaderResult.diaSources,
                                           diffIm)

        # Force photometer on the Difference and Calibrated exposures using
        # the new and updated DiaObject locations.
        diaForcedSources = self.diaForcedSource.run(assocResults.diaObjects,
                                                    ccdExposureIdBits,
                                                    exposure,
                                                    diffIm)

        # Store DiaSources and updated DiaObjects in the Apdb.
        self.apdb.storeDiaSources(assocResults.diaSources)
        self.apdb.storeDiaObjects(
            assocResults.UpdatedDiaObjects,
            exposure.getInfo().getVisitInfo().getDate().toPython())
        self.apdb.storeDiaForcedSources(diaForcedSources)

        return pipeBase.Struct(apdb_maker=self.config.apdb.value)
