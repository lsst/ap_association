#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__all__ = ["APDataMapperConfig", "APDataMapperTask",
           "DiaSourceMapperConfig", "DiaSourceMapperTask"]

import lsst.afw.table as afwTable
import lsst.pipe.base as pipeBase
import lsst.pex.config as pexConfig
import lsst.afw.image as afwImage


class APDataMapperConfig(pexConfig.Config):

    copyColumns = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Mapping of input SciencePipelines columns to output DPDD "
            "columns.",
        default={"id": "id",
                 "parent": "parent",
                 "coord_ra": "coord_ra",
                 "coord_dec": "coord_dec"}
    )


class APDataMapperTask(pipeBase.Task):

    ConfigClass = APDataMapperConfig
    _DefaultName = "apDataMapperTask"

    def __init__(self, inputSchema, outputSchema, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.inputSchema = inputSchema
        self.outputSchema = outputSchema

        self.mapper = afwTable.SchemaMapper(inputSchema, outputSchema)

        for inputName, outputName in self.config.copyColumns.items():
            self.mapper.addMapping(
                self.inputSchema.find(inputName).key,
                outputName,
                True)

    def run(self, inputCatalog, exposure=None):
        """Copy data from the inputCatalog into an output catalog with
        requested columns.

        Parameters
        ----------
        inputCatalog: `lsst.afw.table.SourceCatalog``
           Input catalog with data to be copied into new output catalog.

        Returns
        -------
        outputCatalog: `lsst.afw.table.SourceCatalog`
            Output catalog with data copied from input and new column names.
        """
        outputCatalog = afwTable.SourceCatalog(self.outputSchema)
        outputCatalog.preallocate(len(inputCatalog))

        for inputRecord in inputCatalog:
            outputRecord = outputCatalog.addNew()
            outputRecord.assign(inputRecord, self.mapper)

        if outputCatalog.isContiguous():
            return outputCatalog

        return outputCatalog.copy(deep=True)


class DiaSourceMapperConfig(pexConfig.Config):

    copyColumns = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Mapping of input SciencePipelines columns to output DPDD "
            "columns.",
        default={"id": "id",
                 "parent": "parent",
                 "coord_ra": "coord_ra",
                 "coord_dec": "coord_dec",
                 "slot_ApFlux_instFlux": "apFlux",
                 "slot_ApFlux_instFluxErr": "apFluxErr",
                 "slot_PsfFlux_instFlux": "psFlux",
                 "slot_PsfFlux_instFluxErr": "psFluxErr"}
    )
    calibrateColumns = pexConfig.ListField(
        dtype=str,
        doc="Flux columns in the input catalog to calibrate.",
        default=["slot_ApFlux", "slot_PsfFlux", "slot_dipMeanFlux",
                 "slot_dipFluxDiff", "slot_totFlux"]
    )


class DiaSourceMapperTask(APDataMapperTask):

    ConfigClass = DiaSourceMapperConfig
    _DefaultName = "diaSourceMapper"

    def __init__(self, inputSchema, **kwargs):
        APDataMapperTask.__init__(self,
                                  inputSchema,
                                  make_dia_source_schema(),
                                  **kwargs)      

    
    def run(self, inputCatalog, exposure):
        """Copy data from the inputCatalog into an output catalog with
        requested columns.

        Parameters
        ----------
        inputCatalog: `lsst.afw.table.SourceCatalog``
           Input catalog with data to be copied into new output catalog.

        Returns
        -------
        outputCatalog: `lsst.afw.table.SourceCatalog`
            Output catalog with data copied from input and new column names.
        """
        visit_info = exposure.getInfo().getVisitInfo()
        ccdVisitId =  visit_info.getExposureId()
        midPointTaiMJD = visit_info.getDate().get(system=DateTime.MJD)

        flux0, flux0Err = exposure.getCalib().getFluxMag0()
        photoCalib = afwImage.PhotoCalib(1 / flux0, 1 / flux0Err)

        outputCatalog = afwTable.SourceCatalog(self.outputSchema)
        outputCatalog.preallocate(len(inputCatalog))

        for inputRecord in inputCatalog:
            outputRecord = outputCatalog.addNew()
            outputRecord.assign(inputRecord, self.mapper)
            self.calibrateFluxes(inputRecord, outputRecord, photoCalib)
     
            outputRecord.set("ccdVisitId", ccdVisitId)
            outputRecord.set("midPointTai", midPointTaiMJD)

        if outputCatalog.isContigious():
            return outputCatalog

        return outputCatalog.copy(deep=True)

    def calibrateFluxes(self, inputRecord, outputRecord, photoCalib):
        """Copy flux values into an output record and calibrate them.

        Parameters
        ----------
        inputRecord: `lsst.afw.table.SourceRecord`
            Record to copy flux values from.
        outputRecord: `lsst.afw.table.SourceRecord`
            Record to copy and calibrate values into.
        photoCalib: `lsst.afw.image.PhotoCalib`
            Calibration object from the difference exposure.
        """
        for col_name in self.config.calibrateColumns:
            val, valErr = photoCalib.instFluxToMaggies(inputRecord, col_name)
            outputRecord.set(self.config.copyColumns[col_name + "_instFlux"],
                             val)
            outputRecord.set(
                self.config.copyColumns[col_name + "_instFluxErr"] + "Err",
                valErr)
