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

"""Classes for taking science pipeline outputs and creating data products for
use in ap_association and the prompt-products database (PPDB).
"""

__all__ = ["MapApDataConfig", "MapApDataTask",
           "MapDiaSourceConfig", "MapDiaSourceTask"]

import os
import yaml

import lsst.afw.table as afwTable
from lsst.daf.base import DateTime
import lsst.pipe.base as pipeBase
import lsst.pex.config as pexConfig
from lsst.pex.exceptions import RuntimeError
import lsst.afw.image as afwImage
from lsst.utils import getPackageDir
from .afwUtils import make_dia_source_schema


class MapApDataConfig(pexConfig.Config):
    """Configuration for the generic MapApDataTask class.
    """
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


class MapApDataTask(pipeBase.Task):
    """Generic mapper class for copying values from a science pipelines catalog
    into a product for use in ap_association or the PPDB.
    """
    ConfigClass = MapApDataConfig
    _DefaultName = "mapApDataTask"

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
        inputCatalog: `lsst.afw.table.SourceCatalog`
           Input catalog with data to be copied into new output catalog.

        Returns
        -------
        outputCatalog: `lsst.afw.table.SourceCatalog`
            Output catalog with data copied from input and new column names.
        """
        outputCatalog = afwTable.SourceCatalog(self.outputSchema)
        outputCatalog.extend(inputCatalog, self.mapper)

        if not outputCatalog.isContiguous():
            raise RuntimeError("Output catalogs must be contiguous.")

        return outputCatalog


class MapDiaSourceConfig(pexConfig.Config):
    """Config for the DiaSourceMapperTask
    """
    copyColumns = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Mapping of input SciencePipelines columns to output DPDD "
            "columns.",
        default={"id": "id",
                 "parent": "parent",
                 "coord_ra": "coord_ra",
                 "coord_dec": "coord_dec",
                 "slot_Centroid_x": "x",
                 "slot_Centroid_xErr": "xErr",
                 "slot_Centroid_y": "y",
                 "slot_Centroid_yErr": "yErr",
                 "slot_ApFlux_instFlux": "apFlux",
                 "slot_ApFlux_instFluxErr": "apFluxErr",
                 "slot_PsfFlux_instFlux": "psFlux",
                 "slot_PsfFlux_instFluxErr": "psFluxErr",
                 "ip_diffim_forced_PsfFlux_instFlux": "totFlux",
                 "ip_diffim_forced_PsfFlux_instFluxErr": "totFluxErr"}
    )
    calibrateColumns = pexConfig.ListField(
        dtype=str,
        doc="Flux columns in the input catalog to calibrate.",
        default=["slot_ApFlux", "slot_PsfFlux", "ip_diffim_forced_PsfFlux"]
    )
    flagMap = pexConfig.Field(
        dtype=str,
        doc="Yaml file specifying SciencePipelines flag fields to bit packs.",
        default=os.path.join(getPackageDir("ap_association"),
                             "data",
                             "association-flag-map.yaml"),
    )


class MapDiaSourceTask(MapApDataTask):
    """Task specific for copying columns from science pipelines catalogs,
    calibrating them, for use in ap_association and the PPDB.

    This task also copies information from the exposure such as the ExpsoureId
    and the exposure date as specified in the DPDD.
    """

    ConfigClass = MapDiaSourceConfig
    _DefaultName = "mapDiaSourceTask"

    def __init__(self, inputSchema, **kwargs):
        MapApDataTask.__init__(self,
                               inputSchema=inputSchema,
                               outputSchema=make_dia_source_schema(),
                               **kwargs)
        self._create_bit_pack_mappings()

    def _create_bit_pack_mappings(self):
        """Setup all flag bit packings.
        """
        output_columns = []
        with open(self.config.flagMap) as yaml_stream:
            table_list = list(yaml.load_all(yaml_stream))
            for table in table_list:
                if table['tableName'] == 'DiaSource':
                    output_columns = table['columns']
                    break
        self.bit_pack_columns = output_columns

        for outputFlag in self.bit_pack_columns:
            try:
                self.outputSchema.find(outputFlag['columnName'])
            except KeyError:
                raise KeyError(
                    "Requested column %s not found in MapDiaSourceTask output "
                    "schema. Please check that the requested output column "
                    "exists." % outputFlag['columnName'])
            bitList = outputFlag['bitList']
            for bit in bitList:
                try:
                    self.inputSchema.find(bit['name'])
                except KeyError:
                    raise KeyError(
                        "Requested column %s not found in MapDiaSourceTask input "
                        "schema. Please check that the requested input column "
                        "exists." % outputFlag['columnName'])

    def run(self, inputCatalog, exposure):
        """Copy data from the inputCatalog into an output catalog with
        requested columns.

        Parameters
        ----------
        inputCatalog: `lsst.afw.table.SourceCatalog`
            Input catalog with data to be copied into new output catalog.
        exposure: `lsst.afw.image.Exposure`
            Exposure with containing the calib object relevant to this catalog.
        Returns
        -------
        outputCatalog: `lsst.afw.table.SourceCatalog`
            Output catalog with data copied from input and new column names.
        """
        visit_info = exposure.getInfo().getVisitInfo()
        ccdVisitId = visit_info.getExposureId()
        midPointTaiMJD = visit_info.getDate().get(system=DateTime.MJD)
        filterId = exposure.getFilter().getId()
        filterName = exposure.getFilter().getName()

        flux0, flux0Err = exposure.getCalib().getFluxMag0()
        photoCalib = afwImage.PhotoCalib(1 / flux0, flux0Err / flux0 ** 2)

        outputCatalog = afwTable.SourceCatalog(self.outputSchema)
        outputCatalog.reserve(len(inputCatalog))

        for inputRecord in inputCatalog:
            outputRecord = outputCatalog.addNew()
            outputRecord.assign(inputRecord, self.mapper)
            self.calibrateFluxes(inputRecord, outputRecord, photoCalib)
            self.bitPackFlags(inputRecord, outputRecord)

            outputRecord.set("ccdVisitId", ccdVisitId)
            outputRecord.set("midPointTai", midPointTaiMJD)
            outputRecord.set("filterId", filterId)
            outputRecord.set("filterName", filterName)

        if not outputCatalog.isContiguous():
            raise RuntimeError("Output catalogs must be contiguous.")

        return outputCatalog

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
            meas = photoCalib.instFluxToMaggies(inputRecord, col_name)
            outputRecord.set(self.config.copyColumns[col_name + "_instFlux"],
                             meas.value)
            outputRecord.set(
                self.config.copyColumns[col_name + "_instFluxErr"],
                meas.err)

    def bitPackFlags(self, inputRecord, outputRecord):
        """.

        Parameters
        ----------
        inputRecord: `lsst.afw.table.SourceRecord`
            Record to copy flux values from.
        outputRecord: `lsst.afw.table.SourceRecord`
            Record to copy and calibrate values into.
        """
        for outputFlag in self.bit_pack_columns:
            bitList = outputFlag['bitList']
            value = 0
            for bit in bitList:
                value += inputRecord[bit['name']] * 2 ** bit['bit']
            outputRecord.set(outputFlag['columnName'], value)
