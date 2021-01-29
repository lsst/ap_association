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
use in ap_association and the alert production database (APDB).
"""

__all__ = ["MapApDataConfig", "MapApDataTask",
           "MapDiaSourceConfig", "MapDiaSourceTask",
           "UnpackApdbFlags"]

import numpy as np
import os
import yaml

import lsst.afw.table as afwTable
from lsst.daf.base import DateTime
import lsst.pipe.base as pipeBase
import lsst.pex.config as pexConfig
from lsst.pex.exceptions import RuntimeError
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
    into a product for use in ap_association or the APDB.
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

    @pipeBase.timeMethod
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
                 "ip_diffim_DipoleFit_orientation": "dipAngle",
                 "ip_diffim_DipoleFit_chi2dof": "dipChi2",
                 "ip_diffim_forced_PsfFlux_instFlux": "totFlux",
                 "ip_diffim_forced_PsfFlux_instFluxErr": "totFluxErr",
                 "ip_diffim_DipoleFit_flag_classification": "isDipole",
                 "slot_Shape_xx": "ixx",
                 "slot_Shape_xxErr": "ixxErr",
                 "slot_Shape_yy": "iyy",
                 "slot_Shape_yyErr": "iyyErr",
                 "slot_Shape_xy": "ixy",
                 "slot_Shape_xyErr": "ixyErr",
                 "slot_PsfShape_xx": "ixxPSF",
                 "slot_PsfShape_yy": "iyyPSF",
                 "slot_PsfShape_xy": "ixyPSF"}
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
    dipFluxPrefix = pexConfig.Field(
        dtype=str,
        doc="Prefix of the Dipole measurement column containing negative and "
            "positive flux lobes.",
        default="ip_diffim_DipoleFit",
    )
    dipSepColumn = pexConfig.Field(
        dtype=str,
        doc="Column of the separation of the negative and positive poles of "
            "the dipole.",
        default="ip_diffim_DipoleFit_separation"
    )


class MapDiaSourceTask(MapApDataTask):
    """Task specific for copying columns from science pipelines catalogs,
    calibrating them, for use in ap_association and the APDB.

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
        self.bit_pack_columns = []
        with open(self.config.flagMap) as yaml_stream:
            table_list = list(yaml.safe_load_all(yaml_stream))
            for table in table_list:
                if table['tableName'] == 'DiaSource':
                    self.bit_pack_columns = table['columns']
                    break

        # Test that all flags requested are present in both the input and
        # output schemas.
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

    @pipeBase.timeMethod
    def run(self, inputCatalog, exposure, return_pandas=False):
        """Copy data from the inputCatalog into an output catalog with
        requested columns.

        Parameters
        ----------
        inputCatalog : `lsst.afw.table.SourceCatalog`
            Input catalog with data to be copied into new output catalog.
        exposure: `lsst.afw.image.Exposure`
            Exposure with containing the PhotoCalib object relevant to this
            catalog.
        return_pandas : `bool`
            Return `pandas.DataFrame` instead of `lsst.afw.table.SourceCatalog`

        Returns
        -------
        outputCatalog: `lsst.afw.table.SourceCatalog` or `pandas.DataFrame`
            Output catalog with data copied from input and new column names.
        """
        visit_info = exposure.getInfo().getVisitInfo()
        ccdVisitId = visit_info.getExposureId()
        midPointTaiMJD = visit_info.getDate().get(system=DateTime.MJD)
        filterName = exposure.getFilterLabel().bandLabel
        wcs = exposure.getWcs()

        photoCalib = exposure.getPhotoCalib()

        outputCatalog = afwTable.SourceCatalog(self.outputSchema)
        outputCatalog.reserve(len(inputCatalog))

        for inputRecord in inputCatalog:
            outputRecord = outputCatalog.addNew()
            outputRecord.assign(inputRecord, self.mapper)
            self.calibrateFluxes(inputRecord, outputRecord, photoCalib)
            self.computeDipoleFluxes(inputRecord, outputRecord, photoCalib)
            self.computeDipoleSep(inputRecord, outputRecord, wcs)
            self.bitPackFlags(inputRecord, outputRecord)
            self.computeBBoxSize(inputRecord, outputRecord)
            outputRecord.set("ccdVisitId", ccdVisitId)
            outputRecord.set("midPointTai", midPointTaiMJD)
            outputRecord.set("filterName", filterName)

        if not outputCatalog.isContiguous():
            raise RuntimeError("Output catalogs must be contiguous.")

        if return_pandas:
            return self._convert_to_pandas(outputCatalog)
        return outputCatalog

    def calibrateFluxes(self, inputRecord, outputRecord, photoCalib):
        """Copy flux values into an output record and calibrate them.

        Parameters
        ----------
        inputRecord : `lsst.afw.table.SourceRecord`
            Record to copy flux values from.
        outputRecord : `lsst.afw.table.SourceRecord`
            Record to copy and calibrate values into.
        photoCalib : `lsst.afw.image.PhotoCalib`
            Calibration object from the difference exposure.
        """
        for col_name in self.config.calibrateColumns:
            meas = photoCalib.instFluxToNanojansky(inputRecord, col_name)
            outputRecord.set(self.config.copyColumns[col_name + "_instFlux"],
                             meas.value)
            outputRecord.set(
                self.config.copyColumns[col_name + "_instFluxErr"],
                meas.error)

    def computeDipoleFluxes(self, inputRecord, outputRecord, photoCalib):
        """Calibrate and compute dipole mean flux and diff flux.

        Parameters
        ----------
        inputRecord : `lsst.afw.table.SourceRecord`
            Record to copy flux values from.
        outputRecord : `lsst.afw.table.SourceRecord`
            Record to copy and calibrate values into.
        photoCalib  `lsst.afw.image.PhotoCalib`
            Calibration object from the difference exposure.
        """

        neg_meas = photoCalib.instFluxToNanojansky(
            inputRecord, self.config.dipFluxPrefix + "_neg")
        pos_meas = photoCalib.instFluxToNanojansky(
            inputRecord, self.config.dipFluxPrefix + "_pos")
        outputRecord.set(
            "dipMeanFlux",
            0.5 * (np.abs(neg_meas.value) + np.abs(pos_meas.value)))
        outputRecord.set(
            "dipMeanFluxErr",
            0.5 * np.sqrt(neg_meas.error ** 2 + pos_meas.error ** 2))
        outputRecord.set(
            "dipFluxDiff",
            np.abs(pos_meas.value) - np.abs(neg_meas.value))
        outputRecord.set(
            "dipFluxDiffErr",
            np.sqrt(neg_meas.error ** 2 + pos_meas.error ** 2))

    def computeDipoleSep(self, inputRecord, outputRecord, wcs):
        """Convert the dipole separation from pixels to arcseconds.

        Parameters
        ----------
        inputRecord : `lsst.afw.table.SourceRecord`
            Record to copy flux values from.
        outputRecord : `lsst.afw.table.SourceRecord`
            Record to copy and calibrate values into.
        wcs : `lsst.afw.geom.SkyWcs`
            Wcs of image inputRecords was observed.
        """
        pixScale = wcs.getPixelScale(inputRecord.getCentroid())
        dipSep = pixScale * inputRecord.get(self.config.dipSepColumn)
        outputRecord.set("dipLength", dipSep.asArcseconds())

    def bitPackFlags(self, inputRecord, outputRecord):
        """Pack requested flag columns in inputRecord into single columns in
        outputRecord.

        Parameters
        ----------
        inputRecord : `lsst.afw.table.SourceRecord`
            Record to copy flux values from.
        outputRecord : `lsst.afw.table.SourceRecord`
            Record to copy and calibrate values into.
        """
        for outputFlag in self.bit_pack_columns:
            bitList = outputFlag['bitList']
            value = 0
            for bit in bitList:
                value += inputRecord[bit['name']] * 2 ** bit['bit']
            outputRecord.set(outputFlag['columnName'], value)

    def computeBBoxSize(self, inputRecord, outputRecord):
        """Compute the size of a square bbox that fully contains the detection
        footprint.

        Parameters
        ----------
        inputRecord : `lsst.afw.table.SourceRecord`
            Record to copy flux values from.
        outputRecord : `lsst.afw.table.SourceRecord`
            Record to copy and calibrate values into.
        """
        footprintBBox = inputRecord.getFootprint().getBBox()
        # Compute twice the size of the largest dimension of the footprint
        # bounding box. This is the largest footprint we should need to cover
        # the complete DiaSource assuming the centroid is withing the bounding
        # box.
        maxSize = 2 * np.max([footprintBBox.getWidth(),
                              footprintBBox.getHeight()])
        recX = inputRecord.getCentroid().x
        recY = inputRecord.getCentroid().y
        bboxSize = int(
            np.ceil(2 * np.max(np.fabs([footprintBBox.maxX - recX,
                                        footprintBBox.minX - recX,
                                        footprintBBox.maxY - recY,
                                        footprintBBox.minY - recY]))))
        if bboxSize > maxSize:
            bboxSize = maxSize
        outputRecord.set("bboxSize", bboxSize)

    def _convert_to_pandas(self, inputCatalog):
        """Convert input afw table to pandas.

        Using afwTable.toAstropy().to_pandas() alone is not sufficient to
        properly store data in the Apdb. We must also convert the RA/DEC values
        from radians to degrees and rename several columns.

        Parameters
        ----------
        inputCatalog : `lsst.afw.table.SourceCatalog`
            Catalog to convert to panads and rename columns.

        Returns
        -------
        catalog : `pandas.DataFrame`
        """
        catalog = inputCatalog.asAstropy().to_pandas()
        catalog.rename(columns={"coord_ra": "ra",
                                "coord_dec": "decl",
                                "id": "diaSourceId",
                                "parent": "parentDiaSourceId"},
                       inplace=True)
        catalog["ra"] = np.degrees(catalog["ra"])
        catalog["decl"] = np.degrees(catalog["decl"])

        return catalog


class UnpackApdbFlags:
    """Class for unpacking bits from integer flag fields stored in the Apdb.

    Attributes
    ----------
    flag_map_file : `str`
        Absolute or relative path to a yaml file specifiying mappings of flags
        to integer bits.
    table_name : `str`
        Name of the Apdb table the integer bit data are coming from.
    """

    def __init__(self, flag_map_file, table_name):
        self.bit_pack_columns = []
        with open(flag_map_file) as yaml_stream:
            table_list = list(yaml.safe_load_all(yaml_stream))
            for table in table_list:
                if table['tableName'] == table_name:
                    self.bit_pack_columns = table['columns']
                    break

        self.output_flag_columns = {}

        for column in self.bit_pack_columns:
            names = []
            for bit in column["bitList"]:
                names.append((bit["name"], np.bool))
            self.output_flag_columns[column["columnName"]] = names

    def unpack(self, input_flag_values, flag_name):
        """Determine individual boolean flags from an input array of unsigned
        ints.

        Parameters
        ----------
        input_flag_values : array-like of type uint
            Input integer flags to unpack.
        flag_name : `str`
            Apdb column name of integer flags to unpack. Names of packed int
            flags are given by the flag_map_file.

        Returns
        -------
        output_flags : `numpy.ndarray`
            Numpy named tuple of booleans.
        """
        bit_names_types = self.output_flag_columns[flag_name]
        output_flags = np.zeros(len(input_flag_values), dtype=bit_names_types)

        for bit_idx, (bit_name, dtypes) in enumerate(bit_names_types):
            masked_bits = np.bitwise_and(input_flag_values, 2 ** bit_idx)
            output_flags[bit_name] = masked_bits

        return output_flags
