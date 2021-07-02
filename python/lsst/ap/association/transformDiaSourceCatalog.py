# This file is part of ap_association
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

__all__ = ("TransformDiaSourceCatalogConnections",
           "TransformDiaSourceCatalogConfig",
           "TransformDiaSourceCatalogTask",
           "UnpackApdbFlags")

import numpy as np
import os
import yaml

from lsst.daf.base import DateTime
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as connTypes
from lsst.pipe.tasks.postprocess import TransformCatalogBaseTask
from lsst.pipe.tasks.parquetTable import ParquetTable
from lsst.utils import getPackageDir


class TransformDiaSourceCatalogConnections(pipeBase.PipelineTaskConnections,
                                           dimensions=("instrument", "visit", "detector"),
                                           defaultTemplates={"coaddName": "deep", "fakesType": ""}):
    """Butler connections for TransformDiaSourceCatalogTask.
    """
    diaSourceSchema = connTypes.InitInput(
        doc="Schema for DIASource catalog output by ImageDifference.",
        storageClass="SourceCatalog",
        name="{fakesType}{coaddName}Diff_diaSrc_schema",
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
    diaSourceTable = connTypes.Output(
        doc=".",
        name="{fakesType}{coaddName}Diff_diaSrcTable",
        storageClass="DataFrame",
        dimensions=("instrument", "visit", "detector"),
    )


class TransformDiaSourceCatalogConfig(pipeBase.PipelineTaskConfig,
                                      pipelineConnections=TransformDiaSourceCatalogConnections):
    """
    """
    flagMap = pexConfig.Field(
        dtype=str,
        doc="Yaml file specifying SciencePipelines flag fields to bit packs.",
        default=os.path.join(getPackageDir("ap_association"),
                             "data",
                             "association-flag-map.yaml"),
    )
    functorFile = pexConfig.Field(
        dtype=str,
        doc='Path to YAML file specifying Science DataModel functors to use '
            'when copying columns and computing calibrated values.',
        default=os.path.join("${AP_ASSOCIATION_DIR}",
                             "data",
                             "DiaSource.yaml")
    )


class TransformDiaSourceCatalogTask(TransformCatalogBaseTask):
    """Apply Science DataModel-ification on the DiaSource afw table.

    This task calibrates and renames columns in the DiaSource catalog
    to ready the catalog for insertion into the Apdb.

    This is a Gen3 Butler only task. It will not run in Gen2.
    """

    ConfigClass = TransformDiaSourceCatalogConfig
    _DefaultName = "transformDiaSourceCatalog"
    RunnerClass = pipeBase.ButlerInitializedTaskRunner

    def __init__(self, initInputs, **kwargs):
        super().__init__(**kwargs)
        self.funcs = self.getFunctors()
        self.inputSchema = initInputs['diaSourceSchema'].schema
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

        # Test that all flags requested are present in the input schemas.
        # Output schemas are flexible, however if names are not specified in
        # the Apdb schema, flag columns will not be persisted.
        for outputFlag in self.bit_pack_columns:
            bitList = outputFlag['bitList']
            for bit in bitList:
                try:
                    self.inputSchema.find(bit['name'])
                except KeyError:
                    raise KeyError(
                        "Requested column %s not found in input DiaSource "
                        "schema. Please check that the requested input "
                        "column exists." % bit['name'])

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        expId, expBits = butlerQC.quantum.dataId.pack("visit_detector",
                                                      returnMaxBits=True)
        inputs["ccdVisitId"] = expId
        inputs["band"] = butlerQC.quantum.dataId["band"]

        outputs = self.run(**inputs)

        butlerQC.put(outputs, outputRefs)

    @pipeBase.timeMethod
    def run(self,
            diaSourceCat,
            diffIm,
            band,
            ccdVisitId,
            funcs=None):
        """Convert input catalog to ParquetTable/Pandas and run functors.

        Additionally, add new columns for stripping information from the
        exposure and into the DiaSource catalog.

        Parameters
        ----------

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            Results struct with components.

            - ``diaSourceTable`` : Catalog of DiaSources with calibrated values
              and renamed columns.
              (`lsst.pipe.tasks.ParquetTable` or `pandas.DataFrame`)
        """
        self.log.info(
            "Transforming/standardizing the DiaSource table ccdVisitId: %i",
            ccdVisitId)

        diaSourceDf = diaSourceCat.asAstropy().to_pandas()
        diaSourceDf["bboxSize"] = self.computeBBoxSizes(diaSourceCat)
        diaSourceDf["ccdVisitId"] = ccdVisitId
        diaSourceDf["filterName"] = band
        diaSourceDf["midPointTai"] = diffIm.getInfo().getVisitInfo().getDate().get(system=DateTime.MJD)
        diaSourceDf["diaObjectId"] = 0
        diaSourceDf["pixelId"] = 0
        self.bitPackFlags(diaSourceDf)

        df = self.transform(band,
                            ParquetTable(dataFrame=diaSourceDf),
                            self.funcs,
                            dataId=None).df
        # The Ra/DecColumn functors preserve the coord_ra/dec original columns.
        # Since we don't need these and keeping them causes a DB insert crash
        # we drop them from the DataFrame before returning output catalog.
        return pipeBase.Struct(
            diaSourceTable=df.drop(columns=["coord_ra", "coord_dec"]),
        )

    def computeBBoxSizes(self, inputCatalog):
        """Compute the size of a square bbox that fully contains the detection
        footprint.

        Parameters
        ----------
        inputCatalog : `lsst.afw.table.SourceCatalog`
            Catalog containing detected footprints.

        Returns
        -------
        outputBBoxSizes : `numpy.ndarray`, (N,)
            Array of bbox sizes.
        """
        outputBBoxSizes = np.empty(len(inputCatalog), dtype=int)
        for idx, record in enumerate(inputCatalog):
            footprintBBox = record.getFootprint().getBBox()
            # Compute twice the size of the largest dimension of the footprint
            # bounding box. This is the largest footprint we should need to cover
            # the complete DiaSource assuming the centroid is withing the bounding
            # box.
            maxSize = 2 * np.max([footprintBBox.getWidth(),
                                  footprintBBox.getHeight()])
            recX = record.getCentroid().x
            recY = record.getCentroid().y
            bboxSize = int(
                np.ceil(2 * np.max(np.fabs([footprintBBox.maxX - recX,
                                            footprintBBox.minX - recX,
                                            footprintBBox.maxY - recY,
                                            footprintBBox.minY - recY]))))
            if bboxSize > maxSize:
                bboxSize = maxSize
            outputBBoxSizes[idx] = bboxSize

        return outputBBoxSizes

    def bitPackFlags(self, df):
        """Pack requested flag columns in inputRecord into single columns in
        outputRecord.

        Parameters
        ----------
        df : `pandas.DataFrame`
            DataFrame to read bits from and pack them into.
        """
        for outputFlag in self.bit_pack_columns:
            bitList = outputFlag['bitList']
            value = np.zeros(len(df), dtype=np.uint64)
            for bit in bitList:
                # Hard type the bit arrays.
                value += (df[bit['name']]*2**bit['bit']).to_numpy().astype(np.uint64)
            df[outputFlag['columnName']] = value


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
                names.append((bit["name"], bool))
            self.output_flag_columns[column["columnName"]] = names

    def unpack(self, input_flag_values, flag_name):
        """Determine individual boolean flags from an input array of unsigned
        ints.

        Parameters
        ----------
        input_flag_values : array-like of type uint
            Array of integer flags to unpack.
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
            masked_bits = np.bitwise_and(input_flag_values, 2**bit_idx)
            output_flags[bit_name] = masked_bits

        return output_flags
