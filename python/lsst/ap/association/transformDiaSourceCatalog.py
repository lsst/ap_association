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
           "TransformDiaSourceCatalogTask")

import numpy as np
import os

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
    functorFile = pexConfig.Field(
        dtype=str,
        doc='Path to YAML file specifying Science DataModel functors to use '
            'when copying columns and computing calibrated values.',
        default=os.path.join(getPackageDir("ap_association"),
                             "data",
                             "DiaSource.yaml")
    )


class TransformDiaSourceCatalogTask(pipeBase.PipelineTask,
                                    TransformCatalogBaseTask):
    """Apply Science DataModel-ification on the DiaSource afw table.

    This task calibrates and renames columns in the DiaSource catalog
    to ready the catalog for insertion into the Apdb.

    This is a Gen3 Butler only task. It will not run in Gen2.
    """

    ConfigClass = TransformDiaSourceCatalogConfig
    _DefaultName = "transformDiaSourceCatalog"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.funcs = self.getFunctors()

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        expId, expBits = butlerQC.quantum.dataId.pack("visit_detector",
                                                      returnMaxBits=True)
        inputs["ccdVisitId"] = expId
        inputs["band"] = butlerQC.quantum.dataId["band"]

        outputs = self.run(**inputs)

        butlerQC.put(outputs, outputRefs)

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

        df = self.transform(band,
                            ParquetTable(dataFrame=diaSourceDf),
                            self.funcs,
                            dataId=None).df
        return pipeBase.Struct(
            diaSourceTable=df
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
