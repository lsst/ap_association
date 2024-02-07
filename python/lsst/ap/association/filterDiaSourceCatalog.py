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

__all__ = (
    "FilterDiaSourceCatalogConfig",
    "FilterDiaSourceCatalogTask",
)

import numpy as np

from lsst.afw.table import SourceCatalog
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as connTypes
from lsst.utils.timer import timeMethod


class FilterDiaSourceCatalogConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("instrument", "visit", "detector"),
    defaultTemplates={"coaddName": "deep", "fakesType": ""},
):
    """Connections class for FilterDiaSourceCatalogTask."""

    diaSourceCat = connTypes.Input(
        doc="Catalog of DiaSources produced during image differencing.",
        name="{fakesType}{coaddName}Diff_diaSrc",
        storageClass="SourceCatalog",
        dimensions=("instrument", "visit", "detector"),
    )

    filteredDiaSourceCat = connTypes.Output(
        doc="Output catalog of DiaSources after filtering.",
        name="{fakesType}{coaddName}Diff_candidateDiaSrc",
        storageClass="SourceCatalog",
        dimensions=("instrument", "visit", "detector"),
    )

    rejectedDiaSources = connTypes.Output(
        doc="Optional output storing all the rejected DiaSources.",
        name="{fakesType}{coaddName}Diff_rejectedDiaSrc",
        storageClass="SourceCatalog",
        dimensions={"instrument", "visit", "detector"},
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)
        if not self.config.doWriteRejectedSources:
            self.outputs.remove("rejectedDiaSources")


class FilterDiaSourceCatalogConfig(
    pipeBase.PipelineTaskConfig, pipelineConnections=FilterDiaSourceCatalogConnections
):
    """Config class for FilterDiaSourceCatalogTask."""

    doRemoveSkySources = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Input DiaSource catalog contains SkySources that should be "
        "removed before storing the output DiaSource catalog.",
    )

    doWriteRejectedSources = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Store the output DiaSource catalog containing all the rejected "
        "sky sources."
    )


class FilterDiaSourceCatalogTask(pipeBase.PipelineTask):
    """Filter out sky sources from a DiaSource catalog."""

    ConfigClass = FilterDiaSourceCatalogConfig
    _DefaultName = "filterDiaSourceCatalog"

    @timeMethod
    def run(self, diaSourceCat):
        """Filter sky sources from the supplied DiaSource catalog.

        Parameters
        ----------
        diaSourceCat : `lsst.afw.table.SourceCatalog`
            Catalog of sources measured on the difference image.

        Returns
        -------
        filterResults : `lsst.pipe.base.Struct`

            ``filteredDiaSourceCat`` : `lsst.afw.table.SourceCatalog`
                The catalog of filtered sources.
            ``rejectedDiaSources`` : `lsst.afw.table.SourceCatalog`
                The catalog of rejected sources.
        """
        rejectedSkySources = None
        if self.config.doRemoveSkySources:
            sky_source_column = diaSourceCat["sky_source"]
            num_sky_sources = np.sum(sky_source_column)
            rejectedSkySources = diaSourceCat[sky_source_column].copy(deep=True)
            diaSourceCat = diaSourceCat[~sky_source_column].copy(deep=True)
            self.log.info(f"Filtered {num_sky_sources} sky sources.")
        if not rejectedSkySources:
            rejectedSkySources = SourceCatalog(diaSourceCat.getSchema())
        filterResults = pipeBase.Struct(filteredDiaSourceCat=diaSourceCat,
                                        rejectedDiaSources=rejectedSkySources)
        return filterResults
