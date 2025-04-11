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

from astropy.utils.iers import conf
conf.auto_max_age = None


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

    diffImVisitInfo = connTypes.Input(
        doc="VisitInfo of diffIm.",
        name="{fakesType}{coaddName}Diff_differenceExp.visitInfo",
        storageClass="VisitInfo",
        dimensions=("instrument", "visit", "detector"),
    )

    longTrailedSources = connTypes.Output(
        doc="Optional output temporarily storing long trailed diaSources.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="ArrowAstropy",
        name="{fakesType}{coaddName}Diff_longTrailedSrc",
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)
        if not self.config.doWriteRejectedSkySources:
            self.outputs.remove("rejectedDiaSources")
        if not self.config.doTrailedSourceFilter:
            self.outputs.remove("longTrailedSources")
        if not self.config.doWriteTrailedSources:
            self.outputs.remove("longTrailedSources")


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

    doWriteRejectedSkySources = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Store the output DiaSource catalog containing all the rejected "
        "sky sources."
    )

    doRemoveNegativeDirectImageSources = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Remove DIASources with negative scienceFlux/scienceFluxErr "
        "according to a configurable threshold.",
    )

    minAllowedDirectSnr = pexConfig.Field(
        dtype=float,
        doc="Minimum allowed ratio of scienceFlux/scienceFluxErr.",
        default=-2.0,
    )

    doTrailedSourceFilter = pexConfig.Field(
        doc="Run trailedSourceFilter to remove long trailed sources from the"
            "diaSource output catalog.",
        dtype=bool,
        default=True,
    )

    doWriteTrailedSources = pexConfig.Field(
        doc="Write trailed diaSources sources to a table.",
        dtype=bool,
        default=True,
        deprecated="Trailed sources will not be written out during production."
    )

    max_trail_length = pexConfig.Field(
        dtype=float,
        doc="Length of long trailed sources to remove from the input catalog, "
            "in arcseconds per second. Default comes from DMTN-199, which "
            "requires removal of sources with trails longer than 10 "
            "degrees/day, which is 36000/3600/24 arcsec/second, or roughly"
            "0.416 arcseconds per second.",
        default=36000/3600.0/24.0,
    )


class FilterDiaSourceCatalogTask(pipeBase.PipelineTask):
    """Filter sources from a DiaSource catalog."""

    ConfigClass = FilterDiaSourceCatalogConfig
    _DefaultName = "filterDiaSourceCatalog"

    @timeMethod
    def run(self, diaSourceCat, diffImVisitInfo):
        """Filter sources from the supplied DiaSource catalog.

        Parameters
        ----------
        diaSourceCat : `lsst.afw.table.SourceCatalog`
            Catalog of sources measured on the difference image.
        diffImVisitInfo:  `lsst.afw.image.VisitInfo`
            VisitInfo for the difference image corresponding to diaSourceCat.

        Returns
        -------
        filterResults : `lsst.pipe.base.Struct`

            ``filteredDiaSourceCat`` : `lsst.afw.table.SourceCatalog`
                The catalog of filtered sources.
            ``rejectedDiaSources`` : `lsst.afw.table.SourceCatalog`
                The catalog of rejected sources.
            ``longTrailedDiaSources`` : `astropy.table.Table`
                DiaSources which have trail lengths greater than
                max_trail_length*exposure_time.
        """
        rejectedSources = None
        exposure_time = diffImVisitInfo.exposureTime
        if self.config.doRemoveSkySources:
            sky_source_column = diaSourceCat["sky_source"]
            num_sky_sources = np.sum(sky_source_column)
            rejectedSources = diaSourceCat[sky_source_column].copy(deep=True)
            diaSourceCat = diaSourceCat[~sky_source_column].copy(deep=True)
            self.log.info(f"Filtered {num_sky_sources} sky sources.")

        if self.config.doRemoveNegativeDirectImageSources:
            direct_snr = (diaSourceCat["ip_diffim_forced_PsfFlux_instFlux"]
                          / diaSourceCat["ip_diffim_forced_PsfFlux_instFluxErr"])
            too_negative = direct_snr < self.config.minAllowedDirectSnr
            rejectedNegative = diaSourceCat[too_negative].copy(deep=True)
            diaSourceCat = diaSourceCat[~too_negative].copy(deep=True)
            self.log.info(f"Filtered {np.sum(too_negative)} negative direct sources.")
            if rejectedSources is None:
                rejectedSources = rejectedNegative
            else:
                rejectedSkySources = rejectedSources
                rejectedSources = SourceCatalog(diaSourceCat.getSchema())
                rejectedSources.reserve(len(rejectedSkySources) + len(rejectedNegative))
                rejectedSources.extend(rejectedSkySources, deep=True)
                rejectedSources.extend(rejectedNegative, deep=True)

        if rejectedSources is None:
            rejectedSources = SourceCatalog(diaSourceCat.getSchema())

        if self.config.doTrailedSourceFilter:
            trail_mask = self._check_dia_source_trail(diaSourceCat, exposure_time)
            longTrailedDiaSources = diaSourceCat[trail_mask].copy(deep=True)
            diaSourceCat = diaSourceCat[~trail_mask]

            self.log.info("%i DiaSources exceed max_trail_length %f arcseconds per second, "
                          "dropping from source catalog."
                          % (self.config.max_trail_length, len(diaSourceCat)))
            self.metadata.add("num_filtered", len(longTrailedDiaSources))

            if self.config.doWriteTrailedSources:
                filterResults = pipeBase.Struct(filteredDiaSourceCat=diaSourceCat,
                                                rejectedDiaSources=rejectedSources,
                                                longTrailedSources=longTrailedDiaSources.asAstropy())
            else:
                filterResults = pipeBase.Struct(filteredDiaSourceCat=diaSourceCat,
                                                rejectedDiaSources=rejectedSources)

        else:
            filterResults = pipeBase.Struct(filteredDiaSourceCat=diaSourceCat,
                                            rejectedDiaSources=rejectedSources)

        return filterResults

    def _check_dia_source_trail(self, dia_sources, exposure_time):
        """Find DiaSources that have long trails or trails with indeterminant
        end points.

        Return a mask of sources with lengths greater than
        (``config.max_trail_length`` multiplied by the exposure time)
        arcseconds.
        Additionally, set mask if
        ``ext_trailedSources_Naive_flag_off_image`` is set or if
        ``ext_trailedSources_Naive_flag_suspect_long_trail`` and
        ``ext_trailedSources_Naive_flag_edge`` are both set.

        Parameters
        ----------
        dia_sources : `pandas.DataFrame`
            Input diaSources to check for trail lengths.
        exposure_time : `float`
            Exposure time from difference image.

        Returns
        -------
        trail_mask : `pandas.DataFrame`
            Boolean mask for diaSources which are greater than the
            Boolean mask for diaSources which are greater than the
            cutoff length or have trails which extend beyond the edge of the
            detector (off_image set). Also checks if both
            suspect_long_trail and edge are set and masks those sources out.
        """
        trail_mask = (dia_sources["ext_trailedSources_Naive_length"]
                      >= (self.config.max_trail_length*exposure_time))
        trail_mask |= dia_sources['ext_trailedSources_Naive_flag_off_image']
        trail_mask |= (dia_sources['ext_trailedSources_Naive_flag_suspect_long_trail']
                       & dia_sources['ext_trailedSources_Naive_flag_edge'])

        return trail_mask
