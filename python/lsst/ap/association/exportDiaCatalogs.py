# This file is part of ap_association.
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

"""Load DiaSources and DiaObjects by patch and tract for reassociation.
"""


import astropy.units

from lsst.daf.base import DateTime
import lsst.daf.butler as dafButler
import lsst.dax.apdb as daxApdb
import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as connTypes
from lsst.pipe.base.utils import RegionTimeInfo
from lsst.skymap import BaseSkyMap
import lsst.sphgeom

from lsst.ap.association.loadDiaCatalogs import LoadDiaCatalogsTask

__all__ = ("ExportDiaCatalogsTask", "ExportDiaCatalogsConfig")


class ExportDiaCatalogsConnections(pipeBase.PipelineTaskConnections,
                                   dimensions=("tract", "patch", "skymap")):

    skyMap = pipeBase.connectionTypes.Input(
        doc="Geometry of the tracts and patches that the coadds are defined on.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        dimensions=("skymap",),
        storageClass="SkyMap",
    )
    coaddExposures = pipeBase.connectionTypes.Input(
        doc="Coadds for the current patch and tract. There may be multiple "
            "bands but we will only use the dataId of the first that is found."
            "Used to constrain the quantum graph and only attempt to load data "
            "from the APDB for patches with templates.",
        dimensions=("tract", "patch", "skymap", "band"),
        storageClass="ExposureF",
        name="template_coadd",
        multiple=True,
        deferLoad=True,
    )
    diaObjects = connTypes.Output(
        doc="DiaObjects preloaded from the APDB.",
        name="apdb_export_diaObjects",
        storageClass="ArrowAstropy",
        dimensions=("tract", "patch"),
    )
    diaSources = connTypes.Output(
        doc="DiaSources preloaded from the APDB.",
        name="apdb_export_diaSources",
        storageClass="ArrowAstropy",
        dimensions=("tract", "patch"),
    )
    diaForcedSources = connTypes.Output(
        doc="DiaForcedSources preloaded from the APDB.",
        name="apdb_export_diaForcedSources",
        storageClass="ArrowAstropy",
        dimensions=("tract", "patch"),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)
        if not config.doLoadForcedSources:
            del self.diaForcedSources


class ExportDiaCatalogsConfig(pipeBase.PipelineTaskConfig,
                              pipelineConnections=ExportDiaCatalogsConnections):
    """Config class for ExportDiaCatalogsConfig.
    """
    apdb_config_url = pexConfig.Field(
        dtype=str,
        default=None,
        optional=False,
        doc="A config file specifying the APDB and its connection parameters, "
            "typically written by the apdb-cli command-line utility. "
            "The database must already be initialized.",
    )
    angleMargin = pexConfig.RangeField(
        doc="Padding to add to the radius of the bounding circle (arcseconds)",
        dtype=float,
        default=2,
        min=0,
    )
    doLoadForcedSources = pexConfig.Field(
        dtype=bool,
        default=True,
        deprecated="Added to allow disabling forced sources for performance "
                   "reasons during the ops rehearsal. "
                   "It is expected to be removed.",
        doc="Load forced DiaSource history from the APDB? "
            "This should only be turned off for debugging purposes.",
    )


class ExportDiaCatalogsTask(LoadDiaCatalogsTask):
    """Retrieve DiaObjects and associated DiaSources from the Apdb given an
    input exposure.
    """
    ConfigClass = ExportDiaCatalogsConfig
    _DefaultName = "exportDiaCatalogs"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.apdb = daxApdb.Apdb.from_uri(self.config.apdb_config_url)
        self.apdb.read_sources_months = 12

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        skymap = inputs.pop("skyMap")
        coaddExposures = inputs.pop("coaddExposures")
        dataId = coaddExposures[0].dataId

        # region = lsst.sphgeom.ConvexPolygon([pp.getVector() for pp in wcs.pixelToSky(bbox.getCorners())])
        regionTime = self._makeRegionTime(skymap, dataId["tract"], dataId["patch"])
        outputs = self.run(regionTime)
        butlerQC.put(outputs, outputRefs)

    @staticmethod
    def _makeRegionTime(skymap, tract, patch):
        """Construct a region and timespan to load catalogs from the APDB.

        Parameters
        ----------
        skymap : `lsst.skymap.SkyMap`
            Geometry of the tracts and patches the coadds are defined on.
        tract : `int`
            Tract id number.
        patch : `int`
            Patch id number.

        Returns
        -------
        regionTime : `lsst.pipe.base.utils.RegionTimeInfo`
            A serializable container for a sky region and timespan.
        """
        wcs = skymap[tract].getWcs()
        bbox = skymap[tract][patch].getOuterBBox()
        patchCorners = wcs.pixelToSky(lsst.geom.Box2D(bbox).getCorners())
        region = lsst.sphgeom.ConvexPolygon([pp.getVector() for pp in patchCorners])

        begin = DateTime.now().toAstropy()
        end = begin + 30*astropy.units.second
        timespan = dafButler.Timespan(begin=begin, end=end)
        return RegionTimeInfo(region=region, timespan=timespan)
