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

"""Methods for force photometering direct and difference images at DiaObject
locations.
"""

__all__ = ["DiaForcedSourceTask", "DiaForcedSourcedConfig"]

from lsst.meas.base.pluginRegistry import register
from lsst.meas.base import (
    ForcedTransformedCentroidConfig,
    ForcedTransformedCentroidPlugin)
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
from lsst.meas.base import ForcedMeasurementTask

from .afwUtils import make_dia_object_schema, make_dia_forced_source_schema


class ForcedTransformedCentroidFromCoordConfig(ForcedTransformedCentroidConfig):
    """Configuration for the forced transformed coord algorithm.
    """
    pass


@register("ap_assoc_TransformedCentroid")
class ForcedTransformedCentroidFromCoordPlugin(ForcedTransformedCentroidPlugin):
    """Record the transformation of the reference catalog coord.
    The coord recorded in the reference catalog is tranformed to the
    measurement coordinate system and stored.

    Parameters
    ----------
    config : `ForcedTransformedCentroidFromCoordConfig`
        Plugin configuration
    name : `str`
        Plugin name
    schemaMapper : `lsst.afw.table.SchemaMapper`
        A mapping from reference catalog fields to output
        catalog fields. Output fields are added to the output schema.
    metadata : `lsst.daf.base.PropertySet`
        Plugin metadata that will be attached to the output catalog.

    Notes
    -----
    This can be used as the slot centroid in forced measurement when only a
    reference coord exits, allowing subsequent measurements to simply refer to
    the slot value just as they would in single-frame measurement.
    """

    ConfigClass = ForcedTransformedCentroidFromCoordConfig

    def measure(self, measRecord, exposure, refRecord, refWcs):
        targetWcs = exposure.getWcs()

        targetPos = targetWcs.skyToPixel(refRecord.getCoord())
        measRecord.set(self.centroidKey, targetPos)

        if self.flagKey is not None:
            measRecord.set(self.flagKey, refRecord.getCentroidFlag())


class DiaForcedSourcedConfig(pexConfig.Config):
    """Configuration for the generic DiaForcedSourcedTask class.
    """
    forcedMeasurement = pexConfig.ConfigurableField(
        target=ForcedMeasurementTask,
        doc="Subtask to force photometer DiaObjects in the direct and "
            "difference images.",
    )

    def setDefaults(self):
        self.forcedMeasurement.plugins = ["ap_assoc_TransformedCentroid",
                                          "base_PsfFlux"]
        self.forcedMeasurement.doReplaceWithNoise = False
        self.forcedMeasurement.copyColumns = {
            "id": "objectId",
            "coord_ra": "coord_ra",
            "coord_dec": "coord_dec"}
        self.forcedMeasurement.slots.centroid = "ap_assoc_TransformedCentroid"
        self.forcedMeasurement.slots.psfFlux = "base_PsfFlux"
        self.forcedMeasurement.slots.shape = None


class DiaForcedSourceTask(pipeBase.Task):
    """Task for measuring and storing forced sources at DiaObject locations
    in both difference and direct images.
    """
    ConfigClass = DiaForcedSourcedConfig
    _DefaultName = "diaForcedSource"

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.dia_forced_source_schema = make_dia_forced_source_schema()
        self.makeSubtask("forcedMeasurement",
                         refSchema=make_dia_object_schema())

    def run(self, dia_objects, expIdBits, exposure, diffim):
        """Measure forced sources on the direct and different images,
        calibrate, and store them in the Ppdb.

        Parameters
        ----------
        dia_objects : `lsst.afw.table.SourceCatalog`
            Catalog of previously observed and newly created DiaObjects
            contained within the difference and direct images.
        expIdBits : `int`
            Bit length of the exposure id.
        exposure : `lsst.afw.image.Exposure`
            Direct image exposure.
        diffim : `lsst.afw.image.Exposure`
            Difference image.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components.

            - ``diffForcedSources`` : Psf forced photometry on the difference
              image at the DiaObject location (`lsst.afw.table.SourceCatalog`)
            - ``directForcedSources`` : Psf forced photometry on the direct
              image at the DiaObject location (`lsst.afw.table.SourceCatalog`)
        """

        idFactoryDiff = afwTable.IdFactory.makeSource(
            diffim.getInfo().getVisitInfo().getExposureId(),
            afwTable.IdFactory.computeReservedFromMaxBits(int(expIdBits)))
        idFactoryDirect = afwTable.IdFactory.makeSource(
            diffim.getInfo().getVisitInfo().getExposureId(),
            afwTable.IdFactory.computeReservedFromMaxBits(int(expIdBits)))

        diffForcedSources = self.forcedMeasurement.generateMeasCat(
            diffim,
            dia_objects,
            diffim.getWcs(),
            idFactory=idFactoryDiff)
        self.forcedMeasurement.run(
            diffForcedSources, diffim, dia_objects, diffim.getWcs())

        directForcedSources = self.forcedMeasurement.generateMeasCat(
            exposure,
            dia_objects,
            exposure.getWcs(),
            idFactory=idFactoryDirect)
        self.forcedMeasurement.run(
            directForcedSources, exposure, dia_objects, exposure.getWcs())

        return pipeBase.Struct(
            diffForcedSources=diffForcedSources,
            directForcedSources=directForcedSources,
        )
