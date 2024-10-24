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

import numpy as np
import pandas as pd

import lsst.afw.table as afwTable
from lsst.daf.base import DateTime
import lsst.geom as geom
from lsst.meas.base import ForcedMeasurementTask
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.utils.timer import timeMethod


class DiaForcedSourcedConfig(pexConfig.Config):
    """Configuration for the generic DiaForcedSourcedTask class.
    """
    forcedMeasurement = pexConfig.ConfigurableField(
        target=ForcedMeasurementTask,
        doc="Subtask to force photometer DiaObjects in the direct and "
            "difference images.",
    )
    dropColumns = pexConfig.ListField(
        dtype=str,
        doc="Columns produced in forced measurement that can be dropped upon "
            "creation and storage of the final pandas data.",
    )
    historyThreshold = pexConfig.Field(
        dtype=int,
        doc="Minimum number of detections of a diaObject required "
            "to run forced photometry. Set to 1 to include all diaObjects.",
        default=2,
    )

    def setDefaults(self):
        self.forcedMeasurement.plugins = ["base_TransformedCentroidFromCoord",
                                          "base_PsfFlux"]
        self.forcedMeasurement.doReplaceWithNoise = False
        self.forcedMeasurement.copyColumns = {
            "id": "diaObjectId",
            "coord_ra": "ra",
            "coord_dec": "dec"}
        self.forcedMeasurement.slots.centroid = "base_TransformedCentroidFromCoord"
        self.forcedMeasurement.slots.psfFlux = "base_PsfFlux"
        self.forcedMeasurement.slots.shape = None
        self.dropColumns = ['coord_ra', 'coord_dec', 'x', 'y', 'parent',
                            'base_TransformedCentroidFromCoord_x',
                            'base_TransformedCentroidFromCoord_y',
                            'base_PsfFlux_instFlux',
                            'base_PsfFlux_instFluxErr', 'base_PsfFlux_area',
                            'slot_PsfFlux_area', 'base_PsfFlux_flag',
                            'slot_PsfFlux_flag',
                            'base_PsfFlux_flag_noGoodPixels',
                            'slot_PsfFlux_flag_noGoodPixels',
                            'base_PsfFlux_flag_edge', 'slot_PsfFlux_flag_edge',
                            'base_PsfFlux_chi2', 'slot_PsfFlux_chi2',
                            'base_PsfFlux_npixels', 'slot_PsfFlux_npixels',
                            ]


class DiaForcedSourceTask(pipeBase.Task):
    """Task for measuring and storing forced sources at DiaObject locations
    in both difference and direct images.
    """
    ConfigClass = DiaForcedSourcedConfig
    _DefaultName = "diaForcedSource"

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.makeSubtask("forcedMeasurement",
                         refSchema=afwTable.SourceTable.makeMinimalSchema())

    @timeMethod
    def run(self,
            dia_objects,
            updatedDiaObjectIds,
            exposure,
            diffim,
            idGenerator):
        """Measure forced sources on the direct and difference images.

        Parameters
        ----------
        dia_objects : `pandas.DataFrame`
            Catalog of previously observed and newly created DiaObjects
            contained within the difference and direct images. DiaObjects
            must be indexed on the ``diaObjectId`` column.
        updatedDiaObjectIds : `numpy.ndarray`
            Array of diaObjectIds that were updated during this dia processing.
            Used to assure that the pipeline includes all diaObjects that were
            updated in case one falls on the edge of the CCD.
        exposure : `lsst.afw.image.Exposure`
            Direct image exposure.
        diffim : `lsst.afw.image.Exposure`
            Difference image.
        idGenerator : `lsst.meas.base.IdGenerator`
            Object that generates source IDs and random number generator seeds.

        Returns
        -------
        output_forced_sources : `pandas.DataFrame`
            Catalog of calibrated forced photometered fluxes on both the
            difference and direct images at DiaObject locations.
        """
        # Restrict forced source measurement to objects with sufficient history to be reliable.
        objectTable = dia_objects.query(f'nDiaSources >= {self.config.historyThreshold}')
        if objectTable.empty:
            # The dataframe will be coerced to the correct (empty) format in diaPipe.
            return pd.DataFrame()

        afw_dia_objects = self._convert_from_pandas(objectTable)

        idFactoryDiff = idGenerator.make_table_id_factory()

        diffForcedSources = self.forcedMeasurement.generateMeasCat(
            diffim,
            afw_dia_objects,
            diffim.getWcs(),
            idFactory=idFactoryDiff)
        self.forcedMeasurement.run(
            diffForcedSources, diffim, afw_dia_objects, diffim.getWcs())

        directForcedSources = self.forcedMeasurement.generateMeasCat(
            exposure,
            afw_dia_objects,
            exposure.getWcs(),
            idFactory=idFactoryDiff)
        self.forcedMeasurement.run(
            directForcedSources, exposure, afw_dia_objects, exposure.getWcs())

        output_forced_sources = self._calibrate_and_merge(diffForcedSources,
                                                          directForcedSources,
                                                          diffim,
                                                          exposure)

        output_forced_sources = self._trim_to_exposure(output_forced_sources,
                                                       updatedDiaObjectIds,
                                                       exposure)
        # Drop superfluous columns from output DataFrame.
        output_forced_sources.drop(columns=self.config.dropColumns, inplace=True)
        return output_forced_sources.set_index(
            ["diaObjectId", "diaForcedSourceId"],
            drop=False)

    def _convert_from_pandas(self, input_objects):
        """Create minimal schema SourceCatalog from a pandas DataFrame.

        We need a catalog of this type to run within the forced measurement
        subtask.

        Parameters
        ----------
        input_objects : `pandas.DataFrame`
            DiaObjects with locations and ids. ``

        Returns
        -------
        outputCatalog : `lsst.afw.table.SourceTable`
            Output catalog with minimal schema.
        """
        schema = afwTable.SourceTable.makeMinimalSchema()

        outputCatalog = afwTable.SourceCatalog(schema)
        outputCatalog.reserve(len(input_objects))

        for obj_id, df_row in input_objects.iterrows():
            outputRecord = outputCatalog.addNew()
            outputRecord.setId(obj_id)
            outputRecord.setCoord(
                geom.SpherePoint(df_row["ra"],
                                 df_row["dec"],
                                 geom.degrees))
        return outputCatalog

    def _calibrate_and_merge(self,
                             diff_sources,
                             direct_sources,
                             diff_exp,
                             direct_exp):
        """Take the two output catalogs from the ForcedMeasurementTasks and
        calibrate, combine, and convert them to Pandas.

        Parameters
        ----------
        diff_sources : `lsst.afw.table.SourceTable`
            Catalog with PsFluxes measured on the difference image.
        direct_sources : `lsst.afw.table.SourceTable`
            Catalog with PsfFluxes measured on the direct (calexp) image.
        diff_exp : `lsst.afw.image.Exposure`
            Difference exposure ``diff_sources`` were measured on.
        direct_exp : `lsst.afw.image.Exposure`
            Direct (calexp) exposure ``direct_sources`` were measured on.

        Returns
        -------
        output_catalog : `pandas.DataFrame`
            Catalog calibrated diaForcedSources.
        """
        diff_calib = diff_exp.getPhotoCalib()
        direct_calib = direct_exp.getPhotoCalib()

        diff_fluxes = diff_calib.instFluxToNanojansky(diff_sources,
                                                      "slot_PsfFlux")
        direct_fluxes = direct_calib.instFluxToNanojansky(direct_sources,
                                                          "slot_PsfFlux")

        output_catalog = diff_sources.asAstropy().to_pandas()
        # afwTable source catalogs store coordinates as radians, but the
        # output must be in degrees
        output_catalog.loc[:, "ra"] = np.rad2deg(output_catalog.loc[:, "ra"])
        output_catalog.loc[:, "dec"] = np.rad2deg(output_catalog.loc[:, "dec"])
        output_catalog.rename(columns={"id": "diaForcedSourceId",
                                       "slot_PsfFlux_instFlux": "psfFlux",
                                       "slot_PsfFlux_instFluxErr": "psfFluxErr",
                                       "slot_Centroid_x": "x",
                                       "slot_Centroid_y": "y"},
                              inplace=True)
        output_catalog.loc[:, "psfFlux"] = diff_fluxes[:, 0]
        output_catalog.loc[:, "psfFluxErr"] = diff_fluxes[:, 1]

        output_catalog["scienceFlux"] = direct_fluxes[:, 0]
        output_catalog["scienceFluxErr"] = direct_fluxes[:, 1]

        midpointMjdTai = direct_exp.visitInfo.date.get(system=DateTime.MJD)
        output_catalog["visit"] = direct_exp.visitInfo.id
        output_catalog["detector"] = direct_exp.detector.getId()
        output_catalog["midpointMjdTai"] = midpointMjdTai
        output_catalog["band"] = diff_exp.getFilter().bandLabel
        output_catalog["time_processed"] = DateTime.now().toPython()
        # TODO: propagate actual flags (DM-42355)

        return output_catalog

    def _trim_to_exposure(self, catalog, updatedDiaObjectIds, exposure):
        """Remove DiaForcedSources that are outside of the bounding box region.

        Paramters
        ---------
        catalog : `pandas.DataFrame`
            DiaForcedSources to check against the exposure bounding box.
        updatedDiaObjectIds : `numpy.ndarray`
            Array of diaObjectIds that were updated during this dia processing.
            Used to assure that the pipeline includes all diaObjects that were
            updated in case one falls on the edge of the CCD.
        exposure : `lsst.afw.image.Exposure`
            Exposure to check against.

        Returns
        -------
        output : `pandas.DataFrame`
            DataFrame trimmed to only the objects within the exposure bounding
            box.
        """
        bbox = geom.Box2D(exposure.getBBox())

        xS = catalog.loc[:, "x"]
        yS = catalog.loc[:, "y"]

        return catalog[
            np.logical_or(bbox.contains(xS, yS),
                          np.isin(catalog.loc[:, "diaObjectId"],
                                  updatedDiaObjectIds))]
