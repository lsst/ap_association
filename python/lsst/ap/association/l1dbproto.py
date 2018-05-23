#
# LSST Data Management System
# Copyright 2017 LSST/AURA.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

"""Simple sqlite3 interface for storing and retrieving DIAObjects and
DIASources. This class is mainly for testing purposes in the context of
ap_pipe/ap_verify.
"""
import os

from lsst.meas.algorithms.indexerRegistry import IndexerRegistry
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
from lsst.ap.association import \
    make_minimal_dia_object_schema, \
    make_minimal_dia_source_schema
from .afwUtils import \
    convert_dia_source_to_asssoc_schema
import lsst.l1dbproto.l1db as l1db
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

__all__ = ["AssociationL1DBProtoConfig",
           "AssociationL1DBProtoTask"]


def _cfg_file_name(basename):
    """Return path name of a data file.
    """
    return os.path.join(os.environ.get("L1DBPROTO_DIR"), "cfg", basename)


class AssociationL1DBProtoConfig(pexConfig.Config):
    """Configuration parameters for the AssociationL1DBProtoTask.
    """
    indexer = IndexerRegistry.makeField(
        doc='Select the spatial indexer to use within the database.',
        default='HTM'
    )
    db_config_file = pexConfig.Field(
        dtype=str,
        doc='',
        default=_cfg_file_name("l1db-ap-pipe-sqlite.py")
    )


class AssociationL1DBProtoTask(pipeBase.Task):
    """
    """

    ConfigClass = AssociationL1DBProtoConfig
    _DefaultName = "association_l1dbroto_task"

    def __init__(self, **kwargs):

        pipeBase.Task.__init__(self, **kwargs)
        self.indexer = IndexerRegistry[self.config.indexer.name](
            self.config.indexer.active)
        self.l1db_config = l1db.L1dbConfig()
        self.l1db_config.load(self.config.db_config_file)

        self._dia_object_afw_schema = make_minimal_dia_object_schema(
            ['u', 'g', 'r', 'i', 'z', 'y'])
        self._dia_source_afw_schema = make_minimal_dia_source_schema()

        afw_schema = dict(
            DiaObject=self._dia_object_afw_schema,
            DiaSource=self._dia_source_afw_schema)
        self.db = l1db.L1db(config=self.l1db_config, afw_schemas=afw_schema)

    def load_dia_objects(self, exposure):
        """Load all DIAObjects within the exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            An exposure with a solved WCS representing the area on the sky to
            load DIAObjects.

        Returns
        -------
        dia_objects : `lsst.afw.table.SourceCatalog`
            Catalog of DIAObjects that are contained with the the bounding
            box defined by expMd.
        """
        bbox = afwGeom.Box2D(exposure.getBBox())
        wcs = exposure.getWcs()
        expMd = pipeBase.Struct(
            bbox=bbox,
            wcs=wcs,)

        ctr_coord = expMd.wcs.pixelToSky(expMd.bbox.getCenter())
        max_radius = max(
            ctr_coord.separation(expMd.wcs.pixelToSky(pp))
            for pp in expMd.bbox.getCorners())

        indexer_indices, on_boundry = self.indexer.get_pixel_ids(
            ctr_coord, max_radius)
        # Index types must be cast to int to work with l1dbproto.
        index_ranges = [[int(idexer_idx), int(idexer_idx) + 1] for idexer_idx in indexer_indices]
        covering_dia_objects = self.db.getDiaObjects(index_ranges)

        output_dia_objects = afwTable.SourceCatalog(
            covering_dia_objects.getSchema())
        for cov_dia_object in covering_dia_objects:
            if self._check_dia_object_position(cov_dia_object, expMd):
                output_dia_objects.append(cov_dia_object)

        return output_dia_objects.copy(deep=True)

    def _check_dia_object_position(self, dia_object_record, expMd):
        """Check the RA, DEC position of the current dia_object_record against
        the bounding box of the exposure.

        Parameters
        ----------
        dia_object_record : `lsst.afw.table.SourceRecord`
            A SourceRecord object containing the DIAObject we would like to
            test against our bounding box.
        expMd : `lsst.pipe.base.Struct` (optional)
            Results struct with components:

            - ``bbox``: Bounding box of exposure (`lsst.afw.geom.Box2D`).
            - ``wcs``: WCS of exposure (`lsst.afw.geom.SkyWcs`).

        Return
        ------
        is_contained : `bool`
            Object position is contained within the bounding box of expMd.
        """
        if expMd is None:
            return True
        point = expMd.wcs.skyToPixel(dia_object_record.getCoord())
        return expMd.bbox.contains(point)

    def load_dia_sources(self, dia_obj_ids):
        """Retrieve all DIASources associated with this collection of DIAObject
        ids.

        Parameters
        ----------
        dia_obj_ids : array-like of `int`s
            Id of the DIAObject that is associated with the DIASources
            of interest.

        Returns
        -------
        dia_sources : `lsst.afw.table.SourceCatalog`
            SourceCatalog of DIASources
        """

        # l1db proto does not currently use the dt (dateTime) variables for
        # this function.
        return self.db.getDiaSources(dia_obj_ids, None)

    def store_dia_objects(self,
                          dia_objects,
                          compute_spatial_index=False,
                          exposure=None):
        """Store all DIAObjects in this SourceCatalog.

        Parameters
        ----------
        dia_objects : `lsst.afw.table.SourceCatalog`
            Catalog of DIAObjects to store.
        compute_spatial_index : `bool`
            If True, compute the spatial search indices using the
            indexer specified at class instantiation.
        exposure: `lsst.afw.image.Exposure` (optional)
            CcdExposure associated with these DIAObjects being inserted.
            Inserts the CcdVisitInfo for this exposure in the CcdVisitTable.
        """
        if compute_spatial_index:
            for dia_object in dia_objects:
                pixelId = self.compute_indexer_id(dia_object.getCoord())
                dia_object['pixelId'] = pixelId

        dt = exposure.getInfo().getVisitInfo().getDate().toPython()
        self.db.storeDiaObjects(dia_objects, dt)

    def compute_indexer_id(self, sphere_point):
        """Compute the pixel index of the given point.

        Parameters
        ----------
        sphere_point : `lsst.afw.geom.SpherePoint`
            Point to compute pixel index for.

        Returns
        -------
        index : `int`
            Index of the pixel the point is contained in.
        """
        return self.indexer.index_points(
            [sphere_point.getRa().asDegrees()],
            [sphere_point.getDec().asDegrees()])[0]

    def store_dia_sources(self,
                          dia_sources,
                          associated_ids=None,
                          exposure=None):
        """Store all DIASources in this SourceCatalog.

        Parameters
        ----------
        dia_sources : `lsst.afw.table.SourceCatalog`
            Catalog of DIASources to store.
        associated_ids : array-like of `int`s (optional)
            DIAObject ids that have been associated with these DIASources
        exposure : `lsst.afw.image.Exposure`
            Exposure object the DIASources were detected in.
        """
        if dia_sources.getSchema().contains(make_minimal_dia_source_schema()) and \
           associated_ids is None and exposure is None:
                self.db.storeDiaSources(dia_sources)
        else:
            dia_sources_to_store = \
                convert_dia_source_to_asssoc_schema(
                    dia_sources, associated_ids, exposure)

            self.db.storeDiaSources(dia_sources_to_store)

    def get_dia_object_schema(self):
        """Retrieve the Schema of the DIAObjects in this database.

        Returns
        -------
        schema : `lsst.afw.table.Schema`
            Schema of the DIAObjects in this database.
        """
        return self._dia_object_afw_schema

    def get_dia_source_schema(self):
        """Retrieve the Schema of the DIASources in this database.

        Returns
        -------
        schema : `lsst.afw.table.Schema`
            Schema of the DIASources in this database.
        """
        return self._dia_source_afw_schema
