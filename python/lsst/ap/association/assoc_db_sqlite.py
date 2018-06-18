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

from __future__ import absolute_import, division, print_function

import numpy as np
import sqlite3

from lsst.meas.algorithms.indexerRegistry import IndexerRegistry
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from .afwUtils import \
    make_minimal_dia_object_schema, \
    make_minimal_dia_source_schema, \
    getCcdVisitSchemaSql, \
    add_dia_source_aliases_to_catalog, \
    get_ccd_visit_info_from_exposure, \
    make_overwrite_dict

__all__ = ["AssociationDBSqliteConfig",
           "AssociationDBSqliteTask",
           "SqliteDBConverter"]


class SqliteDBConverter(object):
    """Class for defining conversions to and from an sqlite database and afw
    SourceRecord objects.

    Parameters
    ----------
    schema : `lsst.afw.table.Schema`
        Schema defining the SourceRecord objects to be converted.
    table_name : `str`
        Name of the sqlite table this converter is to be used for.
    """

    def __init__(self, schema, table_name):
        self._schema = schema
        self._table_name = table_name
        self._afw_to_db_types = {
            "Angle": "REAL",
            "D": "REAL",
            "L": "INTEGER",
            "String": "TEXT",
        }

    @property
    def table_name(self):
        """Return name of the sqlite table this catalog is for
        """
        return self._table_name

    @property
    def schema(self):
        """Return the internal catalog schema.
        """
        return self._schema

    def make_table_from_afw_schema(self, table_name):
        """Convert the schema into a sqlite CREATE TABLE command.

        Parameters
        ----------
        table_name : `str`
            Name of the new table to create

        Returns
        -------
        sql_query : `str`
            A string of the query command to create the new table in sqlite.
        """
        name_type_string = ""
        for sub_schema in self._schema:
            tmp_name = sub_schema.getField().getName()
            tmp_type = self._afw_to_db_types[
                sub_schema.getField().getTypeString()]
            if tmp_name == 'id':
                tmp_type += " PRIMARY KEY"
            name_type_string += "%s %s," % (tmp_name, tmp_type)
        name_type_string = name_type_string[:-1]

        return "CREATE TABLE %s (%s)" % (table_name, name_type_string)

    def source_record_from_db_row(self, db_row):
        """Create a source record from the values stored in a database row.

        Parameters
        ----------
        db_row : `list` of ``values``
            Retrieved values from the database to convert into a SourceRecord.

        Returns
        -------
        record : `lsst.afw.table.SourceRecord`
            Converted source record.
        """

        output_source_record = afwTable.SourceTable.makeRecord(
            afwTable.SourceTable.make(self._schema))

        for sub_schema, value in zip(self._schema, db_row):
            if value is None:
                value = np.nan
            if sub_schema.getField().getTypeString() == 'Angle':
                output_source_record.set(
                    sub_schema.getKey(), value * afwGeom.degrees)
            else:
                output_source_record.set(sub_schema.getKey(), value)
        return output_source_record

    def source_record_to_value_list(self, source_record, overwrite_dict={}):
        """Convert a source record object into a list of its internal values.

        Parameters
        ----------
        source_record : `lsst.afw.table.SourceRecord`
            SourceRecord to convert.
        overwrite_dict : `dict` (optional)
            Mapping specifying the names of columns to overwrite with
            specified values.

        Returns
        -------
        source_list : `list` of ``values``
            Extracted values from ``source_record`` in `list` form.
        """
        values = []
        for sub_schema in self._schema:
            field_name = sub_schema.getField().getName()
            if field_name in overwrite_dict:
                values.append(overwrite_dict[field_name])
            else:
                if sub_schema.getField().getTypeString() == 'Angle':
                    values.append(
                        source_record.get(sub_schema.getKey()).asDegrees())
                else:
                    values.append(source_record.get(sub_schema.getKey()))
        return values


class AssociationDBSqliteConfig(pexConfig.Config):
    """Configuration parameters for the AssociationDBSqliteTask
    """
    db_name = pexConfig.Field(
        dtype=str,
        doc='Location on disk and name of the sqlite3 database for storing '
        'and loading DIASources and DIAObjects.',
        default=':memory:'
    )
    filter_names = pexConfig.ListField(
        dtype=str,
        doc='List of filter names to store and expect from in this DB.',
        default=[],
    )
    indexer = IndexerRegistry.makeField(
        doc='Select the spatial indexer to use within the database.',
        default='HTM'
    )


class AssociationDBSqliteTask(pipeBase.Task):
    """Enable storage of and reading of DIAObjects and DIASources from a
    sqlite database.

    Create a simple sqlite database and implement wrappers to store and
    retrieve DIAObjects and DIASources from within that database. This task
    functions as a testing ground for the L1 database and should mimic this
    database's eventual functionality. This specific database implementation is
    useful for the verification packages which may not be run with access to
    L1 database.
    """

    ConfigClass = AssociationDBSqliteConfig
    _DefaultName = "association_db_sqlite"

    def __init__(self, **kwargs):

        pipeBase.Task.__init__(self, **kwargs)
        self.indexer = IndexerRegistry[self.config.indexer.name](
            self.config.indexer.active)
        self._db_connection = sqlite3.connect(self.config.db_name)
        self._db_cursor = self._db_connection.cursor()

        self._dia_object_converter = SqliteDBConverter(
            make_minimal_dia_object_schema(self.config.filter_names),
            "dia_objects")
        self._dia_source_converter = SqliteDBConverter(
            make_minimal_dia_source_schema(),
            "dia_sources")
        self._ccd_visit_schema = getCcdVisitSchemaSql()

    def _commit(self):
        """Save changes to the sqlite database.
        """
        self._db_connection.commit()

    def close(self):
        """Close the connection to the sqlite database.
        """
        self._db_connection.close()

    def create_tables(self):
        """If no sqlite database with the correct tables exists we can create
        one using this method.

        Returns
        -------
        succeeded : `bool`
            Successfully created a new database with specified tables.
        """

        self._db_cursor.execute(
            'select name from sqlite_master where type = "table"')
        db_tables = self._db_cursor.fetchall()

        # If this database currently contains any tables exit and do not
        # create tables.
        if db_tables:
            return False
        else:
            # Create tables to store the individual DIAObjects and DIASources
            self._db_cursor.execute(
                self._dia_object_converter.make_table_from_afw_schema(
                    "dia_objects"))
            self._db_cursor.execute(
                "CREATE INDEX indexer_id_index ON dia_objects(pixelId)")
            self._commit()
            self._db_cursor.execute(
                self._dia_source_converter.make_table_from_afw_schema(
                    "dia_sources"))
            self._db_cursor.execute(
                "CREATE INDEX diaObjectId_index ON dia_sources(diaObjectId)")
            self._commit()

            table_schema = ",".join(
                "%s %s" % (key, self._ccd_visit_schema[key])
                for key in self._ccd_visit_schema.keys())
            self._db_cursor.execute(
                "CREATE TABLE CcdVisit (%s)" % table_schema)
            self._commit()

        return True

    @pipeBase.timeMethod
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

        dia_objects = self._get_dia_object_catalog(indexer_indices, expMd)

        return dia_objects

    @pipeBase.timeMethod
    def load_dia_sources(self, dia_obj_ids):
        """Retrieve all DIASources associated with this DIAObject id.

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

        dia_source_schema = make_minimal_dia_source_schema()
        output_dia_sources = afwTable.SourceCatalog(dia_source_schema)

        rows = self._query_dia_sources(dia_obj_ids)
        for row in rows:
            output_dia_sources.append(
                self._dia_source_converter.source_record_from_db_row(row))

        return output_dia_sources.copy(deep=True)

    @pipeBase.timeMethod
    def store_dia_objects(self, dia_objects, compute_spatial_index=False):
        """Store all DIAObjects in this SourceCatalog.

        Parameters
        ----------
        dia_objects : `lsst.afw.table.SourceCatalog`
            Catalog of DIAObjects to store.
        compute_spatial_index : `bool`
            If True, compute the spatial search indices using the
            indexer specified at class instantiation.
        """
        if compute_spatial_index:
            for dia_object in dia_objects:
                sphPoint = dia_object.getCoord()
                dia_object.set('pixelId',
                               self.compute_indexer_id(sphPoint))
        self._store_catalog(dia_objects, self._dia_object_converter, None)
        self._commit()

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

    @pipeBase.timeMethod
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
        self._store_catalog(dia_sources,
                            self._dia_source_converter,
                            associated_ids,
                            exposure)
        self._commit()

    def store_ccd_visit_info(self, exposure):
        """Store information describing the exposure for this ccd, visit.

        Paramters
        ---------
        exposure : `lsst.afw.image.Exposure`
            Exposure to store information from.
        """

        values = get_ccd_visit_info_from_exposure(exposure)
        self._db_cursor.execute(
            "INSERT OR REPLACE INTO CcdVisit VALUES (%s)" %
            ",".join("?" for idx in range(len(self._ccd_visit_schema))),
            [values[key] for key in self._ccd_visit_schema.keys()])

    def _get_dia_object_catalog(self, indexer_indices, expMd=None):
        """Retrieve the DIAObjects from the database whose indexer indices
        are with the specified list of indices.

        Retrieves a list of DIAObjects that are covered by the pixels with
        indices, indexer_indices. Use this to retrieve complete DIAObjects.

        Parameters
        ----------
        indexer_indices : array-like of `int`s
            Pixelized indexer indices from which to load.
        expMd : `lsst.pipe.base.Struct` (optional)
            Results struct with components:

            - ``bbox``: Bounding box of exposure (`lsst.afw.geom.Box2D`).
            - ``wcs``: WCS of exposure (`lsst.afw.geom.SkyWcs`).

        Returns
        -------
        dia_objects : `lsst.afw.table.SourceCatalog`
            Catalog of DIAObjects with the specified indexer index and
            contained within the expMd bounding box.
        """
        dia_object_rows = self._query_dia_objects(indexer_indices)

        output_dia_objects = afwTable.SourceCatalog(
            self._dia_object_converter.schema)

        for row in dia_object_rows:
            dia_object_record = \
                self._dia_object_converter.source_record_from_db_row(row)
            if self._check_dia_object_position(dia_object_record, expMd):
                output_dia_objects.append(dia_object_record)

        return output_dia_objects.copy(deep=True)

    def _query_dia_objects(self, indexer_indices):
        """Query the database for the stored DIAObjects given a set of
        indices in the indexer.

        Parameters
        ----------
        indexer_indices : array-like of `int`s
            Spatial indices in the indexer specifying the area on the sky
            to load DIAObjects for.

        Returns
        -------
        dia_objects : `list` of `tuples` containing ``dia_object`` ``values``
            Query result containing the catalog values of the DIAObject.
        """
        self._db_cursor.execute(
            "CREATE TEMPORARY TABLE tmp_indexer_indices "
            "(pixelId INTEGER PRIMARY KEY)")

        self._db_cursor.executemany(
            "INSERT OR REPLACE INTO tmp_indexer_indices VALUES (?)",
            [(int(indexer_index),) for indexer_index in indexer_indices])

        self._db_cursor.execute(
            "SELECT o.* FROM dia_objects AS o "
            "INNER JOIN tmp_indexer_indices AS i "
            "ON o.pixelId = i.pixelId")

        output_rows = self._db_cursor.fetchall()

        self._db_cursor.execute("DROP TABLE tmp_indexer_indices")
        self._commit()

        return output_rows

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

    def _query_dia_sources(self, dia_object_ids):
        """Query the database for the stored DIASources given a set of
        DIAObject ids.

        Parameters
        ----------
        dia_object_ids : array-like of `int`s
            Spatial indices in the indexer specifying the area on the sky
            to load DIAObjects for.

        Return
        ------
        dia_objects : `list` of `tuples`
            Query result containing the values representing DIASources
        """
        self._db_cursor.execute(
            "CREATE TEMPORARY TABLE tmp_object_ids "
            "(diaObjectId INTEGER PRIMARY KEY)")

        self._db_cursor.executemany(
            "INSERT OR REPLACE INTO tmp_object_ids VALUES (?)",
            [(int(dia_object_id),) for dia_object_id in dia_object_ids])

        self._db_cursor.execute(
            "SELECT s.* FROM dia_sources AS s "
            "INNER JOIN tmp_object_ids AS o "
            "ON s.diaObjectId = o.diaObjectId")

        output_rows = self._db_cursor.fetchall()

        self._db_cursor.execute("DROP TABLE tmp_object_ids")
        self._commit()

        return output_rows

    def _store_catalog(self, source_catalog, converter, obj_ids=None, exposure=None):
        """ Store a SourceCatalog into the database.

        Parameters
        ----------
        source_catalog : `lsst.afw.table.SourceCatalog`
            SourceCatalog to store in the database table specified by
            converter.
        converter : `lsst.ap.association.SqliteDBConverter`
            A converter object specifying the correct database table to write
            into.
        obj_id : array-like of `int`s (optional)
            Ids of the DIAObjects these objects are associated with. Use only
            when storing DIASources.
        exposure : `lsst.afw.image.Exposure` (optional)
            Exposure that the sources in source_catalog were detected in. If
            set the fluxes are calibrated and stored using the exposure Calib
            object. The filter the exposure was taken in as well as the
            ccdVisitId are also stored.
        """
        values = []

        if converter.table_name == 'dia_sources':
            # Create aliases to appropriate flux fields if they exist.
            add_dia_source_aliases_to_catalog(source_catalog)

        if exposure is None:
            exp_dict = None
        else:
            exp_dict = get_ccd_visit_info_from_exposure(exposure)

        for src_idx, source_record in enumerate(source_catalog):
            if obj_ids is None:
                obj_id = None
            else:
                obj_id = obj_ids[src_idx]

            overwrite_dict = make_overwrite_dict(source_record,
                                                 obj_id,
                                                 exp_dict)

            values.append(converter.source_record_to_value_list(
                source_record, overwrite_dict))

        insert_string = ("?," * len(values[0]))[:-1]

        self._db_cursor.executemany(
            "INSERT OR REPLACE INTO %s VALUES (%s)" %
            (converter.table_name, insert_string), values)

    def get_dia_object_schema(self):
        """Retrieve the Schema of the DIAObjects in this database.

        Returns
        -------
        schema : `lsst.afw.table.Schema`
            Schema of the DIAObjects in this database.
        """
        return self._dia_object_converter.schema

    def get_dia_source_schema(self):
        """Retrieve the Schema of the DIASources in this database.

        Returns
        -------
        schema : `lsst.afw.table.Schema`
            Schema of the DIASources in this database.
        """
        return self._dia_source_converter.schema
