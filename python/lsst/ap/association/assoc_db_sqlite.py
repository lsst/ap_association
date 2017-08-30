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

from __future__ import absolute_import, division, print_function

import sqlite3

from lsst.meas.algorithms.indexerRegistry import IndexerRegistry
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from .dia_collection import DIAObjectCollection
from .dia_object import DIAObject, \
                        make_minimal_dia_object_schema, \
                        make_minimal_dia_source_schema


__all__ = ["AssociationDBSqliteConfig", "AssociationDBSqliteTask"]


afw_to_db_types = {
    "L": "INTEGER",
    "Angle": "REAL"
}


class AssociationDBSqliteConfig(pexConfig.Config):
    """!
    \anchor AssociationDBSqliteConfig_

    \brief Configuration parameters for the AssociationDBSqliteTask
    """
    db_name = pexConfig.Field(
        dtype=str,
        doc='Location on disk and name of the sqlite3 database for storing '
        'and loading DIASources and DIAObjects.',
        default=':memory:'
    )
    indexer = IndexerRegistry.makeField(
        doc='Select the spatial indexer to use within the database.',
        default='HTM'
    )


class AssociationDBSqliteTask(pipeBase.Task):
    """!
    \anchor AssociationDBSqliteTask_

    \brief Enable storage of reading of DIAObjects and DIASources from a
    sqlite database.

    \section ap_association_AssociationDBSqliteTask_Contents Contents
      - \ref ap_association_AssociationDBSqliteTask_Purpose
      - \ref ap_association_AssociationDBSqliteTask_Config
      - \ref ap_association_AssociationDBSqliteTask_Run
      - \ref ap_association_AssociationDBSqliteTask_Debug
      - \ref ap_association_AssociationDBSqliteTask_Example

    \section ap_association_AssociationDBSqliteTask_Purpose   Description
    Create a simple sqlite database and implement wrappers to store and
    retrieve DIAObjects and DIASources within that database. This task
    functions as a testing ground for the L1 database and should mimic this
    database's eventual functionaly. This specific database implementation is
    useful for the verification packages which may not be run with access to
    L1 database.

    \section ap_association_AssociationDBSqliteTaskInitialize
        Task initialization
    \copydoc \_\_init\_\_

    \section ap_association_AssociationDBSqliteTask_Run
        Invoking the Task
    \copydoc run

    \section ap_association_AssociationDBSqliteTask_Config
        Configuration parameters
    See \ref AssociationDBSqliteConfig

    \section ap_association_AssociationDBSqliteTask_Debug
        Debug variables
    This task has no debug variables

    \section ap_association_AssociationDBSqliteTask_Example
        Example of using AssociationDBSqliteTask
    This task has no standalone example, however it is applied as a subtask of
    ap.association.AssociationTask
    """

    ConfigClass = AssociationDBSqliteConfig
    _DefaultName = "association_db_sqlite"

    def __init__(self, **kwargs):

        pipeBase.Task.__init__(self, **kwargs)
        self.indexer = IndexerRegistry[self.config.indexer.name](
            self.config.indexer.active)
        self._db_connection = sqlite3.connect(self.config.db_name)
        self._db_cursor = self._db_connection.cursor()

    def run(self):
        """ This subtask does not make use of the run method for Tasks.

        Please use the store or load methods.
        """
        return None

    def _commit(self):
        """ Save changes to the sqlite database.
        """
        self._db_connection.commit()

    def close(self):
        """ Close the connection to the sqlite database.
        """
        self._db_connection.close()

    def create_tables(self):
        """ If no sqlite database exists with the correct tables we can create
        one using this method.

        Returns
        -------
        bool
            Successfully created a new database with specified tables.
        """

        obj_schema = make_minimal_dia_object_schema()
        src_schema = make_minimal_dia_source_schema()

        self._db_cursor.execute(
            'select name from sqlite_master where type = "table"')
        db_tables = self._db_cursor.fetchall()

        if db_tables:
            return False
        else:
            # Create databases to store the individual DIAObjects and DIASources
            self._table_from_afw_shcema('dia_objects', obj_schema)
            self._table_from_afw_shcema('dia_sources', src_schema)

            # Create linkage database between associated dia_objects and
            # dia_sources.
            self._db_cursor.execute(
                "CREATE TABLE dia_objects_to_dia_sources ("
                "FOREIGN KEY(src_id) REFERENCES dia_sources(id), "
                "FOREIGN KEY(obj_id) REFERENCES dia_objects(id)"
                ")")
            self._db_connection.commit()

        return True

    def _table_from_afw_shcema(table_name, schema):
        """ Create a new table from an afw.table.Schema

        The primary key of the table is assumed to have the name id.

        Parameters
        ----------
        table_name : str
            Name of the table to create
        schema : afw.table.Schema
            Schema of the objects that his database will store.
        """
        name_type_string = ""
        for sub_schema in obj_schema:
            tmp_name = sub_schema.getField().getName()
            tmp_type = afw_to_db_types[
                sub_schema.getField().getTypeString()]
            if tmp_name == 'id':
                tmp_type += " PRIMARY KEY"
            name_type_string += "%s %s," % (tmp_name, tmp_type)
        name_type_string = name_type_string[:-1]
        self._db_cursor.execute(
            "CREATE TABLE %s (%s)" % (table_name, name_type_string))
        self._db_connection.commit()

    @pipeBase.timeMethod
    def load(self, ctr_coord, radius):
        """ Load all DIAObjects and associated DIASources.

        Parameters
        ----------
        ctr_coord : lsst.afw.geom.SpherePoint
            Center position of the circle on the sky to load.
        radius : lsst.afw.geom.Angle


        Returns
        -------
        lsst.ap.association.DIAObjectCollection
        """
        indexer_indices, on_boundry = self.indexer.get_pixel_ids(
            ctr_coord, radius)

        dia_objects = self.get_dia_objects(indexer_indices)

        dia_collection = DIAObjectCollection(dia_objects)

        return dia_collection

    @pipeBase.timeMethod
    def store(self, dia_collection, compute_spatial_index=False):
        """ Store all DIAObjects and DIASources in this dia_collection.

        Parameters
        ----------
        dia_collection : lsst.ap.association.DIAObjectCollection
            Collection of DIAObjects to store. Also stores the DIASources
            associated with these DIAObjects.
        compute_spatial_index : bool
            If True, compute the spatial search indices using the
            indexer specified at class creation.
        """
        for dia_object in dia_collection.dia_objects:
            if compute_spatial_index:
                dia_object.dia_object_record.set(
                    'indexer_id', self.indexer.index_points(
                        [dia_object.ra], [dia_object.dec])[0])
            self.store_dia_object(dia_object)
            for dia_source in dia_object.dia_source_catalog:
                self.store_dia_source(dia_source)
                self.store_dia_object_source_pair(
                    dia_object.id, dia_source.getId())
        self.commit()

    @pipeBase.timeMethod
    def store_updated(self, dia_collection, updated_indices):
        """ Store new DIAObjects and sources in the sqlite database.

        Parameters
        ----------
        dia_collection : lsst.ap.association.DIAObjectCollection
            A collection of DIAOjbects containing newly created or updated
            DIAObjects.
        updated_indices : int ndarray
            Indices within the set DIAObjctCollection that should be stored as
            updated DIAObjects in the database.
        """
        for updated_collection_index in updated_indices:
            dia_object = dia_collection.dia_objects[updated_collection_index]
            if dia_object.n_dia_sources == 1:
                dia_object.dia_object_record.set(
                    'indexer_id', self.indexer.index_points(
                        [dia_object.ra], [dia_object.dec])[0])
            self.store_dia_object(dia_object)
            self.store_dia_source(dia_object.dia_source_catalog[-1])
            self.store_dia_object_source_pair(
                dia_object.id,
                dia_object.dia_source_catalog[-1].getId())

        self.commit()

    def get_dia_objects(self, indexer_indices):
        """ Retrive the DIAObjects from the database whose HTM indices
        are within the specified range.

        Parameters
        ----------
        indexer_indices : array like of ints
            Pixelized indexer indices from which to load.

        Returns
        -------
        list of lsst.ap.association.DIAObjects
        """
        output_dia_objects = []

        dia_object_schema = make_minimal_dia_object_schema()

        n_indices = len(indexer_indices)
        for batch_idx in range(int(np.ceil(n_indices / 127))):
            indices_left = (n_indices - batch_idx * 127)
            if indices_left < 127:
                db_string = ("?," * (indices_left))[:-1]
            else:
                db_string = ("?," * 127)[:-1]
            self._db_cursor.execute(
                "SELECT * FROM dia_objects WHERE indexer_id IN (%s)" %
                db_string, indexer_indices[
                    batch_idx * 127: (batch_idx + 1) * 127])
            for row in self._db_cursor.fetchall():
                dia_object_record = afwTable.SourceTable.makeRecord(
                    afwTable.SourceTable.make(dia_object_schema))
                self._edit_source_record(row, dia_object_record,
                                         dia_object_schema)
                dia_sources = self.get_associated_dia_sources(
                    dia_object_record.getId())
                output_dia_objects.append(
                    DIAObject(dia_sources, dia_object_record))

        return output_dia_objects

    def get_dia_object_records(self, indexer_indices):
        """ Retrive the SourceRecord objects representing the DIAObjects
        in the catalog.

        Retrives the SourceRecords that are covered by the pixels with
        indices, indexer_indices. Use this to retrive the summary statistics
        of the DIAObjects themselves rather than taking the extra overhead
        of loading the associated DIASources as well.

        Parameters
        ----------
        indexer_indices : list of ints
            Pixelized indexer indices from which to load.

        Returns
        -------
        a lsst.afw.table.SourceCatalog
        """
        dia_object_schema = make_minimal_dia_object_schema()
        output_dia_objects = afwTable.SourceCatalog(dia_object_schema)

        dia_object_schema = make_minimal_dia_object_schema()
        for indexer_idx in indexer_indices:
            self._db_cursor.execute(
                "SELECT * FROM dia_objects WHERE indexer_id = ?",
                (indexer_idx,))
            rows = self._db_cursor.fetchall()
            for row in rows:
                dia_object_record = afwTable.SourceTable.makeRecord(
                    afwTable.SourceTable.make(dia_object_schema))
                self._edit_source_record(row, dia_object_record,
                                         dia_object_schema)
                output_dia_objects.append(dia_object_record)

        return output_dia_objects

    def get_associated_dia_sources(self, dia_obj_id):
        """ Retrive all DIASources associated with this DIAObject.

        Parameters
        ----------
        dia_obj_id : int
            Id of the DIAObject that is asssociated with the DIASources
            of interest.

        Returns
        -------
        lsst.afw.table.SourceCatalog
            SourceCatalog of DIASources associated with the DIAObject
        """
        self._db_cursor.execute(
            "SELECT src_id FROM dia_objects_to_dia_sources "
            "WHERE obj_id = ?",
            (dia_obj_id,))
        src_ids = self._db_cursor.fetchall()

        return self.get_dia_sources(np.array(src_ids, dtype=np.int).flatten())

    def get_dia_sources(self, dia_source_ids):
        """ Retrive the DIASources from the database

        Parameters
        ----------
        dia_source_id : list of ints
            Ids of the DIASources to load.

        Returns
        -------
        lsst.afw.table.SourceCatalog
            SourceCatalog of DIASources associated with the DIAObject
        """

        dia_source_schema = make_minimal_dia_source_schema()
        output_dia_sources = afwTable.SourceCatalog(dia_source_schema)
        n_sources = len(dia_source_ids)
        output_dia_sources.reserve(len(dia_source_ids))

        for batch_idx in range(int(np.ceil(n_sources / 127))):
            n_sources_left = (n_sources - batch_idx * 127)
            if n_sources_left < 127:
                db_string = ("?," * n_sources_left)[:-1]
            else:
                db_string = ("?," * 127)[:-1]
            self._db_cursor.execute(
                "SELECT * FROM dia_sources WHERE id IN (%s)" %
                db_string, dia_source_ids[
                    batch_idx * 127: (batch_idx + 1) * 127])
            for src_idx, row in enumerate(self._db_cursor.fetchall()):
                output_dia_sources.addNew()
                self._edit_source_record(row, output_dia_sources[src_idx],
                                         dia_source_schema)

        return output_dia_sources

    def _edit_source_record(self, db_row, src_record, schema):
        """ Edit the contents of a DIASource record.

        Parameters
        ----------
        db_row : list of values
            Database values retived to write into the sr_record
        src_record : lsst.afw.table.SourceRecord
            SourceRecord object to write values into
        schema : lsst.afw.table.Schema
            Schema defining the columns in src_record.
        """
        for sub_schema, value in zip(schema, db_row):
            if sub_schema.getField().getTypeString() == 'Angle':
                src_record.set(
                    sub_schema.getKey(),
                    afwGeom.Angle(value, units=afwGeom.degrees))
            else:
                src_record.set(
                    sub_schema.getKey(), value)

    def store_dia_object_source_pair(self, obj_id, src_id):
        """ Store a link between a DIAObject id and a DIASource.

        Parameters
        ----------
        obj_id : int
            Id of DIAObject
        src_id : int
            Id of DIASource
        """
        self._db_cursor.execute(
            "INSERT OR REPLACE INTO dia_objects_to_dia_sources "
            "VALUES (?, ?)", (src_id, obj_id))

    def store_dia_object(self, dia_object):
        """ Store an individual DIAObject into the database.

        Also stores the mapping of parent DIAObject id to DIASource id.

        Parameters
        ----------
        dia_object : lsst.ap.association.DIAObject
            DIAObject to store in the database. If the DIAObject already
            exists it is replaced in the database with the updated centroids,
            fluxes, etc.
        """
        insert_string = ""
        values = []
        for sub_schema in dia_object.dia_object_record.getSchema():
            if sub_schema.getField().getTypeString() == 'Angle':
                values.append(
                    dia_object.get(sub_schema.getKey()).asDegrees())
            else:
                values.append(dia_object.get(sub_schema.getKey()))
            insert_string += '?,'
        insert_string = insert_string[:-1]

        self._db_cursor.execute(
            "INSERT OR REPLACE INTO dia_objects VALUES (%s)" % insert_string,
            values)

    def store_dia_source(self, dia_source):
        """ Store this DIASource in the database.

        Parameters
        ----------
        dia_soruce : lsst.afw.table.SourceRecord
            DIASource object whose values we would like to store.
        """
        insert_string = ""
        values = []
        for sub_schema in dia_source.getSchema():
            if sub_schema.getField().getTypeString() == 'Angle':
                values.append(
                    dia_source.get(sub_schema.getKey()).asDegrees())
            else:
                values.append(dia_source.get(sub_schema.getKey()))
            insert_string += '?,'
        insert_string = insert_string[:-1]

        self._db_cursor.execute(
            "INSERT OR REPLACE INTO dia_sources VALUES (%s)" % insert_string,
            values)
