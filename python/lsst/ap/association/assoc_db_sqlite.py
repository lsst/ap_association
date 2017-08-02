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

__all__ = ["AssociationDBSqliteConfig", "AssociationDBSqliteTask"]

import sqlite3

from lsst.meas.algorithms.indexerRegistry import IndexerRegistry
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from .dia_collection import *
from .dia_object import *


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
        doc = 'Select the spatial indexer to use within the database.',
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
        """ Create a new connection to the sqlite database specified in
        the config.
        """
        pipeBase.Task.__init__(self, **kwargs)
        self.indexer = IndexerRegistry[self.config.indexer.name](self.config.indexer.active)
        self._db_connection = sqlite3.connect(self.config.db_name)
        self._db_cursor = self._db_connection.cursor()

    def run(self):
        """ Not implemented in this subtask. Please use the store and load
        methods.
        """
        Raise(NotImplementedError)

    def commit(self):
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
            name_type_string = ""
            for sub_schema in obj_schema:
                tmp_name = sub_schema.getField().getName()
                tmp_type = afw_to_db_types[
                    sub_schema.getField().getTypeString()]
                if tmp_name == 'id':
                    tmp_type += " PRIMARY KEY"  
                name_type_string += "%s %s, " % (tmp_name,tmp_type)
            name_type_string = name_type_string[:-2]
            self._db_cursor.execute(
                "CREATE TABLE dia_objects (%s)" % name_type_string)

            name_type_string = ""
            for sub_schema in src_schema:
                tmp_name = sub_schema.getField().getName()
                tmp_type = afw_to_db_types[
                    sub_schema.getField().getTypeString()]
                if tmp_name == 'id':
                    tmp_type += " PRIMARY KEY"
                name_type_string += "%s %s, " % (tmp_name,tmp_type)
            name_type_string = name_type_string[:-2]
            self._db_cursor.execute(
                "CREATE TABLE dia_sources (%s)" % name_type_string)

            self._db_cursor.execute(
                "CREATE TABLE dia_objects_to_dia_soruces "
                "(obj_id INTEGER, src_id INTEGER PRIMARY KEY)")

            self._db_connection.commit()

        return True

    @pipeBase.timeMethod
    def load(self, ctr_coord, radius):
        """ Load all DIAObjects and associated DIASources.

        Parameters
        ----------
        ctr_coord : lsst.afw.geom.Coord
            Center position of on on sky circle to load.
        radius : lsst.afw.geom.Angle


        Returns
        -------
        lsst.ap.association.DIAObjectCollection
        """
        indexer_indices, on_boundry = self.indexer.get_pixel_ids(
            ctr_coord, radius)

        dia_object_records = get_dia_object_records(indexer_indices)
        dia_object_list = []

        for dia_object_record in dia_objects:
            dia_sources = self.get_dia_sources(dia_object_record.getId())
            dia_object_list.append(DIAObject(dia_sources, dia_object_record))

        dia_colletion = DIAObjectCollection(dia_object_list)

        return dia_collection

    @pipeBase.timeMethod
    def store(self, dia_collection, updated_indices):
        """ Store new DIAObjects and sources in the sqlite database.

        Parameters
        ----------
        dia_collection : lsst.ap.association.DIAObjectCollection
        updated_indices : int ndarray
            Indices within the DIAObjctCollection that contain DIAObjects
            that were updated by AssociationTask.
        """
        for updated_index in updated_indices:
            dia_object = dia_collection.dia_objects[updated_index]
            if dia_object.n_dia_sources == 1:
                dia_object.dia_object_record.set(
                    'indexer_id', self.indexer.index_points(
                        [dia_object.ra], [dia_object.dec])[0])
            self.store_dia_object(dia_object)

        self.commit()

    def get_dia_object_records(self, indexer_indices):
        """ Retrive the DIAObjects from the database whose HTM indices
        are within the specified range.

        Parameters
        ----------
        indexer_indices : list of ints
            Pixelized indexer indices from which to load.

        Returns
        -------
        list of lsst.ap.association.DIAObjects
        """
        output_dia_objects = []

        self._db_cursor.excecutemany(
            "SELECT * FROM dia_objects WHERE indexer_id = ?",
            indexer_indices)

        rows = self._db_cursor.fetchall()

        dia_object_schema = make_minimal_dia_object_schema()
        for row in rows:
            dia_object_record = afwTable.SourceTable.makeRecord(
                afwTable.SourceTable.make(dia_object_schema))
            self._edit_source_record(row, dia_object_record,
                                     dia_object_schema)
            dia_sources = self.get_dia_sources(dia_object_record.getId())
            output_dia_objects.append(
                DIAObject(dia_sources, dia_object_record))

        return output_dia_objects

    def get_dia_sources(self, dia_obj_id):
        """ Retrive the DIASources association with the specified DIAObject
        id.

        Parameters
        ----------
        dia_obj_id : int
            Load the DIASource association with the DIAObject with id
            dia_obj_id

        Returns
        -------
        lsst.afw.table.SourceCatalog
            SourceCatalog of DIASources associated with the DIAObject
        """
        self._db_cursor.excecute(
            "SELECT dia_source_id FROM dia_objects_to_dia_soruces "
            "WHERE dia_object_id = ?",
            dia_obj_id)
        src_ids = self._db_cursor.fetchall()

        dia_source_schema =  make_minimal_dia_source_schema()
        output_dia_sources = afwTable.SourceCatalog(dia_source_schema)
        output_dia_sources.reserve(len(src_ids))

        self._db_cursor.excecutemany(
            "SELECT * FROM dia_soucres WHERE id = ?", src_ids)
        for src_idx, row in self._db_cursor.fetchall():
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
            if sub_schema.getField().getTypeString() == 'Anlge':
                src_record.set(
                    sub_schema.getKey(),
                    afwGeom.Angle(value, units=afwGeom.degrees))
            else:
                src_record.set(
                    sub_schema.getKey(), value)


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

        self._db_cursor.execute(
            "INSERT INTO dia_objects_to_dia_soruces VALUES (?, ?)",
            (dia_object.get('id'), dia_object.dia_source_catalog[-1].getId()))

        insert_string = ""
        values = []
        for sub_schema in dia_object.dia_object_record.getSchema():
            if sub_schema.getField().getTypeString() == 'Angle':
                values.append(
                    dia_object.get(sub_schema.getKey()).asDegrees())
            else:
                values.append(dia_object.get(sub_schema.getKey()))
            insert_string += '?, '
        insert_string = insert_string[:-2]

        self._db_cursor.execute(
            "INSERT OR REPLACE INTO dia_objects VALUES (%s)" % insert_string,
            values)

        self.store_dia_source(dia_object.dia_source_catalog[-1])

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
            insert_string += '?, '
        insert_string = insert_string[:-2]

        self._db_cursor.execute(
            "INSERT INTO dia_sources VALUES (%s)" % insert_string,
            values)
