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

"""Utilities for working with the APDB.
"""
__all__ = ("convertTableToSdmSchema", "readSdmSchemaFile", "readSchemaFromApdb",
           "dropEmptyColumns", "make_empty_catalog", "getMidpointFromTimespan",
           "makeEmptyForcedSourceTable")

from collections.abc import Mapping
import os

from lsst.dax.apdb import Apdb, ApdbTables, schema_model
import felis.datamodel
import numpy as np
import pandas as pd
import yaml

from lsst.daf.butler import Timespan
import lsst.dax.apdb as daxApdb


_dtype_map: Mapping[felis.datamodel.DataType | daxApdb.schema_model.ExtraDataTypes, type | str] = {
    felis.datamodel.DataType.double: np.float64,
    felis.datamodel.DataType.float: np.float32,
    felis.datamodel.DataType.timestamp: "datetime64[ms]",
    felis.datamodel.DataType.long: np.int64,
    felis.datamodel.DataType.int: np.int32,
    felis.datamodel.DataType.short: np.int16,
    felis.datamodel.DataType.byte: np.int8,
    felis.datamodel.DataType.binary: object,
    felis.datamodel.DataType.char: object,
    felis.datamodel.DataType.text: object,
    felis.datamodel.DataType.string: object,
    felis.datamodel.DataType.unicode: object,
    felis.datamodel.DataType.boolean: bool,
}


def column_dtype(felis_type: felis.datamodel.DataType | daxApdb.schema_model.ExtraDataTypes) -> type | str:
    """Return Pandas data type for a given Felis column type.

    Parameters
    ----------
    felis_type : `felis.datamodel.DataType`
        Felis type, on of the enums defined in `felis.datamodel` module.

    Returns
    -------
    column_dtype : `type` or `str`
        Type that can be used for columns in Pandas.

    Raises
    ------
    TypeError
        Raised if type is cannot be handled.
    """
    try:
        return _dtype_map[felis_type]
    except KeyError:
        raise TypeError(f"Unexpected Felis type: {felis_type}")


def readSdmSchemaFile(schemaFile: str,
                      schemaName: str = "ApdbSchema",
                      ):
    """Read a schema file in YAML format.

    Parameters
    ----------
    schemaFile : `str`
        Fully specified path to the file to be read.
    schemaName : `str`, optional
        Name of the table of schemas to read from the file.

    Returns
    -------
    schemaTable : dict[str, schema_model.Table]
        A dict of the schemas in the given table defined in the specified file.

    Raises
    ------
    ValueError
        If the schema file can't be parsed.
    """
    schemaFile = os.path.expandvars(schemaFile)
    with open(schemaFile) as yaml_stream:
        schemas_list = list(yaml.load_all(yaml_stream, Loader=yaml.SafeLoader))
        schemas_list = [schema for schema in schemas_list if schema.get("name") == schemaName]
        if not schemas_list:
            raise ValueError(f"Schema file {schemaFile!r} does not define schema {schemaName!r}")
        elif len(schemas_list) > 1:
            raise ValueError(f"Schema file {schemaFile!r} defines multiple schemas {schemaName!r}")
        felis_schema = felis.datamodel.Schema.model_validate(schemas_list[0])
        schema = daxApdb.schema_model.Schema.from_felis(felis_schema)
    schemaTable = {}

    for singleTable in schema.tables:
        schemaTable[singleTable.name] = singleTable
    return schemaTable


def readSchemaFromApdb(apdb: Apdb) -> dict[str, schema_model.Table | None]:
    """Extract the schema from an APDB instance.

    Parameters
    ----------
    apdb : `lsst.dax.apdb.Apdb`
        AP database connection object.

    Returns
    -------
    schemaTable : dict[str, schema_model.Table | None]
        A dict of the schemas in the given table defined in the specified file.
    """
    return {table.table_name(): apdb.tableDef(table) for table in ApdbTables}


def convertTableToSdmSchema(apdbSchema, sourceTable, tableName):
    """Force a table to conform to the schema defined by the APDB.

    This method uses the table definitions in ``sdm_schemas`` to
    load the schema of the APDB, and does not actually connect to the APDB.

    Parameters
    ----------
    apdbSchema : `dict` [`str`, `lsst.dax.apdb.schema_model.Table`]
        Schema from ``sdm_schemas`` containing the table definition to use.
    sourceTable : `pandas.DataFrame`
        The input table to convert.
    tableName : `str`
        Name of the table in the schema to use.

    Returns
    -------
    `pandas.DataFrame`
        A table with the correct schema for the APDB and data copied from
        the input ``sourceTable``.
    """
    table = apdbSchema[tableName]

    data = {}
    nSrc = len(sourceTable)

    for columnDef in table.columns:
        if columnDef.name in sourceTable.columns:
            dtype = column_dtype(columnDef.datatype)
            data[columnDef.name] = pd.Series(sourceTable[columnDef.name], dtype=dtype)
        else:
            dataInit = np.zeros(nSrc, dtype=column_dtype(columnDef.datatype))
            if columnDef.nullable:
                try:
                    dataInit.fill(None)
                except (TypeError, ValueError):
                    pass
            data[columnDef.name] = pd.Series(dataInit, index=sourceTable.index)
    return pd.DataFrame(data)


def dropEmptyColumns(apdbSchema, sourceTable, tableName):
    """Drop empty columns that are nullable.

    This method uses the table definitions in ``sdm_schemas`` to
    load the schema of the APDB, and does not actually connect to the APDB.

    Parameters
    ----------
    apdbSchema : `dict` [`str`, `lsst.dax.apdb.schema_model.Table`]
        Schema from ``sdm_schemas`` containing the table definition to use.
    sourceTable : `pandas.DataFrame`
        The input table to remove missing data columns from.
    tableName : `str`
        Name of the table in the schema to use.
    """
    table = apdbSchema[tableName]

    nullableList = [columnDef.name for columnDef in table.columns if columnDef.nullable]
    nullColumns = sourceTable.isnull().all()
    nullColNames = nullColumns[nullColumns].index.tolist()
    dropColumns = list(set(nullColNames) & set(nullableList))
    return sourceTable.drop(columns=dropColumns)


def make_empty_catalog(apdbSchema, tableName):
    """Make an empty catalog for a table with a given name.

    Parameters
    ----------
    apdbSchema : `dict` [`str`, `lsst.dax.apdb.schema_model.Table`]
        Schema from ``sdm_schemas`` containing the table definition to use.
    tableName : `str`
        Name of the table in the schema to use.

    Returns
    -------
    catalog : `pandas.DataFrame`
        An empty catalog.
    """
    table = apdbSchema[tableName]

    data = {
        columnDef.name: pd.Series(dtype=column_dtype(columnDef.datatype))
        for columnDef in table.columns
    }
    return pd.DataFrame(data)


def getMidpointFromTimespan(timespan, allowUnbounded=True):
    """Safely retrieve the midpoint in TAI from a Timespan.

    Parameters
    ----------
    timespan : `lsst.daf.butler.Timespan` of `astropy.time.Time`
        A Timespan centered on the midpoint of a visit.
    allowUnbounded : `bool`, optional
        If set, return the start or end of an unbounded timespan.

    Returns
    -------
    midpoint : `astropy.time.Time`
        The midpoint of the timespan.

    Raises
    ------
    ValueError
        Raised if either the start or end of the timespan is None, and
        ``allowUnbounded`` is not set.
    ValueError
        Raised if the timespan is empty.
    """
    if (timespan.begin == Timespan.EMPTY) or (timespan.begin == Timespan.EMPTY):
        raise ValueError("Cannot compute midpoint: EMPTY Timespan.")

    try:
        interval = timespan.end - timespan.begin
        return (timespan.begin + interval/2).tai
    except TypeError as e:
        if allowUnbounded:
            if timespan.end is not None:
                return timespan.end.tai
            elif timespan.begin is not None:
                return timespan.begin.tai
            else:
                raise ValueError("Cannot compute midpoint: unbounded timespan.") from e
        else:
            raise ValueError("Cannot compute midpoint: unbounded timespan.") from e


def makeEmptyForcedSourceTable(schema):
    """Return a dataframe with the correct columns for diaForcedSources table.

    Returns
    -------
    diaForcedSources : `pandas.DataFrame`
        Empty dataframe.
    """
    diaForcedSources = convertTableToSdmSchema(schema, pd.DataFrame(), tableName="DiaForcedSource")
    return diaForcedSources
