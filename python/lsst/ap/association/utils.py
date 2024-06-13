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
__all__ = ("convertTableToSdmSchema", "readSdmSchemaFile", "readSchemaFromApdb", "dropEmptyColumns")

from collections.abc import Mapping
import os

import felis.datamodel
import numpy as np
import pandas as pd
import yaml

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
    schemaTable : 'dict' of `lsst.dax.apdb.apdbSchema.ApdbSchema`
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


def readSchemaFromApdb(apdb):
    """Extract the schema from an APDB instance.

    Parameters
    ----------
    apdb : `lsst.dax.apdb.Apdb`
        AP database connection object.

    Returns
    -------
    schemaTable : 'dict' of `lsst.dax.apdb.apdbSchema.ApdbSchema`
        A dict of the schemas in the given table defined in the specified file.
    """
    schemaTable = {}
    for singleTable in apdb._schema.tableSchemas:
        schemaTable[singleTable.name] = apdb._schema.tableSchemas[singleTable]
    return schemaTable


def convertTableToSdmSchema(apdbSchema, sourceTable, tableName):
    """Force a table to conform to the schema defined by the APDB.

    This method uses the table definitions in ``sdm_schemas`` to
    load the schema of the APDB, and does not actually connect to the APDB.

    Parameters
    ----------
    apdbSchema : `lsst.dax.apdb.apdbSchema.ApdbSchema`
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
    apdbSchema : `lsst.dax.apdb.apdbSchema.ApdbSchema`
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
