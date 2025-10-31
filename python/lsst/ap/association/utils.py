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
__all__ = ("convertDataFrameToSdmSchema", "readSdmSchemaFile", "readSchemaFromApdb",
           "dropEmptyColumns", "make_empty_catalog", "getMidpointFromTimespan",
           "makeEmptyForcedSourceTable", "checkSdmSchemaColumns", "getRegion")

from collections.abc import Mapping
import os

from lsst.dax.apdb import Apdb, ApdbTables, schema_model
import felis.datamodel
import numpy as np
import pandas as pd
from astropy.table import Table
import yaml

from lsst.daf.butler import Timespan
import lsst.dax.apdb as daxApdb
import lsst.geom
import lsst.sphgeom


# The first entry in the returned mapping is for nullable columns,
# the second entry is for non-nullable columns.
_dtype_map: Mapping[felis.datamodel.DataType, tuple[str, str]] = {
    felis.datamodel.DataType.double: ("float64", "float64"),  # Cassandra utilities need np.nan not pd.NA
    felis.datamodel.DataType.float: ("float32", "float32"),  # Cassandra utilities need np.nan not pd.NA
    felis.datamodel.DataType.timestamp: ("datetime64[ms]", "datetime64[ms]"),
    felis.datamodel.DataType.long: ("Int64", "int64"),
    felis.datamodel.DataType.int: ("Int32", "int32"),
    felis.datamodel.DataType.short: ("Int16", "int16"),
    felis.datamodel.DataType.byte: ("Int8", "int8"),
    felis.datamodel.DataType.binary: ("object", "object"),
    felis.datamodel.DataType.char: ("object", "object"),
    felis.datamodel.DataType.text: ("object", "object"),
    felis.datamodel.DataType.string: ("object", "object"),
    felis.datamodel.DataType.unicode: ("object", "object"),
    felis.datamodel.DataType.boolean: ("boolean", "bool"),
}


def column_dtype(felis_type: felis.datamodel.DataType | daxApdb.schema_model.ExtraDataTypes,
                 nullable=False) -> str:
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
        return _dtype_map[felis_type][0] if nullable else _dtype_map[felis_type][1]
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


def checkSdmSchemaColumns(apdbSchema, colNames, tableName):
    """Check if supplied column names exists in the schema.

    Parameters
    ----------
    apdbSchema : `dict` [`str`, `lsst.dax.apdb.schema_model.Table`]
        Schema from ``sdm_schemas`` containing the table definition to use.
    colNames : `list` of ``str`
        Names of the columns to check for in the table.
    tableName : `str`
        Name of the table in the schema to use.

    Returns
    -------
    missing : `list` of `str`
        All column names that are not in the schema
    """
    table = apdbSchema[tableName]
    missing = []

    names = [columnDef.name for columnDef in table.columns]
    for col in colNames:
        if col not in names:
            missing.append(col)
    return missing


def convertDataFrameToSdmSchema(apdbSchema, sourceTable, tableName):
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
    if sourceTable.empty:
        make_empty_catalog(apdbSchema, tableName)
    table = apdbSchema[tableName]

    data = {}
    nSrc = len(sourceTable)

    for columnDef in table.columns:
        dtype = column_dtype(columnDef.datatype, nullable=columnDef.nullable)
        if columnDef.name in sourceTable.columns:
            data[columnDef.name] = pd.Series(sourceTable[columnDef.name], dtype=dtype,
                                             index=sourceTable.index)
        else:
            if columnDef.nullable:
                try:
                    data[columnDef.name] = pd.Series([pd.NA]*nSrc, dtype=dtype, index=sourceTable.index)
                except TypeError:
                    data[columnDef.name] = pd.Series([np.nan]*nSrc, dtype=dtype, index=sourceTable.index)
            else:
                data[columnDef.name] = pd.Series([0]*nSrc, dtype=dtype, index=sourceTable.index)
    return pd.DataFrame(data)


def convertTableToSdmSchema(apdbSchema, sourceTable, tableName):
    """Force a table to conform to the schema defined by the APDB.

    This method uses the table definitions in ``sdm_schemas`` to
    load the schema of the APDB, and does not actually connect to the APDB.

    Parameters
    ----------
    apdbSchema : `dict` [`str`, `lsst.dax.apdb.schema_model.Table`]
        Schema from ``sdm_schemas`` containing the table definition to use.
    sourceTable : `astropy.table.Table`
        The input table to convert.
    tableName : `str`
        Name of the table in the schema to use.

    Returns
    -------
    `astropy.table.Table`
        A table with the correct schema for the APDB and data copied from
        the input ``sourceTable``.
    """
    table = apdbSchema[tableName]

    data = {}
    nSrc = len(sourceTable)

    for columnDef in table.columns:
        dtype = column_dtype(columnDef.datatype, nullable=columnDef.nullable)
        if columnDef.name in sourceTable.columns:
            data[columnDef.name] = Table.Column(sourceTable[columnDef.name], dtype=dtype.lower())
        else:
            if columnDef.nullable:
                try:
                    data[columnDef.name] = Table.Column([pd.NA]*nSrc, dtype=object)
                except TypeError:
                    data[columnDef.name] = Table.Column([pd.nan]*nSrc, dtype=dtype)
            else:
                data[columnDef.name] = Table.Column([0]*nSrc, dtype=dtype)
    return Table(data)


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
        columnDef.name: pd.Series(dtype=column_dtype(columnDef.datatype, nullable=columnDef.nullable))
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
    diaForcedSources = convertDataFrameToSdmSchema(schema, pd.DataFrame(), tableName="DiaForcedSource")
    return diaForcedSources


def getRegion(exposure):
    """Calculate an enveloping region for an exposure.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure object with calibrated WCS.

    Returns
    -------
    region : `lsst.sphgeom.Region`
        Region enveloping an exposure.
    """
    # Bounding box needs to be a `Box2D` not a `Box2I` for `wcs.pixelToSky()`
    bbox = lsst.geom.Box2D(exposure.getBBox())
    wcs = exposure.getWcs()

    region = lsst.sphgeom.ConvexPolygon([pp.getVector() for pp in wcs.pixelToSky(bbox.getCorners())])

    return region
