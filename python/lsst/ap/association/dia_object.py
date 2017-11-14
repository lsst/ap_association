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

""" Definitions for DIAObject and a minimal schema for them.
"""

from __future__ import absolute_import, division, print_function

import numpy as np

import lsst.afw.table as afwTable
from lsst.afw.coord import averageCoord

__all__ = ["DIAObject",
           "make_minimal_dia_object_schema",
           "make_minimal_dia_source_schema"]


def make_minimal_dia_object_schema():
    """ Define and create the minimal schema required for a DIAObject.

    Return
    ------
    lsst.afw.table.schema.schema.Schema
    """

    schema = afwTable.SourceTable.makeMinimalSchema()
    # For the MVP/S we currently only care about the position though
    # in the future we will add summary computations for fluxes etc.
    # as well as their errors.

    # In the future we would like to store a covariance of the coordinate.
    # This functionality is not defined in currently in the stack, so we will
    # hold off until it is implemented. This is to be addressed in DM-7101.
    schema.addField("indexer_id", type=np.int64)
    schema.addField("n_dia_sources", type=np.int64)

    return schema


def make_minimal_dia_source_schema():
    """ Define and create the minimal schema required for a DIASource.

    Return
    ------
    lsst.afw.table.schema.schema.Schema
    """

    schema = afwTable.SourceTable.makeMinimalSchema()

    return schema


class DIAObject(object):
    """ A class specifying a collection of single frame difference image
    sources and statistics on these collections.

    """
    def __init__(self, dia_source_catalog, object_source_record=None):
        """  Create a DIAObject given an input SourceCatalog of
        DIASources.

        Takes as input an lsst.afw.table.SourceCatlog object specifying a
        collection of DIASources that make up this DIAObject. The optional
        input object_source_record should contain summary statistics on the
        SourceCatalog of DIASources. Using this optional input escapes the
        need to recompute the summary statistics when not necessary.

        Parameters
        ----------
        dia_source_catalog : lsst.afw.table.SourceCatalog
            SourceCatalog of DIASources associated to this DIAObject
        object_source_record : lsst.afw.table.SourceRecord, optional
            Optional input SourceRecord containing summary statistics on
            the input SourceCatalog.

        Returns
        -------
        A DIAObject instance
        """

        self._dia_source_catalog = dia_source_catalog
        self._dia_source_schema = self._dia_source_catalog.getSchema()
        self._updated = False

        if object_source_record is None:
            self._dia_object_record = afwTable.SourceTable.makeRecord(
                afwTable.SourceTable.make(
                    make_minimal_dia_object_schema()))
            # For now we copy the first DIASource's object id as our DIAObject
            # id.
            self._dia_object_record.set(
                'id', self._dia_source_catalog[0].getId())
            self.update()
        else:
            self._dia_object_record = object_source_record
            self._updated = True

    def get(self, name):
        """ Retrieve a specific summary statistic from this DIAObject

        Parameters
        ----------
        name : str or lsst.afw.table.Key

        Return
        ------
        A SourceRecord column value
        """

        # This will in the future be replaced with a overwriting of __getattr
        # and __dir__ for this class.
        return self._dia_object_record.get(name)

    def update(self):
        """ Compute all summary statistics given the current catalog of
        DIASources assigned to this DIAObject.

        Store these summaries (e.g. median RA/DEC position, fluxes...) in
        the object_source_record attribute and set the class variable
        updated to True.
        """

        self._updated = False

        # To quickly compute the summary statistics we check if the catalog
        # is currently continuous and if not we make a deep copy.
        if not self._dia_source_catalog.isContiguous():
            tmp_dia_source_catalog = self._dia_source_catalog.copy(deep=True)
            self._dia_source_catalog = tmp_dia_source_catalog

        self._compute_summary_statistics()

        self._updated = True

    def _compute_summary_statistics(self):
        """ Retrieve properties from DIASourceCatalog attribute and update the
        summary statistics that represent this DIAObject
        """

        # Loop through DIASources, compute summary statistics (TBD) and store
        # them in dia_object_record attribute.
        self._store_n_associated_sources()
        self._compute_mean_coordinate()

        # In the future we will calculate covariances on this centroid,
        # however generalized coordinate covariances are not defined (DM-7101)
        # we also do not need them yet for the MVP/S

    def _store_n_associated_sources(self):
        """ Store the number of DIASources associated with the DIAObject.
        """

        self._dia_object_record.set("n_dia_sources", self.n_dia_sources)

    def _compute_mean_coordinate(self):
        """ Compute the mean coordinate of this DIAObject given the current
        DIASources associated with it.
        """

        coord_list = [src.getCoord() for src in self._dia_source_catalog]
        ave_coord = averageCoord(coord_list)
        self._dia_object_record.setCoord(ave_coord)

    def append_dia_source(self, input_dia_source_record):
        """ Append the input_dia_source to the dia_source_catalog attribute.

        Additionally set update boolean to False.

        Parameters
        ----------
        input_dia_source : lsst.afw.table.SourceRecord
            Single DIASource object to append to this DIAObject's source
            catalog.
        """

        # Since we are adding to the SourceCatalog our summary statistics are
        # no longer valid. We set this to false and hold off on recomputing
        # them until we are finished adding sources.
        self._updated = False

        # We need to test that schema of the source to appended lines up with
        # the schema specified in the catalog of dia_sources associated with
        # this dia_object. If not we attempt to access the those that overlap
        # assuming that the schema defined in this dia_object is a subset of
        # of the input source's schema.
        input_schema = input_dia_source_record.getSchema()
        if input_schema == self._dia_source_schema:
            self._dia_source_catalog.append(
                self._dia_source_catalog.getTable().copyRecord(
                    input_dia_source_record))
        else:
            tmp_source_record = afwTable.SourceTable.makeRecord(
                afwTable.SourceTable.make(self._dia_source_schema))
            for name in self._dia_source_schema.getNames():
                tmp_source_record.set(
                    self._dia_source_schema[name].asKey(),
                    input_dia_source_record.getSchema()[name].asKey())

            self._dia_source_catalog.append(
                self._dia_source_catalog.getTable().copyRecord(
                    tmp_source_record))

    def get_light_curve(self):
        """ Retrieve the light curve of fluxes for the DIASources that make up
        this DIAObject.

        Returns
        -------
        An array like object specifying the light curve for this object.
        """

        # Loop through DIASources and return the "light curve"
        # Right now I'm making this the same as returning the
        # dia_source_catalog.

        raise NotImplementedError(
            "Light curves not yet implemented. Use dia_source_catalog property"
            "instead.")

    @property
    def is_updated(self):
        """ Return the current state of this DIAObject.

        If True is returned this DIAObject has computed all of its summary
        statistics for the current collection of DIASources that make it up.

        Return
        ------
        bool
        """

        return self._updated

    @property
    def dia_object_record(self):
        """ Retrieve the SourceRecord that represents the summary statistics on
        this DIAObject's set of DIASources.

        Return
        ------
        A lsst.afw.table.SourceRecord
        """
        return self._dia_object_record

    @property
    def dia_source_catalog(self):
        """ Retrieve the SourceCatalog that represents the DIASources that make
        up this DIAObject.

        Return
        ------
        A lsst.afw.table.SourceCatalog
        """
        return self._dia_source_catalog

    @property
    def dia_source_schema(self):
        """ Retrieve the SourceCatalog that represents the DIASources that make
        up this DIAObject.

        Return
        ------
        A lsst.afw.table.SourceCatalog
        """
        return self._dia_source_schema

    @property
    def n_dia_sources(self):
        """ Return the number of DIASources currently associated with this
        object.

        Return
        ------
        """
        return len(self._dia_source_catalog)

    @property
    def schema(self):
        """ Return the schema of the DIAObject record.

        Returns
        -------
        lsst.afw.table.schema.schema.Schema
        """
        return self._dia_object_record.schema

    @property
    def source_catalog_schema(self):
        """ Return the schema of the DIASourceCatalog associated with this
        DIAObject.

        Returns
        -------
        lsst.afw.table.schema.schema.Schema
        """
        return self._dia_source_catalog.schema

    @property
    def ra(self):
        """ Get the RA of this DIAObject.

        Return
        ------
        A lsst.afw.geom.angle.angle.Angle
        """
        return self._dia_object_record.getRa()

    @property
    def dec(self):
        """ Get the DEC of this DIAObject.

        Return
        ------
        A lsst.afw.geom.angle.angle.Angle
        """
        return self._dia_object_record.getDec()

    @property
    def coord(self):
        """ Get the Coordinate of this DIAObject.

        Return
        ------
        A lsst.afw.coord._coord.IcrsCoord
        """
        return self._dia_object_record.getCoord()

    @property
    def id(self):
        """ Get the unique catalog identifier for this object.

        Return
        ------
        int
        """
        return self._dia_object_record.getId()
