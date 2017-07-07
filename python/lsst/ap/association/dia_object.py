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

import numpy as np

import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom

__all__ = ["DIAObject", "make_minimal_dia_object_schema"]


def make_minimal_dia_object_schema():
    """ Define and create the minimal schema required for a DIAObject.

    Return
    ------
    lsst.afw.table.schema.schema.Schema
    """

    schema = afwTable.SourceTable.makeMinimalSchema()

    # For the MVP/S we currently only care about the scatter on the possition
    # so we add those.
    schema.addField("coord_ra_rms", doc="rms position in ra/dec",
                    type="Angle")
    schema.addField("coord_dec_rms", doc="rms position in ra/dec",
                    type="Angle")

    return schema


class DIAObject(object):
    """ A class specifying a collection of single frame difference image
    sources and statistics on these collections.

    Attributes
    ----------
    _dia_object_record : lsst.afw.table.SourceRecord
        A SourceRecord object containing the summary statistics for the
        collection of DIASources this DIAObject represents (e.g. median
        RA/DEC position).
    _dia_source_catalog : lsst.afw.table.SourceCatalog
        A set of SourceRecords specifying the DIASources that make up
        this DIAObject.
    _updated : bool
        boolean specifying if the summary statistics for this DIAObject have
        been updated with the current set of DIASources in the SourceCatalog.
        This variable should be set to false whenever the SourceCatalog
        of DIASources changes and set to true when the initialize method is
        run.
    """
    def __init__(self, dia_source_catalog, object_source_record=None):
        """  Create a DIAObject given an input SourceCatalog of
        DIASources.

        Takes as input an lsst.afw.table.SourceCaatlog object specifying a
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
        """ Return the data stored in column name within the internal
        dia_object_record.

        Parameters
        ----------
        name : str or lsst.afw.table.Key

        Return
        ------
        An lsst.afw data type
        """

        # This will in the future be replaced with a overwritting of __getattr
        # and __dir__ for this class.
        return self._dia_object_record.get(name)

    def update(self):
        """ Compute all summary statistics given the current catalog of
        DIASources asigned to this DIAObject.

        Store these summaries (e.g. median RA/DEC position, fluxes...) in
        the object_source_record attribute and set the class variable
        updated to True

        Returns
        -------
        None
        """

        self._updated = False

        # To quickly compute the summary statistics we check if the catalog
        # is currently contious and if not we make a deep copy.
        if not self._dia_source_catalog.isContiguous():
            tmp_dia_source_catalog = self._dia_source_catalog.copy(deep=True)
            del self._dia_source_catalog
            self._dia_source_catalog = tmp_dia_source_catalog

        self._compute_summary_statistics()

        self._updated = True

        return None

    def _compute_summary_statistics(self):
        """ Retrive properties from DIASourceCatalog attribute and update the
        summary statistics that represent this DIAObject

        Returns
        -------
        None
        """

        # Loop through DIASources, compute summary statistics (TBD) and store
        # them in dia_object_record attribute.

        for name in self.schema.getNames():
            # For the MVP/S we are only dealing with angles so we skip over
            # everything that isn't coord_ra or coord_dec.
            if name != 'coord_ra' and name != 'coord_dec':
                continue

            mean_value = np.mean(self._dia_source_catalog[name])
            # Currently hard coded to work with angles.
            self._dia_object_record.set(name, afwGeom.Angle(mean_value))
            # If we only have 1 source we can't compute an rms and hence we
            # set the the rms variables to NaN.
            if self.n_dia_sources > 1:
                self._dia_object_record[name + '_rms'] = afwGeom.Angle(
                    np.std(self._dia_source_catalog[name]))
            else:
                self._dia_object_record[name + '_rms'] = afwGeom.Angle(np.nan)

        return None

    def append_dia_source(self, input_dia_source_record):
        """ Append the input_dia_source to the dia_source_catalog attribute.

        Additionally set update boolean to False.

        Parameters
        ----------
        input_dia_source : lsst.afw.table.SourceRecord
            Single DIASource object to append to this DIAObject's source
            catalog.

        Return
        ------
        None
        """

        # Since we are adding to the SourceCatalog our summary statistics are
        # no longer valid. We set this to false and hold off on recomputing
        # them until we are finished adding sources.
        self._updated = False

        self._dia_source_catalog.append(input_dia_source_record)

        return None

    def get_light_curve(self):
        """ Retreve the light curve of fluxes for the DIASources that make up
        this DIAObject.

        Returns
        -------
        An array like object specifying the light curve for this object.
        """

        # Loop through DIASources and return the "light curve"
        # Right now I'm making this the same as returning the
        # dia_source_catalog.

        return self.dia_source_catalog()

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
        """ Retrive the SourceRecord that represents the summary statistics on
        this DIAObject's set of DIASources.

        Return
        ------
        A lsst.afw.table.SourceRecord
        """
        return self._dia_object_record

    @property
    def dia_source_catalog(self):
        """ Retrive the SourceCatalog that represents the DIASources that make
        up this DIAObject.

        Return
        ------
        A lsst.afw.table.SourceCatalog
        """
        return self._dia_source_catalog

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
