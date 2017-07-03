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


class DIAObject(object):
    """ A class specifying a collection of single frame difference image
    sources and statistics on these collections.

    Attributes
    ----------
    dia_object_record : lsst.afw.table.SourceRecord
        A SourceRecord object containing the summary statistics for the
        collection of DIASources this DIAObject represents (e.g. median
        RA/DEC position).
    source_catalog : lsst.afw.table.SourceCatalog
        A set of SourceRecords specifying the DIA sources that make up
        this DIAObject.
    initialized : bool
        boolean specifying if the summary statistics for this DIAObject have
        been updated with the current set of DIASources in the SourceCatalog.
        This variable should be set to false whenever a the SourceCatalog
        of DIASources changes and set to true when the initialize method is
        run.
    """

    def __init__(self, source_catalog, object_source_record=None):
        """  Create a DIAObject given an input SourceCatalog of
        DIASources.

        Takes as input an lsst.afw.table.SourceCatalog object specifying a
        collection of DIASources that make up this DIAObject. The optional
        input object_source_record should contain summary statictics on the
        SourceCatalog of DIASources. Using this optional input escapes the
        need to recompute the summary statictis when not nessisary.

        Parameters
        ----------
        source_catalog : lsst.afw.table.SourceCatalog
            SourceCatalog of DIASource associated to this DIAObject
        object_source_record : lsst.afw.table.SourceRecord, optional
            Optional input SourceRecord containing summary statistics on
            the input SourceCatalog.

        Returns
        -------
        A DIAObject instance
        """

        self.source_catalog = source_catalog
        self.initialized = False

        if object_source_record is None:
            self.initialize()
        else:
            self.dia_object_record = object_source_record
            self.initialized = True

    def initilize(self):
        """ Compute all summary statistics given the current catalog of
        DIASources asigned to this DIAObject.

        Store these summaries (e.g. median RA/DEC position, fluxes...) in
        the object_source_record attribute and set the class variable
        intialized to True

        Returns
        -------
        None
        """

        self.initialized = False

        # compute all summary statistics on the catalog of DIASources.

        self.compute_summary_statistics()

        self.initialized = True

        return None

    def isInitialized(self):
        """ Return the current state of this DIAObject.

        If True is returned this DIAObject has computed all of it's summary
        statistics for the current collection of DIAObjects that make it up.

        Return
        ------
        bool
        """

        return self.initialized

    def append_dia_source(self, input_dia_source):
        """ Append the input_dia_source to the source_catalog attribute

        Additionally set initialized to False.

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
        self.initialized = False

        # Do stuff to append this SourceRecord.

        return None

    def compute_light_curve(self):
        """ Retreve the precomputed light curve of  fluxes for the DIASources
        that make up this DIAObject.

        Returns
        -------
        An array like object specifying the light curve for this object.
        """

        # Loop through DIASources and return the "light curve"
        # I'm keeping this separate from coupte summary statistics for the
        # moment.

        Raise(NotImplementedError)

    def compute_summary_statistics(self):
        """ Retrive properties from DIASourceCatalog attribute and update the
        summary statistics that represent this DIAObject

        Returns
        -------
        None
        """

        # Loop through DIASources, compute summary statistics (TBD) and store
        # them in dia_object_record attribute.

        Raise(NotImplementedError)

    def get_object_record(self):
        """ Retrive the SourceRecord that represents the summary statistics on
        this DIAObject's set of DIASources.

        Return
        ------
        A lsst.afw.table.SourceRecord
        """
        return self.dia_object_record

    def get_dia_source_catalog(self):
        """ Retrive the SourceCatalog that represents the DIASources that make
        up this DIAObject.

        A lsst.afw.table.SourceCatalog
        """
        return self.source_catalog
