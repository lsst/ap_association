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

"""Defines afw schemas and conversions for use in ap_association tasks.

Until a more finalized interface between the prompt products database (PPD) can
be established, we put many utility functions for converting `lsst.afw.table` and
`lsst.afw.image` objects to a PPD. This includes both mapping schemas between
different catalogs and the DB.
"""

__all__ = ["make_dia_object_schema",
           "make_dia_source_schema",
           "getCcdVisitSchemaSql"]

from collections import OrderedDict as oDict

import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
from lsst.daf.base import DateTime


def make_dia_object_schema(filter_names=None):
    """Define and create the minimal schema required for a DIAObject.

    Parameters
    ----------
    filter_names : `list` of `str`
        Names of the filters expect and compute means for.

    Returns
    -------
    schema : `lsst.afw.table.Schema`
        Minimal schema for DIAObjects.
    """
    schema = afwTable.SourceTable.makeMinimalSchema()
    # For the MVP/S we currently only care about the position though
    # in the future we will add summary computations for fluxes etc.
    # as well as their errors.

    # TODO: In the future we would like to store a covariance of the coordinate.
    # This functionality is not defined currently in the stack, so we will hold
    # off until it is implemented. This is to be addressed in DM-7101.
    schema.addField("pixelId", type='L',
                    doc='Unique spherical pixelization identifier.')
    schema.addField("nDiaSources", type='L',
                    doc="Number of total DiaSources associated with this "
                        "object.")
    if filter_names is None:
        filter_names = []
    for filter_name in filter_names:
        schema.addField("psFluxMean_%s" % filter_name, type='D',
                        doc='Calibrated mean flux in %s band.' % filter_name)
        schema.addField("psFluxMeanErr_%s" % filter_name, type='D',
                        doc='Calibrated error on mean flux in %s band.' %
                        filter_name)
        schema.addField("psFluxSigma_%s" % filter_name, type='D',
                        doc='Calibrated scatter in flux in %s band.' %
                        filter_name)
    return schema


def make_dia_source_schema():
    """ Define and create the minimal schema required for a DIASource.

    Returns
    -------
    schema : `lsst.afw.table.Schema`
        Minimal schema for DIAObjects.
    """
    schema = afwTable.SourceTable.makeMinimalSchema()
    schema.addField("diaObjectId", type='L',
                    doc='Unique identifier of the DIAObject this source is '
                        'associated to.')
    schema.addField("ccdVisitId", type='L',
                    doc='Id of the exposure and ccd this object was detected '
                        'in.')
    schema.addField("midPointTai", type="D",
                    doc="Time of mid-exposure for this DIASource")
    schema.addField("psFlux", type='D',
                    doc='Calibrated PSF flux of this source.')
    schema.addField("psFluxErr", type='D',
                    doc='Calibrated PSF flux err of this source.')
    schema.addField("filterName", type='String', size=10,
                    doc='String name of the filter this source was observed '
                        'in.')
    schema.addField("filterId", type='L',
                    doc='Obs package id of the filter this source was '
                        'observed in.')
    schema.addField("flags", type='L',
                    doc='Quality flags for this DIASource.')
    schema.addField("pixelId", type='L',
                    doc='Unique spherical pixelization identifier.')
    return schema


def getCcdVisitSchemaSql():
    """Define the schema for the CcdVisit table.

    Returns
    -------
    ccdVisitNames : `collections.OrderedDict`
       Names of columns in the ccdVisit table.
    """
    return oDict([("ccdVisitId", "INTEGER PRIMARY KEY"),
                  ("ccdNum", "INTEGER"),
                  ("filterName", "TEXT"),
                  ("filterId", "INTEGER"),
                  ("ra", "REAL"),
                  ("decl", "REAL"),
                  ("expTime", "REAL"),
                  ("expMidptMJD", "REAL"),
                  ("fluxZeroPoint", "REAL"),
                  ("fluxZeroPointErr", "REAL")])


def get_ccd_visit_info_from_exposure(exposure):
    """
    Extract info on the ccd and visit from the exposure.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure to store information from.

    Returns
    -------
    values : `dict` of ``values``
        List values representing info taken from the exposure.
    """
    visit_info = exposure.getInfo().getVisitInfo()
    date = visit_info.getDate()
    sphPoint = exposure.getWcs().getSkyOrigin()
    # TODO: Calib is going away and being replaced with photoCalib as in
    # DM-10153.
    flux0, flux0_err = exposure.getCalib().getFluxMag0()
    filter_obj = exposure.getFilter()
    # Values list is:
    # [CcdVisitId ``int``,
    #  ccdNum ``int``,
    #  filterName ``str``,
    #  RA WCS center ``degrees``,
    #  DEC WCS center ``degrees``,
    #  exposure time ``seconds``,
    #  dateTimeMJD ``days``,
    #  flux zero point ``counts``,
    #  flux zero point error ``counts``]
    values = {'ccdVisitId': visit_info.getExposureId(),
              'ccdNum': exposure.getDetector().getId(),
              'filterName': filter_obj.getName(),
              'filterId': filter_obj.getId(),
              'ra': sphPoint.getRa().asDegrees(),
              'decl': sphPoint.getDec().asDegrees(),
              'expTime': visit_info.getExposureTime(),
              'expMidptMJD': date.get(system=DateTime.MJD),
              'fluxZeroPoint': flux0,
              'fluxZeroPointErr': flux0_err,
              'photoCalib': afwImage.PhotoCalib(1 / flux0,
                                                flux0_err / flux0 ** 2)}
    return values
