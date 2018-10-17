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
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
from lsst.daf.base import DateTime


def make_dia_object_schema(
        filter_names=['u', 'g', 'r', 'i', 'z', 'y']):
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
    # schema.addField("raSigma", type="Angle",
    #                 doc='Scatter of DiaObject position in ra.')
    # schema.addField("declSigma", type="Angle",
    #                 doc='Scatter of DiaObject position in dec.')
    schema.addField("ra_decl_Cov", type="F",
                    doc='Covariance of ra and dec.')
    schema.addField("radecTai", type="D",
                    doc="MJD time DiaObject was at position ra/dec.")
    schema.addField("pmRa", type="F",
                    doc='Proper motion in ra.')
    schema.addField("pmDec", type="F",
                    doc='Proper motion in dec.')
    if filter_names is None:
        filter_names = []
    for filter_name in filter_names:
        schema.addField("%sPSFluxMean" % filter_name, type='F',
                        doc='Calibrated mean flux in %s band.' % filter_name)
        schema.addField("%sPSFluxMeanErr" % filter_name, type='F',
                        doc='Calibrated error on mean flux in %s band.' %
                        filter_name)
        schema.addField("%sPSFluxSigma" % filter_name, type='F',
                        doc='Calibrated scatter in flux in %s band.' %
                        filter_name)
        schema.addField("%sPSFluxChi2" % filter_name, type="F",
                        doc='Chi2 statistic for the scatter of PSFFlux around '
                            'PSFluxMean.')
        schema.addField("%sPSFluxNdata" % filter_name, type="F",
                        doc='Calibrated mean flux in %s band.' % filter_name)
        schema.addField("%sFSFluxMean" % filter_name, type="F",
                        doc='Calibrated mean flux in %s band.' % filter_name)
        schema.addField("%sFSFluxMeanErr" % filter_name, type="F",
                        doc='Calibrated error on mean flux in %s band.' %
                        filter_name)
        schema.addField("%sFSFluxSigma" % filter_name, type="F",
                        doc='Calibrated scatter in flux in %s band.' %
                        filter_name)
        """
        schema.addField("%sLcPeriodic" % filter_name, type="ArrayF", size=32,
                        doc="Periodic features extracted from light-curves "
                        "using generalized Lomb-Scargle periodogram.")
        schema.addField("%sLcNonPeriodic" % filter_name, type="ArrayF",
                        size=20, doc="Non-periodic features extracted from "
                        "light-curves using generalized Lomb-Scargle "
                        "periodogram.")
        """
    for obj_idx in range(6):
        schema.addField("nearbyObj%i" % obj_idx, type="L",
                        doc="Id of the %ith closest object." % obj_idx)
        schema.addField("nearbyObj%iDist" % obj_idx, type="F",
                        doc="Distance to the %ith closest object." % obj_idx)
        schema.addField("nearbyObj%iDistLnP" % obj_idx, type="F",
                        doc="Natural log of the probability that the observed "
                        "diaObject is the same as the nearbyObj%i." % obj_idx)
    schema.addField("flags", type="L",
                    doc="Flags, bitwise.")
    schema.addField("pixelId", type='L',
                    doc='Unique spherical pixelization identifier.')
    return schema


def make_dia_source_schema():
    """ Define and create the minimal schema required for a DIASource.

    Returns
    -------
    schema : `lsst.afw.table.Schema`
        Minimal schema for DIAObjects.
    """
    schema = afwTable.SourceTable.makeMinimalSchema()
    schema.addField("ccdVisitId", type='L',
                    doc='Id of the exposure and ccd this object was detected '
                    'in.')
    schema.addField("diaObjectId", type='L',
                    doc='Unique identifier of the DIAObject this source is '
                    'associated to.')
    schema.addField("ssObjectId", type='L',
                    doc='Unique identifier of the SolarSystemObjects this '
                    'source is associated to.')
    schema.addField("midPointTai", type="D",
                    doc="Midpoint of time of the exposure.")
    schema.addField("x", type="F",
                    doc="x position of source centroid in the diffim.")
    schema.addField("y", type="F",
                    doc="y position of source centroid in the diffim.")
    schema.addField("xErr", type="F",
                    doc="x position of source centroid in the diffim.")
    schema.addField("yErr", type="F",
                    doc="y position error of source centroid in the diffim.")
    schema.addField("xyCov", type="F",
                    doc="xy Covariance of source centroid in the diffim.")
    schema.addField("apFlux", type="F",
                    doc='Calibrated PSF flux of this source.')
    schema.addField("apFluxErr", type="F",
                    doc='Calibrated PSF flux err of this source.')
    schema.addField("snr", type="F",
                    doc='Signal to noise ratio of ')
    schema.addField("psFlux", type="F",
                    doc='Calibrated PSF flux of this source.')
    schema.addField("psFluxErr", type="F",
                    doc='Calibrated PSF flux err of this source.')
    """ 
    schema.addField("psRa", type="Angle",
                    doc="RA of point source model.")
    schema.addField("psDec", type="Angle",
                    doc="DEC of point source model.")
    schema.addField("psRaErr", type="Angle",
                    doc="Error on RA of point source model.")
    schema.addField("psDecErr", type="Angle",
                    doc="Error on DEC of point source model.")
    schema.addField("psRaDecCov", type="Angle",
                    doc="Covariance of RA/DEC in the point source model.")
    """
    schema.addField("psLnL", type="F",
                    doc="Log likelihood of the point source model.")
    schema.addField("psChi2", type="F",
                    doc="Chi2 value of the point source model fit.")
    schema.addField("psNdata", type="F",
                    doc="Number of pixels used in the point source model fit.")
    schema.addField("trailFlux", type="F",
                    doc="Flux of a trailed source model.")
    """
    schema.addField("trailRa", type="Angle",
                    doc="RA of the trailed flux model.")
    schema.addField("trailDec", type="Angle",
                    doc="DEC of the trailed flux model.")
    """
    schema.addField("trailLength", type="F",
                    doc="Best fit trail length of the source.")
    """
    schema.addField("trailAngle", type="F",
                    doc="Angle of the trail relative to the meridian through "
                    "the centroid.")
    schema.addField("trailCov", type="ArrayF", size=15,
                    doc="Covariance of the trailed model fit parameters.")
    """
    schema.addField("trailLnL", type="F",
                    doc="Log likelihood of the trailed source model.")
    schema.addField("trailChi2", type="F",
                    doc="Chi2 of the trailed source model.")
    schema.addField("trailNdata", type="F",
                    doc="Number of pixels used to fit the model.")
    schema.addField("dipMeanFlux", type="F",
                    doc=" Maximum likelihood value for the mean absolute flux "
                    "of the two lobes for a dipole model.")
    schema.addField("dipFluxDiff", type="F",
                    doc="Maximum likelihood value for the difference of "
                    "absolute fluxes of the two lobes for a dipole model.")
    """
    schema.addField("dipRa", type="Angle",
                    doc="RA of the dipole flux model.")
    schema.addField("dipDec", type="Angle",
                    doc="DEC of the dipole flux model.")
    """
    schema.addField("dipLength", type="F",
                    doc="Best fit length between the poles.")
    """
    schema.addField("dipAngle", type="F",
                    doc="Angle of the dipole relative to the meridian through "
                    "the centroid.")
    """
    # schema.addField("dipCov", type="ArrayF", size=21,
    #                 doc="Covariance of the dipole model fit parameters.")
    schema.addField("dipLnL", type="F",
                    doc="Log likelihood of the dipole source model.")
    schema.addField("dipChi2", type="F",
                    doc="Chi2 of the dipole source model.")
    schema.addField("dipNdata", type="F",
                    doc="Number of pixels used to fit the model.")
    schema.addField("totFlux", type="F",
                    doc="Point source calibrated flux in the visit image, "
                    "measured at the location of the source.")
    schema.addField("totFluxErr", type="F",
                    doc="Error on totFlux.")
    schema.addField("diffFlux", type="F",
                    doc="Calibrated flux in the difference of the snaps at "
                    "the position of the detected source.")
    schema.addField("diffFluxErr", type="F",
                    doc="Error on diffFlux.")
    schema.addField("fpBkgd", type="F",
                    doc="Estimated background at the centroid position in the "
                    "template image.")
    schema.addField("fpBkgdErr", type="F",
                    doc="Estimated uncertainty of fpBkgd.")
    schema.addField("ixx", type="F",
                    doc="Adaptive second moment for the source.")
    schema.addField("iyy", type="F",
                    doc="Adaptive second moment for the source.")
    schema.addField("ixy", type="F",
                    doc="Adaptive second moment for the source.")
    # schema.addField("icov", type="ArrayF", size=6,
    #                 doc="Covariance of the second moments.")
    schema.addField("ixxPSF", type="F",
                    doc="Adaptive second moment for the PSF.")
    schema.addField("iyyPSF", type="F",
                    doc="Adaptive second moment for the PSF.")
    schema.addField("ixyPSF", type="F",
                    doc="Adaptive second moment for the PSF.")
    schema.addField("extendedness", type="F",
                    doc="Measure of how extended the source is.")
    schema.addField("spuriousness", type="F",
                    doc="Measure of spuriousness given telescope defect and "
                    "ghosting maps.")
    schema.addField("flags", type='L',
                    doc='Quality flags for this DIASource.')
    schema.addField("filterName", type='String', size=10,
                    doc='String name of the filter this source was observed '
                        'in.')
    schema.addField("filterId", type='I',
                    doc='Obs package id of the filter this source was '
                        'observed in.')
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
