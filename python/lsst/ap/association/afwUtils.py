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
           "make_dia_forced_source_schema",
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

    # Generated automatically from ppdb-schema.yaml in dax_ppdb/data.
    schema.addField('validityStart', type='L',
                    doc='Time when validity of this diaObject starts.')
    schema.addField('validityEnd', type='L',
                    doc='Time when validity of this diaObject ends.')
    schema.addField('raErr', type='Angle',
                    doc='Uncertainty of ra.')
    schema.addField('declErr', type='Angle',
                    doc='Uncertainty of decl.')
    schema.addField('ra_decl_Cov', type='F',
                    doc='Covariance between ra and decl.')
    schema.addField('radecTai', type='D',
                    doc='Time at which the object was at a position ra/decl.')
    schema.addField('pmRa', type='F',
                    doc='Proper motion (ra).')
    schema.addField('pmRaErr', type='F',
                    doc='Uncertainty of pmRa.')
    schema.addField('pmDecl', type='F',
                    doc='Proper motion (decl).')
    schema.addField('pmDeclErr', type='F',
                    doc='Uncertainty of pmDecl.')
    schema.addField('parallax', type='F',
                    doc='Parallax.')
    schema.addField('parallaxErr', type='F',
                    doc='Uncertainty of parallax.')
    schema.addField('pmRa_pmDecl_Cov', type='F',
                    doc='Covariance of pmRa and pmDecl.')
    schema.addField('pmRa_parallax_Cov', type='F',
                    doc='Covariance of pmRa and parallax.')
    schema.addField('pmDecl_parallax_Cov', type='F',
                    doc='Covariance of pmDecl and parallax.')
    schema.addField('pmParallaxLnL', type='F',
                    doc='Natural log of the likelihood of the linear proper motion parallax fit.')
    schema.addField('pmParallaxChi2', type='F',
                    doc='Chi^2 static of the model fit.')
    schema.addField('pmParallaxNdata', type='I',
                    doc='The number of data points used to fit the model.')
    schema.addField('uPSFluxMean', type='F',
                    doc='Weighted mean point-source model magnitude for u filter.')
    schema.addField('uPSFluxMeanErr', type='F',
                    doc='Standard error of uPSFluxMean.')
    schema.addField('uPSFluxSigma', type='F',
                    doc='Standard deviation of the distribution of uPSFlux.')
    schema.addField('uPSFluxChi2', type='F',
                    doc='Chi^2 statistic for the scatter of uPSFlux around uPSFluxMean.')
    schema.addField('uPSFluxNdata', type='I',
                    doc='The number of data points used to compute uPSFluxChi2.')
    schema.addField('uFPFluxMean', type='F',
                    doc='Weighted mean forced photometry flux for u filter.')
    schema.addField('uFPFluxMeanErr', type='F',
                    doc='Standard error of uFPFluxMean.')
    schema.addField('uFPFluxSigma', type='F',
                    doc='Standard deviation of the distribution of uFPFlux.')
    schema.addField('gPSFluxMean', type='F',
                    doc='Weighted mean point-source model magnitude for g filter.')
    schema.addField('gPSFluxMeanErr', type='F',
                    doc='Standard error of gPSFluxMean.')
    schema.addField('gPSFluxSigma', type='F',
                    doc='Standard deviation of the distribution of gPSFlux.')
    schema.addField('gPSFluxChi2', type='F',
                    doc='Chi^2 statistic for the scatter of gPSFlux around gPSFluxMean.')
    schema.addField('gPSFluxNdata', type='I',
                    doc='The number of data points used to compute gPSFluxChi2.')
    schema.addField('gFPFluxMean', type='F',
                    doc='Weighted mean forced photometry flux for g filter.')
    schema.addField('gFPFluxMeanErr', type='F',
                    doc='Standard error of gFPFluxMean.')
    schema.addField('gFPFluxSigma', type='F',
                    doc='Standard deviation of the distribution of gFPFlux.')
    schema.addField('rPSFluxMean', type='F',
                    doc='Weighted mean point-source model magnitude for r filter.')
    schema.addField('rPSFluxMeanErr', type='F',
                    doc='Standard error of rPSFluxMean.')
    schema.addField('rPSFluxSigma', type='F',
                    doc='Standard deviation of the distribution of rPSFlux.')
    schema.addField('rPSFluxChi2', type='F',
                    doc='Chi^2 statistic for the scatter of rPSFlux around rPSFluxMean.')
    schema.addField('rPSFluxNdata', type='I',
                    doc='The number of data points used to compute rPSFluxChi2.')
    schema.addField('rFPFluxMean', type='F',
                    doc='Weighted mean forced photometry flux for r filter.')
    schema.addField('rFPFluxMeanErr', type='F',
                    doc='Standard error of rFPFluxMean.')
    schema.addField('rFPFluxSigma', type='F',
                    doc='Standard deviation of the distribution of rFPFlux.')
    schema.addField('iPSFluxMean', type='F',
                    doc='Weighted mean point-source model magnitude for i filter.')
    schema.addField('iPSFluxMeanErr', type='F',
                    doc='Standard error of iPSFluxMean.')
    schema.addField('iPSFluxSigma', type='F',
                    doc='Standard deviation of the distribution of iPSFlux.')
    schema.addField('iPSFluxChi2', type='F',
                    doc='Chi^2 statistic for the scatter of iPSFlux around iPSFluxMean.')
    schema.addField('iPSFluxNdata', type='I',
                    doc='The number of data points used to compute iPSFluxChi2.')
    schema.addField('iFPFluxMean', type='F',
                    doc='Weighted mean forced photometry flux for i filter.')
    schema.addField('iFPFluxMeanErr', type='F',
                    doc='Standard error of iFPFluxMean.')
    schema.addField('iFPFluxSigma', type='F',
                    doc='Standard deviation of the distribution of iFPFlux.')
    schema.addField('zPSFluxMean', type='F',
                    doc='Weighted mean point-source model magnitude for z filter.')
    schema.addField('zPSFluxMeanErr', type='F',
                    doc='Standard error of zPSFluxMean.')
    schema.addField('zPSFluxSigma', type='F',
                    doc='Standard deviation of the distribution of zPSFlux.')
    schema.addField('zPSFluxChi2', type='F',
                    doc='Chi^2 statistic for the scatter of zPSFlux around zPSFluxMean.')
    schema.addField('zPSFluxNdata', type='I',
                    doc='The number of data points used to compute zPSFluxChi2.')
    schema.addField('zFPFluxMean', type='F',
                    doc='Weighted mean forced photometry flux for z filter.')
    schema.addField('zFPFluxMeanErr', type='F',
                    doc='Standard error of zFPFluxMean.')
    schema.addField('zFPFluxSigma', type='F',
                    doc='Standard deviation of the distribution of zFPFlux.')
    schema.addField('yPSFluxMean', type='F',
                    doc='Weighted mean point-source model magnitude for y filter.')
    schema.addField('yPSFluxMeanErr', type='F',
                    doc='Standard error of yPSFluxMean.')
    schema.addField('yPSFluxSigma', type='F',
                    doc='Standard deviation of the distribution of yPSFlux.')
    schema.addField('yPSFluxChi2', type='F',
                    doc='Chi^2 statistic for the scatter of yPSFlux around yPSFluxMean.')
    schema.addField('yPSFluxNdata', type='I',
                    doc='The number of data points used to compute yPSFluxChi2.')
    schema.addField('yFPFluxMean', type='F',
                    doc='Weighted mean forced photometry flux for y filter.')
    schema.addField('yFPFluxMeanErr', type='F',
                    doc='Standard error of yFPFluxMean.')
    schema.addField('yFPFluxSigma', type='F',
                    doc='Standard deviation of the distribution of yFPFlux.')
    # Mapping of arrays and BLOBs is not currently supported by the PPDB so we
    # these columns out.
    # schema.addField('uLcPeriodic', type='ArrayF',
    #                 doc='Periodic features extracted from light-curves using generalized Lomb-Scargle '
    #                     'periodogram for u filter. [32 FLOAT].')
    # schema.addField('gLcPeriodic', type='ArrayF',
    #                 doc='Periodic features extracted from light-curves using generalized Lomb-Scargle '
    #                     'periodogram for g filter. [32 FLOAT].')
    # schema.addField('rLcPeriodic', type='ArrayF',
    #                 doc='Periodic features extracted from light-curves using generalized Lomb-Scargle '
    #                     'periodogram for r filter. [32 FLOAT].')
    # schema.addField('iLcPeriodic', type='ArrayF',
    #                 doc='Periodic features extracted from light-curves using generalized Lomb-Scargle '
    #                     'periodogram for i filter. [32 FLOAT].')
    # schema.addField('zLcPeriodic', type='ArrayF',
    #                 doc='Periodic features extracted from light-curves using generalized Lomb-Scargle '
    #                     'periodogram for z filter. [32 FLOAT].')
    # schema.addField('yLcPeriodic', type='ArrayF',
    #                 doc='Periodic features extracted from light-curves using generalized Lomb-Scargle '
    #                     'periodogram for y filter. [32 FLOAT].')
    # schema.addField('uLcNonPeriodic', type='ArrayF',
    #                 doc='Non-periodic features extracted from light-curves using generalized Lomb-Scargle '
    #                     'periodogram for u filter. [20 FLOAT].')
    # schema.addField('gLcNonPeriodic', type='ArrayF',
    #                 doc='Non-periodic features extracted from light-curves using generalized Lomb-Scargle '
    #                     'periodogram for g filter. [20 FLOAT].')
    # schema.addField('rLcNonPeriodic', type='ArrayF',
    #                 doc='Non-periodic features extracted from light-curves using generalized Lomb-Scargle '
    #                     'periodogram for r filter. [20 FLOAT].')
    # schema.addField('iLcNonPeriodic', type='ArrayF',
    #                 doc='Non-periodic features extracted from light-curves using generalized Lomb-Scargle '
    #                     'periodogram for i filter. [20 FLOAT].')
    # schema.addField('zLcNonPeriodic', type='ArrayF',
    #                 doc='Non-periodic features extracted from light-curves using generalized Lomb-Scargle '
    #                     'periodogram for z filter. [20 FLOAT].')
    # schema.addField('yLcNonPeriodic', type='ArrayF',
    #                 doc='Non-periodic features extracted from light-curves using generalized Lomb-Scargle '
    #                     'periodogram for y filter. [20 FLOAT].')
    schema.addField('nearbyObj1', type='L',
                    doc='Id of the closest nearby object.')
    schema.addField('nearbyObj1Dist', type='F',
                    doc='Distance to nearbyObj1.')
    schema.addField('nearbyObj1LnP', type='F',
                    doc='Natural log of the probability that the observed diaObject is the same as the '
                        'nearbyObj1.')
    schema.addField('nearbyObj2', type='L',
                    doc='Id of the second-closest nearby object.')
    schema.addField('nearbyObj2Dist', type='F',
                    doc='Distance to nearbyObj2.')
    schema.addField('nearbyObj2LnP', type='F',
                    doc='Natural log of the probability that the observed diaObject is the same as the '
                        'nearbyObj2.')
    schema.addField('nearbyObj3', type='L',
                    doc='Id of the third-closest nearby object.')
    schema.addField('nearbyObj3Dist', type='F',
                    doc='Distance to nearbyObj3.')
    schema.addField('nearbyObj3LnP', type='F',
                    doc='Natural log of the probability that the observed diaObject is the same as the '
                        'nearbyObj3.')
    schema.addField('flags', type='L',
                    doc='Flags, bitwise OR tbd.')
    schema.addField('pixelId', type='L',
                    doc='HTM index.')
    schema.addField('lastNonForcedSource', type='L',
                    doc='Last time when non-forced DIASource was seen for this object')
    schema.addField('nDiaSources', type='I',
                    doc='Total number of DiaSources associated with this DiaObject.')

    return schema


def make_dia_source_schema():
    """ Define and create the minimal schema required for a DIASource.

    Returns
    -------
    schema : `lsst.afw.table.Schema`
        Minimal schema for DiaSources.
    """

    # Generated automatically from ppdb-schema.yaml in dax_ppdb/data.
    schema = afwTable.SourceTable.makeMinimalSchema()
    schema.addField('ccdVisitId', type='L',
                    doc='Id of the ccdVisit where this diaSource was measured. Note that we are allowing a '
                        'diaSource to belong to multiple amplifiers, but it may not span multiple ccds.')
    schema.addField('diaObjectId', type='L',
                    doc='Id of the diaObject this source was associated with, if any. If not, it is set to '
                        'NULL (each diaSource will be associated with either a diaObject or ssObject).')
    schema.addField('ssObjectId', type='L',
                    doc='Id of the ssObject this source was associated with, if any. If not, it is set to '
                        'NULL (each diaSource will be associated with either a diaObject or ssObject).')
    schema.addField('prv_procOrder', type='I',
                    doc='Position of this diaSource in the processing order relative to other diaSources '
                        'within a given diaObjectId or ssObjectId.')
    schema.addField('ssObjectReassocTime', type='D',
                    doc='Time when this diaSource was reassociated from diaObject to ssObject (if such '
                        'reassociation happens, otherwise NULL).')
    schema.addField('midPointTai', type='D',
                    doc='Effective mid-exposure time for this diaSource.')
    schema.addField('raErr', type='D',
                    doc='Uncertainty of ra.')
    schema.addField('declErr', type='D',
                    doc='Uncertainty of decl.')
    schema.addField('ra_decl_Cov', type='D',
                    doc='Covariance between ra and decl.')
    schema.addField('x', type='D',
                    doc='x position computed by a centroiding algorithm.')
    schema.addField('xErr', type='F',
                    doc='Uncertainty of x.')
    schema.addField('y', type='D',
                    doc='y position computed by a centroiding algorithm.')
    schema.addField('yErr', type='F',
                    doc='Uncertainty of y.')
    schema.addField('x_y_Cov', type='D',
                    doc='Covariance between x and y.')
    schema.addField('apFlux', type='D',
                    doc='Calibrated aperture flux. Note that this actually measures the difference between '
                        'the template and the visit image.')
    schema.addField('apFluxErr', type='D',
                    doc='Estimated uncertainty of apFlux.')
    schema.addField('snr', type='D',
                    doc='The signal-to-noise ratio at which this source was detected in the difference '
                        'image.')
    schema.addField('psFlux', type='D',
                    doc='Calibrated flux for Point Source model. Note this actually measures the flux '
                        'difference between the template and the visit image.')
    schema.addField('psFluxErr', type='D',
                    doc='Uncertainty of psFlux.')
    schema.addField('psRa', type='D',
                    doc=' RA-coordinate of centroid for point source model.')
    schema.addField('psRaErr', type='D',
                    doc='Uncertainty of psRa.')
    schema.addField('psDecl', type='D',
                    doc=' Decl-coordinate of centroid for point source model.')
    schema.addField('psDeclErr', type='D',
                    doc='Uncertainty of psDecl.')
    schema.addField('psFlux_psRa_Cov', type='D',
                    doc='Covariance between psFlux and psRa.')
    schema.addField('psFlux_psDecl_Cov', type='D',
                    doc='Covariance between psFlux and psDecl.')
    schema.addField('psRa_psDecl_Cov', type='D',
                    doc='Covariance between psRa and psDecl.')
    schema.addField('psLnL', type='D',
                    doc='Natural log likelihood of the observed data given the Point Source model.')
    schema.addField('psChi2', type='D',
                    doc='Chi^2 statistic of the model fit.')
    schema.addField('psNdata', type='I',
                    doc='The number of data points (pixels) used to fit the model.')
    schema.addField('trailFlux', type='D',
                    doc='Calibrated flux for a trailed source model. Note this actually measures the flux '
                        'difference between the template and the visit image.')
    schema.addField('trailFluxErr', type='D',
                    doc='Uncertainty of trailFlux.')
    schema.addField('trailRa', type='D',
                    doc=' RA-coordinate of centroid for trailed source model.')
    schema.addField('trailRaErr', type='D',
                    doc='Uncertainty of trailRa.')
    schema.addField('trailDecl', type='D',
                    doc=' Decl-coordinate of centroid for trailed source model.')
    schema.addField('trailDeclErr', type='D',
                    doc='Uncertainty of trailDecl.')
    schema.addField('trailLength', type='D',
                    doc='Maximum likelihood fit of trail length.')
    schema.addField('trailLengthErr', type='D',
                    doc='Uncertainty of trailLength.')
    schema.addField('trailAngle', type='D',
                    doc='Maximum likelihood fit of the angle between the meridian through the centroid and '
                        'the trail direction (bearing).')
    schema.addField('trailAngleErr', type='D',
                    doc='Uncertainty of trailAngle.')
    schema.addField('trailFlux_trailRa_Cov', type='D',
                    doc='Covariance of trailFlux and trailRa.')
    schema.addField('trailFlux_trailDecl_Cov', type='D',
                    doc='Covariance of trailFlux and trailDecl.')
    schema.addField('trailFlux_trailLength_Cov', type='D',
                    doc='Covariance of trailFlux and trailLength')
    schema.addField('trailFlux_trailAngle_Cov', type='D',
                    doc='Covariance of trailFlux and trailAngle')
    schema.addField('trailRa_trailDecl_Cov', type='D',
                    doc='Covariance of trailRa and trailDecl.')
    schema.addField('trailRa_trailLength_Cov', type='D',
                    doc='Covariance of trailRa and trailLength.')
    schema.addField('trailRa_trailAngle_Cov', type='D',
                    doc='Covariance of trailRa and trailAngle.')
    schema.addField('trailDecl_trailLength_Cov', type='D',
                    doc='Covariance of trailDecl and trailLength.')
    schema.addField('trailDecl_trailAngle_Cov', type='D',
                    doc='Covariance of trailDecl and trailAngle.')
    schema.addField('trailLength_trailAngle_Cov', type='D',
                    doc='Covariance of trailLength and trailAngle')
    schema.addField('trailLnL', type='D',
                    doc='Natural log likelihood of the observed data given the trailed source model.')
    schema.addField('trailChi2', type='D',
                    doc='Chi^2 statistic of the model fit.')
    schema.addField('trailNdata', type='I',
                    doc='The number of data points (pixels) used to fit the model.')
    schema.addField('dipMeanFlux', type='D',
                    doc='Maximum likelihood value for the mean absolute flux of the two lobes for a dipole '
                        'model.')
    schema.addField('dipMeanFluxErr', type='D',
                    doc='Uncertainty of dipMeanFlux.')
    schema.addField('dipFluxDiff', type='D',
                    doc='Maximum likelihood value for the difference of absolute fluxes of the two lobes for '
                        'a dipole model.')
    schema.addField('dipFluxDiffErr', type='D',
                    doc='Uncertainty of dipFluxDiff.')
    schema.addField('dipRa', type='D',
                    doc=' RA-coordinate of centroid for dipole model.')
    schema.addField('dipRaErr', type='D',
                    doc='Uncertainty of trailRa.')
    schema.addField('dipDecl', type='D',
                    doc=' Decl-coordinate of centroid for dipole model.')
    schema.addField('dipDeclErr', type='D',
                    doc='Uncertainty of dipDecl.')
    schema.addField('dipLength', type='D',
                    doc='Maximum likelihood value for the lobe separation in dipole model.')
    schema.addField('dipLengthErr', type='D',
                    doc='Uncertainty of dipLength.')
    schema.addField('dipAngle', type='D',
                    doc='Maximum likelihood fit of the angle between the meridian through the centroid and '
                        'the dipole direction (bearing, from negative to positive lobe).')
    schema.addField('dipAngleErr', type='D',
                    doc='Uncertainty of dipAngle.')
    schema.addField('dipMeanFlux_dipFluxDiff_Cov', type='D',
                    doc='Covariance of dipMeanFlux and dipFluxDiff.')
    schema.addField('dipMeanFlux_dipRa_Cov', type='D',
                    doc='Covariance of dipMeanFlux and dipRa.')
    schema.addField('dipMeanFlux_dipDecl_Cov', type='D',
                    doc='Covariance of dipMeanFlux and dipDecl.')
    schema.addField('dipMeanFlux_dipLength_Cov', type='D',
                    doc='Covariance of dipMeanFlux and dipLength.')
    schema.addField('dipMeanFlux_dipAngle_Cov', type='D',
                    doc='Covariance of dipMeanFlux and dipAngle.')
    schema.addField('dipFluxDiff_dipRa_Cov', type='D',
                    doc='Covariance of dipFluxDiff and dipRa.')
    schema.addField('dipFluxDiff_dipDecl_Cov', type='D',
                    doc='Covariance of dipFluxDiff and dipDecl.')
    schema.addField('dipFluxDiff_dipLength_Cov', type='D',
                    doc='Covariance of dipFluxDiff and dipLength.')
    schema.addField('dipFluxDiff_dipAngle_Cov', type='D',
                    doc='Covariance of dipFluxDiff and dipAngle.')
    schema.addField('dipRa_dipDecl_Cov', type='D',
                    doc='Covariance of dipRa and dipDecl.')
    schema.addField('dipRa_dipLength_Cov', type='D',
                    doc='Covariance of dipRa and dipLength.')
    schema.addField('dipRa_dipAngle_Cov', type='D',
                    doc='Covariance of dipRa and dipAngle.')
    schema.addField('dipDecl_dipLength_Cov', type='D',
                    doc='Covariance of dipDecl and dipLength.')
    schema.addField('dipDecl_dipAngle_Cov', type='D',
                    doc='Covariance of dipDecl and dipAngle.')
    schema.addField('dipLength_dipAngle_Cov', type='D',
                    doc='Covariance of dipLength and dipAngle.')
    schema.addField('dipLnL', type='D',
                    doc='Natural log likelihood of the observed data given the dipole source model.')
    schema.addField('dipChi2', type='D',
                    doc='Chi^2 statistic of the model fit.')
    schema.addField('dipNdata', type='I',
                    doc='The number of data points (pixels) used to fit the model.')
    schema.addField('totFlux', type='D',
                    doc='Calibrated flux for Point Source model measured on the visit image centered at the '
                        'centroid measured on the difference image (forced photometry flux).')
    schema.addField('totFluxErr', type='D',
                    doc='Estimated uncertainty of totFlux.')
    schema.addField('diffFlux', type='D',
                    doc='Calibrated flux for Point Source model centered on radec but measured on the '
                        'difference of snaps comprising this visit.')
    schema.addField('diffFluxErr', type='D',
                    doc='Estimated uncertainty of diffFlux.')
    schema.addField('fpBkgd', type='D',
                    doc='Estimated sky background at the position (centroid) of the object.')
    schema.addField('fpBkgdErr', type='D',
                    doc='Estimated uncertainty of fpBkgd.')
    schema.addField('ixx', type='D',
                    doc='Adaptive second moment of the source intensity.')
    schema.addField('ixxErr', type='D',
                    doc='Uncertainty of ixx.')
    schema.addField('iyy', type='D',
                    doc='Adaptive second moment of the source intensity.')
    schema.addField('iyyErr', type='D',
                    doc='Uncertainty of iyy.')
    schema.addField('ixy', type='D',
                    doc='Adaptive second moment of the source intensity.')
    schema.addField('ixyErr', type='D',
                    doc='Uncertainty of ixy.')
    schema.addField('ixx_iyy_Cov', type='D',
                    doc='Covariance of ixx and iyy.')
    schema.addField('ixx_ixy_Cov', type='D',
                    doc='Covariance of ixx and ixy.')
    schema.addField('iyy_ixy_Cov', type='D',
                    doc='Covariance of iyy and ixy.')
    schema.addField('ixxPSF', type='D',
                    doc='Adaptive second moment for the PSF.')
    schema.addField('iyyPSF', type='D',
                    doc='Adaptive second moment for the PSF.')
    schema.addField('ixyPSF', type='D',
                    doc='Adaptive second moment for the PSF.')
    schema.addField('extendedness', type='D',
                    doc='A measure of extendedness, Computed using a combination of available moments and '
                        'model fluxes or from a likelihood ratio of point/trailed source models (exact '
                        'algorithm TBD). extendedness = 1 implies a high degree of confidence that the '
                        'source is extended. extendedness = 0 implies a high degree of confidence that the '
                        'source is point-like.')
    schema.addField('spuriousness', type='D',
                    doc='A measure of spuriousness, computed using information from the source and image '
                        'characterization, as well as the information on the Telescope and Camera system '
                        '(e.g., ghost maps, defect maps, etc.).')
    schema.addField('flags', type='L',
                    doc='Flags, bitwise OR tbd.')
    schema.addField('pixelId', type='L',
                    doc='HTM index.')
    schema.addField("filterName", type='String', size=10,
                    doc='String name of the filter this source was observed '
                        'in.')
    schema.addField("filterId", type='L',
                    doc='Obs package id of the filter this source was '
                        'observed in.')
    return schema


def make_dia_forced_source_schema():
    """ Define and create the minimal schema required for a DiaForcedSource.

    Returns
    -------
    schema : `lsst.afw.table.Schema`
        Minimal schema for DiaForcedSources.
    """

    # Generated automatically from ppdb-schema.yaml in dax_ppdb/data.
    schema = afwTable.SourceTable.makeMinimalSchema()
    schema.addField('ccdVisitId', type='L',
                    doc='Id of the ccdVisit where this diaSource was measured. Note that we are allowing a '
                        'diaSource to belong to multiple amplifiers, but it may not span multiple ccds.')
    schema.addField('psFlux', type='D',
                    doc='Calibrated flux for Point Source model. Note this actually measures the flux '
                        'difference between the template and the visit image.')
    schema.addField('psFluxErr', type='D',
                    doc='Uncertainty of psFlux.')
    schema.addField('totFlux', type='D',
                    doc='Calibrated flux measured in direct image.')
    schema.addField('totFluxErr', type='D',
                    doc='Uncertainty of totFlux.')
    schema.addField('x', type='D',
                    doc='x position at which psFlux has been measured.')
    schema.addField('y', type='D',
                    doc='y position at which psFlux has been measured.')
    schema.addField('flags', type='L',
                    doc='Flags from measurement on the difference image, '
                        'bitwise OR tbd')
    schema.addField('flagsDirectIm', type='L',
                    doc='Flags from measurement on the direct image, bitwise '
                        'OR tbd')
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
    # TODO: need to scale these until DM-10153 is completed and PhotoCalib has replaced Calib entirely
    referenceFlux = 1e23 * 10**(48.6 / -2.5) * 1e9
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
              'photoCalib': afwImage.PhotoCalib(referenceFlux / flux0,
                                                referenceFlux*flux0_err / flux0 ** 2)}
    return values
