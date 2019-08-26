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

"""A simple implementation of source association task for ap_verify.
"""

__all__ = ["AssociationConfig", "AssociationTask"]

from astropy.stats import median_absolute_deviation
import numpy as np
import pandas
from scipy.optimize import lsq_linear
from scipy.spatial import cKDTree
from scipy.stats import skew

import lsst.geom as geom
from lsst.daf.base import DateTime
from lsst.meas.algorithms.indexerRegistry import IndexerRegistry
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from .afwUtils import make_dia_object_schema

# Enforce an error for unsafe column/array value setting in pandas.
pandas.options.mode.chained_assignment = 'raise'


def _set_mean_position(dia_object_record, dia_sources):
    """Compute and set the mean position of the input dia_object_record using
    the positions of the input catalog of DIASources.

    Parameters
    ----------
    dia_object_record : `dict` or `pandas.Series`
        SourceRecord of the DIAObject to edit.
    dia_sources : `pandas.DataFrame`
        Catalog of DIASources to compute a mean position from.
    """
    coord_list = [geom.SpherePoint(src["ra"], src["decl"], geom.degrees)
                  for idx, src in dia_sources.iterrows()]
    ave_coord = geom.averageSpherePoint(coord_list)
    dia_object_record["ra"] = ave_coord.getRa().asDegrees()
    dia_object_record["decl"] = ave_coord.getDec().asDegrees()


def _set_flux_stats(dia_object_record, dia_sources, filter_name, filter_id):
    """Compute the mean, standard error, and variance of a DIAObject for
    a given band.

    Parameters
    ----------
    dia_object_record : `dict` or `pandas.Series`
        SourceRecord of the DIAObject to edit.
    dia_sources : `pandas.DataFrame`
        Catalog of DIASources to compute a mean position from.
    filter_name : `str`
        Name of the band pass filter to update.
    filter_id : `int`
        id of the filter in the AssociationDB.
    """
    currentFluxMask = dia_sources.loc[:, "filterId"] == filter_id
    fluxes = dia_sources.loc[currentFluxMask, "psFlux"].replace([None], np.nan)
    fluxErrors = dia_sources.loc[currentFluxMask, "psFluxErr"].replace([None], np.nan)
    midpointTais = dia_sources.loc[currentFluxMask, "midPointTai"]

    noNanMask = np.logical_and(np.isfinite(fluxes), np.isfinite(fluxErrors))
    fluxes = fluxes[noNanMask]
    fluxErrors = fluxErrors[noNanMask]
    midpointTais = midpointTais[noNanMask]

    if len(fluxes) == 1:
        dia_object_record['%sPSFluxMean' % filter_name] = fluxes.iloc[0]
        dia_object_record['%sPSFluxNdata' % filter_name] = 1
    elif len(fluxes > 1):
        fluxMean = np.average(fluxes, weights=1 / fluxErrors ** 2)

        # Standard, DDPD defined columns.
        dia_object_record['%sPSFluxMean' % filter_name] = fluxMean
        dia_object_record['%sPSFluxMeanErr' % filter_name] = np.sqrt(
            1 / np.sum(1 / fluxErrors ** 2))
        dia_object_record['%sPSFluxSigma' % filter_name] = np.std(fluxes,
                                                                  ddof=1)
        dia_object_record["%sPSFluxChi2" % filter_name] = np.sum(
            ((fluxMean - fluxes) / fluxErrors) ** 2)
        dia_object_record['%sPSFluxNdata' % filter_name] = len(fluxes)

        # Columns below are created in DM-18316 for use in ap_pipe/verify
        # testing.
        ptiles = np.percentile(fluxes, [5, 25, 50, 75, 95])
        dia_object_record['%sPSFluxPercentile05' % filter_name] = ptiles[0]
        dia_object_record['%sPSFluxPercentile25' % filter_name] = ptiles[1]
        dia_object_record['%sPSFluxMedian' % filter_name] = ptiles[2]
        dia_object_record['%sPSFluxPercentile75' % filter_name] = ptiles[3]
        dia_object_record['%sPSFluxPercentile95' % filter_name] = ptiles[4]

        dia_object_record['%sPSFluxMAD' % filter_name] = \
            median_absolute_deviation(fluxes)

        dia_object_record['%sPSFluxSkew' % filter_name] = skew(fluxes)

        dia_object_record['%sPSFluxMin' % filter_name] = fluxes.min()
        dia_object_record['%sPSFluxMax' % filter_name] = fluxes.max()

        deltaFluxes = fluxes[1:] - fluxes[:-1]
        deltaTimes = midpointTais[1:] - midpointTais[:-1]
        dia_object_record['%sPSFluxMaxSlope' % filter_name] = np.max(
            deltaFluxes / deltaTimes)

        m, b = _fit_linear_flux_model(fluxes, fluxErrors, midpointTais)
        dia_object_record['%sPSFluxLinearSlope' % filter_name] = m
        dia_object_record['%sPSFluxLinearIntercept' % filter_name] = b

        dia_object_record['%sPSFluxStetsonJ' % filter_name] = _stetson_J(
            fluxes, fluxErrors)

        dia_object_record['%sPSFluxErrMean' % filter_name] = \
            np.mean(fluxErrors)

    totFluxes = dia_sources.loc[:, "totFlux"].replace([None], np.nan)
    totFluxes = totFluxes[currentFluxMask]
    totFluxErrors = dia_sources.loc[:, "totFluxErr"].replace([None], value=np.nan)
    totFluxErrors = totFluxErrors[currentFluxMask]
    noNanMask = np.logical_and(np.isfinite(totFluxes),
                               np.isfinite(totFluxErrors))
    totFluxes = totFluxes[noNanMask]
    totFluxErrors = totFluxErrors[noNanMask]

    if len(totFluxes) == 1:
        dia_object_record['%sTOTFluxMean' % filter_name] = totFluxes.iloc[0]
    elif len(totFluxes) > 1:
        fluxMean = np.average(totFluxes, weights=1 / totFluxErrors ** 2)
        dia_object_record['%sTOTFluxMean' % filter_name] = fluxMean
        dia_object_record['%sTOTFluxMeanErr' % filter_name] = np.sqrt(
            1 / np.sum(1 / totFluxErrors ** 2))
        dia_object_record['%sTOTFluxSigma' % filter_name] = np.std(totFluxes,
                                                                   ddof=1)


def _fit_linear_flux_model(fluxes, errors, times):
    """Fit a linear model (m*x + b) to flux vs time.

    Parameters
    ----------
    fluxes : `numpy.ndarray`, (N,)
        Input fluxes.
    errors : `numpy.ndarray`, (N,)
        Import errors associated with fluxes.
    times : `numpy.ndarray`, (N,)
        Time of the flux observation.

    Returns
    -------
    ans : tuple, (2,)
        Slope (m) and intercept (b) values fit to the light-curve.
    """
    A = np.array([times / errors, 1 / errors]).transpose()
    ans = lsq_linear(A, fluxes / errors).x

    return ans


def _stetson_J(fluxes, errors):
    """Compute the single band stetsonJ statistic.

    Parameters
    ----------
    fluxes : `numpy.ndarray` (N,)
        Calibrated lightcurve flux values.
    errors : `numpy.ndarray`
        Errors on the calibrated lightcurve fluxes.

    Returns
    -------
    stetsonJ : `float`
        stetsonJ statistic for the input fluxes and errors.

    References
    ----------
    .. [1] Stetson, P. B., "On the Automatic Determination of Light-Curve
       Parameters for Cepheid Variables", PASP, 108, 851S, 1996
    """
    n_points = len(fluxes)
    flux_mean = _stetson_mean(fluxes, errors)
    delta_val = (
        np.sqrt(n_points / (n_points - 1)) * (fluxes - flux_mean) / errors)
    p_k = delta_val ** 2 - 1

    return np.mean(np.sign(p_k) * np.sqrt(np.fabs(p_k)))


def _stetson_mean(values, errors, alpha=2., beta=2., n_iter=20, tol=1e-6):
    """Compute the stetson mean of the fluxes which down-weights outliers.

    Weighted biased on an error weighted difference scaled by a constant
    (1/``a``) and raised to the power beta. Higher betas more harshly penalize
    outliers and ``a`` sets the number of sigma where a weighted difference of
    1 occurs.

    Parameters
    ----------
    values : `numpy.dnarray`, (N,)
        Input values to compute the mean of.
    errors : `numpy.ndarray`, (N,)
        Errors on the input values.
    alpha : `float`
        Scalar downweighting of the fractional difference. lower->more clipping
    beta : `float`
        Power law slope of the used to down-weight outliers. higher->more
        clipping
    n_iter : `int`
        Number of iterations of clipping.
    tol : `float`
        Fractional and absolute tolerance goal on the change in the mean before
        exiting early.

    Returns
    -------
    wmean : `float`
        Weighted stetson mean result.

    References
    ----------
    .. [1] Stetson, P. B., "On the Automatic Determination of Light-Curve
       Parameters for Cepheid Variables", PASP, 108, 851S, 1996
    """
    n_points = len(values)
    n_factor = np.sqrt(n_points / (n_points - 1))

    wmean = np.average(values, weights=1 / errors ** 2)
    for iter_idx in range(n_iter):
        chi = np.fabs(n_factor * (values - wmean) / errors)
        weights = 1 / (1 + (chi / alpha) ** beta)
        tmp_wmean = np.average(values, weights=weights)
        diff = np.fabs(tmp_wmean - wmean)
        wmean = tmp_wmean
        if diff / wmean < tol and diff < tol:
            break
    return wmean


class AssociationConfig(pexConfig.Config):
    """Config class for AssociationTask.
    """
    maxDistArcSeconds = pexConfig.Field(
        dtype=float,
        doc='Maximum distance in arcseconds to test for a DIASource to be a '
        'match to a DIAObject.',
        default=1.0,
    )
    indexer = IndexerRegistry.makeField(
        doc='Select the spatial indexer to use within the database.',
        default='HTM'
    )


class AssociationTask(pipeBase.Task):
    """Associate DIAOSources into existing DIAObjects.

    This task performs the association of detected DIASources in a visit
    with the previous DIAObjects detected over time. It also creates new
    DIAObjects out of DIASources that cannot be associated with previously
    detected DIAObjects.
    """

    ConfigClass = AssociationConfig
    _DefaultName = "association"

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.indexer = IndexerRegistry[self.config.indexer.name](
            self.config.indexer.active)
        self.dia_object_schema = make_dia_object_schema()

    @pipeBase.timeMethod
    def run(self, dia_sources, exposure, ppdb):
        """Load DIAObjects from the database, associate the sources, and
        persist the results into the L1 database.

        Parameters
        ----------
        dia_sources : `pandas.DataFrame`
            DIASources to be associated with existing DIAObjects.
        exposure : `lsst.afw.image`
            Input exposure representing the region of the sky the dia_sources
            were detected on. Should contain both the solved WCS and a bounding
            box of the ccd.
        ppdb : `lsst.dax.ppdb.Ppdb`
            Ppdb connection object to retrieve DIASources/Objects from and
            write to.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components.

            - ``dia_objects`` : Complete set of dia_objects covering the input
              exposure. Catalog contains newly created, updated, and untouched
              diaObjects. (`pandas.DataFrame`)
        """
        # Assure we have a Box2D and can use the getCenter method.
        dia_objects = self.retrieve_dia_objects(exposure, ppdb)

        dia_sources = self.check_dia_souce_radec(dia_sources)

        updated_obj_ids = self.associate_sources(dia_objects, dia_sources)

        # Store newly associated DIASources.
        ppdb.storeDiaSources(dia_sources)
        # Update previously existing DIAObjects with the information from their
        # newly association DIASources and create new DIAObjects from
        # unassociated sources.
        dia_objects = self.update_dia_objects(dia_objects,
                                              updated_obj_ids,
                                              exposure,
                                              ppdb)

        return pipeBase.Struct(
            dia_objects=dia_objects,
        )

    @pipeBase.timeMethod
    def retrieve_dia_objects(self, exposure, ppdb):
        """Convert the exposure object into HTM pixels and retrieve DIAObjects
        contained within the exposure.

        DiaObject DataFrame will be indexed on ``diaObjectId``.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            An exposure specifying a bounding region with a WCS to load
            DIAOjbects within.
        ppdb : `lsst.dax.ppdb.Ppdb`
            Ppdb connection object to retrieve DIAObjects from.

        Returns
        -------
        diaObjects : `pandas.DataFrame`
            DiaObjects within the exposure boundary.
        """
        bbox = geom.Box2D(exposure.getBBox())
        wcs = exposure.getWcs()

        ctr_coord = wcs.pixelToSky(bbox.getCenter())
        max_radius = max(
            ctr_coord.separation(wcs.pixelToSky(pp))
            for pp in bbox.getCorners())

        indexer_indices, on_boundry = self.indexer.getShardIds(
            ctr_coord, max_radius)
        # Index types must be cast to int to work with dax_ppdb.
        index_ranges = [[int(indexer_idx), int(indexer_idx) + 1]
                        for indexer_idx in indexer_indices]
        covering_dia_objects = ppdb.getDiaObjects(index_ranges,
                                                  return_pandas=True)
        ccd_mask = pandas.Series(False, index=covering_dia_objects.index)

        for df_idx, cov_dia_object in covering_dia_objects.iterrows():
            if self._check_dia_object_position(cov_dia_object, bbox, wcs):
                ccd_mask.loc[df_idx] = True

        diaObjects = covering_dia_objects[ccd_mask]
        diaObjects.set_index("diaObjectId", inplace=True)

        return diaObjects

    def _check_dia_object_position(self, dia_object_record, bbox, wcs):
        """Check the RA, DEC position of the current dia_object_record against
        the bounding box of the exposure.

        Parameters
        ----------
        dia_object_record : `pandas.Series`
            A SourceRecord object containing the DIAObject we would like to
            test against our bounding box.
        bbox : `lsst.geom.Box2D`
            Bounding box of exposure.
        wcs : `lsst.afw.geom.SkyWcs`
            WCS of exposure.

        Return
        ------
        is_contained : `bool`
            Object position is contained within the bounding box.
        """
        point = wcs.skyToPixel(geom.SpherePoint(dia_object_record["ra"],
                                                dia_object_record["decl"],
                                                geom.degrees))
        return bbox.contains(point)

    def check_dia_souce_radec(self, dia_sources):
        """Check that all DiaSources have non-NaN values for RA/DEC.

        If one or more DiaSources are found to have NaN values, throw a
        warning to the log with the ids of the offending sources. Drop them
        from the table.

        Parameters
        ----------
        dia_sources : `pandas.DataFrame`
            Input DiaSources to check for NaN values.

        Returns
        -------
        trimmed_sources : `pandas.DataFrame`
            DataFrame of DiaSources trimmed of all entries with NaN values for
            RA/DEC.
        """
        nan_mask = (dia_sources.loc[:, "ra"].isnull() |
                    dia_sources.loc[:, "decl"].isnull())
        if np.any(nan_mask):
            nan_idxs = np.argwhere(nan_mask).flatten()
            for nan_idx in nan_idxs:
                self.log.warning(
                    "DiaSource %i has NaN value for RA/DEC, "
                    "dropping from association." %
                    dia_sources.loc[nan_idx, "diaSourceId"])
            dia_sources = dia_sources[~nan_mask]
        return dia_sources

    @pipeBase.timeMethod
    def associate_sources(self, dia_objects, dia_sources):
        """Associate the input DIASources with the catalog of DIAObjects.

        DiaObject DataFrame must be indexed on ``diaObjectId``.

        Parameters
        ----------
        dia_objects : `pandas.DataFrame`
            Catalog of DIAObjects to attempt to associate the input
            DIASources into.
        dia_sources : `pandas.DataFrame`
            DIASources to associate into the DIAObjectCollection.

        Returns
        -------
        updated_ids : array-like of `int`
            Ids of the DIAObjects that the DIASources associated to including
            the ids of newly created DIAObjects.
        """

        scores = self.score(
            dia_objects, dia_sources,
            self.config.maxDistArcSeconds * geom.arcseconds)
        match_result = self.match(dia_objects, dia_sources, scores)

        self._add_association_meta_data(match_result)

        return match_result.associated_dia_object_ids

    @pipeBase.timeMethod
    def update_dia_objects(self, dia_objects, updated_obj_ids, exposure, ppdb):
        """Update select dia_objects currently stored within the database or
        create new ones.

        Modify the dia_object catalog in place to post-pend newly created
        DiaObjects.

        Parameters
        ----------
        dia_objects : `pandas.DataFrame`
            Pre-existing/loaded DIAObjects to copy values that are not updated
            from.
        updated_obj_ids : array-like of `int`
            Ids of the dia_objects that should be updated.
        exposure : `lsst.afw.image.Exposure`
            Input exposure representing the region of the sky the dia_sources
            were detected on. Should contain both the solved WCS and a bounding
            box of the ccd.
        ppdb : `lsst.dax.ppdb.Ppdb`
            Ppdb connection object to retrieve DIASources from and
            write DIAObjects to.

        Returns
        -------
        outputDiaObjects : `pandas.DataFrame`
            Union of updated and un-touched DiaObjects indexed on
            ``diaObjectId``.
        """
        filter_name = exposure.getFilter().getName()
        filter_id = exposure.getFilter().getId()

        updated_obj_ids.sort()
        dia_object_used = pandas.DataFrame(
            False,
            index=dia_objects.index,
            columns=["used"])

        dateTime = exposure.getInfo().getVisitInfo().getDate()

        dia_sources = ppdb.getDiaSources(updated_obj_ids, dateTime.toPython(),
                                         return_pandas=True)
        dia_sources.set_index(["diaObjectId", "diaSourceId"], inplace=True)
        dia_sources["psFlux"].replace(None, np.nan, inplace=True)
        dia_sources["psFluxErr"].replace(None, np.nan, inplace=True)
        dia_sources["totFlux"].replace(None, np.nan, inplace=True)
        dia_sources["totFluxErr"].replace(None, np.nan, inplace=True)

        updated_dia_objects = []

        for obj_id in updated_obj_ids:

            try:
                updated_dia_obj_df = dia_objects.loc[obj_id]
                updated_dia_object = updated_dia_obj_df.to_dict()
                updated_dia_object["diaObjectId"] = obj_id
                dia_object_used.loc[obj_id] = True
            except KeyError:
                updated_dia_object = self._initialize_dia_object(obj_id)

            # Select the dia_sources associated with this DIAObject id and
            # copy the subcatalog for fast slicing.
            obj_dia_sources = dia_sources.loc[obj_id]

            _set_mean_position(updated_dia_object, obj_dia_sources)
            indexer_id = self.indexer.indexPoints(
                [updated_dia_object["ra"]],
                [updated_dia_object["decl"]])[0]
            updated_dia_object['pixelId'] = indexer_id
            updated_dia_object["radecTai"] = dateTime.get(system=DateTime.MJD)

            updated_dia_object["nDiaSources"] = len(obj_dia_sources)
            _set_flux_stats(updated_dia_object,
                            obj_dia_sources,
                            filter_name,
                            filter_id)

            # TODO: DM-15930
            #     Define and improve flagging for DiaObjects
            # Set a flag on this DiaObject if any DiaSource that makes up this
            # object has a flag set for any reason.
            if np.any(obj_dia_sources["flags"] > 0):
                updated_dia_object["flags"] = 1
            else:
                updated_dia_object["flags"] = 0

            updated_dia_objects.append(updated_dia_object)

        # Initialize DataFrame with full set of columns to properly set
        # NaN/Null values.
        updated_dia_objects = pandas.DataFrame(
            data=updated_dia_objects,
            columns=list(ppdb._schema.getColumnMap('DiaObject').keys()))

        ppdb.storeDiaObjects(updated_dia_objects, dateTime.toPython())

        updated_dia_objects.set_index(["diaObjectId"], inplace=True)

        return dia_objects[
            np.logical_not(dia_object_used.loc[:, "used"])].append(updated_dia_objects,
                                                                   sort=True)

    def _initialize_dia_object(self, obj_id):
        """Input default values for non-Nullable DiaObject columns.

        Parameters
        ----------
        obj_id : `int`
            Id to asign to the new object.

        Returns
        -------
        new_dia_object : `dict`
            Dictionary with default values.
        """
        new_dia_object = {"diaObjectId": obj_id,
                          "pmParallaxNdata": 0,
                          "nearbyObj1": 0,
                          "nearbyObj2": 0,
                          "nearbyObj3": 0}
        for f in ["u", "g", "r", "i", "z", "y"]:
            new_dia_object["%sPSFluxNdata" % f] = 0
        return new_dia_object

    @pipeBase.timeMethod
    def score(self, dia_objects, dia_sources, max_dist):
        """Compute a quality score for each dia_source/dia_object pair
        between this catalog of DIAObjects and the input DIASource catalog.

        ``max_dist`` sets maximum separation in arcseconds to consider a
        dia_source a possible match to a dia_object. If the pair is
        beyond this distance no score is computed.

        Parameters
        ----------
        dia_objects : `pandas.DataFrame`
            A contiguous catalog of DIAObjects to score against dia_sources.
        dia_sources : `pandas.DataFrame`
            A contiguous catalog of dia_sources to "score" based on distance
            and (in the future) other metrics.
        max_dist : `lsst.geom.Angle`
            Maximum allowed distance to compute a score for a given DIAObject
            DIASource pair.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results struct with components:

            - ``scores``: array of floats of match quality updated DIAObjects
                (array-like of `float`).
            - ``obj_idxs``: indexes of the matched DIAObjects in the catalog.
                (array-like of `int`)
            - ``obj_ids``: array of floats of match quality updated DIAObjects
                (array-like of `int`).

            Default values for these arrays are
            INF, -1, and -1 respectively for unassociated sources.
        """
        scores = np.full(len(dia_sources), np.inf, dtype=np.float64)
        obj_idxs = np.full(len(dia_sources), -1, dtype=np.int)
        obj_ids = np.full(len(dia_sources), -1, dtype=np.int)

        if len(dia_objects) == 0:
            return pipeBase.Struct(
                scores=scores,
                obj_idxs=obj_idxs,
                obj_ids=obj_ids)

        spatial_tree = self._make_spatial_tree(dia_objects)

        max_dist_rad = max_dist.asRadians()

        vectors = self._radec_to_xyz(dia_sources)

        scores, obj_idxs = spatial_tree.query(
            vectors,
            distance_upper_bound=max_dist_rad)
        matched_src_idxs = np.argwhere(np.isfinite(scores))
        obj_ids[matched_src_idxs] = dia_objects.index[
            obj_idxs[matched_src_idxs]]

        return pipeBase.Struct(
            scores=scores,
            obj_idxs=obj_idxs,
            obj_ids=obj_ids)

    def _make_spatial_tree(self, dia_objects):
        """Create a searchable kd-tree the input dia_object positions.

        Parameters
        ----------
        dia_objects : `pandas.DataFrame`
            A catalog of DIAObjects to create the tree from.

        Returns
        -------
        kd_tree : `scipy.spatical.cKDTree`
            Searchable kd-tree created from the positions of the DIAObjects.
        """
        vectors = self._radec_to_xyz(dia_objects)
        return cKDTree(vectors)

    def _radec_to_xyz(self, catalog):
        """Convert input ra/dec coordinates to spherical unit-vectors.

        Parameters
        ----------
        catalog : `pandas.DataFrame`
            Catalog to produce spherical unit-vector from.

        Returns
        -------
        vectors : `numpy.ndarray`, (N, 3)
            Output unit-vectors
        """
        ras = np.radians(catalog["ra"])
        decs = np.radians(catalog["decl"])
        vectors = np.empty((len(ras), 3))

        sin_dec = np.sin(np.pi / 2 - decs)
        vectors[:, 0] = sin_dec * np.cos(ras)
        vectors[:, 1] = sin_dec * np.sin(ras)
        vectors[:, 2] = np.cos(np.pi / 2 - decs)

        return vectors

    @pipeBase.timeMethod
    def match(self, dia_objects, dia_sources, score_struct):
        """Match DIAsources to DIAObjects given a score and create new
        DIAObject Ids for new unassociated DIASources.

        Parameters
        ----------
        dia_objects : `pandas.DataFrame`
            A SourceCatalog of DIAObjects to associate to DIASources.
        dia_sources : `pandas.DataFrame`
            A contiguous catalog of dia_sources for which the set of scores
            has been computed on with DIAObjectCollection.score.
        score_struct : `lsst.pipe.base.Struct`
            Results struct with components:

            - ``scores``: array of floats of match quality
                updated DIAObjects (array-like of `float`).
            - ``obj_ids``: array of floats of match quality
                updated DIAObjects (array-like of `int`).
            - ``obj_idxs``: indexes of the matched DIAObjects in the catalog.
                (array-like of `int`)

            Default values for these arrays are
            INF, -1 and -1 respectively for unassociated sources.

        Returns
        -------
        result : `lsst.pipeBase.Struct`
            Results struct with components:

            - ``updated_and_new_dia_object_ids`` : ids of new and updated
              dia_objects as the result of association. (`list` of `int`).
            - ``n_updated_dia_objects`` : Number of previously know dia_objects
              with newly associated DIASources. (`int`).
            - ``n_new_dia_objects`` : Number of newly created DIAObjects from
              unassociated DIASources (`int`).
            - ``n_unupdated_dia_objects`` : Number of previous DIAObjects that
              were not associated to a new DIASource (`int`).
        """

        n_previous_dia_objects = len(dia_objects)
        used_dia_object = np.zeros(n_previous_dia_objects, dtype=np.bool)
        used_dia_source = np.zeros(len(dia_sources), dtype=np.bool)
        associated_dia_object_ids = np.zeros(len(dia_sources),
                                             dtype=np.uint64)

        n_updated_dia_objects = 0
        n_new_dia_objects = 0

        # We sort from best match to worst to effectively perform a
        # "handshake" match where both the DIASources and DIAObjects agree
        # their the best match. By sorting this way, scores with NaN (those
        # sources that have no match and will create new DIAObjects) will be
        # placed at the end of the array.
        score_args = score_struct.scores.argsort(axis=None)
        for score_idx in score_args:
            if not np.isfinite(score_struct.scores[score_idx]):
                # Thanks to the sorting the rest of the sources will be
                # NaN for their score. We therefore exit the loop to append
                # sources to a existing DIAObject, leaving these for
                # the loop creating new objects.
                break
            dia_obj_idx = score_struct.obj_idxs[score_idx]
            if used_dia_object[dia_obj_idx]:
                continue
            used_dia_object[dia_obj_idx] = True
            used_dia_source[score_idx] = True
            obj_id = score_struct.obj_ids[score_idx]
            associated_dia_object_ids[score_idx] = obj_id
            n_updated_dia_objects += 1
            dia_sources.loc[score_idx, "diaObjectId"] = obj_id

        # Argwhere returns a array shape (N, 1) so we access the index
        # thusly to retrieve the value rather than the tuple
        for (src_idx,) in np.argwhere(np.logical_not(used_dia_source)):
            src_id = dia_sources.loc[src_idx, "diaSourceId"]
            associated_dia_object_ids[src_idx] = src_id
            dia_sources.loc[src_idx, "diaObjectId"] = src_id
            n_new_dia_objects += 1

        # Return the ids of the DIAObjects in this DIAObjectCollection that
        # were updated or newly created.
        n_unassociated_dia_objects = \
            n_previous_dia_objects - n_updated_dia_objects
        return pipeBase.Struct(
            associated_dia_object_ids=associated_dia_object_ids,
            n_updated_dia_objects=n_updated_dia_objects,
            n_new_dia_objects=n_new_dia_objects,
            n_unassociated_dia_objects=n_unassociated_dia_objects,)

    def _add_association_meta_data(self, match_result):
        """Store summaries of the association step in the task metadata.

        Parameters
        ----------
        match_result : `lsst.pipeBase.Struct`
            Results struct with components:

            - ``updated_and_new_dia_object_ids`` : ids new and updated
              dia_objects in the collection (`list` of `int`).
            - ``n_updated_dia_objects`` : Number of previously know dia_objects
              with newly associated DIASources. (`int`).
            - ``n_new_dia_objects`` : Number of newly created DIAObjects from
              unassociated DIASources (`int`).
            - ``n_unupdated_dia_objects`` : Number of previous DIAObjects that
              were not associated to a new DIASource (`int`).
        """
        self.metadata.add('numUpdatedDiaObjects',
                          match_result.n_updated_dia_objects)
        self.metadata.add('numNewDiaObjects',
                          match_result.n_new_dia_objects)
        self.metadata.add('numUnassociatedDiaObjects',
                          match_result.n_unassociated_dia_objects)
