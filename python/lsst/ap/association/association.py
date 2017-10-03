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

""" A simple implementation of source association task for ap_verify.
"""

from __future__ import absolute_import, division, print_function

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.geom as afwGeom
from .assoc_db_sqlite import AssociationDBSqliteTask

__all__ = ["AssociationConfig", "AssociationTask"]


class AssociationConfig(pexConfig.Config):
    """ Config class for AssociationTask.
    """
    level1_db = pexConfig.ConfigurableField(
        target=AssociationDBSqliteTask,
        doc='Specify where and how to load and store DIAObjects and '
        'DIASources.',
    )
    maxDistArcSeconds = pexConfig.Field(
        dtype=float,
        doc='Maximum distance in arcseconds to test for a DIASource to be a '
        'match to a DIAObject.',
        default=1.0,
    )


class AssociationTask(pipeBase.Task):
    """!
    Associate DIAOSources into existing DIAObjects.

    This task performs the assocation of detected DIASources in a visit
    with the previous DIAObjects detected over time. It also creates new
    DIAObjects out of DIASources that cannot be associated with previously
    detected DIAObjects.

    Attributes
    ----------
    level1_db : lsst.ap.assoiation.AssoiationDBSqlite
        A wrapper class for handling persitence of DIAObjects and DIASources.
    """

    ConfigClass = AssociationConfig
    _DefaultName = "association"

    def __init__(self, **kwargs):
        """ Initialize the the association task and create the database link.
        """
        pipeBase.Task.__init__(self, **kwargs)
        self.makeSubtask('level1_db')

    @pipeBase.timeMethod
    def run(self, dia_sources, exposure):
        """ Load DIAObjects from the database, associate the sources, and
        persist the results into the L1 database.

        Parameters
        ----------
        dia_sources : lsst.afw.table.SourceCatalog
            DIASources to be associated with existing DIAObjects.
        expMd : lsst.daf.base.PropertySet
            Input exposure metadata containing the bounding box for this
            exposure.
        """
        # Assure we have a Box2D and can use the getCenter method.
        bbox = afwGeom.Box2D(exposure.getBBox())
        wcs = exposure.getWcs()
        ctr_coord = wcs.pixelToSky(bbox.getCenter())
        max_radius = max(
            ctr_coord.angularSeparation(wcs.pixelToSky(pp))
            for pp in bbox.getCorners())

        dia_collection = self.level1_db.load(ctr_coord, max_radius)

        association_result = self.associate_sources(dia_collection,
                                                    dia_sources)

        dia_collection = association_result.dia_collection
        updated_obj_indices = association_result.updated_indices

        self.level1_db.store_updated(dia_collection, updated_obj_indices)

    @pipeBase.timeMethod
    def associate_sources(self, dia_collection, dia_sources):
        """ Associate the input DIASources in to the collection of DIAObjects.

        Parameters
        ----------
        dia_collection : lsst.ap.association.DIAObjectCollection
            Collection of DIAObjects to atempt to associate the input
            DIASources into.
        dia_sources : lsst.afw.table.SourceCatalog
            DIASources to associate into the DIAObjectCollection.
        """
        scores = dia_collection.score(
            dia_sources, self.config.maxDistArcSeconds * afwGeom.arcseconds)
        updated_dia_objects = dia_collection.match(dia_sources, scores)

        return pipeBase.Struct(
            dia_collection=dia_collection,
            updated_indices=updated_dia_objects
        )
