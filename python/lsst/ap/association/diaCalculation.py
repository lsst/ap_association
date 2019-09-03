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

from collections import namedtuple
import numpy as np
import pandas as pd

from lsst.meas.base import (
    BasePlugin,
    CatalogCalculationPluginConfig,
    CatalogCalculationPlugin,
    CatalogCalculationConfig,
    CatalogCalculationTask,
    PluginRegistry,
    PluginMap)
from lsst.meas.base.catalogCalculation import CCContext
import lsst.pipe.base

# Enforce an error for unsafe column/array value setting in pandas.
pd.options.mode.chained_assignment = 'raise'

__all__ = ("DiaObjectCalculationPlugin", "DiaObjectCalculationPluginConfig",
           "DiaObjectCalculationTask", "DiaObjectCalculationConfig")


class DiaObjectCalculationPluginConfig(CatalogCalculationPluginConfig):
    """Default configuration class DIA catalog calculation plugins.
    """
    pass


class DiaObjectCalculationPlugin(CatalogCalculationPlugin):
    """Base class for DIA catalog calculation plugins.

    Task follows CatalogCalculationPlugin with modifications for use in AP.

    Parameters
    ----------
    config : `DiaObjectCalculationPlugin.ConfigClass`
        Plugin configuration.
    name : `str`
        The string the plugin was registered with.
    metadata : `lsst.daf.base.PropertySet`
        Plugin metadata that will be attached to the output catalog
    """

    ConfigClass = DiaObjectCalculationPluginConfig  # documentation inherited

    registry = PluginRegistry(DiaObjectCalculationPluginConfig)
    """List of available plugins (`lsst.meas.base.PluginRegistry`).
    """

    """Add order after flux means and stds are calculated.
    """
    FLUX_MOMMENTS_CALCULATED = 5.0

    plugType = 'single'
    """Does the plugin operate on a single source or the whole catalog (`str`)?
    If the plugin operates on a single source at a time, this should be set to
    ``"single"``; if it expects the whoe catalog, to ``"multi"``.  If the
    plugin is of type ``"multi"``, the `fail` method must be implemented to
    accept the whole catalog. If the plugin is of type ``"single"``, `fail`
    should accept a single source record.
    """

    def __init__(self, config, name, metadata):
        BasePlugin.__init__(self, config, name)

    @classmethod
    def getExecutionOrder(cls):
        r"""Used to set the relative order of plugin execution.
        The values returned by `getExecutionOrder` are compared across all
        plugins, and smaller numbers run first.

        Notes
        -----
        `DiaObjectCalculationPlugin`\s must run with
        `BasePlugin.DEFAULT_CATALOGCALCULATION` or higher.
        All plugins must implement this method with an appropriate run level
        """
        raise NotImplementedError()

    def calculate(self, cat, **kwargs):
        """Perform the calculation specified by this plugin.
        This method can either be used to operate on a single catalog record
        or a whole catalog, populating it with the output defined by this
        plugin.

        Note that results may be added to catalog records as new columns, or
        may result in changes to existing values.

        Parameters
        ----------
        cat : `lsst.afw.table.SourceCatalog` or `lsst.afw.table.SourceRecord`
            May either be a `~lsst.afw.table.SourceCatalog` or a single
            `~lsst.afw.table.SourceRecord`, depending on the plugin type. Will
            be updated in place to contain the results of plugin execution.
        **kwargs
            Any additional keyword arguments that may be passed to the plugin.
        """
        raise NotImplementedError()


class DiaObjectCalculationConfig(CatalogCalculationConfig):
    """Config class for the catalog calculation driver task.
    Specifies which plugins will execute when the `CatalogCalculationTask`
    associated with this configuration is run.
    """

    plugins = DiaObjectCalculationPlugin.registry.makeField(
        multi=True,
        default=["ap_meanPosition",
                 "ap_meanFlux"],
        doc="Plugins to be run and their configuration")


class DiaObjectCalculationTask(CatalogCalculationTask):
    """Run plugins which operate on a catalog of sources.
    This task facilitates running plugins which will operate on a source
    catalog. These plugins may do things such as classifynig an object based
    on source record entries inserted during a measurement task.

    Parameters
    ----------
    plugMetaData : `lsst.daf.base.PropertyList` or `None`
        Will be modified in-place to contain metadata about the plugins being
        run. If `None`, an empty `~lsst.daf.base.PropertyList` will be
        created.
    **kwargs
        Additional arguments passed to the superclass constructor.

    Notes
    -----
    Plugins may either take an entire catalog to work on at a time, or work on
    individual records.
    """
    ConfigClass = DiaObjectCalculationConfig
    _DefaultName = "diaObjectCalculation"

    def __init__(self, plugMetadata=None, **kwargs):
        lsst.pipe.base.Task.__init__(self, **kwargs)
        if plugMetadata is None:
            plugMetadata = lsst.daf.base.PropertyList()
        self.plugMetadata = plugMetadata
        self.plugins = PluginMap()

        self.initializePlugins()

    def initializePlugins(self):
        """Initialize the plugins according to the configuration.
        """

        pluginType = namedtuple('pluginType', 'single multi')
        self.executionDict = {}
        # Read the properties for each plugin. Allocate a dictionary entry for
        # each run level. Verify that the plugins are above the minimum run
        # level for an catalogCalculation plugin. For each run level, the
        # plugins are sorted into either single record, or multi record groups
        # to later be run appropriately
        for executionOrder, name, config, PluginClass in sorted(self.config.plugins.apply()):
            if executionOrder not in self.executionDict:
                self.executionDict[executionOrder] = pluginType(single=[], multi=[])
            if PluginClass.getExecutionOrder() >= BasePlugin.DEFAULT_CATALOGCALCULATION:
                plug = PluginClass(config, name, metadata=self.plugMetadata)
                self.plugins[name] = plug
                if plug.plugType == 'single':
                    self.executionDict[executionOrder].single.append(plug)
                elif plug.plugType == 'multi':
                    self.executionDict[executionOrder].multi.append(plug)
            else:
                errorTuple = (PluginClass, PluginClass.getExecutionOrder(),
                              BasePlugin.DEFAULT_CATALOGCALCULATION)
                raise ValueError("{} has an execution order less than the minimum for an catalogCalculation "
                                 "plugin. Value {} : Minimum {}".format(*errorTuple))

    @lsst.pipe.base.timeMethod
    def run(self, diaObjectCat, diaSourceCat, updatedDiaObjectIds, filterName):
        """The entry point for the DIA catalog calculation task.

        Run method both updates the values in the diaObjectCat and appends
        newly created DiaObjects to the catalog.

        Parameters
        ----------
        diaObjectCat : `pandas.DataFrame`
            DiaObjects to update values of and append new objects to. DataFrame
            should be indexed on "diaObjectId"
        diaSourceCat : `pandas.DataFrame`
            DiaSources associated with the DiaObjects in diaObjectCat. DataFrame
            should be indexed on `["diaObjectId", "filterName", "diaSourceId"]`
        updatedDiaObjectIds : `numpy.ndarray`
            Integer ids of the DiaObjects to update and create.
        filterName : `str`
            String name of the filter being processed.

        Returns
        -------
        returnStruct : `lsst.pipe.base.Struct`
        """
        return self.callCompute(diaObjectCat,
                                diaSourceCat,
                                updatedDiaObjectIds,
                                filterName)

    def callCompute(self,
                    diaObjectCat,
                    diaSourceCat,
                    updatedDiaObjectIds,
                    filterName):
        """Run each of the plugins on the catalog.

        Parameters
        ----------
        catalog : `lsst.afw.table.SourceCatalog`
            The catalog on which the plugins will operate.
        """

        diaObjectUsed = pd.DataFrame(
            False,
            index=diaObjectCat.index,
            columns=["used"])

        updatedDiaObjects = []

        for objId in updatedDiaObjectIds:
            try:
                updatedDiaObjDF = diaObjectCat.loc[objId]
                updatedDiaObject = updatedDiaObjDF.to_dict()
                updatedDiaObject["diaObjectId"] = objId
                diaObjectUsed.loc[objId] = True
            except KeyError:
                updatedDiaObject = self._initialize_dia_object(objId)

            # Sub-select diaSources associated with this diaObject.
            objDiaSources = diaSourceCat.loc[objId]
            # Currently needed as dataFrames loaded from sql do not currently
            # map Null to NaN for custom queries. This can either stay here
            # or move to dax_ppdb or ap_association.
            objDiaSources.replace(to_replace=[None],
                                  value=np.nan)

            # Sub-select on diaSources observed in the current filter.
            filterObjDiaSources = objDiaSources.loc[filterName]

            for runlevel in sorted(self.executionDict):
                for plug in self.executionDict[runlevel].single:
                    with CCContext(plug, updatedDiaObject, self.log):
                        plug.calculate(diaObject=updatedDiaObject,
                                       diaSources=objDiaSources,
                                       filterDiaSources=filterObjDiaSources,
                                       filterName=filterName)

            updatedDiaObjects.append(updatedDiaObject)

        updatedDiaObjects = pd.DataFrame(data=updatedDiaObjects)

        return lsst.pipe.base.Struct(
            diaObjectCat=diaObjectCat[~diaObjectUsed["used"]].append(
                updatedDiaObjects.set_index("diaObjectId"),
                sort=False),
            updatedDiaObjects=updatedDiaObjects)

    def _initialize_dia_object(self, objId):
        """
        """
        new_dia_object = {"diaObjectId": objId,
                          "pmParallaxNdata": 0,
                          "nearbyObj1": 0,
                          "nearbyObj2": 0,
                          "nearbyObj3": 0}
        for f in ["u", "g", "r", "i", "z", "y"]:
            new_dia_object["%sPSFluxNdata" % f] = 0
        return new_dia_object
