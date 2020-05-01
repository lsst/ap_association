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

__all__ = ("PackageAlertsConfig", "PackageAlertsTask")

import os
import tempfile

import lsst.alert.packet as alertPack
import lsst.geom as geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.utils import getPackageDir


class PackageAlertsConfig(pexConfig.Config):
    """Config class for AssociationTask.
    """
    schemaFile = pexConfig.Field(
        dtype=str,
        doc="Location to write alerts to.",
        default=os.path.join(getPackageDir("sample_avro_alert"),
                             "schema/2/1/lsst.alert.avsc"),
    )
    cutoutSize = pexConfig.RangeField(
        dtype=int,
        min=0,
        max=1000,
        default=30,
        doc="Dimension of the square image cutouts to package in the alert."
    )
    alertWriteLocation = pexConfig.Field(
        dtype=str,
        doc="Location to write alerts to.",
        default=os.path.join(os.getcwd(), "alerts"),
    )


class PackageAlertsTask(pipeBase.Task):
    """Associate DIAOSources into existing DIAObjects.

    This task performs the association of detected DIASources in a visit
    with the previous DIAObjects detected over time. It also creates new
    DIAObjects out of DIASources that cannot be associated with previously
    detected DIAObjects.
    """
    ConfigClass = PackageAlertsConfig
    _DefaultName = "packageAlerts"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.alertSchema = alertPack.Schema.from_file(self.config.schemaFile)
        self.cutoutBBox = geom.Extent2I(self.config.cutoutSize,
                                        self.config.cutoutSize)
        os.makedirs(self.config.alertWriteLocation, exist_ok=True)

    def run(self,
            diaSourceCat,
            diaObjectCat,
            diaSrcHistory,
            diffIm,
            template,
            ccdExposureIdBits):
        """
        """
        alerts = []
        self._patchDiaSources(diaSourceCat)
        self._patchDiaSources(diaSrcHistory)
        self._patchDiaObjects(diaObjectCat)
        ccdVisitId = diffIm.getInfo().getVisitInfo().getExposureId()
        for srcIndex, diaSource in diaSourceCat.iterrows():
            # Get all diaSources for the associated diaObject.
            diaObject = diaObjectCat.loc[srcIndex[0]]
            if  diaObject["nDiaSources"] > 1:
                objSourceHistory = diaSrcHistory.loc[srcIndex[0]]
            else:
                objSourceHistory = None
            sphPoint = geom.SpherePoint(diaSource["ra"],
                                        diaSource["decl"],
                                        geom.degrees)
            diffImCutout = diffIm.getCutout(sphPoint, self.cutoutBBox)
            templateCutout = None
            alertId = diaSource["diaSourceId"]
            alerts.append(
                self.makeJsonAlert(alertId,
                                   diaSource,
                                   diaObject,
                                   objSourceHistory,
                                   diffImCutout,
                                   templateCutout))
        try:
            with open(os.path.join(self.config.alertWriteLocation,
                                   f"{ccdVisitId}.avro"),
                      "ab+") as f:
                self.alertSchema.store_alerts(f, alerts)
        except:
            print("EXCEPTION THROWN")
            import pdb; pdb.set_trace()

    def _patchDiaSources(self, diaSources):
        """
        """
        diaSources["programId"] = 0
        diaSources["ra_decl_Cov"] = None
        diaSources["x_y_Cov"] = None
        diaSources["ps_Cov"] = None
        diaSources["trail_Cov"] = None
        diaSources["dip_Cov"] = None
        diaSources["i_cov"] = None
        
    def _patchDiaObjects(self, diaObjects):
        """
        """
        diaObjects["ra_decl_Cov"] = None
        diaObjects["pm_parallax_Cov"] = None
        diaObjects["uLcPeriodic"] = None
        diaObjects["gLcPeriodic"] = None
        diaObjects["rLcPeriodic"] = None
        diaObjects["iLcPeriodic"] = None
        diaObjects["zLcPeriodic"] = None
        diaObjects["yLcPeriodic"] = None
        diaObjects["uLcNonPeriodic"] = None
        diaObjects["gLcNonPeriodic"] = None
        diaObjects["rLcNonPeriodic"] = None
        diaObjects["iLcNonPeriodic"] = None
        diaObjects["zLcNonPeriodic"] = None
        diaObjects["yLcNonPeriodic"] = None

    def makeJsonAlert(self,
                      alertId,
                      diaSource,
                      diaObject,
                      objDiaSrcHistory,
                      diffImCutout,
                      templateCutout):
        """
        """
        alert = dict()
        alert['alertId'] = alertId
        alert['diaSource'] = diaSource.to_dict()

        if objDiaSrcHistory is None:
            alert['prvDiaSources'] = objDiaSrcHistory
        else:
            alert['prvDiaSources'] = objDiaSrcHistory.to_dict("records")

        alert['prvDiaForcedSources'] = None
        alert['prvDiaNondetectionLimits'] = None

        alert['diaObject'] = diaObject.to_dict()

        alert['ssObject'] = None

        alert['cutoutDifference'] = {
            'fileName': '',
            'stampData': self.makeCutoutBytes(diffImCutout),
        }
        alert["cutoutTemplate"] = None

        return alert

    def makeCutoutBytes(self, cutout):
        """
        """
        # writeFits seems to want filenames, not file pointers
        temp = tempfile.NamedTemporaryFile()
        cutout.writeFits(temp.name)
        with open(temp.name, 'rb') as f:
            cutout_bytes = f.read()
        return cutout_bytes
