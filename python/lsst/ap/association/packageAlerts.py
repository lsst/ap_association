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


class PackageAlertsConfig(pexConfig.Config):
    """Config class for AssociationTask.
    """
    schemaFile = pexConfig.Field(
        dtype=str,
        doc="Location to write alerts to.",
        default="/project/ebellm/sample-avro-alert/schema/2/1/lsst.alert.avsc",
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
        ccdVisitId = diffIm.getInfo().getVisitInfo().getExposureId()
        for srcIndex, diaSource in diaSourceCat.iterrows():
            # Get all diaSources for the associated diaObject.
            objSourceHistory = diaSrcHistory.loc[srcIndex[0]]
            sphPoint = geom.SpherePoint(diaSource["ra"],
                                        diaSource["decl"],
                                        geom.degrees)
            diffImCutout = diffIm.getCutout(sphPoint, self.cutoutBBox)
            templateCutout = None
            alertId = diaSource["diaSourceId"]
            alerts.append(
                self.makeJsonAlert(alertId,
                                   diaSource,
                                   diaObjectCat.loc[srcIndex[0]],
                                   objSourceHistory,
                                   diffImCutout,
                                   templateCutout))

        with open(os.path.join(self.config.alertWriteLocation,
                               f"{ccdVisitId}.avro"),
                  "wb") as f:
            self.alertSchema.store_alerts(f, alerts)

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
