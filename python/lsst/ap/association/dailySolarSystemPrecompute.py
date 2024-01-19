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



import lsst.pex.config
import lsst.pipe.base 

class DailySolarSystemPrecomputeConfig(lsst.pex.config.Config):
    pass

class DailySolarSystemPrecomputeTask(lsst.pipe.base.Task):
    _DefaultName = "dailySolarSystemPrecompute"
    ConfigClass = DailySolarSystemPrecomputeConfig

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        #Make the database connection
        self.mpcorb_db = None

    def run(self, date):
        orbits = self._load_mpcorb()
        state_vectors = self._compute_positions(orbits, date)
        self._write_apdb(state_vectors)

    def _compute_positions(self, orbits, date):
        """Summary

        Parameters
        ----------
        orbits : `pandas.DataFrame`
            MPCORB 6-parameter orbits, plus epoch and ObjID.
        date : `datetime.datetime`
            Time at which to compute states.
        """

        return orbits
    
