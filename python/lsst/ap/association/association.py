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

__all__ = ["AssociationConfig", "AssociationTask"]

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

class AssociationConfig(pexConfig.Config):
    pass


class AssociationTask(pipeBase.Task):
    
    ConfigClass = AssociationConfig
    _DefaultName = "associate_sources"

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)

    def associate_sources(self):
        pass
