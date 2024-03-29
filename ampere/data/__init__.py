"""
Welcome to AMPERE

This sub-package provides classes for observational data.
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("pgmuvi")
except PackageNotFoundError:
    # package is not installed
    __version__ = "unknown"


# __version__ = "0.1.1"
__copyright__ = """ Copyright (C) 2017  P. Scicluna, F. Kemper, S. Srinivasan
J.P. Marshall, L. Fanciullo, T. Dharmawardena, A. Trejo, S. Hony

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
"""

__license__ = "GNU Public License v3" #

#from .data import Data #Base class should not be exposed, power users should need to go looking for it explicitly
from .photometry import Photometry
from .spectrum import Spectrum


__all__=["Photometry", "Spectrum"]
