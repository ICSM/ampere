"""
Welcome to AMPERE

"""

__version__ = "0.1.1"
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
""" #Add a statement like this to each file/module/subpackage we include. Exactly how it should be included is a matter of debate... but including it in variables makes it both part of the source and part of the program itself, meaning it is accessible in python.

__license__ = "GNU Public License v3" #

from . import models
#from . import extinction
from . import basesearch
from . import emceesearch
#from . import PowerLawAGN #ClassRT
from . import data

__all__=["models","extinction","emceesearch","PowerLawAGN","data"] #Nothing to import just yet


#input - data format(s)?
      # n * photometry (filter name, mag/flux, err, limits, JD, grouping) & n * spectroscopy (lam, fnu, sigma_fnu, limits, JD, grouping) & dust properties (optional) & model details (choices, rt code, geometry etc) & classification () 

# Dust library
# Parameter search   - try to use existing methods
# RT                 - interfaces ((MoDust, 2Dust, RadMC-3D, any more as necessary) & (pre-processors))
# goodness of fit    - flexible likelihood function including construction of covariances

#output - file format(s)?
      # summary - parameters, range of parameters
      # samples from search - parameters
                            # spectrum
                            # synthetic photometry
                            # GOF


#Stopping criteria?? 

