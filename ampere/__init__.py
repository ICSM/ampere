"""
Welcome to AMPERE

Ampere is an attempt to produce a fitting environment that can natively handle multiple different kinds of astronomical data with differing information content, even when the model being applied to them might be missing some of the key processes and be unable to actually reproduce all aspects of the observations properly.
This is motivated by the need within out team for a tool to simultaneously fit the SEDs and spectra of dusty objects to constrain, among other things, the dust properties and mineralogy. 
However, the final product will be more general, such that it can be applied to a wide range of astronomical questions.

To achieve this, we include a simple parametric model of correlated noise for each dataset which is marginalised over when fitting the parameters of the models. 
This approach was chosen because, typically, deficiencies in the model result in structured residuals, and structured residuals are equivalent to correlated noise.
This effectively downweights parts of the data which the model doesn't represent well without having to manually identify these regions.

At present, ampere is in the alpha testing phase, but we anticipate a beta release in the near future. If you are interested, please get in touch with us!


## Installation

git clone repository and install with pip (ideally in a new conda environment):

```bash
> git clone git@github.com:ICSM/ampere.git
> cd ampere
> pip install -e .
```

This will install the minimal package, but leave out some extra features. For all features, instead do
```bash
> pip install -e .[all]
```
"""

__version__ = "0.1.2"
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
from . import infer
#from . import emceesearch
#from . import PowerLawAGN #ClassRT
from . import data
from . import utils

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

