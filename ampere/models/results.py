import numpy as np
from functools import cached_property

import logging


class ModelResults(object):
    """
    A class for holding, processing and passing the results of an Ampere Model

    
    """

    def __init__(self, spectrum=None, profile = None, image = None, cube = None, visibilities = None, closure_phase = None):
        if isinstance(spectrum, dict) and "wavelength" in spectrum and "flux" in spectrum:
            self.spectrum=spectrum
        elif spectrum is not None:
            logging.warn("Invalid model output, spectrum not used.")

        #pass


    #@cachedproperty
    #def spectrum(self):
    #    #This method will automatically compute the spectrum from other results if they exist, but only when asked for
    #    pass

    
#####
#Possible model outputs:
#F(lambda) - spectrum
#F(lambda, IQUV) -Spectropolarimetry
#F(r) - radial profile
#F(r, lambda)
#F(r, IQUV)
#F(r, lambda, IQUV)
#F(s) - profile in some other direction, like a cut or along a slit
#F(s, lambda)
#F(s, IQUV)
#F(s, lambda, IQUV)
#F(xy) - image
#F(xy, IQUV) - imaging polarimetry
#F(xy, lambda) - datacube
#F(xy, lambda, IQUV)
#F(uv) - visibilities
#F(uv, lambda)
#F(uv, IQUV)
#F(uv, lambda, IQUV)
#phi(uv) - (closure) phases
#phi(uv, lambda)
#phi(uv, IQUV)
#phi(uv, lambda, IQUV)

#Do we want to include astrometric observables? Would be nice to be able to fit parallax/proper motion/etc at the same time as fluxes
#Need to think about how best to represent all these and the connections between them
