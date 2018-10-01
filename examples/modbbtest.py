import numpy as np
import ampere
from ampere.data import Spectrum, Photometry
from ampere.emceesearch import EmceeSearch
from ampere.PowerLawAGN import SingleModifiedBlackBody, PowerLawAGN
import corner
import matplotlib as mpl
#mpl.use("Agg")
import matplotlib.pyplot as plt
import os
from astropy import constants as const
from astropy import units as u
from spectres import spectres
import pyphot



if __name__=="__main__":
    """ Set up the inputs for the model """
    """ wavelength grid """
    wavelengths = 10**np.linspace(0.,1.9, 2000)

    """ Choose some model parameters  """
    t_in = 200. #units = K (or 300K)
    scale_in = 1.
    index_in = 0.
    distance = 1.
    lims_in=np.array([ [100., 300.], [0., 2.], [-5, 5], [0.9, 1.1]
                  ])


    """ Initialise the model  """
    model = SingleModifiedBlackBody(wavelengths, lims=lims_in)


    """ Get a test spectrum out of the model """
    #model(t_in, scale_in, index_in, distance)
    #model_flux = model.modelFlux
    model_flux = model(t_in, scale_in, index_in, distance).modelFlux #Do the two lines above with just one line here. 


    """ get synthetic photometry and spectra """
    


    """ (optionally) add some noise """


    """ now set up ampere to try and fit the same stuff """
