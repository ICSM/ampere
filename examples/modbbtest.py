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


    """ Initialise the model  """


    """ Get a test spectrum out of the model """


    """ get synthetic photometry and spectra """


    """ (optionally) add some noise """


    """ now set up ampere to try and fit the same stuff """
