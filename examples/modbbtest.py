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
    
