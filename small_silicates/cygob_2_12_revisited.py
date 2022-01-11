import sys
sys.path.insert(1, '/home/zeegers/git_ampere/ampere/')
import numpy as np
import os
import ampere
from ampere.data import Spectrum, Photometry
from ampere.emceesearch import EmceeSearch
from ampere.models import Model
from spectres import spectres
import pyphot
from emcee import moves
from astropy.modeling import models
from astropy import units as u
from astropy.modeling.models import BlackBody
import matplotlib.pyplot as plt
from scipy import interpolate
