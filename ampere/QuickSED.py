#
# Ampere Stuff
#

import numpy as np
from astropy import constants as const
from astropy import units
try:
    from astropy.modeling import blackbody
except ImportError:
    from astropy.modeling.physical_models import BlackBody
from .models import AnalyticalModel

from scipy.stats import dirichlet

# 
# QuickSED Stuff
# 

import QuickSED

#
# Other Stuff
#

import subprocess

# Silly string to link Ampere and Hyperion together
class QuickSEDModel(Model):
    '''
    Class to link Ampere and Hyperion.

    '''

    def __init__(self,directory = '',
                        prefix = 'new_qsed_model',
                        obsv = False,
        #wavelength grid
                        wave_min = 3e-1,
                        wave_max = 3e3,
                        nwave = int(1001),
        #stellar model
                        dstar = 10.0,
                        tstar = 5770.0,
                        rstar = 1.0,
        #disc model
                        fdisc = [1e-4],
                        tdisc = [300.0],
                        lam0  = [200.0],
                        beta  = [1.5],
        #outputs
                        sed_star  = [],
                        sed_disc  = [],
                        sed_total = [],
                        **kwargs):
        
        #Assign keyword variables to the object
        #Assign all the inputs to __init__ to instance variables with the same name as above
        #this is equivalent to many lines of self.niter = niter etc
        l = locals()
        for key, value in l.items():
            if key == "self": #skip self to avoid possible recursion
                continue
            setattr(self, key, value)

        self.model = QuickSED.QSED(self)
        
        print("QuickSED Model setup complete.")
        #self.dictionary of fixed parameters that are needed for modelling in __call__

    #Only give __call__ the numbers that emcee is going to change.
    def __call__(self,**kwargs):

        #Same as in __init__!
        l = locals()
        for key, value in l.items():
            if key == "self": #skip self to avoid possible recursion
                continue
            setattr(self, key, value)
        
        self.npars += len(self.lam0) + len(self.beta) + len(self.tdisc) + 3 - 1 #factor 3 comes from stellar parameters

        self.model = QuickSED.QSED(self)

        self.model.wave()
        self.model.star()
        self.model.disc()

        self.wave = self.model.sed_wave
        self.flux = self.model.sed_total

    def lnprior(self, theta, **kwargs):
        if self.flatprior:
            if (self.lims[0,0] < theta[0] < self.lims[0,1]) and \
               (self.lims[1,0] < theta[1] < self.lims[1,1]) and \
               (self.lims[2,0] < theta[2] < self.lims[2,1]) and \
               (self.lims[3,0] < theta[3] < self.lims[3,1]) and \
                np.sum(10**theta[4:]) <= 1. and np.all(theta[4:] < 0.): 
                return 0
            else:
                return -np.inf
        else:
            raise NotImplementedError()

    def __str__(self, **kwargs):
        return  str(self.__class__) + '\n'+ '\n'.join(('{} = {}'.format(item, type(self.__dict__[item])) for item in self.__dict__))
    #This will provide a summary of the properties of the model
    def __repr__(self, **kwargs):
        from pprint import pformat
        return "<" + type(self).__name__ + "> " + pformat(vars(self), indent=4, width=1)
