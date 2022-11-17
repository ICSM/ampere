#
# Ampere Stuff
#
import ampere
from ampere.models import Model
import numpy as np
from astropy import constants as const
from astropy import units as u
try:
    from astropy.modeling import blackbody
except ImportError:
    from astropy.modeling.physical_models import BlackBody
from .models import AnalyticalModel

from scipy.stats import dirichlet
from scipy.stats import uniform

#
# Other Stuff
#

import subprocess

#constants
h = 6.626e-34
c = 299792458.0 # m/s
k = 1.38e-23
sb = 5.67e-8 #
au     = 1.495978707e11 # m 
pc     = 3.0857e16 # m
lsol   = 3.828e26 # W
rsol   = 6.96342e8 # m
MEarth = 5.97237e24 # kg

um = 1e-6 #for wavelengths in microns

class QuickSEDModel(Model):
    '''This is a very simple model in which the flux is linear in wavelength.

    This model shows you the basics of writing a model for ampere

    '''
    def __init__(self, wavelengths, flatprior = True,
                 lims = np.array([[0.,1e30], #pc
                                  [1.0,1e5], #K
                                  [0.0,1e3], #Rsun
                                  [0.0,1.0], #Ldust/Lstar
                                  [1.0,1e5], #K
                                  [0.0,1e3], #micron 
                                  [-3.0,3.0]]),**kwargs):

        '''The model constructor, which will set everything up

        This method does essential setup actions, primarily things that 
        may change from one fit to another, but will stay constant throughout 
        the fit. This may be things like the grid of wavelengths to calculate
        the model output on, or establishing the dust opacities if involved.
        There are also several important variables it *MUST* define here
        '''
        self.wave = wavelengths
        self.freq = (wavelengths * u.micron).to(u.Hz, equivalencies=u.spectral()).value
        self.npars = 7 #Number of free parameters for the model (__call__()). For some models this can be determined through introspection, but it is still strongly recommended to define this explicitly here. Introspection will only be attempted if self.npars is not defined.
        self.npars_ptform = 7 #Sometimes the number of free parameters is different when using the prior transform instead of the prior. In that case, self.npars_ptform should also be defined.
        #You can do any other set up you need in this method.
        #For example, we could define some cases to set up different priors
        #But that's for a slightly more complex example.
        #Here we'll just use a simple flat prior
        self.lims = lims
        self.flatprior = flatprior
        self.parLabels = ['dstar','tstar','rstar','fdust','tdust','lam0','beta']

        '''Assign any other keywords.'''
        l = locals()
        for key, value in l.items():
            if key == "self": #skip self to avoid possible recursion
                continue
            setattr(self, key, value)

    def __call__(self,dstar,tstar,rstar,fdust,tdust,lam0,beta,**kwargs):
        '''The model itself, using the callable class functionality of python.

        This is an essential method. It should do any steps required to 
        calculate the output fluxes. Once done, it should stop the output fluxes
        in self.modelFlux.
        '''

        l = locals()
        for key, value in l.items():
            if key == "self": #skip self to avoid possible recursion
                continue
            setattr(self, key, value)

        self.lstar = 4.*np.pi*(self.rstar*rsol)**2*sb*self.tstar**4 / lsol

        sc = 1e26*(1./(self.dstar*pc)**2)*(self.wave*um)**2 /c

        photosphere = (rstar*rsol)**2*np.pi*BlackBody().evaluate(self.freq, self.tstar, 1*u.Jy/u.sr)

        photopshere = photosphere * ( self.lstar / np.trapz(photosphere) ) * sc

        self.star = photosphere

        modified = np.where(self.wave >= self.lam0) 
        
        emission = BlackBody().evaluate(self.freq, self.tdust, 1*u.Jy/u.sr)    

        emission[modified] = emission[modified]*(self.lam0/self.wave[modified])**beta

        emission = emission*(self.fdust*np.trapz(self.star)/np.trapz(emission))*sc
        
        self.dust = emission

        self.modelFlux = self.star + self.dust

        return {"spectrum":{"wavelength":self.wave, "flux": self.modelFlux}}

    def lnprior(self, theta, **kwargs):
        if self.flatprior:
            if (self.lims[0,0] < theta[0] < self.lims[0,1]) and \
               (self.lims[1,0] < theta[1] < self.lims[1,1]) and \
               (self.lims[2,0] < theta[2] < self.lims[2,1]) and \
               (self.lims[3,0] < theta[3] < self.lims[3,1]) and \
               (self.lims[4,0] < theta[4] < self.lims[4,1]) and \
               (self.lims[5,0] < theta[5] < self.lims[5,1]) and \
               (self.lims[6,0] < theta[6] < self.lims[6,1]): 
                return 0
            else:
                return -np.inf
        else:
            raise NotImplementedError()

    def prior_transform(self, u, **kwargs):
        '''The prior transform, which takes samples from the Uniform(0,1) distribution to the desired distribution.

        This is only included for completeness and to demonstrate how a prior 
        transform function should look. This example only uses emcee for 
        fitting, which uses the lnprior function instead. Prior transforms are 
        required by nested-sampling codes and similar approaches.
        '''
        if self.flatprior:
            theta = np.zeros_like(u)
            return (self.lims[:,1] - self.lims[:,0]) * u + self.lims[:,0]
        else:
            raise NotImplementedError()

    def __str__(self, **kwargs):
        return  str(self.__class__) + '\n'+ '\n'.join(('{} = {}'.format(item, type(self.__dict__[item])) for item in self.__dict__))
    #This will provide a summary of the properties of the model
    def __repr__(self, **kwargs):
        from pprint import pformat
        return "<" + type(self).__name__ + "> " + pformat(vars(self), indent=4, width=1)


if __name__ == "__main__":
    """ Set up the inputs for the model """
    """ wavelength grid """
    raise NotImplementedError()
