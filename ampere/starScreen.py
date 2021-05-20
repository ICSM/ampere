
# Input: wavelengths from obs., dust model (q_ij), parameters A, B, C_j
# j = 0, nmaterial
# Wavelengths need include both pivotal wavelengths for photometry and wavelengths for (multiple) spectra
# Regrid the q values onto the wavelengths
# Execution: calculate the model flux; return the result

import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.modeling import blackbody
from .models import AnalyticalModel

class PolynomialSource(AnalyticalModel):
    '''Input: fit parameters (multiplicationFactor, powerLawIndex, dustAbundances), 
              opacities (opacity_array): q_ij where i = wavelength, j = species
              wavelength grid from data (wavelengths)
    Output: model fluxes (modelFlux)'''
    
    def __init__(self, wavelengths, flatprior=True, redshift=None, **kwargs):
        import os
        from scipy import interpolate
        self.flatprior=flatprior
        #print(os.path.dirname(__file__))
        opacityDirectory = os.path.dirname(__file__)+'/Opacities/'
        opacityFileList = os.listdir(opacityDirectory)
        opacityFileList = np.array(opacityFileList)[['sub.q' in zio for zio in opacityFileList]] # Only files ending in .q are valid (for now)
        nSpecies = opacityFileList.__len__()
        #wavelengths = np.logspace(np.log10(8), np.log10(35), 100) # For testing purposes
        opacity_array = np.zeros((wavelengths.__len__(), nSpecies))
        if redshift is not None:
            self.restwaves = wavelengths / (1+redshift)
        else:
            self.restwaves = wavelengths
        for j in range(nSpecies):
            tempData = np.loadtxt(opacityDirectory + opacityFileList[j], comments = '#')
            print(opacityFileList[j])
            tempWl = tempData[:, 0]
            tempOpac = tempData[:, 1]
            #I think we need to put in the continuum subtraction here as well, in case the data isn't continuum subtracted already. These ones are though, so let's see how it goes.
            #Hopefully Sundar can help us out with this.
            f = interpolate.interp1d(tempWl, tempOpac, assume_sorted = False)
            opacity_array[:,j] = f(self.restwaves)#wavelengths)
        self.wavelength = wavelengths
        self.opacity_array = opacity_array
        self.nSpecies = nSpecies
        self.redshift = redshift
        print(self.restwaves)
        print(self.wavelength)
        print(self.redshift)
        #self.npars = nSpecies - 1 + 2 #there are two parameters for the continuum, then we need n-1 abundances
        self.npars = nSpecies + 2 #there are two parameters for the continuum
        
    def __call__(self,
                 multiplicationFactor, # = log(1),
                 powerLawIndex, # = log(2),
                 *args, # = (np.ones(self.nSpecies)/self.nSpecies),
                 **kwargs):
        #relativeAbundances = np.append(10**np.array(args),1.-np.sum(10**np.array(args)))
        dustAbundances = 10**np.array(args)
        #moved from using relative abundances to absolute abundances. # species/parameters had
        #to be changed throughout the program too
        if self.redshift is not None:
            waves = self.restwaves
        else:
            waves = self.wavelength
        fModel = (np.matmul(self.opacity_array, dustAbundances)+1)
        fModel = fModel*(waves**powerLawIndex) #*(10**multiplicationFactor) --> not needed when using absolute abundances
        self.modelFlux = fModel

    def lnprior(self, theta, **kwargs):
        if self.flatprior:
            if np.sum(10**theta[2:]) <= 1. and np.all(theta[2:] < 0.) and -10 < theta[0] < 10. and 1.5 > theta[1] > 0.: #basic physical checks first
                return 0
            else:
                return -np.inf
        else:
            raise NotImplementedError()

    def prior_transform(self, u, **kwargs):
        if self.flatprior:
            theta = np.zeros_like(u)
            theta[0] = 20. * u[0] - 10
            theta[1] = 1.5 * u[1]
        pass



