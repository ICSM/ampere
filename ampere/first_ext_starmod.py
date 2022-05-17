
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
import matplotlib.pyplot as plt


class PolynomialSource(AnalyticalModel):
    '''Input: fit parameters (multiplicationFactor, powerLawIndex, dustAbundances), 
              opacities (opacity_array): q_ij where i = wavelength, j = species
              wavelength grid from data (wavelengths)
    Output: model fluxes (modelFlux)'''
    
    def __init__(self, wavelengths, flatprior=True, **kwargs):
        import os
        from scipy import interpolate
        self.flatprior=flatprior
        opacityDirectory = os.path.dirname(__file__)+'/Opacities/'
        opacityFileList = os.listdir(opacityDirectory)
        opacityFileList = np.array(opacityFileList)[['sub.q' in zio for zio in opacityFileList]] # Only files ending in sub.q are valid (for now). At the moment there are 6 files that meet this criteria
        nSpecies = opacityFileList.__len__()
        #wavelengths = np.logspace(np.log10(8), np.log10(35), 100) # For testing purposes
        opacity_array = np.zeros((wavelengths.__len__(), nSpecies))
        
        for j in range(nSpecies):
            tempData = np.loadtxt(opacityDirectory + opacityFileList[j], comments = '#')
            print(opacityFileList[j])
            tempWl = tempData[:, 0]
            tempOpac = tempData[:, 1]

            f = interpolate.interp1d(tempWl, tempOpac, assume_sorted = False)
            opacity_array[:,j] = f(wavelengths)#wavelengths)
        self.wavelength = wavelengths
        self.opacity_array = opacity_array
        self.nSpecies = nSpecies
#        print(self.wavelength)
        #self.npars = nSpecies - 1 + 2 #there are two parameters for the continuum, then we need n-1 abundances
        self.npars = nSpecies + 3 #there are three parameters for the continuum, as it is a second order polynomial
        
    def __call__(self,
                 secondOrderConstant, # a in a x^2
                 firstOrderConstant, # b in b x
                 constant, # c
                 *args, # = (np.ones(self.nSpecies)/self.nSpecies),
                 **kwargs):

        dustAbundances = np.array(args) # instead of 10**np.array(args)
        #print('dustAbundances = ',dustAbundances)
        #print('constants a, b, c = ', secondOrderConstant, firstOrderConstant, constant)
        waves = self.wavelength
        #print('opacity_array[0] = ', self.opacity_array[:,0])
        fModel = (np.matmul(self.opacity_array, dustAbundances))
        #plt.plot(self.opacity_array[:,0], label = "0")
        #plt.plot(self.opacity_array[:,1], label = "1")
        #plt.plot(self.opacity_array[:,2], label = "2")
        #plt.plot(self.opacity_array[:,3], label = "3")
        #plt.plot(self.opacity_array[:,4], label = "4")
        #plt.plot(self.opacity_array[:,5], label = "5")
        #plt.plot(fModel, label = "opacities")
        #plt.legend()
        #plt.show()

        #print('fModel.shape = ',fModel.shape)
        fModel = 10**(firstOrderConstant*waves**1 + constant)+10**np.exp(-fModel) 
        #plt.plot(fModel, label = "fModel")
        #plt.legend()
        #plt.show()

        self.modelFlux = fModel

    def lnprior(self, theta, **kwargs):
        if self.flatprior:
            if np.all(theta[3:] > 0.) and -1.0 < theta[0] < 1. and -2.0 < theta[1] < 2.0 and -1.0 < theta[2] < 1.0: #basic physical checks first
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
