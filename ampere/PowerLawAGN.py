
# Input: wavelengths from obs., dust model (q_ij), parameters A, B, C_j
# j = 0, nmaterial-1
# Wavelengths need include both pivotal wavelengths for photometry and wavelengths for (multiple) spectra
# Regrid the q values onto the wavelengths
# Execution: calculate the model flux; return the result

import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.analytic_functions import blackbody
from models import AnalyticalModel

class SingleModifiedBlackBody(AnalyticalModel):
    def __init__(self, wavelengths, flatprior=True,
                 normWave = 1., sigmaNormWave = 1.,
                 dist = 1., **kwargs):
        self.wavelengths = wavelengths #grid of wavelengths to calculate BB for
        self.freq = const.c / wavelengths #unit conversions will be required...
        self.flatprior = flatprior #whether to assume flat priors
        self.normWave = normWave #wavelength at which opacity is normalised
        self.sigmaNormWave = sigmaNormWave #value to which the opacity is normalised at wavelength normWave
        self.dist = dist

    def __call__(self, t = 1., scale = 1., index = 1., **kwargs):
        
        #return (blackbody.blackbody_nu(const.c.value*1e6/self.wavelengths,t).to(u.Jy / u.sr).value / (dist_lum.value)**2 * kappa230.value * ((wave/230.)**betaf) * massf) #*M_sun.cgs.value
        pass

class PowerLawAGN(AnalyticalModel):
    '''Input: fit parameters (multiplicationFactor, powerLawIndex, relativeAbundances), 
              opacities (opacity_array): q_ij where i = wavelength, j = species
              wavelength grid from data (wavelengths)
    Output: model fluxes (modelFlux)'''

	def __init__(self, wavelengths, flatprior=True, **kwargs):
            import os
            from scipy import interpolate
            self.flatprior=flatprior
            opacityDirectory = './Opacities/'
            opacityFileList = os.listdir(opacityDirectory)
            opacityFileList = np.array(opacityFileList)[['.q' in zio for zio in opacityFileList]] # Only files ending in .q are valid (for now)
            nSpecies = opacityFileList.__len__()
            #wavelengths = np.logspace(np.log10(8), np.log10(35), 100) # For testing purposes
            opacity_array = np.zeros((wavelengths.__len__(), nSpecies))
            for j in range(nSpecies):
                tempData = np.loadtxt(opacityDirectory + opacityFileList[j], comments = '#')
                tempWl = tempData[:, 0]
                tempOpac = tempData[:, 1]
                f = interpolate.interp1d(tempWl, tempOpac, assume_sorted = False)
                opacity_array[:,j] = f(wavelengths)
            self.wavelengths = wavelengths
            self.opacity_array = opacity_array
            self.nSpecies = nSpecies
        
        def __call__(self, multiplicationFactor = 1, powerLawIndex = 2, relativeAbundances = np.ones(self.nSpecies)/self.nSpecies, **kwargs):
            fModel = (np.matmul(self.opacity_array, relativeAbundances)+1)*self.wavelengths**powerLawIndex*multiplicationFactor
            self.modelFlux = fModel

        def lnprior(self, theta, **kwargs):
            if self.flatprior:
                return 0
            else:
                raise NotImplementedError()
