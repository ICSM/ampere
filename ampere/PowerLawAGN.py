
# Input: wavelengths from obs., dust model (q_ij), parameters A, B, C_j
# j = 0, nmaterial-1
# Wavelengths need include both pivotal wavelengths for photometry and wavelengths for (multiple) spectra
# Regrid the q values onto the wavelengths
# Execution: calculate the model flux; return the result

import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.modeling import blackbody
from .models import AnalyticalModel

class SingleModifiedBlackBody(AnalyticalModel):
    def __init__(self, wavelengths, flatprior=True,
                 normWave = 1., sigmaNormWave = 1.,
                 redshift = False, lims=np.array([[0,1e6],[-100,100],[-10,10],[0,np.inf]]),
                 **kwargs):
        self.wavelength = wavelengths #grid of observed wavelengths to calculate BB for
        #self.freq = const.c.value / (wavelengths*1e-6) #unit conversions will be required...
        self.flatprior = flatprior #whether to assume flat priors
        self.normWave = normWave #wavelength at which opacity is normalised
        self.sigmaNormWave = sigmaNormWave #value to which the opacity is normalised at wavelength normWave
        self.redshift = redshift
        self.lims = lims
        print(self.lims)
        if redshift:
            from astropy.cosmology import FlatLambdaCDM
            self.cosmo=FlatLambdaCDM(H0=70, Om0=0.3)

    def __call__(self, t = 1., scale = 1., index = 1., dist=1., **kwargs):
        if self.redshift:
            z = dist
            dist = cosmo.luminosity_distance(z).to(u.m)
            freq = const.c.value / ((self.wavelength/(1.+z))*1e-6)
        else:
            dist = dist*u.pc.to(u.m)
            freq = const.c.value / (self.wavelength*1e-6)
        #Simple bug catches for inputs to blackbody model
        #if scale <= 0.0:
        #    print('Scale factor must be positive and non-zero.') #not relevant - scale is a log
        if t <= 0.0:
            print('Temperature must be positive and in Kelvin.')
        bb = blackbody.blackbody_nu(freq,t).to(u.Jy / u.sr).value
        bb = bb / dist**2
        bb = bb * 10**(scale) * self.sigmaNormWave * ((self.wavelength / self.normWave)**index)
        self.modelFlux = bb
        #return (blackbody.blackbody_nu(const.c.value*1e6/self.wavelengths,t).to(u.Jy / u.sr).value / (dist_lum.value)**2 * kappa230.value * ((wave/230.)**betaf) * massf) #*M_sun.cgs.value

    def lnprior(self, theta, **kwargs):
        if self.flatprior:
            if (self.lims[0,0] < theta[0] < self.lims[0,1]) and (self.lims[1,0] < theta[1] < self.lims[1,1]) and (self.lims[2,0] < theta[2] < self.lims[2,1]) and (self.lims[3,0] < theta [3] < self.lims[3,1]): 
                return 0
            else:
                return -np.inf
        else:
            raise NotImplementedError()

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
        self.wavelength = wavelengths
        self.opacity_array = opacity_array
        self.nSpecies = nSpecies
        
    def __call__(self,
                 multiplicationFactor, # = 1,
                 powerLawIndex, # = 2,
                 relativeAbundances, # = np.ones(self.nSpecies)/self.nSpecies,
                 **kwargs):
        fModel = (np.matmul(self.opacity_array, relativeAbundances)+1)*self.wavelength**powerLawIndex*multiplicationFactor
        self.modelFlux = fModel

    def lnprior(self, theta, **kwargs):
        if self.flatprior:
            return 0
        else:
            raise NotImplementedError()
