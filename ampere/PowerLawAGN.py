
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

class DualBlackBodyDust(AnalyticalModel):
    def __init__(self, wavelengths, flatprior=True,
                 normWave = 1., sigmaNormWave = 1.,
                 opacityFileList=[], #opacities, #turning this into an empty list so the code will import correctly
                 lims=np.array([[100,1000],[0,np.inf],
                                [100,1000],[0,np.inf],
                                [100,1000],[0,np.inf],
                                [100,1000],[0,np.inf]]),
                               **kwargs): #arguments are T1c, F1c, T2c, F2c, T1f, F1f, T2f, F2f; with c for continuum and f for features
        self.wavelength = wavelengths #grid of observed wavelengths to calculate BB for
        #self.freq = const.c.value / (wavelengths*1e-6) #unit conversions will be required...
        self.flatprior = flatprior #whether to assume flat priors
        self.normWave = normWave #wavelength at which opacity is normalised
        self.sigmaNormWave = sigmaNormWave #value to which the opacity is normalised at wavelength normWave
        self.lims = lims
                 
        #Define opacities for use in model calculation
        import os
        from scipy import interpolate
        opacityDirectory = os.path.dirname(__file__)+'/Opacities/'
        opacityFileList = np.array(opacityFileList)
        #opacityFileList = np.array(opacityFileList)[['.q' in zio for zio in opacityFileList]] # Only files ending in .q are valid (for now)
        nSpecies = opacityFileList.__len__()
        #print(opacityFileList,nSpecies)
        opacity_array = np.zeros((wavelengths.__len__(), nSpecies))
        for j in range(nSpecies):
            tempData = np.loadtxt(opacityDirectory + opacityFileList[j], comments = '#')
            tempWl = tempData[:, 0]
            tempOpac = tempData[:, 1]
            f = interpolate.interp1d(tempWl, tempOpac, assume_sorted = False)
            opacity_array[:,j] = f(self.wavelength)
        self.opacity_array = opacity_array
        self.nSpecies = nSpecies
        self.npars = nSpecies - 1 + 8 # 8: two temperatures and scaling factors
                 # for the continuum and features each. 

                 

    def __call__(self, T1c = 200., F1c = 0.5, T2c = 800., F2c =0.5, T1f = 200., F1f = 0.5, T2f = 800., F2f = 0.5, *args, **kwargs): #default values for temperatures and fractions will most likely never be used. *args will be a list of (7) relative abundances of the dust species. 
        relativeAbundances = np.append(10**np.array(args),1.-np.sum(10**np.array(args))) # The 8th abundance is 1 minus the sum of the other 7. We are doing the parameter space exploration in the log space to ensure proper sampling of small fractions.

        freq = const.c.value / (self.wavelength*1e-6)
        
        fModel = (np.matmul(self.opacity_array, relativeAbundances))
        fModel = fModel*(F1f*blackbody.blackbody_nu(freq,T1f).to(u.Jy / u.sr).value + F2f*blackbody.blackbody_nu(freq,T2f).to(u.Jy / u.sr).value) + (F1c*blackbody.blackbody_nu(freq,T1c).to(u.Jy / u.sr).value + F2c*blackbody.blackbody_nu(freq,T2c).to(u.Jy / u.sr).value)           
        self.modelFlux = fModel

                 

    def lnprior(self, theta, **kwargs):
        if self.flatprior:
            if (self.lims[0,0] < theta[0] < self.lims[0,1]) and (self.lims[1,0] < theta[1] < self.lims[1,1]) and (self.lims[2,0] < theta[2] < self.lims[2,1]) and (self.lims[3,0] < theta [3] < self.lims[3,1]): 
                return 0
            else:
                return -np.inf
        else:
            raise NotImplementedError()

class SingleModifiedBlackBody(AnalyticalModel):
    def __init__(self, wavelengths, flatprior=True,
                 normWave = 1., sigmaNormWave = 1.,
                 redshift = False, lims=np.array([[0,1e6],
                                                  [-100,100],
                                                  [-10,10],
                                                  [0,np.inf]]),
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
        
    def prior_transform(self, u, **kwargs):
        if self.flatprior:
            theta = np.zeros_like(u)
            return (lims[:,1] - lims[:,0]) * u - lims[:,0]
        else:
            raise NotImplementedError()

class PowerLawAGN(AnalyticalModel):
    '''Input: fit parameters (multiplicationFactor, powerLawIndex, relativeAbundances), 
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
        fModel = fModel*(waves**powerLawIndex)*(10**multiplicationFactor)
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

class OpacitySpectrum(AnalyticalModel):
    '''Input: fit parameters (multiplicationFactor, powerLawIndex, relativeAbundances), 
              opacities (opacity_array): q_ij where i = wavelength, j = species
              wavelength grid from data (wavelengths)
    Output: model fluxes (modelFlux)'''

    def __init__(self, wavelengths, flatprior=True,
                 normWave = 1., sigmaNormWave = 1.,opacityFileList=None,
                 redshift = False, lims=np.array([[0,1e6],[-100,100],[-10,10],[0,np.inf]]),
                 **kwargs):
        #self.wavelength = wavelengths #grid of observed wavelengths to calculate BB for
        self.redshift = redshift
        if redshift:
            from astropy.cosmology import FlatLambdaCDM
            self.cosmo=FlatLambdaCDM(H0=70, Om0=0.3)
        if redshift is not None:
            self.wavelength = wavelengths / (1+redshift)
        else:
            self.wavelength = wavelengths
        #self.freq = const.c.value / (wavelengths*1e-6) #unit conversions will be required...
        self.flatprior = flatprior #whether to assume flat priors
        self.normWave = normWave #wavelength at which opacity is normalised
        self.sigmaNormWave = sigmaNormWave #value to which the opacity is normalised at wavelength normWave
        self.lims = lims #in same order as inputs to __call__
        #print(self.lims)
        #Define opacities for use in model calculation
        import os
        from scipy import interpolate
        opacityDirectory = os.path.dirname(__file__)+'/Opacities/'
        opacityFileList = np.array(opacityFileList)
        #opacityFileList = np.array(opacityFileList)[['.q' in zio for zio in opacityFileList]] # Only files ending in .q are valid (for now)
        nSpecies = opacityFileList.__len__()
        #print(opacityFileList,nSpecies)
        opacity_array = np.zeros((wavelengths.__len__(), nSpecies))
        for j in range(nSpecies):
            tempData = np.loadtxt(opacityDirectory + opacityFileList[j], comments = '#')
            tempWl = tempData[:, 0]
            tempOpac = tempData[:, 1]
            f = interpolate.interp1d(tempWl, tempOpac, assume_sorted = False)
            opacity_array[:,j] = f(self.wavelength)
        self.opacity_array = opacity_array
        self.nSpecies = nSpecies
        self.npars = nSpecies - 1 + 4

    def __call__(self, t = 1., scale = 1., index = 1., dist=1.,
                 *args,
                 **kwargs):
        print(args)
        relativeAbundances = np.append(np.array(args),1.-np.sum(np.array(args))) #np.append(10**np.array(args),1.-np.sum(10**np.array(args)))
        #print(relativeAbundances)
        if self.redshift:
            z = dist
            dist = cosmo.luminosity_distance(z).to(u.m)
            freq = const.c.value / ((self.wavelength/(1.+z))*1e-6)
        else:
            dist = dist*u.pc.to(u.m)
            freq = const.c.value / (self.wavelength*1e-6)
        #Simple bug catches for inputs to opacity spectrum
        if len(relativeAbundances) != self.nSpecies :
            print('Number of weights must be same as number of species')
        #if scale <= 0.0:
        #    print('Scale factor must be positive and non-zero.') #not relevant - scale is a log
        if t <= 0.0:
            print('Temperature must be positive and in Kelvin.')
        bb = blackbody.blackbody_nu(freq,t).to(u.Jy / u.sr).value
        bb = bb / dist**2
        bb = bb * 10**(scale) * self.sigmaNormWave * ((self.wavelength / self.normWave)**index)
        #Subtract the sum of opacities from the blackbody continuum to calculate model spectrum
        fModel = bb * (1.0 - (np.matmul(self.opacity_array, relativeAbundances)))
        self.modelFlux = fModel
        #return (blackbody.blackbody_nu(const.c.value*1e6/self.wavelengths,t).to(u.Jy / u.sr).value / (dist_lum.value)**2 * kappa230.value * ((wave/230.)**betaf) * massf) #*M_sun.cgs.value

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
