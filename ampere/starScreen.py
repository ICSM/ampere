
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
        self.npars = nSpecies - 1 + 3 #there are three parameters for the continuum, as it is a second order polynomial
        
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
        fModel = np.exp((secondOrderConstant*np.log(waves)**2 + firstOrderConstant*np.log(waves)**1 + constant)-fModel)
        #plt.plot(fModel, label = "fModel")
        #plt.legend()
        #plt.show()

        self.modelFlux = fModel

    def lnprior(self, theta, **kwargs):
        if self.flatprior:
            if np.all(theta[3:] > 0.) and -2.0 < theta[0] < 2. and -2.0 < theta[1] < 3.0 and -2.0 < theta[2] < 2.0: #basic physical checks first                
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


    
class PolynomialSource2(AnalyticalModel):
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
        fModel = (secondOrderConstant*waves**2 + firstOrderConstant*waves**1 + constant)*np.exp(-fModel) 
        #plt.plot(fModel, label = "fModel")
        #plt.legend()
        #plt.show()

        self.modelFlux = fModel

    def lnprior(self, theta, **kwargs):
        if self.flatprior:
            if np.all(theta[3:] > 0.) and -10.0 < theta[0] < 10. and -2.0 < theta[1] < 2.0 and -1.0 < theta[2] < 1.0: #basic physical checks first
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


# creating blackbody model with dust opacities or try to use the black

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

                 

    def __call__(self, T1c = 200., F1c = 0.5, T2c = 800., F2c =0.5, T1f = 200., F1f = 0.5, T2f = 800., F2f = 0.5, *args, **kwargs): #default values for temperatures and fractions will most likely never be used. *args will be a list of (8) absolute abundances of the dust species. 
        #relativeAbundances = np.append(10**np.array(args),1.-np.sum(10**np.array(args))) # The 8th abundance is 1 minus the sum of the other 7. We are doing the parameter space exploration in the log space to ensure proper sampling of small fractions.
        dustAbundances = 10**np.array(args) #We are doing the parameter space exploration in the log space to ensure proper sampling of small fractions.

        freq = const.c.value / (self.wavelength*1e-6)
        
        fModel = (np.matmul(self.opacity_array, dustAbundances))
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



