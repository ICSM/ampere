

import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.modeling import models
from .models import AnalyticalModel
from scipy.stats import dirichlet


class PowerLawContinuumAbsoluteAbundances(AnalyticalModel):
    """ A power-law continuum with features on top of it

    Computes the flux as a function of wavelength from a continuum source whose emission 
    follows a power law. This continuum is multiplied with (continuum-subtracted) 
    opacities for species which produce features, scaled by abundances (defined in 
    absolute values).
    

    Parameters
    ----------
    wavelengths : float, array-like
        The wavelengths at which to calculate the spectrum.
    redshift : optional Float, default None
        Known redshift of the target. If None, source is treated as having z = 0, i.e. 
        rest frame is equivalent to observed frame.
    lims : float, array-like (4x2)
        The minimum and maximum limits of the free parameters, when assuming flat priors 
        on all parameters.

    Attributes
    ----------
    None

    Methods
    -------
    

    Notes
    -----
    
    
    Examples
    --------
    >>> 
    """

    
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
        self.npars = nSpecies + 1 #there are two parameters for the continuum - only ONE because of absolute abundances
        
    def __call__(self,
                 #multiplicationFactor, # = log(1),
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
            #if np.sum(10**theta[1:]) <= 1. and np.all(theta[1:] < 0.) and -10 < theta[0] < 10.:# and 1.5 > theta[1] > 0.: #basic physical checks first
            if 1.5 > theta[0] > 0. and np.all(theta[1:]) > -20. and np.all(theta[1:]) < 20.: #basic physical checks first
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



class PowerLawContinuumRelativeAbundances(AnalyticalModel):
    """ A power-law continuum with features on top of it

    Computes the flux as a function of wavelength from a continuum source whose emission 
    follows a power law. This continuum is multiplied with (continuum-subtracted) 
    opacities for species which produce features, scaled by relative abundances.
    

    Parameters
    ----------
    wavelengths : float, array-like
        The wavelengths at which to calculate the spectrum.
    redshift : optional Float, default None
        Known redshift of the target. If None, source is treated as having z = 0, i.e. 
        rest frame is equivalent to observed frame.
    lims : float, array-like (4x2)
        The minimum and maximum limits of the free parameters, when assuming flat priors 
        on all parameters.

    Attributes
    ----------
    None

    Methods
    -------
    

    Notes
    -----
    
    
    Examples
    --------
    >>> 
    """
    '''Input: fit parameters (multiplicationFactor, powerLawIndex, dustAbundances), 
              opacities (opacity_array): q_ij where i = wavelength, j = species
              wavelength grid from data (wavelengths)
    Output: model fluxes (modelFlux)'''
    
    def __init__(self, wavelengths, flatprior=True, redshift=None, **kwargs):
        import os
        from scipy import interpolate
        self.flatprior=flatprior
        #print(os.path.dirname(__file__))
        opacityDirectory = os.path.dirname(__file__)+'/../Opacities/'
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
        opacity_array[opacity_array < 0] = 0 #Eliminate negative continuum-subtracted opacities, because this results in absurd fluxes and NaN likelihoods sometimes.
        print(opacity_array)
        self.wavelength = wavelengths
        self.opacity_array = opacity_array
        self.nSpecies = nSpecies
        self.redshift = redshift
        print(self.restwaves)
        print(self.wavelength)
        print(self.redshift)
        self.npars = nSpecies - 1 + 2 #there are two parameters for the continuum, then we need n-1 abundances
        self.npars_ptform = nSpecies + 2 #If the code uses a prior transform, it needs to 
        #self.npars = nSpecies + 1 #there are two parameters for the continuum - only ONE because of absolute abundances
        
    def __call__(self,
                 multiplicationFactor, # = log(1),
                 powerLawIndex, # = log(2),
                 *args, # = (np.ones(self.nSpecies)/self.nSpecies),
                 **kwargs):
        if len(args) == self.nSpecies-1: #emcee-like codes need to end up in this pathway
            relativeAbundances = np.append(10**np.array(args),1.-np.sum(10**np.array(args)))
        elif len(args) == self.nSpecies: #codes which use a prior transform need to end up in this branch
            relativeAbundances = args
        else: #something is wrong, raise an error
            raise ValueError("The number of input abundances does not match. The number of input abundances must have either the same number of entries as, or one fewer than, the number of dust species, but the number of abundances is %s and the number of species is %s."%(len(args), self.nSpecies))
        if np.any(np.array(relativeAbundances) < 0):
            print(relativeAbundances)
            print(self.nSpecies)
            print(len(args), len(relativeAbundances))
            print("Some abundances are negative!")
            raise ValueError()
        #dustAbundances = 10**np.array(args)
        #moved from using relative abundances to absolute abundances. # species/parameters had
        #to be changed throughout the program too
        if self.redshift is not None:
            waves = self.restwaves
        else:
            waves = self.wavelength
        fModel = (np.matmul(self.opacity_array, relativeAbundances)+1)
        fModel = fModel*(waves**powerLawIndex)*(10**multiplicationFactor) #--> not needed when using absolute abundances
        self.modelFlux = fModel

    def lnprior(self, theta, **kwargs):
        if self.flatprior:
            if np.sum(10**theta[2:]) <= 1. and np.all(theta[2:] < 0.) and -10 < theta[0] < 10. and 1.5 > theta[1] > 0.: #basic physical checks first
            #if 1.5 > theta[0] > 0. and np.all(theta[1:]) > -20. and np.all(theta[1:]) < 20.: #basic physical checks first
                try:
                    return dirichlet((1.,)*self.nSpecies).logpdf(10**theta[2:]) #When you have a set of random variates whose sum = 1 by definition, the correct prior is a dirichlet distribution (which is a generalisation of the beta distribution).
                    #The above line returns the log of the probability of getting those n-1 numbers from an n-dimensional dirichlet distribution, so that we can correctly get the last value by taking the sum of the n-1 variates and subtract it from 1.
                except ValueError: #something went wrong with the checks above and the abundances are not within the simplex
                    #print("ValueError!")
                    return -np.inf
                #return 0
            else:
                #print("prior 0")
                #print(theta)
                #print(np.sum(10**theta[2:]))
                #exit()
                return -np.inf
        else:
            raise NotImplementedError()

    def prior_transform(self, u, **kwargs):
        #For a prior transform, u is a set of uniform random variates
        if self.flatprior:
            theta = np.zeros_like(u)
            theta[0] = 20. * u[0] - 10
            theta[1] = 1.5 * u[1]
            #Now we encode the transformation from N uniform random variates U(0,1) to N dirichlet random variates:
            gamma_quantiles = -np.log(u[2:])
            theta[2:] = gamma_quantiles/gamma_quantiles.sum()
        else:
            raise NotImplementedError()
            
        return theta
