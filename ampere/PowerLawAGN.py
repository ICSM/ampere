
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
from scipy.stats import dirichlet

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



class PowerLawAGNRelativeAbundances(AnalyticalModel):
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

class OpacitySpectrum(AnalyticalModel):
    '''This model fits a modfied blackbody multiplied by sum of opacties, consisting of a warm and cold component, not necessarily of the same composition, over two temperature ranges, following Kemper et al. 2002.
We use the following parameters:
       wavelengths : the array of wavelengths considered
       opacityFileList
       acold : (relative) abundances for the cold component
       awarm : (relative) abundances for the warm component
       Tcold : temperature range of the cold component
       Twarm : temperature range of the warm component
       indexp : index p of the density distribution
       indexq : index q of the temperature distribution
       multfact : multiplication factor that the sum of the modified black bodies has to be multiplied with to fit the spectrum
       Output is an array of model fluxes (fluxes), to match wavelengths'''

    def __init__(self, wavelengths, flatprior=True,
                 opacityFileList=None,
                 **kwargs):
        self.wavelength = wavelengths #grid of observed wavelengths to calculate BB for

        #self.freq = const.c.value / (wavelengths*1e-6) #unit conversions will be required...
        self.flatprior = flatprior #whether to assume flat priors

        import os
        from scipy import interpolate
        opacityDirectory = os.path.dirname(__file__)+'/Opacities/'
        opacityFileList = np.array(opacityFileList)
        nSpecies = opacityFileList.__len__()
        opacity_array = np.zeros((wavelengths.__len__(), nSpecies))
        for j in range(nSpecies):
            tempData = np.loadtxt(opacityDirectory + opacityFileList[j], comments = '#')
            tempWl = tempData[:, 0]
            tempOpac = tempData[:, 1]
            f = interpolate.interp1d(tempWl, tempOpac, assume_sorted = False)
            opacity_array[:,j] = f(self.wavelength)
        self.opacity_array = opacity_array
        self.nSpecies = nSpecies
        self.npars = nSpecies + 7

    def __call__(self,
                 *args,
                 **kwargs):
        # print out *args
        print(args)
        # print out **kwargs
        for key, value in kwargs.items():
            print(key + " : " + value)

        #number of dust species in opacityDirectory are all going to be fitted
        #for each we will consider a hot and a cold component. We will allow only 2
        #temperature ranges, the same for all dust species.

        #the module translated from ck_modbb calculates the flux for each component at a single temperature

        #ckmodbb and shbb are now written. Now I just need to figure out how to do the
        #fit over 8 dust species and 2 temperature ranges.

        #let's first do this over fixed temperature ranges, instead of allowing them to be free
        #the temperature ranges are 118-100 K and 60-30 K, following Kemper et al. 2002

        coldcomponent = np.zeros((2, wavelengths.__len__()))
        coldcomponent[0,:] = self.wavelengths
        warmcomponent = np.zeros((2, wavelengths.__len__()))
        warmcomponent[0,:] = self.wavelengths

        counter = 0
        for ac in abundances[0,:]
            onespeciescold = ckmodbb(opacity_array[:,counter], tin = Tcold[0], tout = Tcold[1], index = indexp, n0 = abundances[0,counter])
            coldcomponent[1,:] = coldcomponent[1,:] + onespeciescold
            counter = counter + 1

        counter = 0
        for aw in abundances[0,:]
            onespecieswarm = ckmodbb(opacity_array[:,counter], tin = Twarm[0], tout = Twarm[1], index = indexp, n0 = abundances[0,counter])
            warmcomponent[1,:] = warmcomponent[1,:] + onespecieswarm
            counter = counter + 1

        fModel = np.like(coldcomponent)
        fModel[1,:] = fModel[1,:] + warmcomponent[1,:]

        return fModel[1,:]

    def ckmodbb(self, q, tin, tout, index = 0.5, n0, r0 = 1e15, distance = 910., grainsize = 0.1, steps = 10)
        d = distance * 3.0857e18 #convert distance from pc to cm
        a = grainsize * 1e-4 #convert grainsize from micron to cm
 
        fnu = np.zeros((2, wavelengths.__len__()))
        fnu[1,] = self.wavelengths

        for i in range(steps - 1):
            t = tin - i * (tin-tout)/steps
            power = (t/tin)^(2*index - 6)
            bb = shbb(fnu[:,], t, 0.) 
            fnu[,:] = fnu[,:] + q[,:]*bb[,:]*power*((tin-tout)/steps)
        extra = r0/d
        factor = 4 * math.pi * a * a * r0 * n0 * extra * extra / (3-index)
        fnu[,:] = fnu[,:] * factor
        return, fnu

    def shbb(self, aar, temp, pinda):
        a1 = 3.97296e19
        a2 = 1.43875e4
        mbb = np.copy(aar)                                 
        bbflux = a1/(mbb[:,]^3)/(exp(a2/(mbb[:,]*temp))-1)
        mbb[,:] = bbflux *[mbb:,]^pinda
        return mbb                                 

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
