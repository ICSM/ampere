from __future__ import print_function

import numpy as np
from models import Model
from extinction import apply as redden #This is the generic extinction package, which is Cython optimised and so should be faster. the dust_extinction package is astropy affiliated and can also be used in many cases but may be slower
#from dust_exctinction import

class BaseExtinctionLaws(Model):
    """

    """
    
    def __init__(self, **kwargs):
        raise NotImplementedError("")
    
    def __str__(self, **kwargs):
        raise NotImplementedError("")
    
    def __repr__(self, **kwargs):
        raise NotImplementedError("")

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def extinctionModel(self, **kwargs): 
        raise NotImplementedError("")
    

class F99Extinction(BaseExtinctionLaws):
    """

    """
    def __init__(self, wavelengths, av = None, rv = None,
                 lims = np.array([[0,100],[0,10]]), #limits on av and rv respectively
                 **kwargs):
        from extinction import fitzpatrick99 as f99
        self.f_ext = f99 #define the extinction function, this is more to provide convenience functions for when people want to evaluate the model themselves
        self.wavelengths = wavelengths #Wavelengths are assumed to be in micron!!
        self.invwaves = 1/wavelengths
        self.npars = 2 #Up to two free parameters for this type of model
        if rv is not None: #R_v will be fixed in the fitting
            self.rv = rv
            self.npars -= 1

        if av is not None: #A_v will be fixed in the fitting
            self.av = av
            self.npars -= 1

        #raise NotImplementedError("")
    
    def __str__(self, **kwargs):
        raise NotImplementedError("")
    
    def __repr__(self, **kwargs):
        raise NotImplementedError("")

    def __call__(self, *args, **kwargs):
        if self.npars == 2:
            #both av and rv are free parameters
            av = args[0]
            rv = args[1]
            self.modelFlux = redden(self.f_ext(1/self.wavelengths, av, rv, unit='invum'), np.ones_like(self.wavelengths)) #f99 gives the extinction in magnitudes for a set of wavelengths for given av, rv. We then redden a set of ones to get the multiplier for the fluxes (i.e. we get A = exp(-tau(lambda)))
        elif self.npars == 1:
            #one of av or rv is fixed
            if self.rv is not None: #Rv is fixed
                av = args[0]
                rv = self.rv
                self.modelFlux = redden(self.f_ext(1/self.wavelengths, av, rv, unit='invum'), np.ones_like(self.wavelengths))
            elif self.av is not None:
                av = self.av
                rv = args[0]
                self.modelFlux = redden(self.f_ext(1/self.wavelengths, av, rv, unit='invum'), np.ones_like(self.wavelengths))
            pass
        elif self.npars == 0:
            #both av and rv are fixed, we just want to return extinction values
            self.modelFlux = redden(self.f_ext(1/self.wavelengths, self.av, self.rv, unit='invum'), np.ones_like(self.wavelengths))
            pass
        #raise NotImplementedError()


    def lnprior(self, theta, **kwargs):
        if self.npars == 2:
            #both av and rv are free parameters
            if self.flatprior:
                if (self.lims[0,0] <= theta[0] <= self.lims[0,1]) and (self.lims[1,0] <= theta[1] <= self.lims[1,1]):
                    return 0
                else:
                    return -np.inf
            pass
        elif self.npars == 1:
            #one of av or rv is fixed
            if self.av is not None:
                if (self.lims[0,0] <= theta[0] <= self.lims[0,1]):
                    return 0
                else:
                    return -np.inf
            elif self.rv is not None:
                if (self.lims[1,0] <= theta[0] <= self.lims[1,1]):
                    return 0
                else:
                    return -np.inf
            pass
        elif self.npars == 0:
            #both av and rv are fixed, we just want to return a constant because the prior is meaningless
            return 0
        pass


    def prior_transform(self, u, **kwargs):
        pass
    
        
class CCMExtinctionLaw(BaseExtinctionLaws):
    """
    """
    
    def __init__(self, **kwargs):
        raise NotImplementedError("")
    
    def extinctionModel(self, extinctionModel, **kwargs):
        return CCMModel(self, wavelength, R_A, A_V, **kwargs)
        
    def a(self, wavelength, type_Wavelength, **kwargs):
        
        inverse_wavelength = 1 / wavelength
        
        if type_Wavelength == "IR":
                
                a = 0.574 * (inverse_wavelength**1.61)
                
        elif type_Wavelength == "NIR" or "Optical":
                
                y = inverse_wavelength - 1.82
                
                a = 1 + (0.17699*y) - (0.5044*(y**2)) - (0.02427*(y**3)) + (0.72085*(y**4)) + (0.01979*(y**5)) - (0.77530*(y**6)) + (0.32999*(y**7))
                
        elif type_Wavelength == "UV":
                
                if 8 >= inverse_wavelength >= 5.9: 
                        
                        P_a = -(0.04473((inverse_wavelength - 5.9)**2)) - (0.009779((inverse_wavelength - 5.9)**3))
                        
                elif inverse_wavelength < 5.9:
                        
                        P_a = 0
                        
                a = 1.752 - (0.316*inverse_wavelength) - (0.104 / ( ((inverse_wavelength - 4.67)**2) + 0.341) ) + P_a
                
        elif type_Wavelength == "FUV":
                
                a = -1.073 - (0.628*(inverse_wavelength - 8)) + (0.137 * ((inverse_wavelength - 8)**2)) - (0.070 * ((inverse_wavelength - 8)**3))

        else: 
                raise NotImplementedError("Incorrect type_Wavelength string") 
        
        
        return a


        
    def b(self, wavelength, type_Wavelength, **kwargs):
            
            inverse_wavelength = 1 / wavelength
            
            if type_Wavelength == "IR":
                    
                    b = 0.527 * (inverse_wavelength**1.61)
                    
            elif type_Wavelength == "NIR" or "Optical":
                    
                    y = inverse_wavelength - 1.82
                    
                    b = (1.41338*y) - (2.28305*(y**2)) - (1.07233*(y**3)) - (5.38434*(y**4)) - (0.62251*(y**5)) + (5.30260*(y**6)) - (2.09002*(y**7))
                    
            elif type_Wavelength == "UV":
                    
                    if 8 >= inverse_wavelength >= 5.9: 
                            
                            P_b = (0.2130((inverse_wavelength - 5.9)**2)) + (0.1207((inverse_wavelength - 5.9)**3))
                                
                    elif inverse_wavelength < 5.9:
                            
                            P_b = 0
                                
                    b = (-3.090) + (1.825*inverse_wavelength) + (1.206 / ( ((inverse_wavelength - 4.62)**2) + 0.263) ) + P_b
                        
            elif type_Wavelength == "FUV":
                    
                    b = 13.670 + (4.257*(inverse_wavelength - 8)) - (0.420 * ((inverse_wavelength - 8)**2)) - (0.374 * ((inverse_wavelength - 8)**3))
                        
            else: 
                raise NotImplementedError("Incorrect type_Wavelength string")
            
            return b

    def CCMModel(self, wavelength, R_A, A_V, **kwargs):
            
        a_lambda = self.a(wavelength, type_Wavelength, **kwargs)
        b_lambda = self.b(wavelength, type_Wavelength, **kwargs)
        
        A_lambda = ( a_lambda + (b_lambda / R_A) ) * A_V
        return A_Lambda

