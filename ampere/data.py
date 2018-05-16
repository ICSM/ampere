from __future__ import print_function

import numpy as np
from astropy.table import Table
import pyphot
from astropy.io import fits

class Data(object):
    """


    """

    def __init__(**kwargs):
        pass

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()

    def lnlike(self, synWave, synFlux, **kwargs):
        pass

    def fromFile(self, filename, format, **kwargs):
        ''' 
        Routine to generate data object from a file containing said data
        '''
        pass

    def fromTable(self, table, format, **kwargs):
        ''' 
        Routine to generate data object from an astropy Table object or a file containing data in a format that can be read in as an astropy Table
        '''
        pass

    
    
#1. Should all the photometry be stored in one object

class Photometry(Data):

    """


    """

    def __init__(self, filterName, value, uncertainty, **kwargs):
        self.filterName = filterName
        
#        self.uncertainty = uncertainty #Error bars may be asymmetric!
        self.fluxUnits = photUnits #May be different over wavelength; mag, Jy  
        self.bandUnits = bandUnits #Should be A or um if taken from pyPhot        
        self.type = 'Photometry'
        
        #identify values in magnitudes, convert to Jy
        np.array(photUnits)
        
        mags = (photUnits == 'mag')
        
        zeropoints = filters[mags].Vega_zero_Jy().magnitude
        value[mags] = zeropoints*10^(-0.4*value[mags])
        uncertainty[mags] = value[mags] - zeropoints*10^(-0.4*(value[mags]+uncertainty[mags]))
        
        #identify values in milliJansky, convert to Jy
        mjy = (photUnits == 'mjy')
        
        value[mjy] = 1000.*value[mjy]
        uncertainty[mjy] = 1000.*uncertainty[mjy]
        self.uncertainty = uncertainty
        self.value = value

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()

    def synPhot(self, **kwargs):
        pass

    def lnlike(self, modWave, modFlux, **kwargs):
        ''' docstring goes here '''
        
        ''' First take the model values (passed in) and compute synthetic photometry '''
        ''' I assume that the filter library etc is already setup '''
        filts, modSed = pyphot.extractPhotometry(modWave,
                                                 modFlux,
                                                 self.filterName,
                                                 Fnu = True,
                                                 absFlux = False
            )

        ''' then update the covariance matrix for the parameters passed in '''
        #skip this for now
        self.covMat
        
        ''' then compute the likelihood for each photometric point in a vectorised statement '''
        a = self.value - modSed

        b = np.log(1./((2*np.pi)**(len(self.value)) * np.linalg.det(self.covMat)
            ) 
        #pass
        probFlux = b + ( -0.5 * ( np.matmul ( a.T, np.matmul(self.covMat, a) ) ) )

    def cov(self, **kwargs):
        pass

class Spectrum(Data):

    """


    """

    def __init__(self, bandpass, value, uncertainty, **kwargs):
        self.frequency = freqSpec #True/False? Wavelength or Frequency spectrum
        self.filterName = filterName #Telescope/Instrument cf photometry
        self.bandpass = bandpass #Will be the grid of wavelength/frequency
        self.value = value #Should always be a flux unless someone is peverse
        self.uncertainty = uncertainty #Ditto
        self.fluxUnits = specFluxUnits #lamFlam/nuFnu, Fnu, Flam, again always be same
        self.bandUnits = specBandUnits #wavelength, frequency (A -> GHz)
        self.type = 'Spectrum'

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()
    
    def cov(self, **kwargs):
        pass

    def lnlike(self, **kwargs):
        pass

class Image(Data):
    #Need separate subclasses for images and radial profiles
    def __init__(**kwargs):
        pass

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()
    
    def cov(self, **kwargs):
        pass

    def lnlike(self, **kwargs):
        pass


class Interferometry(Data):

    def __init__(**kwargs):
        pass

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()
    
    def cov(self, **kwargs):
        pass

    def lnlike(self, **kwargs):
        pass


class Cube(Data):
    #IFU or Channel Maps
    def __init__(**kwargs):
        pass

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()
    
    def cov(self, **kwargs):
        pass

    def lnlike(self, **kwargs):
        pass
