from __future__ import print_function

import numpy as np
from astropy.table import Table

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
        #JPM:
        #Expecting a structure with something like the following:
        #Filter name (or Wavelength/Frequency for spectra)
        #Values (Flux,Magnitude)
        #Uncertainty (Flux,Magnitude,Fractional? -- e.g. IRAS)
        #Units (Either a separate column (phot) or as a header/keyword (spec)?)
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
        self.filterName = filterName #For consistency get lam value from pyPhot
        import pyphot
        filters = filterLibrary.load_filters(filterName)
        self.wavelength = filters.lpivot() #Pivot wavelengths from pyPhot        
        self.value = value
        self.uncertainty = uncertainty #Error bars may be asymmetric!
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
        
        pass

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()

    def synPhot(self, **kwargs):
        pass

    def lnlike(self, synWave, synFlux, **kwargs):
        pass

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