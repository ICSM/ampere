from __future__ import print_function

import numpy as np
from astropy.table import Table
from astropy import constants as const
import pyphot
from astropy.io import fits
# for reading VO table formats into a single table
from astropy.io.votable import parse_single_table

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
    Routine to take input from Data(), convert all fluxes into Jy

    """

    def __init__(self, filterName, value, uncertainty, photUnits, bandUnits, **kwargs):
        self.filterName = filterName
        
        #Create wavelength array for photometry based on pivot wavelengths of
        #filters
        filters = filterLibrary.load_filters(filterName)
        self.wavelength = filters.lpivot.magnitude
        
#        self.uncertainty = uncertainty #Error bars may be asymmetric!
        self.fluxUnits = photUnits #May be different over wavelength; mag, Jy  
        self.bandUnits = bandUnits #Should be A or um if taken from pyPhot        
        self.type = 'Photometry'
        
        #identify values in magnitudes, convert to Jy
        np.array(photUnits)
        
        mags = (photUnits == 'mag')
        
        zeropoints = filters[mags].Vega_zero_Jy.magnitude
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

        return probFlux
        

    def cov(self, **kwargs):
        pass

    def fromFile(self, filename, format, **kwargs):
        ''' 
        Routine to generate photometry data object from a file containing said data
        '''
                # First type votable as take from vizier/sed http://vizier.u-strasbg.fr/vizier/sed/
                # following the astropy docs: http://docs.astropy.org/en/stable/io/votable/index.html
                   
                   # of course in the more elaborate version this will be a case statement with file types
                   #switch
                   #VO
                   # this command reads a VOTable into a single table which is then converted in to astropy table format
                   table = parse_single_table(filename).to_table()
                   ## convert this table with standardised column names?

                   #other formats

                   #endswitch
                   ## pass control to fromTable to define the variables
                   fromTable(table)

    def fromTable(self, table, format, **kwargs):
        ''' 
        Routine to generate data object from an astropy Table object or a file containing data in a format that can be read in as an astropy Table
        '''
                   # extract the variables that we are intrested in from the table
                   # for the moment we use the columns from the VO Table
                   value = table['sed_flux'].data
                   photUnits = table['sed_flux'].unit
                   uncertainty = table['sed_eflux'].data
                   filterName = table['sed_filter'].data

                   # We don't see a need for bandUnits because this info should be part of the filterName in pyphot
                   # The SED VO table has frequency units in GHz.
                   # How does this map on the needs of pyphot?
                   
class Spectrum(Data):

    """


    """

    def __init__(self, wavelength, value, uncertainty, **kwargs):
        self.filterName = filterName #Telescope/Instrument cf photometry
        self.type = 'Spectrum'

        #Wavelength conversion
        self.frequency = freqSpec #True/False? Wavelength or Frequency spectrum
        
        if freqSpec == 'True': #Assume given in GHz if freq, convert to microns
            wavelength = 1.0e6*(const.c.value/(wavelength*1.0e9))
        
        self.bandUnits = bandUnits
        
        if bandUnits == 'A':
            wavelength = 1.0e-4*wavelength
        if bandUnits == 'nm' :
            wavelength = 1.0e-3*wavelength
        if bandUnits == 'um':
            wavelength = wavelength
        
        self.wavelength = wavelength #Will be the grid of wavelength/frequency
        
        self.fluxUnits = specFluxUnits #lamFlam/nuFnu, Fnu, Flam, again always be same

        if fluxUnits == 'Fnu':
            value = CONVERSIONS FROM WHATEVER INTO JANSKYS
        if fluxUnits == 'lamFlam' or fluxUnits = 'nuFnu':
            value = CONVERSIONS FROM WHATEVER INTO JANSKYS
        if fluxUnits == 'Flam':
            value = CONVERSIONS FROM WHATEVER INTO JANSKYS

        self.value = value #Should always be a flux unless someone is peverse
        self.uncertainty = uncertainty #Ditto


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

    def fromFile(self, filename, format, **kwargs):
        ''' 
        Routine to generate photometry data object from a file containing said data
        '''
                # First type votable as take from vizier/sed http://vizier.u-strasbg.fr/vizier/sed/
                # following the astropy docs: http://docs.astropy.org/en/stable/io/votable/index.html
                   
                   # of course in the more elaborate version this will be a case statement with file types
                   #switch
                   #VO
                   # this command reads a CASSIS fits file *yaaar*.fits

                   filename = 'cassis_yaaar_spcfw_14203136t.fits' 
                   hdul = fits.open(filename)
                   hdu = hdul[0]
                   header=hdu.header
                   data = hdu.data
                   table = Table(data,names=[header['COL01DEF'],header['COL02DEF'],header['COL03DEF'],header['COL04DEF'],header['COL05DEF'],header['COL06DEF'],header['COL07DEF'],header['COL08DEF'],header['COL09DEF'],header['COL10DEF'],header['COL11DEF'],header['COL12DEF'],header['COL13DEF'],header['COL14DEF'],header['COL15DEF'],'DUMMY'])
                   table['wavelength'].unit='um'
                   table['flux'].unit='Jy'
                   table['error (RMS+SYS)'].unit='Jy'
                   table['error (RMS)'].unit='Jy'
                   table['error (SYS)'].unit='Jy'
                   table['offset uncertainty (CAL)'].unit='Jy'
                   table['sky'].unit='Jy'
                   table['sky error'].unit='Jy'
                   # note that there is a column called module. this has values 0.0 1.0 2.0 and 3.0 corresponding to SL1, SL2, LL1 and LL2 resp.
                   # we may want consider the independent
                   
                   #other formats

                   #endswitch
                   ## pass control to fromTable to define the variables
                   fromTable(table)

    def fromTable(self, table, format, **kwargs):
        ''' 
        Routine to generate data object from an astropy Table object or a file containing data in a format that can be read in as an astropy Table
        '''
                   # extract the variables that we are intrested in from the table
                   # for the moment we use the columns from the cassis
                   value = table['flux'].data
                   # we should read the paper to see which uncertainties to include
                   photUnits = table['flux'].unit
                   uncertainty = table['error (RMS+SYS)'].data

                   ## here we assign the wavelength unit to the bandUnit. is this correct?
                   bandUnits = table['wavelength'].unit

                   ## what about the modules?

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
