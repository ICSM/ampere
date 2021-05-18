from __future__ import print_function

import numpy as np
from astropy.table import Table
from astropy.table import vstack
from astropy import constants as const
import pyphot
from astropy.io import fits
# for reading VO table formats into a single table
from astropy.io.votable import parse_single_table
import astropy.units as u
from astropy.io import ascii
from spectres import spectres
from scipy.stats import norm, halfnorm
from scipy.linalg import inv

class Data(object):
    """A base class to represent data objects and their properties

    This is intended purely as a base class to define the interface. When creating 
    your own Data types you must reimplement all methods except: 
    selectWaves
    which you only need to reimplement if your type of data needs to handle them 
    differently. 


    Parameters
    ----------

    None

    Attributes
    ----------

    None

    Notes
    -----

    Data objects are intended to encapsulate both measurements and their covariances, 
    and provide the means to calculate the likelihood of the encapsulated data given
    some model. If the type of data has some nuisance parameters associated with it
    (e.g. a normalisation term) it must also define the prior for those parameters.
    
    Examples
    --------

    None, since this is the base class
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
        '''Routine to generate data object from a file containing said data
        '''
        pass

    def fromTable(self, table, format, **kwargs):
        '''Routine to generate data object from an astropy Table object or a file containing data in a format that can be read in as an astropy Table
        '''
        pass

    def selectWaves(self, low = 0, up = np.inf, interval = "closed", **kwargs):
        '''Method to generate a mask for the data by wavelength. Uses interval terminology to determine how the limits are treated.
        '''

        if interval == "closed": #Both arguments will be treated with less/greater-than-or-equal-to

            mask = np.logical_and(self.wavelength >= low, self.wavelength <= up)

        elif interval == "left-open": #only the upper limit will be treated with less-than-or-equal-to
            mask = np.logical_and(self.wavelength > low, self.wavelength <= up)
        elif interval == "right-open": #only the lower limit will be treated with less-than-or-equal-to
            mask = np.logical_and(self.wavelength >= low, self.wavelength < up)
        elif interval == "open": #neither limit will be treated with less-than-or-equal-to
            mask = np.logical_and(self.wavelength > low, self.wavelength < up)
            
        #Now we add check to make sure that if masks have previously been defined we don't overwrite them, and only accept values 
        #that pass both masks. Otherwise, we define a mask.
        try:
            self.mask = np.logical_and(mask, self.mask)
        except NameError:
            self.mask = mask
            
        #now we need to create a mask for the covariance matrix
        #The outer product does what we want, producing a matrix which has elements such that cov_mask[i,j] = mask[i] * mask[j]
        #This produces the right answer because boolean multiplication is treated as an AND operation in python
        self.cov_mask = np.outer(self.mask, self.mask)

    def maskNaNs(self, **kwargs):
        '''Method to generate a mask which blocks NaNs in the data.
        '''

        mask = np.logical_and(np.isfinite(self.value), np.isfinite(self.uncertainty))
        
        #Now we add check to make sure that if masks have previously been defined we don't overwrite them, and only accept values 
        #that pass both masks. Otherwise, we define a mask.
        try:
            self.mask = np.logical_and(mask, self.mask)
        except NameError:
            self.mask = mask
        
        #now we need to create a mask for the covariance matrix
        #The outer product does what we want, producing a matrix which has elements such that cov_mask[i,j] = mask[i] * mask[j]
        #This produces the right answer because boolean multiplication is treated as an AND operation in python
        self.cov_mask = np.outer(self.mask, self.mask)
        
    def unmask(self, **kwargs):
        '''A method to overwrite previous masks with True in case something goes wrong
        
        '''
        try:
            mask = np.ones_like(self.mask, dtype=np.bool)
        except NameError:
            mask = np.ones_like(self.value, dtype=np.bool)
            
        self.mask = mask
        self.cov_mask = np.outer(self.mask, self.mask)

    
    
#1. Should all the photometry be stored in one object

class Photometry(Data):
    """A class to represent photometric data objects and their properties

    This is a class to hold photometric data points, and their covariances, 
    along with details about the filters that they were observed in. Given the
    covariances and a spectrum from a model, it computes the photometry that 
    would be observed given the model, and computes the likelihood of the model.
    It also contains routines to read information from a file ready to fit.


    Parameters
    ----------
    filterName : string, array-like
        The names of the filters that this object will hold
    value : float, array-like
        The fluxes or magnitudes corresponding to each filter. 
        Fluxes and magnitudes can be mixed, see `photUnits`.
    uncertainty : float, array-like
        The uncertainty on the fluxes or magnitudes.
    photUnits : {'Jy', 'mJy', 'mag'}, array-like
        The units of the photometry. Should be an array-like of the same length as filterName.
    bandUnits : optional, string, scalar or array-like
        Currently assumes micron ('um') as pyphot converts internally. May be updated in future.
    libName : string
        the path to the pyphot library that holds the relevant filter curves.

    Attributes
    ----------
    None

    Methods
    --------
    lnlike : Calculate the likelihood of the data given the model
    
    Examples
    --------
    create an instance by directly instantiating it

    >>> phot = Photometry(['2MASS_J', '2MASS_K'], [10., 5.], [0.2, 0.1], ['Jy', 'Jy'], libname='path/to/filter/library')

    Or create it from a file containing a table

    >>> phot = Photometry.fromFile('filename.vot')

    """

    def __init__(self, filterName, value, uncertainty, photUnits, bandUnits=None, libName = None, **kwargs):
        self.filterName = filterName

        ''' setup pyphot for this set of photometry '''
        self.pyphotSetup(libName)
        self.filterNamesToPyphot()

        #print(self.filterMask)
        if np.all(self.filterMask):
            self.filterName = filterName
        else:
            self.filterName = filterName[self.filterMask]
        #Create wavelength array for photometry based on pivot wavelengths of
        #filters
        filters = self.filterLibrary.load_filters(self.filterName)#[self.filterMask])
        self.wavelength = np.array([filt.lpivot.magnitude for filt in filters])
        
        #        self.uncertainty = uncertainty #Error bars may be asymmetric!
        try:
            self.fluxUnits = photUnits[self.filterMask] #May be different over wavelength; mag, Jy
        except:
            self.fluxUnits = np.repeat(photUnits,len(filters))
        self.bandUnits = 'um' #Should be um if taken from pyPhot        
        self.type = 'Photometry'
        value = value[self.filterMask]
        uncertainty = uncertainty[self.filterMask]
        
        #identify values in magnitudes, convert to Jy
        np.array(photUnits)
        
        mags = self.fluxUnits == 'mag'
        #print(mags.__repr__())
        #pyphot returns a list of filters, this means this nice boolean masking doesn't work :(
        zeropoints = np.zeros_like(value)
        for i in range(len(mags)):
            if mags[i]:
                zeropoints[i] = filters[i].Vega_zero_Jy.magnitude
                value[i] = zeropoints[i]*10^(-0.4*value[i])
                uncertainty[i] = value[i] - zeropoints*10^(-0.4*(value[i]+uncertainty[i]))
        
        #identify values in milliJansky, convert to Jy
        mjy = (photUnits == 'mjy')
        
        value[mjy] = 1000.*value[mjy]
        uncertainty[mjy] = 1000.*uncertainty[mjy]
        self.uncertainty = uncertainty
        self.value = value

        self.selectWaves()

        self.cov()

        self.npars = 0

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        ''' 
        
        '''
        
        return '\n'.join(self.pformat())
        #raise NotImplementedError()

    def pformat(self, **kwargs):
        '''
        Return the instances as a list of formatted strings
        '''
        
        l=[] #define an empty list for the strings
        #now we need to do a little pre-processing to determine the maximum length of each field
        # not sure exactly how to do this yet, so here are some placeholders
        nFilt=max([len(max(np.array(self.filterName).astype(str), key=len)),6])#6
        nWave=max([len(max(np.array(self.wavelength).astype(str), key=len)),18])
        nVal=max([len(max(np.array(self.value).astype(str), key=len)),8])
        nUnc=max([len(max(np.array(self.uncertainty).astype(str), key=len)),12])
        
        ''' first comes header info '''


        ''' then a table of data '''
        ''' this consists of a few header rows '''
        l.append(
            '{:^{nFilt}} {:^{nWave}} {:^{nVal}} {:^{nUnc}} '.format(
                'Filter','Pivot wavelength','Flux','Uncertainty',
                nFilt = nFilt, nWave = nWave, nVal = nVal, nUnc = nUnc
            )
        )
        l.append('{} {} {} {}'.format('-'*(nFilt),'-'*(nWave),'-'*(nVal),'-'*(nUnc)))

        ''' then a table of values '''
        for i in range(len(self.filterName)):
            l.append(
                '{:<{nFilt}} {:>{nWave}.2e} {:>{nVal}.2e} {:>{nUnc}.2e}'.format(
                    self.filterName[i],self.wavelength[i],self.value[i],self.uncertainty[i],
                nFilt = nFilt, nWave = nWave, nVal = nVal, nUnc = nUnc
                )
            )

        ''' finally we might want a way to output the coviariance matrix, although that might only be useful for plotting '''
        return l
    
    def __repr__(self, **kwargs):
        raise NotImplementedError()

    def pyphotSetup(self, libName = None, **kwargs):
        ''' Given the data, read in the pyphot filter library and make sure we have the right list of filters in memory 

        Parameters
        ----------
        libName : str, optional
            The name of the filter library to use

        Notes
        ------
        Future work: go through multiple libraries from different (user-defined) locations and import htem all
        '''
        
        if libName is None:
            print("No library given, using default pyphot filters")
            libDir = pyphot.__file__.strip('__init__.py')+'libs/'
            libName = libDir + 'synphot_nonhst.hd5' #PhIReSSTARTer.hd5'
        
        self.filterLibrary = pyphot.get_library(fname=libName)

    def filterNamesToPyphot(self, **kwargs):
        """Attempt to convert the set of filter names that the objects was instantiated with so that they match the contents of the pyphot library """
        pyphotFilts = self.filterLibrary.get_library_content()
        #print(pyphotFilts)
        filtsOrig = self.filterName
        l = []
        for filt in self.filterName: ## Problem: filt is not defined at this point.
            l.append(filt in pyphotFilts)
        #try replacing colons and / with _
        #print(l)
        try:
            newTry = [filt.astype(str).replace(':','_').replace('/','_').replace('WISE','WISE_RSR').replace('Spitzer','SPITZER') for filt in self.filterName]
        except AttributeError:
            newTry = [filt.replace(':','_').replace('/','_').replace('WISE','WISE_RSR').replace('Spitzer','SPITZER') for filt in self.filterName]
        #change type to str from byte for filt to make it run <CK>
        newTry = [filt.replace('RSR_RSR','RSR') for filt in newTry]
        print(newTry)
        #newTry = [filt.replace('WISE','WISE_RSR').replace( for filt in newTry]
        for i in range(len(l)):
            l[i] = (newTry[i] in pyphotFilts)
            if l[i]:
                self.filterName[i] = newTry[i]
        self.filterMask = np.array(l)
        print(l,self.filterMask.__repr__())
        #print(np.logical_not(np.array(l)))
        if not np.all(l):
            print("Some filters were not recognised by pyphot. The following filters will be ignored:")
            print(filtsOrig[np.logical_not(np.array(l))])

    def reloadFilters(self, modwaves):
        ''' Use this method to reload the filters after your model has got a defined wavelength grid. 

        This method reloads the filters and interpolates their definitions onto the given 
        wavelength grid. This should be used to make sure the filter definitions are ready 
        to compute synthetic photometry in the likelihood calls.
        '''
        #Create wavelength array for photometry based on pivot wavelengths of
        #filters
        #if modwaves is not None:
        #    pass
        #else:
        filters = self.filterLibrary.load_filters(self.filterName[self.filterMask],
                                                  interp = True,
                                                  lamb = modwaves*pyphot.unit['micron'])
        self.wavelength = np.array([filt.lpivot.magnitude for filt in filters])
        self.filters=filters

    def lnprior(self, theta, **kwargs):
        """Return the prior of any nuisance parameters. 

        Since this implementation has no nuisance parameters, it does nothing."""
        return 0

    def lnlike(self, theta, model, **kwargs):
        r'''Compute the likelihood of the photometry given the model. 

        The likelihood is computed as 
        .. math::
        
            \frac{1}{2} N \ln\left(2\pi\right) - \frac{1}{2}\ln\left(\mathrm{det}C\right) - \frac{1}{2} \left(F_\mathrm{obs} - F_\mathrm{mod})^T C^{-1} \left(F_\mathrm{obs} - F_\mathrm{mod}\right)  

        where N is the number of photometric points, C is the covariance matrix, and F_obs and F_mod are the observed and predicted photmetry, respectively.

        Parameters
        ----------
        theta: None,
            included for compatibility reasons
        model: an instance of Model or a subclass,
            The model for which the likelihood will be computed

        Returns
        -------
        probFlux: float
            The natural logarithm of the likelihood of the data given the model
        '''
        
        ''' First take the model values (passed in) and compute synthetic photometry '''
        ''' I assume that the filter library etc is already setup '''
        filts, modSed = pyphot.extractPhotometry(model.wavelength,
                                                 model.modelFlux,
                                                 self.filters,
                                                 Fnu = True,
                                                 absFlux = False,
                                                 progress=False
            )

        ''' then update the covariance matrix for the parameters passed in '''
        #skip this for now
        self.covMat = self.cov()
        
        ''' then compute the likelihood for each photometric point in a vectorised statement '''
        a = self.value[self.mask] - modSed

        b = -0.5*len(self.value[self.mask]) * np.log(2*np.pi) - (0.5*self.logDetCovMat)
            #np.log(1./((2*np.pi)**(len(self.value)) * np.linalg.det(self.covMat))
            #) 
        probFlux = b + ( -0.5 * ( np.matmul ( a.T, np.matmul(inv(self.covMat[self.cov_mask]), a) ) ) )
        return probFlux
        

    def cov(self, **kwargs):
        '''This routine populates a covariance matrix given some methods to call and parameters for them.

        At present, photometric points are assumed to be uncorrelated, so this 
        routine is provided to ensure consistency of the interface between 
        different data types. The resulting covariance matrix is therefore
        always diagonal.

        Parameters
        ----------
        None

        Returns
        -------
        covMat: float, N times N array-like
            The covariance matrix

        Raises
        ------
        Nothing yet!
        '''

        ''' inititalise covariance matrix as a diagonal matrix '''
        #self.covMat = np.diag(uncertainty**2)#np.diag(np.ones(len(uncertainty)))
        self.varMat = np.outer(self.uncertainty, self.uncertainty) #self.uncertainty[self.mask] * self.uncertainty[self.mask][:,np.newaxis] #make a square array of sigma_i * sigma_j, i.e. variances
        self.covMat = np.diag(np.ones_like(self.uncertainty))#[self.mask]))
        a = self.covMat > 0
        self.covMat[a] = self.covMat[a] * self.varMat[a]# = np.diag(uncertainty**2)
        self.logDetCovMat = np.linalg.slogdet(self.covMat[self.cov_mask])[1]# / np.log(10.)
        print(self.logDetCovMat)
        if self.logDetCovMat == -np.inf: #This needs to be updated to raise an error!
            print("""The determinant of the covariance matrix for this dataset is 0.
            Please check that the uncertainties are positive real numbers """)
            print(self)
            exit()
        return self.covMat

    @classmethod
    def fromFile(cls, filename, format=None, **kwargs):
        '''Create Photometry object from a file containing data

        Parameters
        ----------
        filename: str
            The name of the file to load photometry from
        format: {'VO', }, optional, default 'VO'
            The format of the file.

        Returns
        -------
        Photometry
            A new Photometry instance
        '''
        self=cls.__new__(Photometry)
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
        self.fromTable(table, **kwargs)
        #now return the instance so that it works as intended
        return self

    def fromTable(self, table, format=None, **kwargs):
        '''Populates the photometry instance from an Astropy Table instance.

        This is used to populate an existing Photometry instance from an astropy table, called from fromFile. 

        Parameters
        ----------
        table: astropy.table.Table
            The table containing the photometry data


        Notes
        -----
        This routine is not intended to be used as a standalone routine to 
        instantiate Photometry instances. If you wish to do so, you must first 
        create the instance using p = Photometry.__new__(Photometry), then 
        instantiate it with p.fromTable(table). Since this is basically what 
        Photometry.fromFile() does,we recommend using that unless you have to 
        build the table in memory. We hope to revise this in future.

        '''
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
        #print(photUnits)

                   
        self.__init__(filterName, value, uncertainty, photUnits, **kwargs)
                   
class Spectrum(Data):
    """A class to represent 1D spectra data objects and their properties

    


    Parameters
    ----------
    wavelength : float, array-like
        The wavelengths of the spectrum.
    value : float, array-like
        The fluxes of the spectrum.
    uncertainty : float, array-like
        The uncertainty on the fluxes.
    fluxUnits : {'Jy', 'mJy', 'W/m^2/Hz', 'W/m^2/Angstrom', 'W/cm^2/um', 'erg/cm^2/s/Hz', 'erg/cm^2/s/Angstrom'}, scalar
        The units of the fluxes..
    bandUnits : {'um', 'A', 'nm'}, optional, string, scalar or array-like
        The units of the wavelengths.
    freqSpec : bool, default False
        If the input is wavelength (False, default) of Frequncy (True). Currently only works for frequncies in GHz, otherwise input must be wavelength
    calUnc : float, optional, default 0.01
        Parameter controlling width of prior on the scale factor for absolute calibration uncertainty. Larger values imply larger calibration uncertainty. The default value corresponds to roughly 2.5% uncertainty. 
    covarianceTruncation : float, optional, default 1e-3
        Threshold below which the covariance matrix is truncated to zero.

    Attributes
    ----------
    None

    Methods
    -------
    lnlike : Calculate the likelihood of the data given the model

    Notes
    -----
    
    
    Examples
    --------
    >>> spec = Spectrum(wavelength_array, flux_array, uncertainty_array, "um", "Jy")

    or alternatively
    >>> spec = Spectrum.fromFile('filename', 'SPITZER-YAAAR')
    """

    """


    """

    def __init__(self, wavelength, value, uncertainty, bandUnits, fluxUnits, freqSpec = False, calUnc = None, covarianceTruncation = 1e-3, **kwargs):
        #self.filterName = filterName #Telescope/Instrument cf photometry
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
        
        self.fluxUnits = fluxUnits#specFluxUnits
        
        if fluxUnits == 'Jy':
            value = value
            uncertainty = uncertainty
        elif fluxUnits == 'mJy':
            value = 1.0E+3*value
            uncertainty = 1.0E+3*uncertainty
        elif fluxUnits == 'W/m^2/Hz':
            value = 1.0E+26*value
            uncertainty = 1.0E+26*uncertainty
        elif fluxUnits == 'W/m^2/Angstrom':
            value = 2.99792458E-12*value/(1.0.E+4*wavelength)^2
            uncertainty = 2.99792458E-12*uncertainty/(1.0E+4*wavelength)^2
        elif fluxUnits == 'W/cm^2/um':
            value = 2.99792458E-16*value/wavelength^2
            uncertainty = 2.99792458E-16*uncertainty/wavelength^2
        elif fluxUnits == 'erg/cm^2/s/Hz':
            value = 1.0E+23*value
            uncertainty = 1.0E+23*uncertainty
        elif fluxUnits == 'erg/cm^2/s/Angstrom':
            value = 3.33564095E+04*value*(wavelength*1e4)^2
            uncertainty = 3.33564095E+04*uncertainty*(wavelength*1e4)^2
        else:
            raise NotImplementedError()

        self.value = value #Should always be a flux unless someone is peverse
        self.uncertainty = uncertainty #Ditto

        ''' inititalise covariance matrix as a diagonal matrix '''
        self.varMat = uncertainty * uncertainty[:,np.newaxis] #make a square array of sigma_i * sigma_j, i.e. variances
        self.covMat = np.diag(np.ones_like(uncertainty))
        a = self.covMat > 0
        self.covMat[a] = self.covMat[a] * self.varMat[a]# = np.diag(uncertainty**2)
        self.logDetCovMat = np.linalg.slogdet(self.covMat)[1]# / np.log(10.)
        print(self.logDetCovMat)

        ''' Assume default of 10% calibration uncertainty unless otherwise specified by the user '''
        if calUnc is None:
            self.calUnc = 0.010 
        else:
            self.calUnc = calUnc

        self.npars = 3 #this will have to change in future to 1 + number of parameters required by GPs for covMat

        self.covarianceTruncation = covarianceTruncation

        self.selectWaves()

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        ''' 
        
        '''
        
        return '\n'.join(self.pformat())
        #raise NotImplementedError()

    def pformat(self, **kwargs):
        '''
        Return the instances as a list of formatted strings
        '''
        
        l=[] #define an empty list for the strings
        #now we need to do a little pre-processing to determine the maximum length of each field
        # not sure exactly how to do this yet, so here are some placeholders
        #nFilt=max([len(max(np.array(self.filterName).astype(str), key=len)),6])#6
        nWave=max([len(max(np.array(self.wavelength).astype(str), key=len)),18])
        nVal=max([len(max(np.array(self.value).astype(str), key=len)),8])
        nUnc=max([len(max(np.array(self.uncertainty).astype(str), key=len)),12])
        
        ''' first comes header info '''


        ''' then a table of data '''
        ''' this consists of a few header rows '''
        l.append(
            ' {:^{nWave}} {:^{nVal}} {:^{nUnc}} '.format(
                'Wavelength','Flux','Uncertainty',
                nWave = nWave, nVal = nVal, nUnc = nUnc
            )
        )
        l.append('{} {} {}'.format('-'*(nWave),'-'*(nVal),'-'*(nUnc)))

        ''' then a table of values '''
        for i in range(len(self.wavelength)):
            l.append(
                '{:>{nWave}.2e} {:>{nVal}.2e} {:>{nUnc}.2e}'.format(
                    self.wavelength[i],self.value[i],self.uncertainty[i],
                 nWave = nWave, nVal = nVal, nUnc = nUnc
                )
            )

        l.append("Number of wavelength points: "+str(len(self.wavelength)))
        
        ''' finally we might want a way to output the coviariance matrix, although that might only be useful for plotting '''
        return l
        #raise NotImplementedError()

    def __repr__(self, **kwargs):
        """ hacky short-term solution so that we get some sort of print out so people can inspect spectra """
        return self.__str__()
        #raise NotImplementedError()
    
    def cov(self, theta, **kwargs):
        '''This routine populates a covariance matrix given some methods to call and parameters for them.

        It populates the covariance matrix assuming the noise consists of two 
        components: 1) an uncorrelated component; and 2) a single squared-exponential
        with two parameters.

        Parameters
        ----------
        theta : float, length-2 sequence
            Theta contains the parameters for the covariance kernel. The first 
            value is a scale factor, and the second is the scale length (in wavelength)
            of the covariance kernel (the standard deviation of the gaussian).
        '''
        ''' Gradually build this up - for the moment we will assume a single squared-exponential kernel + an uncorrelated component '''

        ''' create a grid of positions/distances on which to build the matrix '''
        #we might be able to offload this step to an instance variable, as it shouldn't change between iterations...
        i, j = np.mgrid[:len(self.wavelength),:len(self.wavelength)]# np.mgrid[:len(self.wavelength[self.mask]),:len(self.wavelength[self.mask])]
        d = i - j

        ''' hardcode a squared-exponential kernel for the moment '''
        m = np.exp(-d**2. / (2.* theta[1]**2.))
        a = m < self.covarianceTruncation
        m[np.logical_not(a)] = 0. #overwrite small values with 0 to speed up some of the algebra

        #self.varMat = #uncertainty[self.mask] * uncertainty[self.mask][:,np.newaxis] #make a square array of sigma_i * sigma_j, i.e. variances
        self.varMat = np.outer(self.uncertainty, self.uncertainty)
        #covMat = (1-theta[0])*np.diag(np.ones_like(self.uncertainty[self.mask])) + theta[0]*m
        covMat = (1-theta[0])*np.diag(np.ones_like(self.uncertainty)) + theta[0]*m
        self.covMat = covMat * self.varMat
        
        self.logDetCovMat = np.linalg.slogdet(self.covMat[self.cov_mask])[1]# / np.log(10.)
        #return self.covMat

    def lnprior(self, theta, **kwargs):
        """Return the prior of any nuisance parameters. 

        This implementation has 3 nuisance parameters. These are 1) a scale factor, 
        which multiplies the observed spectrum to include the uncertainty on 
        absolute calibration, 2) the strength of the correlated component of the
        covariances and 3) the scale length (standard deviation) of the Gaussian
        kernel for the correlated component of the covariances. All three parameters
        are passed in with a single sequence.
        """
        try:
            scaleFac = theta[0]
        except IndexError: #Only possible if theta is scalar or can't be indexed
            scaleFac = theta
        if scaleFac > 0 and 0. < theta[1] < 1. and theta[2] > 0.:
            #print(scalefac)
            return norm.logpdf(np.log10(scaleFac), loc=0., scale = self.calUnc) + halfnorm.logpdf(theta[2], 0., 1.)
        return -np.inf

    def lnlike(self, theta, model, **kwargs):
        '''Compute the likelihood of the photometry given the model. 

        The likelihood is computed as:
        .. math:: \frac{1}{2} N \ln\left(2\pi\right) - \frac{1}{2}\ln\left(\mathrm{det}C\right) - \frac{1}{2} \left(\alpha F_\mathrm{obs} - F_\mathrm{mod})^T C^{-1} \left(\alpha F_\mathrm{obs} - F_\mathrm{mod}\right)
        where N is the number of points in the spectrum, C is the covariance matrix, and F_obs and F_mod are the observed and predicted photmetry, respectively.

        Parameters
        ----------
        theta: empty, included for compatibility reasons
        model: an instance of Model or a subclass

        Returns
        -------
        probFlux: float
            The natural logarithm of the likelihood of the data given the model
        '''
        
        ''' First take the model values (passed in) and compute synthetic Spectrum '''
        
        try:
            scaleFac = theta[0]
        except IndexError: #Only possible if theta is scalar or can't be indexed
            scaleFac = theta
            
        #print(self)
        #wavelength = self.wavelength
        #modSpec = model.modelFlux #
        modSpec = spectres(self.wavelength[self.mask], model.wavelength, model.modelFlux) #For some reason spectres isn't cooperating :/ actually, looks like it was just a stupid mistake

        ''' then update the covariance matrix for the parameters passed in '''
        #skip this for now
        #self.covMat =
        self.cov(theta[1:])
        #import matplotlib.pyplot as plt
        #plt.imshow(self.covMat)
        #plt.show()
        
        ''' then compute the likelihood for each photometric point in a vectorised statement '''
        a = scaleFac*self.value[self.mask] - modSpec
        #if not np.all(np.isfinite(a)):
        #    print(a)
        #    print(modSpec)

        #make this a try: except OverflowError to protect against large spectra (which will be most astronomical ones...)?
        b = 0#np.log10(1./((np.float128(2.)*np.pi)**(len(self.value)) * np.linalg.det(self.covMat))
            #)

        b = -0.5*len(self.value[self.mask]) * np.log(2*np.pi) - (0.5*self.logDetCovMat) #less computationally intensive version of above
        #pass
        probFlux = b + ( -0.5 * ( np.matmul ( a.T, np.matmul(inv(self.covMat[self.cov_mask]), a) ) ) )
        #print(((np.float128(2.)*np.pi)**(len(self.value))), np.linalg.det(self.covMat))
        #print(((np.float128(2.)*np.pi)**(len(self.value)) * np.linalg.det(self.covMat)))
        #print(b, probFlux)
        #exit()
        return probFlux 

    @classmethod
    def fromFile(cls, filename, format, filetype=None, **kwargs):
        '''Create Spectrum object from a file containing data

        Parameters
        ----------
        filename: str
            The name of the file to load photometry from
        format: {'SPITZER-YAAAR', 'SWS-AAR', 'IRAS-LRS', 'User-Defined'}
            The format of the file.
        filetype: {'fits', 'text'}, required if format=="User-Defined" 
            The type of file for user-defined data formats, required to ensure
            that the correct routines are used to open the file.

        Returns
        -------
        Spectrum
            A new Spectrum instance
        '''
        # this command reads a CASSIS fits file *yaaar*.fits
        #self=cls.__new__(Spectrum)

        # we use the format argument to deal with different kind of inputs
        # note that the might want to normalise the names of the columns coming out
        # out of the different files here. 
        if format == 'SPITZER-YAAAR':
            #filename = 'Testdata/cassis_yaaar_spcfw_14203136t.fits' 
            hdul = fits.open(filename)
            hdu = hdul[0]
            header=hdu.header
            data = hdu.data
            table = Table(data,names=[header['COL01DEF'],header['COL02DEF'],header['COL03DEF'],header['COL04DEF'],header['COL05DEF'],header['COL06DEF'],header['COL07DEF'],header['COL08DEF'],header['COL09DEF'],header['COL10DEF'],header['COL11DEF'],header['COL12DEF'],header['COL13DEF'],header['COL14DEF'],header['COL15DEF'],'DUMMY'])
            #table.pprint(max_lines = -1)
            #exit()
            table['wavelength'].unit='um'
            table['flux'].unit='Jy'
            table['error (RMS+SYS)'].unit='Jy'
            table['error (RMS)'].unit='Jy'
            table['error (SYS)'].unit='Jy'
            table['offset uncertainty (CAL)'].unit='Jy'
            table['sky'].unit='Jy'
            table['sky error'].unit='Jy'
            mask = np.logical_and.reduce([np.isfinite(c) for c in table.columns.values()]) #require the elements to be non-NaNs
            table = table[mask]
            chunks = np.zeros_like(table['module'].data)
            sl = np.logical_or(table['module'] == 0.0, table['module'] == 1.0) #SL
            ll = np.logical_or(table['module'] == 2.0, table['module'] == 3.0) #LL
            chunks[sl] = 1.
            chunks[ll] = 2.
            table['chunk'] = chunks

            #normalise the column names (wavelength,flux,uncertainty)
            table.rename_column('error (RMS+SYS)','uncertainty')
            
            #TEMPORARY HACK - just rearrange data into ascending order of wavelength
            #In future we need to extract SL and LL into separate objects, then re-arrange the orders and stitch them together so we end up with one spectrum
            ''' If I'm interpreting the CASSIS data correctly, Module(SL) = 0, 1; Module(LL) = 2, 3 '''
            #a = table['module'] > 1.
            #tablell = table[a]#.sort(keys='wavelength')
            #tablell.sort(keys='wavelength')
            #tablell.pprint()
            #table = tablell
            table.sort(keys='wavelength')
            
        # ISO SWS AAR fits files
        if format == 'SWS-AAR':
            #filename = 'Testdata/sws_ngc6790.fits'
            hdul = fits.open(filename)
            # data is actually in extension 1
            hdu = hdul[1]
            header=hdu.header
            data = hdu.data
            table = Table(data,names=['wavelength','flux','uncertainty','integration count','detector number','time key','universal time key','rpid','spare','line number','scan direction','counts','status','flag'])
            table['wavelength'].unit='um'
            table['flux'].unit='Jy'
            table['uncertainty'].unit='Jy'

        # as found on http://isc83.astro.unc.edu/iraslrs/getlrs_test.html
        # corrected raw text
        # very basic files with the following columns
        # 1 line header with name
        # two columns with wavelength (um) and flux (lambda*Flambda in Watts/meter^2)
        # the two orders of the spectrum are separated by a line with
        # -------------------
        if format == 'IRAS-LRS':
            #filename='Testdata/IRAS06266-0905_lrs.dat'
            lookup='-------------------'
            with open(filename,'r') as f:
                for (num, line) in enumerate(f):
                    if lookup in line:
                        splitnum=num
            table1=ascii.read(filename,data_start=1,data_end=splitnum,header_start=None)
            table1['module']=1
            table2=ascii.read(filename,data_start=splitnum+1,header_start=None)
            table2['module']=2

            table=vstack([table1,table2])

            table.rename_column('col1','wavelength')
            table.rename_column('col2','flux')
            table['wavelength'].unit='um'
            table['flux'].unit='W/m^2'
            # We probably want to convert here to Jy
            # Also we have to define uncertainties

        #create a user-defined class where the relevant keywords and column
        #names are provided to the read-in 
        #user will need to provide:
        #file type (fits/text)
        #column names for flux,unc,wavelength
        #keywords for units of flux+unc, wavelength if in header:
        #waveCol = wavelength column name
        #fluxCol = flux column name
        #uncsCol = uncertainty column name
        #waveUnit = wavelengths unit
        #fluxUnit = flux/uncertainty unit
        #filterName = instrument 
        #
        if format == 'User-Defined':
            filename=filename
            
            if filetype == 'fits':
                hdul = fits.open(filename)
                hdu = hdul[1]
                header=hdu.header
                data = hdu.data
                table = Table(data,names=columns_names_list)
                table.rename_column(keywords['waveCol'],'wavelength')
                table.rename_column(keywords['fluxCol'],'flux')
                table.rename_column(keywords['uncsCol'],'uncertainty')
                table['wavelength'].unit=header[keywords['waveUnit']]
                table['flux'].unit=header[keywords['fluxUnit']]
                table['uncertainty'].unit=header[keywords['fluxUnit']]
                # We probably want to convert here to Jy, um
                # Need to add optional columns for separate orders/modules
                     
            if filetype == 'text':
                
                table=ascii.read(filename,data_start=1)
    
                table.rename_column(keywords['waveCol'],'wavelength')
                table.rename_column(keywords['fluxCol'],'flux')

                table['wavelength'].unit=keywords['waveUnit']
                table['flux'].unit=keywords['fluxUnit']
                # Need to add optional columns for separate orders/modules

                # We probably want to convert here to Jy, um
                if keywords['waveUnit'] != 'um':
                    pass #CONVERT TO MICRONS HERE
                if keywords['fluxUnit'] != 'Jy':
                    pass #CONVERT TO JANSKYS HERE
                # Also we have to define uncertainties
                if keywords['uncsCol'] != '':
                    table.rename_column(keywords['uncsCol'],'uncertainty')
                    table['uncertainty'].unit=keywords['fluxUnit']  
                else:
                    table.add_column('uncertainty')
                    table['uncertainty'].data[:] = 0.0
        ## pass control to fromTable to define the variables
        specList = cls.fromTable(table)
        return specList#elf

    @classmethod
    def fromTable(cls, table, format=None, **kwargs):
        ''' 
            Routine to generate data object from an astropy Table object or a file containing data in a format that can be read in as an astropy Table
        '''
        # extract the variables that we are intrested in from the table
        # for the moment we use the columns from the cassis
        value = table['flux'].data
        # we should read the paper to see which uncertainties to include
        photUnits = table['flux'].unit
        uncertainty = table['uncertainty'].data
        
        ## here we assign the wavelength unit to the bandUnit. is this correct?
        ## Yes, that is correct.
        bandUnits = table['wavelength'].unit
        specList = []
        ## what about the modules?
        # the loop below breaks the spectrum up into chunks based on a column in the table. We can decide whether to group SL and LL, or split all the orders separately
        for chunk in np.unique(table['chunk'].data):
            selection = table['chunk'] == chunk
            self=cls.__new__(Spectrum)
            self.__init__(table['wavelength'][selection].data, value[selection], uncertainty[selection], bandUnits, photUnits) #Also pass in flux units
            specList.append(self)
        return specList #self

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

        ''' docstring goes here '''

        raise NotImplementedError()
        
        ''' First take the model values (passed in) and compute synthetic Image '''
        modImage = 1 #Needs to be added here. Take synthetic photometry on the image plane and convolve with the PSF. 

        ''' then update the covariance matrix for the parameters passed in '''
        #skip this for now
        self.covMat = self.cov()
        
        ''' then compute the likelihood for each photometric point in a vectorised statement '''
        a = self.value - modImage

        b = np.log(1./((2*np.pi)**(len(self.value)) * np.linalg.det(self.covMat))
            ) 
        #pass
        probFlux = b + ( -0.5 * ( np.matmul ( a.T, np.matmul(inv(self.covMat), a) ) ) )

        return probFlux

    @classmethod
    def fromFile(cls, filename, format, **kwargs):
        ''' 
        Routine to generate image data object from a file containing said data
        '''
        # Assume that all image files are fits files in the first instance
        # Further assume that we may be given either a list of images, or a
        # multiple extension fits file.
        self=cls.__new__(Image)
        
        for fname in filename:
            #Look for 'image' extension name to identify image extensions
            #Use hdu[1+].data as the image in the absence of other knowledge
            hdul = fits.open(fname)
            images = np.empty()
            wavelength = np.empty()
            filterName = np.empty()
            imageScale = np.empty() 
            imageSize  = np.empty()
            try:
                for i in range(0,len(hdul)):
                    hdutype = hdul[i]._summary()[0]
                    if hdutype == 'image':
                        header=hdul[i].header
                        #Take keywords here i.e. wavelength, pixel scale, etc.
                        wavelength = np.append(hdul[i].header['WAVELENG'])
                        fluxUnits  = np.append(hdul[i].header['BUNIT'])
                        filterName = np.append(hdul[i].header['INSTRMNT'])
                        imageScale = np.append(hdul[i].header['PIXSCALE'])
                        imageSize  = np.append([hdul[i].header['NAXIS1'],
                                           hdul[i].header['NAXIS2']])
                        images = np.append(images,hdul[i].data)
                        #Get PSF from image extension or internal library?
            except:
                print('No HDU marked "image" found...')
            
            self.wavelength = wavelength
            self.fluxUnits  = fluxUnits
            self.filterName = filterName
            self.imageScale = imageScale
            self.imageSize  = imageSize
            self.images     = images
        #now return the instance so that it works as intended
        return self
            

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

        ''' docstring goes here '''

        raise NotImplementedError()
        
        ''' First take the model values (passed in) and compute synthetic Interferometry image '''
        modInterferImage = 1 #Needs to be added here. Take synthetic photometry on the Interferometry image plane.

        ''' then update the covariance matrix for the parameters passed in '''
        #skip this for now
        self.covMat = self.cov()
        
        ''' then compute the likelihood for each photometric point in a vectorised statement '''
        a = self.value - modInterferImage

        b = np.log(1./((2*np.pi)**(len(self.value)) * np.linalg.det(self.covMat))
            ) 
        #pass
        probFlux = b + ( -0.5 * ( np.matmul ( a.T, np.matmul(inv(self.covMat), a) ) ) )

        return probFlux


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

        ''' docstring goes here '''

        raise NotImplementedError()
        
        ''' First take the model values (passed in) and compute synthetic IFU data '''
        modCube = 1 #Needs to be added here. Generate model IFU data.

        ''' then update the covariance matrix for the parameters passed in '''
        #skip this for now
        self.covMat = self.cov()
        
        ''' then compute the likelihood for each photometric point in a vectorised statement '''
        a = self.value - modCube

        b = np.log(1./((2*np.pi)**(len(self.value)) * np.linalg.det(self.covMat))
            ) 
        #pass
        probFlux = b + ( -0.5 * ( np.matmul ( a.T, np.matmul(inv(self.covMat), a) ) ) )

        return probFlux
