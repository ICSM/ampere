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

        ''' inititalise covariance matrix as a diagonal matrix '''
        #self.covMat = np.diag(uncertainty**2)#np.diag(np.ones(len(uncertainty)))
        self.varMat = uncertainty * uncertainty[:,np.newaxis] #make a square array of sigma_i * sigma_j, i.e. variances
        self.covMat = np.diag(np.ones_like(uncertainty))
        a = self.covMat > 0
        self.covMat[a] = self.covMat[a] * self.varMat[a]# = np.diag(uncertainty**2)
        self.logDetCovMat = np.linalg.slogdet(self.covMat)[1]# / np.log(10.)
        print(self.logDetCovMat)
        if self.logDetCovMat == -np.inf:
            print("""The determinant of the covariance matrix for this dataset is 0.
            Please check that the uncertainties are positive real numbers """)
            print(self)
            exit()

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
        
        Future work: go through multiple libraries from different (user-defined) locations and import htem all
        '''
        
        if libName is None:
            print("No library given, using default pyphot filters")
            libDir = pyphot.__file__.strip('__init__.py')+'libs/'
            libName = libDir + 'synphot_nonhst.hd5' #PhIReSSTARTer.hd5'
        
        self.filterLibrary = pyphot.get_library(fname=libName)

    def filterNamesToPyphot(self, **kwargs):
        pyphotFilts = self.filterLibrary.get_library_content()
        #print(pyphotFilts)
        filtsOrig = self.filterName
        l = []
        for filt in self.filterName:
            l.append(filt in pyphotFilts)
        #try replacing colons and / with _
        #print(l)
        newTry = [filt.replace(':','_').replace('/','_').replace('WISE','WISE_RSR').replace('Spitzer','SPITZER') for filt in self.filterName]
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
        return 0

    def lnlike(self, theta, model, **kwargs):
        ''' docstring goes here '''
        
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
        a = self.value - modSed

        b = -0.5*len(self.value) * np.log(2*np.pi) - (0.5*self.logDetCovMat)
            #np.log(1./((2*np.pi)**(len(self.value)) * np.linalg.det(self.covMat))
            #) 
        probFlux = b + ( -0.5 * ( np.matmul ( a.T, np.matmul(inv(self.covMat), a) ) ) )
        return probFlux
        

    def cov(self, **kwargs):
        ''' 
        This routine populates a covariance matrix given some methods to call and parameters for them.

        For the moment, however, it does nothing.
        '''
        return self.covMat

    @classmethod
    def fromFile(cls, filename, format=None, **kwargs):
        ''' 
        Routine to generate photometry data object from a file containing said data
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
        raise NotImplementedError()
    
    def cov(self, theta, **kwargs):
        ''' 
        This routine populates a covariance matrix given some methods to call and parameters for them.

        For the moment, however, it does nothing.
        '''
        ''' Gradually build this up - for the moment we will assume a single squared-exponential kernel + an uncorrelated component '''

        ''' create a grid of positions/distances on which to build the matrix '''
        #we might be able to offload this step to an instance variable, as it shouldn't change between iterations...
        i, j = np.mgrid[:len(self.wavelength),:len(self.wavelength)]
        d = i - j

        ''' hardcode a squared-exponential kernel for the moment '''
        m = np.exp(-d**2. / (2.* theta[1]**2.))
        a = m < self.covarianceTruncation
        m[np.logical_not(a)] = 0. #overwrite small values with 0 to speed up some of the algebra

        covMat = (1-theta[0])*np.diag(np.ones_like(self.uncertainty)) + theta[0]*m
        self.covMat = covMat * self.varMat
        
        self.logDetCovMat = np.linalg.slogdet(self.covMat)[1]# / np.log(10.)
        #return self.covMat

    def lnprior(self, theta, **kwargs):
        try:
            scaleFac = theta[0]
        except IndexError: #Only possible if theta is scalar or can't be indexed
            scaleFac = theta
        if scaleFac > 0 and 0. < theta[1] < 1. and theta[2] > 0.:
            #print(scalefac)
            return norm.logpdf(np.log10(scaleFac), loc=0., scale = self.calUnc) + halfnorm.logpdf(theta[2], 0., 1.)
        return -np.inf

    def lnlike(self, theta, model, **kwargs):
        ''' docstring goes here '''
        
        ''' First take the model values (passed in) and compute synthetic Spectrum '''
        
        try:
            scaleFac = theta[0]
        except IndexError: #Only possible if theta is scalar or can't be indexed
            scaleFac = theta
            
        #print(self)
        #wavelength = self.wavelength
        #modSpec = model.modelFlux #
        modSpec = spectres(self.wavelength, model.wavelength, model.modelFlux) #For some reason spectres isn't cooperating :/ actually, looks like it was just a stupid mistake

        ''' then update the covariance matrix for the parameters passed in '''
        #skip this for now
        #self.covMat =
        self.cov(theta[1:])
        #import matplotlib.pyplot as plt
        #plt.imshow(self.covMat)
        #plt.show()
        
        ''' then compute the likelihood for each photometric point in a vectorised statement '''
        a = scaleFac*self.value - modSpec
        #if not np.all(np.isfinite(a)):
        #    print(a)
        #    print(modSpec)

        #make this a try: except OverflowError to protect against large spectra (which will be most astronomical ones...)?
        b = 0#np.log10(1./((np.float128(2.)*np.pi)**(len(self.value)) * np.linalg.det(self.covMat))
            #)

        b = -0.5*len(self.value) * np.log(2*np.pi) - (0.5*self.logDetCovMat) #less computationally intensive version of above
        #pass
        probFlux = b + ( -0.5 * ( np.matmul ( a.T, np.matmul(inv(self.covMat), a) ) ) )
        #print(((np.float128(2.)*np.pi)**(len(self.value))), np.linalg.det(self.covMat))
        #print(((np.float128(2.)*np.pi)**(len(self.value)) * np.linalg.det(self.covMat)))
        #print(b, probFlux)
        #exit()
        return probFlux 

    @classmethod
    def fromFile(cls, filename, format, **kwargs):
        ''' 
        Routine to generate spectrum data object from a file containing said data
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
            # note that there is a column called module. this has values 0.0 1.0 2.0 and 3.0 corresponding to SL1, SL2, LL1 and LL2 resp.
            # we may want consider them independent?
            #print(np.unique(table['module']))
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
            ''' If I'm interpreting the CASSIS data correectly, Module(SL) = 0, 1; Module(LL) = 2, 3 '''
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
