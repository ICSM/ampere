from __future__ import print_function

import numpy as np
from numpy import ma
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
from scipy.stats import rv_continuous, norm, halfnorm
#from scipy.linalg import inv
from numpy.linalg import inv

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

    _ismasked = False

    def __init__(**kwargs):
        pass

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
       raise NotImplementedError()  # switched off by sascha 19/04/2021

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

    def setPlotParameters(self):
        pass

    def plot(self):
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
        if self._ismasked:
            self.mask = np.logical_and(mask, self.mask)
        else:
            self.mask = mask
            self._ismasked = True
#        try:
#            self.mask = np.logical_and(mask, self.mask)
#        except NameError:
#            self.mask = mask
            
        #now we need to create a mask for the covariance matrix
        #The outer product does what we want, producing a matrix which has elements such that cov_mask[i,j] = mask[i] * mask[j]
        #This produces the right answer because boolean multiplication is treated as an AND operation in python

        #self.cov_mask = np.outer(self.mask, self.mask)

        #now we need to update the covariance matrix by extracting the unmasked elements from the original one.
        #However, it has to be reshaped because numpy always returns 1D arrays when using boolean masks
        #self.covMat = self.covMat_orig[self.cov_mask].reshape((np.sum(mask), np.sum(mask)))
        self.cov(None)

        pass

#    def maskNaNs(self, **kwargs):
#        '''Method to generate a mask which blocks NaNs in the data.
#        '''
#
#        mask = np.logical_and(np.isfinite(self.value), np.isfinite(self.uncertainty))
        
        #Now we add check to make sure that if masks have previously been defined we don't overwrite th\em, and only accept values 
        #that pass both masks. Otherwise, we define a mask.
#        try:
#            self.mask = np.logical_and(mask, self.mask)
#        except NameError:
#            self.mask = mask
        
        #now we need to create a mask for the covariance matrix
        #The outer product does what we want, producing a matrix which has elements such that cov_mask[i,j] = mask[i] * mask[j]
        #This produces the right answer because boolean multiplication is treated as an AND operation in python
#        self.cov_mask = np.outer(self.mask, self.mask)
        
    def unmask(self, **kwargs):
        '''A method to overwrite previous masks with True in case something goes wrong
        
        '''
        if self._ismasked:
            mask = np.ones_like(self.mask, dtype=np.bool)
        else:
            mask = np.ones_like(self.value, dtype=np.bool)
            
        self.mask = mask
        self.cov_mask = np.outer(self.mask, self.mask)

        self.cov()

    
    
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

    #Create a dictionary of default plotting parameters for Photometry objects
    plotParams={"alpha": 1.0,
                     "label": "Photometry",
                     "linestyle":'None',
                     "marker": 'o',
                     "markeredgecolor": 'black',
                     "markerfacecolor": 'orange',
                     "markersize": 2.,
                     "zorder": 10
        }

    _hasNoiseModel = False

    def __init__(self, filterName, value, uncertainty, photUnits,
                 bandUnits=None, libName = None, label = "Photometry", **kwargs):
        self.filterName = filterName

        ''' setup pyphot for this set of photometry '''
        self.pyphotSetup(libName)
        print(self.filterName.astype('str'))
        #newTry = [str(filt).replace(':','_').replace('/','_').replace('WISE','WISE_RSR').replace('Spitzer','SPITZER') for filt in self.filterName]
        self.filterNamesToPyphot()

        self.label = label
        self.plotParams["label"] = label

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
                
        mags = self.fluxUnits == 'mag'
        #print(mags.__repr__())
        #pyphot returns a list of filters, this means this nice boolean masking doesn't work :(
        zeropoints = np.zeros_like(value)
        for i in range(len(mags)):
            if mags[i]:
                zeropoints[i] = filters[i].Vega_zero_Jy.magnitude
                value[i] = zeropoints[i]*10^(-0.4*value[i])
                uncertainty[i] = value[i] - zeropoints*10^(-0.4*(value[i]+uncertainty[i]))
        
        try:
            assert len(photUnits) == len(value)
        except AssertionError: #We have more than one unit entry, but not one per flux entry, raise an error and force the user to do something about it:
            if isinstance(photUnits, str):
                photUnits = [photUnits] * len(value)
            else:
                raise RunTimeError("The wrong number of unit entries appear to have been provided. Please check this and try again. You provided {0} units, but {1} fluxes. \n The fluxes are \n {2} \nand the units are \n {3}".format(len(photUnits), len(value), photunits, values))
        except TypeError: #only one unit was provided, let's forcibly turn it into an iterable
            photUnits = len(value) * (photUnits,)
                       
        #identify values in milliJansky, convert to Jy
        uconv = np.array([u.Jy.to(pU) for pU in photUnits])
        self.uncertainty = uncertainty / uconv
        self.value = value / uconv

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
        for i in range(len(self.filterName[self.mask])):
            l.append(
                '{:<{nFilt}} {:>{nWave}.2e} {:>{nVal}.2e} {:>{nUnc}.2e}'.format(
                    self.filterName[self.mask][i],self.wavelength[self.mask][i],self.value[self.mask][i],self.uncertainty[self.mask][i],
                nFilt = nFilt, nWave = nWave, nVal = nVal, nUnc = nUnc
                )
            )

        ''' finally we might want a way to output the coviariance matrix, although that might only be useful for plotting '''
        return l
    
    def __repr__(self, **kwargs):
    #    raise NotImplementedError()   modified on 19/04/2021 by sascha
        return self.__str__()
        
    
    def geef_data(self, **kwargs):
        return(self.filterName,self.wavelength,self.value, self.uncertainty)

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
            newTry = [filt.decode("utf-8").replace(':','_').replace('/','_').replace('WISE','WISE_RSR').replace('Spitzer','SPITZER') for filt in self.filterName]
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

    def prior_transform(self, u, **kwargs):
        """Transform from uniform RVs to the prior of any nuisance parameters. 

        Since this implementation has no nuisance parameters, it does nothing."""
        return None
        
            
    def lnprior(self, theta, **kwargs):
        """Return the prior of any nuisance parameters. 

        Since this implementation has no nuisance parameters, it does nothing."""
        return 0

    def lnlike(self, theta, model, **kwargs):
        '''Compute the likelihood of the photometry given the model. 

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
        if self.logDetCovMat == -np.inf or self.signDetCovMat < 0:
            return -np.inf
        
        ''' then compute the likelihood for each photometric point in a vectorised statement '''
        a = self.value[self.mask] - modSed[self.mask] 

        b = -0.5*len(self.value[self.mask]) * np.log(2*np.pi) - (0.5*self.logDetCovMat)
            #np.log(1./((2*np.pi)**(len(self.value)) * np.linalg.det(self.covMat))

            #)

        #covMatmask = np.reshape(self.covMat[self.cov_mask], np.shape(self.covMat))
        probFlux = b + ( -0.5 * ( np.matmul ( a.T, np.matmul(inv(self.covMat), a) ) ) )

        if np.isnan(probFlux):
            print("NaN probability")
            print(theta)
            print(b)
            print(inv(self.covMat))
            print(self.value[self.mask])
            print(modSed[self.mask])
            print(a)
            print(np.matmul(inv(self.covMat), a))
            print(np.matmul ( a.T, np.matmul(inv(self.covMat), a) ) )
            return -np.inf #hack for now so we can see how often this occurs and hopefully troubleshoot it!

        return probFlux
        

    def cov(self, *args, **kwargs):
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
        self.varMat = np.outer(self.uncertainty[self.mask], self.uncertainty[self.mask]) #self.uncertainty[self.mask] * self.uncertainty[self.mask][:,np.newaxis] #make a square array of sigma_i * sigma_j, i.e. variances
        self.covMat = np.diag(np.ones_like(self.uncertainty[self.mask]))
        a = self.covMat > 0
        self.covMat[a] = self.covMat[a] * self.varMat[a]# = np.diag(uncertainty**2)
        #self.covMat = np.reshape(#self.covMat[self.cov_mask], np.shape(self.covMat))
        self.signDetCovMat, self.logDetCovMat = np.linalg.slogdet(self.covMat)# / np.log(10.)
#        self.logDetCovMat = np.linalg.slogdet(self.covMat[self.cov_mask])[1]# / np.log(10.)
        #print(self.logDetCovMat)
        if self.logDetCovMat == -np.inf: #This needs to be updated to raise an error!

            print("""The determinant of the covariance matrix for this dataset is 0.
            Please check that the uncertainties are positive real numbers """)
            print(self)
            exit()
        elif self.signDetCovMat < 0:
            print("""The determinant of the covariance matrix for this dataset is negative.
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


    def setPlotParams(self,**kwargs):
        '''
            Routine to set up the plotting parameters for any photometry data set, and optionally plot (and show and/or save) the data. 
            A limited set of keywords can be set by passing **kwargs to the function.
        '''

        print("Setting plotting parameters for photometry data set.")

        #update plotParams with any keywords passed via **kwargs
        #This only changes things that were passed in, so anything else will retain its default
        for key, value in kwargs.items():
            self.plotParams[key] = value
        
        #if doPlot == True:
        #    import matplotlib.pyplot as plt
        #    
        #    fig, ax = plt.subplots(1,1,figsize=(6, 8))
        #    ax.set_xtitle(r"Wavelength ($\mu$m)")
        #    ax.set_ytitle(r"Flux density (mJy)")
        #    ax.errorbar(self.wavelength,self.value,xerr=None,yerr=self.uncertainty,\
        #                linestyle=self.linestyle,marker=self.marker,mec=self.mec,mfc=self.mfc,\
        #                color=self.mec,ecolor=self.mec,alpha=self.alpha,legend=self.label)
        #    plt.legend("lower right",fontsize="small")
        #    plt.tight_layout()

        #    if showPlot == True:
        #        plt.show()
        #    if savePlot == True: 
        #        fig.savefig("dataset.png",dpi=200,overwrite=True)
        #    plt.close()
        #    plt.clf()


    def plot(self, fig = None, ax = None, unmask=False,
             doPlot=True,savePlot=False,showPlot=False, **kwargs):
        
        self.setPlotParams(**kwargs)

        if doPlot:
        
            if ax is not None:
                #Easiest case, we just have to use ax.errorbar to plot onto existing axes
                ax.errorbar(self.wavelength[self.mask], self.value[self.mask],
                            yerr=self.uncertainty[self.mask],
                            **self.plotParams, **kwargs)
            else:
                if fig is not None:
                    #Now we have a figure but no axes.
                    #Therefore we will add some new axes to the figure and proceed
                    if len(fig.get_axes()) == 0:
                        ax = fig.add_subplot(111)
                    else:
                        #For now, assume that this is supposed to be added to the last axes if some exist
                        ax = fig.get_axes()[-1]
                    ax.errorbar(self.wavelength[self.mask], self.value[self.mask],
                                yerr=self.uncertainty[self.mask],
                                **self.plotParams, **kwargs)
                else: #no axis or figure, let's create everything
                    import matplotlib.pyplot as plt
                    fig, ax = plt.subplots(1,1,figsize=(6, 8))
                    ax.set_xtitle(r"Wavelength ($\mu$m)")
                    ax.set_ytitle(r"Flux density (mJy)")
                    ax.errorbar(self.wavelength[self.mask],self.value[self.mask],
                                yerr=self.uncertainty[self.mask],
                                **self.plotParams, **kwargs)
                    plt.legend("lower right",fontsize="small")
                    plt.tight_layout()
                if savePlot:
                    fig.savefig(self.label+"plot.png", dpi = 200, overwrite=True)
                if showPlot:
                    plt.show()
        pass
                   

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
    bandUnits : {'um', 'AA', 'nm'}, optional, string, scalar or array-like
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

    plotParams={"alpha": 1.0,
                     "color": "blue",
                     "drawstyle":"steps-mid",
                     "label": "Spectrum",
                     "linestyle":'-',
                     "linewidth":1.5,
                     "zorder": 9,
                }
    plotParamsSample={"sample color": "cyan",
                     "sample alpha": 0.6,
                     "sample linewidth": 1.,
                     "sample zorder": 5
        }

    _hasNoiseModel = True

    def __init__(self, wavelength, value, uncertainty, bandUnits, fluxUnits,
                 freqSpec = False, calUnc = None, covarianceTruncation = 1e-3,
                 scaleLengthPrior = None, scaleFacPrior = None, label="Spectrum",
                 resampleMethod="exact",
                 **kwargs):
        #self.filterName = filterName #Telescope/Instrument cf photometry
        self.type = 'Spectrum'
        self.label = label
        self.plotParams["label"] = label

        #Wavelength conversion
        self.frequency = freqSpec #True/False? Wavelength or Frequency spectrum
        
        if freqSpec == 'True': #If freq, convert to wavelength in micron
            wavelength = (const.c / wavelength).to('um')
        
        self.bandUnits = bandUnits
        
        # as long as bandUnits has the dimensions of length, the wavelength
        # is computed in um.
        if u.m.is_equivalent(bandUnits):
            wavelength = wavelength / u.um.to(bandUnits)
        else:
            mesg = """The value for bandUnits must have dimensions of length!
            Allowed values: 'AA', 'um', 'nm', etc."""
            raise ValueError(mesg)
        
        self.wavelength = wavelength #Will be the grid of wavelength/frequency
                                     #using numpy.ma.asarray()? to make it a MaskedArray so that we can easily
                                     #mask out unwanted parts of the spectrum along with its covariances
        
        self.fluxUnits = fluxUnits#specFluxUnits

        if u.Jy.is_equivalent(fluxUnits):
            uconv = u.Jy.to(fluxUnits)
            value = value / uconv
            uncertainty = uncertainty / uconv
        elif (u.W / u.m**3).is_equivalent(fluxUnits):
            # first convert flux and uncertainty to W/m^3
            uconv = (u.W / u.m**3).to(fluxUnits)
            value = value / uconv
            uncertainty = uncertainty / uconv
            # now convert flux and uncertainty to Jy
            conv_fac = (wavelength * u.um)**2 / const.c
            value = (value * conv_fac).to('Jy')
            uncertainty = (uncertainty * conv_fac).to('Jy')
        else:
            raise NotImplementedError()

        self.value = value #Should always be a flux unless someone is peverse, as with wavelengths using MaskedArrays?
        self.uncertainty = uncertainty #Ditto #ma.asarray()?

        ''' inititalise covariance matrix as a diagonal matrix '''
        self.varMat = uncertainty * uncertainty[:,np.newaxis] #make a square array of sigma_i * sigma_j, i.e. variances
        self.covMat = np.diag(np.ones_like(uncertainty)) #ma.asarray()
        #self.covMat.mask = np.logical_not(self.covMat > 0) #The mask for a MaskedArray is backwards compared to boolean indexing!
        #                                                   #By applying this mask we minimise the operations that need to be done, but this might really be overkill
        a = self.covMat > 0

        self.covMat = self.covMat * self.varMat # = np.diag(uncertainty**2)
        self.logDetCovMat = np.linalg.slogdet(self.covMat)[1]# / np.log(10.)
        #print(self.logDetCovMat)


        ''' Assume default of 10% calibration uncertainty unless otherwise specified by the user '''
        if calUnc is None:
            self.calUnc = 0.010 #beware! This is a variance!
        else:
            self.calUnc = calUnc

        self.npars = 3 #this will have to change in future to 1 + number of parameters required by GPs for covMat

        if resampleMethod == "exact":
            self.resampler = spectres
        elif resampleMethod=="fast":
            self.resampler = np.interp

        #Set up the objects for some of the priors:
        if isinstance(scaleFacPrior, rv_continuous):
            #The user has provided a prior, use it!
            self.scaleFacPrior = scaleFacPrior
        elif isinstance(scaleFacPrior, float) and calUnc is None:
            #The user has provided the calibration uncertainty here rather than in calUnc
            self.scaleFacPrior = norm(loc=0., scale = scaleFacPrior)
        elif isinstance(scaleFacPrior, float):
            #The user has provided a float here as well as a number of calUnc, we'll assume this is the logarithm of the bias on the scale factor
            self.scaleFacPrior = norm(loc=scaleFacPrior, scale = self.calUnc)
        elif isinstance(scaleFacPrior, (list, tuple, np.ndarray)) and len(scaleFacPrior) == 2: #Consider changing the first condition to isinstance(scaleFacPrior, (collections.Sequence, np.ndarray)) ?
            #Two numbers, they should be the loc and scale parameters
            self.scaleFacPrior = norm(loc=scaleFacPrior[0], scale = scaleFacPrior[1])
        else:
            self.scaleFacPrior = norm(loc=0., scale = self.calUnc)

        #Now the scale length of the noise model
        if isinstance(scaleLengthPrior, rv_continuous):
            #The user has provided a prior, use it!
            self.scaleLengthPrior = scaleLengthPrior
        elif isinstance(scaleLengthPrior, float):
            #The user has provided the width of the prior
            self.scaleLengthPrior = halfnorm(loc=0., scale = scaleLengthPrior) #Using the half normal results in infinite bounds on the length scale, and hence there is always room to improve the likelihood by increasing the length scale. It may be necessary to alter this to a truncated normal distribution.
        else:
            self.scaleLengthPrior = halfnorm(loc=0., scale = 1)
        
        #self.scaleLengthPrior = halfnorm(loc = 0, scale= 1.) #.logpdf(theta[2], 0., 1.)

        self.covarianceTruncation = covarianceTruncation #if defined, covariance matrices will be truncated when their values are less than this number, to minimise the number of operation

        self.theta_last = [0.,1.] #Set some default theta values that correspond to uncorrelated noise so that selectWaves() will work without unwanted side effects

        self.selectWaves()

    def __call__(self, **kwargs):
        raise NotImplementedError()

    def setResampler(self, resampleMethod = "exact", **kwargs):
        """Change the method used to resample the model spectra to the observed wavelengths

        Parameters
        ----------
        resampleMethod : {'exact', 'fast', callable}, str or callable

        Notes
        -----
        'exact' uses spectres for exact, flux-conserving resampling. This is 
        important if your model contains features that are unresolved at the 
        sampling of the observed data, so that the flux is correct. As a result, 
        this is the default method, although it is relatively slow. We are 
        working on making an optimised version available which uses numba for 
        just-in-time compilation.

        When all features in the model are resolved in the observations, however, 
        this is overkill. For those cases, interpolation of the model to the 
        observed wavelengths is usualy sufficient, and so the much faster 
        numpy.interp can be used by selecting the "fast" method.

        As we are sure this will not cover all potential use cases, we provide 
        the option for a user-supplied callable. This must have inputs in the 
        same format as spectres and numpy.interp, so it must be called as
        new_fluxes = f(observed_wavelengths, model_wavelengths, model_fluxes)
    
    
        Examples
        --------
        To switch on fast resampling, simply call
        >>> spec = Spectrum(wavelength_array, flux_array, uncertainty_array, "um", "Jy")
        >>> spec.setResampler(resampleMethod = "fast")

        and to switch back to exact
        >>> spec.setResampler(resampleMethod = "exact")
        
        if you have a specific function you want to use pass it as
        >>> from package import resampler
        >>> spec.setResampler(resampleMethod = resampler)
        
        """

        if resampleMethod == "exact":
            self.resampler = spectres
            print("Using exact resampling")
        elif resampleMethod=="fast":
            self.resampler = np.interp
            print("Using fast resampling")
        elif callable(resampleMethod):
            self.resampler = resampleMethod
            print("Using user-specified resampling")
        #elif resampleMethod=="exactcompiled":
        #    #Raise a warning here, this is going to be experimental when it's implemented
        #    self.resampler = spectres_numba
        else:
            raise ValueError("resampleMethod must be a callable or one of \"exact\" or \"fast\", but \"{0}\" was used".format(resampleMethod))
    
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
        for i in range(len(self.wavelength[self.mask])):
            l.append(
                '{:>{nWave}.2e} {:>{nVal}.2e} {:>{nUnc}.2e}'.format(
                    self.wavelength[self.mask][i],self.value[self.mask][i],self.uncertainty[self.mask][i],
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

        if theta is None:
            theta = self.theta_last
        #we might be able to offload this step to an instance variable, as it shouldn't change between iterations...
        i, j = np.mgrid[:len(self.wavelength),:len(self.wavelength)]# np.mgrid[:len(self.wavelength[self.mask]),:len(self.wavelength[self.mask])]

        i, j = np.meshgrid(self.wavelength[self.mask], self.wavelength[self.mask]) #Now it's wavelengths, not pixels!
        d = i - j #These lines need to be switched to *wavelength* space, not *pixel* space



        ''' hardcode a squared-exponential kernel for the moment '''
        m = np.exp(-(d**2.) / (2.* theta[1]**2.)) 
        #m = m / np.max(m) #at this point we have a correlation matrix, which should always have ones on the diagonal unless something strange is going on
        a = m < self.covarianceTruncation
        m[a] = 0.#[np.logical_not(a)] = 0. #overwrite small values with 0 to speed up some of the algebra

        #self.varMat = #uncertainty[self.mask] * uncertainty[self.mask][:,np.newaxis] #make a square array of sigma_i * sigma_j, i.e. variances
        self.varMat = np.outer(self.uncertainty[self.mask], self.uncertainty[self.mask])
        #covMat = (1-theta[0])*np.diag(np.ones_like(self.uncertainty[self.mask])) + theta[0]*m
        covMat = (1-theta[0])*np.diag(np.ones_like(self.uncertainty[self.mask])) + theta[0]*m
        self.covMat = covMat * self.varMat

        #covMatmask = np.reshape(self.covMat[self.cov_mask], (self.value[mask].shape[0], self.value[mask].shape[0]))

        self.signDetCovMat, self.logDetCovMat = np.linalg.slogdet(self.covMat)
        if self.logDetCovMat == -np.inf: #This needs to be updated to raise an error! Or at least force the likelihood to be zero
            print("""The determinant of the covariance matrix for this dataset is 0.
            Please check that the uncertainties are positive real numbers """)
            print(self)
            print(theta)
            print(self.signDetCovMat,self.logDetCovMat)
            print(self.covMat)
            exit()
        #elif self.signDetCovMat < 0:
        #    print("""The determinant of the covariance matrix for this dataset is negative.
        #    Please check that the uncertainties are positive real numbers """)
        #    print(self)
        #    print(theta)
        #    print(self.signDetCovMat,self.logDetCovMat)
        #    print(self.covMat)
        #    print(d)
        #    print(m)
        #    print(covMat)
        #    exit()
        self.theta_last = theta
        #return self.covMat
        
    def prior_transform(self, u, **kwargs):
        """Transform from uniform RVs to the prior of any nuisance parameters. 

        This implementation has 3 nuisance parameters. These are 1) a scale factor, 
        which multiplies the observed spectrum to include the uncertainty on 
        absolute calibration, 2) the strength of the correlated component of the
        covariances and 3) the scale length (standard deviation) of the Gaussian
        kernel for the correlated component of the covariances. All three parameters
        are passed in with a single sequence.
        """

        theta = np.zeros_like(u)
        theta[0] = 10**self.scaleFacPrior.ppf(u[0])
        theta[1] = u[1]
        theta[2] = self.scaleLengthPrior.ppf(u[2])
        return theta
    
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
        if scaleFac > 0 and 0. <= theta[1] <= 1. and theta[2] > 0.:
            #print(scalefac)
            return self.scaleFacPrior.logpdf(np.log10(scaleFac)) + self.scaleLengthPrior.logpdf(theta[2])#norm.logpdf(np.log10(scaleFac), loc=0., scale = self.calUnc) + halfnorm.logpdf(theta[2], 0., 1.)
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
        #print(model.wavelength)
        modSpec = self.resampler(self.wavelength, model.wavelength, model.modelFlux)[self.mask] #For some reason spectres isn't cooperating :/ actually, looks like it was just a stupid mistake
        ''' then update the covariance matrix for the parameters passed in '''
        #skip this for now
        #self.covMat =
        self.cov(theta[1:])
        #if self.logDetCovMat == -np.inf:
        #    return -np.inf
        if self.signDetCovMat < 0: #Hyperparameters for covariance matrix result in a negative determinant
            return -np.inf         #Therefore this isn't a valid covariance matrix and we can't use this parameter combination
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
        #covMatmask = np.reshape(self.covMat[self.cov_mask], np.shape(self.covMat))
        probFlux = b + ( -0.5 * ( np.matmul ( a.T, np.matmul(inv(self.covMat), a) ) ) )
        if np.isnan(probFlux):
            print("NaN probability")
            print(theta)
            print(b)
            print(inv(self.covMat))
            print(self.value[self.mask])
            print(modSpec)
            print(model.modelFlux)
            print(a)
            print(np.matmul(inv(self.covMat), a))
            print(np.matmul ( a.T, np.matmul(inv(self.covMat), a) ) )
            return -np.inf #hack for now so we can see how often this occurs and hopefully troubleshoot it!
        #print(((np.float128(2.)*np.pi)**(len(self.value))), np.linalg.det(self.covMat))
        #print(((np.float128(2.)*np.pi)**(len(self.value)) * np.linalg.det(self.covMat)))
        #print(b, probFlux)
        #exit()
        return probFlux 

    @classmethod
    def fromFile(cls, filename, format, filetype=None, keywords=None, **kwargs):
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
            #should we divide this in four chunks: module 0 = SL2, module 1 = SL1, module 2 = LL2, module 3 is LL1?
            chunks[sl] = 1.
            chunks[ll] = 2.
            table['chunk'] = chunks

            #normalise the column names (wavelength,flux,uncertainty)
            table.rename_column('error (RMS+SYS)','uncertainty')
            
            ''' If I'm interpreting the CASSIS data correctly, Module(SL) = 0, 1; Module(LL) = 2, 3 '''
            #a = table['module'] > 1.
            #tablell = table[a]#.sort(keys='wavelength')
            #tablell.sort(keys='wavelength')
            #tablell.pprint()
            #table = tablell
            table.sort(keys='wavelength')
            
        if format == 'SPITZER-YAAAR_OPTDIFFHR':
            #filename = 'Testdata/cassis_yaaar_spcfw_14203136t.fits'
            hdul = fits.open(filename)
            hdu = hdul[0]
            header=hdu.header
            data = hdu.data
            table = Table(data,names=[header['COL01DEF'],header['COL02DEF'],header['COL03DEF'],header['COL04DEF'],header['COL05DEF'],header['COL06DEF'],header['COL07DEF'],header['COL08DEF']])
            table['wavelength'].unit='um'
            table['flux'].unit='Jy'
            table['flux_error'].unit='Jy'
            
            #table['error (RMS)'].unit='Jy'
            #table['error (SYS)'].unit='Jy'
            #table['offset uncertainty (CAL)'].unit='Jy'
            #table['sky'].unit='Jy'
            #table['sky error'].unit='Jy'
            mask = np.logical_and.reduce([np.isfinite(c) for c in table.columns.values()]) #require the elements to be non-NaNs

            table = table[mask]
            wavelengths = np.array(table['wavelength'].data)
            orders = np.array(table['IRS_order'].data)
            chunks = np.zeros_like(table['module'].data)
            order_array=np.arange(11,21,1)
            
            for i_order,order in enumerate(order_array):
                wl_sel=wavelengths[order==orders]
                nums=np.where(order==orders)[0]
                breaknum=np.where(np.diff(wl_sel)>5)[0][0]+1
                chunks[nums[0:breaknum]]=1 #sh
                chunks[nums[breaknum::]]=2 #lh

            table['chunk'] = chunks

            #normalise the column names (wavelength,flux,uncertainty)
            table.rename_column('flux_error','uncertainty')
            
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
        specList = cls.fromTable(table, **kwargs)
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
            self.__init__(table['wavelength'][selection].data, value[selection], uncertainty[selection], bandUnits, photUnits, **kwargs) #Also pass in flux units
            specList.append(self)
        return specList #self

    #def setPlotParameters(self,doPlot=False,savePlot=False,showPlot=True,**kwargs):
    #    '''
    #        Routine to set up the plotting parameters for any spectroscopy data set, and optionally plot (and show and/or save) the data. 
    #        A limited set of keywords can be set by passing **kwargs to the function.
    #    '''

    #    print("Setting plotting parameters for spectroscopy data set.")

    #    self.label = "Spectroscopy"

    #    try:
    #        self.marker = marker
    #    except:
    #        self.marker = "."

    #    try:
    #        self.mec = mec
    #    except:
    #        self.mec = "blue"
        
    #    try:
    #        self.mfc = mfc
    #    except:
    #        self.mfc = "dodgerblue"
        
    #    try:
    #        self.linestyle = linestyle
    #    except:
    #        self.linestyle = "-"
        
    #    try:
    #        self.alpha = alpha
    #    except:
    #        self.alpha = 0.1
        
    #    if doPlot == True:
    #        import matplotlib.pyplot as plt
            
    #        fig, ax = plt.subplots(1,1,figsize=(6, 8))
    #        ax.set_xtitle(r"Wavelength ($\mu$m)")
    #        ax.set_ytitle(r"Flux density (mJy)")

    #        nk = len(self.specList)
    #        k = 1
    #        for spec in self.specList:
    #            ax.errorbar(spec.wavelength,spec.value,xerr=None,yerr=spec.uncertainty,\
    #                        linestyle=self.linestyle,marker=self.marker,mec=self.mec,mfc=self.mfc,\
    #                        color=self.mec,ecolor=self.mec,alpha=self.alpha,legend=self.label+"_chunk_"+k)
    #            k += 1
    #        plt.legend("lower left",fontsize="small")
    #        plt.tight_layout()

    #        if showPlot == True:
    #            plt.show()
    #        if savePlot == True: 
    #            fig.savefig("dataset.png",dpi=200,overwrite=True)
    #        plt.close()
    #        plt.clf()


    def setPlotParams(self,**kwargs):
        '''
            Routine to set up the plotting parameters for any photometry data set, and optionally plot (and show and/or save) the data. 
            A limited set of keywords can be set by passing **kwargs to the function.
        '''

        print("Setting plotting parameters for photometry data set.")

        #update plotParams with any keywords passed via **kwargs
        #This only changes things that were passed in, so anything else will retain its default
        for key, value in kwargs.items():
            self.plotParams[key] = value
        
        #if doPlot == True:
        #    import matplotlib.pyplot as plt
        #    
        #    fig, ax = plt.subplots(1,1,figsize=(6, 8))
        #    ax.set_xtitle(r"Wavelength ($\mu$m)")
        #    ax.set_ytitle(r"Flux density (mJy)")
        #    ax.errorbar(self.wavelength,self.value,xerr=None,yerr=self.uncertainty,\
        #                linestyle=self.linestyle,marker=self.marker,mec=self.mec,mfc=self.mfc,\
        #                color=self.mec,ecolor=self.mec,alpha=self.alpha,legend=self.label)
        #    plt.legend("lower right",fontsize="small")
        #    plt.tight_layout()

        #    if showPlot == True:
        #        plt.show()
        #    if savePlot == True: 
        #        fig.savefig("dataset.png",dpi=200,overwrite=True)
        #    plt.close()
        #    plt.clf()


    def plot(self, fig = None, ax = None, unmask=False,
             doPlot=True,savePlot=False,showPlot=False,
             **kwargs):
        
        self.setPlotParams(**kwargs)
        if "steps-" in self.plotParams["drawstyle"]:
            s = self.plotParams["drawstyle"].replace("steps-","")

        if doPlot:
        
            if ax is not None:
                #Easiest case, we just have to use ax.errorbar to plot onto existing axes
                ax.plot(self.wavelength[self.mask], self.value[self.mask],
                        #yerr=self.uncertainty[self.mask],
                        **self.plotParams, **kwargs)
                ax.fill_between(self.wavelength[self.mask],
                                self.value[self.mask] - self.uncertainty[self.mask],
                                self.value[self.mask] + self.uncertainty[self.mask],
                                alpha = self.plotParams["alpha"]/2,
                                label = self.plotParams["label"]+" 68\% confidence band",
                                step = s,
                                zorder = self.plotParams["zorder"] - 1,
                                **{k:v for k,v in self.plotParams.items() if k not in ["drawstyle", "zorder", "alpha", "label"]},
                                **kwargs)
            else:
                if fig is not None:
                    #Now we have a figure but no axes.
                    #Therefore we will add some new axes to the figure and proceed
                    if len(fig.get_axes()) == 0:
                        ax = fig.add_subplot(111)
                    else:
                        #For now, assume that this is supposed to be added to the last axes if some exist
                        ax = fig.get_axes()[-1]
                    ax.plot(self.wavelength[self.mask], self.value[self.mask],
                            #yerr=self.uncertainty[self.mask],
                            **self.plotParams, **kwargs)
                    ax.fill_between(self.wavelength[self.mask],
                                    self.value[self.mask] - self.uncertainty[self.mask],
                                    self.value[self.mask] + self.uncertainty[self.mask],
                                    alpha = self.plotParams["alpha"]/2,
                                    label = self.plotParams["label"]+" 68\% confidence band",
                                    step = s,
                                    zorder = self.plotParams["zorder"] - 1,
                                    **{k:v for k,v in self.plotParams.items() if k not in ["drawstyle", "zorder", "alpha", "label"]},
                                    **kwargs)
                else: #no axis or figure, let's create everything
                    import matplotlib.pyplot as plt
                    fig, ax = plt.subplots(1,1,figsize=(6, 8))
                    ax.set_xtitle(r"Wavelength ($\mu$m)")
                    ax.set_ytitle(r"Flux density (mJy)")
                    ax.plot(self.wavelength[self.mask],self.value[self.mask],
                            #yerr=self.uncertainty[self.mask],
                            **self.plotParams, **kwargs)
                    ax.fill_between(self.wavelength[self.mask],
                                    self.value[self.mask] - self.uncertainty[self.mask],
                                    self.value[self.mask] + self.uncertainty[self.mask],
                                    alpha = self.plotParams["alpha"]/2,
                                    label = self.plotParams["label"]+" 68\% confidence band",
                                    step = s,
                                    zorder = self.plotParams["zorder"] - 1,
                                    **{k:v for k,v in self.plotParams.items() if k not in ["drawstyle", "zorder", "alpha", "label"]},
                                    **kwargs)
                    plt.legend("lower right",fontsize="small")
                    plt.tight_layout()
                if savePlot:
                    fig.savefig(self.label+"plot.png", dpi = 200, overwrite=True)
                if showPlot:
                    plt.show()

    def plotRealisation(self, s, ax=None, **kwargs):

        localplotParams = {k:v for k,v in self.plotParams.items() if k not in ["color", "linewidth", "zorder", "alpha", "label"]}
        waves = self.wavelength[self.mask]
        values = s[0]*self.value[self.mask] #only update for scaling, not sure how to include covariance info or if it's even worth it.
        ax.plot(waves, values,
                #yerr=self.uncertainty[self.mask],
                linewidth = self.plotParamsSample["sample linewidth"],
                zorder = self.plotParamsSample["sample zorder"],
                label = self.plotParams["label"]+" posterior samples",
                color = self.plotParamsSample["sample color"],
                alpha = self.plotParamsSample["sample alpha"],
                **localplotParams, **kwargs)

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
