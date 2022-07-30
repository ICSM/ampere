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

from .data import Data


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

        print(photUnits)
        print(type(photUnits))

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
            print(len(photUnits))
        except TypeError:
            print(photUnits)
        try:
            print(len(value))
        except TypeError:
            print(value)
        

        #print(len(photUnits))
        #print(len(value))
        try:
            assert len(photUnits) == len(value)
        except AssertionError: #We have more than one unit entry, but not one per flux entry, raise an error and force the user to do something about it:
            if isinstance(photUnits, str):
                print("photUnits is a string")
                photUnits = [photUnits] * len(value)
            else:
                print("photUnits is weird")
                raise RunTimeError("The wrong number of unit entries appear to have been provided. Please check this and try again. You provided {0} units, but {1} fluxes. \n The fluxes are \n {2} \nand the units are \n {3}".format(len(photUnits), len(value), photunits, values))
        except TypeError: #only one unit was provided, let's forcibly turn it into an iterable
            print("photunits is not iterable")
            photUnits = len(value) * (photUnits,)
        else:
            if isinstance(photUnits, str):
                print("photUnits is a string")
                photUnits = [photUnits] * len(value)
            else:
                print("photUnits is very weird")

        print(photUnits)
        print(type(photUnits))
                       
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
        
    #This should be updated to only get recalculated if anything ever changes
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

    def simulate(self, theta, model, **kwargs):
        """ Simulate photometry, given a model result"""

        #First, we compute the model photometry
        filts, modSed = pyphot.extractPhotometry(model.spectrum["wavelength"],
                                                 model.spectrum["flux"],
                                                 self.filters,
                                                 Fnu = True,
                                                 absFlux = False,
                                                 progress=False
            )

        #Then, we add some noise
        rng = np.random.default_rng() #Optimise this, no need to re-create object each time, plus need seed to be stored so it can be saved
        simulated_data = rng.multivariate_normal(modSed[self.mask], self.covMat) #We draw one sample from the multivariate normal distribution with mean = model photometry and covariance = data covariances
        
        #now we can return the simulated data
        return simulated_data

    #@property
    #def 

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

class LineStrengths(Photometry):
    """A class for the special case of line strengths (e.g. EWs or integrated intensities)

    At present, this is a placeholder.
    """
