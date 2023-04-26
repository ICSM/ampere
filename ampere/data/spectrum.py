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
from astropy import units
from astropy.io import ascii
from spectres import spectres
from scipy.stats import rv_continuous, norm, halfnorm
#from scipy.linalg import inv
from numpy.linalg import inv

from .data import Data

class Spectrum(Data):
    '''A class to represent 1D spectra data objects and their properties

    


    Parameters
    ----------
    wavelength : float, array-like or Quantity
        The wavelengths of the spectrum. Required if value is *not* a Spectrum1D
    value : float, array-like, Quantity or specutils.Spectrum1D
        The fluxes of the spectrum. If it is a Spectrum1D, it is assumed to contain the spectral axis and uncertainties along with their units
    uncertainty : float, array-like or Quantity
        The uncertainty on the fluxes. Required if value is *not* a Spectrum1D
    fluxUnits : string or astropy.units.core.Unit, scalar
        The units of the fluxes. If a string is passed, it will be interpreted with astropy.units.Unit. Required if either value is not a Spectrum1D and either value or uncertainty are not a Quantity
    bandUnits : string or astropy.units.core.Unit, scalar
        The units of the wavelengths. If a string is passed, it will be interpreted with astropy.units.Unit. Required if value is not a Spectrum1D and if wavelength is not a Quantity
    freqSpec : bool, default False
        If the input is wavelength (False, default) of Frequncy (True). Currently only works for frequncies in GHz, otherwise input must be wavelength (deprecated)
    calUnc : float, optional, default 0.01
        Parameter controlling width of prior on the scale factor for absolute calibration uncertainty. Larger values imply larger calibration uncertainty. The default value corresponds to roughly 2.5% uncertainty. This is essentially an alias of scaleFacPrior (see below) when a float is passed. 
    covarianceTruncation : float, optional, default 1e-3
        Threshold below which the covariance matrix is truncated to zero.
    scaleFacPrior : float, scipy.stats.rv_continuous, or None 
        Sets the prior on the (base-10 logarithm of the) scale factor which re-normalises the spectrum. If a float is passed, it is assumed to be the sigma for a LogNormal distribution (the same as calUnc). If a float is passed to both calUnc and scaleFacPrior, calUnc is treated as the sigma, and scaleFacPrior as mu (a bias term). Finally, you can pass an instance of rv_continuous, in which case that will be treated as the prior on the logarithm of the scale factor.
    scaleLengthPrior : float, scipy.stats.rv_continuous, or None
        Sets the prior on the scale length of the RBF kernel of the additional correlated noise. You can pass an instance of rv_continuous, in which case that will be treated as the prior on the length scale - beware, distributions which allow negative values will cause errors. If a float is passed, it is treated as the sigma (in micron) for a half-normal distribution used as the prior. If None (default) then a half-normal with sigma = 0.1 micron is used.
    covWeightPrior : float, scipy.stats.rv_continuous, or None
        Sets the prior on the scale factor for the additional correlated noise. As above, if an rv_continuous is passed it will be used directly as the prior on this scale factor - negative values will result in errors as the covariance matrix may not be positive semi-definite. If a float is passed, it will be used as the sigma for a half-normal distribution prior on this weighting factor. If None (default) then a half-normal with sigma = 0.2 is used. 
    label : str, default "Spectrum", 
    resampleMethod : callable or str, default "exact"
    
    
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
    create an instance by directly instantiating it
    
    >>> spec = Spectrum(wavelength_array, flux_array, uncertainty_array, "um", "Jy")

    or alternatively
    
    >>> spec = Spectrum.fromFile('filename', 'SPITZER-YAAAR')
    
    creates an instance from a Spitzer spectrum downloaded from the CASSIS database.
    '''

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

    def __init__(self, wavelength=None, value=None, uncertainty=None,
                 bandUnits=None, fluxUnits=None,
                 freqSpec = False, calUnc = None, covarianceTruncation = 1e-3,
                 scaleLengthPrior = None, scaleFacPrior = None, covWeightPrior = None,
                 label="Spectrum", resampleMethod="exact",
                 **kwargs):
        #self.filterName = filterName #Telescope/Instrument cf photometry
        self.type = 'Spectrum'
        self.label = label
        self.plotParams["label"] = label

        if value is None:
            raise ValueError("Value is a required input!")
        elif "Spectrum1D" in [t.__name__ for t in type(value).__mro__]:
            #this checks if it is (or inherets from) Spectrum1D without us having to have specutils as a dependence!
            #Having the input as a Spectrum1D object is nice, because they support all sorts of automatic transformations between different systems and units
            self.wavelength = value.wavelength.to(units.micron, equivalencies = units.spectral()).value
            self.value = value.flux.to(units.Jy,
                                      equivalencies = units.spectral_density(self.wavelength*units.micron)).value
            self.uncertainty = value.uncertainty.to(units.Jy,
                                      equivalencies = units.spectral_density(self.wavelength*units.micron)).value
        else:
            #It's not a Spectrum1D, so we have to unpack flux, uncertainty and wavelengths separately
            #wavelength must come first, since it is required to unpack the equivalencies for value and uncertainty
            if isinstance(wavelength, units.Quantity):
                self.wavelength = wavelength.to(units.micron, equivalencies = units.spectral()).value #this only works if things are defined in absolute units, not doppler (i.e. velocity). This should be accomodated somehow.
            elif isinstance(bandUnits, units.core.Unit) and np.ndim(wavelength) > 0:
                self.wavelength = (wavelength*bandUnits).to(units.micron, equivalencies = units.spectral()).value
            elif isinstance(bandUnits, str) and np.ndim(wavelength) > 0:
                bandUnits = units.Unit(bandUnits)
                self.wavelength = (wavelength*bandUnits).to(units.micron, equivalencies = units.spectral()).value
            else:
                raise TypeError("If wavelength is not a Quantity, bandUnits must be either a Unit or a unit string. \n However, wavelength is {0} and bandUnits is {1}".format(type(wavelength), type(bandUnits)))
            
            if isinstance(value, units.Quantity):
                self.value = value.to(units.Jy,
                                      equivalencies = units.spectral_density(self.wavelength*units.micron)).value
            elif isinstance(fluxUnits, units.core.Unit) and np.ndim(value) > 0:
                self.value = (value*fluxUnits).to(units.Jy,
                                                  equivalencies = units.spectral_density(self.wavelength*units.micron)).value
            elif isinstance(fluxUnits, str) and np.ndim(value) > 0:
                fluxUnits = units.Unit(fluxUnits)
                self.value = (value*fluxUnits).to(units.Jy,
                                                  equivalencies = units.spectral_density(self.wavelength*units.micron)).value
            else:
                raise TypeError("If value is not a Quantity, fluxUnits must be either a Unit or a unit string\n However, wavelength is {0} and bandUnits is {1}".format(type(value), type(fluxUnits)))

            if isinstance(uncertainty, units.Quantity):
                self.uncertainty = uncertainty.to(units.Jy,
                                                  equivalencies = units.spectral_density(self.wavelength*units.micron)).value
            elif isinstance(fluxUnits, units.core.Unit) and np.ndim(uncertainty) > 0:
                self.uncertainty = (np.array(uncertainty)*fluxUnits).to(units.Jy,
                                                              equivalencies = units.spectral_density(self.wavelength*units.micron)).value
            elif isinstance(fluxUnits, str) and np.ndim(uncertainty) > 0:
                fluxUnits = units.Unit(fluxUnits)
                self.uncertainty = (uncertainty*fluxUnits).to(units.Jy,
                                                              equivalencies = units.spectral_density(self.wavelength*units.micron)).value
            else:
                raise TypeError("If uncertainty is not a Quantity, fluxUnits must be either a Quantity or a unit string\n However, wavelength is {0} and bandUnits is {1}".format(type(uncertainty), type(fluxUnits)))

            
        """ Commenting this segment out, as above unit conversions should now take care of everything.
        #Wavelength conversion
        self.frequency = freqSpec #True/False? Wavelength or Frequency spectrum
        
        if freqSpec == 'True': #If freq, convert to wavelength in micron
            wavelength = (const.c / wavelength).to('um')
        
        self.bandUnits = bandUnits
        
        # as long as bandUnits has the dimensions of length, the wavelength
        # is computed in um.
        if units.m.is_equivalent(bandUnits):
            wavelength = wavelength / units.um.to(bandUnits)
        else:
            mesg = \"""The value for bandUnits must have dimensions of length!
            Allowed values: 'AA', 'um', 'nm', etc.""\"
            raise ValueError(mesg)
        
        self.wavelength = wavelength #Will be the grid of wavelength/frequency
                                     #using numpy.ma.asarray()? to make it a MaskedArray so that we can easily
                                     #mask out unwanted parts of the spectrum along with its covariances
        
        self.fluxUnits = fluxUnits#specFluxUnits

        if units.Jy.is_equivalent(fluxUnits):
            uconv = units.Jy.to(fluxUnits)
            value = value / uconv
            uncertainty = uncertainty / uconv
        elif (units.W / units.m**3).is_equivalent(fluxUnits):
            # first convert flux and uncertainty to W/m^3
            uconv = (units.W / units.m**3).to(fluxUnits)
            value = value / uconv
            uncertainty = uncertainty / uconv
            # now convert flux and uncertainty to Jy
            conv_fac = (wavelength * units.um)**2 / const.c
            value = (value * conv_fac).to('Jy')
            uncertainty = (uncertainty * conv_fac).to('Jy')
        else:
            raise NotImplementedError()

        self.value = value #Should always be a flux unless someone is peverse, as with wavelengths using MaskedArrays?
        self.uncertainty = uncertainty #Ditto #ma.asarray()?
        """

        ''' inititalise covariance matrix as a diagonal matrix '''
        self.varMat = uncertainty * uncertainty[:,np.newaxis] #make a square array of sigma_i * sigma_j, i.e. variances
        self.covMat = np.diag(np.ones_like(uncertainty)) #ma.asarray()
        #self.covMat.mask = np.logical_not(self.covMat > 0) #The mask for a MaskedArray is backwards compared to boolean indexing!
        #                                                   #By applying this mask we minimise the operations that need to be done, but this might really be overkill
        a = self.covMat > 0

        self.covMat = self.covMat * self.varMat # = np.diag(uncertainty**2)
        self.logDetCovMat = np.linalg.slogdet(self.covMat)[1]# / np.log(10.)
        #print(self.logDetCovMat)

        self.parLabels = ["calVar", "cov scale factor", "cov scale length"]


        #These bits need to be moved to @propertys to make getting and setting much cleaner!!
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
            self.scaleLengthPrior = halfnorm(loc=0., scale = 0.1) #The scale is in the same units as the wavelengths, so micron. For IR spectra, this is basically assuming that correlated noise (or residuals) on the scale of 0.1 micron is typical, but this would need to be chosen more carefully for spectra at shorter wavelenghts.
        
        #self.scaleLengthPrior = halfnorm(loc = 0, scale= 1.) #.logpdf(theta[2], 0., 1.)

        if isinstance(covWeightPrior, rv_continuous):
            self.covWeightPrior = covWeightPrior
        elif isinstance(covWeightPrior, float):
            self.covWeightPrior = halfnorm(loc=0., scale = covWeightPrior)
        else:
            self.covWeightPrior = halfnorm(loc=0., scale = 0.2) #choosing 0.2 as a default so that the uncertainty on the data dominates at the 5-sigma level

        self.covarianceTruncation = covarianceTruncation #if defined, covariance matrices will be truncated when their values are less than this number, to minimise the number of operation

        self.theta_last = [0.,1.] #Set some default theta values that correspond to uncorrelated noise so that selectWaves() will work without unwanted side effects

        self.selectWaves()

    def __call__(self, **kwargs):
        raise NotImplementedError()

    def setResampler(self, resampleMethod = "exact", **kwargs):
        '''Change the method used to resample the model spectra to the observed wavelengths

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
        
        '''

        if resampleMethod == "exact":
            self.resampler = spectres
            #l.append("Using exact resampling")
            print("Using exact resampling")
        elif resampleMethod=="fast":
            self.resampler = np.interp
            #l.append("Using fast resampling")
            print("Using fast resampling")
        elif callable(resampleMethod):
            self.resampler = resampleMethod
            #l.append("Using user-specified resampling")
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
        covMat = np.diag(np.ones_like(self.uncertainty[self.mask])) + theta[0]*m
        self.covMat = covMat * self.varMat

        #covMatmask = np.reshape(self.covMat[self.cov_mask], (self.value[mask].shape[0], self.value[mask].shape[0]))

        self.signDetCovMat, self.logDetCovMat = np.linalg.slogdet(self.covMat)
        if self.logDetCovMat == -np.inf: #This needs to be updated to raise an error! Or at least force the likelihood to be zero
            l.error("""The determinant of the covariance matrix for this dataset is 0.
            Please check that the uncertainties are positive real numbers """)
            l.error(self)
            l.error(theta)
            l.error(self.signDetCovMat,self.logDetCovMat)
            l.error(self.covMat)
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
        '''Transform from uniform RVs to the prior of any nuisance parameters. 

        This implementation has 3 nuisance parameters. These are 1) a scale factor, 
        which multiplies the observed spectrum to include the uncertainty on 
        absolute calibration, 2) the strength of the correlated component of the
        covariances and 3) the scale length (standard deviation) of the Gaussian
        kernel for the correlated component of the covariances. All three parameters
        are passed in with a single sequence.

        Parameters
        ----------
            u : array
            Nuisance parameters

        Returns
        -------
            theta : array
            Prior
        '''

        theta = np.zeros_like(u)
        theta[0] = 10**self.scaleFacPrior.ppf(u[0])
        theta[1] = self.covWeightPrior.ppf(u[1])
        theta[2] = self.scaleLengthPrior.ppf(u[2])
        return theta
    
    def lnprior(self, theta, **kwargs):
        '''Return the prior of any nuisance parameters. 

        This implementation has 3 nuisance parameters. These are 1) a scale factor, 
        which multiplies the observed spectrum to include the uncertainty on 
        absolute calibration, 2) the strength of the correlated component of the
        covariances and 3) the scale length (standard deviation) of the Gaussian
        kernel for the correlated component of the covariances. All three parameters
        are passed in with a single sequence.

        Parameters
        ----------
            theta : array
            Nuisance parameters

        Returns
        -------
            self.xxxPrior : array
            Prior of the nuisance parameters
        '''
        try:
            scaleFac = theta[0]
            #return self.scaleFacPrior.logpdf(np.log10(scaleFac)) + self.covWeightPrior.logpdf(theta[1]) + self.scaleLengthPrior.logpdf(theta[2])#norm.logpdf(np.log10(scaleFac), loc=0., scale = self.calUnc) + halfnorm.logpdf(theta[2], 0., 1.)
        except IndexError: #Only possible if theta is scalar or can't be indexed, in which case only the normalisation factor is in use
            scaleFac = theta
            return self.scaleFacPrior.logpdf(np.log10(scaleFac))
        #if scaleFac > 0 and 0. <= theta[1] <= 1. and theta[2] > 0.:
            #print(scalefac)
        return self.scaleFacPrior.logpdf(np.log10(scaleFac)) + self.covWeightPrior.logpdf(theta[1]) + self.scaleLengthPrior.logpdf(theta[2])#norm.logpdf(np.log10(scaleFac), loc=0., scale = self.calUnc) + halfnorm.logpdf(theta[2], 0., 1.)
        #return -np.inf

    def lnlike(self, theta, model, **kwargs):
        '''Compute the likelihood of the photometry given the model. 

        The likelihood is computed as
        
        .. math:: 
            \\frac{1}{2} N \\ln\\left(2\\pi\\right) - \\frac{1}{2}\\ln\\left(\\mathrm{det}C\\right) - \\frac{1}{2} \\left(\\alpha F_\\mathrm{obs} - F_\\mathrm{mod}\\right)^T C^{-1} \\left(\\alpha F_\\mathrm{obs} - F_\\mathrm{mod}\\right)
        
        where :math:`N` is the number of points in the spectrum, :math:`C` is the covariance matrix, :math:`\\alpha` is a free parameter to rescale the spectrum, and :math:`F_\\mathrm{obs}` and :math:`F_\\mathrm{mod}` are the observed and predicted spectrum, respectively.

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
            raise IndexError("theta is not scalar or can't be indexed.")
            scaleFac = theta

        #print(self)
        #wavelength = self.wavelength
        #modSpec = model.modelFlux #
        #print(model.wavelength)
        modSpec = self.resampler(self.wavelength, model.spectrum["wavelength"], model.spectrum["flux"])[self.mask] #For some reason spectres isn't cooperating :/ actually, looks like it was just a stupid mistake
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
            # l.append("NaN probability")
            # l.append(theta)
            # l.append(b)
            # l.append(inv(self.covMat))
            # l.append(self.value[self.mask])
            # l.append(modSpec)
            # l.append(model.modelFlux)
            # l.append(a)
            # l.append(np.matmul(inv(self.covMat), a))
            # l.append(np.matmul( a.T, np.matmul(inv(self.covMat), a) ) )

            print("NaN probability")
            print(theta)
            print(b)
            print(inv(self.covMat))
            print(self.value[self.mask])
            print(modSpec)
            print(model.spectrum["flux"])
            print(a)
            print(np.matmul(inv(self.covMat), a))
            print(np.matmul( a.T, np.matmul(inv(self.covMat), a) ) )
            return -np.inf #hack for now so we can see how often this occurs and hopefully troubleshoot it!
        #print(((np.float128(2.)*np.pi)**(len(self.value))), np.linalg.det(self.covMat))
        #print(((np.float128(2.)*np.pi)**(len(self.value)) * np.linalg.det(self.covMat)))
        #print(b, probFlux)
        #exit()
        return probFlux

    def simulate(self, theta, model, **kwargs):
        '''This function simulates a spectrum, given a model. The model spectrum is resampled to match the wavelengths of
        the data and then scaled by a factor theta[0]. Some noise is added to the simulated spectrum using a
        multivariate normal distribution with mean equal to the scaled model spectrum and covariance equal to the object's
        covariance matrix. The covariance matrix is updated using the remaining elements of theta.

        The simulated spectrum is returned as an array.

        Parameters
        ----------
            theta : array-like
                The scale factor for the model spectrum and the parameters for the covariance matrix.
            model : Model
                A model object containing a spectrum.

        Return
        ------
            simulated_data : array
                The simulated data of a model spectrum with noise added
        '''

        try:
            scaleFac = theta[0]
        except IndexError: #Only possible if theta is scalar or can't be indexed
            scaleFac = theta
        #First we resample the spectrum
        modSpec = self.resampler(self.wavelength, model.spectrum["wavelength"], model.spectrum["flux"])[self.mask]

        #now we rescale it with the scale factor
        #For the likelihood, we multiply the data by the scale factor
        #To maintain the same logic for the simulation system, we have to *divide* the simulated data by the scale factor
        modSpec = modSpec/scaleFac
        

        #Then, we add some noise - first we need an RNG
        rng = np.random.default_rng() #Optimise this, no need to re-create object each time, plus need seed to be stored so it can be saved
        #Now we update the covariance matrix
        self.cov(theta[1:])
        simulated_data = rng.multivariate_normal(modSpec, self.covMat, method='eigh') #We draw one sample from the multivariate normal distribution with mean = model spectrum and covariance = data covariances
        
        #now we can return the simulated data
        return simulated_data

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
        try:
            print("Spectrum has ",np.unique(table['chunk'].data)," chunks, breaking up the spectrum into separate chunks.")
            for chunk in np.unique(table['chunk'].data):
                selection = table['chunk'] == chunk
                self=cls.__new__(Spectrum)
                self.__init__(table['wavelength'][selection].data, value[selection], uncertainty[selection], bandUnits, photUnits, **kwargs) #Also pass in flux units
                specList.append(self)
            #l.append("Table['chunk'] column found : spectrum has",np.unique(table['chunk'].data," chunks.")) #I don't know where this line comes from or what it is for - will have to investigate...
            #Looks like it was a botched attempt to add logging here - needs to be fixed properly but that requires more extensive updates to the ampere logger!
        except:
            print("Spectrum has no chunks, assuming one continuous order.")
            self=cls.__new__(Spectrum)
            self.__init__(table['wavelength'][selection].data, value[selection], uncertainty[selection], bandUnits, photUnits, **kwargs) #Also pass in flux units
            specList.append(self)
            #l.append("No table['chunk'] column found : spectrum assumed to be a single chunk.")

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

        '''
        Plot the spectrum. Figure is either presented on screen or saved to a file.
    
        Parameters
        ----------
        fig : matplotlib.figure.Figure, optional
            Figure to plot on. If not provided, a new figure will be created.
        ax : matplotlib.axes.Axes, optional
            Axes object. If not provided, new axes will be created on the provided figure or a new figure.
        unmask : bool, optional
            If True, plot all wavelengths, ignoring any previously defined mask. Default is False.
        doPlot : bool, optional
            If False, do not perform the plot. Default is True.
        savePlot : bool, optional
            If True, save the plot to a file. Default is False.
        showPlot : bool, optional
            If True, show the plot using `matplotlib.pyplot.show`. Default is False.
        **kwargs : optional
            Additional arguments to pass to `matplotlib.pyplot.plot` or `matplotlib.pyplot.fill_between`.
        '''
        
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
        """
        Plot a realisation of the data model.

        Parameters
        ----------
        s : array
            The realisation of the data model to plot. The first element of `s` should be the scaling factor applied to the model.
        ax : matplotlib.axes.Axes, optional
            The Axes to plot on. If not provided, a new Axes will be created.
        **kwargs
            Additional keyword arguments to be passed to the plot function.
        """

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
