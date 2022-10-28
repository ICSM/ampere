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

    def simulate(self, results, **kwargs):
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
