import numpy as np
import os
import ampere
from ampere.data import Spectrum, Photometry
from ampere.infer.zeussearch import ZeusSearch
from ampere.models import Model
from spectres import spectres
import pyphot
from zeus import moves

#First we will define a rather simple model

class ASimpleModel(Model):
    '''This is a very simple model in which the flux is linear in wavelength.

    This model shows you the basics of writing a model for ampere

    '''
    def __init__(self, wavelenghts, flatprior=True,
                 lims=np.array([[-10, 10],
                                [-10, 10]])):
        '''The model constructor, which will set everything up

        This method does essential setup actions, primarily things that 
        may change from one fit to another, but will stay constant throughout 
        the fit. This may be things like the grid of wavelengths to calculate
        the model output on, or establishing the dust opacities if involved.
        There are also several important variables it *MUST* define here
        '''
        self.wavelength = wavelengths
        self.npars = 2 #Number of free parameters for the model (__call__()). For some models this can be determined through introspection, but it is still strongly recommended to define this explicitly here. Introspection will only be attempted if self.npars is not defined.
        self.npars_ptform = 2 #Sometimes the number of free parameters is different when using the prior transform instead of the prior. In that case, self.npars_ptform should also be defined.
        #You can do any other set up you need in this method.
        #For example, we could define some cases to set up different priors
        #But that's for a slightly more complex example.
        #Here we'll just use a simple flat prior
        self.lims = lims
        self.flatprior = flatprior
        self.parLabels = ["slope", "intercept"]

    def __call__(self, slope, intercept, **kwargs):
        '''The model itself, using the callable class functionality of python.

        This is an essential method. It should do any steps required to 
        calculate the output fluxes. Once done, it should stop the output fluxes
        in self.modelFlux.
        '''
        self.modelFlux = slope*self.wavelength + intercept
        return {"spectrum":{"wavelength":self.wavelength, "flux": self.modelFlux}}

    def lnprior(self, theta, **kwargs):
        slope = theta[0]
        intercept = theta[1]
        if self.flatprior:
            if (self.lims[0,0] < theta[0] < self.lims[0,1]) and (self.lims[1,0] < theta[1] < self.lims[1,1]):
                return 0
            else:
                return -np.inf
        else:
            raise NotImplementedError()
        

    def prior_transform(self, u, **kwargs):
        '''The prior transform, which takes samples from the Uniform(0,1) distribution to the desired distribution.

        This is only included for completeness and to demonstrate how a prior 
        transform function should look. This example only uses emcee for 
        fitting, which uses the lnprior function instead. Prior transforms are 
        required by nested-sampling codes and similar approaches.
        '''
        if self.flatprior:
            theta = np.zeros_like(u)
            return (self.lims[:,1] - self.lims[:,0]) * u + self.lims[:,0]
        else:
            raise NotImplementedError()

if __name__ == "__main__":
    """ Set up the inputs for the model """
    """ wavelength grid """
    wavelengths = 10**np.linspace(0.,1.9, 2000)

    """ Choose some model parameters """
    slope = 1. #Keep it super simple for now
    intercept = 1.

    #Now init the model:
    model = ASimpleModel(wavelengths)
    #And call it to produce the fluxes for our chosen parameters
    model(slope, intercept)
    model_flux = model.modelFlux

    #Now we create synthetic data:
    #First photometry
    filterName = np.array(['WISE_RSR_W1', 'SPITZER_MIPS_70']) #This is minimal, so we'll just have two bands well separated

    libDir = ampere.__file__.strip('__init__.py') # '/home/peter/pythonlibs/ampere/ampere/'
    libname = libDir + 'ampere_allfilters.hd5'
    filterLibrary = pyphot.get_library(fname=libname)
    filters = filterLibrary.load_filters(filterName, interp=True, lamb = wavelengths*pyphot.unit['micron'])
    #Now we need to extract the photometry with pyphot
    #first we need to convert the flux from Fnu to Flambda
    flam = model_flux / wavelengths**2
    modSed = []
    for i, f in enumerate(filters):
        lp = f.lpivot.to("micron").value
        fphot = f.get_flux(wavelengths*pyphot.unit['micron'], flam*pyphot.unit['flam'], axis=-1).value
        print(fphot)
        modSed.append(fphot*lp**2)

    modSed = np.array(modSed)

    input_noise_phot = 0.1 #Fractional uncertainty
    photunc = input_noise_phot * modSed #Absolute uncertainty
    modSed = modSed + np.random.randn(len(filterName)) * photunc #Now perturb data by drawing from a Gaussian distribution

    
    #now we'll create a synthetic spectrum from the model fluxes, using a Spitzer IRS observation to get the wavelength sampling
    dataDir = f'{os.getcwd()}/test_data/'
    specFileExample = 'cassis_yaaar_spcfw_14191360t.fits'
    irsEx = Spectrum.fromFile(os.path.normpath(dataDir+specFileExample),format='SPITZER-YAAAR')
    spec0 = spectres(irsEx[0].wavelength,wavelengths,model_flux)
    spec1 = spectres(irsEx[1].wavelength,wavelengths,model_flux)

    #And again, add some noise to it
    input_noise_spec = 0.1
    unc0 = input_noise_spec*spec0
    unc1 = input_noise_spec*spec1
    spec0 = spec0 + np.random.randn(len(spec0))*unc0
    spec1 = spec1 + np.random.randn(len(spec1))*unc1
    
    spec0 = Spectrum(irsEx[0].wavelength, spec0, unc0,"um", "Jy",calUnc=0.0025, scaleLengthPrior = 0.01) #, resampleMethod=resmethod)
    spec1 = Spectrum(irsEx[1].wavelength, spec1, unc1,"um", "Jy",calUnc=0.0025, scaleLengthPrior = 0.01) #, resampleMethod=resmethod)

    #Now let's try changing the resampling method so it's faster
    #This model is very simple so exact flux conservation is not important
    resmethod = "fast" #"exact"#"fast"#
    spec0.setResampler(resampleMethod=resmethod)
    spec1.setResampler(resampleMethod=resmethod)

    """ now set up ampere to try and fit the same stuff """
    photometry = Photometry(filterName=filterName, value=modSed, uncertainty=photunc, photunits='Jy', libName=libname)
    #print(photometry.filterMask)
    photometry.reloadFilters(wavelengths)

    dataset = [photometry,
               #spec0, #Fitting spectra is slow because it needs to do a lot of resampling
               spec1   #As a result, we're leaving some of them out
               ]


    #Ampere exposes acces to Zeus' moves interface. This can be useful if the posterior turns out to not be well behaved. This shows you how to define moves, but doesn't actually use them in the inference
    m = [(moves.GlobalMove(), 0.8),
        (moves.DifferentialMove(), 0.2),
         ]
    #For example, You might want to burn in your sampler, then re-build it with new moves if the posterior is very multimodal

    #Now we set up the optimizer object:
    optimizer = ZeusSearch(model=model, data=dataset, nwalkers=100, vectorize = False)
    guess = [
        [1, 1, #The parameters of the model
         #1.0, 0.1, 0.1, #Each Spectrum object contains a noise model with three free parameters
         #The first one is a calibration factor which the observed spectrum will be multiplied by
         #The second is the fraction of correlated noise assumed
         #And the third is the scale length (in microns) of the correlated component of the noise
         1.0 ,0.1, 0.1
       ] #
        + np.random.rand(optimizer.npars)*[1,1,
                                           #1,1,1,
                                           1,1,1
                                           ]
        for i in range(optimizer.nwalkers)]

    #guess = "None"

    #Then we tell it to explore the parameter space
    optimizer.optimise(nsamples = 150, burnin=100, guess=guess
                       )


    optimizer.postProcess() #now we call the postprocessing to produce some figures
    # Save the optimizer object to a file so we can read it back in later
    import pickle
    with open("test_pickle_zeus.pkl", 'wb') as f:
        pickle.dump(optimizer, f)

    #Now try reading it back in:
    with open("test_pickle_zeus.pkl", 'rb') as f2:
        opt2 = pickle.load(f2)
