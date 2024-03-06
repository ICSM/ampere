# %% [markdown]
# # Truncated Marginal Neural Ratio Estimation with _Swyft_
# This notebook demonstrates how to use _Swyft_ to perform Truncated Marginal Neural
# Ratio Estimation (TMNRE) in Ampere. We will use the same simple example to demonstrate
# the basic steps. The example is a simple linear model of flux as a function of
# wavelength.

# %% [markdown]
# ## Imports
#
# Before we do anything else, we need to import a few key packages.

# %%
import numpy as np
import os
import ampere
from ampere.data import Spectrum, Photometry
from ampere.infer.sbi import Swyft_TMNRE
from ampere.models import Model
from spectres import spectres
import pyphot
from astropy import units
from pathlib import Path

# %% [markdown]
# ##Define the model
#
# Now we can define our model. The model is a simple linear model of flux as a function
# of wavelength. We will use this to generate some synthetic data in a minute.

# %%

class ASimpleModel(Model):
    '''This is a very simple model in which the flux is linear in wavelength.

    This model shows you the basics of writing a model for ampere

    '''
    def __init__(self, wavelengths, flatprior=True, lims=None):
        '''The model constructor, which will set everything up

        This method does essential setup actions, primarily things that 
        may change from one fit to another, but will stay constant throughout 
        the fit. This may be things like the grid of wavelengths to calculate
        the model output on, or establishing the dust opacities if involved.
        There are also several important variables it *MUST* define here
        '''
        if lims is None:
            lims = np.array([[-10, 10], [-10, 10]])
        self.wavelength = wavelengths
        self.npars = 2 #Number of free parameters for the model (__call__()). For some 
        # models this can be determined through introspection, but it is still strongly 
        # recommended to define this explicitly here. Introspection will only be 
        # attempted if self.npars is not defined.
        self.npars_ptform = 2 #Sometimes the number of free parameters is different 
        # when using the prior transform instead of the prior. In that case, 
        # self.npars_ptform should also be defined.

        # You can do any other set up you need in this method.
        # For example, we could define some cases to set up different priors
        # But that's for a slightly more complex example.
        # Here we'll just use a simple flat prior
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
        return {"spectrum": {"wavelength": self.wavelength, "flux": self.modelFlux}}

    def lnprior(self, theta, **kwargs):
        """The model prior probability distribution

        The prior is essential to most types of inference with ampere. The
        prior describes the relative weights (or probabilities if normalised)
        of different parameter combinations. Using a normalised prior with
        SNPE is strongly recommended, otherwise ampere will attempt
        approximate normalisation using Monte Carlo integration.
        """
        if not self.flatprior:
            raise NotImplementedError()
        slope = theta[0]
        #print(slope)
        intercept = theta[1]
        return (
            0
            if self.lims[0, 0] < slope < self.lims[0, 1]
            and self.lims[1, 0] < intercept < self.lims[1, 1]
            else -np.inf
        )

    def prior_transform(self, u, **kwargs):
        '''The prior transform, which takes samples from the Uniform(0,1)
        distribution to the desired distribution.

        Prior transforms are essential for SNPE. SNPE needs to be able to
        generate samples from the prior, and this method is integral to doing 
        so. Therefore, unlike other inference methods, if you want to use 
        SNPE (or other SBI approaches) you need to define *both* lnprior and 
        prior_transform.
        '''
        if self.flatprior:
            return (self.lims[:,1] - self.lims[:,0]) * u + self.lims[:,0]
        else:
            raise NotImplementedError()

# %% [markdown]
# ## Set up the data
#
# Now we need to set up the data. We will use the model to generate some synthetic
# data with intercept and slope set to 1.0. We will then use this synthetic data to
# perform inference.
        
# %%
""" Set up the inputs for the model """
""" wavelength grid """
wavelengths = 10**np.linspace(0.,1.9, 2000)

""" Choose some model parameters """
slope = 1. #Keep it super simple for now
intercept = 1.

# Now we can initialise the model
model = ASimpleModel(wavelengths)
# And call it to produce the fluxes for our chosen parameters
model(slope, intercept)
model_flux = model.modelFlux

# Now we create synthetic data:
# First photometry
# This is minimal, so we'll just have two bands well separated in 
# wavelength
filterName = np.array(['WISE_RSR_W1', 'SPITZER_MIPS_70'])  

libDir = Path(ampere.__file__.strip('__init__.py')) # '/home/peter/pythonlibs/ampere/ampere/'
print(f'{libDir=}')
libname = libDir/'ampere_allfilters.hd5'
print(f'{libname=}')
filterLibrary = pyphot.get_library(fname=libname)
filters = filterLibrary.load_filters(filterName, interp=True,
                                     lamb=wavelengths*pyphot.unit['micron'])
# Now we need to extract the photometry with pyphot
# first we need to convert the flux from Fnu to Flambda
flam = model_flux / wavelengths**2
modSed = []
for i, f in enumerate(filters):
    lp = f.lpivot.to("micron").value
    fphot = f.get_flux(wavelengths*pyphot.unit['micron'], flam*pyphot.unit['flam'], axis=-1).value
    print(fphot)
    modSed.append(fphot*lp**2)
print(modSed)
modSed = np.array(modSed)
print(modSed)
print(units.Quantity(modSed, 'Jy'))
print(modSed*units.Jy)
print(modSed*units.Unit('Jy'))

#exit()

input_noise_phot = 0.1  # Fractional uncertainty
photunc = input_noise_phot * modSed  # Absolute uncertainty
# Now perturb data by drawing from a Gaussian distribution
modSed = modSed + np.random.randn(len(filterName)) * photunc 

# now we'll create a synthetic spectrum from the model fluxes, using a Spitzer IRS 
# observation to get the wavelength sampling
#print(ampere.__file__.split('/'))
dataDir = Path(__file__).parent / 'test_data'
#dataDir = f"{'/'.join(ampere.__file__.split('/')[:-2])}/examples/test_data/"
print(f'{dataDir=}')
specFileExample = dataDir/'cassis_yaaar_spcfw_14191360t.fits'
irsEx = Spectrum.fromFile(specFileExample,  #os.path.normpath(dataDir+specFileExample),
                          format='SPITZER-YAAAR')
spec0 = spectres(irsEx[0].wavelength,wavelengths,model_flux)
spec1 = spectres(irsEx[1].wavelength,wavelengths,model_flux)

# And again, add some noise to it
input_noise_spec = 0.1
unc0 = input_noise_spec*spec0
unc1 = input_noise_spec*spec1
spec0 = spec0 + np.random.randn(len(spec0))*unc0
spec1 = spec1 + np.random.randn(len(spec1))*unc1

spec0 = Spectrum(irsEx[0].wavelength, spec0, unc0, "um", "Jy", 
                    calUnc=0.0025, 
                    scaleLengthPrior=0.01)  # , resampleMethod=resmethod)
spec1 = Spectrum(irsEx[1].wavelength, spec1, unc1, "um", "Jy", 
                    calUnc=0.0025, 
                    scaleLengthPrior=0.01)  # , resampleMethod=resmethod)

# Now let's try changing the resampling method so it's faster
# This model is very simple so exact flux conservation is not important
resmethod = "fast"  # "exact"#"fast"#
spec0.setResampler(resampleMethod=resmethod)
spec1.setResampler(resampleMethod=resmethod)

""" now set up ampere to try and fit the same stuff """
photometry = Photometry(filterName=filterName, value=modSed, 
                        uncertainty=photunc, photunits="Jy", 
                        libName=libname)
# print(photometry.filterMask)
photometry.reloadFilters(wavelengths)

dataset = [photometry,
            # spec0, #Fitting spectra is slow because it needs to do a lot of resampling
            spec1   # As a result, we're leaving some of them out
            ]

# %% [markdown]
# ## Set up the optimizer
#
# Now we have a model and some data, we can set up the optimizer. We will use
# Swyft_TMNRE, which wraps the official implementation of Truncated Marginal
# Neural Ratio Estimation (TMNRE) in _Swyft_. This is a powerful method for
# likelihood-free inference. It is particularly useful when the likelihood is
# expensive to evaluate, and when the data or model is high-dimensional.

# %%
# Now we set up the optimizer object:
optimizer = Swyft_TMNRE(model=model, data=dataset)

# Now we can do the inference
optimizer.optimise()

# %% [markdown]
# ## Post-processing
#
# Now we can do some post-processing. This will produce some figures and save the
# optimizer object to a file so we can read it back in later.

# %%
optimizer.postProcess() #now we call the postprocessing to produce some figures
# Save the optimizer object to a file so we can read it back in later
import pickle
with open("test_pickle_swyft.pkl", 'wb') as f:
    pickle.dump(optimizer, f)

# Now try reading it back in:
with open("test_pickle_swyft.pkl", 'rb') as f2:
    opt2 = pickle.load(f2)

# %% [markdown]
# ## Conclusion
#
# This notebook demonstrates how to use _Swyft_ to perform Truncated Marginal Neural
# Ratio Estimation (TMNRE) in Ampere. We used a simple example to demonstrate the
# basic steps. The example is a simple linear model of flux as a function of wavelength.
# We generated some synthetic data and used it to perform inference. We then did some
# post-processing, producing some figures and saving the optimizer object to a file so
# we can read it back in later.
