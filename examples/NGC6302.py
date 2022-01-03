import numpy as np
import os
import ampere
from ampere.data import Spectrum, Photometry
from ampere.emceesearch import EmceeSearch
from ampere.models import Model
from spectres import spectres
import pyphot
from emcee import moves

# First we will define a rather simple model


class OpacitySpectrum(Model):
    '''This model fits a modfied blackbody multiplied by sum of opacties,
        consisting of a warm and cold component, not necessarily of the same
        composition, over two temperature ranges, following Kemper et al. 2002.
        We use the following parameters:
       wavelengths : the array of wavelengths considered
       opacityFileList
       acold : (relative) abundances for the cold component
       awarm : (relative) abundances for the warm component
       Tcold : temperature range of the cold component
       Twarm : temperature range of the warm component
       indexp : index p of the density distribution
       indexq : index q of the temperature distribution
       multfact : multiplication factor that the sum of the modified black
                  bodies has to be multiplied with to fit the spectrum
       Output is an array of model fluxes (fluxes), to match wavelengths
    '''
    def __init__(self, wavelengths, flatprior=True,
                 opacityFileList=None):
        '''The model constructor, which will set everything up

        This method does essential setup actions, primarily things that
        may change from one fit to another, but will stay constant throughout
        the fit. This may be things like the grid of wavelengths to calculate
        the model output on, or establishing the dust opacities if involved.
        There are also several important variables it *MUST* define here
        '''
        self.wavelength = wavelengths
        import os
        from scipy import interpolate
        opacityDirectory = os.path.dirname(__file__)+'/Opacities/'
        opacityFileList = np.array(opacityFileList)
        nSpecies = opacityFileList.__len__()
        opacity_array = np.zeros((wavelengths.__len__(), nSpecies))
        for j in range(nSpecies):
            tempData = np.loadtxt(opacityDirectory + opacityFileList[j],
                                  comments='#')
            tempWl = tempData[:, 0]
            tempOpac = tempData[:, 1]
            f = interpolate.interp1d(tempWl, tempOpac, assume_sorted=False)
            opacity_array[:, j] = f(self.wavelength)
        self.opacity_array = opacity_array
        self.nSpecies = nSpecies
        self.npars = nSpecies + 7
        # 7 is the number of free parameters for the model (__call__()). For
        # some models this can be determined through introspection, but it is
        # still strongly recommended to define this explicitly here.
        # Introspection will only be attempted if self.npars is not defined.
        self.npars_ptform = 2
        # Sometimes the number of free parameters is different when using the
        # prior transform instead of the prior. In that case, self.npars_ptform
        # should also be defined.
        # You can do any other set up you need in this method.
        # For example, we could define some cases to set up different priors
        # But that's for a slightly more complex example.
        # Here we'll just use a simple flat prior
        self.lims = lims  # This one is still to be defined based on npars
        self.flatprior = flatprior

    def __call__(self, *args, **kwargs):
        '''The model itself, using the callable class functionality of python.

        This is an essential method. It should do any steps required to
        calculate the output fluxes. Once done, it should stop the output fluxes
        in self.modelFlux.
        '''
        # print out *args
        print(args)
        # print out **kwargs
        for key, value in kwargs.items():
            print(key + " : " + value)

        # number of dust species in opacityDirectory are all going to be fitted
        # for each we will consider a hot and a cold component. We will allow
        # only 2 temperature ranges, the same for all dust species.
        # the module translated from ck_modbb calculates the flux for each
        # component at a single temperature
        # ckmodbb and shbb are now written. Now I just need to figure out how
        # to do the fit over 8 dust species and 2 temperature ranges.
        # let's first do this over fixed temperature ranges, instead of
        # allowing them to be free the temperature ranges are 118-100 K and
        # 60-30 K, following Kemper et al. 2002

        coldcomponent = np.zeros((2, wavelengths.__len__()))
        coldcomponent[0, :] = self.wavelengths
        warmcomponent = np.zeros((2, wavelengths.__len__()))
        warmcomponent[0, :] = self.wavelengths
        for i, ac in enumerate(abundances[0,:]):
            onespeciescold = ckmodbb(opacity_array[:,i], tin = Tcold[0], tout = Tcold[1], n0 = abundances[0, i], index = indexp)
            coldcomponent[1, :] = coldcomponent[1,:] + onespeciescold
        for i, aw in enumerate(abundances[0,:]):
            onespecieswarm = ckmodbb(opacity_array[:,i], tin = Twarm[0], tout = Twarm[1], n0 = abundances[0,i], index = indexp)
            warmcomponent[1,:] = warmcomponent[1,:] + onespecieswarm
        fModel = np.like(coldcomponent)
        fModel[1,:] = fModel[1,:] + warmcomponent[1,:]
        self.modelFlux = fModel[1,:]

# hiero, I am getting all my indices of fnu etc. wrong. Check the functions
# below, especially whether I am using wavelength and flux in the correct place    
    def ckmodbb(self, q, tin, tout, n0, index = 0.5, r0 = 1e15, distance = 910., grainsize = 0.1, steps = 10):
        d = distance * 3.0857e18 #convert distance from pc to cm
        a = grainsize * 1e-4 #convert grainsize from micron to cm
 
        fnu = np.zeros((2, wavelengths.__len__()))
        fnu[0,:] = self.wavelengths

        for i in range(steps - 1):
            t = tin - i * (tin-tout)/steps
            power = (t/tin)^(2*index - 6)
            bb = shbb(fnu, t, 0.) 
            fnu[1,:] = np.add(fnu[1,:],np.multiply(q[1,:],bb[1,:])*(power*((tin-tout)/steps)))

        extra = r0/d
        factor = 4 * math.pi * a * a * r0 * n0 * extra * extra / (3-index)
        fnu[1,:] = fnu[1,:] * factor
        return fnu

    def shbb(self, aar, temp, pinda):
        a1 = 3.97296e19
        a2 = 1.43875e4
        mbb = np.copy(aar)                                 
        bbflux = a1/(mbb[0,:]^3)/(exp(a2/(mbb[0,:]*temp))-1)
        mbb[1,:] = bbflux * mbb[0,:]^pinda
        return mbb                                 

    def lnprior(self, theta, **kwargs): #hiero: too be done still
        slope = theta[0]
        intercept = theta[1]
        if self.flatprior:
            if (self.lims[0,0] < theta[0] < self.lims[0,1]) and \
               (self.lims[1,0] < theta[1] < self.lims[1,1]) and \
               (self.lims[2,0] < theta[2] < self.lims[2,1]) and \
               (self.lims[3,0] < theta[3] < self.lims[3,1]) and \
                np.sum(10**theta[4:]) <= 1. and np.all(theta[4:] < 0.):
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

if __name__ == "__main__": #hiero: too be done still, check all below
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
    filts, modSed = pyphot.extractPhotometry(wavelengths,
                                             model_flux,
                                             filters,
                                             Fnu = True,
                                             absFlux = False,
                                             progress=False
            )

    input_noise_phot = 0.1 #Fractional uncertainty
    photunc = input_noise_phot * modSed #Absolute uncertainty
    modSed = modSed + np.random.randn(len(filterName)) * photunc #Now perturb data by drawing from a Gaussian distribution

    
    #now we'll create a synthetic spectrum from the model fluxes, using a Spitzer IRS observation to get the wavelength sampling
    dataDir = os.getcwd() + '/PGQuasars/PG1011-040/'
    specFileExample = 'cassis_yaaar_spcfw_14191360t.fits'
    irsEx = Spectrum.fromFile(dataDir+specFileExample,format='SPITZER-YAAAR')
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
    photometry = Photometry(filterName=filterName, value=modSed, uncertainty=photunc, photUnits='Jy', libName=libname)
    #print(photometry.filterMask)
    photometry.reloadFilters(wavelengths)

    dataset = [photometry,
               #spec0, #Fitting spectra is slow because it needs to do a lot of resampling
               spec1   #As a result, we're leaving some of them out
               ]


    #Ampere exposes acces to emcee's moves interface. This can be useful if the posterior turns out to not be well behaved - the default move only deals well with posteriors that are monomodal and approximately Gaussian. Here's an example that usually deals a bit better with posteriors that don't meet these criteria:
    m = [(moves.DEMove(), 0.8),
        (moves.DESnookerMove(), 0.2),
         ]

    #Now we set up the optimizer object:
    optimizer = EmceeSearch(model=model, data=dataset, nwalkers=100, moves=m)
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
    optimizer.optimise(nsamples = 1500, burnin=1000, guess=guess
                       )


    optimizer.postProcess() #now we call the postprocessing to produce some figures
