import numpy as np
import ampere
from ampere.models import Model
from astropy.modeling.physical_models import BlackBody
from astropy import units as u
#from scipy.stats import norm, halfnorm, truncnorm
from scipy.stats import uniform


class ModifiedBlackBody(Model):
    def __init__(self, wavelength, 
                 kappa=10, kappawave=250,
                 parallax = 10, sigmaparallax = 1,
                 lims = np.array([[10, 50],
                                  [-3,3],
                                  [-3,0],
                                  [0.05,0.15]]),
                 **kwargs):
        self.freq = (wavelength * u.micron).to(u.Hz, equivalencies=u.spectral()).value
        self.wavelength = wavelength
        self.kappa = kappa
        self.kappawave = kappawave
        self.npars = 4
        self.npars_ptform = 4
        self.parallax = parallax
        self.sigmaparallax = sigmaparallax
        self.parallaxsnr = parallax/sigmaparallax
        self.parLabels = ['temperature', 'log10 mass', 'beta', 'distance']
        self.priors = [uniform(lim[0], lim[1]-lim[0]) for lim in lims]

    def __call__(self, t, logm, beta, d, **kwargs):
        modelFlux = BlackBody().evaluate(self.freq, t, 1*u.Jy/u.sr)
        modelFlux = modelFlux/(d*u.pc.to(u.cm))**2
        modelFlux = modelFlux * 10**(logm) * (self.kappa) * (self.wavelength/self.kappawave)**beta
        self.modelFlux = modelFlux
        return {"spectrum":{"wavelength":self.wavelength, "flux": modelFlux}}

    def lnprior(self, theta, **kwargs):
        #t_lnprior = halfnorm.logpdf(theta[0], loc=0, scale = 100) #Very weakly-informative prior suggesting we're looking for cold dust
        #logm_lnprior = norm.logpdf(theta[1], loc = 0, scale = 3)  #weakly informative prior on log(mass)
        #beta_lnprior = norm.logpdf(theta[2], loc = 2, scale = 0.5) #We expect to get beta values close to 2
        #dist_lnprior = truncnorm.logpdf(((1/theta[3]) - self.parallax)/self.sigmaparallax, -self.parallaxsnr, np.inf) #The 
        return np.sum([self.priors[i].logpdf(theta[i]) for i in range(len(theta))])#t_lnprior+logm_lnprior+beta_lnprior+dist_lnprior

    def prior_transform(self, u, **kwargs):
        theta = np.zeros_like(u)
        #theta[0] = halfnorm.ppf(u[0], loc=0, scale=100)
        #theta[1] = norm.ppf(u[1], loc=0, scale=3)
        #theta[2] = norm.ppf(u[2], loc=-2, scale=0.5)
        #theta[3] = 1/truncnorm.ppf(u[3], -self.parallaxsnr, np.inf) * self.sigmaparallax + self.parallax #transform from a uniform RV to an RV sampled from a truncated normal described by the parallax (postive values only) and then inversed
        theta = np.array([
                            self.priors[i].ppf(u[i]) for i in range(len(u))
                         ])
        return theta



if __name__=="__main__":
    wavelengths = 10**np.linspace(0.,2.9, 2000)
    t_true = 30.
    logm_true = 1
    beta_true = -2
    d_true = 0.11 #kpc
    
    model = ModifiedBlackBody(wavelengths)
    
    res = model(t_true, logm_true, beta_true, d_true)
    f_true = res['spectrum']['flux']
    
    import pyphot
    #Now we create synthetic photometry
    filterName = np.array(['WISE_RSR_W4', 'SPITZER_MIPS_70']) #This is minimal, so we'll just have two bands well separated
    filterName = np.array(['AKARI_FIS_N160', 'AKARI_FIS_N60', 'AKARI_FIS_WIDEL', 'AKARI_FIS_WIDES', 'HERSCHEL_PACS_100', 'HERSCHEL_PACS_160', 'HERSCHEL_PACS_70', 'HERSCHEL_SPIRE_250', 'HERSCHEL_SPIRE_350', 'HERSCHEL_SPIRE_500'])

    libDir = ampere.__file__.strip('__init__.py') # '/home/peter/pythonlibs/ampere/ampere/'
    libname = libDir + 'ampere_allfilters.hd5'
    filterLibrary = pyphot.get_library(fname=libname)
    filters = filterLibrary.load_filters(filterName, interp=True, lamb = wavelengths*pyphot.unit['micron'])
    filts, modSed = pyphot.extractPhotometry(wavelengths,
                                             f_true,
                                             filters,
                                             Fnu = True,
                                             absFlux = False,
                                             progress=False
                                             )

    input_noise_phot = 0.1 #Fractional uncertainty
    photunc = input_noise_phot * modSed #Absolute uncertainty
    modSed = modSed + np.random.randn(len(filterName)) * photunc #Now perturb data by drawing from a Gaussian distribution

    from ampere.data import Photometry
    photometry = Photometry(filterName=filterName, value=modSed, uncertainty=photunc, photUnits='Jy', libName=libname)
    #print(photometry.filterMask)
    photometry.reloadFilters(wavelengths)
    print(photometry)
    dataset = [photometry]

    from ampere.infer.emceesearch import EmceeSearch
    from emcee import moves

    m = [(moves.DEMove(), 0.8),
        (moves.DESnookerMove(), 0.2),
         ]

    optimizer = EmceeSearch(model=model, data=dataset, nwalkers=100, moves=m, vectorize = False)
    optimizer.optimise(nsamples = 150, burnin=100, guess='None'
                       )
    optimizer.postProcess()


    #from ampere.infer.zeussearch import ZeusSearch
    #from zeus import moves
    #m = [(moves.GlobalMove(), 0.8),
    #    (moves.DifferentialMove(), 0.2),
    #     ]
    
    #optimizer = ZeusSearch(model=model, data=dataset, nwalkers=100, vectorize = False)
    #optimizer.optimise(nsamples = 150, burnin=100, guess=guess
    #                   )
    #optimizer.postProcess()

    #from ampere.infer.dynestysearch import DynestyNestedSampler
    #optimizer = DynestyNestedSampler(model=model, data=dataset)
    #optimizer.optimise(dlogz = 5.)
    #optimizer.postProcess()


    from ampere.infer.sbi import SBI_SNPE
    optimizer = SBI_SNPE(model=model, data=dataset)
    optimizer.optimise(nsamples = 10000, nsamples_post = 10000
                       )

    optimizer.postProcess()
