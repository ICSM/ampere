import numpy as np
import os
from spectres import spectres
from astropy import units, constants
from astropy.io import fits, ascii
from scipy.interpolate import interp1d
from scipy.stats import uniform
import Starfish
from Starfish.emulator import Emulator
from Starfish.grid_tools import GridInterface
from Starfish.grid_tools import PHOENIXGridInterfaceNoAlpha as PHOENIX #NoAlpha
#from Starfish.grid_tools import BTSettlGridInterface as BTSettl
import Starfish.constants as C
from Starfish.utils import create_log_lam_grid
from Starfish.grid_tools.utils import vacuum_to_air, idl_float
from Starfish.grid_tools import download_PHOENIX_models as dpm
from Starfish.grid_tools import HDF5Creator
import ampere
from ampere.models import Model
import matplotlib.pyplot as plt
from extinction import apply, ccm89

#calculate the bolometric flux of a one solar luminosity star at a distance of 1 pc:
fbol_1l1p = 1*units.Lsun
fbol_1l1p = ((fbol_1l1p/(4*np.pi*(1*units.pc)**2)).to(units.erg/units.s/units.cm/units.cm)).value
c = constants.c.to(units.AA/units.s)

class StarfishStellarModel(Model):
    def __init__(self,emulator, wavelength, distance = 100.,
                 lims = np.array([[-1, 1],
                                  [6000,8000],
                                  [4,5],
                                  [0,3],
                                  [2.9, 4.0]]),
                 **kwargs):
        self.emulator = emulator
        self.wavelength = wavelength
        self.distance = distance
        self.parLabels=["log(L)", "T_eff", "log(g)", #"log(Fe/H)",
                        "A(V)", "R(V)"]
        self.npars = 5
        self.npars_ptform = 5
        self.priors = [uniform(lim[0], lim[1]-lim[0]) for lim in lims]

    def __call__(self, logl, teff, logg, av, rv):
        #First get flux in starfish's strange system, which we have decided is erg/s/cm^2/AA
        fl = self.emulator.load_flux([teff, logg, 0.0])
        #Now to convert to a more useful unit, we calculate the bolometric flux of the model, and we will rescale it to the desired luminosity
        fbol_init = np.trapz(fl, self.emulator.wl)
        nu = c/self.emulator.wl
        scale = (fbol_1l1p/fbol_init)*(10**logl)
        fl = fl*scale * (3.34e5 * self.emulator.wl**2)
        fl = fl/(4*np.pi*self.distance**2)
        #Now we have to extrapolate to a wider wavelength range, exactly how is not too important
        #print(fl)
        ext = ccm89(self.wavelength*1e4, av, rv)
        #print(ext)
        ff = interp1d(np.log10(self.emulator.wl/1e4), np.log10(fl), fill_value=np.nan, bounds_error=False) # A simple interpolator works well within the data, but extrapolates badly. we will instead take a polynomial that passes through the data and has the form 2lambda+c at the end
        fint = ff(np.log10(self.wavelength))#*1e4))
        
        #fint[np.isnan(fint)] = self.wavelength[np.isnan(fint)] - self.emulator.wl[-5]
        #print(fint)
        i = 1
        fint[np.isnan(fint)] = fint[~np.isnan(fint)][-i] + np.log10((self.wavelength[np.isnan(fint)]/(self.emulator.wl[-i]/1e4))**-2)
        fint = 10**fint
        
        #print(fint)
        
        modelFlux = apply(ext,fint)
        return {"spectrum":{"wavelength":self.wavelength, "flux": modelFlux}}

    def lnprior(self, theta, **kwargs):
        return np.sum([self.priors[i].logpdf(theta[i]) for i in range(len(theta))])#t_lnprior+logm_lnprior+beta_lnprior+dist_lnprior

    def prior_transform(self, u, **kwargs):
        theta = np.zeros_like(u)
        theta = np.array([
                            self.priors[i].ppf(u[i]) for i in range(len(u))
                         ])
        return theta



if __name__=="__main__":
    path = './phoenix_models/'
    ranges = [
        [5700, 8000], #T
        [4.0, 5.0], #log g
        [0., 0.5] #We'll stick with solar metallicity here
    ]
    download_models = False #Should be True if this is the first time you're running this script, others False
    if download_models:
        #First, we download a subsection of the PHOENIX grid using the tools that STARFISH provides.
        
        dpm(path, ranges=ranges)
    
    grid = PHOENIX(path)
    #grid.load_flux(header=True)

    recreate_hdf5 = False
    if recreate_hdf5:
        creator = HDF5Creator(grid, "ampere_test_grid.hdf5", ranges = ranges)
        creator.process_grid()
    
    retrain_emulator = False
    if retrain_emulator:
        emu = Emulator.from_grid("ampere_test_grid.hdf5")
        emu.train()
        emu.save("ampere_test_emulator.hdf5")
    else:
        emu = Emulator.load("ampere_test_emulator.hdf5")
    
    #Now we have a spectral emulator that can be called to interpolate to arbitrary stellar parameters within the bounds of the grid it was trained on.
    #print(emu)
    test_params = [7054, 4.0324, 0.01]
    #print(emu.load_flux(test_params, norm=True))
    #print(emu.wl)
    #print(emu.norm_factor(test_params))
    #print(np.trapz(emu.load_flux(test_params, norm=True), emu.wl))
    
    #print(fbol_1l1p)#.to(units.erg/units.s/units.cm/units.cm))#a.to(units.erg/units.s/units.cm/units.cm))
    theta = [np.log10(4.68), 6750, 4., 0.0]
    fl = emu.load_flux(theta[1:])
    #Now to convert to a more useful unit, we calculate the bolometric flux of the model, and we will rescale it to the desired luminosity
    fbol_init = np.trapz(fl, emu.wl)
    nu = c/emu.wl
    scale = (fbol_1l1p/fbol_init)*(10**theta[0])
    #print(scale)
    fl = fl*scale * (3.34e5 * emu.wl**2)
    fl /= (4*np.pi*454**2)
    #plt.plot(np.log10(emu.wl), np.log10(fl))
    #print(fbol_init)
    #print(scale)
    print(fl)
    
    #print(np.trapz(emu.load_flux([7999, 4.0324, 0.01]), emu.wl))
    #print(np.trapz(emu.load_flux([5701, 4.0324, 0.01]), emu.wl))
    
    #plt.loglog(emu.wl, emu.load_flux(test_params))
    #plt.show()
    wavelengths = emu.wl/10000

    

    import pyphot
    #Now we create synthetic photometry
    filterName = np.array(['Gaia_BP', 'Gaia_G', 'Gaia_RP', "2MASS_J", "2MASS_H", "2MASS_Ks", 'WISE_RSR_W1', 'WISE_RSR_W2',])

    libDir = ampere.__file__.strip('__init__.py') 
    libname = libDir + 'ampere_allfilters.hd5'
    filterLibrary = pyphot.get_library(fname=libname)
    filters = filterLibrary.load_filters(filterName, interp=True, lamb = wavelengths*pyphot.unit['micron'])
    fwaves = np.array([filt.lpivot.to("micron").value for filt in filters])#*u.micron
    modSed = np.zeros_like(filters)
    flam = fl / wavelengths**2#res["spectrum"]["flux"] / res["spectrum"]["wavelength"]**2
    for i, (f, lp) in enumerate(zip(filters, fwaves)):
        fphot = f.get_flux(wavelengths*pyphot.unit['micron'],
                           flam*pyphot.unit['flam'], axis=-1).value
        modSed[i] = (fphot*lp**2)
        print(filterName[i], modSed[i])

    input_noise_phot = 0.05 #Fractional uncertainty
    photunc = input_noise_phot * modSed #Absolute uncertainty
    modSed = modSed + np.random.randn(len(filterName)) * photunc #Now perturb data by drawing from a Gaussian distribution

    wavelengths = 10**np.linspace(-.5,1., 2000)

    model = StarfishStellarModel(emu, wavelengths)
    theta = [np.log10(4.68), 6750, 4.37, 0.0 ]
    theta_true = theta[:-1]
    theta_true.extend([1.0, 3.2])
    print(theta_true)
    res = model(*theta_true)
    print(res)
    #exit()
    f_true = res['spectrum']['flux']

    modSed = np.zeros_like(filters)
    flam = res["spectrum"]["flux"] / res["spectrum"]["wavelength"]**2
    for i, (f, lp) in enumerate(zip(filters, fwaves)):
        fphot = f.get_flux(wavelengths*pyphot.unit['micron'],
                           flam*pyphot.unit['flam'], axis=-1).value
        modSed[i] = (fphot*lp**2)

    input_noise_phot = 0.05 #Fractional uncertainty
    photunc = input_noise_phot * modSed #Absolute uncertainty
    modSed = modSed + np.random.randn(len(filterName)) * photunc #Now perturb data by drawing from a Gaussian distribution

    from ampere.data import Photometry
    photometry = Photometry(filterName=filterName, value=modSed, uncertainty=photunc, photunits='Jy', libName=libname)
    #print(photometry.filterMask)
    photometry.reloadFilters(wavelengths)
    print(photometry)
    dataset = [photometry]

    #from ampere.infer.emceesearch import EmceeSearch
    #from emcee import moves

    #m = [(moves.DEMove(), 0.8),
    #    (moves.DESnookerMove(), 0.2),
    #     ]

    #optimizer = EmceeSearch(model=model, data=dataset, nwalkers=100, moves=m, vectorize = False, name = "star_test_emcee", namestyle="short")
    #optimizer.optimise(nsamples = 150, burnin=100, guess='None'
    #                   )
    #optimizer.postProcess()

    from ampere.infer.sbi import SBI_SNPE
    #optimizer = SBI_SNPE(model=model, data=dataset, name = "star_test_sbi", namestyle="short")
    #optimizer.optimise(nsamples = 10000, nsamples_post = 10000
    #                   )

    #optimizer.postProcess(logx=True, logy=True)

    # exit()

    #Now we attempt to generate a synthetic optical spectrum to see how that changes the outcome
    #Let's decide we have something similar to a Gaia RVS spectrum (R~11000, 842-872 nm)
    r = 11000
    lam_start = 0.842
    lam_end = 0.872

    specwaves = [lam_start]
    while specwaves[-1] < lam_end:
        specwaves.append(specwaves[-1]*(1+(1/r)))

    specwaves = np.array(specwaves)
    spec0 = spectres(specwaves,wavelengths,res['spectrum']['flux'])
    #And again, add some noise to it
    input_noise_spec = 0.01
    unc0 = input_noise_spec*spec0
    spec0 = spec0 + np.random.randn(len(spec0))*unc0
    from ampere.data import Spectrum
    spec0 = Spectrum(specwaves, spec0, unc0,"um", "Jy",calUnc=0.0025, scaleLengthPrior = 0.01) #, resampleMethod=resmethod)

    #Now let's try changing the resampling method so it's faster
    #This model is very simple so exact flux conservation is not important
    resmethod = "fast" #"exact"#"fast"#
    spec0.setResampler(resampleMethod=resmethod)

    dataset = [photometry,
               spec0, 
               ]

    optimizer = SBI_SNPE(model=model, data=dataset, name = "star_test_sbi_spec", namestyle="short")
    optimizer.optimise(nsamples = 10000, nsamples_post = 10000
                       )

    optimizer.postProcess(logx=True, logy=True)

    import pickle
    with open("pickle_pheonix_star_spec.pkl", 'wb') as f:
        pickle.dump(optimizer, f) #This way we can read it back in later and make some more convincing plots if necessary!
