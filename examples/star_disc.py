import ampere
from ampere import data
from ampere.infer.emceesearch import EmceeSearch
from ampere import QuickSED
from astropy import units, constants
import numpy as np
import scipy

libdir= ampere.__file__.strip('__init__.py')


#Load data for fitting process - photometry, irs spectrum, stellar spectrum


#Photometry
photometry = data.Photometry.fromFile('./star_disc/HD105_SED.vot', format = 'votable',libName = libdir+ 'ampere_allfilters.hd5')
photometry.reloadFilters(np.logspace(np.log10(0.2),np.log10(500),1000,base=10))

#IRS spectrum
irs_spec = data.Spectrum.fromFile('./star_disc/HD105_IRS_Spectrum.csv', 'User-Defined', filetype = 'text', keywords = {'waveCol': 'wavelength', 'fluxCol': 'fluxdensity', 'uncsCol': 'uncertainty', 'waveUnit': 'micron', 'fluxUnit': 'Jy'})

#Stellar spectrum - we'll also need some of this for the fitting
import Starfish 
from Starfish.emulator import Emulator
from Starfish.grid_tools import GridInterface
from Starfish.grid_tools import PHOENIXGridInterfaceNoAlpha as PHOENIX #NoAlpha
import Starfish.constants as C
from Starfish.utils import create_log_lam_grid
from Starfish.grid_tools.utils import vacuum_to_air, idl_float
from Starfish.grid_tools import download_PHOENIX_models as dpm
from Starfish.grid_tools import HDF5Creator

#Conversions for Starfish
fbol_1l1p = 1*units.Lsun
fbol_1l1p = ((fbol_1l1p/(4*np.pi*(1*units.pc)**2)).to(units.erg/units.s/units.cm/units.cm)).value
c = constants.c.to(units.AA/units.s)

path = libdir + '../examples/examples_paper/phoenix_models/'
print(path)
download_models = False #Should be True if this is the first time you're running this script, others False
if download_models:
    #First, we download a subsection of the PHOENIX grid using the tools that STARFISH provides.
    ranges = [
        [5000, 8000], #T
        [4.0, 5.0], #log g
        [0., 0.5] #We'll stick with solar metallicity here
    ]

    dpm(path, ranges=ranges)

grid = PHOENIX(path)
#grid.load_flux(header=True)

recreate_hdf5 = False
if recreate_hdf5:
    creator = HDF5Creator(grid, "ampere_test_grid.hdf5", ranges = ranges)
    creator.process_grid()

retrain_emulator = True
if retrain_emulator:
    emu = Emulator.from_grid(path + "../" + "ampere_test_grid.hdf5")
    emu.train()
    emu.save(path + "../" + "ampere_test_emulator.hdf5")
else:
    emu = Emulator.load(path + "../" + "ampere_test_emulator.hdf5")

#Now we have a spectral emulator that can be called to interpolate to arbitrary stellar parameters within the bounds of the grid it was trained on.

#Use the Gaia DR2 distance to HD 105 as a fixed value, which should constrain the degeneracy between radius and luminosity
parallax = 25.7534
parallax_error = 0.0211

#Create a synthetic spectrum for HD 105 based on the fitted parameters from Marshall et al. 2018
hd105_params = [6034, 4.478, 0.02]
print(emu.load_flux(hd105_params, norm=True))
print(emu.wl)
print(emu.norm_factor(hd105_params))
print(np.trapz(emu.load_flux(hd105_params, norm=True), emu.wl))

theta = [np.log10(1.216), 6034, 4.478, 0.02]
fl = emu.load_flux(theta[1:])
#Now to convert to a more useful unit, we calculate the bolometric flux of the model, and we will rescale it to the desired luminosity
fbol_init = np.trapz(fl, emu.wl)
nu = c/emu.wl
scale = (fbol_1l1p/fbol_init)*(10**theta[0])
fl = fl*scale * (3.34e5 * emu.wl**2)
fl /= (4*np.pi*(1000/parallax)**2)
#print(fl)

wfl = np.linspace(4e3,1.2e4,int(8e2),endpoint=True)

ffl = scipy.interpolate.interp1d(emu.wl,fl)

fl = ffl(wfl)

wav = wfl/10000 #A -> um

#Add Gaussian noise to calculated spectrum assuming s/n 100 spectrum
unc = np.zeros(len(fl))
for i in range(0,len(fl)):
	unc[i] = fl[i]*0.01*np.random.normal()

star_spec = data.Spectrum(wfl,fl,unc,"um","Jy",calUnc=0.0025, scaleLengthPrior = 0.01)

dataSet = [photometry,star_spec,irs_spec]

#Now we get to the searching
#Create a wavelength grid
wavelengths = np.logspace(np.log10(0.1),np.log10(500),1000,endpoint=True)

model = ampere.QuickSED.QuickSEDModel(wavelengths,lims=np.array([[0.5,1.5],[5500.,6500.],[4.0,5.0],[0.0,0.1],[-30,30],[30.,300.],[50.,500.],[-2.,2.]]),dstar=1000./parallax)

print(model.__repr__)

optimizer = EmceeSearch(model=model,data=dataSet,nwalkers=20,vectorize=False)

#self.parLabels = ['lstar','tstar','log_g','[fe/h]','Adust','tdust','lam0','beta'] + noise parameters
guess = [1.0,5784.0,4.5,0.01,-3,50.0,200.0,1.0]#,1.0,0.5,1.0,1.0,0.5,1.0,1.0,0.5,1.0,1.0,0.5,1.0,1.0,0.5,1.0,1.0,0.5,1.0]

optimizer.optimise(nsamples=1000,nburnin=250,guess=guess + np.random.rand(10)*[10,1,1,0.1, 1,1,1,1,1,1])

optimizer.summary()

optimizer.print_summary(outfile='jonty_quicksed_test_summary.txt')

optimizer.postProcess(logx=True,logy=True)

#save object before crashing.
import dill
with open('circumstellar_disc.pkl', 'wb') as f:
	dill.dump(optimizer,f)