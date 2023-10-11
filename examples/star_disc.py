import ampere
import pyphot
from pyphot import unit, Filter
from spectres import spectres
from ampere import data
from ampere.infer.emceesearch import EmceeSearch
from ampere.models import QuickSED
from astropy import units, constants
import numpy as np
import scipy

libdir= ampere.__file__.strip('__init__.py')


#Load data for fitting process - photometry, irs spectrum, stellar spectrum

#Create the ATCA photometry to the data set manually, as they don't have PyPhot filters
#Here we approximate them as simple tophat functions rather than taking into account the full 
#details of the spectral windows used in the original observations

#Photometry
photometry = data.Photometry.fromFile('./star_disc/HD105_SED.vot', format = 'votable',libName = libdir+ 'ampere_allfilters_extended.hd5')
#photometry.reloadFilters(np.logspace(np.log10(0.2),np.log10(10000),1000,base=10))




#IRS spectrum
irs_spec = data.Spectrum.fromFile('./star_disc/cassis_yaaar_spcfw_5295616t.fits', format='SPITZER-YAAAR', filetype = 'fits')

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
#print(path)
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

retrain_emulator = False
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
#print(emu.load_flux(hd105_params, norm=True))
#print(emu.wl)
#print(emu.norm_factor(hd105_params))
#print(np.trapz(emu.load_flux(hd105_params, norm=True), emu.wl))

theta = [np.log10(1.216), 6034, 4.478, 0.02]
fl = emu.load_flux(theta[1:])
#Now to convert to a more useful unit, we calculate the bolometric flux of the model, and we will rescale it to the desired luminosity
fbol_init = np.trapz(fl, emu.wl)
nu = c/emu.wl
scale = (fbol_1l1p/fbol_init)*(10**theta[0])
fl = fl*scale * (3.34e5 * emu.wl**2)
fl /= (4*np.pi*(1000/parallax)**2)
#print(fl)

#Create s/n 100 synthetic stellar spectrum that looks like Gaia RV spec
r = 11000
lam_start = 0.847
lam_end = 0.871

specwaves = [lam_start]
while specwaves[-1] < lam_end:
    specwaves.append(specwaves[-1]*(1+(1/r)))

specwaves = np.array(specwaves)
spec0 = spectres(specwaves,emu.wl*1e-4,1.216*fl)
#And again, add some noise to it
input_noise_spec = 0.05
unc0 = input_noise_spec*spec0
spec0 = spec0 + np.random.randn(len(spec0))*unc0
star_spec = data.Spectrum(specwaves, spec0, unc0,"um", "Jy",calUnc=0.10, scaleLengthPrior = 0.01) #, resampleMethod=resmethod)

dataSet = [photometry,star_spec,
		   ]

for irs in irs_spec:
	dataSet.append(irs)

#Now we get to the searching
#Create a wavelength grid
wavelengths = np.logspace(np.log10(0.31),np.log10(10000),1000,endpoint=True)

model = QuickSED.QuickSEDModel(wavelengths,
							   lims=np.array([[0.7,1.3],[5500.,6500.],[4.1,4.9],[0.01,0.1],[-3,3],[30.,100.],[50.,500.],[0.,4.]]),
							   dstar=1000./parallax,
							   starfish_dir=path+'../',
							   fbol_1l1p=fbol_1l1p)

#print(model.__repr__)

optimizer = EmceeSearch(model=model,data=dataSet,nwalkers=40,vectorize=False)

#self.parLabels = ['lstar','tstar','log_g','[fe/h]','Adust','tdust','lam0','beta'] + noise parameters
guess = [1.0,5900.0,4.3,0.05,0.5,50.0,100.0,1.0,1.0,0.5,1.0,1.0,0.5,1.0]#,1.0,0.5,1.0,1.0,0.5,1.0,1.0,0.5,1.0,1.0,0.5,1.0,1.0,0.5,1.0,1.0,0.5,1.0]

optimizer.optimise(nsamples=4000,nburnin=1000,guess=guess)

optimizer.summary()

optimizer.print_summary(outfile='jonty_quicksed_test_summary.txt')

optimizer.postProcess(logx=True,logy=True)

#save object before crashing.
import dill
with open('star_disc_PRODUCTION_v1.pkl', 'wb') as f:
	dill.dump(optimizer,f)