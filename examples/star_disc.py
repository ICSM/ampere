import ampere
from ampere import data
from ampere.infer.emceesearch import EmceeSearch
import ampere.QuickSED
import numpy as np
import dill

libdir= ampere.__file__.strip('__init__.py')

p = data.Photometry.fromFile('HD105_SED.vot', format = 'votable',libName = libdir+ 'ampere_allfilters.hd5')

p.reloadFilters(np.logspace(np.log10(0.2),np.log10(500),1000,base=10))

s = data.Spectrum.fromFile('HD105_IRS_Spectrum.csv', 'User-Defined', filetype = 'text', keywords = {'waveCol': 'wavelength', 'fluxCol': 'fluxdensity', 'uncsCol': 'uncertainty', 'waveUnit': 'micron', 'fluxUnit': 'Jy'})

dataset = [p]
#for sp in s:
#	dataset.append(sp)

wavelengths = np.logspace(np.log10(0.1),np.log10(500),1000,endpoint=True)


parallax = 25.7534
parallax_error = 0.0211

model = ampere.QuickSED.QuickSEDModel(wavelengths,lims=np.array([[(1000./(parallax+5*parallax_error)),1000./(parallax-5*parallax_error)],[5500.,6500.],[-1,1],[-30,30],[30.,300.],[50.,500.],[-2.,2.]]))

print(model.__repr__)


optimizer = EmceeSearch(model=model,data=dataset,nwalkers=20,vectorize=False)

#self.parLabels = ['dstar','tstar','rstar','Adust','tdust','lam0','beta'] + noise parameters
guess = [1000./parallax,6000.0,0.0,-3,50.0,200.0,1.0]#,1.0,0.5,1.0,1.0,0.5,1.0,1.0,0.5,1.0,1.0,0.5,1.0,1.0,0.5,1.0,1.0,0.5,1.0]

optimizer.optimise(nsamples=10000,nburnin=2500,guess=guess)

optimizer.summary()

optimizer.print_summary(outfile='jonty_quicksed_test_summary.txt')

optimizer.postProcess(logx=True,logy=True)

#save object before crashing.
with open('jonty_quicksed_test.pkl', 'wb') as file:
	dill.dump(optimizer.samples,file)