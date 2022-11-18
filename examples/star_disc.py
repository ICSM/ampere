import ampere
from ampere import data
from ampere.infer.emceesearch import EmceeSearch
from ampere.infer.sbi import SBI_SNPE
import ampere.QuickSED
import numpy as np
import dill

libdir= ampere.__file__.strip('__init__.py')

p = data.Photometry.fromFile('./star_disc/HD105_SED.vot', format = 'votable',libName = libdir+ 'ampere_allfilters.hd5')

p.reloadFilters(np.logspace(np.log10(0.2),np.log10(500),1000,base=10))

s = data.Spectrum.fromFile('./star_disc/HD105_IRS_Spectrum.csv', 'User-Defined', filetype = 'text', keywords = {'waveCol': 'wavelength', 'fluxCol': 'fluxdensity', 'uncsCol': 'uncertainty', 'waveUnit': 'micron', 'fluxUnit': 'Jy'})

dataset = [p]
#for sp in s:
#	dataset.append(sp)

wavelengths = np.logspace(np.log10(0.1),np.log10(500),1000,endpoint=True)

model = ampere.QuickSED.QuickSEDModel(wavelengths,dstar=40.0,tstar=5770.0,rstar=1.0,fdisc=1e-4,\
                        tdisc=50.0,lam0=200.0,beta=1.5)

print(model.__repr__)

optimizer = EmceeSearch(model=model,data=dataset,nwalkers=20,vectorize=False)

#self.parLabels = ['dstar','tstar','rstar','fdust','tdust','lam0','beta'] + noise parameters
guess = [20.0,6000.0,1.0,1e-4,100.0,200.0,1.0]#,1.0,0.5,1.0,1.0,0.5,1.0,1.0,0.5,1.0,1.0,0.5,1.0,1.0,0.5,1.0,1.0,0.5,1.0]

optimizer.optimise(nsamples=10000,nburnin=2500,guess=guess)

optimizer.summary()

optimizer.print_summary(outfile='jonty_quicksed_test_summary.txt')

optimizer.postProcess()

#save object before crashing.
with open('jonty_quicksed_test.pkl', 'wb') as file:
	dill.dump(optimiser.samples,file)
