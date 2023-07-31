import ampere
from ampere import data
from ampere.infer.emceesearch import EmceeSearch
from ampere.infer.sbi import SBI_SNPE
import ampere.models.Hyperion
import numpy as np
import dill
from astropy import units

limits = np.asarray([[-10,-6],
          [-2,2],
          [2,4],
          [-2,2],
          [1,3],
          [1e3,1e4],
          [0.0,1.0],
          [0.0,1.0]])

labels = ['envelope_mass','envelope_rin','envelope_rout','envelope_r0','stellar_mass','stellar_luminosity','abundance_1','abundance_2']

model = ampere.models.Hyperion.HyperionCStarRTModel(np.logspace(np.log10(0.2),np.log10(200),1000,base=10),nproc=70,lims=limits,parLabels=labels)

libdir= ampere.__file__.strip('__init__.py')

p = data.Photometry.fromFile('Observed_SED.vot', format = 'votable', libName = libdir+ 'ampere_allfilters.hd5')

p.reloadFilters(np.logspace(np.log10(0.2),np.log10(200),1000,base=10))

s = data.Spectrum.fromFile('SPEC_OGLE_CAGB_IRS.csv', 'User-Defined', filetype = 'text', keywords = {'waveCol': 'wavelength', 'fluxCol': 'fluxdensity', 'uncsCol': 'uncertainty', 'waveUnit': 'micron', 'fluxUnit': 'Jy'})

#t = data.Spectrum.fromFile('OGLE_LMC_LPV_28579_photosphere_interpolated.vot', 'User-Defined', filetype = 'votable', keywords = {'waveCol': 'LAMINTERP', 'fluxCol': 'SPECINTERP', 'waveUnit': 'micron', 'fluxUnit': 'Jy'})

dataset = [p]
for sp in s:
    dataset.append(sp)

embedding = {'type': 'FC', 'num_hiddens': 100, 'n_layers': 3, "output_dim": 28}
optimizer_sbi = SBI_SNPE(model=model, data=dataset, name="Cstar_test_sbi_v3",
                         namestyle="short", n_prior_norm_samples=100000,
                         get_prior_bounds=False, n_rounds=2,
                         embedding_net=embedding)

optimizer_sbi.optimise(nsamples = 10000, nsamples_post = 10000)

model.HyperionRTSED = 0.0
model.sed = 0.0
model.model = 0.0
model.envelope_shell = 0.0

optimizer_sbi.postProcess(logx=True,logy=True)

#save object before crashing.
with open('sundar_hyperion_test_sbi_alternate_comp_and_photosphere_OUTPUT.pkl', 'wb') as file:
        dill.dump(optimizer_sbi,file)

exit()

#emcee
optimizer = EmceeSearch(model=model,data=dataset,nwalkers=40,vectorize=False)

guess = [0.9,1,0.1,0.1,1,0.1,0.1]

optimizer.optimise(nsamples=1000,nburnin=250,guess=guess)

optimizer.summary()

optimizer.print_summary(outfile='sundar_hyperion_test_summary.txt')

optimizer.postProcess(logx=True,logy=True)

#save object before crashing.
with open('sundar_hyperion_test_sbi_whole_new_comp_and_photosphere_OUTPUT.pkl', 'wb') as file:
    dill.dump(optimiser_emcee,file)
