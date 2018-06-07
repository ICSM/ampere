import numpy as np
import ampere
from ampere.data import Spectrum, Photometry
from ampere.emceesearch import EmceeSearch
from ampere.PowerLawAGN import PowerLawAGN
import corner
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import os
from astropy import constants as const
from astropy import units as u
from spectres import spectres
import pyphot


if __name__=="__main__":
    
    ''' Let's start with PG 1011-040 as an example '''
    dataDir = os.getcwd() + '/PGQuasars/PG1011-040/'
    specFile = 'cassis_yaaar_spcfw_14191360t.fits'
    photFile = 'vizier_votable_pruned.vot'
    irs = Spectrum.fromFile(dataDir+specFile,format='SPITZER-YAAAR')
    libDir = '/home/peter/pythonlibs/ampere/ampere/'
    libname = libDir + 'ampere_allfilters.hd5'
    phot = Photometry.fromFile(dataDir+photFile, libName = libname)

    modwaves = 10**np.linspace(0.,1.9, 2000)

    model = PowerLawAGN(modwaves, redshift=0.058)
    phot.reloadFilters(modwaves)
    dataSet = [phot]
    for s in irs:
        dataSet.append(s)
    for s in dataSet:
        print(s)

    opt = EmceeSearch(model = model, data = dataSet, nwalkers = 200)

    print(opt.npars)

    pos = [
           [
               0., 2., 0.1, 0.1, 0.1, 0.1, 0.1, 1., 0.5, 1., 1., 0.5, 1.
               #20 + np.random.randn() for i in range(np.int(opt.npars))
           ]
           + np.random.randn(np.int(opt.npars)) for j in range(opt.nwalkers)
          ]
    for i in range(len(pos)):
        pos[i][0] = pos[i][0] / 1000.
        pos[i][2:7] = pos[i][2:7] / 30. + 0.1
    print(pos[0])
    print(np.max(pos, axis=0))
    print(np.min(pos, axis=0))
    opt.optimise(nsamples = 20000, burnin=10000,guess=pos)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    print(opt.samples.shape)
    print(np.max(opt.samples, axis=0))
    print(np.min(opt.samples, axis=0))
    nneg=0
    opt.postProcess()
    for i in range(0,opt.samples.shape[0],100):
        #print(opt.samples[i,:])
        if opt.samples[i,0] > 0.:
            opt.model(opt.samples[i,0],opt.samples[i,1],*opt.samples[i,2:7])
            ax.plot(modwaves,opt.model.modelFlux, '-', color='black', alpha=0.02) #irs.wavelength
        else:
            nneg += 1
            
    #    ax.plot(irs.wavelength, model(opt.samples[i,:]), '-', alpha = 0.1)
    print(nneg)
    for i in irs:
        ax.plot(i.wavelength, i.value, '-',color='blue')
    ax.plot(phot.wavelength, phot.value, 'o',color='blue')
    ax.set_ylim(0., 1.5*np.max([np.max(i.value) for i in dataSet]))
    fig.savefig("seds.png")
    #fig2 = corner.corner(opt.samples)#,labels=opt.labels)
    #plt.show()
