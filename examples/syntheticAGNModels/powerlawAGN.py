import numpy as np
import ampere
from ampere.data import Spectrum, Photometry
from ampere.emceesearch import EmceeSearch
from ampere.PowerLawAGN import SingleModifiedBlackBody, PowerLawAGN
import corner
import matplotlib as mpl
#mpl.use("Agg")
import matplotlib.pyplot as plt
import os
from astropy import constants as const
from astropy import units as u
from spectres import spectres
import pyphot



if __name__=="__main__":
    """ Set up the inputs for the model """
    """ wavelength grid """
    wavelengths = 10**np.linspace(0.,1.9, 2000)

    """ Choose some model parameters  """
    redshift=None
    multiplicationFactor=-2.#0.0
    powerLawIndex=0.5 #1.0
    relativeAbundances=[-1.30,-9.99,-0.29,-0.39,-9.99]#[-10.,-10.,-0.5,-0.16509,-10.]

    model = PowerLawAGN(wavelengths,redshift=0.5)
    #exit()
    model(multiplicationFactor, powerLawIndex, *relativeAbundances)
    model_flux = model.modelFlux

    """ get synthetic photometry and spectra """
    filterName = np.array(['WISE_RSR_W4', 'WISE_RSR_W3', 'WISE_RSR_W2', 'WISE_RSR_W1', 'SPITZER_MIPS_24', 'SPITZER_MIPS_70']) #

    #libDir = pyphot.__file__.strip('__init__.py')+'libs/'
    #libName = libDir + 'synphot_nonhst.hd5' #PhIReSSTARTer.hd5'

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

    print(filts,modSed)

    """ (optionally) add some noise """
    #print(modSed)
    photunc = 0.1 * modSed
    modSed = modSed + np.random.randn(6) * photunc
    #print(modSed)
    #exit()

    """ now we'll create a synthetic spectrum from the model"""
    dataDir = os.getcwd() + '/../PGQuasars/PG1011-040/'
    specFileExample = 'cassis_yaaar_spcfw_14191360t.fits'
    irsEx = Spectrum.fromFile(dataDir+specFileExample,format='SPITZER-YAAAR')
    spec0 = spectres(irsEx[0].wavelength,wavelengths,model_flux)
    spec1 = spectres(irsEx[1].wavelength,wavelengths,model_flux)
    unc0 = 0.01*spec0
    unc1 = 0.01*spec1
    spec0 = spec0 + np.random.randn(len(spec0))*unc0
    spec1 = spec1 + np.random.randn(len(spec1))*unc1
    spec0 = Spectrum(irsEx[0].wavelength, spec0, unc0,"um", "Jy",calUnc=0.0025)
    spec1 = Spectrum(irsEx[1].wavelength, spec1, unc1,"um", "Jy",calUnc=0.0025)

    """ now set up ampere to try and fit the same stuff """
    photometry = Photometry(filterName=filterName, value=modSed, uncertainty=photunc, photUnits='Jy', libName=libname)
    #print(photometry.filterMask)
    photometry.reloadFilters(wavelengths)
    
    optimizer = EmceeSearch(model=model, data=[photometry,spec0,spec1], nwalkers=100)

    optimizer.optimise(nsamples=300000, burnin=290000, guess=[[-2.5, .5, -1., -10, -0.5, -0.3, -10., 1.0, 1.0, 1.0, 1.0 ,1.0, 1.0] + np.random.rand(13)*[1,1,1,1,0.2,0.2, 1,1,1,1,1,1,1] for i in range(optimizer.nwalkers)])

    optimizer.postProcess()

    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #print(optimizer.samples.shape)
    #print(np.max(optimizer.samples, axis=0))
    #print(np.min(optimizer.samples, axis=0))
    nneg=0
    #optimizer.postProcess()
    for i in range(0,optimizer.samples.shape[0],1000):
        #print(opt.samples[i,:])
        #if opt.samples[i,0] > 0.:
        optimizer.model(optimizer.samples[i,0],optimizer.samples[i,1],*optimizer.samples[i,2:7])
        ax.plot(wavelengths,optimizer.model.modelFlux, '-', color='black', alpha=0.05) #irs.wavelength
        #else:
        #    nneg += 1
            
    #    ax.plot(irs.wavelength, model(opt.samples[i,:]), '-', alpha = 0.1)
    #print(nneg)
    for i in [spec0,spec1]:
        ax.plot(i.wavelength, i.value, '-',color='blue')
    ax.plot(photometry.wavelength, photometry.value, 'o',color='blue')
    ax.set_ylim(0., 1.5*np.max([np.max(i.value) for i in [photometry, spec0, spec1]]))
    fig.savefig("seds.png")

    print("Acceptance fractions: ",optimizer.sampler.acceptance_fraction)
    try:
        print("Estimates of the autocorelation lengths: ",optimizer.sampler.acor)
    except Exception as e:
        print(str(e))
        print("Try using more samples")
        print("Current settings - nwalkers = ",opt.nwalkers,",  nsamples = ",opt.nsamp)
        
