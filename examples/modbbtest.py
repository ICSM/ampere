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
    t_in = 200. #units = K (or 300K)
    scale_in = 1.
    index_in = 0.
    distance = 1.
    lims_in=np.array([ [100., 300.], [0., 2.], [-5, 5], [0.9, 1.1]
                  ])


    """ Initialise the model  """
    model = SingleModifiedBlackBody(wavelengths, lims=lims_in)


    """ Get a test spectrum out of the model """
    model(t_in, scale_in, index_in, distance)
    model_flux = model.modelFlux
    #model_flux = model(t_in, scale_in, index_in, distance).modelFlux #Do the two lines above with just one line here. 


    """ get synthetic photometry and spectra """
    filterName = np.array(['WISE_RSR_W4', 'WISE_RSR_W3', 'WISE_RSR_W2', 'WISE_RSR_W1', 'SPITZER_MIPS_24', 'SPITZER_MIPS_70']) #

    #libDir = pyphot.__file__.strip('__init__.py')+'libs/'
    #libName = libDir + 'synphot_nonhst.hd5' #PhIReSSTARTer.hd5'

    libDir = '/home/peter/pythonlibs/ampere/ampere/'
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
    print(modSed)
    modSed = modSed + np.random.randn(6) * 0.1 * modSed
    print(modSed)
    #exit()

    """ now we'll create a synthetic spectrum from the model"""
    dataDir = os.getcwd() + '/PGQuasars/PG1011-040/'
    specFileExample = 'cassis_yaaar_spcfw_14191360t.fits'
    irsEx = Spectrum.fromFile(dataDir+specFileExample,format='SPITZER-YAAAR')
    spec0 = spectres(irsEx[0].wavelength,wavelengths,model_flux)
    spec1 = spectres(irsEx[1].wavelength,wavelengths,model_flux)
    spec0 = Spectrum(irsEx[0].wavelength, spec0, 0.1*spec0,"um", "Jy")
    spec1 = Spectrum(irsEx[1].wavelength, spec1, 0.1*spec1,"um", "Jy")

    """ now set up ampere to try and fit the same stuff """
    photometry = Photometry(filterName=filterName, value=modSed, uncertainty=0.1*modSed, photUnits='Jy', libName=libname)
    print(photometry.filterMask)
    photometry.reloadFilters(wavelengths)
    
    optimizer = EmceeSearch(model=model, data=[photometry,spec0,spec1], nwalkers=100)

    optimizer.optimise(nsamples=1000, burnin=500, guess=[[200, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0 ,1.0, 1.0] + np.random.rand(10)*[10,1,1,0.1, 1,1,1,1,1,1] for i in range(optimizer.nwalkers)])

    optimizer.postProcess()

    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    print(optimizer.samples.shape)
    print(np.max(optimizer.samples, axis=0))
    print(np.min(optimizer.samples, axis=0))
    nneg=0
    optimizer.postProcess()
    for i in range(0,optimizer.samples.shape[0],1000):
        #print(opt.samples[i,:])
        #if opt.samples[i,0] > 0.:
        optimizer.model(optimizer.samples[i,0],optimizer.samples[i,1],*optimizer.samples[i,2:7])
        ax.plot(modwaves,optimizer.model.modelFlux, '-', color='black', alpha=0.05) #irs.wavelength
        #else:
        #    nneg += 1
            
    #    ax.plot(irs.wavelength, model(opt.samples[i,:]), '-', alpha = 0.1)
    print(nneg)
    for i in data[1:]:
        ax.plot(i.wavelength, i.value, '-',color='blue')
    ax.plot(photometry.wavelength, photometry.value, 'o',color='blue')
    ax.set_ylim(0., 1.5*np.max([np.max(i.value) for i in dataSet]))
    fig.savefig("seds.png")



































