import numpy as np
import sys
sys.path.insert(1, '/home/zeegers/git_ampere/ampere/')
import ampere
from ampere.data import Spectrum, Photometry
from ampere.emceesearch import EmceeSearch
from ampere.PowerLawAGN import PowerLawAGN, SingleModifiedBlackBody
from ampere.extinction import CCMExtinctionLaw
import corner
import matplotlib.pyplot as plt
import os

from astropy import constants as const
from astropy import units as u
from spectres import spectres
import pyphot

if __name__=="__main__":
    """
    This script is a demonstration of how ampere (should) work(s). 
    """


    """ Read in some data """
    filename = ampere.__file__.strip('__init__.py')+'Testdata/all_low_res2_nolines_hd.fits' # path to data of cyg ob 2 12
    print(filename)
    irs = Spectrum.fromFile(filename,format='SPITZER-YAAAR')
    #for s in irs:
    #    print(s)
    #exit()
    #If we have more than one dataset we need multiple data objects
    #print(irs.wavelength)
    
    """ Define the model """
    modwaves = 10**np.linspace(0.,2., 2000)
    #print(modwaves)
    #print(10**modwaves)
    #exit()
    Tbb = 250.0
    area = (10.*u.au.to(u.m)**2)#.value)^2
    print(Tbb, np.log10(area), 0., 3.)
    #exit()
    TestData = SingleModifiedBlackBody(modwaves, # irs.wavelength, #modwaves,flatprior=True,
                                       normWave = 20., sigmaNormWave = 100.,
                                       redshift = False)
    #TestData(250,23,0,3)
    #print(TestData.modelFlux)
    #exit()
    print(np.min(modwaves),np.max(modwaves))
    TestData(t = Tbb, scale = np.log10(area), index = 0., dist=3.)
    for i in irs:
        #print(i)
        i.value = spectres(i.wavelength, modwaves, TestData.modelFlux)
        #print(i)
        i.uncertainty = 0.11*i.value
        #print(i)
        #irs.value = TestData.modelFlux #spectres(irs.wavelength, modwaves, TestData.modelFlux) #['flux'].data = TestData.modelFlux
    #irs.uncertainty = 0.11*irs.value#['flux'].data
    #print(irs)
    #exit()

    ''' Now we create some synthetic photometry too *** look for the SED information of Schulte 12 ''' 
    filters=[
             #'SPITZER_IRAC_36',
             'WISE_RSR_W1',
             'WISE_RSR_W2',
             'WISE_RSR_W3',
             'WISE_RSR_W4',
             'HERSCHEL_PACS_BLUE'
             ]
    libDir = os.getcwd() + '/ampere/'
#    libname = libDir + 'ampere_allfilters.hd5'
    libname = '/home/zeegers/ampere/ampere/ampere_allfilters.hd5'
    print(libname)
    filterLib = pyphot.Library.from_hd5(libname)
    filts = filterLib.load_filters(filters, interp=True, lamb = modwaves*pyphot.unit['micron'])
    f, sed =  pyphot.extractPhotometry(modwaves, TestData.modelFlux, filts, Fnu = True, absFlux = False)
    lamCen = np.array([b.magnitude for b in f])
    print(lamCen, sed)

    ''' now turn the numbers into a Photometry object '''
    phot = Photometry(np.array(filters), np.array(sed), np.array(0.1*sed), np.array(['Jy', 'Jy','Jy','Jy','Jy']), bandUnits='um',libName = libname)
    phot.reloadFilters(modwaves)
    print(phot)
    irs.append(phot)
    #exit()
    #model = PowerLawAGN(irs.wavelengths) #We should probably define a wavelength grid separate from the data wavelengths for a number of reasons, primarily because the photometry needs an oversampled spectrum to compute synthetic photometry

    model = SingleModifiedBlackBody(modwaves, #irs.wavelength,flatprior=True,
                                    normWave = 20., sigmaNormWave = 100.,
                                    redshift = False, lims = np.array([[150.,400.],
                                                                       [15.,35.],
                                                                       [-2.,2.],
                                                                       [1.,5.]]
                                                                      )
                                    )
    
    """ and the extinction law (if necessary) (could be moved inside model as only some models will need this) """
    #ext = CCMExtinctionLaw()

    #model.extinction = ext #combination of multiple different models is something we need to consider how to implement. 

    """ Connect the dots: """
    """ Hook it up to an optimiser """
    opt = EmceeSearch(model = model, data = irs, nwalkers = 30) #introspection is used inside the optimiser to determine the number of parameters that each part of the model has

    """ if you want, overload the default priors with a new function """
    def lnprior(self, inputs, **kwargs):
        return 0 #flat prior
    
    #model.lnprior = lnprior

    print(opt.npars,np.int(opt.npars))
    """ Run it """
    pos = [
           [
               250., 23., 0., 3., 1., 0.5, 1., 1., 0.5, 1.
               #20 + np.random.randn() for i in range(np.int(opt.npars))
           ]
           + np.random.randn(np.int(opt.npars)) for j in range(opt.nwalkers)
          ]
    print(pos[0])
    print(np.max(pos, axis=0))
    print(np.min(pos, axis=0))
    opt.optimise(nsamples = 1000, burnin=100,guess=pos)

    """save optimiser state for later """
    #opt.save(filename="output_file",pickle=True) #Save as a python object
    
    """ Sort out the results """
    """ First produce some terminal output for the numbers and confidence intervals """

    """ Then produce a couple of plots """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    print(opt.samples.shape)
    print(np.max(opt.samples, axis=0))
    print(np.min(opt.samples, axis=0))
    nneg=0
    for i in range(0,opt.samples.shape[0],100):
        #print(opt.samples[i,:])
        if opt.samples[i,0] > 0.:
            opt.model(opt.samples[i,0],opt.samples[i,1],opt.samples[i,2],opt.samples[i,3])
            ax.plot(modwaves,opt.model.modelFlux, '-', 'k', alpha=0.02) #irs.wavelength
        else:
            nneg += 1
            
    #    ax.plot(irs.wavelength, model(opt.samples[i,:]), '-', alpha = 0.1)
    print(nneg)
    for i in irs:
        ax.plot(i.wavelength, i.value, '-','b')
    fig2 = corner.corner(opt.samples)#,labels=opt.labels)
    plt.show()

    """ Figure out what it all means """
