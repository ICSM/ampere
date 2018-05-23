import numpy as np
import ampere
from ampere.data import Spectrum, Photometry
from ampere.emceesearch import EmceeSearch
from ampere.PowerLawAGN import PowerLawAGN, SingleModifiedBlackBody
from ampere.extinction import CCMExtinctionLaw
import corner
import matplotlib.pyplot as plt

from astropy import constants as const
from astropy import units as u
from spectres import spectres

if __name__=="__main__":
    """
    This script is a demonstration of how ampere (should) work(s). 
    """


    """ Read in some data """
    filename = ampere.__file__.strip('__init__.py')+'Testdata/cassis_yaaar_spcfw_14203136t.fits'
    print(filename)
    irs = Spectrum.fromFile(filename,format='SPITZER-YAAAR')
    #If we have more than one dataset we need multiple data objects
    #print(irs.wavelength)
    
    """ Define the model """
    #modwaves = np.linspace(1.,1.6, 10000)
    #print(10**modwaves)
    #exit()
    Tbb = 250.0
    area = (10.*u.au.to(u.m)**2)#.value)^2
    print(Tbb, np.log10(area), 0., 3.)
    #exit()
    TestData = SingleModifiedBlackBody(irs.wavelength, #modwaves,flatprior=True,
                                       normWave = 20., sigmaNormWave = 100.,
                                       redshift = False)
    #TestData(250,23,0,3)
    #print(TestData.modelFlux)
    #exit()
    TestData(t = Tbb, scale = np.log10(area), index = 0., dist=3.)
    irs.value = TestData.modelFlux #spectres(irs.wavelength, modwaves, TestData.modelFlux) #['flux'].data = TestData.modelFlux
    irs.uncertainty = 0.11*irs.value#['flux'].data
    print(irs)
    
    #model = PowerLawAGN(irs.wavelengths) #We should probably define a wavelength grid separate from the data wavelengths for a number of reasons, primarily because the photometry needs an oversampled spectrum to compute synthetic photometry

    model = SingleModifiedBlackBody(irs.wavelength,flatprior=True,
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
    opt = EmceeSearch(model = model, data = [irs], nwalkers = 300) #introspection is used inside the optimiser to determine the number of parameters that each part of the model has

    """ if you want, overload the default priors with a new function """
    def lnprior(self, inputs, **kwargs):
        return 0 #flat prior
    
    #model.lnprior = lnprior

    print(opt.npars,np.int(opt.npars))
    """ Run it """
    pos = [
           [
               250., 23., 0., 3.
               #20 + np.random.randn() for i in range(np.int(opt.npars))
           ]
           + np.random.randn(np.int(opt.npars)) for j in range(opt.nwalkers)
          ]
    print(pos[0])
    print(np.max(pos, axis=0))
    print(np.min(pos, axis=0))
    opt.optimise(nsamples = 2000, burnin=1000,guess=pos)

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
            ax.plot(irs.wavelength,opt.model.modelFlux, '-', alpha=0.02)
        else:
            nneg += 1
            
    #    ax.plot(irs.wavelength, model(opt.samples[i,:]), '-', alpha = 0.1)
    print(nneg)
    ax.plot(irs.wavelength, irs.value, '-')
    fig2 = corner.corner(opt.samples)#,labels=opt.labels)
    plt.show()

    """ Figure out what it all means """
