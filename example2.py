#F(T) = (1/D^2) (4 pi a^2 r0^3 n0) / (3 - p) (T/T0)^(- (3-p)/q)  Q Bnu(T)
#2x temperature ranges T1 - T2, T3 - T4 (30-60K, 100-118K)


import numpy as np
import ampere
from ampere.data import Spectrum, Photometry
from ampere.emceesearch import EmceeSearch
from ampere.PowerLawAGN import OpacitySpectrum
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
    This script will take one blackbody and modify it by a wavelength dependent opacity.
    """

    """ Read in some data """
    
    """ Read in spectrum """
    filename = 'Users/jonty/Documents/GitHub/ampere/examples/NGC6302/NGC6302_100.tab'
    print(filename)
    spec = Spectrum.fromFile(filename,format='User-defined')

    """ Define the model """

    modwaves = 10**np.linspace(2.,200., 2000) #2-200 um spectrum w/ 2000 data points
    Tbb = 45.0 #cold component from Kemper et al. 2002 is 30-60K, warm component is 100-118K
    distance = 140.0*u.pc #assume 140pc distance
    area = (10.*u.au.to(u.m)**2)#.value)^2
	
    """ Generate model """
    #modfied blackbody multiplied by sum of opacties
    opacities = ['ss_Dorschneretal1995_Olivine_0.10.q',
    			 'ss_Jaegeretal1998_Enstatite_0.10.q',
    			 'ss_Jaegeretal1998_Forsterite_0.10.q'] # list of opacity spectra to be used in modelling
    relativeAbundances=[0.01,0.01,0.01]#initial guess
    model = OpacitySpectrum(modwaves,
                            normWave = 1., sigmaNormWave = 1.,
                            opacityFileList=opacities,
                            weights=relativeAbundances,
                            redshift = False, lims = np.array([[100.,200.],
                                                               [30.,60.],
                                                               [-2.,2.],
                                                               [1.,5.]]
                                                              )
                            )

    """ Connect the dots: """
    """ Hook it up to an optimiser """
    opt = EmceeSearch(model = model, data = spec, nwalkers = 100) #introspection is used inside the optimiser to determine the number of parameters that each part of the model has

    """ if you want, overload the default priors with a new function """
    def lnprior(self, inputs, **kwargs):
        return 0 #flat prior
    
    #model.lnprior = lnprior

    print(opt.npars,np.int(opt.npars))
    """ Run it """
    pos = [
           [
               150., 45., 0., 3., 1., 0.5, 1., 1., 0.5, 1.
               #20 + np.random.randn() for i in range(np.int(opt.npars))
           ]
           + np.random.randn(np.int(opt.npars)) for j in range(opt.nwalkers)
          ]
    print(pos[0])
    print(np.max(pos, axis=0))
    print(np.min(pos, axis=0))
    opt.optimise(nsamples = 1000, burnin=500,guess=pos)

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
    for i in spec:
        ax.plot(i.wavelength, i.value, '-','b')
    fig2 = corner.corner(opt.samples)#,labels=opt.labels)
    plt.show()

    """ Figure out what it all means """