# This routine is written to reproduce the model used by Clio Gielen.
# The model is described by Gielen et al. (2008, A&A 490, 725) for post-AGB
# stars. It has also been used by Gielen et al. (2011, A&A 533, A99) to
# model and compare Galactic and Magellanic post-AGB stars.
# An erratum was issued (Gielen et al. 2010, A&A 515, C2), but this seems to
# be only for the code which contained a bug.
#
# According to Gielen et al. (2008) the model emission is given by:
#
# $F_\lambda \propto ( \sum_i \alpha_i \kappa_i ) \times ( \sum_j \beta_j B_{\lambda}(T_j) )$
#
# with
# kappa_i = mass absorption coefficient of dust component i
# alpha_i = the (mass?) fraction of dust component i
# B_lambda (T_j) = the Planck function at temperature T_j
# beta_j = the (mass?) fraction of dust at temperature T_j
#
# but inspecting their code and the paper in detail reveals that the model
# implemented actually reads:
# $F_\lambda \propto \sum_{j=1}^{\{2,3\}} (( \sum_{i=1}^8 \alpha_i \kappa_i ) \times \beta_j B_{\lambda} (T_j)) + \sum_{k=1}^{\{2,3\}} \beta_k B_{\lambda}(T_k)$ 
# with the left hand term giving rise to the features, and the right hand term
# contributing to a continuum only. The temperatures of the two terms are not
# the same, nor are the distribution over the two temperatures.
# The mass fractions in each of the dust components is the same for each of the
# 2 or 3 temperature components on the left hand side.
#
# Constraints on the parameters used by Gielen et al. 2008
# Only two temperature components: j=1 and j=2, between 100 and 1000 K, with
# their scaling factors: 2 * 2 free parameters for the features and
# 2 * 2 for the continuum term. 
# Four silicate species (amorphous and crystalline olivine and pyroxene),
# each with two dust sizes, and their relative fractions (7 free parameters).
#
# Thus, we arrive at a total of 2*2 * 2 +7 = 15 free parameters, as claimed
# by Gielen et al. 2008. 


import numpy as np
import ampere
from ampere.data import Spectrum, Photometry
from ampere.emceesearch import EmceeSearch
from ampere.PowerLawAGN import TripleBlackBodyDust
from ampere.extinction import CCMExtinctionLaw
import corner
import matplotlib.pyplot as plt
import os

from astropy import constants as const
from astropy import units as u
from astropy.io import ascii
from astropy.table import Table, Column
from spectres import spectres
import pyphot

if __name__=="__main__":
    """
    This script applies the model from Gielen et al. (2008)
    """

    """ Read in some data """
    
    """ Read in spectrum """
#    filename = './NGC6302/NGC6302_100.tab'
#    print(filename)
    #spec = Spectrum.fromFile(filename,format='User-Defined',filetype='text')
#    sws = ascii.read(filename,header_start=0,data_start=2)
#    unc = Column(0.05*np.abs(sws['Flux'].data),name='Uncertainty')
#    sws.add_column(unc,index=2)

#    spec = Spectrum(wavelength=sws['Wavelength'].data,value=sws['Flux'].data,uncertainty=sws['Uncertainty'].data,fluxUnits='Jy',bandUnits='um')
    
    """ Define the model """

    modwaves = 10**np.linspace(0.3,2.3, 2000) #2-200 um spectrum w/ 2000 data points

    Tmin = 100.  #minimum temperature allowed in fit
    Tmax = 1000. #maximum temperature allowed in fit
    
    """ Generate model """

    opacities = ['oliv_c_s.q',
    		 'oliv_c_b.q',
    		 'oliv_a_s.q',
                 'oliv_a_b.q',
                 'pyr_c_s.q',
                 'pyr_c_b.q',
                 'pyr_a_s.q',
                 'pyr_a_s.q']
    # opacities to be calculated separately. We consider 8 dust components:
    # olivine and pyroxene, either in crystalline or amorphous form, and in
    # the form of small or large grains. The q values are to be calculated
    # separately from n,k values.
    # Gielen et al. (2008) considered grain sizes of 0.1 and 2.0 micron.
    # The grains are non-spherical, and should be modelled using either
    # DHS or GRF. CDE does not allow for a distinction between grain sizes.
    # The optical properties are taken from:
    # Servoin & Pirou (1973), Dorschner et al. (1995),
    # Henning & Stognienko (1996), Jaeger et al. (1998).
    

    # list of opacity spectra to be used in modelling
    #relativeAbundances=np.array([0.01,0.01])#initial guess
    #nSpecies = len(opacities)-1

    mdl = DualBlackBodyDust(modwaves,
                          normWave = 1., sigmaNormWave = 1.,
                          opacityFileList=opacities,
                          lims = np.array([[0,1e6],
                                                             [-100,100],
                                                             [-10,10],
                                                             [0,np.inf]])
                            )
    
    """ Connect the dots: """
    """ Hook it up to an optimiser """

    opt = EmceeSearch(model = mdl, data = [spec], nwalkers = 100) #introspection is used inside the optimiser to determine the number of parameters that each part of the model has

    """ if you want, overload the default priors with a new function """
    def lnprior(self, inputs, **kwargs):
        return 0 #flat prior
    
    #model.lnprior = lnprior

    print(opt.npars,int(opt.npars))
    """ Run it """
    pos = [
           [
               45., 25., 0., 910., 1., 0.5, 1., 0.5, 1.
               #20 + np.random.randn() for i in range(int(opt.npars))
           ]
           + np.random.randn(int(opt.npars)) for j in range(opt.nwalkers)
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
    #print(opt.samples.shape)
    #print(np.max(opt.samples, axis=0))
    #print(np.min(opt.samples, axis=0))
    nneg=0

    ax.set_ylim(0.1,1000.0)
    ax.set_yscale('log')
    ax.set_ylabel('Flux (Jy)')
    ax.set_xlim(2.0,200.0)
    ax.set_xlabel('Wavelength (microns)')

    for i in range(0,opt.samples.shape[0],100):
        #print(opt.samples[i,:])
        if opt.samples[i,0] > 0.:
            opt.model(opt.samples[i,0],opt.samples[i,1],opt.samples[i,2],opt.samples[i,3],*opt.samples[i,4:6])
            ax.plot(modwaves,opt.model.modelFlux, '-', 'k', alpha=0.02) #irs.wavelength
        else:
            nneg += 1
            
    #    ax.plot(irs.wavelength, model(opt.samples[i,:]), '-', alpha = 0.1)
    print(nneg)
    #for i in range(0,len(spec.value)):
    ax.plot(spec.wavelength, spec.value, '-','b')
    fig2 = corner.corner(opt.samples)#,labels=opt.labels)
    plt.show()

    """ Figure out what it all means """
