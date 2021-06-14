import numpy as np
import sys
sys.path.insert(1, '/home/zeegers/git_ampere/ampere/')
import ampere
from ampere.data import Spectrum, Photometry
from ampere.emceesearch import EmceeSearch
from ampere.PowerLawAGN import PowerLawAGN, SingleModifiedBlackBody
from ampere.PowerLawAGN import PowerLawAGN, PowerLawAGN
from ampere.starScreen import PolynomialSource
#from extinction import CCMExtinctionLaw
#from extinction import apply, fitzpatrick99
import corner
import matplotlib.pyplot as plt
import os

from astropy import constants as const
from astropy import units as u
from spectres import spectres
import pyphot
from astropy.io.votable import parse_single_table

import pdb

if __name__=="__main__":
    """
    This script should contain
    """
    #ardila_2008=readfits('/home/zeegers/small_silicates_crystal/small_silicates/infrared_spectra/cygob2_12/cassis_yaaar_spcfw_27570176t.fits',ardila_header) 
    #evans_2004_HR=readfits('/home/zeegers/small_silicates_crystal/small_silicates/infrared_spectra/cygob2_12/cassis_yaaar_optdiff_9834496.fits',evans_header) 
    #evans_2004_LR=readfits('/home/zeegers/small_silicates_crystal/small_silicates/infrared_spectra/cygob2_12/cassis_yaaar_spcfw_9834496t.fits',evans_header2)

    """ Read in some data """
    filename = ampere.__file__.strip('__init__.py')+'Testdata/cassis_yaaar_spcfw_27570176t.fits' # path to data of cyg ob 2 12
    print(filename)
    irs = Spectrum.fromFile(filename,format='SPITZER-YAAAR')
    #irs[0].selectWaves(low = 5.4, up = 22.) #following Srinivasan et al. 2017, check if limits are ok
    #irs[1].selectWaves(low = 5.4, up = 22.) #following Srinivasan et al. 2017, check if limiets are ok
    
    #for s in irs:
    #    print(s)
    #exit()
    #If we have more than one dataset we need multiple data objects
    #print(irs.wavelength)

    filename2 = ampere.__file__.strip('__init__.py')+'Testdata/cassis_yaaar_spcfw_9834496t.fits' # path to data of cyg ob 2 12
    
    irs2 = Spectrum.fromFile(filename2,format='SPITZER-YAAAR')
    
    modwaves = 10**np.linspace(0.7,1.6, 1000) #setting up a wavelength grid from 5.25 to 37.4 um. 

    model = PolynomialSource(modwaves)
    dataSet = [s for s in irs] 
    
    for s in irs2:             #include the next two lines when appending spectroscopy to photometry
        dataSet.append(s)
    

    model = PolynomialSource(modwaves)
                                    
                                                                        
    print(model.npars)                                
    
    """ and the extinction law (if necessary) (could be moved inside model as only some models will need this) """
    #ext = CCMExtinctionLaw()

    #model.extinction = ext #combination of multiple different models is something we need to consider how to implement. 

    """ Connect the dots: """
    """ Hook it up to an optimiser """
    opt = EmceeSearch(model = model, data = dataSet, nwalkers = 50) #introspection is used inside the optimiser to determine the number of parameters that each part of the model has
    
    print(opt.npars)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in irs:
        ax.loglog(i.wavelength, i.value, '-',color='blue')
    for j in irs2:
        ax.loglog(j.wavelength, j.value, '-',color='blue')        
    #ax.set_ylim(0., 1.5*np.max([np.max(i.value) for i in dataSet]))
    #fig.savefig("sed_test.png")

    model(0.0,-1.2,2.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1)
    
    modspec = model.modelFlux
    print(modspec)
    ax.plot(modwaves,modspec, color="red")
    plt.show() #this plots the spectrum and photometry plus the shape of the model SED using the input parameters


    """ if you want, overload the default priors with a new function """
#    def lnprior(self, inputs, **kwargs):
#        return 0 #flat prior
    
    #model.lnprior = lnprior
    
    #20 + np.random.randn() for i in range(np.int(opt.npars))
    
    #pdb.set_trace()
    print(opt.npars,np.int(opt.npars))
    print("hier ben ik geweest")
    """ Run it """
          
    pos = [
           [
               1.0, -1.2, 2.0, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 1., 0.5, 1., 1., 0.5, 1., 1., 0.5, 1. # the last six parameters represent the three parameters for the noise model, for the two chunks in the data set.
           ]
           + np.random.randn(int(opt.npars))/100 for j in range(opt.nwalkers)
          ]
           
    print(pos[0])
    print(np.max(pos, axis=0))
    print(np.min(pos, axis=0))
    #pdb.set_trace()
    
 #**************************************************************************************************8   
    
    '''Probability space seems to be very complicated, so we're going to try a bit of a trick'''
    ''' First, we do a very short run with the initial guess that we specified before '''
    opt.optimise(nsamples = 60, burnin=0, guess=pos)
    plt.hist(opt.samples[np.isfinite(opt.samples)].flatten(), bins =100)
    plt.show()
    
    repeat = True
    acrate=2.0
    lnprob = opt.sampler.lnprobability
    print(np.max(lnprob), np.min(lnprob))
    ''' Now we extract the best 200 samples from the chain based on their lnprob, and use those as initial guesses for another short run '''
    ''' After this, we check the acceptance fraction, and if it's at least 0.2, we go to a production run, if not, we loop back to extracting the best 200 samples and repeat the process '''
    i=0
    while repeat:
        ''' Find 200 best samples from previous run'''
        i+=1
        print("starting loop iteration ", i, "with a = ",acrate)
        print(lnprob.shape)
        n = opt.nwalkers
        flat_indices = np.argpartition(lnprob.ravel(), -n-1)[-n:]
        print(np.max(lnprob), np.min(lnprob))
        plt.hist(lnprob[np.isfinite(lnprob)].flatten(), bins =100)
        plt.show()
        
        row_indices, col_indices = np.unravel_index(flat_indices, lnprob.shape)
        ''' put these into new short run as initial guess'''
        print('k = ',row_indices)
        print('iterations = ', col_indices)
        newGuess = opt.sampler.chain[row_indices, col_indices, :].reshape((-1, opt.npars))
        print(newGuess.shape)
        #print(newGuess)
        print(np.min(newGuess,axis=0))
        print(np.max(newGuess, axis=0))
        #opt.sampler.reset()
        if i > 1:
            opt.rebuildSampler(acceptRate=acrate)
        print(newGuess)
        opt.optimise(nsamples = 30, burnin=0,guess=newGuess)
        lnprob = opt.sampler.lnprobability
        #b += 50
        
        ''' Check acceptance '''
        acceptanceFraction = opt.sampler.acceptance_fraction
        acmin = np.min(acceptanceFraction)
        acmax = np.max(acceptanceFraction)
        print(acmin, acmax)
        #acmid = np.median(acceptanceFraction)

        if  acmin > 0.10: # and acmax < 0.40:
            repeat = False
        else:
            acrate *= 1.05
        #if afrac > 0.2:
        #    repeat = False
    print("exiting loop, moving to production run with a = ",acrate)
    #exit()
    opt.optimise(nsamples = 500, burnin=100, guess = opt.sampler.chain[:, -1, :]) #burnin should discard all steps taken before exiting the loop
    
    a = 1. - np.sum(10**opt.samples[2:8],axis=0)
    b = np.percentile(a, [16, 50, 84])
    print(a.shape)
    print(b[1], b[2]-b[1], b[1]-b[0])

#***********************************************************

    
    """ Sort out the results """
    """ First produce some terminal output for the numbers and confidence intervals """

    """ Then produce a couple of plots """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    print(opt.samples.shape)
    print(np.max(opt.samples, axis=0))
    print(np.min(opt.samples, axis=0))
    nneg=0
    print("hier ben ik geweest #2")
    for i in range(0,opt.samples.shape[0],100):
        #print(opt.samples[i,:])
        if opt.samples[i,0] > 0.:
            opt.model(opt.samples[i,0],opt.samples[i,1],opt.samples[i,2],opt.samples[i,3:9])
            ax.plot(modwaves,opt.model.modelFlux, '.', color='k', alpha=0.02) #irs.wavelength
        else:
            nneg += 1
            
    #    ax.plot(irs.wavelength, model(opt.samples[i,:]), '-', alpha = 0.1)
    print(nneg)
    for i in irs:
        ax.plot(i.wavelength, i.value, '.',color='b')
    for i in irs2:
        ax.plot(i.wavelength, i.value, '.',color='b')
    fig2 = corner.corner(opt.samples)#,labels=opt.labels)
    plt.show()

    """ Figure out what it all means """
