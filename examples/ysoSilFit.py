import numpy as np
import ampere
from ampere.data import Spectrum, Photometry
from ampere.emceesearch import EmceeSearch
from ampere.starScreen import PolynomialScreen
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
    
# S140 is an example of an SWS spectrum
# HOPS-68 is an example of a Spitzer spectrum
    dataDir = os.getcwd() + '/YSOsils/'
    hops68Data = 'cassis_yaaar_spcfw_20838656t.fits'
    s140Dats = '06401081_sws.fit'
    hops68irs = Spectrum.fromFile(dataDir+specFile,format='SPITZER-YAAAR')
    #irs = [chunk.selectWaves(low = 8., up=35.) for chunk in irs] # alternative to next two lines, doesn't work
    #irs[0].selectWaves(low = 8., up = 35.) #following Srinivasan et al. 2017
    #irs[1].selectWaves(low = 8., up = 35.) #following Srinivasan et al. 2017

    #libDir = '../ampere/'
    #libname = libDir + 'ampere_allfilters.hd5'
    #phot = Photometry.fromFile(dataDir+photFile, libName = libname)
    #phot.selectWaves(low = 35., interval = "right-open") #using only MIPS-70 and PACS, following Srinivasan et al. 2017
    #print(phot.mask)
    #print("hello")
    modwaves = 10**np.linspace(0.,1.9, 2000)

    model = PolynomialScreen(modwaves)
#    phot.reloadFilters(modwaves)
#    dataSet = [phot]          #use this line when using photometry
    dataSet = [s for s in irs] #comment out when using photometry
    for s in irs:             #include the next two lines when appending spectroscopy to photometry
        dataSet.append(s)
    for s in dataSet:
        print(s)

    opt = EmceeSearch(model = model, data = dataSet, nwalkers = 200)

    print(opt.npars)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in irs:
        ax.plot(i.wavelength, i.value, '-',color='blue')
    ax.plot(phot.wavelength, phot.value, 'o',color='blue')
    ax.set_ylim(0., 1.5*np.max([np.max(i.value) for i in dataSet]))
    fig.savefig("sed_test.png")

    model(-2.,0.3,-0.5,
                    -0.5,-0.5,-0.5,
                    -0.5,-0.5)
    modspec = model.modelFlux
    print(modspec)
    ax.plot(modwaves,modspec)
    plt.show() #this plots the spectrum and photometry plus the shape of the model SED using the input parameters
    #exit()
    
    pos = [
           [
               -2., 0.3, -0.78, -0.78, -0.78, -0.78, -0.78, -0.78, 1., 0.5, 1., 1., 0.5, 1. #spectroscopy included
               #-2., 0.3, -0.78, -0.78, -0.78, -0.78, -0.78, -0.78 # only photometry
               #20 + np.random.randn() for i in range(np.int(opt.npars))
           ]
           + np.random.randn(np.int(opt.npars))/1e3 for j in range(opt.nwalkers)
          ]
    #for i in range(len(pos)):
        #pos[i][0] = pos[i][0] / 1000.
        #pos[i][1] = pos[i][1] / 3. + 1.
        #pos[i][2:7] = pos[i][2:7] / 30. + 0.1
    print(pos[0])
    print(np.max(pos, axis=0))
    print(np.min(pos, axis=0))
    #print(np.min(newGuess,axis=0))
    #print(np.max(newGuess, axis=0))
    '''Probability space seems to be very complicated, so we're going to try a bit of a trick'''
    ''' First, we do a very short run with the initial guess that we specified before '''
    opt.optimise(nsamples = 20, burnin=0,guess=pos)
    repeat = True
    acrate=2.0
    lnprob = opt.sampler.lnprobability
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
        row_indices, col_indices = np.unravel_index(flat_indices, lnprob.shape)
        ''' put these into new short run as initial guess'''
        newGuess = opt.sampler.chain[row_indices, col_indices, :].reshape((-1, opt.npars))
        print(newGuess.shape)
        #print(newGuess)
        print(np.min(newGuess,axis=0))
        print(np.max(newGuess, axis=0))
        #opt.sampler.reset()
        if i > 1:
            opt.rebuildSampler(acceptRate=acrate)
        opt.optimise(nsamples = 50, burnin=0,guess=newGuess)
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

    fig = plt.figure()
    ax = fig.add_subplot(111)
    print(opt.samples.shape)
    print(np.max(opt.samples, axis=0))
    print(np.min(opt.samples, axis=0))
    nneg=0
    opt.postProcess()
    for i in range(0,opt.samples.shape[0],1000):
        #print(opt.samples[i,:])
        #if opt.samples[i,0] > 0.:
        opt.model(opt.samples[i,0],opt.samples[i,1],*opt.samples[i,2:8])
        ax.plot(modwaves,opt.model.modelFlux, '-', color='black', alpha=0.05) #irs.wavelength
        #else:
        #    nneg += 1
            
    #    ax.plot(irs.wavelength, model(opt.samples[i,:]), '-', alpha = 0.1)
    print(nneg)
    for i in irs:
        ax.plot(i.wavelength, i.value, '-',color='blue')
    ax.plot(phot.wavelength, phot.value, 'o',color='blue')
    ax.set_ylim(0., 1.5*np.max([np.max(i.value) for i in dataSet]))
    fig.savefig("seds.png")
    #fig2 = corner.corner(opt.samples)#,labels=opt.labels)
    #plt.show()


    ''' Now save the chain in case further analysis is interesting later '''
    opt.save('hops68_optimser_dump.pkl')
