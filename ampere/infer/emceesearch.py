from __future__ import print_function

import numpy as np
import emcee
import logging
from .basesearch import BaseSearch
from .mcmcsearch import EnsembleSampler
from inspect import signature
from ..data import Photometry, Spectrum
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from .mixins import SimpleMCMCPostProcessor, ScipyMinMixin

class EmceeSearch(EnsembleSampler,
                  SimpleMCMCPostProcessor, #):#,
                  ScipyMinMixin):
    """
    A class to use the emcee affine-invariant ensemble sampler
    """
    _inference_method = "Affine-invariance ensemble MCMC with emcee"

    def __init__(self, nwalkers = None, model = None, verbose = False,
                 data = None, lnprior = None, vectorize = True, 
                 parameter_labels = None, acceptRate = 2.0, moves=None,
                 name='', namestyle="full",
                 **kwargs):

        #Call super with everything required:
        super().__init__(nwalkers = nwalkers, model = model, data= data, verbose = verbose,
                         parameter_labels = parameter_labels, name = name,
                         namestyle=namestyle, **kwargs)
        self.moves = moves
        self.acceptRate = acceptRate
        ''' then set up the sampler '''
        if vectorize:
            self.logger.info("Using vectorised posterior")
            #self.lnprob = self.lnprob_vector
            self.sampler = emcee.EnsembleSampler(self.nwalkers, np.int(self.npars), self.lnprob_vector, a = self.acceptRate, moves=self.moves, vectorize=True)#, args=self.dataSet)
        else:
            self.sampler = emcee.EnsembleSampler(self.nwalkers, np.int(self.npars), self.lnprob, a = self.acceptRate, moves=self.moves)#, args=self.dataSet)


    def __call__(self, guess=None, **kwargs):
        if guess is None:
            guess = [np.random.randn(np.int(self.npars)) for i in range(self.nwalkers)]
        self.sampler.run_mcmc(guess, self.nsamp)
        self.allSamples = self.sampler.chain
        self.samples = self.sampler.chain[:, self.burnin, :].reshape((-1, self.npars))
        pass

    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()

    def rebuildSampler(self, nwalkers = None, model = None,
                       data = None, lnprior = None,
                       labels = None, acceptRate = None,
                       moves=None,
                       **kwargs):
        ''' This method will replace parts of the optimiser object if it is passed them. It can 
        be used to update the model, data, sampler setup, or prior part-way through a run '''

        #if np.any([model, data]):
        #    if model is not None:
        #        self.model=model
        #        try:
        #            self.nparsMod = self.model.npars
        #        except:
        #            sig = signature(model.__call__)
        #            self.nparsMod = len(sig.parameters) - 1 #Always subtract **kwargs from the parameters, but don't need to worry about self once it is bound to an instance
        #    if data is not None:
        #        self.dataSet = data
        #        self.nparsData = [data.npars for data in self.dataSet] #number of parameters to be passed into each set of data

        #    self.npars = np.int(self.nparsMod + np.sum(self.nparsData))
        
        #if lnprior is not None: 
       #    self.lnprior = lnprior
        ##print(self.npars, self.nparsMod, self.nparsData)
        #if nwalkers is not None:
        #        self.nwalkers=nwalkers
        ''' then set up the sampler '''
        
        if np.any([nwalkers, acceptRate, data, model, lnprior, labels, moves, kwargs]):
            super().rebuildSampler(nwalkers = nwalkers, model = model, data = data, lnprior = lnprior, labels = labels, **kwargs)
            if moves:
                self.moves = moves
            if acceptRate:
                self.acceptRate = acceptRate
            ''' if anything has been changed, update the sampler and re-intialise it '''
            ''' first destroy the sampler '''
            self.sampler=None
            ''' now rebuild it '''
            self.sampler = emcee.EnsembleSampler(self.nwalkers, np.int(self.npars), self.lnprob, a = self.acceptRate, moves=self.moves, **kwargs)
        else:
            self.logger.info("No arguments passed, proceeding with sampler as is")

    def lnprior(self, theta, **kwargs):
        """
        We delegate the priors to the models and the data
        """
        #print(theta)
        lp = self.model.lnprior(theta[:self.nparsMod])
        i = self.nparsMod
        for data in self.dataSet:
            #print(theta,theta[:i],theta[i:i+data.npars],i,data.npars)
            #print(lp)
            lp+= data.lnprior(theta[i:i+data.npars])
            i+=data.npars
        #print(lp)
        return lp #return 0

    #functionality moved to mixins.ScipyMinMixin.minimize()
    #def preopt(self, start, **kwargs):
    #    neglnprob = lambda *args: -self.lnprob(*args)
    #    from scipy.optimize import minimize
    #    self.logger.info("Using scipy.minimize to find approximate MAP solution")
    #    self.logger.info("starting from position: %s",start)
    #    solution = minimize(neglnprob, start)
    #    guess = [p for p in solution.x]
    #    self.logger.info("Minimization complete, final position: %s ",guess)
    #    return guess
        

    def optimise(self, nsamples = None, burnin = None, guess = None,
                 preopt = True, guessscale = 1e-3, noguess=False, progress=True, **kwargs):
        from collections.abc import Sequence, Iterable
        self.logger.info("Preparing to sample")
        if guess == 'None': # and not noguess:
            self.logger.info("Setting initial guess randomly")
            try: #Attempting to use the ptform
                self.logger.info("No guess specified, attempting to draw from prior transform")
                rng = np.random.default_rng() #.uniform(-1,0,1000)
                guess = [self.prior_transform(rng.uniform(size=self.npars)) for i in range(self.nwalkers)]
                #print(guess)
            except AttributeError: #not all components have a ptform, so we'll do this the simple way
                self.logger.info("Drawing from prior transform failed, drawing randomly")
                guess = [np.random.randn(np.int(self.npars)) for i in range(self.nwalkers)]
        #print('dimensions of guess = ', np.shape(guess))
        #print('nsamples = ', nsamples)
        #print('burnin = ', burnin)
        #print("np.max(guess, axis=0) = ", np.max(guess, axis=0))
        #print("np.min(guess, axis=0) = ", np.min(guess, axis=0))
        #print(self.dataSet)
        #print(self.model.lims)

        if preopt:
            self.logger.info("Searching for approximate MAP solution as starting point for MCMC")
            self.logger.info("Selecting start point for scipy.minimize")
            if isinstance(guess, Sequence):
                if not isinstance(guess[0], (Sequence, Iterable)):#Only a single entry
                    self.logger.debug("Only one entry in guess")
                    self.logger.debug(guess[0])
                    start = guess
                else: #guess contains many entries, randomly select one
                    self.logger.debug("Multiple entries in guess, selecting at random!")
                    rng = np.random.default_rng()
                    i = rng.integers(0, len(guess))
                    start = guess[i]
            else:
                #Guess is not appropriate, need to do something here!
                raise ValueError("guess does not conform to expectations")
            newguess = self.preopt(start)
            guess = [newguess + guessscale*np.random.randn(np.int(self.npars)) for i in range(self.nwalkers)]
            
        self.nsamp = nsamples
        self.burnin = burnin
        self.logger.info("Starting to sample: ")
        self.logger.info("Each walker will produce %d samples, of which the first %d will be discarded as burn-in", self.nsamp, self.burnin)
        if noguess:
            self.sampler.run_mcmc(nsamples, progress=progress)
        else:
            self.sampler.run_mcmc(guess, nsamples, progress=progress)
        self.allSamples = self.sampler.chain
        #print('do we get here (no): ',np.max(self.allSamples), np.min(self.allSamples))
        self.samples = self.sampler.chain[:, self.burnin:, :].reshape((-1, self.npars))
        #print(np.max(self.samples), np.min(self.samples))
        #pass



    def postProcess(self, show=False, textfile=None, **kwargs):
        ''' 
        A method to post-process the sampler results 
        '''

        #Compute things like autocorrelation time.
        #ESS?
        #Acceptance fraction?
        self.print_summary() #outfile=textfile)

        ''' Now produce some diagnostic plots '''

        ''' A plot of the walkers' sampling history '''
        self.plot_trace()
        
        ''' And a corner plot of the post-burnin results '''
        self.plot_corner()

        ''' plot the evolution of the likelihood function '''
        self.plot_lnprob()

        '''finally, we want to look at the correlation/covariance matrices for the spectra, if any '''
        self.plot_covmats()


        self.plot_posteriorpredictive(**kwargs)
        
    #def plot_walkers(self):
    #    #MUST USE autocorrelation time and burnin info on plots!
    #    #Should probably have quiet set 'False' to pick up too short emcee runs
    #    tauto = self.sampler.get_autocorr_time(quiet=True)

    #    fig, axes = plt.subplots(self.npars, 1, sharex=True, figsize=(8, 9))
    #    for i in range(self.npars):
    #        axes[i].plot(self.sampler.chain[:, :, i].T, color="k", alpha=0.4) #,label=self.parLabels[i])
    #        axes[i].set_xlim(0, self.nsamp)
    #        ylims = axes[i].get_ylim()
    #        xlims = axes[i].get_xlim()
    #        axes[i].plot([10*tauto[i],10*tauto[i]],[ylims[0],ylims[1]],color="red",marker="",linestyle="-")#,label=r"10$\times t_{\rm autocorr}$")
    #        axes[i].add_patch(Polygon([[xlims[0], ylims[0]], [xlims[0], ylims[1]], [self.burnin, ylims[1]], [self.burnin, ylims[0]]], closed=True,
    #                  fill=True, color='darkgrey'))
    #        axes[i].set_xlabel("Steps")
    #        axes[i].set_ylabel(self.parLabels[i])
    #        axes[i].label_outer()
    #        axes[i].set_xlim(0, self.nsamp)

    #    #plt.legend(fontsize="small")
    #    plt.tight_layout()

    #        #    axes[0].yaxis.set_major_locator(MaxNLocator(5))
    #        #axes[0].axhline(m_true, color="#888888", lw=2)
    #        #axes[i].set_ylabel("$m$")
    #    
    #
    #    #fig.tight_layout(h_pad=0.0)
    #    fig.savefig("walkers.png")
    
    #def plot_corner(self):
    #    import corner
    #    fig2 = corner.corner(self.samples,labels=self.parLabels)
    #    fig2.savefig("corner.png")

    def plot_lnprob(self):
        #USE autocorrelation time and burnin info on plots?
        tauto = self.sampler.get_autocorr_time(quiet=True)

        fig = plt.figure()
        axes = fig.add_subplot(111) #plt.subplots(1, 1, sharex=True, figsize=(6, 8))
        for i in range(self.nwalkers):
            axes.plot(self.sampler.lnprobability[i,:])
        try:
            axes.set_ylim(np.min(self.sampler.lnprobability),np.max(self.sampler.lnprobability))
        except ValueError:
            try:
                axes.set_ylim(-2000.,np.max(self.sampler.lnprobability))
            except ValueError:
                axes.set_ylim(-2000.,2000)

        ylims = axes.get_ylim()
        xlims = axes.get_xlim()
        axes.add_patch(Polygon([[xlims[0], ylims[0]], [xlims[0], ylims[1]], [self.burnin, ylims[1]], [self.burnin, ylims[0]]], closed=True,
                      fill=True, color='grey'))
        axes.plot([10*tauto,10*tauto],[ylims[0],ylims[1]],color="red",marker="",linestyle="-",label=r"10$\times t_{\rm autocorr}$")
        fig.savefig(self.name+"_lnprob.png")        
        
    #def plot_covmats(self):
    #    istart = self.nparsMod
    #    for i, d in enumerate(self.dataSet):
    #        if isinstance(d, Photometry):
    #            continue
    #        elif isinstance(d, Spectrum):
    #            fig, ax = plt.subplots(1,1)
    #            #for d in self.dataSet[1:]:
    #            #d=self.dataSet[1]
    #            #print("Using these parameters for the covariance matrix:")
    #            #print(self.parLabels[istart+1], self.parLabels[istart+2])
    #            d.cov([self.res[istart+1][0],self.res[istart+2][0]])
    #            #ax0.set_title('Covariance matrix')
    #            im = ax.imshow(np.log10(d.covMat))
    #            istart+=d.npars
    #            fig.savefig("covMat_"+str(i)+".png")

    #def plot_posteriorpredictive(self, n_post_samples = 1000,plotfile="posteriorpredictive.png", logx = False, logy = False, alpha = 0.1,**kwargs):
    #    '''
    #    Function to plot the data and samples drawn from the posterior of the models and data.
    #    '''
    #    fig,axes = plt.subplots(1,1,figsize=(8,6))
    #    axes.set_xlabel(r"Wavelength ($\mu$m)")
    #    axes.set_ylabel(r"Flux density (mJy)")

    #    #observations
    #    for d in self.dataSet:
    #        if isinstance(d,(Photometry,Spectrum)):
    #            d.plot(ax = axes)
    #    for s in self.samples[np.random.randint(len(self.samples), size=n_post_samples)]:
    #        try:
    #            self.model(*s[:self.nparsMod])
    #            axes.plot(self.model.wavelength,self.model.modelFlux, '-', color='k', alpha=alpha,label='Samples', zorder = 0)
    #        except ValueError:
    #            pass
    #        i = self.nparsMod
    #        for d in self.dataSet:
    #            if isinstance(d,(Photometry,Spectrum)):
    #                if d._hasNoiseModel:
    #                    d.plotRealisation(s[i:i+d.npars], ax=axes)
    #                i+= d.npars


    #    #best fit model
    #    try:
    #        self.model(*self.bestPars[:self.nparsMod])
    #        axes.plot(self.model.wavelength,self.model.modelFlux, '-', color='k', alpha=1.0,label='MAP', zorder=8)
    #    except ValueError:
    #        print("Error in MAP solution \n Skipping MAP in plot")

    #    #These plots end up with too many labels for the legend, so we clobber the label information so that only one of each one is plotted
    #    handles, labels = plt.gca().get_legend_handles_labels()
    #    # labels will be the keys of the dict, handles will be values
    #    temp = {k:v for k,v in zip(labels, handles)}
    #    plt.legend(temp.values(), temp.keys(), loc='best')

    #    #plt.legend()
    #    plt.tight_layout()
    #    fig.savefig(plotfile,dpi=200)
    #    plt.close(fig)
    #    plt.clf()

    def summary(self):
        super().summary()
        self.maf = np.mean(self.sampler.acceptance_fraction)

    def print_summary(self, outfile=None):

        ''' First calculate some diagnostic information, like the autocorrelation time and acceptance fraction '''

        super().print_summary()
        self.logger.info("Mean Acceptance Fraction : %.5f", self.maf)
        
       # print("Mean Acceptance Fraction : {0:.5f}",np.mean(self.sampler.acceptance_fraction))
       # #try:
       # tauto = self.sampler.get_autocorr_time(quiet=True)
       # print("Mean Autocorrelation Time: {0:.5f}",np.mean(tauto))
       # 
       # print("Autocorrelation times for each parameter:")
       # #tauto = self.sampler.get_autocorr_time()
       # for i in range(self.npars):
       #     print("{0}  = {1:.0f} steps".format(self.parLabels[i],tauto[i]))
       # #except emcee.autocorr.AutocorrError as e:
       # #    print(repr(e))
       # #    print("Skipping autocorrelation time")
       #     

       # #Compute things like autocorrelation time.
       # #ESS?
       # #Acceptance fraction?

       # '''Now find the median and 68% interval '''
       # self.res=[]
       # print("Median and confidence intervals for parameters in order:")
       # for i in range(self.npars):
       #     a = np.percentile(self.samples[:,i], [16, 50, 84])
       #     self.res.append([a[1], a[2]-a[1], a[1]-a[0]])
       #     #res.append((lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
       #     #                     *zip(np.percentile(self.samples[:,i], [16, 50, 84])))#,
       #                                             #axis=0)))
       #     #      )
       #     print("{0}  = {1[0]:.5f} + {1[1]:.5f} - {1[2]:.5f}".format(self.parLabels[i],self.res[i])) #make this look prettier
            

        #''' Then check what the "best fit" was '''
        #print("Range of lnprob values in model from sampler: {0:.5f} to {1:.5f}",np.min(self.sampler.lnprobability),np.max(self.sampler.lnprobability))
        #row_ind, col_ind = np.unravel_index(np.argmax(self.sampler.lnprobability), self.sampler.lnprobability.shape)
        #self.bestPars = self.sampler.chain[row_ind, col_ind, :]

        #print("MAP Solution: ")
        #for i in range(self.npars):
        #    print("{0}  = {1:.5f}".format(self.parLabels[i],self.bestPars[i])) #

        #if outfile is not None:
        #    with open(outfile, 'w') as f:
        #        f.write("Posterior means and 1-sigma confidence intervals of the parameters marginalising over all other parameters: \n ")
        #        for i in range(self.npars):
        #            f.write("{0}  = {1[0]:.5f} + {1[1]:.5f} - {1[2]:.5f}".format(self.parLabels[i],self.res[i])
        #            )

        #        f.write("\n")
        #        f.write("MAP Solution: \n")
        #        for i in range(self.npars):
        #            f.write("{0}  = {1:.5f}".format(self.parLabels[i],self.bestPars[i])) #
        #            #f.write("with Posterior probability ln(P*) = {0:.5f}\n".format(self.results.logwt[-1]))
        #            #f.write("and likelihood ln(L*) = {0:.5f}\n".format(self.results.logl[-1]))
