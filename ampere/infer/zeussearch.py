from __future__ import print_function

import numpy as np
import zeus
import logging
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from .basesearch import BaseSearch
from .mcmcsearch import EnsembleSampler
from inspect import signature
from ..data import Photometry, Spectrum
from .mixins import SimpleMCMCPostProcessor, ScipyMinMixin

class ZeusSearch(EnsembleSampler,
                  SimpleMCMCPostProcessor, #):#,
                  ScipyMinMixin):
    """
    A class to use the Zeus ensemble slice sampler
    """

    def __init__(self, nwalkers = None, model = None, verbose = False,
                 data = None, lnprior = None, vectorize = True,
                 parameter_labels = None, moves = None,
                 **kwargs):

        super().__init__(nwalkers = nwalkers, model = model, data= data,
                         verbose = verbose,
                         parameter_labels = parameter_labels, **kwargs)
        ''' then set up the sampler '''
        self.moves = moves
        ''' then set up the sampler '''
        if vectorize:
            logging.info("Using vectorised posterior")
            #self.lnprob = self.lnprob_vector
            self.sampler = zeus.EnsembleSampler(self.nwalkers, np.int(self.npars), self.lnprob_vector, moves=self.moves, vectorize=True)#, args=self.dataSet)
        else:
            self.sampler = zeus.EnsembleSampler(self.nwalkers, np.int(self.npars), self.lnprob, moves=self.moves)#, args=self.dataSet)


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
                       labels = None,
                       **kwargs):
        ''' This method will replace parts of the optimiser object if it is passed them. It can 
        be used to update the model, data, sampler setup, or prior part-way through a run '''

        if np.any([nwalkers, data, model, lnprior, labels, moves, kwargs]):
            ''' if anything has been changed, update the sampler and re-intialise it '''
            super().rebuildSampler(nwalkers = nwalkers, model = model, data = data, lnprior = lnprior, labels = labels, **kwargs)
            if moves:
                self.moves = moves
            ''' first destroy the sampler '''
            self.sampler=None
            ''' now rebuild it '''
            self.sampler = zeus.EnsembleSampler(self.nwalkers, np.int(self.npars), self.lnprob, moves = self.moves, **kwargs)

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

    def optimise(self, nsamples = None, burnin = None, guess = None,
                 preopt = True, guessscale = 1e-3, noguess=False, progress=True, **kwargs):
        from collections.abc import Sequence, Iterable
        logging.info("Preparing to sample")
        if guess == 'None': # and not noguess:
            logging.info("Setting initial guess randomly")
            try: #Attempting to use the ptform
                logging.info("No guess specified, attempting to draw from prior transform")
                rng = np.random.default_rng() #.uniform(-1,0,1000)
                guess = [self.prior_transform(rng.uniform(size=self.npars)) for i in range(self.nwalkers)]
                print(guess)
            except AttributeError: #not all components have a ptform, so we'll do this the simple way
                logging.info("Drawing from prior transform failed, drawing randomly")
                guess = [np.random.randn(np.int(self.npars)) for i in range(self.nwalkers)]

        if preopt:
            logging.info("Searching for approximate MAP solution as starting point for MCMC")
            logging.info("Selecting start point for scipy.minimize")
            if isinstance(guess, Sequence):
                if not isinstance(guess[0], (Sequence, Iterable)):#Only a single entry
                    logging.debug("Only one entry in guess")
                    logging.debug(guess[0])
                    start = guess
                else: #guess contains many entries, randomly select one
                    logging.debug("Multiple entries in guess, selecting at random!")
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
        logging.info("Starting to sample: ")
        logging.info("Each walker will produce %d samples, of which the first %d will be discarded as burn-in", self.nsamp, self.burnin)
        if noguess:
            self.sampler.run_mcmc(nsamples, progress=progress)
        else:
            self.sampler.run_mcmc(guess, nsamples, progress=progress)
        self.allSamples = self.sampler.chain
        self.samples = self.sampler.chain[self.burnin:, :,  :].reshape((-1, self.npars))


    def postProcess(self, show=False, **kwargs):
        ''' 
        A method to post-process the sampler results 
        '''

        '''First find the median and 68% interval '''
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
        
    #Now we need to overload some of the SimpleMCMCPostProcessor methods
    def get_map(self, **kwargs):
        self.lnprobability = self.sampler.get_log_prob()
        row_ind, col_ind = np.unravel_index(np.argmax(self.lnprobability), self.lnprobability.shape)
        self.bestPars = self.sampler.chain[row_ind, col_ind, :]
        return self.bestPars

    def get_autocorr(self, chain = None, quiet=True, **kwargs):
        self.tauto = self.sampler.act
        return self.tauto

    def plot_trace(self, **kwargs):
        try:
            tauto = self.tauto
        except AttributeError:
            tauto = self.get_autocorr()
        fig, axes = plt.subplots(self.npars, 1, sharex=True, figsize=(8, 9))
        for i in range(self.npars):
            axes[i].plot(self.sampler.chain[:, :, i], color="k", alpha=0.4) #No need to transpose the chain here, because zeus has steps in the first axis, not second (like emcee)
            axes[i].set_xlim(0, self.nsamp)
            ylims = axes[i].get_ylim()
            xlims = axes[i].get_xlim()
            axes[i].plot([10*tauto[i],10*tauto[i]],[ylims[0],ylims[1]],color="red",marker="",linestyle="-")#,label=r"10$\times t_{\rm autocorr}$")
            axes[i].add_patch(Polygon([[xlims[0], ylims[0]], [xlims[0], ylims[1]], [self.burnin, ylims[1]], [self.burnin, ylims[0]]], closed=True,
                      fill=True, color='darkgrey'))
            axes[i].set_xlabel("Steps")
            axes[i].set_ylabel(self.parLabels[i])
            axes[i].label_outer()
            axes[i].set_xlim(0, self.nsamp)

        plt.tight_layout()
        fig.savefig("walkers.png")

    def plot_lnprob(self):
        #USE autocorrelation time and burnin info on plots?
        tauto = self.tauto #self.sampler.get_autocorr_time(quiet=True)

        fig = plt.figure()
        axes = fig.add_subplot(111) #plt.subplots(1, 1, sharex=True, figsize=(6, 8))
        for i in range(self.nwalkers):
            axes.plot(self.lnprobability[:, i])
        try:
            axes.set_ylim(np.min(self.lnprobability),np.max(self.lnprobability))
        except ValueError:
            try:
                axes.set_ylim(-2000.,np.max(self.lnprobability))
            except ValueError:
                axes.set_ylim(-2000.,2000)

        ylims = axes.get_ylim()
        xlims = axes.get_xlim()
        axes.add_patch(Polygon([[xlims[0], ylims[0]], [xlims[0], ylims[1]], [self.burnin, ylims[1]], [self.burnin, ylims[0]]], closed=True,
                      fill=True, color='grey'))
        axes.plot([10*tauto,10*tauto],[ylims[0],ylims[1]],color="red",marker="",linestyle="-",label=r"10$\times t_{\rm autocorr}$")
        fig.savefig("lnprob.png")     


