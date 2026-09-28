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

    Parameters
    ----------
nwalkers : int, optional
The number of Markov Chain Monte Carlo (MCMC) walkers. The default is None.
nsteps : int, optional
The number of MCMC steps to take. The default is None.
model : object, optional
The model to use for the MCMC sampling. The default is None.
verbose : bool, optional
Whether to print progress statements during the MCMC sampling. The default is False.
data : list, optional
The data to use for the MCMC sampling. The default is None.
lnprior : function, optional
The log prior distribution to use for the MCMC sampling. The default is None.
vectorize : bool, optional
Whether to use a vectorized version of the log posterior distribution. The default is True.
parameter_labels : list, optional
The labels for the parameters being sampled. The default is None.
acceptRate : float, optional
The acceptance rate for the MCMC sampling. The default is 2.0.
moves : list, optional
The moves to use during the MCMC sampling. The default is None.
name : str, optional
The name for the MCMC sampler. The default is an empty string.
namestyle : str, optional
The style for the name of the MCMC sampler. The default is 'full'.
**kwargs : additional keyword arguments
Additional keyword arguments to pass to the super class.

    Generated with Chat-GPT
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
            self.sampler = emcee.EnsembleSampler(self.nwalkers, int(self.npars), self.lnprob_vector, a = self.acceptRate, moves=self.moves, vectorize=True)#, args=self.dataSet)
        else:
            self.sampler = emcee.EnsembleSampler(self.nwalkers, int(self.npars), self.lnprob, a = self.acceptRate, moves=self.moves)#, args=self.dataSet)


    def __call__(self, guess=None, **kwargs):
        if guess is None:
            guess = [np.random.randn(int(self.npars)) for i in range(self.nwalkers)]
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
        be used to update the model, data, sampler setup, or prior part-way through a run 

        Parameters
        ----------
        nwalkers : int, optional
        The number of walkers to use in the MCMC sampling. If not provided, the current number of walkers will be used.
        model : object, optional
        A new model to use in the MCMC sampling. If not provided, the current model will be used.
        data : object, optional
        New data to use in the MCMC sampling. If not provided, the current data will be used.
        lnprior : function, optional
        A new prior distribution to use in the MCMC sampling. If not provided, the current prior will be used.
        labels : list of str, optional
        A list of parameter labels to use in post-processing. If not provided, the 
        current labels will be used.
        **kwargs : additional keyword arguments
        Additional keyword arguments to pass to the sampler.

        Generated with Chat-GPT

        '''

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

        #    self.npars = int(self.nparsMod + np.sum(self.nparsData))
        
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
            self.sampler = emcee.EnsembleSampler(self.nwalkers, int(self.npars), self.lnprob, a = self.acceptRate, moves=self.moves, **kwargs)
        else:
            self.logger.info("No arguments passed, proceeding with sampler as is")

    def lnprior(self, theta, **kwargs):
        """
        We delegate the priors to the models and the data

        Parameters
        ----------
        theta : array_like
            The parameter values for which the prior probability is to be calculated.
        kwargs : optional
            Any additional keyword arguments required for the prior probability calculation.
        
        Returns
        -------
        lp : float
            The prior probability for the given parameter values.

        Generated with Chat-GPT
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
        """
        Optimise the model using the Markov Chain Monte Carlo (MCMC) method with emcee
        
        Parameters
        ----------
        nsamples : int, optional
        Number of samples to produce. The default is None.
        burnin : int, optional
        Number of samples to discard as burn-in. The default is None.
        guess : list or array-like, optional
        Initial guess for the model parameters. The default is None.
        preopt : bool, optional
        Whether to use a pre-optimization step to find an approximate maximum a posteriori (MAP) solution. The default is True.
        guessscale : float, optional
        Scaling factor for the initial guess when using pre-optimization. The default is 1e-3.
        noguess : bool, optional
        Whether to use initial guess for sampling. The default is False.
        progress : bool, optional
        Whether to display progress bar while sampling. The default is True.
        **kwargs : optional
        Additional keyword arguments to pass to scipy.optimize.minimize when using pre-optimization.
        
        Returns
        -------
        None
        
        Raises
        --------
        ValueError
        If the input guess does not conform to expectations.

        Generated with Chat-GPT
        """
        
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
                guess = [np.random.randn(int(self.npars)) for i in range(self.nwalkers)]
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
            guess = [newguess + guessscale*np.random.randn(int(self.npars)) for i in range(self.nwalkers)]
            
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

        Parameters
        ----------
        show : bool, optional
        Whether to show the plots, by default False
        textfile : str, optional
        Path to a text file to save a summary of the results, by default None
        **kwargs
        Additional keyword arguments to be passed to plot_posteriorpredictive()
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

    def plot_lnprob(self):
        """ Plot the log probability of each walker

        Plots the evolution of the log-probability of each walker throughout the sampling process. 
        A grey patch is added to indicate the burn-in period and a red line is added to indicate 
        10 times the autocorrelation time.

        """
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

    def summary(self):
        super().summary()
        self.maf = np.mean(self.sampler.acceptance_fraction)

    def print_summary(self, outfile=None):

        ''' First calculate some diagnostic information, like the autocorrelation time and acceptance fraction '''

        super().print_summary()
        self.logger.info("Mean Acceptance Fraction : %.5f", self.maf)
