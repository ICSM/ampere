from __future__ import print_function

import numpy as np
import emcee
from .basesearch import BaseSearch
from inspect import signature


class EmceeSearch(BaseSearch):
    """
    A class to use the emcee affine-invariant ensemble sampler
    """

    def __init__(self, nwalkers = None, model = None,
                 data = None, lnprior = None,
                 **kwargs):

        #self.npars=npars
        #self.nsamp=nsamp
        self.nwalker=nwalkers
        #self.burnin=burnin
        self.model=model
        self.dataSet = data
        if lnprior is not None:
            self.lnprior = lnprior
        ''' now do some introspection on the various bits of model to 
        understand how many parameters there are for each compponent '''
        sig = signature(model.__call__)
        self.nparsMod = len(sig.parameters) - 2 #Always subtract self and **kwargs from the parameters
        self.nparsData = np.zeros(len(self.dataSet))
        #self.npars = something # total number of parameters
        #self.nparsMod = something #number of parameters for the model
        #self.nparsData = [something for data in self.dataSet] #number of parameters to be passed into each set of data
        self.npars = self.nparsMod + np.sum(self.nparsData)
        ''' then set up the sampler '''
        self.sampler = emcee.EnsembleSampler(self.nwalkers, self.npars, self.lnprob, args=dataSet)

 
        
        #raise NotImplementedError()

    def __call__(self, guess=None, **kwargs):
        if guess is None:
            guess = [np.random.randn(npars) for i in range(self.nwalkers)]
        self.sampler.run_mcmc(guess, self.nsamp)
        self.allSamples = self.sampler.chain
        self.samples = self.sampler.chain[:, self.burnin, :].reshape((-1, self.npars))
        
        pass

    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()

    def lnprior(self, theta, **kwargs):
        """
        Simple uniform prior over all parameter space
        """
        return 0

    def optimise(self, nsamples = None, burnin = None, guess = None, **kwargs):
        if guess is None:
            guess = [np.random.randn(npars) for i in range(self.nwalkers)]
        self.nsamp = nsamples
        self.burnin = burnin
        self.sampler.run_mcmc(guess, nsamples)
        self.allSamples = self.sampler.chain
        self.samples = self.sampler.chain[:, self.burnin, :].reshape((-1, self.npars))
        pass


#    def lnlike(self, theta, **kwargs):
#        model = self.model(theta)
#        like = 0.
#        for data in self.dataSet:
#            """ do something for each bit of data """
#            deltaLike = something
#            like = like + deltaLike
#        return like #-0.5 * np.sum((((y - model)**2)/(yerr**2)))


