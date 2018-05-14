from __future__ import print_function

import numpy as np
import emcee
from basesearch import BaseSearch


class EmceeSearch(BaseSearch):
    """
    A class to use the emcee affine-invariant ensemble sampler
    """

    def __init__(self, npars, nsamp,
                 nwalkers, burnin, model,
                 dataSet, lnprior,
                 **kwargs):

        self.npars=npars
        self.nsamp=nsamp
        self.nwalker=nwalkers
        self.burnin=burnin
        self.model=model
        self.dataSet = dataSet
        if lnprior is not None:
            self.lnprior = lnprior
        self.sampler = emcee.EnsembleSampler(nwalkers, npars, self.lnprob, args=dataSet)
        
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

#    def lnlike(self, theta, **kwargs):
#        model = self.model(theta)
#        like = 0.
#        for data in self.dataSet:
#            """ do something for each bit of data """
#            deltaLike = something
#            like = like + deltaLike
#        return like #-0.5 * np.sum((((y - model)**2)/(yerr**2)))


