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
                 labels = None,
                 **kwargs):

        #self.npars=npars
        #self.nsamp=nsamp
        self.nwalkers=nwalkers
        #self.burnin=burnin
        self.model=model
        self.dataSet = data
        if lnprior is not None:
            self.lnprior = lnprior
        ''' now do some introspection on the various bits of model to 
        understand how many parameters there are for each compponent '''
        try:
            self.nparsMod = self.model.npars
        except:
            sig = signature(model.__call__)
            self.nparsMod = len(sig.parameters) - 1 #Always subtract **kwargs from the parameters, but don't need to worry about self once it is bound to an instance
        #print(self.nparsMod, len(sig.parameters), sig.parameters)
        #self.nparsData = #np.zeros(len(self.dataSet))
        #self.npars = something # total number of parameters
        #self.nparsMod = something #number of parameters for the model
        self.nparsData = [data.npars for data in self.dataSet] #number of parameters to be passed into each set of data
        self.npars = np.int(self.nparsMod + np.sum(self.nparsData))
        print(self.npars, self.nparsMod, self.nparsData)
        #exit()
        ''' then set up the sampler '''
        self.sampler = emcee.EnsembleSampler(self.nwalkers, np.int(self.npars), self.lnprob)#, args=self.dataSet)

 
        
        #raise NotImplementedError()

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

    def optimise(self, nsamples = None, burnin = None, guess = None, **kwargs):
        if guess is None:
            guess = [np.random.randn(np.int(self.npars)) for i in range(self.nwalkers)]
        self.nsamp = nsamples
        self.burnin = burnin
        self.sampler.run_mcmc(guess, nsamples)
        self.allSamples = self.sampler.chain
        self.samples = self.sampler.chain[:, self.burnin:, :].reshape((-1, self.npars))
        pass



    def postProcess(self, show=False, **kwargs):
        ''' 
        A method to post-process the sampler results 
        '''

        '''First find the median and 68% interval '''
        for i in range(self.npars):
            print((lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                 np.percentile(self.samples[:,i], [16, 50, 84]))#,
                                                    #axis=0)))
                  )

        ''' Then check what the "best fit" was '''
        print(np.min(self.sampler.lnprobability))
        print(np.max(self.sampler.lnprobability))
        ''' Now produce some diagnostic plots '''

        ''' A plot of the walkers' sampling history '''
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(self.npars, 1, sharex=True, figsize=(8, 9))
        for i in range(self.npars):
            axes[i].plot(self.sampler.chain[:, :, i].T, color="k", alpha=0.4)
            #    axes[0].yaxis.set_major_locator(MaxNLocator(5))
            #axes[0].axhline(m_true, color="#888888", lw=2)
            #axes[i].set_ylabel("$m$")
        
    
        #fig.tight_layout(h_pad=0.0)
        fig.savefig("walkers.png")
        #plt.close(fig)
        ''' And a corner plot of the post-burnin results '''
        import corner
        fig2 = corner.corner(self.samples)#,labels=opt.labels)
        fig2.savefig("corner.png")

        ''' plot the evolution of the likelihood function '''
        fig3 = plt.figure()
        axes3 = fig3.add_subplot(111) #plt.subplots(1, 1, sharex=True, figsize=(6, 8))
        for i in range(self.nwalkers):
            axes3.plot(self.sampler.lnprobability[i,:])
        try:
            axes3.set_ylim(np.min(self.sampler.lnprobability),np.max(self.sampler.lnprobability))
        except ValueError:
            axes3.set_ylim(-2000.,np.max(self.sampler.lnprobability))
        fig3.savefig("lnprob.png")
        pass


#    def lnlike(self, theta, **kwargs):
#        model = self.model(theta)
#        like = 0.
#        for data in self.dataSet:
#            """ do something for each bit of data """
#            deltaLike = something
#            like = like + deltaLike
#        return like #-0.5 * np.sum((((y - model)**2)/(yerr**2)))


