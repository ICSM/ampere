from __future__ import print_function

import numpy as np
import emcee
from .basesearch import BaseSearch
from inspect import signature
from .data import Photometry, Spectrum
import matplotlib.pyplot as plt

class EmceeSearch(BaseSearch):
    """
    A class to use the emcee affine-invariant ensemble sampler
    """

    def __init__(self, nwalkers = None, model = None,
                 data = None, lnprior = None,
                 labels = None, acceptRate = 2.0,
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
        #print(self.npars, self.nparsMod, self.nparsData)
        #exit()
        ''' then set up the sampler '''
        self.sampler = emcee.EnsembleSampler(self.nwalkers, np.int(self.npars), self.lnprob, a = acceptRate)#, args=self.dataSet)

 
        
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

    def rebuildSampler(self, nwalkers = None, model = None,
                       data = None, lnprior = None,
                       labels = None, acceptRate = None,
                       **kwargs):
        ''' This method will replace parts of the optimiser object if it is passed them. It can 
        be used to update the model, data, sampler setup, or prior part-way through a run '''

        if np.any([model, data]):
            if model is not None:
                self.model=model
                try:
                    self.nparsMod = self.model.npars
                except:
                    sig = signature(model.__call__)
                    self.nparsMod = len(sig.parameters) - 1 #Always subtract **kwargs from the parameters, but don't need to worry about self once it is bound to an instance
            if data is not None:
                self.dataSet = data
                self.nparsData = [data.npars for data in self.dataSet] #number of parameters to be passed into each set of data

            self.npars = np.int(self.nparsMod + np.sum(self.nparsData))
        
        if lnprior is not None:
            self.lnprior = lnprior
        #print(self.npars, self.nparsMod, self.nparsData)
        if nwalkers is not None:
                self.nwalkers=nwalkers
        ''' then set up the sampler '''
        
        if np.any([nwalkers, acceptRate, data, model, lnprior, labels]):
            ''' if anything has been changed, update the sampler and re-intialise it '''
            ''' first destroy the sampler '''
            self.sampler=None
            ''' now rebuild it '''
            self.sampler = emcee.EnsembleSampler(self.nwalkers, np.int(self.npars), self.lnprob, a = acceptRate)

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

    def optimise(self, nsamples = None, burnin = None, guess = None, noguess=False, **kwargs):
        if guess == 'None': # and not noguess:
            guess = [np.random.randn(np.int(self.npars)) for i in range(self.nwalkers)]
        #print('dimensions of guess = ', np.shape(guess))
        #print('nsamples = ', nsamples)
        #print('burnin = ', burnin)
        #print("np.max(guess, axis=0) = ", np.max(guess, axis=0))
        #print("np.min(guess, axis=0) = ", np.min(guess, axis=0))
        self.nsamp = nsamples
        self.burnin = burnin
        if noguess:
            self.sampler.run_mcmc(nsamples)
        else:
            self.sampler.run_mcmc(guess, nsamples)
        self.allSamples = self.sampler.chain
        #print('do we get here (no): ',np.max(self.allSamples), np.min(self.allSamples))
        self.samples = self.sampler.chain[:, self.burnin:, :].reshape((-1, self.npars))
        #print(np.max(self.samples), np.min(self.samples))
        #pass



    def postProcess(self, show=False, **kwargs):
        ''' 
        A method to post-process the sampler results 
        '''

        '''First find the median and 68% interval '''
        res=[]
        print("Median and confidence intervals for parameters in order:")
        for i in range(self.npars):
            a = np.percentile(self.samples[:,i], [16, 50, 84])
            res.append([a[1], a[2]-a[1], a[1]-a[0]])
            #res.append((lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
            #                     *zip(np.percentile(self.samples[:,i], [16, 50, 84])))#,
                                                    #axis=0)))
            #      )
            print(res[i])
            

        ''' Then check what the "best fit" was '''
        print(np.min(self.sampler.lnprobability))
        print(np.max(self.sampler.lnprobability))
        row_ind, col_ind = np.unravel_index(np.argmax(self.sampler.lnprobability.ravel), self.sampler.lnprobability.shape)
        bestPars = self.sampler.chain[row_ind, col_ind, :]
        print("MAP Solution: ", bestPars)
        ''' Now produce some diagnostic plots '''

        ''' A plot of the walkers' sampling history '''
        self.plot_walkers()
        
        #plt.close(fig)
        ''' And a corner plot of the post-burnin results '''
        self.plot_corner()

        ''' plot the evolution of the likelihood function '''
        self.plot_lnprob()

        '''finally, we want to look at the correlation/covariance matrices for the spectra, if any '''
        self.plot_covmats()
        #fig4,(ax0,ax1) = plt.subplots(1,2)
        #ax=[ax0, ax1]
        #i=0
        #for d in self.dataSet:
        #    if isinstance(d, Photometry):
        #        continue
        #    elif isinstance(d, Spectrum):
        #        i+=1
        #        fig, ax = plt.subplots(1,1)
        #        #for d in self.dataSet[1:]:
        #        #d=self.dataSet[1]
        #        d.cov([res[8][0],res[9][0]])
        #        #ax0.set_title('Covariance matrix')
        #        im = ax.imshow(np.log10(d.covMat))
        #        fig.savefig("covMat_"+str(i)+".png")

        
        #cax0 = divider4.append_axes("right", size="20%", pad=0.05)
        #cbar0 = plt.colorbar(im0, cax=cax0)
        #d=self.dataSet[2]
        #d.cov([res[11][0],res[12][0]])
            #ax0.set_title('Covariance matrix')
        #im1 = ax1.imshow(np.log10(d.covMat))
        #cax1 = divider4.append_axes("right", size="20%", pad=0.05)
        #cbar1 = plt.colorbar(im1, cax=cax1)
        #fig4.savefig("covMat.png")
        
        pass
    def plot_walkers(self):
        fig, axes = plt.subplots(self.npars, 1, sharex=True, figsize=(8, 9))
        for i in range(self.npars):
            axes[i].plot(self.sampler.chain[:, :, i].T, color="k", alpha=0.4)
            #    axes[0].yaxis.set_major_locator(MaxNLocator(5))
            #axes[0].axhline(m_true, color="#888888", lw=2)
            #axes[i].set_ylabel("$m$")
        
    
        #fig.tight_layout(h_pad=0.0)
        fig.savefig("walkers.png")
    
    def plot_corner(self):
        import corner
        fig2 = corner.corner(self.samples)#,labels=opt.labels)
        fig2.savefig("corner.png")

    def plot_lnprob(self):
        fig3 = plt.figure()
        axes3 = fig3.add_subplot(111) #plt.subplots(1, 1, sharex=True, figsize=(6, 8))
        for i in range(self.nwalkers):
            axes3.plot(self.sampler.lnprobability[i,:])
        try:
            axes3.set_ylim(np.min(self.sampler.lnprobability),np.max(self.sampler.lnprobability))
        except ValueError:
            try:
                axes3.set_ylim(-2000.,np.max(self.sampler.lnprobability))
            except ValueError:
                axes3.set_ylim(-2000.,2000)
        fig3.savefig("lnprob.png")        
        
    def plot_covmats(self):
        for i, d in enumerate(self.dataSet):
            if isinstance(d, Photometry):
                continue
            elif isinstance(d, Spectrum):
                fig, ax = plt.subplots(1,1)
                #for d in self.dataSet[1:]:
                #d=self.dataSet[1]
                d.cov([res[8][0],res[9][0]])
                #ax0.set_title('Covariance matrix')
                im = ax.imshow(np.log10(d.covMat))
                fig.savefig("covMat_"+str(i)+".png")

    def plot_posteriorpredictive(self):
        pass

#    def lnlike(self, theta, **kwargs):
#        model = self.model(theta)
#        like = 0.
#        for data in self.dataSet:
#            """ do something for each bit of data """
#            deltaLike = something
#            like = like + deltaLike
#        return like #-0.5 * np.sum((((y - model)**2)/(yerr**2)))


