"""
A file for simulation-based inference approaches, starting with the sbi package

"""

from __future__ import print_function

import numpy as np
import torch
import pickle
import logging
from tqdm import tqdm
from datetime import datetime
from .basesearch import BaseSearch
from .mixins import SBIPostProcessor
from ..logger import Logger
from ..models.results import ModelResults
from sbi.inference import SNPE
from sbi import utils as utils
from sbi import analysis as analysis

from sbi.utils import process_prior #CustomPriorWrapper

class CustomPriorWrapper:
    pass

class LFIBase(BaseSearch, Logger):
    """ A base class for simulation-based (or likelihood-free) inference 
        approaches
    
    
    """

    _inference_method="LFI"


    def __init__(self, model=None, data=None, verbose = False,
                 parameter_labels = None, cache_models = True,
                 name='', namestyle="full",
                 **kwargs):
        self.model = model
        self.dataSet = data
        self.cache_models = cache_models
        if cache_models:
            self.cache = {}

        try:
            modelname = model.name
        except AttributeError:
            modelname = model.__class__.__name__
        if namestyle=='full':
            self.name = "ampere_"+str(datetime.now()).replace(' ','_').replace(":","-")[:-7] + "_" + modelname + "_" + name
        elif namestyle=="short":
            self.name = name
        elif namestyle=="stamp":
            self.name = "ampere_"+str(datetime.now()).replace(' ','_').replace(":","-")[:-7] + "_" + name
        elif namestyle=="model":
            self.name = "ampere_"+modelname + "_" + name

        self.setup_logging(verbose=verbose) #For now we will exclusively use default logging settings, this will be modified once logging is tested.
        self.logger.info("Welcome to ampere")
        self.logger.info("Setting up your inference problem:")
        self.logger.info("You are using %s", self._inference_method)
        self.logger.info("You have %s items in your dataset", str(len(data)))

        ''' now do some introspection on the various bits of model to 
        understand how many parameters there are for each compponent '''
        try:
            self.nparsMod = self.model.npars_ptform
        except:
            sig = signature(model.__call__)
            self.nparsMod = len(sig.parameters) - 1 #Always subtract **kwargs from the parameters, but don't need to worry about self once it is bound to an instance
        self.nparsData = [data.npars for data in self.dataSet] #number of parameters to be passed into each set of data
        self.npars = np.int(self.nparsMod + np.sum(self.nparsData))

        if parameter_labels is None:
            #The user hasn't specified parameter labels, let's see if the models and data have instead
            try: #First the model parameters
                self.parLabels = self.model.parLabels
            except AttributeError:
                self.parLabels = ['x'+str(i) for i in range(self.nparsMod)]
            i = self.nparsMod
            for data in self.dataSet:
                try:
                    self.parLabels.extend(data.parLabels)
                except AttributeError:
                    self.parLabels.extend(['x'+str(i) for i in range(i, i+data.npars)])
                finally:
                    i+=data.npars
        else: #User isn't really supposed to use this interface, however...
            self.parLabels = parameter_labels
        #self.parLabels = ['x'+str(i) for i in range(self.npars)] #Parameter for parameter names (labels) to associate with output in post processing - emcee does this internally with parameter_names, but we want a universal system across all the search methods. Called parLabels to distinguish it from Data Class labels.

        self.verbose = verbose
        #if self.verbose:
        self.logger.info("This model has %d parameters.", self.nparsMod)
        self.logger.info("There are also %d parameters for the noise model", self.npars - self.nparsMod)
        self.logger.info("Hence, there are a total of %d parameters to sample", self.npars)
        self.logger.info("The parameter names are:")
        for l in self.parLabels:
            self.logger.info("%s", l)

        #Now we should check whether the number of parameters matches with the parameter labels!
        if len(self.parLabels) != self.npars:
            self.logger.critical("You have %d free parameters but %d parameter labels", self.npars, len(self.parLabels))
            self.logger.critical("Please check the number of parameters and labels")
            raise ValueError("Mismatch between number of parameters and labels")

    def simulate(self, theta): #, theta_noise):
        """A method to simulate an observation """

        #First we call the model
        if self.cache_models: #There should be a better way of updating the cache than this, but for now this is fine...
            if (key := tuple(theta[:self.nparsMod])) in self.cache:
                result = self.cache[key]
            else:
                result = self.model(*theta[:self.nparsMod])
        else:
            result = self.model(*theta[:self.nparsMod])
        if not isinstance(result, ModelResults):
            result = ModelResults(**result) #Try unpacking the results here if the user didn't already define their model with it
        #Cache the result
        if self.cache_models:
            self.cache[tuple(theta[:self.nparsMod])] = result #There's some overhead here, because it will overwrite the cache if it's already cached...

        #Then we make each Data object produce a sample synthetic observation
        sims = []
        i=self.nparsMod
        for data in self.dataSet:
            sim = data.simulate(theta[i:i+data.npars],result)
            sims.extend(sim)
            i+=data.npars

        return sims #Now we return the simulated data. Inidividual implementations should modify this to accomodate their idiosyncracies


    def simulate_vector(self, thetas):
        """ A method to compute a batch of simulations """
        #This method takes the same approach as simulate(), but iterates over many thetas instead.
        n_sims = thetas.shape[0]
        for j, theta in enumerate(tqdm(thetas)):

            if self.cache_models:
                if (key := tuple(theta[:self.nparsMod])) in self.cache:
                    result = self.cache[key]
                else:
                    result = self.model(*theta[:self.nparsMod])
            #print(theta)
            #print(theta[:self.nparsMod])
            #First we call the model
            else:
                result = self.model(*theta[:self.nparsMod])
            if not isinstance(result, ModelResults):
                result = ModelResults(**result) #Try unpacking the results here if the user didn't already define their model with it
            #Cache the result
            if self.cache_models:
                self.cache[tuple(theta[:self.nparsMod])] = result #There's some overhead here, because it will overwrite the cache if it's already cached...
            #Then we make each Data object produce a sample synthetic observation
            sims = []
            i=self.nparsMod
            for data in self.dataSet:
                sim = data.simulate(theta[i:i+data.npars],result)
                sims.extend(sim)
                i+=data.npars
            if j == 0: #first time though we need to create the array for the results
                all_sims = np.zeros((n_sims, len(sims)))
            all_sims[j,:] = sims

        return all_sims

    def simulate_from_cached_models(self, thetas):
        pass

    def save_cached_models(self, format="pickle"):
        #first we're going to get all the keys for the models so we know what we're outputting
        keys = self.cache.keys()

        #Then we have to open output files in the requested formats
        #For now we're implementing pickling, but we should also offer:
        #Dill
        #HDF5
        #JSON
        #NetCDF?
        pass


    def load_cached_models(self):
        pass

            
    
class SBI_SNPE(LFIBase,SBIPostProcessor):
    _inference_method = "Sequential Neural Posterior Estimation"

    def __init__(self, model=None, data=None, verbose = False,
                 parameter_labels = None,
                 name='', namestyle="full",
                 check_prior_normalisation = True,
                 n_prior_norm_samples = 100000,
                 prior_norm_thres = 0.01,
                 cache_models = True,
                 **kwargs):
        super().__init__(model = model, data= data, verbose = verbose,
                         parameter_labels = parameter_labels,
                         cache_models = cache_models,
                         name = name, namestyle=namestyle, **kwargs)

        if check_prior_normalisation:
            self.logger.info("Checking prior normalisation")
            #Now we need to jump through a few hoops to get the prior into a format that SBI understands
            #We need to check if lnprior is normalised
            #We will do this with MC integration
            #First, we generate samples (default 100000) from the prior transform.
            #However, since the prior may have infinite support in some or all variables, we restrict ourselves to
            #a volume 0.5**npars so that we are computing a *fixed* and *known* fraction of the integral
            #which should have a value of 1 * (0.5**npars) if the distribution is normalised
            self.logger.info("Drawing prior samples to test normalisation")
            u = np.random.default_rng().uniform(low = 0.25, high=0.75, size=(n_prior_norm_samples, self.npars))
            samples = np.zeros_like(u)
            #Now we call the prior transform. For safety, we will manually iterate over the batch dimension

            for i in tqdm(range(n_prior_norm_samples)):
                samples[i] = self.prior_transform(u[i,:])
            self.logger.info("Computing integral")
            lp = self.lnprior_vector(samples)
            from scipy.special import logsumexp
            ranges = np.max(samples, axis=0) - np.min(samples, axis=0)
            log_int_est = logsumexp(lp) - np.log(n_prior_norm_samples) + np.sum(np.log(ranges))
            if log_int_est - self.npars*np.log(0.5) < prior_norm_thres: #We can consider this "close enough", since it corresponds to just over 1% for the default input
                self.logger.debug("Prior is normalised")
                self._prior_is_normalised = True
            else:
                #If the deviation is more than 1%, we're going to assume the prior is actually not normalised,
                self.logger.warning("Prior does not appear to be normalised! Inference will proceed with an approximate normalisation, but it would be better to ensure it is normalised in future. To prevent this warning in future, pass ``check_prior_normalisation = False'' to the inference object")
                self._prior_is_normalised = True
                self.logprior_norm = log_int_est
        
        prior, *_ = process_prior(self)#not sure if this kind of recursion will work :/ In principle
                                        #it should, since and SBI_SNPE instance provides both sample() and logprob()
                                        #methods which meet the requirements.
        self.sampler = SNPE(prior=prior)
        pass


    def optimise(self, nsamples = None, nsamples_post = None,
                 **kwargs):

        self.logger.info("Preparing to sample")

        #first we draw a batch of samples from the prior
        self.logger.info("Generating samples from the prior")
        thetas = self.sample(nsamples)
        #Now we simulate models for all of these samples
        self.logger.info("Simulating all prior samples")
        sims = self.simulate_vector(thetas)

        #now that we have simulations and parameter values, we're going to make sure they're torch Tensors
        thetas = torch.Tensor(thetas)#, dtype=torch.float32)
        sims = torch.Tensor(sims)#, dtype=torch.float32)

        #now we give the sampler the simulations ...
        self.inference = self.sampler.append_simulations(thetas, sims)
        #...and start doing inference
        self.logger.info("Training posterior")
        self.density_estimator = self.inference.train()
        self.posterior = self.inference.build_posterior(self.density_estimator)

        #now that the posterior is trained, we can apply it to our observation
        #first we need to check how many observations we're dealing with
        obs = []
        for d in self.dataSet:
            #print(d.value[d.mask])
            obs.extend(d.value[d.mask]) #for now, we'll do this the ugly way, and improve it later.
        #print(obs)
        #print(sims[-1])
        obs = torch.Tensor(obs)
        self.logger.info("Sampling trained posterior")
        self.posterior.set_default_x(obs)
        #self.samples = self.posterior.sample((nsamples_post,), x = obs)
        self.samples = self.posterior.sample((nsamples_post,))
        #for posterity (and post-processing) we're going to save the parameter values for the prior samples
        self.thetas = thetas
        check_cache=False#True
        if check_cache:
            self.logger.info("Checking model cache")
            self.logger.info("%s models were cached in this run", len(self.cache))
            self.logger.info("The keys are:")
            self.logger.info("%s", self.cache.keys())
        pass

    def postProcess(self, show=False, textfile=None, **kwargs):
        ''' 
        A method to post-process the sampler results 
        '''

        self.print_summary(**kwargs)

        self.plot_corner(**kwargs)

        self.plot_posteriorpredictive(**kwargs)


    def sample(self, nsamples = 1, **kwargs):
        """ 
        Produce samples from the prior

        The sbi package requires priors that have sample() and logprob() methods
        to generate samples and evaluate the prior. Seeing as we have effectively 
        already designed such functions for other purposes, we 
        """
        
        if isinstance(nsamples, torch.Size):
            nsamples = 1
        elif isinstance(nsamples, tuple):
            #if len
            #print(nsamples)
            #print(len(nsamples))
            nsamples = nsamples[0]
            
        #First, we generate nsamples samples of self.npars uniform random variates
        u = np.single(np.random.default_rng().uniform(size=(nsamples, self.npars)))
        samples = np.zeros((nsamples, self.npars))

        #Now we call the prior transform. For safety, we will manually iterate over the batch dimension
        for i in range(nsamples):
            samples[i] = self.prior_transform(u[i,:])
        return samples

    def log_prob(self, theta, **kwargs):
        """
        Calculate the logprobability of the prior
        """

        if self._prior_is_normalised:
            return self.lnprior_vector(theta)
        else:
            return self.lnprior_vector(theta) - self.logprior_norm

        

        
