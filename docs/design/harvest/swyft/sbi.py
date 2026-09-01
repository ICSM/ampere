"""
A file for simulation-based inference approaches, starting with the sbi package

"""

from __future__ import print_function
from typing import Optional, Any, Tuple, List


import numpy as np
import torch
import pickle
import logging
import copy
from tqdm import tqdm
import matplotlib.pyplot as plt
from datetime import datetime
from .basesearch import BaseSearch
from .mixins import SBIPostProcessor
from ..logger import Logger
from ..models.results import ModelResults
from sbi.inference import SNPE, DirectPosterior
from sbi import utils as utils
from sbi import analysis as analysis

from sbi.utils import process_prior  # CustomPriorWrapper

try:
    import swyft
    is_swift_installed = True
except (ModuleNotFoundError, ImportError):
    is_swift_installed = False


class CustomPriorWrapper:
    pass


class LFIBase(BaseSearch, Logger):
    """ A base class for simulation-based (or likelihood-free) inference 
        approaches

    This should be subclassed for specific approaches
    
    Parameters
    ---------
    model : function
        The model that is used to fit the data.
    data : list of Data objects
        The data that the model is fit to.
    verbose : bool, optional
        Flag to print output messages. Default is False.
    parameter_labels : list of str, optional
        Labels to give to each parameter. Default is None.
    cache_models : bool, optional
        Flag to cache models. Default is True.
    name : str, optional
        Name to give to the Inference object. Default is empty string.
    namestyle : str, optional
        Style to give to the name of the Inference object. Options are 
        "full", "short", "stamp", and "model". Default is "full".
    **kwargs : dict, optional
        Presently provided only for compatibility
    """

    _inference_method = "LFI"


    def __init__(self, model=None, data=None, verbose=False,
                 parameter_labels=None, cache_models=True,
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
            self.name = "ampere_"+str(datetime.now()).replace(' ', '_').replace(":", "-")[:-7] + "_" + modelname + "_" + name
        elif namestyle=="short":
            self.name = name
        elif namestyle=="stamp":
            self.name = "ampere_"+str(datetime.now()).replace(' ', '_').replace(":", "-")[:-7] + "_" + name
        elif namestyle=="model":
            self.name = f"ampere_{modelname}_{name}"

        self.setup_logging(verbose=verbose) #For now we will exclusively use default logging settings, this will be modified once logging is tested.
        self.logger.info("Welcome to ampere")
        self.logger.info("Setting up your inference problem:")
        self.logger.info("You are using %s", self._inference_method)
        self.logger.info("You have %s items in your dataset", str(len(data)))

        ''' now do some introspection on the various bits of model to 
        understand how many parameters there are for each compponent '''
        try:
            self.nparsMod = self.model.npars_ptform
        except AttributeError:
            sig = signature(model.__call__)
            self.nparsMod = len(sig.parameters) - 1 #Always subtract **kwargs from the parameters, but don't need to worry about self once it is bound to an instance
        self.nparsData = [data.npars for data in self.dataSet] #number of parameters to be passed into each set of data
        self.npars = int(self.nparsMod + np.sum(self.nparsData))

        if parameter_labels is None:
            # The user hasn't specified parameter labels, let's see if the models and data have instead
            try: #First the model parameters
                self.parLabels = self.model.parLabels
            except AttributeError:
                self.parLabels = [f'x{str(i)}' for i in range(self.nparsMod)]
            i = self.nparsMod
            for data in self.dataSet:
                try:
                    self.parLabels.extend(data.parLabels)
                except AttributeError:
                    self.parLabels.extend([f'x{str(i)}' for i in range(i, i+data.npars)])
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
        """A method to simulate an observation 

        Simulate an observation using the current model and noise model parameters.

        Parameters
        ----------
        theta : array-like
            Model parameter values to be simulated
        
        Returns
        ---------
        sims : list
            A list of simulated observations corresponding to the input data objects.
        """

        #First we call the model
        if self.cache_models: 
            # There should be a better way of updating the cache than this, but for 
            # now this is fine...
            if (key := tuple(theta[:self.nparsMod])) in self.cache:
                result = self.cache[key]
            else:
                result = self.model(*theta[:self.nparsMod])
        else:
            result = self.model(*theta[:self.nparsMod])
        if not isinstance(result, ModelResults):
            result = ModelResults(**result) # Try unpacking the results here if the user
                                            # didn't already define their model with it
        #Cache the result
        if self.cache_models:
            # There's some overhead here, because it will overwrite the cache if it's 
            # already cached...
            self.cache[tuple(theta[:self.nparsMod])] = result 

        # Then we make each Data object produce a sample synthetic observation
        sims = []
        i=self.nparsMod
        for data in self.dataSet:
            sim = data.simulate(theta[i:i+data.npars],result)
            sims.extend(sim)
            i+=data.npars

        # Now we return the simulated data. Inidividual implementations should modify 
        # this to accomodate their idiosyncracies
        return sims 


    def simulate_vector(self, thetas):
        """ A method to compute a batch of simulations 

        Simulate a batch of observations using many parameter values.

        At present, this method simply loops over all sets of parameter values

        Parameters
        ----------
        theta : array-like
        Model parameter values to be simulated
        
        Returns
        ---------
        sims : list
        A list of simulated observations corresponding to the input data objects.

        Generated by Chat-GPT
        """
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


class SBI_SNPE(LFIBase, SBIPostProcessor):
    _inference_method = "Sequential Neural Posterior Estimation"

    def __init__(self, model=None, data=None, verbose=False,
                 parameter_labels=None,
                 name='', namestyle="full",
                 check_prior_normalisation=True,
                 get_prior_bounds=True,
                 n_prior_norm_samples=100000,
                 prior_norm_thres=0.01,
                 cache_models=True,
                 embedding_net=None,
                 **kwargs):
        super().__init__(model=model, data= data, verbose=verbose,
                         parameter_labels=parameter_labels,
                         cache_models=cache_models,
                         name=name, namestyle=namestyle, **kwargs)

        if check_prior_normalisation:
            self.logger.info("Checking prior normalisation")
            self._check_prior_normalisation(n_prior_norm_samples, 
                                            prior_norm_thres)

        if get_prior_bounds:
            lower_bounds, upper_bounds=self._check_prior_bounds()

            prior, *_ = process_prior(self, custom_prior_wrapper_kwargs=dict(lower_bound=torch.from_numpy(lower_bounds.astype(np.float32)), upper_bound=torch.from_numpy(upper_bounds.astype(np.float32))))
        else:
            prior, *_ = process_prior(self)

        self.prior = prior
        if embedding_net:
            from torch import nn
            data_shape = (np.sum([data.value.shape[0] for data in self.dataSet]) ,)#this highlights an important issue - data objects need to define their __len__ in an appropriate way, but this will do for now
            if embedding_net in [True, "default", "Conv", "CNN"]:
                from sbi.neural_nets.embedding_nets import CNNEmbedding
                embedding_net = CNNEmbedding(input_shape=data_shape, output_dim=2*self.npars)
                self.logger.info("Using default CNN embedding net")
            elif embedding_net in ["FC"]:
                from sbi.neural_nets.embedding_nets import FCEmbedding
                embedding_net = FCEmbedding(input_dim=data_shape[0], output_dim=2*self.npars)
                self.logger.info("Using default fully-connected embedding net")
            elif isinstance(embedding_net, nn.Module):
                self.logger.info("Using user-defined embedding net")
            elif isinstance(embedding_net, dict):
                #embedding_net is a dictionary of parameters for one of the default embedding nets
                net = embedding_net.get("type", "CNN")
                output_dim = embedding_net.get("output_dim", 2*self.npars)
                if net in ["CNN", "Conv"]:
                    from sbi.neural_nets.embedding_nets import CNNEmbedding
                    n_conv_layers = embedding_net.get("n_conv_layers", 2)
                    out_channels_per_layer = embedding_net.get("out_channels_per_layer", [6, 12])
                    kernel_size_per_layer = embedding_net.get("kernel_size_per_layer", 5)
                    n_linear_layers = embedding_net.get("n_linear_layers", 2)
                    num_linear_units = embedding_net.get("num_linear_units", 50)
                    pool_kernel_size = embedding_net.get("pool_kernel_size", 2)
                    embedding_net = CNNEmbedding(input_shape=data_shape, 
                                                 output_dim=output_dim, 
                                                 num_conv_layers = n_conv_layers, 
                                                 out_channels_per_layer=out_channels_per_layer, 
                                                 kernel_size_per_layer=kernel_size_per_layer, 
                                                 num_linear_layers=n_linear_layers, 
                                                 num_linear_units=num_linear_units, 
                                                 pool_kernel_size=pool_kernel_size)
                    self.logger.info("Using default CNN embedding net with user-defined parameters")
                    self.logger.info("Network:")
                    self.logger.info(embedding_net)
                elif net in ["FC", "FullyConnected"]:
                    from sbi.neural_nets.embedding_nets import FCEmbedding
                    n_layers = embedding_net.get("n_layers", 2)
                    num_hiddens = embedding_net.get("num_hiddens", 40)
                    embedding_net = FCEmbedding(input_dim=data_shape[0],
                                                output_dim=output_dim,
                                                num_layers=n_layers,
                                                num_hiddens=num_hiddens)
                    self.logger.info("Using default fully-connected embedding net with user-defined parameters")
                    self.logger.info("Network:")
                    self.logger.info(embedding_net)
            elif isinstance(embedding_net, str):
                raise ValueError(f"Unknown embedding net type {embedding_net}")
            else:
                raise ValueError(f"Unknown embedding net type {embedding_net}")
            neural_posterior = utils.posterior_nn(model="maf", 
                                                  embedding_net=embedding_net, 
                                                  hidden_features=10, 
                                                  num_transforms=2
                                                  )
            self.sampler = SNPE(prior=prior, density_estimator=neural_posterior)
        else:
            self.sampler = SNPE(prior=prior)

    def _check_prior_normalisation(self, n_prior_norm_samples=100000, 
                                   prior_norm_thres=0.01):
        """ Check the normalisation of the prior by integrating it over the prior volume.

        This method uses MC integration to check the normalisation of the prior. It 
        draws a number of samples from the prior (default 100000) and then integrates 
        the prior over the central 50% of each parameter of the prior - this guarantees 
        that we integrate over a *fixed* and *known* fraction of the prior volume (and 
        hence integral). As a result, the integral should be approximately equal to 0.5 
        to the power of the number of parameters. For numerical stability, we work with 
        the log of the integral, which hence should be equal to n_pars*log(0.5). If the 
        absolute value of the difference between the log of the integral and this value 
        is less than prior_norm_thres, then the prior is considered to be normalised and
        will be used as is in inference. Otherwise, the normalisation factor will be 
        computed and is applied to the prior logprob function during inference.

        Parameters
        ----------
        n_prior_norm_samples : int, optional
            Number of samples to draw for MC integration of the prior, by default 100000
        prior_norm_thres : float, optional
            _description_, by default 0.01
        """
        self.logger.info("Drawing prior samples to test normalisation")
        u = np.random.default_rng().uniform(low=0.25, high=0.75, size=(n_prior_norm_samples, self.npars))
        samples = np.zeros_like(u)
        # Now we call the prior transform. For safety, we will manually iterate over the
        # batch dimension
        for i in tqdm(range(n_prior_norm_samples)):
            samples[i] = self.prior_transform(u[i,:])
        self.logger.info("Computing integral")
        lp = self.lnprior_vector(samples)
        from scipy.special import logsumexp
        ranges = np.max(samples, axis=0) - np.min(samples, axis=0)
        log_int_est = logsumexp(lp) - np.log(n_prior_norm_samples) + np.sum(np.log(ranges))
        if np.abs(log_int_est - self.npars*np.log(0.5)) < prior_norm_thres: 
            # We can consider this "close enough", since it corresponds to just over 1% 
            # for the default input
            self.logger.debug("Prior is normalised")
            self._prior_is_normalised = True
        else:
            # If the deviation is more than prior_norm_thres, we're going to assume the 
            # prior is actually not normalised,
            self.logger.warning("Prior does not appear to be normalised! Inference will proceed with an approximate normalisation, but it would be better to ensure it is normalised in future. To prevent this warning in future, pass ``check_prior_normalisation=False'' to the inference object")
            self._prior_is_normalised = False
            self.logprior_norm = log_int_est

    def _check_prior_bounds(self):
        """Check the bounds of the prior.

        _extended_summary_

        Returns
        -------
        lower_bounds, np.ndarray
            _description_
        """
        # now we should try to establish what the bounds of the prior are so we can pass
        # it to SBI to limit the number of samples outside the range.
        self.logger.info("Trying to determine prior bounds from prior transforms")
        bounds_test = self.prior_transform(np.zeros(self.npars))
        for i in range(self.npars): 
            # Now we can iterate over all the parameters and get the maximum value of 
            # each one in turn
            u = np.zeros(self.npars)
            u[i] = 1 #
            bounds_test = np.vstack((bounds_test, self.prior_transform(u)))
        lower_bounds = np.nanmin(bounds_test, axis=0)
        upper_bounds = np.nanmax(bounds_test, axis=0)

        self.logger.info("Inferred prior bounds:")
        for i in range(self.npars):
            self.logger.info("%s: Lower bounds =  %.5f, upper bounds = %.5f", self.parLabels[i], lower_bounds[i], upper_bounds[i])
        return lower_bounds, upper_bounds

    def optimise(self, nsamples=None, nsamples_post=None, n_rounds=1,
                 **kwargs):

        """Simulate samples from a model and optimise the posterior.

        Parameters
        ----------
        nsamples : int, optional
            The number of samples to draw from the proposal each round.
        nsamples_post : int, optional
            The number of samples to draw from the final trained posterior.
        n_rounds : int, optional
            The number of rounds of sampling to do. If n_rounds > 1, the prior is used
            as the proposal in the first round, and the previous round's posterior is 
            used as the new proposal in each subsequent round.
        **kwargs : optional
            Keyword arguments provided for compatibility and consistency

        Returns
        -------
        None
        """

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
        #self.posterior = self.inference.build_posterior(self.density_estimator)
        self.posterior = DirectPosterior(self.density_estimator, self.prior, enable_transform=False) #Doing it this way avoids problems with re-transformation of samples by stopping SBI from applying any of its own transformations.

        #now that the posterior is trained, we can apply it to our observation
        #first we need to check how many observations we're dealing with
        obs = []
        for d in self.dataSet:
            #print(d.value[d.mask])
            obs.extend(d.value[d.mask]) #for now, we'll do this the ugly way, and improve it later.
        #print(obs)
        #print(sims[-1])
        obs = torch.Tensor(obs)
        #self.logger.info("Sampling trained posterior")
        self.posterior.set_default_x(obs)
        if n_rounds > 1:
            self.posteriors = []
            #self.caches = []
            self.thetas_hist = []
            self.sims_hist = []
            self.posteriors.append(copy.deepcopy(self.posterior))
            #self.caches.append(copy.deepcopy(self.cache))
            self.thetas_hist.append(copy.deepcopy(thetas))
            self.sims_hist.append(copy.deepcopy(sims))
            for i in range(n_rounds-1):
                self.logger.info("Starting inference round %i", i+2)
                #self.cache = {}
                proposal = self.posterior.set_default_x(obs)
                theta = proposal.sample((nsamples,))
                self.logger.info("Generating round %i samples", i+2)
                sims = self.simulate_vector(thetas.detach().numpy())
                thetas = torch.Tensor(thetas)#, dtype=torch.float32)
                sims = torch.Tensor(sims)#, dtype=torch.float32)
                self.thetas_hist.append(copy.deepcopy(thetas))
                self.sims_hist.append(copy.deepcopy(sims))
                self.logger.info("Training round %i posterior", i+2)
                density_estimator = self.sampler.append_simulations(thetas, sims, proposal=proposal).train()
                posterior = DirectPosterior(density_estimator, self.prior, enable_transform=False)
                self.posteriors.append(posterior)
                #self.caches.append(self.cache)
            self.posterior = posterior
            self.posterior.set_default_x(obs)
            
        #self.samples = self.posterior.sample((nsamples_post,), x=obs)
        self.logger.info("Sampling trained posterior")
        self.samples = self.posterior.sample((nsamples_post,))
        #for posterity (and post-processing) we're going to save the parameter values for the prior samples
        self.thetas = thetas
        check_cache=False#True
        if check_cache:
            self.logger.info("Checking model cache")
            self.logger.info("%s models were cached in this run", len(self.cache))
            self.logger.info("The keys are:")
            self.logger.info("%s", self.cache.keys())


    def postProcess(self, show=False, textfile=None, **kwargs):
        ''' 
        A method to post-process the sampler results 

        Post-process the results of the sampler.

        Parameters
        ----------
        show : bool, optional
            Flag to show the plots after creation.
        textfile : str, deprecated
            Filename for a text file to save the results to.
        **kwargs : additional keyword arguments
            Additional arguments to pass to the print_summary, plot_corner, and 
            plot_posteriorpredictive functions.
        '''

        self.print_summary(**kwargs)

        self.plot_corner(**kwargs)

        self.plot_posteriorpredictive(**kwargs)


    def sample(self, nsamples=1, **kwargs):
        """ 
        Produce samples from the prior

        The sbi package requires priors that have sample() and logprob() methods
        to generate samples and evaluate the prior. Seeing as we have effectively 
        already designed such functions for other purposes, we simply define a wrapper
        to them here

        Parameters
        ----------
        nsamples : int
            number of prior samples to generate
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

        The sbi package requires priors that have sample() and logprob() methods
        to generate samples and evaluate the prior. Seeing as we have effectively 
        already designed such functions for other purposes, we simply define a wrapper
        to them here

        Parameters
        ----------
        theta : array-like
            parameter values to evaluate the prior for
        """

        if self._prior_is_normalised:
            try:
                return self.lnprior_vector(theta)
            except TypeError:
                return self.lnprior_vector(theta.detach().numpy())
        else:
            try:
                return self.lnprior_vector(theta) - self.logprior_norm
            except TypeError:
                return self.lnprior_vector(theta.detach().numpy()) - self.logprior_norm



### TMNRE with Swyft Stars here:

class SwyftNetwork(swyft.SwyftModule):
    def __init__(self, npars, num_features, parLabels, marginals=None):
        super().__init__()
        #self.embedding = torch.nn.Linear(self.num_features, self.npars)
        # If marginals is not defined, use all combinations of parameters
        # create a tuple of tuples that contains the indices of all the marginals
        # This needs all the combinations in the lower triangle of the matrix
        if marginals is None:
            marginals = tuple((i, j) for i in range(npars) for j in range(i+1, npars))
        self.logratios1 = swyft.LogRatioEstimator_1dim(num_features = num_features,
                                                        num_params = npars,
                                                        varnames = parLabels)
        self.logratios2 = swyft.LogRatioEstimator_Ndim(num_features = num_features,
                                                        marginals = marginals,
                                                        varnames = self.parLabels)

    def forward(self, A, B):
        logratios1 = self.logratios1(A['data'], B['pars'])
        logratios2 = self.logratios2(A['data'], B['pars'])
        return logratios1, logratios2


class SwyftNetworkLinearEmbedding(swyft.SwyftModule):
    def __init__(self, npars, num_features, parLabels, marginals=None):
        super().__init__()
        # If marginals is not defined, use all combinations of parameters
        if marginals is None:
            marginals = tuple((i, j) for i in range(npars) for j in range(i+1, npars))
        self.embedding = torch.nn.Linear(num_features, npars)
        self.logratios1 = swyft.LogRatioEstimator_1dim(num_features =npars,
                                                        num_params = npars,
                                                        varnames = parLabels)
        self.logratios2 = swyft.LogRatioEstimator_Ndim(num_features = npars,
                                                        marginals = marginals,
                                                        varnames = parLabels)

    def forward(self, A, B):
        embedding = self.embedding(A['data'])
        logratios1 = self.logratios1(embedding, B['pars'])
        logratios2 = self.logratios2(embedding, B['pars'])
        return logratios1, logratios2

class Swyft_TMNRE(LFIBase):
    """A class to perform likelihood-free inference using the Swyft package.

    Swyft is a package that provides a framework for likelihood-free inference using
    Truncated Marginal Neural Ratio Estimation (TMNRE). TMNRE is a method for
    likelihood-free inference that targets the marginal likelihood ratio rather than
    the joint likelihood ratio. As a result, it typically requires fewer simulations
    than other methods, and is able to scale to high-dimensional problems. However, it
    is not amortised, and so it is most effective for working with small samples, high-
    dimensional parameter spaces, and particularly expensive models.

    Parameters
    ----------
    model : Callable
        The model to be fit to the data.
    data : list of Data objects
        The data to be fit to the model.
    verbose : bool, optional
        Flag to print output messages. Default is False.
    parameter_labels : list of str, optional
        Labels to give to each parameter. Default is None.
    """

    _inference_method = "Truncated Marginal Neural Ratio Estimation"

    def __init__(self, model=None, data=None, verbose=False,
                 parameter_labels=None,
                 name='', namestyle="full", accelerator=None,
                 **kwargs):
        if not is_swift_installed:
            raise ImportError("""The Swyft package is not installed.
                              To enable TMNRE, please install it using
                              `>>> pip install swyft`""")
        super().__init__(model=model, data= data, verbose=verbose,
                         parameter_labels=parameter_labels,
                         name=name, namestyle=namestyle, **kwargs)

        # Set the device to use for training
        # This rather long conditional assignment is to ensure that the device is set
        # while allowing the user to override it if they want to
        self.DEVICE = (
            'gpu'
            if torch.cuda.is_available() and accelerator == 'cuda'
            else accelerator
            if accelerator is not None
            else 'cpu'
        )

        # TODO: Add a check for the device to ensure that it is valid
        # TODO: Enable using multiple GPUs if available

    def optimise(self,
                 nsamples: Optional[int] = 10000,
                 nsamples_post: Optional[int] = 10000,
                 n_rounds: int = 1,
                 embedding_net: bool = True,
                 marginals: Optional[Tuple[int, Tuple[int], Any]] = None,
                 num_workers: Optional[int] = 0,
                 **kwargs: Any) -> None:

        #first we draw a batch of samples from the prior
        self.logger.info("Generating samples from the prior")
        thetas = self.sample(nsamples)
        #Now we simulate models for all of these samples
        self.logger.info("Simulating all prior samples")
        sims = self.simulate_vector(thetas)
        # Now translate them to swyft format
        sims_swyft = swyft.Samples(pars = thetas, data = sims)

        # Now we define the SwyftModule that will be used to train the TMNRE
        self.logger.info("Defining SwyftModule")
        self.num_features = (np.sum([data.value.shape[0] for data in self.dataSet]) ,)[0]


        trainer = swyft.SwyftTrainer(accelerator = self.DEVICE, precision = 64)
        dm = swyft.SwyftDataModule(sims_swyft, num_workers=num_workers)
        if embedding_net is True: # Being specific so that there can be more cases in
                                  # the future
            self.logger.info("Using default embedding net")
            network = SwyftNetworkLinearEmbedding(self.npars,
                                                  self.num_features,
                                                  self.parLabels,
                                                  marginals)
        else:
            network = SwyftNetwork(self.npars,
                                   self.num_features,
                                   self.parLabels,
                                   marginals)
        self.logger.info("Training TMNRE")
        trainer.fit(network, dm)

        # Now we can use the trained network to sample from the posterior
        # First we need to generate samples from the prior
        self.logger.info("Generating samples from the prior for inference")
        prior_samples = self.sample(nsamples_post)
        prior_samples = swyft.Samples(pars = prior_samples)
        obs = []
        for d in self.dataSet:
            obs.extend(d.value[d.mask])
        obs_swyft = swyft.Samples(data = torch.Tensor(obs))
        print(f'obs_swyft: {obs_swyft}')
        print(f'obs_swyft["data"].shape: {obs_swyft["data"].shape}')
        self.logger.info("Drawing samples from trained TMNRE")
        self.predictions = trainer.infer(network, obs_swyft, prior_samples)
        self.logger.info("TMNRE inference complete")

    def sample(self, nsamples: int) -> np.ndarray:
        """ 
        Produce samples from the prior

        The sbi package requires priors that have sample() and logprob() methods
        to generate samples and evaluate the prior. Seeing as we have effectively 
        already designed such functions for other purposes, we simply define a wrapper
        to them here

        Parameters
        ----------
        nsamples : int
            number of prior samples to generate
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

    def postProcess(self, **kwargs):
        pass

    def print_summary(self, **kwargs):
        """ Print a summary of the results"""
        self.summary(**kwargs)
        if self.bestPars is not None:
            self.logger.info("MAP solution:")
            #print("MAP Solution: ")
            for i in range(self.npars):
                self.logger.info("%s  = %.5f", self.parLabels[i],self.bestPars[i]) #


        if self.res is not None:
            self.logger.info("Median and confidence intervals for parameters in order:")
            for i in range(self.npars):
                self.logger.info("%s  = %.5f \u00B1 %.5f", self.parLabels[i],self.res[i][0],self.res[i][1])


    def summary(self, **kwargs):
        self.get_credible_interval(**kwargs)
        self.bestPars = self.get_map(**kwargs)

    def _mean_and_cov(self, samples, weights):
        """Calculate the mean and covariance of the samples

        Parameters
        ----------
        samples : np.ndarray
            The samples to calculate the mean and covariance of
        weights : np.ndarray
            The weights to use for the calculation

        Returns
        -------
        mean : np.ndarray
            The mean of the samples
        cov : np.ndarray
            The covariance of the samples
        """
        mean = np.average(samples, axis=0, weights=weights)

        # Compute the weighted covariance.
        dx = samples - mean
        wsum = np.sum(weights)
        w2sum = np.sum(weights**2)
        cov = wsum / (wsum**2 - w2sum) * np.einsum('i,ij,ik', weights, dx, dx)
        return mean, cov

    def get_credible_interval(self, **kwargs):
        samples, weights = swyft.utils.get_weighted_samples(self.predictions)
        self.mean, self.cov = self._mean_and_cov(samples, weights)
        self.res = np.array([[self.mean[i], np.diag(self.cov)[i]] 
                             for i in range(self.npars)])


    def get_map(self, **kwargs):
        samples, weights = swyft.utils.get_weighted_samples(self.predictions)

        #the map is the sample with the highest weight
        map_index = np.argmax(weights)
        self.map = samples[map_index]
        return self.map

    def plot_corner(self, **kwargs):
        # We will make plots of the 2-d marginals using the swyft package
        # First we will plot all the parameters
        swyft.plot_corner(self.predictions, self.parLabels, **kwargs)
        fig = plt.gcf()
        fig.savefig(f"{self.name}_corner_all.png")

        # Then we will make separate plots of the physical parameters
        swyft.plot_corner(self.predictions, self.parLabels[:self.nparsMod], **kwargs)
        fig = plt.gcf()
        fig.savefig(f"{self.name}_corner_model.png")

        # And finally we will make separate plots of the noise parameters
        swyft.plot_corner(self.predictions, self.parLabels[self.nparsMod:], **kwargs)
        fig = plt.gcf()
        fig.savefig(f"{self.name}_corner_data.png")

    def plot_marginals(self, **kwargs):
        # We will make plots of the marginals using the swyft package
        # First we will plot all the parameters
        fig = swyft.plot_posterior(self.predictions, **kwargs)
        fig.savefig(f"{self.name}_marginals_all.png")
        plt.close(fig)

        # Then we will make separate plots of the physical parameters
        fig = swyft.plot_posterior(self.predictions,
                                   parnames=self.parLabels[:self.nparsMod],
                                   **kwargs)
        fig.savefig(f"{self.name}_marginals_model.png")
        plt.close(fig)

        # And finally we will make separate plots of the noise parameters
        fig = swyft.plot_posterior(self.predictions,
                                   parnames=self.parLabels[self.nparsMod:],
                                   **kwargs)
        fig.savefig(f"{self.name}_marginals_data.png")
        plt.close(fig)

    def plot_posteriorpredictive(self, **kwargs):
        pass
