from __future__ import print_function

import numpy as np
from inspect import signature
from datetime import datetime
from .basesearch import BaseSearch
from ..logger import Logger


class MCMCSampler(BaseSearch, Logger):

    """
    A base class for MCMC Sampling inference approaches.

    Intended to be subclassed for package-specific implementations.

    Parameters
    ----------
    model : ampere.model.Model
        Model whose posterior will be inferred
    data : list or iterable of ampere.data.Data
        the dataset to fit
    verbose : bool, optional
        If True, print verbose output.
    parameter_labels : list of str, optional
        List of strings containing names of each parameter.
    name : str, optional
        Name of the sampler.
    namestyle : str, optional
        String specifying style for naming the sampler.
        Options: 'full', 'short', 'stamp', 'model'.
    """

    _inference_method = "MCMC"

    def __init__(self, model=None, data=None, verbose=False,
                 parameter_labels=None, name='', namestyle="full",
                 **kwargs):
        self.model = model
        self.dataSet = data

        try:
            modelname = model.name
        except AttributeError:
            modelname = model.__class__.__name__
        if namestyle == 'full':
            self.name = ("ampere_"+str(datetime.now()
                                       ).replace(' ',
                                                 '_').replace(":", "-")[:-7]
                         + "_" + modelname + "_" + name
                         )
        elif namestyle == "short":
            self.name = name
        elif namestyle == "stamp":
            self.name = ("ampere_"+str(datetime.now()).replace(' ',
                                                               '_').replace(":",
                                                                            "-")[:-7]
                                  + "_" + name
                         )
        elif namestyle == "model":
            self.name = f"ampere_{modelname}_{name}"

        # For now we will exclusively use default logging settings, this will
        #  be modified once logging is tested.
        self.setup_logging(verbose=verbose)
        self.logger.info("Welcome to ampere")
        self.logger.info("Setting up your inference problem:")
        self.logger.info("You are using %s", self._inference_method)
        self.logger.info("You have %s items in your dataset", str(len(data)))

        # now do some introspection on the various bits of model to
        # understand how many parameters there are for each compponent
        try:
            self.nparsMod = self.model.npars
        except AttributeError:
            sig = signature(model.__call__)
            # Always subtract **kwargs from the parameters, but don't need to
            # worry about self once it is bound to an instance
            self.nparsMod = len(sig.parameters) - 1
        # number of parameters to be passed into each set of data
        self.nparsData = [data.npars for data in self.dataSet]
        self.npars = int(self.nparsMod + np.sum(self.nparsData))

        if parameter_labels is None:
            # The user hasn't specified parameter labels, let's see if the
            # models and data have instead
            self.parLabels = [] # first set parLabels to an empty list, to prevent endlessly 
            # extending the same list if re-creating the same object.
            try:  # First the model parameters
                self.parLabels = self.model.parLabels
            except AttributeError:
                self.parLabels = [f'x{str(i)}' for i in range(self.nparsMod)]
            i = self.nparsMod
            for data in self.dataSet:
                try:
                    self.parLabels.extend(data.parLabels)
                except AttributeError:
                    self.parLabels.extend([f'x{str(i)}'
                                           for i in range(i, i+data.npars)])
                finally:
                    i += data.npars
        else:  # User isn't really supposed to use this interface, however...
            self.parLabels = parameter_labels

        self.verbose = verbose
        # if self.verbose:
        self.logger.info("This model has %d parameters.", self.nparsMod)
        self.logger.info("There are also %d parameters for the noise model",
                         self.npars - self.nparsMod)
        self.logger.info("Hence, there are a total of %d parameters to sample",
                         self.npars)
        self.logger.info("The parameter names are:")
        for label in self.parLabels:
            self.logger.info("%s", label)

        # Now we should check whether the number of parameters matches with
        # the parameter labels!
        if len(self.parLabels) != self.npars:
            self.logger.critical("You have %d free parameters but %d parameter labels",
                                 self.npars, len(self.parLabels))
            self.logger.critical("Please check the number of parameters and labels")
            raise ValueError("Mismatch between number of parameters and labels")


class EnsembleSampler(MCMCSampler):
    """ A base-class for emcee-like samplers (e.g. emcee, zeus, ptemcee)

     Intended to be subclassed for package-specific implementations.

    Parameters
    ----------
    model : ampere.model.Model
        Model which will be used
    data : list or iterable of ampere.data.Data
        the dataset to fit
    nwalkers : int
        number of walkers to use in the ensemble
    verbose : bool, optional
        If True, log verbosely.
    parameter_labels : list of str, optional
        List of strings containing names of each parameter.
    name : str, optional
        Name of the sampler.
    namestyle : str, optional
        String specifying style for naming the sampler.
        Options: 'full', 'short', 'stamp', 'model'.
    """

    _inference_method = "Ensemble MCMC"

    def __init__(self, model=None, data=None, verbose=False,
                 parameter_labels=None, name='', namestyle="full",
                 nwalkers=None,
                 **kwargs):

        super().__init__(model=model, data=data, verbose=verbose,
                         parameter_labels=parameter_labels, name=name,
                         namestyle=namestyle, **kwargs)

        self.nwalkers = nwalkers
        self.logger.info("This problem will be sampled with %d walkers", self.nwalkers)

    def rebuildSampler(self, nwalkers=None, model=None,
                       data=None, lnprior=None,
                       labels=None,
                       **kwargs):

        """
        Rebuilds the MCMC sampler with the specified parameters.

        Parameters
        ----------
        nwalkers : int, optional
            The number of walkers to use in the MCMC sampling. If not provided,
            the current number of walkers will be used.
        model : object, optional
            A new model to use in the MCMC sampling. If not provided, the
            current model will be used.
        data : object, optional
            New data to use in the MCMC sampling. If not provided, the current
            data will be used.
        lnprior : function, optional
            A new prior distribution to use in the MCMC sampling. If not
            provided, the current prior will be used.
        labels : list of str, optional
            A list of parameter labels to use in post-processing. If not
            provided, the current labels will be used.
        **kwargs : additional keyword arguments
            Additional keyword arguments to pass to the sampler.
        """

        self.logger.info("Rebuilding your sampler")
        if np.any([model, data]):
            if model is not None:
                self.logger.info("Your model is being updated")
                self.model = model
                try:
                    self.nparsMod = self.model.npars
                except AttributeError:
                    sig = signature(model.__call__)
                    # Always subtract **kwargs from the parameters, but don't
                    # need to worry about self once it is bound to an instance
                    self.nparsMod = len(sig.parameters) - 1
            if data is not None:
                self.dataSet = data
                self.logger.info("Your data is being updated")
                # number of parameters to be passed into each set of data
                self.nparsData = [data.npars for data in self.dataSet]

            self.npars = int(self.nparsMod + np.sum(self.nparsData))

        if lnprior is not None:
            self.lnprior = lnprior
        # print(self.npars, self.nparsMod, self.nparsData)
        if nwalkers is not None:
            self.nwalkers = nwalkers
            self.logger.info("The number of walkers is now %d", self.nwalkers)
