from __future__ import print_function

import logging
import numpy as np
from inspect import signature
from .basesearch import BaseSearch
from ..logger import Logger

class MCMCSampler(BaseSearch, Logger):

    _inference_method = "MCMC"

    def __init__(self, model = None, data= None, verbose = False,
                 parameter_labels = None,
                 **kwargs):
        self.model = model
        self.dataSet = data

        self.setup_logging(verbose=verbose) #For now we will exclusively use default logging settings, this will be modified once logging is tested.
        self.logger.info("Welcome to ampere")
        self.logger.info("Setting up your inference problem:")
        self.logger.info("You are using %s", self._inference_method)
        self.logger.info("You have %s items in your dataset", str(len(data)))

        ''' now do some introspection on the various bits of model to 
        understand how many parameters there are for each compponent '''
        try:
            self.nparsMod = self.model.npars
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
            


class EnsembleSampler(MCMCSampler):
    """ A base-class for emcee-like samplers (e.g. emcee, zeus, ptemcee)

    """

    _inference_method = "Ensemble MCMC"

    def __init__(self, model = None, data= None, verbose = False,
                 parameter_labels = None,
                 nwalkers = None, 
                 **kwargs):

        super().__init__(model = model, data= data, verbose = verbose,
                         parameter_labels = parameter_labels, **kwargs)

        self.nwalkers = nwalkers
        self.logger.info("This problem will be sampled with %d walkers", self.nwalkers)

    def rebuildSampler(self, nwalkers = None, model = None,
                       data = None, lnprior = None,
                       labels = None, 
                       **kwargs):

        self.logger.info("Rebuilding your sampler")
        if np.any([model, data]):
            if model is not None:
                self.logger.info("Your model is being updated")
                self.model=model
                try:
                    self.nparsMod = self.model.npars
                except:
                    sig = signature(model.__call__)
                    self.nparsMod = len(sig.parameters) - 1 #Always subtract **kwargs from the parameters, but don't need to worry about self once it is bound to an instance
            if data is not None:
                self.dataSet = data
                self.logger.info("Your data is being updated")
                self.nparsData = [data.npars for data in self.dataSet] #number of parameters to be passed into each set of data

            self.npars = np.int(self.nparsMod + np.sum(self.nparsData))
        
        if lnprior is not None:
            self.lnprior = lnprior
        #print(self.npars, self.nparsMod, self.nparsData)
        if nwalkers is not None:
                self.nwalkers=nwalkers
                self.logger.info("The number of walkers is now %d", self.nwalkers)
