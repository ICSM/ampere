from __future__ import print_function

import logging
import numpy as np
from inspect import signature
from .basesearch import BaseSearch
from ..logger import Logger

class MCMCSampler(BaseSearch, Logger):

    def __init__(self, model = None, data= None, verbose = False,
                 parameter_labels = None,
                 **kwargs):
        self.model = model
        self.dataSet = data

        self.setup_logging(verbose=verbose) #For now we will exclusively use default logging settings, this will be modified once logging is tested.
        logging.info("Welcome to ampere")
        logging.info("Setting up your inference problem:")
        logging.info("You have %s items in your dataset", str(len(data)))

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
        logging.info("This model has %s parameters.", str(self.nparsMod))
        logging.info("There are also %s parameters for the noise model", str(self.npars - self.nparsMod))
        logging.info("Hence, there are a total of %s parameters to sample", str(self.npars))
        logging.info("The parameter names are:")
        for l in self.parLabels:
            logging.info(str(l))
            


class EnsembleSampler(MCMCSampler):
    """ A base-class for emcee-like samplers (e.g. emcee, zeus, ptemcee)

    """

    def __init__(self, model = None, data= None, verbose = False,
                 parameter_labels = None,
                 nwalkers = None, 
                 **kwargs):

        super().__init__(model = model, data= data, verbose = verbose,
                         parameter_labels = parameter_labels, **kwargs)

        self.nwalkers = nwalkers
        logging.info("This problem will be sampled with %s walkers", str(self.nwalkers))

    def rebuildSampler(self, nwalkers = None, model = None,
                       data = None, lnprior = None,
                       labels = None, 
                       **kwargs):

        logging.info("Rebuilding your sampler")
        if np.any([model, data]):
            if model is not None:
                logging.info("Your model is being updated")
                self.model=model
                try:
                    self.nparsMod = self.model.npars
                except:
                    sig = signature(model.__call__)
                    self.nparsMod = len(sig.parameters) - 1 #Always subtract **kwargs from the parameters, but don't need to worry about self once it is bound to an instance
            if data is not None:
                self.dataSet = data
                logging.info("Your data is being updated")
                self.nparsData = [data.npars for data in self.dataSet] #number of parameters to be passed into each set of data

            self.npars = np.int(self.nparsMod + np.sum(self.nparsData))
        
        if lnprior is not None:
            self.lnprior = lnprior
        #print(self.npars, self.nparsMod, self.nparsData)
        if nwalkers is not None:
                self.nwalkers=nwalkers
                logging.info("The number of walkers is now %s", str(self.nwalkers))
