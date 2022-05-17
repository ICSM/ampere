from __future__ import print_function

import numpy as np
from inspect import signature
import logging
from ..logger import Logger
from .basesearch import BaseSearch


class BaseNestedSampler(BaseSearch, Logger):
    _inference_method = "Nested Sampling"
    
    def __init__(self, model = None, data= None, verbose = False,
                 parameter_labels = None,
                 **kwargs):
        self.model = model
        self.dataSet = data

        self.setup_logging(verbose=verbose) #For now we will exclusively use default logging settings, this will be modified once logging is tested.
        logging.info("Welcome to ampere")
        logging.info("Setting up your inference problem:")
        logging.info("You are using %s", self._inference_method)
        logging.info("You have %s items in your dataset", str(len(data)))

        ''' now do some introspection on the various bits of model to 
        understand how many parameters there are for each compponent '''
        try: #First we try to see if the number of parameters is documented specifically for methods which use a prior transform
            self.nparsMod = self.model.npars_ptform
        except AttributeError:
            try: #If not, we assume that the number of parameters
                self.nparsMod = self.model.npars
            except AttributeError:
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

        self.verbose = verbose
        #if self.verbose:
        logging.info("This model has %d parameters.", self.nparsMod)
        logging.info("There are also %d parameters for the noise model", self.npars - self.nparsMod)
        logging.info("Hence, there are a total of %d parameters to sample", self.npars)
        logging.info("The parameter names are:")
        for l in self.parLabels:
            logging.info("%s", l)


    #This method is now defined in Basesearch
    def prior_transform(self, u, **kwargs):
        """
        We delegate the prior transforms to the models and the data
        """
        #print(u)
        theta = np.zeros_like(u)
        theta[:self.nparsMod] = self.model.prior_transform(u[:self.nparsMod])
        i = self.nparsMod
        for data in self.dataSet:
            theta[i:i+data.npars] = data.prior_transform(u[i:i+data.npars])
            i+=data.npars
        return theta
