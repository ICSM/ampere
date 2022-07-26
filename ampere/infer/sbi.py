"""
A file for simulation-based inference approaches, starting with the sbi package

"""

from __future__ import print_function

import numpy as np
import pickle
import logging
from .basesearch import BaseSearch
from ..logger import Logger
from ..models.results import ModelResults

class LFIBase(BaseSearch, Logger):
    """ A base class for simulation-based (or likelihood-free) inference 
        approaches
    
    
    """

    _inference_method="LFI"


    def __init__(self, model=None, data=None, verbose = False,
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
        logging.info("This model has %d parameters.", self.nparsMod)
        logging.info("There are also %d parameters for the noise model", self.npars - self.nparsMod)
        logging.info("Hence, there are a total of %d parameters to sample", self.npars)
        logging.info("The parameter names are:")
        for l in self.parLabels:
            logging.info("%s", l)

        #Now we should check whether the number of parameters matches with the parameter labels!
        if len(self.parLabels) != self.npars:
            logging.critical("You have %d free parameters but %d parameter labels", self.npars, len(self.parLabels))
            logging.critical("Please check the number of parameters and labels")
            raise ValueError("Mismatch between number of parameters and labels")

    def simulate(self, theta): #, theta_noise):
        """A method to simulate an observation """

        #First we call the model
        result = self.model(theta[:self.nparsMod])
        if not isinstance(result, ModelResults):
            result = ModelResults(**result) #Try unpacking the results here if the user didn't already define their model with it

        #Then we make each Data object produce a sample synthetic observation
        sims = []
        i=self.nparsMod
        for d in self.dataset:
            sim = d.simulate(theta[i:i+data.npars],result)
            sims.append(sim)
            i+=data.npars

        return sims #Now we return the simulated data. Inidividual implementations should modify this to accomodate their idiosyncracies


    def simulate_vector(self, theta):
        """ A method to compute a batch of simulations """
        pass

            
    
