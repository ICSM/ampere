from __future__ import print_function

import numpy as np
from inspect import signature
from .basesearch import BaseSearch


class MCMCSampler(Basesearch):

    def __init__(self, model = None, data= None, verbose = False,
                 parameter_labels = None,
                 **kwargs):
        self.model = model
        self.dataSet = data

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
        if self.verbose:
            print("This model has ", self.nparsMod," parameters.")
            print("There are also ", self.npars - self.nparsMod, " parameters for the noise model")
            print("Hence, there are a total of ", self.npars, " parameters to sample")
            print("The parameter names are:")
            for l in self.parLabels:
                print(l)
            


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
