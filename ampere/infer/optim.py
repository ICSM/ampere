from __future__ import print_function

import numpy as np
import pickle
from datetime import datetime
from ..logger import Logger
from .basesearch import BaseSearch
from inspect import signature


class OptimBase(BaseSearch, Logger):
    """ A base class for MAP inferece algorithms.
    
    This class will be used to derive subclasses that wrap specific
    optimisation algorithms.
    """

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

    def __init__(self, 
                 model=None,
                 data=None,
                 verbose=False,
                 parameter_labels=None,
                 name='',
                 namestyle="full",
                 **kwargs):
        self.model = model
        self.dataSet = data
        self.namestyle = namestyle
        self.name = name
        self.verbose = verbose

        self.setup_logging(verbose=verbose)
        self.logger.info("Welcome to ampere")
        self.logger.info("Setting up your inference problem:")
        self.logger.info("You are using %s", self._inference_method)
        self.logger.info("You have %s items in your dataset", str(len(data)))

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

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name=''):
        try:
            modelname = self.model.name
        except AttributeError:
            modelname = self.model.__class__.__name__
        if self.namestyle == 'full':
            self._name = ("ampere_"+str(datetime.now()
                                       ).replace(' ',
                                                 '_').replace(":", "-")[:-7]
                          + "_" + modelname + "_" + name
                          )
        elif self.namestyle == "model":
            self._name = f"ampere_{modelname}_{name}"
        elif self.namestyle == "short":
            self._name = name
        elif self.namestyle == "stamp":
            self._name = ("ampere_"+str(datetime.now()).replace(' ',
                                                               '_').replace(":",
                                                                            "-")[:-7]
                                  + "_" + name
                         )


class ScipyMinOpt(OptimBase):
    """A wrapper for scipy.optimize.minimize

    This class wraps the scipy.optimize.minimize function to provide
    a common interface for optimisation algorithms. It is intended to
    provide an easy starting point for quickly attempting different
    optimisation algorithms. It can also be used to set initial guesses
    for the more sophisticated inference algorithms provided by Ampere.
    Since it is only an optimiser, not a sampler, it gives rather poor
    estimates of the uncertainty in the parameters. It is therefore
    not advisable to trust the results of this algorithm for publication.
    In addition, this algorithm performs rather poorly when the posterior
    is multimodal or complex, and other approaches are likely to give better
    results.

    Parameters
    ----------
    OptimBase : _type_
        _description_
    """

    def __init__(self,
                 **kwargs):
        
        pass

    def optimize(self):
        pass

    def postProcess(self):
        pass

    pass


class ScipyBasinOpt(OptimBase):
    """A wrapper for scipy.optimize.basinhopping

    _extended_summary_

    Parameters
    ----------
    OptimBase : _type_
        _description_
    """

    pass

class ScipyDE(OptimBase):
    """A wrapper for scipy.optimize.differential_evolution

    _extended_summary_

    Parameters
    ----------
    OptimBase : _type_
        _description_
    """

    pass


class ScipyDualAnneal(OptimBase):
    """A wrapper for scipy.optimize.dual_annealing

    _extended_summary_

    Parameters
    ----------
    OptimBase : _type_
        _description_
    """

    pass


class AxOpt(OptimBase):
    """A wrapper for Ax

    _extended_summary_

    Parameters
    ----------
    OptimBase : _type_
        _description_
    """

    pass

class AxBO(OptimBase):
    """A wrapper for Bayesian Optimisation with Ax

    _extended_summary_

    Parameters
    ----------
    OptimBase : _type_
        _description_
    """

    pass
