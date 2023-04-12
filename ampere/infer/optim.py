from __future__ import print_function

import numpy as np
import pickle
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

    pass


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
