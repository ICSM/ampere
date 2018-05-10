from __future__ import print_function

import numpy as np



class BaseSearch(object):
    """
    A base class for parameter search algorithms
    """
    
    def __init__(self, **kwargs):
        raise NotImplementedError()

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()

    def lnprior(self, **kwargs):
        raise NotImplementedError()

    def lnlike(self, **kwargs):
        raise NotImplementedError()

    def lnprob(self, **kwargs):
        p = self.lnprior()
        if p == -np.inf:
            return p
        raise p + self.lnlike()

    def sampler(self, **kwargs):
        raise NotImplementedError()

    
