from __future__ import print_function

import numpy as np



class BaseSearch(object):
    """
    A base class for parameter search algorithms
    """
    
    def __init__(self, model=None, data=None, **kwargs):
        '''
        Initialise the optimiser with a model and some data, then use introspection on both to establish a few important numbers 
        '''
        raise NotImplementedError()

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()

    def lnprior(self, theta,**kwargs):
        return self.model.prior(theta)
        #raise NotImplementedError()

    def lnlike(self, theta,**kwargs):
        model = self.model(theta)
        l=np.array([])
        for data in self.datasets:
            l = np.r_[l,np.sum(data.lnlike(model))]
        return np.sum(l)
            
        #raise NotImplementedError()

    def lnprob(self, **kwargs):
        p = self.lnprior()
        if p == -np.inf:
            return p
        raise p + self.lnlike()

    def sampler(self, **kwargs):
        raise NotImplementedError()

    
