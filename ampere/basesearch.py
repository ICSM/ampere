from __future__ import print_function

import numpy as np
import pickle


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

    def lnprior(self, theta, **kwargs):
        return self.model.prior(theta)
        #raise NotImplementedError()

    def lnlike(self, theta, **kwargs):
        model = self.model(*theta[:self.nparsMod])
        l=np.array([])
        i = self.nparsMod
        for data in self.dataSet:
            #print(data)
            lnew = data.lnlike(theta[i:i+data.npars],self.model)
            i+=data.npars
            l = np.r_[l,lnew]
        return np.sum(l)
            
        #raise NotImplementedError()

    def lnprob(self, theta):#, dataSet, **kwargs):
        p = self.lnprior(theta)#[:self.nparsMod])
        if p == -np.inf:
            return p
        return p + self.lnlike(theta)

    def sampler(self, **kwargs):
        raise NotImplementedError()

    def save(self, filename, pickle_it=True,**kwargs):
        ''' A method to save the object to a file. For the moment, this only supports pickling '''
        if not pickle_it:
            raise NotImplementedError('Only pickling is supported at this time')
        with open(filename, 'wb') as f:
            pickle.dump(self.__dict__, f)

    
