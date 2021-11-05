from __future__ import print_function

import numpy as np
import pickle


class BaseSearch(object):
    """ A base class for parameter search algorithms.

    This is the base class from which all other search classes will be derived.
    It doesn't implement much, other than the prior, posterior and likelihood 
    functions. It also implements a general method for serialising (pickling)
    the search objects, so that they can be saved for later or to make parallel
    execution easier. __getstate__ and __setstate__ should be implemented as I 
    expect non of ampere's classes will be easily pickle/dill-able. The methods 
    defined here only need to be re-implemented if you expect to change their 
    behaviour significantly, and the NotImplemented methods must be defined in 
    subclasses to make sure everything is useful.
    
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
    
    def lnprior(self, theta, **kwargs):
        """ Calculate the probability of the model (the prior)
        
        This just calls the prior for the model, but should also call the priors
        for the data too.
        """
        return self.model.prior(theta)
        #raise NotImplementedError()

    def lnlike(self, theta, **kwargs):
        """ Calculate the likelihood of the model given the data

        This routine takes a simple approach to calculating the log-likelihood
        of the model given the data. It calculates the model, then loops over 
        all items in the dataset and calculates their likelihoods.

        """
        model = self.model(*theta[:self.nparsMod])
        l=np.array([])
        i = self.nparsMod
        for data in self.dataSet:
            lnew = data.lnlike(theta[i:i+data.npars],self.model)
            i+=data.npars
            l = np.r_[l,lnew]
        return np.sum(l)
            
        #raise NotImplementedError()

    def lnprob(self, theta):#, dataSet, **kwargs):
        """ Calculate the probability of the model given the data
        
        This just calls the prior, followed by the likelihood, and returns the 
        sum of the two.
        """
        
        p = self.lnprior(theta)#[:self.nparsMod])
        if p == -np.inf:
            return p
        return p + self.lnlike(theta)

    def lnprob_vector(self, theta, pool):
        ''' This method is a drop-in replacement for the above method but where theta is a vector of arguments. 

        It assumes you want to parallelise execution of models across a pool
        '''

        lnprobs = pool.map(lnprob_pool, theta, worker_init=self.init_workers)
        return lnprobs
        
        pass

    def mpire_setup(self):
        ''' Do some setup steps to initialise the WorkerPool and have it ready for parallel processing. 
        '''
        from mpire import WorkerPool
        
        pool = WorkerPool(njobs = self.njobs,
                          use_worker_state=True,
                          keep_alive=True,
                          use_dill=self.use_dill
        )

        self.init_workers = create_worker_init(self)
        return pool

    def sampler(self, **kwargs):
        raise NotImplementedError()

    def save(self, filename, pickle_it=True,**kwargs):
        ''' A method to save the object to a file. For the moment, this only supports pickling '''
        if not pickle_it:
            raise NotImplementedError('Only pickling is supported at this time')
        with open(filename, 'wb') as f:
            pickle.dump(self.__dict__, f)

def create_worker_init(searcher):    
    def init_workers(worker_state):
        worker_state['searcher'] = searcher
    return init_workers


def lnlike_pool(worker_state, theta):
    return worker_state['searcher'].lnlike(theta)

def lnprob_pool(worker_state, theta):
    return worker_state['searcher'].lnprob(theta)


    
