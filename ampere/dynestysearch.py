from __future__ import print_function

import numpy as np
#import dynesty
from .basesearch import BaseSearch
from inspect import signature
#from .data import Photometry, Spectrum
from dynesty import NestedSampler, DynamicNestedSampler
from dynesty import plotting as dyplot
from dynesty import utils as dyfunc


class DynestySearch(BaseSearch):
    """
    A class to use Dynesty to explore parameters space with nested sampling
    """

    def __init__(self, #Many of the NestedSampler options are exposed here
                 prior_transform = None,
                 nlive = 1500,
                 bound = 'multi',
                 sample = 'auto',
                 dataset = None,
                 model = None,
                 update_interval=None, first_update=None,
                 queue_size=None, pool=None, use_pool=None,
                 enlarge=None, bootstrap=0, vol_dec=0.5,
                 vol_check=2.0, walks=25, facc=0.5, slices=5, **kwargs
                 ):
        self.model=model
        self.dataSet = dataset
        if prior_transform is not None: #This should probably be removed! We don't want the user change the ptform at this point, but by changing it for the Model or Data individually
            self.prior_transform = prior_transform
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
        #print(self.nparsMod, len(sig.parameters), sig.parameters)
        #self.nparsData = #np.zeros(len(self.dataSet))
        #self.npars = something # total number of parameters
        #self.nparsMod = something #number of parameters for the model
        self.nparsData = [data.npars for data in self.dataSet] #number of parameters to be passed into each set of data
        self.npars = np.int(self.nparsMod + np.sum(self.nparsData))
        print(self.npars, self.nparsMod, self.nparsData)
        self.sampler = NestedSampler(self.lnlike, self.prior_transform, self.npars,
                                     nlive=nlive, bound=bound, sample=sample,
                                     update_interval = update_intervale,
                                     first_update = first_update, queue_size = queue_size,
                                     pool = pool, use_pool = use_pool, enlarge = enlarge,
                                     bootstrap = bootstrap, vol_dec = vol_dec,
                                     vol_check = vol_check, walks = walks, facc = facc,
                                     slices = slices, **kwargs
        )


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

    def optimise(self, dlogz=None, maxiter=None, maxcall=None,
                 logl_max=inf, add_live=True,
                 print_progress=True, print_func=None,
                 save_bounds=True **kwargs
                 ):
        self.sampler.run_nested(maxiter=maxiter, maxcall=maxcall,
                 dlogz=dlogz, logl_max=logl_max, add_live=add_live,
                 print_progress=print_progress, print_func=print_func,
                 save_bounds=save_bounds **kwargs)
        self.results = self.sampler.results

    def postProcess(self, maxiter=None, maxcall = None, dlogz = None, **kwargs):
        """ Some simple post-processing and plotting """

        """ first, a summary plot """
        fig, axes = dyplot.runplot(self.results) #Summary plot
        fig.savefig("summary.png")


        """ next, a corner plot """
        fg, ax = dyplot.cornerplot(self.results, color='red', #truths=np.zeros(ndim), truth_color='black',
                                   show_titles=True,
                                   quantiles=None, max_n_ticks=5)
        fg.savefig("corner.png")

        """ Finally, the trace plot """
        fig, axes = dyplot.traceplot(self.results, #truths=np.zeros(ndim),
                                     show_titles=True,
                                     trace_cmap='viridis', connect=True,
                                     connect_highlight=range(5))
        fig.savefig("trace.png")
        """ with the plotting completed, we now compute estimates of the posterior """
        samples, weights = self.results.samples, np.exp(self.results.logwt - self.results.logz[-1])
        self.mean, self.cov = dyfunc.mean_and_cov(samples, weights)
        print(self.mean)
        print(self.cov)
        
        pass

class DynestyDynamicSearch(DynestySearch): #I think this can inheret almost everything except __init__ from the Static sampler
    """
    A class to use Dynesty to explore parameters space with dynamic nested sampling
    """

    def __init__(self, 
                 ):
        pass
    



                 
