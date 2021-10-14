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
                 data = None,
                 model = None,
                 update_interval=None, first_update=None,
                 queue_size=None, pool=None, use_pool=None,
                 enlarge=None, bootstrap=0, vol_dec=0.5,
                 vol_check=2.0, walks=25, facc=0.5, slices=5, **kwargs
                 ):
        self.model=model
        self.dataSet = data
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
                                     update_interval = update_interval,
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

    def optimise(self, dlogz=50., maxiter=None, maxcall=None,
                 logl_max=np.inf, add_live=True,
                 print_progress=True, print_func=None,
                 save_bounds=True, **kwargs
                 ):
        self.sampler.run_nested(maxiter=maxiter, maxcall=maxcall,
                                dlogz=dlogz, logl_max=logl_max,
                                add_live=add_live,print_progress=print_progress,
                                print_func=print_func,
                                save_bounds=save_bounds, **kwargs)
        self.results = self.sampler.results

    def plot_summary(self, plotfile="summary.png"):
        fig, axes = dyplot.runplot(self.results) #Summary plot
        fig.savefig(plotfile)

    def plot_corner(self, plotfile="corner.png"):
        """ next, a corner plot """
        fg, ax = dyplot.cornerplot(self.results, color='red', #truths=np.zeros(ndim), truth_color='black',
                                   show_titles=True,
                                   quantiles=None, max_n_ticks=5)
        fg.savefig(plotfile)

    def plot_trace(self, plotfile="trace.png"):
        fig, axes = dyplot.traceplot(self.results, #truths=np.zeros(ndim),
                                     show_titles=True,
                                     trace_cmap='viridis', connect=True,
                                     connect_highlight=range(5))
        fig.savefig(plotfile)


    def plot_posteriorpredictive(self, n_samples = 1000, plotfile="posteriorpredictive.png", logx = False, logy = False, alpha = 0.05):
        ''' Generate the posterior-predictive plots so that the suitability of the model for the data can be inspected. 
        '''

        #This function needs to get vastly more complex in future, somehow. For now, it assumes we want to plot SEDs/spectra

        #This will probably require some small modifications to Data objects to make this as easy as possible.
        #It should include plotting multiple realisations of the data given the noise model (if appropriate)
        #as well as realisations from the model and the raw data
        #Plot MAP solution model?
        #And model with the median of marginal posterior parameters
        #Must be plotted in reverse order: Model realisations first (i.e. lowest zorder), then data realisations, then data (use zorder keyword to get things in the right place, lowest first as highest end up on top)
        #Also produce individual posterior predictive plots for each Data object.

        #Reweight the samples
        samples = self.results.samples.T
        weights = np.exp(self.results.logwt - self.results.logz[-1])
        samples_unif = resample_equal(samples.T, weights)

        #First set up the plot
        fig = plt.figure()
        ax = fig.add_subplot(111)

        for s in samples_unif[np.random.randint(len(samples_unif), size=n_post_samples)]:
            optimizer.model(s)
            ax.plot(optimizer.model.wavelengths,optimizer.model.modelFlux, '-', color='k', alpha=alpha)

        for d in self.dataSet:
            #We need to create some properties in each Data class that contain some plotting info
            #e.g. use plot for spectra with fill_between for uncertainties
            #use errorbar for photometry with circles
            #for now make bad plots with lines everywhere
            ax.plot(d.wavelength, d.value, '-',color='blue')
        
        fig.savefig(plotfile)
        pass

    def postProcess(self, maxiter=None, maxcall = None, dlogz = None, **kwargs):
        """ Some simple post-processing and plotting """

        """ first, a summary plot """
        try:
            self.plot_summary()
            #fig, axes = dyplot.runplot(self.results) #Summary plot
            #fig.savefig("summary.png")
        except ValueError as e:
            print("summary failed with error",e)
            print("skipping summary plot and moving on to corner plot")


        """ next, a corner plot """
        self.plot_corner()
        #fg, ax = dyplot.cornerplot(self.results, color='red', #truths=np.zeros(ndim), truth_color='black',
        #                           show_titles=True,
        #                           quantiles=None, max_n_ticks=5)
        #fg.savefig("corner.png")

        """ Finally, the trace plot """
        self.plot_trace()
        #fig, axes = dyplot.traceplot(self.results, #truths=np.zeros(ndim),
        #                             show_titles=True,
        #                             trace_cmap='viridis', connect=True,
        #                             connect_highlight=range(5))
        #fig.savefig("trace.png")
        """ with the plotting completed, we now compute estimates of the posterior """
        #This should also be extracted to a standalone method which we can call from here
        samples, weights = self.results.samples, np.exp(self.results.logwt - self.results.logz[-1])
        self.mean, self.cov = dyfunc.mean_and_cov(samples, weights)
        print("Posterior means of the parameters: ",self.mean)
        print("1-sigma confidence intervals of the parameters marginalising over all other parameters:")
        print(np.sqrt(np.diag(self.cov)))
        print("Posterior covariances of the parameters: ", self.cov)


        self.plot_posteriorpredictive()
        
        pass

class DynestyDynamicSearch(DynestySearch): #I think this can inheret almost everything except __init__ from the Static sampler
    """
    A class to use Dynesty to explore parameters space with dynamic nested sampling
    """

    def __init__(self, 
                 ):
        pass
    



                 
