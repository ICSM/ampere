from __future__ import print_function

import numpy as np
#import dynesty
from .basesearch import BaseSearch
from inspect import signature
from .data import Photometry, Spectrum
from dynesty import NestedSampler, DynamicNestedSampler
from dynesty import plotting as dyplot
from dynesty import utils as dyfunc
import matplotlib.pyplot as plt


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
        self.parLabels = ['x'+str(i) for i in range(self.npars)] #Parameter for parameter names (labels) to associate with output in post processing
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
        fg, ax = dyplot.cornerpoints(self.results,
                                     max_n_ticks=5, cmap="plasma", kde=True)
        #fg, ax = dyplot.cornerplot(self.results)#, color='red', #truths=np.zeros(ndim), truth_color='black',
                                   #show_titles=True,
                                   #quantiles=None, max_n_ticks=5)
        fg.savefig(plotfile)

    def plot_trace(self, plotfile="trace.png"):
        fig, axes = dyplot.traceplot(self.results, #truths=np.zeros(ndim),
                                     show_titles=True,
                                     trace_cmap='viridis', connect=True,
                                     connect_highlight=range(5))
        fig.savefig(plotfile)

    def plot_covmats(self):
        istart = self.nparsMod
        for i, d in enumerate(self.dataSet):
            if isinstance(d, Photometry):
                continue
            elif isinstance(d, Spectrum):
                fig, ax = plt.subplots(1,1)
                #for d in self.dataSet[1:]:
                #d=self.dataSet[1]
                #print("Using these parameters for the covariance matrix:")
                #print(self.parLabels[istart+1], self.parLabels[istart+2])
                d.cov([self.res[istart+1][0],self.res[istart+2][0]])
                #ax0.set_title('Covariance matrix')
                im = ax.imshow(np.log10(d.covMat))
                istart+=d.npars
                fig.savefig("covMat_"+str(i)+".png")


    def plot_posteriorpredictive(self, n_post_samples = 1000, plotfile="posteriorpredictive.png", logx = False, logy = False, alpha = 0.1, **kwargs):
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
        samples_unif = dyfunc.resample_equal(samples.T, weights)

        #First set up the plot
        fig = plt.figure()
        axes = fig.add_subplot(111)

        #observations
        for d in self.dataSet:
            #Use each data object's self-knowledge to plot stuff
            d.plot(ax = axes)

        for s in samples_unif[np.random.randint(len(samples_unif), size=n_post_samples)]:
            self.model(*s[:self.nparsMod])
            axes.plot(self.model.wavelength,self.model.modelFlux, '-', color='k', alpha=alpha, label='Samples', zorder=0)

            i = self.nparsMod
            for d in self.dataSet:
                if d._hasNoiseModel:
                    d.plotRealisation(s[i:i+d.npars], ax=axes)
                i+= d.npars

        #best fit model
        try:
            self.model(*self.bestPars[:self.nparsMod])
            axes.plot(self.model.wavelength,self.model.modelFlux, '-', color='k', alpha=1.0,label='MAP', zorder=8)
        except ValueError:
            print("Error in MAP solution \n Skipping MAP in plot")


        #These plots end up with too many labels for the legend, so we clobber the label information so that only one of each one is plotted
        handles, labels = plt.gca().get_legend_handles_labels()
        # labels will be the keys of the dict, handles will be values
        temp = {k:v for k,v in zip(labels, handles)}
        plt.legend(temp.values(), temp.keys(), loc='best')

        plt.tight_layout()
        fig.savefig(plotfile)
        plt.close(fig)
        plt.clf()

    def postProcess(self, maxiter=None, maxcall = None, dlogz = None, testfile=None, **kwargs):
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
        self.print_summary(outfile=textfile)

        self.plot_covmats()

        self.plot_posteriorpredictive()
        
        pass


    def print_summary(self,outfile=None):
        ''' Standalone method to compute and present estimates of the posterior.
        '''
        
        #Calculate
        samples, weights = self.results.samples, np.exp(self.results.logwt - self.results.logz[-1])
        self.mean, self.cov = dyfunc.mean_and_cov(samples, weights)
        self.res = np.array([[self.mean[i], self.cov[i]] for i in range(self.npars)])
        
        #Present
        print("Posterior means and 1-sigma confidence intervals of the parameters marginalising over all other parameters: ")
        for i in range(self.npars):
            print("{0}  = {1:.5f} +/- {2:.5f}".format(
                self.parLabels[i],self.mean[i],np.sqrt(np.diag(self.cov)[i]))
                  )
        #print(np.sqrt(np.diag(self.cov))
        #print("Posterior covariances of the parameters: ", self.cov)

        #Now produce ML and MAP solution
        self.bestPars = self.results.samples[-1]
        print("MAP Solution: ")
        for i in range(self.npars):
            print("{0}  = {1:.5f}".format(self.parLabels[i],self.bestPars[i])) #
        print("with Posterior probability ln(P*) = {0:.5f}".format(self.results.logwt[-1]))
        print("and likelihood ln(L*) = {0:.5f}".format(self.results.logl[-1]))



        #print(self.results)
        #help(self.results)

        #Then print out evidence, uncertainty, and estimate of remaining evidence
        print("Model evidence ln(Z) = {0:.5f} +/- {1:.5f}".format(self.results.logz[-1], self.results.logzerr[-1]))
        #print("Estimated remaining evidence dln(Z) = ",

        if outfile is not None:
            with open(outfile, 'w') as f:
                f.write("Posterior means and 1-sigma confidence intervals of the parameters marginalising over all other parameters: \n ")
                for i in range(self.npars):
                    f.write("{0}  = {1:.5f} +/- {2:.5f} \n".format(
                        self.parLabels[i],self.mean[i],np.sqrt(np.diag(self.cov)[i]))
                    )

                f.write("\n")
                f.write("MAP Solution: \n")
                for i in range(self.npars):
                    f.write("{0}  = {1:.5f} \n".format(self.parLabels[i],self.bestPars[i])) #
                    f.write("with Posterior probability ln(P*) = {0:.5f}\n".format(self.results.logwt[-1]))
                    f.write("and likelihood ln(L*) = {0:.5f}\n".format(self.results.logl[-1]))
                
        
class DynestyDynamicSearch(DynestySearch): #I think this can inheret almost everything except __init__ from the Static sampler
    """
    A class to use Dynesty to explore parameters space with dynamic nested sampling
    """

    def __init__(self, 
                 ):
        pass
    



                 
