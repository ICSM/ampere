""" Mixins for inference

This module contains a set of mixins used by Ampere's inference classes to 
enable smoother composition of samplers. The idea is to maximise the sharing of 
code without ending up with too many layers of single inheretance and avoiding 
problems with multiple inheretance. 

Classes defined here *MUST* be defined as mixins only. That means their only 
parent class should be object (or another mixin whose methods they overload or 
extend) and they must *NOT* define an __init__() method, or any methods defined 
in BaseSearch or any of its children. All class names should have "Mixin" 
appended to them so it is obvious what they are being used for when defining 
other classes

Mixins may define only one method for maximal modularity, or may define several 
if there is no reason to separate them. If the method intended to be public 
(i.e. that users might interact with) depends on other methods, they should 
either be defined within the mixin as well, or inhereted by the mixin from other
mixins for clarity (so that all the methods required by one public method get
inhereted at the same time).
"""


from __future__ import print_function

import contextlib
import numpy as np
import logging
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy.special import logsumexp

from ..data import Photometry, Spectrum

class ScipyMinMixin(object):
    """
    A simple mixin for including Scipy's minimiser in a class

    The idea with this is to provide an interface for simple MAP inference. 
    """
    def preopt(self, start, **kwargs):
        """Approximate the maximum a posteriori (MAP) solution using scipy.optimize.minimize.


        Parameters
        ----------
        start : array_like
            The starting position for the optimization.
        kwargs : optional
            Additional keyword arguments to pass to `scipy.optimize.minimize`.
        
        Returns
        -------
        solution : array_like
            The approximate MAP solution.

        Generated with Chat-GPT
        """
        
        neglnprob = lambda *args: -self.lnprob(*args)
        from scipy.optimize import minimize
        self.logger.info("Using scipy.minimize to find approximate MAP solution")
        self.logger.info("starting from position: %s",start)
        solution = minimize(neglnprob, start)
        solution = [p for p in solution.x]
        self.logger.info("Minimization complete, final position: %s ",solution)
        return solution

class _PostProcessorMixin(object):
    """ A base mixin for post-processing results

    This is a base class for mixins providing methods for computing diagnostic
    statistics and creating plots to explore the output. This includes e.g. producing
    corner plots. Methods that are used by multiple mixins should be defined here.
    """

    def _plot_data(self, axes, **kwargs):
        for d in self.dataSet:
            if isinstance(d,(Photometry,Spectrum)):
                d.plot(ax = axes)

    def _plot_samples(self, axes,
                      n_post_samples=1000,
                      alpha=0.1,
                      sample_label=None,
                      sample_color=None,
                      **kwargs):
        if sample_label is None:
            sample_label = 'Samples'
        if sample_color is None:
            sample_color = 'k'
        for s in self.samples[np.random.randint(len(self.samples),
                                                size=n_post_samples)]:
            with contextlib.suppress(ValueError):
                self.model(*s[:self.nparsMod])
                axes.plot(self.model.wavelength,self.model.modelFlux, '-',
                        color=sample_color,
                        alpha=alpha,label=sample_label, zorder = 0)
            i = self.nparsMod
            for d in self.dataSet:
                if isinstance(d,(Photometry,Spectrum)):
                    if d._hasNoiseModel:
                        d.plotRealisation(s[i:i+d.npars], ax=axes, **kwargs)
                    i+= d.npars


    def plot_posteriorpredictive(self, n_post_samples = 1000,
                                 plotfile=None, save = True,
                                 logx = False, logy = False,
                                 xlims = None, ylims = None,
                                 alpha = 0.1,show=False,
                                 return_figure = False,
                                 close_figure = False,
                                 **kwargs):
        '''
        Function to plot the data and samples drawn from the posterior of the models
        and data.

        Generates a plot of the data and posterior predictive samples for the given
        model and data. The plot is saved to a file specified by plotfile, or
        self.name_posteriorpredictive.png if plotfile is not provided. The x-axis is the
        wavelength (in $\mu$m) and the y-axis is the flux density (in mJy). The x-axis
        can be set to a log scale by setting logx to True. The y-axis can be set to a
        log scale by setting logy to True. The alpha parameter controls the transparency
        of the posterior predictive samples. Additional keyword arguments can be passed
        in **kwargs to specify additional plot properties. The posterior predictive
        samples are drawn from self.samples and the number of samples plotted is
        controlled by n_post_samples. The best-fit model is plotted with solid lines and
        the posterior predictive samples are plotted with transparent lines.

        Parameters
        ----------
        n_post_samples : int, optional
            The number of posterior samples to plot, by default 1000
        plotfile : str, optional
            The filename for the plot, by default None
        logx : bool, optional
            Whether to plot the x-axis on a log scale, by default False
        logy : bool, optional
            Whether to plot the y-axis on a log scale, by default False
        alpha : float, optional
            The transparency of the samples, by default 0.1
        **kwargs : optional
            Additional keyword arguments to pass to matplotlib.pyplot.plot

        Returns
        -------
        None

        Generated with Chat-GPT
        '''

        if plotfile is None:
            plotfile = f"{self.name}_posteriorpredictive.png"

        fig,axes = plt.subplots(1,1,figsize=(8,6))
        axes.set_xlabel(r"Wavelength ($\mu$m)")
        axes.set_ylabel(r"Flux density (mJy)")

        if logx:
            axes.set_xscale('log')
        if logy: 
            axes.set_yscale('log')

        #observations
        for d in self.dataSet:
            if isinstance(d,(Photometry,Spectrum)):
                d.plot(axes = axes)
        # posterior predictive samples
        for s in self.samples[np.random.randint(len(self.samples), size=n_post_samples)]:
            with contextlib.suppress(ValueError):
                self.model(*s[:self.nparsMod])
                axes.plot(self.model.wavelength,self.model.modelFlux, '-', color='k',
                          alpha=alpha,label='Samples', zorder = 0)
            i = self.nparsMod
            for d in self.dataSet:
                if isinstance(d,(Photometry,Spectrum)):
                    if d._hasNoiseModel:
                        d.plotRealisation(s[i:i+d.npars], ax=axes)
                    i+= d.npars


        #best fit model
        try:
            self.model(*self.bestPars[:self.nparsMod])
            axes.plot(self.model.wavelength,self.model.modelFlux, '-', color='magenta',
                      alpha=1.0,label='MAP', zorder=8)
        except ValueError:
            print("Error in MAP solution \n Skipping MAP in plot")

        if xlims is None:
            xmin = np.min([np.min(d.wavelength) for d in self.dataSet])
            xmax = np.max([np.max(d.wavelength) for d in self.dataSet])
            xmin = 0.5*xmin if logx else 0.8*xmin  # a little extra buffer if we're plotting in log space
            xmax = 2*xmax if logx else 1.2*xmax  # a little extra buffer if we're plotting in log space
            xlims = [xmin,xmax]
            axes.set_xlim(xlims)
        elif xlims == -1:
            xlims = None
        else:
            axes.set_xlim(xlims)
        if ylims is None:
            ymin = np.min([np.min(d.value) for d in self.dataSet])
            ymax = np.max([np.max(d.value) for d in self.dataSet])
            ymin = 0.5*ymin if logy else 0.8*ymin  # a little extra buffer if we're plotting in log space
            ymax = 2*ymax if logy else 1.2*ymax  # a little extra buffer if we're plotting in log space
            if ymax/ymin > 1e5 and logy:
                ymin = 1e-5*ymax
            ylims = [ymin,ymax]
            axes.set_ylim(ylims)
        elif ylims == -1:
            ylims = None
        else:
            axes.set_ylim(ylims)


        #These plots end up with too many labels for the legend, so we clobber the label
        # information so that only one of each one is plotted
        handles, labels = plt.gca().get_legend_handles_labels()
        # labels will be the keys of the dict, handles will be values
        temp = {k:v for k,v in zip(labels, handles)}
        plt.legend(temp.values(), temp.keys(), loc='best')

        #plt.legend()
        plt.tight_layout()
        if save:
            fig.savefig(plotfile,dpi=200)
        if show:
            plt.show()
        else:
            if return_figure:
                return fig
            elif close_figure:
                plt.close(fig)
                plt.clf()

class SBIPostProcessor(_PostProcessorMixin):
    """
    A mixin for post-processing SBI results

    This mixin provides functions for computing diagnostic statistics and creating 
    plots to explore the output. This includes e.g. producing corner plots. This 
    assumes structure similar to that of emcee. The `get_` methods are intended to be 
    extended or overloaded by classes that inheret from this mixin if necessary.
    """

    def get_credible_interval(self, percentiles=None, chain=None, **kwargs):
        """
        Compute a credible interval given an iterable of percentiles by the user

        Parameters
        ----------
        percentiles: list, optional
            List of percentiles to compute. The default is [16, 50, 84].
        chain: array, optional
            Samples from the posterior distribution. If not provided, the 
            `samples` attribute of the object will be used.
        kwargs: dict, optional
            Additional keyword arguments to pass to the `np.percentile` function.

        Returns
        -------
        res: list
            List of credible intervals for each parameter.
        
        Generated by Chat-GPT
        """
        if percentiles is None:
            percentiles = [16, 50, 84]
        if chain is None:
            chain = self.samples
        self.res = []
        for i in range(self.npars):
            a = np.percentile(self.samples[:,i], percentiles)
            self.res.append([a[1], a[2]-a[1], a[1]-a[0]])
        return self.res

    def get_map(self, **kwargs):
        """Returns the maximum a posteriori parameters of the posterior distribution.

Parameters
**kwargs:
Additional keyword arguments passed to the scipy.optimize.minimize function.

Returns
bestPars : ndarray
The maximum a posteriori parameters.
        
        Generated with Chat-GPT
        """
        self.logger.info("Estimating MAP parameters")
        self.bestPars = self.posterior.map()
        return self.bestPars

    def plot_covmats(self):
        """
    Plot the covariance matrices for each spectrum in the data set.

    Parameters
    ----------
    self : object
        An instance of the `BayesModel` class.

    Returns
    -------
    None

        Generated with Chat-GPT
    """
        istart = self.nparsMod
        for i, d in enumerate(self.dataSet):
            if isinstance(d, Photometry):
                continue
            elif isinstance(d, Spectrum):
                fig, ax = plt.subplots(1,1)
                d.cov([self.res[istart+1][0],self.res[istart+2][0]])
                #ax0.set_title('Covariance matrix')
                im = ax.imshow(np.log10(d.covMat))
                istart+=d.npars
                fig.savefig(f"{self.name}_covMat_{str(i)}.png")

    def plot_posteriorpredictive(self,
                                 n_post_samples = 100,
                                 plotfile=None, save = True,
                                 logx = False, logy = False,
                                 xlims = None, ylims = None,
                                 alpha = 0.1,show=False,
                                 return_figure = False,
                                 close_figure = False,
                                 fig=None, axes=None,
                                 sample_label=None,
                                 sample_color=None,
                                 **kwargs):
        '''
        Function to plot the data and samples drawn from the posterior of the models
        and data.

        Generates a plot of the data and posterior predictive samples for the given
        model and data. The plot is saved to a file specified by plotfile, or
        self.name_posteriorpredictive.png if plotfile is not provided. The x-axis is the
        wavelength (in $\mu$m) and the y-axis is the flux density (in mJy). The x-axis
        can be set to a log scale by setting logx to True. The y-axis can be set to a log
        scale by setting logy to True. The alpha parameter controls the transparency of
        the posterior predictive samples. Additional keyword arguments can be passed in
        **kwargs to specify additional plot properties. The posterior predictive samples
        are drawn from self.samples and the number of samples plotted is controlled by
        n_post_samples. The best-fit model is plotted with solid lines and the posterior
        predictive samples are plotted with transparent lines.

        Parameters
        ----------
        n_post_samples : int, optional
            The number of posterior samples to plot, by default 1000
        plotfile : str, optional
            The filename for the plot, by default None
        logx : bool, optional
            Whether to plot the x-axis on a log scale, by default False
        logy : bool, optional
            Whether to plot the y-axis on a log scale, by default False
        alpha : float, optional
            The transparency of the samples, by default 0.1
        **kwargs : optional
            Additional keyword arguments to pass to matplotlib.pyplot.plot

        Returns
        -------
        None
        '''
        if plotfile is None:
            plotfile = f"{self.name}_posteriorpredictive.png"
        if fig is None and axes is None:
            fig,axes = plt.subplots(1,1,figsize=(8,6))
        elif fig is not None and axes is None:
            axes = fig.add_subplot(111)
        elif fig is None:
            fig = axes.get_figure()

        if sample_label is None:
            sample_label = 'Samples'
        if sample_color is None:
            sample_color = 'k'

        #observations
        self._plot_data(axes = axes, **kwargs)

        #For SBI, we will plot cached models if we can, because the model is probably too slow to generate new ones
        if self.cache_models:
            #First, we have to calculate the *posterior* logprob of each *prior* sample
            post_probs = self.posterior.log_prob(self.thetas).exp()
            print(self.thetas.size())
            post_probs = post_probs/ post_probs[post_probs.isfinite()].sum() #divide by log(sum(probabilites) to normalise them
            print(np.sum(post_probs[post_probs.isfinite()].numpy()))
            #next we need to sample from these, with the most-probable models being most likely to be selected
            rng = np.random.default_rng()
            chosen_thetas = rng.choice(self.thetas[(post_probs.numpy() > 0) & post_probs.isfinite().numpy()], size = n_post_samples, p = post_probs.numpy()[(post_probs.numpy() > 0) & post_probs.isfinite().numpy()], axis = 0) #, replace=False)
            for t in chosen_thetas:
                key = tuple(t[:self.nparsMod])
                r = self.cache[key]
                axes.plot(r.spectrum["wavelength"], r.spectrum["flux"], '-', color='k', alpha = alpha, label = 'Samples', zorder=0)
        else:
            #posterior predictive samples
            self._plot_samples(axes = axes, n_post_samples=n_post_samples, alpha = alpha,
                               sample_label=sample_label, sample_color=sample_color,
                               **kwargs)
            # for s in self.samples[np.random.randint(len(self.samples), size=n_post_samples)]:
            #     with contextlib.suppress(ValueError):
            #         self.model(*s[:self.nparsMod])
            #         axes.plot(self.model.wavelength,self.model.modelFlux, '-', color='k', alpha=alpha,label='Samples', zorder = 0)
            #     i = self.nparsMod
            #     for d in self.dataSet:
            #         if isinstance(d,(Photometry,Spectrum)):
            #             if d._hasNoiseModel:
            #                 d.plotRealisation(s[i:i+d.npars], ax=axes)
            #             i+= d.npars


        #We won't try plotting the MAP in this case, because it almost certainly hasn't been generated a priori
        #best fit model
        #try:
        #    self.model(*self.bestPars[:self.nparsMod])
        #    axes.plot(self.model.wavelength,self.model.modelFlux, '-', color='k', alpha=1.0,label='MAP', zorder=8)
        #except ValueError:
        #    print("Error in MAP solution \n Skipping MAP in plot")

        #Now just before saving, we're going to clean up the presentation of the figure a bit
        axes.set_xlabel(r"Wavelength ($\mu$m)")
        axes.set_ylabel(r"Flux density (mJy)")

        if logx:
            axes.set_xscale("log")
        if logy:
            axes.set_yscale("log")

        if xlims is None:
            xmin = np.min([np.min(d.wavelength) for d in self.dataSet])
            xmax = np.max([np.max(d.wavelength) for d in self.dataSet])
            xmin = 0.5*xmin if logx else 0.8*xmin  # a little extra buffer if we're plotting in log space
            xmax = 2*xmax if logx else 1.2*xmax  # a little extra buffer if we're plotting in log space
            xlims = [xmin,xmax]
            axes.set_xlim(xlims)
        elif xlims == -1:
            xlims = None
        else:
            axes.set_xlim(xlims)
        if ylims is None:
            ymin = np.min([np.min(d.value) for d in self.dataSet])
            ymax = np.max([np.max(d.value) for d in self.dataSet])
            ymin = 0.5*ymin if logy else 0.8*ymin  # a little extra buffer if we're plotting in log space
            ymax = 2*ymax if logy else 1.2*ymax  # a little extra buffer if we're plotting in log space
            ylims = [ymin,ymax]
            axes.set_ylim(ylims)
        elif ylims == -1:
            ylims = None
        else:
            axes.set_ylim(ylims)



        #These plots end up with too many labels for the legend, so we clobber the label information so that only one of each one is plotted
        handles, labels = plt.gca().get_legend_handles_labels()
        # labels will be the keys of the dict, handles will be values
        temp = {k:v for k,v in zip(labels, handles)}
        plt.legend(temp.values(), temp.keys(), loc='best')

        #plt.legend()
        plt.tight_layout()
        if save:
            fig.savefig(plotfile,dpi=200)
        if show:
            plt.show()
        elif return_figure:
            return fig
        elif close_figure:
            plt.close(fig)
            plt.clf()

    def plot_corner(self, **kwargs):
        """
        Generate a corner plot
        """
        if self.npars > 10:
            self.plot_multiple_corner(**kwargs)
        else:
            import corner
            try:
                fig2 = corner.corner(self.samples,labels=self.parLabels, **kwargs)
            except (ImportError, ValueError):
                fig2 = corner.corner(self.samples.numpy(),labels=self.parLabels, **kwargs)
            fig2.savefig(self.name+"_"+"corner.png")

    def plot_multiple_corner(self, **kwargs):
        '''
        Generate separate corner plots for batches of parameters.
        
        Parameters
        ----------
        **kwargs : optional
            Additional keyword arguments to pass to corner.corner.
        '''
        import corner
        try:
            # first we make the corner plot for the model parameters
            fig2 = corner.corner(self.samples[:,:self.nparsMod],labels=self.parLabels[:self.nparsMod], **kwargs)
            fig2.savefig(self.name+"_"+"corner_model.png")
            # then we make the corner plot for the data parameters
            fig3 = corner.corner(self.samples[:,self.nparsMod:],labels=self.parLabels[self.nparsMod:], **kwargs)
            fig3.savefig(self.name+"_"+"corner_data.png")
        except (ImportError, ValueError):
            # first we make the corner plot for the model parameters
            fig2 = corner.corner(self.samples.numpy()[:,:self.nparsMod],labels=self.parLabels[:self.nparsMod], **kwargs)
            fig2.savefig(self.name+"_"+"corner_model.png")
            # then we make the corner plot for the data parameters
            fig3 = corner.corner(self.samples.numpy()[:,self.nparsMod:],labels=self.parLabels[self.nparsMod:], **kwargs)
            fig3.savefig(self.name+"_"+"corner_data.png")

    def summary(self, interval=None, **kwargs):
        """
        Generate the required summary statistics
        """
        
        if interval is None:
            interval = [16, 50, 84]
        self.res = self.get_credible_interval(interval)
        self.bestPars = self.get_map()

    def print_summary(self, outfile=None, **kwargs):
        """
        Generate a summary of the run

        This is currently overly complicated. In future, this functionality should end up being handled and written
        by a custom logger - summary() should handle this logging
        """
        #First generate the summary statistics
        self.summary()
        file_output = {}
        #Now start printing things
        #with open(outfile, 'w') as f:
        if self.bestPars is not None:
            self.logger.info("MAP solution:")
            #print("MAP Solution: ")
            for i in range(self.npars):
                self.logger.info("%s  = %.5f", self.parLabels[i],self.bestPars[i]) #


        if self.res is not None:
            self.logger.info("Median and confidence intervals for parameters in order:")
            for i in range(self.npars):
                self.logger.info("%s  = %.5f + %.5f - %.5f", self.parLabels[i],self.res[i][0],self.res[i][1],self.res[i][2])




class SimpleMCMCPostProcessor(_PostProcessorMixin):
    """
    A mixin for post-processing MCMC results

    This mixin provides functions for computing diagnostic statistics and creating 
    plots to explore the output. This includes e.g. the autocorrelation time of the
    chains, or producing corner plots. This assumes structure similar to that of 
    emcee. The `get_` methods are intended to be extended or overloaded by classes 
    that inheret from this mixin if necessary.

    
    """

    
    
    def plot_trace(self, **kwargs):
        """
        Plot the trace of each parameter in the MCMC sampling.
    
        Parameters
        ----------
        **kwargs : optional
            Additional keyword arguments to be passed to `matplotlib.pyplot.plot`
        
        Generated with Chat-GPT
        """
        try:
            tauto = self.tauto
        except AttributeError:
            tauto = self.get_autocorr()
        fig, axes = plt.subplots(self.npars, 1, sharex=True, figsize=(8, 9))
        for i in range(self.npars):
            axes[i].plot(self.sampler.chain[:, :, i].T, color="k", alpha=0.4) #,label=self.parLabels[i])
            axes[i].set_xlim(0, self.nsamp)
            ylims = axes[i].get_ylim()
            xlims = axes[i].get_xlim()
            axes[i].plot([10*tauto[i],10*tauto[i]],[ylims[0],ylims[1]],color="red",marker="",linestyle="-")#,label=r"10$\times t_{\rm autocorr}$")
            axes[i].add_patch(Polygon([[xlims[0], ylims[0]], [xlims[0], ylims[1]], [self.burnin, ylims[1]], [self.burnin, ylims[0]]], closed=True,
                      fill=True, color='darkgrey'))
            axes[i].set_xlabel("Steps")
            axes[i].set_ylabel(self.parLabels[i])
            axes[i].label_outer()
            axes[i].set_xlim(0, self.nsamp)

        plt.tight_layout()
        fig.savefig(self.name+"_walkers.png")

    # Geweke convergence test threshold
    geweke_max = 2.0
    
    def get_geweke(self, chain = None, **kwargs):
        """
        Measure the Geweke statistic

        The Geweke statistics is a useful criterion for convergence, based on 
        the drift in the mean of the chains. It takes the first 25% of the 
        burned-in chain, and the last 25%, computes the mean of each parameter 
        for those two sections, and computes the difference between the means 
        divided by the sqrt of the sum of the variances. If this is greater 
        than a threshold (default: 2.0) then the drift is significant and the 
        chain has not converged.

        It is important to note that just because G < 2, it doesn't mean that 
        the chain *HAS* converged.

        Parameters
        ----------
        chain : array-like, optional
            The chain to compute the statistic for. If not passed, we will
            use the full set of samples
        """
        self.converged = True
        if chain is None:
            chain = self.allSamples 
        a = chain.reshape((-1, self.npars))[:self.nsamp//4, :] #sampler.get_chain(flat=True)[:num_steps // 4, i]
        b = chain.reshape((-1, self.npars))[-self.nsamp//4:, :] #sampler.get_chain(flat=True)[-num_steps // 4:, i]
        geweke_z = (a.mean(axis=0) - b.mean(axis=0)) / (np.var(a, axis=0) + np.var(b, axis=0))**0.5
        for i, p in enumerate(self.parLabels):
            if geweke_z[i] > self.geweke_max:
                self.logger.warning("geweke drift (z=%.3f) detected for parameter '%s'", geweke_z[i], p)
                self.converged = False
        return geweke_z

    def get_autocorr(self, chain = None, quiet=True, **kwargs):
        """
        Compute the autocorrelation time of a chain

        This relies on importing the integrated_time function from emcee. Other 
        samplers which define their own methods should be used to overload this.

        Parameters
        ----------
        chain : array-like, optional
            The chain to compute the statistic for. If not passed, we will
            use the full set of samples
        """
        try:
            self.tauto = self.sampler.get_autocorr_time(quiet=True)
        except AttributeError:
            from emcee.autocorr import integrated_time
            if chain is None:
                chain = self.allSamples
            self.tauto = integrated_time(chain, quiet=True)
        return self.tauto

    def get_ess(self, chain = None, **kwargs):
        """
        Placeholder for Effective Sample Size computation
        """
        self.ess = self.nwalkers * self.nsamp / np.mean(self.tauto)
        return self.ess

    def get_rhat(self, chain = None, n_splits = 2, threshold = 0.05, discard = None, **kwargs):
        """
        Placeholder for the Split Rhat statistic

        This uses the post-burn in flattened chain, and splits it into two 
        chunks to simulate having multiple chains (which is what is really 
        required for rhat).

        Parameters
        ----------
        chain : array-like, optional
            The chain to compute the statistic for. If not passed, we will
            use the full set of samples
        n_splits : float, optional, default 2
            The number of chunks to split the chain in to
        threshold : float, default 0.05
            The convergence threshold for the Rhat statistic
            If Rhat > 1 + threshold, a warning will be issued
        """
        #First get the chain
        if chain is None:
            chain = self.samples
        
        #Then split the samples into n_splits chunks
        ndim = chain.shape[-1]
        chunks = np.array_split(chain, n_splits)
        #print(chunks)

        means = []
        var = []
        length = []
        
        #Then compute the mean, variance and length of each chunk
        for chunk in chunks:
            means.append(np.mean(chunk.reshape(-1, ndim), axis = 0))
            var.append(np.var(chunk.reshape(-1, ndim), axis = 0))
            length= len(chunk.reshape(-1, ndim))


        
        #Then we get on with the RHat bit

        #First we need the variance between chains:
        B = length * np.var(means, ddof=1, axis = 0)
        

        #And the variance within chains
        W = np.mean(var)
        

        #Now computed a weighted variances
        weighted_var = (1 - 1/length) * W + B/length
        

        self.rhat = np.mean(np.sqrt(weighted_var / W))
        if self.rhat - 1 > threshold:
            self.logger.warning("Rhat = %.3f > %.3f. Chain is not converged! \n You should run a longer chain. ", self.rhat, 1+threshold)
            self.converged = False
        else:
            self.logger.info("Rhat = %.3f <= %.3f ", self.rhat, 1+threshold)
            self.converged = True
        return self.rhat

    def get_credible_interval(self, percentiles = [16, 50, 84], chain = None, **kwargs):
        """
        Compute a credible interval given an iterable of percentiles by the user

Parameters
----------
percentiles: list, optional
    List of percentiles to compute. The default is [16, 50, 84].
chain: array, optional
    Samples from the posterior distribution. If not provided, the `samples` attribute of the object will be used.
kwargs: dict, optional
    Additional keyword arguments to pass to the `np.percentile` function.

Returns
-------
res: list
    List of credible intervals for each parameter.
        
        Generated by Chat-GPT
        """
        if chain is None:
            chain = self.samples
        self.res = []
        for i in range(self.npars):
            a = np.percentile(self.samples[:,i], percentiles)
            self.res.append([a[1], a[2]-a[1], a[1]-a[0]])
        return self.res

    def get_map(self, **kwargs):
        """Returns the maximum a posteriori parameters of the posterior distribution.

Parameters
**kwargs:
Additional keyword arguments passed to the scipy.optimize.minimize function.

Returns
bestPars : ndarray
The maximum a posteriori parameters.
        
        Generated with Chat-GPT
        """
        row_ind, col_ind = np.unravel_index(np.argmax(self.sampler.lnprobability), self.sampler.lnprobability.shape)
        self.bestPars = self.sampler.chain[row_ind, col_ind, :]
        return self.bestPars

    def summary(self, interval = [16, 50, 84], **kwargs):
        """
        Generate the required summary statistics
        """
        self.tauto = self.get_autocorr()
        self.geweke = self.get_geweke()
        self.ess = self.get_ess()
        self.rhat = self.get_rhat()
        self.res = self.get_credible_interval(interval)
        self.bestPars = self.get_map()

    def print_summary(self, outfile=None, **kwargs):
        """
        Generate a summary of the run

        This is currently overly complicated. In future, this functionality should end up being handled and written
        by a custom logger - summary() should handle this logging
        """
        #First generate the summary statistics
        self.summary()
        file_output = {}
        #Now start printing things
        #with open(outfile, 'w') as f:
        if self.bestPars is not None:
            self.logger.info("MAP solution:")
            #print("MAP Solution: ")
            for i in range(self.npars):
                self.logger.info("%s  = %.5f", self.parLabels[i],self.bestPars[i]) #
                #file_ouput['MAP']+="{0}  = {1:.5f}\n".format(self.parLabels[i],self.bestPars[i])

        if self.res is not None:
            self.logger.info("Median and confidence intervals for parameters in order:")
            for i in range(self.npars):
                self.logger.info("%s  = %.5f + %.5f - %.5f", self.parLabels[i],self.res[i][0],self.res[i][1],self.res[i][2])
                #file_output['CI']+="{0}  = {1[0]:.5f} + {1[1]:.5f} - {1[2]:.5f}\n".format(self.parLabels[i],self.res[i])
        if self.tauto is not None:
            self.logger.info("Mean Autocorrelation Time: %.5f",np.mean(self.tauto))
            self.logger.info("Autocorrelation times for each parameter:")
            #tauto = self.sampler.get_autocorr_time()
            for i in range(self.npars):
                self.logger.info("%s  = %.0f steps",self.parLabels[i],self.tauto[i])

        if self.ess is not None:
            self.logger.info("Effective sample size: %.3f", self.ess)

    def plot_corner(self, **kwargs):
        """
        Generate the triangle plot for the chains using the `corner` package

        Parameters
        ---------
        **kwargs:
        Additional keyword arguments to be passed to corner.corner().

        Generated with Chat-GPT
        """
        import corner
        fig2 = corner.corner(self.samples,labels=self.parLabels, **kwargs)
        fig2.savefig(self.name+"_corner.png")

    def plot_covmats(self):
        """
    Plot the covariance matrices for each spectrum in the data set.

    Parameters
    ----------
    self : object
        An instance of the `BayesModel` class.

    Returns
    -------
    None

        Generated with Chat-GPT
    """
        istart = self.nparsMod
        for i, d in enumerate(self.dataSet):
            if isinstance(d, Photometry):
                continue
            elif isinstance(d, Spectrum):
                fig, ax = plt.subplots(1,1)
                d.cov([self.res[istart+1][0],self.res[istart+2][0]])
                #ax0.set_title('Covariance matrix')
                im = ax.imshow(np.log10(d.covMat))
                istart+=d.npars
                fig.savefig(self.name+"_covMat_"+str(i)+".png")

    def plot_posteriorpredictive(self, n_post_samples = 1000,
                                 plotfile=None, save = True,
                                 logx = False, logy = False,
                                 xlims = None, ylims = None,
                                 alpha = 0.1,show=False,
                                 return_figure = False,
                                 close_figure = False,
                                 fig = None, axes=None,
                                 sample_label=None,
                                 sample_color=None,
                                 **kwargs):
        '''
        Function to plot the data and samples drawn from the posterior of the models
        and data.

        Generates a plot of the data and posterior predictive samples for the given
        model and data. The plot is saved to a file specified by plotfile, or
        self.name_posteriorpredictive.png if plotfile is not provided. The x-axis is the
        wavelength (in $\mu$m) and the y-axis is the flux density (in mJy). The x-axis
        can be set to a log scale by setting logx to True. The y-axis can be set to a
        log scale by setting logy to True. The alpha parameter controls the transparency
        of the posterior predictive samples. Additional keyword arguments can be passed
        in **kwargs to specify additional plot properties. The posterior predictive
        samples are drawn from self.samples and the number of samples plotted is
        controlled by n_post_samples. The best-fit model is plotted with solid lines and
        the posterior predictive samples are plotted with transparent lines.

        Parameters
        ----------
        n_post_samples : int, optional
            The number of posterior samples to plot, by default 1000
        plotfile : str, optional
            The filename for the plot, by default None
        logx : bool, optional
            Whether to plot the x-axis on a log scale, by default False
        logy : bool, optional
            Whether to plot the y-axis on a log scale, by default False
        alpha : float, optional
            The transparency of the samples, by default 0.1
        **kwargs : optional
            Additional keyword arguments to pass to matplotlib.pyplot.plot

        Returns
        -------
        None

        Generated with Chat-GPT
        '''

        if plotfile is None:
            plotfile = f"{self.name}_posteriorpredictive.png"

        if fig is None and axes is None:
            fig,axes = plt.subplots(1,1,figsize=(8,6))
        elif fig is not None and axes is None:
            axes = fig.add_subplot(111)
        elif fig is None:
            fig = axes.get_figure()

        axes.set_xlabel(r"Wavelength ($\mu$m)")
        axes.set_ylabel(r"Flux density (mJy)")

        if logx:
            axes.set_xscale('log')
        if logy: 
            axes.set_yscale('log')

        #observations
        self._plot_data(axes = axes, **kwargs)
        # for d in self.dataSet:
        #     if isinstance(d,(Photometry,Spectrum)):
        #         d.plot(ax = axes)
        # posterior predictive samples
        self._plot_samples(axes = axes, n_post_samples=n_post_samples, alpha=alpha,
                           **kwargs)
        # for s in self.samples[np.random.randint(len(self.samples), size=n_post_samples)]:
        #     with contextlib.suppress(ValueError):
        #         self.model(*s[:self.nparsMod])
        #         axes.plot(self.model.wavelength,self.model.modelFlux, '-', color='k',
        #                   alpha=alpha,label='Samples', zorder = 0)
        #     i = self.nparsMod
        #     for d in self.dataSet:
        #         if isinstance(d,(Photometry,Spectrum)):
        #             if d._hasNoiseModel:
        #                 d.plotRealisation(s[i:i+d.npars], ax=axes)
        #             i+= d.npars


        #best fit model
        try:
            self.model(*self.bestPars[:self.nparsMod])
            axes.plot(self.model.wavelength,self.model.modelFlux, '-', color='magenta',
                      alpha=1.0,label='MAP', zorder=8)
        except ValueError:
            print("Error in MAP solution \n Skipping MAP in plot")

        if xlims is None:
            xmin = np.min([np.min(d.wavelength) for d in self.dataSet])
            xmax = np.max([np.max(d.wavelength) for d in self.dataSet])
            xmin = 0.5*xmin if logx else 0.8*xmin  # a little extra buffer if we're plotting in log space
            xmax = 2*xmax if logx else 1.2*xmax  # a little extra buffer if we're plotting in log space
            xlims = [xmin,xmax]
            axes.set_xlim(xlims)
        elif xlims == -1:
            xlims = None
        else:
            axes.set_xlim(xlims)
        if ylims is None:
            ymin = np.min([np.min(d.value) for d in self.dataSet])
            ymax = np.max([np.max(d.value) for d in self.dataSet])
            ymin = 0.5*ymin if logy else 0.8*ymin  # a little extra buffer if we're plotting in log space
            ymax = 2*ymax if logy else 1.2*ymax  # a little extra buffer if we're plotting in log space
            if ymax/ymin > 1e5 and logy:
                ymin = 1e-5*ymax
            ylims = [ymin,ymax]
            axes.set_ylim(ylims)
        elif ylims == -1:
            ylims = None
        else:
            axes.set_ylim(ylims)


        #These plots end up with too many labels for the legend, so we clobber the label
        # information so that only one of each one is plotted
        handles, labels = plt.gca().get_legend_handles_labels()
        # labels will be the keys of the dict, handles will be values
        temp = {k:v for k,v in zip(labels, handles)}
        plt.legend(temp.values(), temp.keys(), loc='best')

        #plt.legend()
        plt.tight_layout()
        if save:
            fig.savefig(plotfile,dpi=200)
        if show:
            plt.show()
        elif return_figure:
            return fig
        elif close_figure:
            plt.close(fig)
            plt.clf()

        
class ArvizPostProcessor(_PostProcessorMixin):
    """
    A mixin for post-processing MCMC results that uses ArViz instead of rolling 
    its own
    """
    pass

class SimulatorMixin(object):
    """
    A mixin for Simulation-Based Inference (SBI) such as ABC
    """

    def simulate(self, theta, **kwargs):
        # First call the model


        #Then iterate over the dataset to produce synthetic observations (including noise)
        pass

class AnotherMixin(object):
    pass
