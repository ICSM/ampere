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

import numpy as np

class ScipyMinMixin(object):
    """
    A simple mixin for including Scipy's minimiser in a class

    The idea with this is to provide an interface for simple MAP inference. 
    """
    def minimize(self, start, **kwargs):
        neglnprob = lambda *args: -self.lnprob(*args)
        from scipy.optimize import minimize
        print("Using scipy.minimize to find approximate MAP solution")
        print("starting from position: ",start)
        solution = minimize(neglnprob, start)
        solution = [p for p in solution.x]
        print("Minimization complete, final position: ",solution)
        return solution


class SimpleMCMCPostProcessor(object):
    """
    A mixin for post-processing MCMC results

    This mixin provides functions for computing diagnostic statistics and creating 
    plots to explore the output. This includes e.g. the autocorrelation time of the
    chains, or producing corner plots. This assumes structure similar to that of 
    emcee. The get_ methods are intended to be extended or overloaded by classes 
    that inheret from this mixin if necessary.

    
    """
    def plot_trace(self, **kwargs):
        tauto = self.autocorr()
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
        fig.savefig("walkers.png")

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
        """
        self.converged = True
        if chain is None:
            chain = self.allSamples 
        a = chain.reshape((-1, self.npars))[:self.nsamp//4, :] #sampler.get_chain(flat=True)[:num_steps // 4, i]
        b = chain.reshape((-1, self.npars))[-self.nsamp//4:, :] #sampler.get_chain(flat=True)[-num_steps // 4:, i]
        geweke_z = (a.mean(axis=0) - b.mean(axis=0)) / (np.var(a, axis=0) + np.var(b, axis=0))**0.5
        for i, p in enumerate(self.parLabels):
            if geweke_z[i] > self.geweke_max:
                print("geweke drift (z={0:.2f}) detected for parameter '{1}'" % (geweke_z, p))
                self.converged = False
        return geweke_z

    def get_autocorr(self, chain = None, quiet=True, **kwargs):
        """
        Compute the autocorrelation time of a chain

        This relies on importing the integrated_time function from emcee. Other 
        samplers which define their own methods should be used to overload this.
        """
        from emcee.autocorr import integrated_time
        if chain is None:
            chain = self.allSamples
        tauto = integrated_time(chain, quiet=True)
        return tauto

    def get_ess(self, chain = None, **kwargs):
        """
        Placeholder for Effective Sample Size computation
        """
        pass

    def get_rhat(self, chain = None, **kwargs):
        """
        Placeholder for the Split Rhat statistic

        This uses the post-burn in flattened chain, and splits it into two 
        chunks to simulate having multiple chains (which is what is really 
        required for rhat).
        """
        pass

    def get_credible_interval(self, percentiles, chain = None, **kwargs):
        """
        Compute a credible interval given an iterable of percentiless by the user
        """
        if chain is None:
            chain = self.samples
        res = []
        for i in range(self.npars):
            a = np.percentile(self.samples[:,i], [16, 50, 84])
            res.append([a[1], a[2]-a[1], a[1]-a[0]])
        return res

    def get_map(self, **kwargs):
        pass

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
        if self.bestPars:
            file_output['MAP'] = ""
            print("MAP Solution: ")
            for i in range(self.npars):
                print("{0}  = {1:.5f}".format(self.parLabels[i],self.bestPars[i])) #
                file_ouput['MAP']+="{0}  = {1:.5f}\n".format(self.parLabels[i],self.bestPars[i])

        if self.res:
            file_output['CI'] = ""
            for i in range(self.npars):
                print("{0}  = {1[0]:.5f} + {1[1]:.5f} - {1[2]:.5f}".format(self.parLabels[i],self.res[i]))
                file_output['CI']+="{0}  = {1[0]:.5f} + {1[1]:.5f} - {1[2]:.5f}\n".format(self.parLabels[i],self.res[i])
        if self.tauto:
            file_output['tau'] = ""
            print("Mean Autocorrelation Time: {0:.5f}",np.mean(self.tauto))
            file_output['tau']+="Mean Autocorrelation Time: {0:.5f}\n",np.mean(self.tauto)
            print("Autocorrelation times for each parameter:")
            file_output['tau']+="Autocorrelation times for each parameter:\m"
            #tauto = self.sampler.get_autocorr_time()
            for i in range(self.npars):
                print("{0}  = {1:.0f} steps".format(self.parLabels[i],self.tauto[i]))
                file_output['tau']+="{0}  = {1:.0f} steps\n".format(self.parLabels[i],self.tauto[i])


    def plot_corner(self, **kwargs):
        import corner
        fig2 = corner.corner(self.samples,labels=self.parLabels, **kwargs)
        fig2.savefig("corner.png")

    def plot_covmats():
        pass

    def plot_posteriorpredictive(**kwargs):
        pass

        
class ArvizPostProcessor(object):
    """
    A mixin for post-processing MCMC results that uses ArViz instead of rolling 
    its own
    """
    pass

class AnotherMixin(object):
    pass
