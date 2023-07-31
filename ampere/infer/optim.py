from __future__ import print_function

import numpy as np
from datetime import datetime
from ..logger import Logger
from .basesearch import BaseSearch
from inspect import signature
import matplotlib.pyplot as plt
from ..data import Photometry, Spectrum


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

    _inference_method = "Optimisation"

    def __init__(self,
                 model=None,
                 data=None,
                 verbose=False,
                 parameter_labels=None,
                 name='',
                 namestyle="full",
                 **kwargs):
        self.model = model
        self.dataSet = data
        self.namestyle = namestyle
        self.name = name
        self.verbose = verbose

        self.setup_logging(verbose=verbose)
        self.logger.info("Welcome to ampere")
        self.logger.info("Setting up your inference problem:")
        self.logger.info("You are using %s", self._inference_method)
        self.logger.info("You have %s items in your dataset", str(len(data)))

        try:
            self.nparsMod = self.model.npars
        except AttributeError:
            sig = signature(model.__call__)
            # Always subtract **kwargs from the parameters, but don't need to
            # worry about self once it is bound to an instance
            self.nparsMod = len(sig.parameters) - 1
        # number of parameters to be passed into each set of data
        self.nparsData = [data.npars for data in self.dataSet]
        self.npars = int(self.nparsMod + np.sum(self.nparsData))

        if parameter_labels is None:
            # The user hasn't specified parameter labels, let's see if the
            # models and data have instead
            try:  # First the model parameters
                self.parLabels = self.model.parLabels
            except AttributeError:
                self.parLabels = [f'x{str(i)}' for i in range(self.nparsMod)]
            i = self.nparsMod
            for data in self.dataSet:
                try:
                    self.parLabels.extend(data.parLabels)
                except AttributeError:
                    self.parLabels.extend([f'x{str(i)}'
                                           for i in range(i, i+data.npars)])
                finally:
                    i += data.npars
        else:  # User isn't really supposed to use this interface, however...
            self.parLabels = parameter_labels

        self.verbose = verbose
        # if self.verbose:
        self.logger.info("This model has %d parameters.", self.nparsMod)
        self.logger.info("There are also %d parameters for the noise model",
                         self.npars - self.nparsMod)
        self.logger.info("Hence, there are a total of %d parameters to sample",
                         self.npars)
        self.logger.info("The parameter names are:")
        for label in self.parLabels:
            self.logger.info("%s", label)

        # Now we should check whether the number of parameters matches with
        # the parameter labels!
        if len(self.parLabels) != self.npars:
            self.logger.critical("You have %d parameters but %d labels",
                                 self.npars, len(self.parLabels))
            self.logger.critical("Please check parameters and labels")
            raise ValueError("Mismatch between number of parameters and labels"
                             )

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name=''):
        try:
            modelname = self.model.name
        except AttributeError:
            modelname = self.model.__class__.__name__
        if self.namestyle == 'full':
            self._name = ("ampere_"+str(datetime.now()
                                        ).replace(' ',
                                                  '_').replace(":", "-")[:-7]
                          + "_" + modelname + "_" + name
                          )
        elif self.namestyle == "model":
            self._name = f"ampere_{modelname}_{name}"
        elif self.namestyle == "short":
            self._name = name
        elif self.namestyle == "stamp":
            self._name = ("ampere_" +
                          str(datetime.now()).replace(' ',
                                                      '_').replace(":",
                                                                   "-")[:-7]
                          + "_" + name
                          )

    def optimize_with_restarts(self,
                               n_restarts=10,
                               guesses=None,
                               **kwargs):
        for i in range(n_restarts):
            self.logger.info("Optimisation restart %s", i)
            guess = guesses[i] if guesses is not None else None
            solution = self.optimize(guess=guess, return_solution=True,
                                     **kwargs)
            # now we have to check if this is the best solution so far
            # and add the final solution to the chain
            if i == 0:
                self.bestPars = self.solution
                self.bestLnprob = solution.fun
                self.chain = [(self.bestPars, self.bestLnprob)]
                i_best = 0
            else:
                if solution.fun > self.bestLnprob:
                    self.bestPars = self.solution
                    self.bestLnprob = solution.fun
                    i_best = i
                self.chain.append((self.solution, solution.fun))
            self.logger.info("Best solution so far: %s", self.bestPars)
            self.logger.info("Found at restart %s", i_best)
            self.logger.info("With lnprob: %s", self.bestLnprob)
        # now we have to set the solution to the best one
        self.solution = self.bestPars  # so that get_map works
        self.logger.info("Optimisation complete")
        # convert the chain to a numpy array so it is easier to work with
        self.chain = np.array(self.chain)

    def get_map(self, **kwargs):
        """Get the maximum a posteriori (MAP) point.
        """
        self.bestPars = self.solution
        return self.bestPars

    def summary(self, **kwargs):
        self.bestPars = self.get_map(**kwargs)

    def print_summary(self, **kwargs):
        self.summary(**kwargs)
        if self.bestPars is not None:
            self.logger.info("MAP solution:")
            for i in range(self.npars):
                self.logger.info("%s  = %.5f", self.parLabels[i],
                                 self.bestPars[i])

    def plot_posteriorpredictive(self,
                                 plotfile=None,
                                 save=True,
                                 logx=False,
                                 logy=False,
                                 xlims=None,
                                 ylims=None,
                                 show=False,
                                 return_figure=False,
                                 close_figure=False,
                                 **kwargs):
        '''
        Function to plot the data and samples drawn from the posterior of the
        models and data.

        Generates a plot of the data and the model for the MAP parameter
        values. The plot is saved to a file specified by plotfile, or
        self.name_posteriorpredictive.png if plotfile is not provided.
        The x-axis is the wavelength (in $\mu$m) and the y-axis is the flux
        density (in mJy). The x-axis can be set to a log scale by setting logx
        to True. The y-axis can be set to a log scale by setting logy to True.
        Additional keyword arguments can be passed in **kwargs to specify
        additional plot properties.

        Parameters
        ----------
        plotfile : str, optional
            The filename for the plot, by default None
        logx : bool, optional
            Whether to plot the x-axis on a log scale, by default False
        logy : bool, optional
            Whether to plot the y-axis on a log scale, by default False
        **kwargs : optional
            Additional keyword arguments to pass to matplotlib.pyplot.plot

        Returns
        -------
        fig : matplotlib.pyplot.figure
            if return_figure is True, the figure object is returned
        '''

        if plotfile is None:
            plotfile = f"{self.name}_posteriorpredictive.png"

        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        axes.set_xlabel(r"Wavelength ($\mu$m)")
        axes.set_ylabel(r"Flux density (mJy)")

        if logx:
            axes.set_xscale('log')
        if logy:
            axes.set_yscale('log')

        # observations
        for d in self.dataSet:
            if isinstance(d, (Photometry, Spectrum)):
                d.plot(ax=axes)

        # best fit model
        try:
            self.model(*self.bestPars[:self.nparsMod])
            axes.plot(self.model.wavelength, self.model.modelFlux, '-',
                      color='k', alpha=1.0, label='MAP', zorder=8)
        except ValueError:
            self.logger.warning("Error in MAP solution")
            self.logger.warning("Skipping MAP in plot")

        # These plots end up with too many labels for the legend, so we
        # clobber the label information so that only one of each one is plotted
        handles, labels = plt.gca().get_legend_handles_labels()
        # labels will be the keys of the dict, handles will be values
        temp = dict(zip(labels, handles))
        plt.legend(temp.values(), temp.keys(), loc='best')

        if xlims is not None:
            axes.set_xlim(xlims)
        if ylims is not None:
            axes.set_ylim(ylims)

        plt.tight_layout()
        if save:
            fig.savefig(plotfile, dpi=200)
        if show:
            plt.show()
        elif return_figure:
            return fig
        elif close_figure:
            plt.close(fig)
            plt.clf()

    def postProcess(self, **kwargs):
        self.print_summary(**kwargs)
        self.plot_posteriorpredictive(**kwargs)


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

    _inference_method = "Scipy Minimisation"

    def __init__(self,
                 model=None,
                 data=None,
                 verbose=False,
                 parameter_labels=None,
                 name='',
                 namestyle="full",
                 **kwargs):

        super().__init__(model=model, data=data, verbose=verbose,
                         parameter_labels=parameter_labels, name=name,
                         namestyle=namestyle, **kwargs)

    def optimize(self,
                 guess=None,
                 method=None,
                 tol=None,
                 return_solution=False,
                 **kwargs):

        from scipy.optimize import minimize
        self.logger.info("Starting optimisation")

        def neglnprob(*args):
            return -self.lnprob(*args)
        if guess is None:
            # no guess, let's attempt to draw one from the prior and use
            # a value roughly in the middle
            self.logger.info("No inital guess provided")
            try:
                self.logger.info("Attempting to draw from prior transform")
                u = np.array([0.5 for _ in range(self.npars)])
                guess = self.prior_transform(u)
            except AttributeError:
                self.logger.warning("No prior transform, using random guess")
                guess = np.random.rand(self.npars)
        self.logger.info("starting from position: %s", guess)
        solution = minimize(neglnprob, guess, method=method, tol=tol, **kwargs)
        self.solution = list(solution.x)
        self.logger.info("Optimisation completed in %s iterations",
                         solution.nit)
        self.logger.info("after %s function calls",
                         solution.nfev)
        if not solution.success:
            self.logger.warning("Optimisation did not converge")
            self.logger.warning("Message: %s", solution.message)
        if return_solution:
            return solution


class ScipyBasinOpt(OptimBase):
    """A wrapper for scipy.optimize.basinhopping

    _extended_summary_

    Parameters
    ----------
    OptimBase : _type_
        _description_
    """

    _inference_method = "Scipy Basin Hopping"

    def __init__(self,
                 model=None,
                 data=None,
                 verbose=False,
                 parameter_labels=None,
                 name='',
                 namestyle="full",
                 **kwargs):

        super().__init__(model=model, data=data, verbose=verbose,
                         parameter_labels=parameter_labels, name=name,
                         namestyle=namestyle, **kwargs)

    def optimize(self,
                 guess=None,
                 niter=100,
                 minimizer_kwargs=None,
                 return_solution=False,
                 **kwargs):
        from scipy.optimize import basinhopping
        self.logger.info("Starting optimisation")

        def neglnprob(*args):
            return -self.lnprob(*args)
        if guess is None:
            # no guess, let's attempt to draw one from the prior and use
            # a value roughly in the middle
            self.logger.info("No inital guess provided")
            try:
                self.logger.info("Attempting to draw from prior transform")
                u = np.array([0.5 for _ in range(self.npars)])
                guess = self.prior_transform(u)
            except AttributeError:
                self.logger.warning("No prior transform, using random guess")
                guess = np.random.rand(self.npars)
        self.logger.info("starting from position: %s", guess)
        solution = basinhopping(neglnprob, guess, niter=niter,
                                minimizer_kwargs=minimizer_kwargs, **kwargs)
        self.solution = list(solution.x)
        self.logger.info("Optimisation completed in %s iterations",
                         solution.nit)
        self.logger.info("after %s function calls",
                         solution.nfev)
        if not solution.success:
            self.logger.warning("Optimisation did not converge")
            self.logger.warning("Message: %s", solution.message)
        if return_solution:
            return solution


class ScipyDE(OptimBase):
    """A wrapper for scipy.optimize.differential_evolution

    _extended_summary_

    Parameters
    ----------
    OptimBase : _type_
        _description_
    """

    _inference_method = "Scipy Differential Evolution"

    def __init__(self,
                 model=None,
                 data=None,
                 verbose=False,
                 parameter_labels=None,
                 name='',
                 namestyle="full",
                 **kwargs):

        super().__init__(model=model, data=data, verbose=verbose,
                         parameter_labels=parameter_labels, name=name,
                         namestyle=namestyle, **kwargs)

    def optimize(self,
                 popsize=15,
                 maxiter=1000,
                 strategy='best1bin',
                 tol=0.01,
                 mutation=(0.5, 1),
                 recombination=0.7,
                 seed=None,
                 polish=True,
                 atol=0,
                 updating='immediate',
                 guess=None,
                 return_solution=False,
                 **kwargs):
        from scipy.optimize import differential_evolution
        self.logger.info("Starting optimisation")

        # DE requires a bounds argument to be passed
        # but our models may not be bounded
        # hence, we'll let DE run in the unit cube
        # and transform the parameters back to the
        # prior volume with the prior transform
        def neglnprob(*args):
            theta = self.prior_transform(*args)
            return -self.lnprob(theta)

        bounds = [(0, 1) for _ in range(self.npars)]
        if guess is None:
            # Unlike other scipy minimisers, DE requires a population
            # of guesses to start from, but since we've already
            # incorporated the prior transform, we'll just use
            # a population of random numbers in the unit cube
            self.logger.info("No inital guess provided")
            self.logger.info("Using random population")
            guess = np.random.rand(popsize, self.npars)
        elif (isinstance(guess, np.ndarray)
              and guess.shape != (popsize, self.npars)):
            self.logger.info("Guess is not a population, ")
            self.logger.info("adding random purturbations")
            guess = np.array([guess + tol*np.random.randn(self.npars)
                              for _ in range(popsize)])
        self.logger.info("starting from population: %s", guess)

        solution = differential_evolution(neglnprob, bounds,
                                          popsize=popsize, maxiter=maxiter,
                                          strategy=strategy, tol=tol,
                                          mutation=mutation,
                                          recombination=recombination,
                                          seed=seed, polish=polish, atol=atol,
                                          updating=updating,
                                          init=guess,
                                          **kwargs)
        # transform back to the prior volume
        self.solution = self.prior_transform(*list(solution.x))
        self.logger.info("Optimisation completed in %s iterations",
                         solution.nit)
        self.logger.info("after %s function calls",
                         solution.nfev)
        if not solution.success:
            self.logger.warning("Optimisation did not converge")
            self.logger.warning("Message: %s", solution.message)
        if return_solution:
            return solution


class ScipyDualAnneal(OptimBase):
    """A wrapper for scipy.optimize.dual_annealing

    _extended_summary_

    Parameters
    ----------
    OptimBase : _type_
        _description_
    """

    _inference_method = "Scipy Dual Annealing"

    def __init__(self,
                 model=None,
                 data=None,
                 verbose=False,
                 parameter_labels=None,
                 name='',
                 namestyle="full",
                 **kwargs):

        super().__init__(model=model, data=data, verbose=verbose,
                         parameter_labels=parameter_labels, name=name,
                         namestyle=namestyle, **kwargs)

    def optimize(self,
                 guess=None,
                 bounds=None,
                 maxiter=1000,
                 local_search_options=None,
                 initial_temp=5230,
                 restart_temp_ratio=2e-05,
                 visit=2.62,
                 accept=-5.0,
                 maxfun=1e+07,
                 seed=None,
                 no_local_search=False,
                 return_solution=False,
                 **kwargs):
        from scipy.optimize import dual_annealing
        self.logger.info("Starting optimisation")

        # DA requires a bounds argument to be passed
        # but our models may not be bounded
        # hence, we'll let DA run in the unit cube
        # and transform the parameters back to the
        # prior volume with the prior transform
        def neglnprob(*args):
            theta = self.prior_transform(*args)
            return -self.lnprob(theta)

        bounds = [(0, 1) for _ in range(self.npars)]
        if guess is None:
            # DA requires a guess to start from
            # but since we've already incorporated the prior transform,
            # we'll just use a random number in the unit cube
            self.logger.info("No inital guess provided")
            self.logger.info("Using random guess")
            guess = np.random.rand(self.npars)
        self.logger.info("starting from position: %s", guess)

        solution = dual_annealing(neglnprob, bounds,
                                  maxiter=maxiter,
                                  local_search_options=local_search_options,
                                  initial_temp=initial_temp,
                                  restart_temp_ratio=restart_temp_ratio,
                                  visit=visit,
                                  accept=accept,
                                  maxfun=maxfun,
                                  seed=seed,
                                  no_local_search=no_local_search,
                                  **kwargs)
        # transform back to the prior volume
        self.solution = self.prior_transform(*list(solution.x))
        self.logger.info("Optimisation completed in %s iterations",
                         solution.nit)
        self.logger.info("after %s function calls",
                         solution.nfev)
        if not solution.success:
            self.logger.warning("Optimisation did not converge")
            self.logger.warning("Message: %s", solution.message)
        if return_solution:
            return solution


class ScipyShgo(OptimBase):
    """A wrapper for scipy.optimize.shgo

    _extended_summary_

    Parameters
    ----------
    OptimBase : _type_
        _description_
    """

    _inference_method = "Scipy SHGO"

    def __init__(self,
                 model=None,
                 data=None,
                 verbose=False,
                 parameter_labels=None,
                 name='',
                 namestyle="full",
                 **kwargs):

        super().__init__(model=model, data=data, verbose=verbose,
                         parameter_labels=parameter_labels, name=name,
                         namestyle=namestyle, **kwargs)

    def optimize(self,
                 guess=None,
                 bounds=None,
                 n=100,
                 iters=1,
                 sampling_method='simplicial',
                 options=None,
                 minimizer_kwargs=None,
                 callback=None,
                 disp=False,
                 return_solution=False,
                 **kwargs):
        from scipy.optimize import shgo
        self.logger.info("Starting optimisation")

        # DA requires a bounds argument to be passed
        # but our models may not be bounded
        # hence, we'll let DA run in the unit cube
        # and transform the parameters back to the
        # prior volume with the prior transform
        def neglnprob(*args):
            theta = self.prior_transform(*args)
            return -self.lnprob(theta)

        bounds = [(0, 1) for _ in range(self.npars)]
        if guess is None:
            # DA requires a guess to start from
            # but since we've already incorporated the prior transform,
            # we'll just use a random number in the unit cube
            self.logger.info("No inital guess provided")
            self.logger.info("Using random guess")
            guess = np.random.rand(self.npars)
        self.logger.info("starting from position: %s", guess)

        solution = shgo(neglnprob, bounds,
                        n=n,
                        iters=iters,
                        sampling_method=sampling_method,
                        options=options,
                        minimizer_kwargs=minimizer_kwargs,
                        callback=callback,
                        disp=disp,
                        **kwargs)
        # transform back to the prior volume
        self.solution = self.prior_transform(*list(solution.x))
        self.logger.info("Optimisation completed in %s iterations",
                         solution.nit)
        self.logger.info("after %s function calls",
                         solution.nfev)
        if not solution.success:
            self.logger.warning("Optimisation did not converge")
            self.logger.warning("Message: %s", solution.message)
        if return_solution:
            return solution


class ScipyDirect(OptimBase):
    """A wrapper for scipy.optimize.direct

    _extended_summary_

    Parameters
    ----------
    OptimBase : _type_
        _description_
    """

    _inference_method = "Scipy DIRECT"

    def __init__(self,
                 model=None,
                 data=None,
                 verbose=False,
                 parameter_labels=None,
                 name='',
                 namestyle="full",
                 **kwargs):

        super().__init__(model=model, data=data, verbose=verbose,
                         parameter_labels=parameter_labels, name=name,
                         namestyle=namestyle, **kwargs)

    def optimize(self,
                 guess=None,
                 bounds=None,
                 eps=1e-4,
                 maxiter=1000,
                 maxfun=None,
                 locally_biased=True,
                 f_min=-np.inf,
                 f_min_rtol=0.0001,
                 vol_tol=1e-16,
                 len_tol=1e-06,
                 return_solution=False,
                 **kwargs):
        from scipy.optimize import direct
        self.logger.info("Starting optimisation")

        # DA requires a bounds argument to be passed
        # but our models may not be bounded
        # hence, we'll let DA run in the unit cube
        # and transform the parameters back to the
        # prior volume with the prior transform
        def neglnprob(*args):
            theta = self.prior_transform(*args)
            return -self.lnprob(theta)

        bounds = [(0, 1) for _ in range(self.npars)]
        if guess is None:
            # DA requires a guess to start from
            # but since we've already incorporated the prior transform,
            # we'll just use a random number in the unit cube
            self.logger.info("No inital guess provided")
            self.logger.info("Using random guess")
            guess = np.random.rand(self.npars)
        self.logger.info("starting from position: %s", guess)

        solution = direct(neglnprob, bounds,
                          maxiter=maxiter,
                          maxfun=maxfun,
                          eps=eps,
                          locally_biased=locally_biased,
                          f_min=f_min,
                          f_min_rtol=f_min_rtol,
                          vol_tol=vol_tol,
                          len_tol=len_tol,
                          **kwargs)
        # transform back to the prior volume
        self.solution = self.prior_transform(*list(solution.x))
        self.logger.info("Optimisation completed in %s iterations",
                         solution.nit)
        self.logger.info("after %s function calls",
                         solution.nfev)
        if not solution.success:
            self.logger.warning("Optimisation did not converge")
            self.logger.warning("Message: %s", solution.message)
        if return_solution:
            return solution


class AxOpt(OptimBase):
    """A wrapper for Ax

    Ax is a platform for 'Adaptive Experimentation' developed by Facebook.
    It is a Bayesian optimisation framework that can be used to optimise
    any function. It is particularly useful for optimising functions
    that are expensive to evaluate. It is also useful for optimising
    functions that are stochastic, as it can be used to optimise the
    expected value of the function.

    This provides access to the simple Loop API, which handles a lot
    of the boilerplate code for you. For more information, see
    https://ax.dev/tutorials/gpei_hartmann_loop.html. It will make
    automatic choices about which optimisation algorithms to use for
    you, and will also handle the parallelisation of the optimisation
    process.

    Parameters
    ----------
    model : instance of Model
        The model to be optimised
    data : list of instances of Data
        The data to be used in the optimisation
    verbose : bool, optional
        Whether to print verbose output, by default False
    parameter_labels : list of str, optional
        The names of the parameters, by default None
    name : str, optional
        The name of the optimisation, by default ''
    namestyle : str, optional
        The style of the name, by default "full"
    """

    _inference_method = "Ax"

    def __init__(self,
                 model=None,
                 data=None,
                 verbose=False,
                 parameter_labels=None,
                 name='',
                 namestyle="full",
                 **kwargs):

        super().__init__(model=model, data=data, verbose=verbose,
                         parameter_labels=parameter_labels, name=name,
                         namestyle=namestyle, **kwargs)

    def objective(self, parametrisation):
        u = np.array([parametrisation[p] for p in self.parLabels])
        theta = self.prior_transform(u)
        # Ax needs the mean and the standard error on the mean
        # for the objective function. For now, we'll assume
        # we're not working with stocahstic models, so the
        # standard error is zero.
        return (self.lnprob(theta), 0.0)

    def optimize(self,
                 total_trials=100,
                 return_solution=False,
                 ):
        """Optimise the model using Ax
        
        This function uses the Ax package to optimise the model. It
        uses the Bayesian optimisation algorithm to find the maximum
        of the posterior probability distribution. It is intended to
        be used for optimising functions that are expensive to evaluate.
        It is also useful for optimising functions that are stochastic,
        as it can be used to optimise the expected value of the function.
        
        Parameters
        ----------
        total_trials : int, optional
            The total number of trials to run, by default 100
        return_solution : bool, optional
            Whether to return the solution, by default False
            
        Returns
        -------
        dict
            A dictionary containing the solution, the lnprobs and the chain
        """
        self.logger.info("Starting optimisation")
        from ax.service.managed_loop import optimize
        # Ax requires a list of dictionaries for the parameter info
        # this includes a name, type and bounds at least
        # since this is bounded inference, we use the prior transform
        # to map the unit cube to the prior volume again.
        parameters = [{'name': p, 'type': 'range', 'bounds': [0, 1]}
                      for p in self.parLabels]
        best_parameters, values, experiment, model = optimize(
            parameters=parameters,
            evaluation_function=self.objective,
            objective_name='logP',
            experiment_name=self.name,
            total_trials=total_trials,
            minimize=False,
            )
        self.logger.info("Optimisation completed in %s iterations",
                         total_trials)
        self.solution = self.prior_transform(*[best_parameters[p]
                                               for p in self.parLabels])
        self.opt_model = model
        self.opt_experiment = experiment
        self.opt_values = values
        # now lets try to construct something approximating a chain from this:
        trials = experiment.trials
        self.chain = np.array([self.prior_transform(*[trials[i].arm.parameters[pn]
                                                      for pn in self.parLabels
                                                      ])
                               for i in range(len(trials))
                               ])
        self.lnprobs = [trials[i].objective_mean for i in range(len(trials))]
        if return_solution:
            return {'solution': self.solution,
                    'lnprobs': self.lnprobs,
                    'chain': self.chain}

    def optimize_with_restarts(self,
                               n_restarts=10,
                               total_trials=100,
                               **kwargs):
        """Optimise the model using Ax with restarts

        This function repeatedly calls optimize() to perform a number
        of BO runs. It then returns the best solution found. The results are
        stored in the all_chains and all_lnprobs attributes, rather than being
        returned. The best solution is stored in the solution attribute.

        Parameters
        ----------
        n_restarts : int, optional
            The number of restarts to perform, by default 10
        total_trials : int, optional
            The total number of trials to run for each restart, by default 100

        Returns
        -------
        None
        """
        self.all_chains = np.zeros((n_restarts, total_trials, self.npars))
        self.all_lnprobs = np.zeros((n_restarts, total_trials))
        for i in range(n_restarts):
            self.logger.info("Optimisation restart %s", i)
            solution = self.optimize(total_trials=total_trials,
                                     return_solution=True,
                                     **kwargs)
            self.all_chains[i, :, :] = self.chain
            self.all_lnprobs[i, :] = self.lnprobs
            # now we have to check if this is the best solution so far
            # and add the final solution to the chain
            if i == 0:
                self.bestPars = self.solution
                self.bestLnprob = self.opt_values[0]  # solution['lnprobs']
                self.chain = [(self.bestPars, self.bestLnprob)]
                i_best = 0
            else:
                if solution.fun > self.bestLnprob:
                    self.bestPars = self.solution
                    self.bestLnprob = solution.fun
                    i_best = i
                self.chain.append((self.solution, solution.fun))
            self.logger.info("Best solution so far: %s", self.bestPars)
            self.logger.info("Found at restart %s", i_best)
            self.logger.info("With lnprob: %s", self.bestLnprob)
        # now we have to set the solution to the best one
        self.solution = self.bestPars  # so that get_map works
        self.logger.info("Optimisation complete")


class AxBO(AxOpt):
    """A wrapper for Bayesian Optimisation with Ax

    This is actually just an alias for AxOpt, but it is provided for
    ease, and to make it obvious that we are using BO. It is intended 
    to be used for optimising functions that are expensive to evaluate. 
    It is also useful for optimising functions that are stochastic, as 
    it can be used to optimise the expected value of the function.

    This provides access to the simple Loop API, which handles a lot
    of the boilerplate code for you. For more information, see
    https://ax.dev/tutorials/gpei_hartmann_loop.html. It will make
    automatic choices about which optimisation algorithms to use for
    you, and will also handle the parallelisation of the optimisation
    process."""

    _inference_method = "Ax Bayesian Optimisation"


class AxSAASBO(AxOpt):
    """A wrapper for Sparse Axis-Aligned Subspace Bayesian Optimisation with Ax

    This class performs Bayesian Optimisation using the SAAS algorithm.
    It is intended to be used for higher-dimensional problems than
    standard BO, up to a few hundred dimensions. For example, Eriksson
    & Jankowiak (2021) showed success with a 388 dimensional SVM problem.

    Unlike standard BO, SAASBO cannot use the managed loop API, and
    so it is necessary to write a custom loop. This is provided in the
    optimise() method. The loop is based on the example provided in
    https://ax.dev/tutorials/saasbo.html.

    Parameters
    ----------
    OptimBase : _type_
        _description_
    """

    _inference_method = "Ax SAAS Bayesian Optimisation"
    _citation = """
@InProceedings{pmlr-v161-eriksson21a,
  title = 	 {High-dimensional {Bayesian} optimization with sparse axis-aligned subspaces},
  author =       {Eriksson, David and Jankowiak, Martin},
  booktitle = 	 {Proceedings of the Thirty-Seventh Conference on Uncertainty in Artificial Intelligence},
  pages = 	 {493--503},
  year = 	 {2021},
  editor = 	 {de Campos, Cassio and Maathuis, Marloes H.},
  volume = 	 {161},
  series = 	 {Proceedings of Machine Learning Research},
  month = 	 {27--30 Jul},
  publisher =    {PMLR},
  pdf = 	 {https://proceedings.mlr.press/v161/eriksson21a/eriksson21a.pdf},
  url = 	 {https://proceedings.mlr.press/v161/eriksson21a.html},
  abstract = 	 {Bayesian optimization (BO) is a powerful paradigm for efficient optimization of black-box objective functions. High-dimensional BO presents a particular challenge, in part because the curse of dimensionality makes it difficult to define—as well as do inference over—a suitable class of surrogate models. We argue that Gaussian process surrogate models defined on sparse axis-aligned subspaces offer an attractive compromise between flexibility and parsimony. We demonstrate that our approach, which relies on Hamiltonian Monte Carlo for inference, can rapidly identify sparse subspaces relevant to modeling the unknown objective function, enabling sample-efficient high-dimensional BO. In an extensive suite of experiments comparing to existing methods for high-dimensional BO we demonstrate that our algorithm, Sparse Axis-Aligned Subspace BO (SAASBO), achieves excellent performance on several synthetic and real-world problems without the need to set problem-specific hyperparameters.}
}
"""

    def optimize(self, total_trials=100, return_solution=False):
        """optimize _summary_

        _extended_summary_

        Parameters
        ----------
        total_trials : int, optional
            _description_, by default 100
        return_solution : bool, optional
            _description_, by default False
        """

        from ax import Data, Experiment, ParameterType, RangeParameter, SearchSpace
        from ax.modelbridge.generation_strategy import GenerationStep, GenerationStrategy
        from ax.modelbridge.registry import Models
        from ax.runners.synthetic import SyntheticRunner
        from ax.core.metric import Metric
        from ax.core.objective import Objective
        from ax.core.optimization_config import OptimizationConfig

        self.logger.info("Starting optimisation")
        # Ax requires a list of dictionaries for the parameter info



class AxBOParallel(AxOpt):
    """Perform Bayesian Optimisation with Ax in parallel using Ray
    
    This class performs Bayesian Optimisation using the Ax package
    in parallel using the Ray package. It is intended to be used for
    optimising functions that are so expensive to evaluate that you
    need to distribute the evaluations to complete in a reasonable 
    time. As a result, it will make more evaluations in total than
    regular BO would, but will hopefully take less time. It is also
    useful for optimising functions that are stochastic, as it can
    be used to optimise the expected value of the function."""