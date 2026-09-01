from __future__ import print_function

import numpy as np

# NOTE (W0.5 / issue #74): this used to be an absolute `import models` /
# `from models import Model`, which only worked when ampere/models/ was
# importable as a top-level package (pre-restructure layout). It is now a
# subpackage, so the import must be relative.
from .models import Model

# NOTE (W0.5 / issue #75): the 'extinction' package this module used to
# depend on (`from extinction import apply as redden, fitzpatrick99 as f99`)
# is abandoned upstream and is not an ampere dependency. We replace it with
# 'dust_extinction' (astropy-affiliated, actively maintained), imported
# lazily inside F99Extinction.__init__ so that merely importing this module
# never requires it -- see the lazy-import quarantine note there.


class BaseExtinctionLaws(Model):
    """

    """
    
    def __init__(self, **kwargs):
        raise NotImplementedError("")
    
    def __str__(self, **kwargs):
        raise NotImplementedError("")
    
    def __repr__(self, **kwargs):
        raise NotImplementedError("")

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def extinctionModel(self, **kwargs): 
        raise NotImplementedError("")
    

class F99Extinction(BaseExtinctionLaws):
    """Fitzpatrick (1999) extinction law.

    NOTE (W0.5 / issue #75): the abandoned 'extinction' package this class
    used to depend on (`extinction.fitzpatrick99` for the curve,
    `extinction.apply` to turn magnitudes into a multiplicative flux
    factor) is replaced by 'dust_extinction' (astropy-affiliated,
    actively maintained; optional extra: ``ampere[extinction]``), imported
    lazily in __init__ so importing this module never requires it.

    Unit-convention note: the old code computed
    ``extinction.apply(fitzpatrick99(1/wavelengths_um, av, rv, unit='invum'),
    ones)``, which -- because ``apply(mag, flux) == flux * 10**(-0.4*mag)``
    applied to an array of ones -- reduces to the plain transmission
    fraction ``10**(-0.4 * A_lambda)``. ``dust_extinction``'s
    ``<model>.extinguish(x, Av=av)`` is defined as exactly
    ``10**(-0.4 * Av * <model>(x))`` where ``<model>(x)`` is the
    dimensionless ``A(x)/A(V)`` curve, i.e. the same transmission-fraction
    convention -- so ``self.modelFlux`` has the same meaning as before.
    Not independently re-verified against the old 'extinction' package's
    numerical output (it is abandoned/unavailable to compare against in
    this environment); flagged here for review since no test exercises
    this module.
    """

    def __init__(self, wavelengths, av = None, rv = None, flatprior = True,
                 lims = np.array([[0,100],[0,10]]), #limits on av and rv respectively
                 **kwargs):
        try:
            from dust_extinction.parameter_averages import F99
        except ImportError as exc:
            raise ImportError(
                "ampere.models.extinctionModels.F99Extinction requires the "
                "optional 'dust_extinction' package (replacement for the "
                "abandoned 'extinction' package -- see ampere issue #75), "
                "which is not installed in this environment. Install it "
                "via the 'ampere[extinction]' extra. "
                f"Original import error: {exc}"
            ) from exc
        self._F99_cls = F99 #define the extinction model class, this is more to provide convenience for when people want to evaluate the model themselves
        self.wavelengths = wavelengths #Wavelengths are assumed to be in micron!!
        self.invwaves = 1/wavelengths
        self.flatprior = flatprior
        self.lims = lims
        self.npars = 2 #Up to two free parameters for this type of model
        # self.rv/self.av default to None so the "is it fixed?" checks below
        # (here and in __call__/lnprior) are well-defined even when only one
        # of the two is fixed (pre-existing bug fixed alongside #74: these
        # attributes were previously left undefined -- AttributeError -- in
        # that case).
        self.rv = None
        self.av = None
        if rv is not None: #R_v will be fixed in the fitting
            self.rv = rv
            self.npars -= 1

        if av is not None: #A_v will be fixed in the fitting
            self.av = av
            self.npars -= 1

        #raise NotImplementedError("")

    def __str__(self, **kwargs):
        raise NotImplementedError("")

    def __repr__(self, **kwargs):
        raise NotImplementedError("")

    def __call__(self, *args, **kwargs):
        from astropy import units as u
        if self.npars == 2:
            #both av and rv are free parameters
            av = args[0]
            rv = args[1]
        elif self.npars == 1:
            #one of av or rv is fixed
            if self.rv is not None: #Rv is fixed
                av = args[0]
                rv = self.rv
            elif self.av is not None:
                av = self.av
                rv = args[0]
        elif self.npars == 0:
            #both av and rv are fixed, we just want to return extinction values
            av = self.av
            rv = self.rv
        #dust_extinction's extinguish() returns 10**(-0.4*Av*curve(x)), i.e.
        #the same transmission-fraction convention as the old
        #extinction.apply(fitzpatrick99(...), ones) -- see the class
        #docstring.
        self.modelFlux = self._F99_cls(Rv=rv).extinguish(self.wavelengths * u.micron, Av=av)


    def lnprior(self, theta, **kwargs):
        if self.npars == 2:
            #both av and rv are free parameters
            if self.flatprior:
                if (self.lims[0,0] <= theta[0] <= self.lims[0,1]) and (self.lims[1,0] <= theta[1] <= self.lims[1,1]):
                    return 0
                else:
                    return -np.inf
            pass
        elif self.npars == 1:
            #one of av or rv is fixed
            if self.av is not None:
                if (self.lims[0,0] <= theta[0] <= self.lims[0,1]):
                    return 0
                else:
                    return -np.inf
            elif self.rv is not None:
                if (self.lims[1,0] <= theta[0] <= self.lims[1,1]):
                    return 0
                else:
                    return -np.inf
            pass
        elif self.npars == 0:
            #both av and rv are fixed, we just want to return a constant because the prior is meaningless
            return 0
        pass


    def prior_transform(self, u, **kwargs):
        pass
    
        
class CCMExtinctionLaw(BaseExtinctionLaws):
    """Cardelli, Clayton & Mathis (1989) extinction law.

    NOTE (W0.5 / issue #76): this class was left unfinished -- ``__init__``
    unconditionally raised ``NotImplementedError("")`` (never accepted or
    stored wavelength/A_V/R_V state), and ``extinctionModel``/``CCMModel``
    referenced undefined names (``wavelength``, ``R_A``, ``A_V``,
    ``type_Wavelength``) and called ``CCMModel`` as a bare function instead
    of ``self.CCMModel``. That is a design gap (how to pick
    ``type_Wavelength`` per input wavelength, and what the constructor's
    signature should be), not a syntax error, so it is out of scope to
    design here -- the class is quarantined: it raises an informative
    ``NotImplementedError`` on construction instead of the original bare
    one. What *is* fixed here are the actual syntax/logic bugs reported by
    #76 in the ``a()``/``b()`` CCM89 polynomial-coefficient helper methods
    (missing multiplication operators that astropy/CPython's parser had
    been silently reinterpreting as call expressions, an
    ``X == "NIR" or "Optical"`` truthiness bug that made the NIR/Optical
    branch always match, and an ``A_Lambda``/``A_lambda`` case typo in
    ``CCMModel``), so they now correctly implement the published CCM89
    coefficients (verified against Cardelli, Clayton & Mathis 1989, ApJ,
    345, 245) and are usable directly as pure functions of
    (wavelength, type_Wavelength) for anyone completing this model.
    """

    def __init__(self, **kwargs):
        raise NotImplementedError(
            "CCMExtinctionLaw was never completed (ampere issue #76): its "
            "constructor never accepted/stored wavelength, A_V or R_V, and "
            "extinctionModel()/CCMModel() reference undefined names. W0.5 "
            "fixed the reported syntax bugs in the a()/b() CCM89 "
            "coefficient helpers (now correct, usable standalone) but did "
            "not design the missing constructor/orchestration -- that is "
            "out of scope for a quarantine fix. Use F99Extinction for a "
            "working extinction law, or complete this class's __init__, "
            "extinctionModel() and CCMModel() (in particular: pick "
            "type_Wavelength from the input wavelength range and pass it "
            "through, and call self.CCMModel(...) rather than the bare "
            "name)."
        )

    def extinctionModel(self, extinctionModel, **kwargs):
        raise NotImplementedError(
            "CCMExtinctionLaw.extinctionModel is not implemented -- see "
            "CCMExtinctionLaw.__init__'s docstring/error (ampere issue #76)."
        )

    def a(self, wavelength, type_Wavelength, **kwargs):

        inverse_wavelength = 1 / wavelength

        if type_Wavelength == "IR":

                a = 0.574 * (inverse_wavelength**1.61)

        elif type_Wavelength in ("NIR", "Optical"):

                y = inverse_wavelength - 1.82

                a = 1 + (0.17699*y) - (0.5044*(y**2)) - (0.02427*(y**3)) + (0.72085*(y**4)) + (0.01979*(y**5)) - (0.77530*(y**6)) + (0.32999*(y**7))

        elif type_Wavelength == "UV":

                if 8 >= inverse_wavelength >= 5.9:

                        P_a = -(0.04473 * ((inverse_wavelength - 5.9)**2)) - (0.009779 * ((inverse_wavelength - 5.9)**3))

                elif inverse_wavelength < 5.9:

                        P_a = 0

                a = 1.752 - (0.316*inverse_wavelength) - (0.104 / ( ((inverse_wavelength - 4.67)**2) + 0.341) ) + P_a

        elif type_Wavelength == "FUV":

                a = -1.073 - (0.628*(inverse_wavelength - 8)) + (0.137 * ((inverse_wavelength - 8)**2)) - (0.070 * ((inverse_wavelength - 8)**3))

        else:
                raise NotImplementedError("Incorrect type_Wavelength string")


        return a



    def b(self, wavelength, type_Wavelength, **kwargs):

            inverse_wavelength = 1 / wavelength

            if type_Wavelength == "IR":

                    b = 0.527 * (inverse_wavelength**1.61)

            elif type_Wavelength in ("NIR", "Optical"):

                    y = inverse_wavelength - 1.82

                    b = (1.41338*y) - (2.28305*(y**2)) - (1.07233*(y**3)) - (5.38434*(y**4)) - (0.62251*(y**5)) + (5.30260*(y**6)) - (2.09002*(y**7))

            elif type_Wavelength == "UV":

                    if 8 >= inverse_wavelength >= 5.9:

                            P_b = (0.2130 * ((inverse_wavelength - 5.9)**2)) + (0.1207 * ((inverse_wavelength - 5.9)**3))

                    elif inverse_wavelength < 5.9:

                            P_b = 0

                    b = (-3.090) + (1.825*inverse_wavelength) + (1.206 / ( ((inverse_wavelength - 4.62)**2) + 0.263) ) + P_b

            elif type_Wavelength == "FUV":

                    b = 13.670 + (4.257*(inverse_wavelength - 8)) - (0.420 * ((inverse_wavelength - 8)**2)) - (0.374 * ((inverse_wavelength - 8)**3))

            else:
                raise NotImplementedError("Incorrect type_Wavelength string")

            return b

    def CCMModel(self, wavelength, type_Wavelength, R_A, A_V, **kwargs):

        a_lambda = self.a(wavelength, type_Wavelength, **kwargs)
        b_lambda = self.b(wavelength, type_Wavelength, **kwargs)

        A_lambda = ( a_lambda + (b_lambda / R_A) ) * A_V
        return A_lambda

