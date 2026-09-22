"""W5.25: ``halfcauchy`` agrees against ``scipy`` on every registered backend.

``regularised_horseshoe``'s global scale (W5.8) is half-Cauchy under *both*
of its tails (``ampere/core/parameter.py``), and until this item ``halfcauchy``
was in neither modern backend's ``lowering.md`` §3.2 table -- reference-only
on NUTS. This is the row §3 asks for: one declaration, ``loc == 0`` (native on
both backends) and shifted (§3.3's composed route), checked against ``scipy``
on every registered fixture, in ``test_parameters.py``'s
``TestPriorsAgainstScipy`` shape.

Nothing here compares ampere against ampere: every oracle is ``scipy``, as
``tests/conformance/oracles.py``'s module docstring requires.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.stats as st

from ampere.core import Parameter, ParameterSet

from .oracles import analytic_lnprior, analytic_prior_transform, free_slices
from .protocol import ConformanceBackend, ParameterSpace, Tolerances

#: ``loc == 0`` and a shift, in one declaration per case: consistency across
#: the two is exactly what "the chain now lowers" (``lowering.md`` §3.2.1) is
#: a claim about, since a shifted half-Cauchy takes the composed route (§3.3)
#: rather than the native one.
LOCS: tuple[float, ...] = (0.0, 3.0)

CUBE = np.array([0.42])


def declaration(loc: float) -> ParameterSet:
    return ParameterSet([Parameter("scale", st.halfcauchy(loc, 2.0))])


@pytest.fixture(params=LOCS, ids=[f"loc={loc:g}" for loc in LOCS])
def space(request: pytest.FixtureRequest, backend: ConformanceBackend) -> ParameterSpace:
    return backend.parameter_space(declaration(request.param))


class TestHalfCauchyAgreesWithScipy:
    """§3.4's fallback used as intended: an exact construction, not an approximation."""

    def test_prior_transform_is_the_analytic_quantile(
        self, space: ParameterSpace, tolerances: Tolerances
    ) -> None:
        """The sample: the same quantile ``scipy``'s own ``ppf`` gives."""
        expected = analytic_prior_transform(space.declaration, CUBE)
        assert space.prior_transform(CUBE) == pytest.approx(expected, abs=tolerances.analytic)

    def test_lnprior_is_the_scipy_logpdf(
        self, space: ParameterSpace, tolerances: Tolerances
    ) -> None:
        """The log-density: summed ``scipy`` ``logpdf`` at the transformed sample."""
        theta = analytic_prior_transform(space.declaration, CUBE)
        expected = analytic_lnprior(space.declaration, theta)
        assert space.lnprior(theta) == pytest.approx(expected, abs=tolerances.analytic)

    def test_a_point_below_loc_has_no_prior_mass(self, space: ParameterSpace) -> None:
        """The sample's support: nothing below ``loc``, shifted or not.

        Exercises exactly the row ``lowering.md``'s module docstring calls out
        for the composed route -- a lowering that reported ``Real()`` instead
        of ``greater_than(loc)`` would pass every density check above and
        still let a sampler spend proposals below the support.
        """
        [(parameter, _)] = free_slices(space.declaration)
        loc, _ = parameter.prior.support()
        below = np.array([loc - 1.0])
        assert space.lnprior(below) == -np.inf

    def test_the_centre_of_the_cube_is_the_prior_median(
        self, space: ParameterSpace, tolerances: Tolerances
    ) -> None:
        centre = np.full(space.free_size, 0.5)
        [(parameter, _)] = free_slices(space.declaration)
        assert space.prior_transform(centre) == pytest.approx(
            [parameter.prior.median()], abs=tolerances.analytic
        )
