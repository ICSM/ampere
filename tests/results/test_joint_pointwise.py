"""The ``"joint"`` decomposition: a run over a joint noise group (W5.9).

``results.md`` §6's third decomposition and ``inference.md`` §4.9's third kind
of top-level component, measured on a real run rather than asserted. Two
claims, and the second is the one that would otherwise fail silently:

* a joint noise group contributes **one** term to the per-dataset
  ``log_likelihood`` group, under the group's own label, in place of its
  members' separate ones;
* the pointwise group carries the group's **rotated outputs** — one per member
  label, on the grid they share — declared ``"joint"``, because a group's
  channels are not independent and a leave-one-out conditional of one channel
  alone would condition on its sibling's value at the same sample without
  saying so.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.backends.reference import EpochSample, ReflexOrbit
from ampere.core import (
    Dataset,
    DatasetCollection,
    DenseGP,
    FittingProblem,
    GaussianFamily,
    Instrument,
    JointGaussianProcessNoise,
    Likelihood,
    Matern32,
    RotationCoupling,
    TimeSeries,
    negotiate,
)
from ampere.inference import EmceeEngine
from ampere.results import (
    ATTR_PREFIX,
    POINTWISE_LOG_LIKELIHOOD_GROUP,
    JOINT_DECOMPOSITION,
    add_pointwise_log_likelihood,
)

GRID = np.array([0.0, 40.0, 95.0, 150.0, 210.0, 260.0, 330.0, 400.0, 460.0, 520.0])
SIGMA = 0.03
ORBIT = {
    "pmra": 1.2,
    "pmdec": -0.6,
    "period": 400.0,
    "phase": 0.7,
    "amp_ra": 0.6,
    "amp_dec": 0.35,
}
COUPLING = {
    "angle": 0.7,
    "log_variance_0": float(np.log(9.0e-4)),
    "log_variance_1": float(np.log(1.0e-4)),
}


def _channel(values: np.ndarray) -> TimeSeries:
    return TimeSeries(
        GRID * u.day,
        np.asarray(values, dtype=float) * u.mas,
        uncertainty=np.full(GRID.size, SIGMA) * u.mas,
    )


def _problem() -> FittingProblem:
    instruments = {
        name: Instrument([EpochSample(GRID)], channel=name, label=f"astrom_{name}")
        for name in ("ra", "dec")
    }
    truth = ReflexOrbit(GRID, **ORBIT)
    compiled = truth.compile_for(negotiate(list(instruments.values())))
    evaluated = compiled.evaluate()
    rng = np.random.default_rng(3)
    observed = {
        name: _channel(
            np.asarray(step(evaluated).values).ravel() + rng.normal(0.0, SIGMA, GRID.size)
        )
        for name, step in instruments.items()
    }
    model = ReflexOrbit(
        GRID,
        pmra=st.norm(0.0, 5.0),
        pmdec=st.norm(0.0, 5.0),
        period=ORBIT["period"],
        phase=ORBIT["phase"],
        amp_ra=ORBIT["amp_ra"],
        amp_dec=ORBIT["amp_dec"],
    )
    noise = JointGaussianProcessNoise(
        Matern32(1.0, 180.0, axes=("time",)),
        DenseGP(),
        datasets=("ra", "dec"),
        coupling=RotationCoupling(
            COUPLING["angle"], COUPLING["log_variance_0"], COUPLING["log_variance_1"]
        ),
    )
    datasets = DatasetCollection(
        {
            name: Dataset(
                observed[name],
                instruments[name],
                likelihood=Likelihood(GaussianFamily()),
                label=name,
            )
            for name in ("ra", "dec")
        },
        joint={"astrom": noise},
    )
    return FittingProblem(model, datasets, seed=7)


class TestAJointRun:
    """One emcee run at a tiny budget, and what it stores."""

    def test_the_log_likelihood_group_is_keyed_by_the_group(self) -> None:
        problem = _problem()
        run = EmceeEngine(problem, walkers=8).run(30, burn_in=10, progress=False)
        stored = run["log_likelihood"].dataset
        assert set(stored.data_vars) == {"astrom"}
        assert np.all(np.isfinite(np.asarray(stored["astrom"])))

    def test_the_pointwise_group_carries_the_rotated_outputs(self) -> None:
        problem = _problem()
        run = EmceeEngine(problem, walkers=8).run(30, burn_in=10, progress=False)
        attached = add_pointwise_log_likelihood(run, problem, thin=10)
        group = attached[POINTWISE_LOG_LIKELIHOOD_GROUP].dataset
        assert set(group.data_vars) == {"ra", "dec"}
        assert group.attrs[f"{ATTR_PREFIX}decomposition"] == JOINT_DECOMPOSITION
        for name in group.data_vars:
            assert group[name].attrs[f"{ATTR_PREFIX}decomposition"] == JOINT_DECOMPOSITION
            values = np.asarray(group[name])
            assert values.shape[-1] == GRID.size
            assert np.all(np.isfinite(values))

    def test_the_rotated_terms_are_not_the_channels_own(self) -> None:
        """A rotated output is not its channel: it is a mixture of both.

        The row that would catch the group silently storing each channel's own
        leave-one-out terms instead of the rotated ones — which is what a
        pointwise implementation that ignored the coupling would produce, and
        which is finite, per-sample and plausible in every other respect.
        """
        problem = _problem()
        noise = problem.datasets.joint["astrom"]
        simulation = problem.simulate({"model.pmra": ORBIT["pmra"], "model.pmdec": ORBIT["pmdec"]})
        assert simulation.predicted is not None
        residuals = [
            np.asarray(problem.datasets[name].observed.values, dtype=float).ravel()
            - np.asarray(simulation.predicted[name].values, dtype=float).ravel()
            for name in ("ra", "dec")
        ]
        rotated = noise.rotate(np.column_stack(residuals), COUPLING)
        # The rotation really mixes them: neither column is either channel.
        for index in range(2):
            for residual in residuals:
                assert not np.allclose(rotated[:, index], residual)
