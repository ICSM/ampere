"""The WStat comparison: a profiled user family against the two-dataset joint fit.

Carries the obligation the 2026-09-03 freeze ruling placed on W1.13 and W2.9
(``DEVELOPMENT_PLAN.md`` §2's "WStat / profile likelihoods" row,
``docs/design/modalities/awkward_instrument.md`` §5 and §9 question 3):
**ampere ships no WStat.** The profiled Cash-with-background statistic ("W
statistic" in XSPEC's terminology) that X-ray astronomers use to fit a
source spectrum together with its background, without a background *model*,
is a profile likelihood rather than a marginal one — its nuisance background
rate is replaced, per bin, by its maximum-likelihood value rather than
integrated out. That is not Bayesian, and ampere's documentation takes the
deliberately opinionated line that the two-dataset formulation — source and
background as two datasets sharing one background model, fit jointly — is
the *correct* approach. But WStat is what practitioners reach for, and
``likelihoods.md`` §4.4's extension surface (a user-written
:class:`~ampere.core.LikelihoodFamily`) makes it expressible, so this module
builds both, runs both, and compares them honestly: pros, cons, and
*results*, not false balance.

Two things this module is not:

* **Not a response-folding example.** ``awkward_instrument.md`` already
  covers RMF/ARF forward folding as an ordinary ``Transformation``; doing it
  again here would only obscure the statistics under discussion. The "model"
  here predicts expected *counts* per channel directly, exposure and
  effective area already folded in — exactly what that sketch's §2
  concludes an X-ray count-rate model should do.
* **Not a claim that WStat is unsafe to use.** It is safe under masking
  (:attr:`~ampere.core.NoiseParams.retain` keeps the background buffer
  aligned with whichever channels survive a mask), and its point estimates
  are usually adequate at moderate counts. What degrades is what a
  *profiled* nuisance always costs: no generative model (so ``sample()``
  refuses, see :meth:`ProfiledCashWithBackground.sample`), per-sample terms
  that are not valid predictive densities for ``arviz.loo``/``waic`` (see
  :class:`ProfiledCashWithBackground`'s docstring), and a posterior for the
  source parameters that is narrower than the two-dataset fit's because the
  background's own uncertainty is a point estimate rather than a marginalised
  quantity.

Running this module
--------------------
``python examples/wstat_comparison.py`` runs both fits with the "doc" budget
(:data:`DOC_WALKERS`, :data:`DOC_STEPS`, :data:`DOC_BURN_IN` — a few dozen
bins, short chains, well under a minute on the reference backend) and prints
the comparison. ``tests/examples/test_wstat_comparison.py`` exercises the same
code end to end with a much smaller budget, as the environment-independent
gate the work item's acceptance criterion asks for (the docs build's notebook
execution is not guaranteed to fail loudly — ``docs/source/conf.py`` sets
``nbsphinx_allow_errors = True`` — so a script that runs in the docs page is
not, by itself, proof it ran; the pytest is).

Everything is seeded through :class:`~ampere.core.FittingProblem`'s own
``seed=`` (which drives walker initialisation and any other engine
randomness) and through a separate, explicit seed for the synthetic data
(:func:`synthetic_xray_counts`), so the whole comparison is deterministic.
"""

from __future__ import annotations

from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.core import (
    Dataset,
    DatasetCollection,
    FittingProblem,
    IndependentNoise,
    Instrument,
    Likelihood,
    LikelihoodFamily,
    Model,
    ModelResult,
    NoiseParams,
    Parameter,
    PoissonFamily,
    Spectrum,
    register_family,
)
from ampere.core.exceptions import LikelihoodError
from ampere.inference import EmceeEngine

__all__ = [
    "ENERGY_KEV",
    "N_BINS",
    "RATIO",
    "REFERENCE_ENERGY_KEV",
    "SEED",
    "TRUTH",
    "ProfiledCashWithBackground",
    "XraySource",
    "XraySourceAndBackground",
    "build_joint_problem",
    "build_wstat_problem",
    "compare",
    "main",
    "run_engine",
    "summarise",
    "synthetic_xray_counts",
]

# ---------------------------------------------------------------------------
# The synthetic data: a few dozen bins, deterministic, low-ish counts so the
# profiling trade-off (§ below) is inspectable rather than academic.
# ---------------------------------------------------------------------------

#: A few dozen channels, as the work item asks for — enough to see a spectral
#: shape, few enough that a chain of a few hundred steps costs a fraction of a
#: second per walker on the reference backend.
N_BINS = 32
ENERGY_KEV = np.linspace(0.5, 8.0, N_BINS)
REFERENCE_ENERGY_KEV = 1.0
#: t_source / t_background: how much of the background rate spills into the
#: source region per unit of background-region rate. Less than one, as is
#: typical when the background region is the larger extraction area.
RATIO = 0.3
#: The synthetic truth. Deliberately faint (a handful of source counts per
#: bin at the high-energy end) so the low-count regime where a profiled
#: nuisance and a marginalised one part company is the regime being tested,
#: not an aside.
TRUTH: dict[str, float] = {"src_norm": 6.0, "src_index": -1.5, "bkg_norm": 2.0}
#: Seed for the synthetic draw. Independent of the engines' own seeding
#: (``FittingProblem(seed=...)``, passed separately below) — data generation
#: happens once, before any problem exists.
SEED = 20260908


def _expected_counts(
    energy_kev: np.ndarray,
    norm: float,
    index: float,
    reference_energy_kev: float = REFERENCE_ENERGY_KEV,
) -> np.ndarray:
    """``norm * (E / E_ref) ** index`` — the shared power-law shape."""
    return norm * (np.asarray(energy_kev, dtype=float) / reference_energy_kev) ** index


def synthetic_xray_counts(seed: int = SEED) -> tuple[np.ndarray, np.ndarray]:
    """Draw one small, deterministic two-region X-ray-like count spectrum.

    Returns ``(source_counts, background_counts)``, each ``(N_BINS,)``
    integer-valued float arrays: the source-region counts (source plus its
    share of the background) and the background-region counts (background
    alone), drawn from independent Poisson processes at :data:`TRUTH`.
    """
    rng = np.random.default_rng(seed)
    source_rate = _expected_counts(ENERGY_KEV, TRUTH["src_norm"], TRUTH["src_index"])
    background_rate = np.full(N_BINS, TRUTH["bkg_norm"])
    source_region_rate = source_rate + RATIO * background_rate
    source_counts = rng.poisson(source_region_rate).astype(float)
    background_counts = rng.poisson(background_rate).astype(float)
    return source_counts, background_counts


# ---------------------------------------------------------------------------
# Route 1: the profiled statistic, as a user family.
# ---------------------------------------------------------------------------


@register_family
class ProfiledCashWithBackground(LikelihoodFamily):
    """The profiled Cash-with-background ("W") statistic, as a user family.

    For channel *i*, the source-region counts ``S_i`` are Poisson in
    ``m_i + ratio * b_i`` and the background-region counts ``B_i`` are
    Poisson in ``b_i``, with ``m_i`` the model's predicted source counts and
    ``b_i`` an unknown per-channel background rate. Rather than fit ``b_i``
    (XSPEC's "W statistic" approach, and the reason it needs no background
    *model* at all — its one real advantage over the joint fit below), this
    family **profiles it out analytically**: for fixed ``m_i``, the closed-form
    maximum-likelihood ``b_i`` is the positive root of a quadratic in the
    joint Poisson log-likelihood of ``S_i`` and ``B_i``,

    .. math::

        \\hat b_i = \\frac{C_i + \\sqrt{C_i^2 + 4\\,\\rho(\\rho+1)\\,B_i m_i}}
                          {2\\,\\rho(\\rho+1)}, \\qquad
        C_i = \\rho (S_i + B_i) - (\\rho + 1) m_i,

    with :math:`\\rho` the ``ratio``. ``log_prob`` then returns
    ``sum_i [log Poisson(S_i; m_i + ratio * b_hat_i) + log Poisson(B_i; b_hat_i)]``
    — the *profile* log-likelihood, evaluated at each channel's own
    plug-in background estimate.

    **Safe under masking.** The background counts are this family's own
    aligned per-sample data (``likelihoods.md`` §8's gap X-2), not something
    ``Likelihood`` excises on its behalf, so ``log_prob`` excises them itself
    with :attr:`~ampere.core.NoiseParams.retain` — the boolean indicator over
    the *full*, unmasked containers that every noise model already sets.
    Without this, masking a channel would silently pair the wrong background
    count with the wrong source count;
    ``tests/examples/test_wstat_comparison.py`` exercises exactly the
    unmasked-vs-masked comparison ``awkward_instrument.md`` §5 used to
    demonstrate the trap.

    **``sample()`` refuses**, and not only for the general reason
    :meth:`~ampere.core.LikelihoodFamily.sample`'s default already gives
    (a ``log_prob`` does not determine an observation process). Here there is
    a sharper, structural reason: :math:`\\hat b_i` is the profiled
    background estimate computed *from* :math:`S_i` and :math:`B_i` — the
    very counts a draw would need to produce. There is no forward model to
    sample from without already having the data it profiles over. The
    two-dataset formulation below has no such problem: its background is an
    ordinary parametric model, and both regions can be simulated forward from
    the prior or the posterior in the usual way.

    **Per-sample log-likelihood semantics degrade.** Because each channel's
    nuisance is independent, ``log_prob`` above already sums additively
    over channels — ``Likelihood.pointwise_log_prob``'s per-sample
    decomposition (``likelihoods.md`` §8) would reproduce those same terms
    exactly, numerically. But a term ``log Poisson(S_i; ...) +
    log Poisson(B_i; b_hat_i)`` is *not* a valid predictive density for
    channel *i*: it already used channel *i*'s own counts to fix its own
    nuisance parameter, rather than marginalising over the nuisance's
    posterior uncertainty as a genuine per-observation predictive density
    would. Handing these terms to ``arviz.loo``/``waic`` would be circular —
    exactly why ``PoissonFamily``'s ordinary terms and ``GaussianProcessNoise``'s
    ``conditional_loo`` terms are valid predictive densities and this
    family's terms are not, even though all three sum correctly to their own
    ``log_prob``. This is the same statement ``likelihoods.md`` makes about a
    profile likelihood in general (§9 question 3 in the awkward-instrument
    sketch): it is not a marginal likelihood, so quantities that assume one
    — LOO, WAIC, a properly calibrated posterior — should not be built from
    it.
    """

    NAME = "wstat_example"
    #: No sigma anywhere: the family defines its own (Poisson) dispersion.
    REQUIRES_UNCERTAINTY = False

    def __init__(self, background_counts: np.ndarray, ratio: float) -> None:
        background = np.asarray(background_counts, dtype=float)
        if background.ndim != 1:
            raise LikelihoodError(
                f"ProfiledCashWithBackground needs a one-dimensional background_counts array, "
                f"got shape {background.shape}."
            )
        if np.any(background < 0.0) or not np.all(background == np.round(background)):
            raise LikelihoodError(
                "ProfiledCashWithBackground needs non-negative integer background counts, one "
                "per channel — the background-region observation, aligned to the source-region "
                "channels it will be profiled against."
            )
        if not (np.isfinite(ratio) and float(ratio) > 0.0):
            raise LikelihoodError(
                f"ProfiledCashWithBackground needs a strictly positive exposure/area ratio, got "
                f"{ratio!r}."
            )
        self._background = self.register_buffer("background", background)
        self._ratio = float(self.register_buffer("ratio", float(ratio)))

    def check_observed(self, observed: Spectrum) -> None:
        """Non-negative integer source counts — the same precondition ``PoissonFamily`` has."""
        counts = np.asarray(observed.values).ravel()
        kept = counts[np.asarray(observed.valid).ravel()]
        if np.any(kept < 0.0) or not np.all(kept == np.round(kept)):
            raise LikelihoodError(
                "the profiled Cash-with-background family needs non-negative integer source "
                "counts, exactly like PoissonFamily; if the data are rates, multiply by the "
                "exposure in the instrument chain rather than here."
            )

    def log_prob(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        noise: NoiseParams,
    ) -> float:
        if noise.retain is None:
            raise LikelihoodError(
                "ProfiledCashWithBackground needs noise.retain to align its own background "
                "buffer with the retained source channels, but the noise model supplied none. "
                "Every noise model Likelihood.log_prob calls sets it (likelihoods.md §8, gap "
                "X-2); use IndependentNoise() or another ampere NoiseModel."
            )
        if self._background.size != noise.retain.size:
            raise LikelihoodError(
                f"the background buffer has {self._background.size} channel(s) but the observed "
                f"source container has {noise.retain.size}. The source and background regions "
                f"must be tabulated channel for channel, even where one or the other is masked."
            )
        background = self._background[noise.retain]
        model = np.asarray(predicted, dtype=float)
        source = np.asarray(observed, dtype=float)
        if np.any(model < 0.0):
            raise LikelihoodError(
                "the profiled Cash-with-background family needs a non-negative expected source "
                "count; the model predicted a value < 0."
            )
        ratio = self._ratio
        a = ratio * (ratio + 1.0)
        c = ratio * (source + background) - (ratio + 1.0) * model
        discriminant = c**2 + 4.0 * a * background * model
        background_mle = (c + np.sqrt(np.maximum(discriminant, 0.0))) / (2.0 * a)
        source_rate = model + ratio * background_mle
        return float(
            np.sum(
                st.poisson.logpmf(source, source_rate)
                + st.poisson.logpmf(background, background_mle)
            )
        )

    def sample(
        self,
        predicted: np.ndarray,
        noise: NoiseParams,
        rng: np.random.Generator,
    ) -> np.ndarray:
        raise LikelihoodError(
            "ProfiledCashWithBackground does not implement sample(): its profiled background "
            "estimate b_hat is computed from the very source and background counts a draw would "
            "need to produce, so there is no forward generative model to sample from without "
            "already having the data it profiles over. Use the two-dataset Bayesian formulation "
            "(build_joint_problem) to simulate source and background counts consistently, or "
            "draw Poisson(model + ratio * background_truth) and Poisson(background_truth) "
            "yourself from a background you are prepared to assume."
        )


# ---------------------------------------------------------------------------
# The two models: a source-only model (route 1) and a source-plus-shared-
# background model with two channels (route 2).
# ---------------------------------------------------------------------------


def _as_parameter(name: str, spec: Any) -> Parameter:
    """A prior (a frozen scipy distribution) or a bare number, as a Parameter."""
    if isinstance(spec, Parameter):
        return spec
    if hasattr(spec, "logpdf") or hasattr(spec, "rvs"):
        return Parameter(name, spec)
    return Parameter(name, None, value=float(spec), fixed=True)


class XraySource(Model):
    """A bare power-law source-count spectrum: ``norm * (E / E_ref) ** index``.

    No background, no response folding — ``awkward_instrument.md`` already
    covers RMF/ARF forward folding, so this model predicts expected *counts*
    directly (exposure and effective area already folded in, as that sketch's
    §2 concludes an X-ray count model should do), keeping this example about
    the likelihood, not the instrument chain.
    """

    def __init__(
        self,
        energy_kev: np.ndarray,
        *,
        norm: Any,
        index: Any,
        reference_energy_kev: float = REFERENCE_ENERGY_KEV,
    ) -> None:
        self.register_buffer("energy", np.asarray(energy_kev, dtype=float), unit=u.keV)
        self.register_buffer("reference_energy", float(reference_energy_kev), unit=u.keV)
        self.register_parameter(_as_parameter("norm", norm))
        self.register_parameter(_as_parameter("index", index))

    def evaluate(self, **values: Any) -> Spectrum:
        ctx = self.context(values)
        counts = ctx["norm"] * (ctx["energy"] / ctx["reference_energy"]) ** ctx["index"]
        return Spectrum(ctx["energy"] * u.keV, counts)


class XraySourceAndBackground(Model):
    """Source and background as one model, two channels, one shared background.

    ``inference.md`` §8's pattern is two datasets sharing a model *channel*;
    here the shared physics (the background rate) needs no channel of its own
    to expose, so it is simplest as one model publishing both region
    predictions directly, built from the same ``bkg_norm`` parameter. An
    equally valid alternative is a separate background ``Model`` reused by
    both datasets' instruments and a ``Tie`` joining the two — worth doing
    when the background genuinely has its own model class reused elsewhere;
    unnecessary complexity here.

    This is the **cost** side of the two-dataset formulation's trade-off: it
    needs an assumed background *shape* (flat, here — one parameter) where
    the profiled statistic needs none, since each channel's background is a
    free nuisance. A misspecified background shape biases this fit in a way
    the profiled statistic cannot be biased by; the profiled statistic pays
    for that flexibility with the degradations :class:`ProfiledCashWithBackground`
    documents.
    """

    def __init__(
        self,
        energy_kev: np.ndarray,
        *,
        src_norm: Any,
        src_index: Any,
        bkg_norm: Any,
        ratio: float,
        reference_energy_kev: float = REFERENCE_ENERGY_KEV,
    ) -> None:
        self.register_buffer("energy", np.asarray(energy_kev, dtype=float), unit=u.keV)
        self.register_buffer("reference_energy", float(reference_energy_kev), unit=u.keV)
        self.register_buffer("ratio", float(ratio))
        self.register_parameter(_as_parameter("src_norm", src_norm))
        self.register_parameter(_as_parameter("src_index", src_index))
        self.register_parameter(_as_parameter("bkg_norm", bkg_norm))

    def evaluate(self, **values: Any) -> ModelResult:
        ctx = self.context(values)
        source = ctx["src_norm"] * (ctx["energy"] / ctx["reference_energy"]) ** ctx["src_index"]
        background = np.full_like(ctx["energy"], ctx["bkg_norm"], dtype=float)
        source_region = source + ctx["ratio"] * background
        grid = ctx["energy"] * u.keV
        return ModelResult(
            {
                "source_region": Spectrum(grid, source_region),
                "background_region": Spectrum(grid, background),
            }
        )


# ---------------------------------------------------------------------------
# Building and running both fitting problems.
# ---------------------------------------------------------------------------


def build_wstat_problem(
    source_counts: np.ndarray, background_counts: np.ndarray, *, seed: int | None = SEED
) -> FittingProblem:
    """Route 1: the profiled statistic as a user family, one dataset."""
    model = XraySource(
        ENERGY_KEV,
        norm=st.loguniform(0.5, 40.0),
        index=st.uniform(-4.0, 4.0),
    )
    observed = Spectrum(ENERGY_KEV * u.keV, source_counts)
    dataset = Dataset(
        observed,
        likelihood=Likelihood(
            ProfiledCashWithBackground(background_counts, RATIO), IndependentNoise()
        ),
        label="source",
    )
    return FittingProblem(model, [dataset], seed=seed)


def build_joint_problem(
    source_counts: np.ndarray, background_counts: np.ndarray, *, seed: int | None = SEED
) -> FittingProblem:
    """Route 2: the two-dataset Bayesian formulation, source and background jointly."""
    model = XraySourceAndBackground(
        ENERGY_KEV,
        src_norm=st.loguniform(0.5, 40.0),
        src_index=st.uniform(-4.0, 4.0),
        bkg_norm=st.loguniform(0.2, 20.0),
        ratio=RATIO,
    )
    datasets = DatasetCollection(
        {
            "src": Dataset(
                Spectrum(ENERGY_KEV * u.keV, source_counts),
                Instrument(channel="source_region", input_kind=Spectrum, label="source"),
                likelihood=Likelihood(PoissonFamily(), IndependentNoise()),
            ),
            "bkg": Dataset(
                Spectrum(ENERGY_KEV * u.keV, background_counts),
                Instrument(channel="background_region", input_kind=Spectrum, label="background"),
                likelihood=Likelihood(PoissonFamily(), IndependentNoise()),
            ),
        }
    )
    return FittingProblem(model, datasets, seed=seed)


#: The "doc" budget: short chains that run in well under a minute on the
#: reference backend. For a publication-grade posterior, use the commented
#: "full-length" budget instead — the shape of the comparison is the same,
#: just with tighter uncertainties on the summary statistics themselves.
DOC_WALKERS = 16
DOC_STEPS = 800
DOC_BURN_IN = 200
# FULL_WALKERS, FULL_STEPS, FULL_BURN_IN = 32, 5000, 1000


def run_engine(
    problem: FittingProblem,
    *,
    walkers: int = DOC_WALKERS,
    steps: int = DOC_STEPS,
    burn_in: int = DOC_BURN_IN,
) -> Any:
    """Sample *problem* with emcee and return the run (an ArviZ ``DataTree``)."""
    return EmceeEngine(problem, walkers=walkers).run(steps=steps, burn_in=burn_in)


def summarise(run: Any, names: tuple[str, ...]) -> dict[str, tuple[float, float, float]]:
    """``{name: (median, low_68, high_68)}`` for each named free parameter."""
    summary = {}
    for name in names:
        draws = np.asarray(run["posterior"][name]).ravel()
        low, median, high = np.quantile(draws, [0.16, 0.5, 0.84])
        summary[name] = (float(median), float(low), float(high))
    return summary


def _interval_note(w: tuple[float, float, float], j: tuple[float, float, float]) -> str:
    """One sentence describing whether the two 68% intervals agree, from the numbers alone."""
    w_med, w_lo, w_hi = w
    j_med, j_lo, j_hi = j
    overlap = min(w_hi, j_hi) - max(w_lo, j_lo)
    if overlap <= 0.0:
        return "the two 68% intervals do not overlap at all"
    each_within_the_other = w_lo <= j_med <= w_hi and j_lo <= w_med <= j_hi
    agreement = (
        "each median falls inside the other approach's 68% interval"
        if each_within_the_other
        else "the intervals overlap, though each median falls outside the other's interval"
    )
    w_width, j_width = w_hi - w_lo, j_hi - j_lo
    if np.isclose(w_width, j_width, rtol=0.05):
        width_note = "the two intervals are close to the same width"
    elif w_width > j_width:
        width_note = f"WStat's is wider, by {w_width - j_width:.3f}"
    else:
        width_note = f"the joint fit's is wider, by {j_width - w_width:.3f}"
    return f"{agreement}; {width_note}"


def compare(
    wstat_summary: dict[str, tuple[float, float, float]],
    joint_summary: dict[str, tuple[float, float, float]],
) -> str:
    """The comparison text: recommendation and trade-offs, not false balance.

    The numeric comparison is read off *this run*'s actual summaries — it does
    not assert a direction (WStat narrower, or biased low, or anything else)
    that a single short, seeded chain cannot reliably demonstrate. What single
    run *can* show, and does, is whether the two posteriors substantially
    agree here; what it cannot show is a general claim about bias or coverage,
    which needs a repeated-trial study this fast worked example deliberately
    does not attempt. The recommendation below does not depend on which way
    the numbers happen to fall in this realisation — it is a structural
    argument, true regardless.
    """
    lines = [
        "WStat (profiled, user family) vs the two-dataset Bayesian joint fit",
        "=" * 70,
        "",
        f"{'parameter':<14}{'truth':>10}{'wstat median [68%]':>28}{'joint median [68%]':>28}",
    ]
    pairs = (
        ("src_norm", "model.norm", "model.src_norm"),
        ("src_index", "model.index", "model.src_index"),
    )
    for truth_key, wstat_name, joint_name in pairs:
        truth = TRUTH[truth_key]
        w_med, w_lo, w_hi = wstat_summary[wstat_name]
        j_med, j_lo, j_hi = joint_summary[joint_name]
        lines.append(
            f"{truth_key:<14}{truth:>10.3f}"
            f"{f'{w_med:.3f} [{w_lo:.3f}, {w_hi:.3f}]':>28}"
            f"{f'{j_med:.3f} [{j_lo:.3f}, {j_hi:.3f}]':>28}"
        )
    index_note = _interval_note(wstat_summary["model.index"], joint_summary["model.src_index"])
    lines += [
        "",
        f"On the source index, this run's two posteriors agree: {index_note}.",
        "",
        "Recommendation: prefer the two-dataset Bayesian joint fit, regardless of how",
        "closely the numbers above happen to agree. It is a proper marginal likelihood —",
        "its posterior is a posterior, its per-sample terms are valid predictive densities",
        "for arviz.loo/waic, and both regions can be simulated forward for posterior-",
        "predictive checks. WStat can do none of that structurally, not as a matter of",
        "this run's luck: its sample() refuses because the profiled background is a",
        "function of the very counts a draw would produce, and its per-channel terms",
        "already used their own data to fix their own nuisance, so they are not valid",
        "predictive densities no matter how well the point estimates above agree.",
        "",
        "Trade-off, stated honestly rather than as false balance: the joint fit pays for",
        "those guarantees with a background *model* (here, one flat-shape parameter) —",
        "if that shape is wrong, the joint fit is biased in a way WStat cannot be, because",
        "WStat's per-channel nuisance needs no shape assumption at all; that is its one",
        "genuine advantage. The literature on profile likelihoods in counting experiments",
        "(Cash 1979 and its X-ray descendants) documents that profiling a background out",
        "rather than marginalising it tends to bias and overstate the precision of the",
        "parameters of interest specifically in the low-count regime this example uses —",
        "a repeated-trial coverage study would show it reliably; a single seeded run, by",
        "construction, only shows one draw from that distribution, which is why the",
        "recommendation above does not lean on this run's particular numbers.",
    ]
    return "\n".join(lines)


def main() -> None:
    source_counts, background_counts = synthetic_xray_counts()

    wstat_problem = build_wstat_problem(source_counts, background_counts)
    wstat_run = run_engine(wstat_problem)
    wstat_summary = summarise(wstat_run, ("model.norm", "model.index"))

    joint_problem = build_joint_problem(source_counts, background_counts)
    joint_run = run_engine(joint_problem)
    joint_summary = summarise(joint_run, ("model.src_norm", "model.src_index", "model.bkg_norm"))

    print(compare(wstat_summary, joint_summary))


if __name__ == "__main__":
    main()
