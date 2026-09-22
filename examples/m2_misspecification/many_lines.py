"""W5.8's extension: one band of lines, one smooth error, three kernels.

M2's four scenarios each inject **one** scale of deviation, and that is why one
stationary length scale copes with all four: a Matérn-3/2 with a single
``length_scale`` has one scale to spend, and each of those scenarios gives it
one scale to spend it on. The plan's Phase-5 "non-stationary flexible
likelihood" bullet asks the question that follows — what happens when the
deviation has *two*, and only one of them is in one part of the band? — and
names the two answers to compare: W5.7's warped Matérn, which bends the
coordinate so that one length scale covers two, and W4.5's ``Sum`` of two
kernels, which simply adds a second.

:data:`~examples.m2_misspecification.generators.MANY_LINES` is that data: a
forest of five narrow lines between 0.860 and 0.870 µm, and a smooth continuum
error across the whole range. This module fits it four ways — the standard
likelihood, and the flexible likelihood with each of the three kernels — and
reports bias, calibration and localisation for each.

**Unlike** W4.5's fringing comparison, this one is *not* only informational.
The prediction is sharp enough to pin, and ``tests/m2/test_many_lines.py`` pins
it as a margin rather than as a number:

* the stationary Matérn is **not enough** here. It stays far better than the
  standard likelihood, but it leaves the truth further from its median than
  M2's own :data:`~examples.m2_misspecification.study.FLEXIBLE_MAX_BIAS_WIDTHS`
  threshold allows — which is the first time in this study a flexible
  likelihood has failed an M2 assertion, and it is failing it for the reason
  the plan predicted;
* the warped Matérn and the ``Sum`` both recover it, each inside that same
  threshold, each covering the truth on all four parameters;
* all three localise the deviation **inside the line band**, so the improvement
  is not the GP having been given somewhere else to put the residual.

The sparsity guard
------------------
The same freedom that makes the ``Sum`` arm possible is the freedom to add a
component the data do not need. :func:`shrinkage` is the demonstration of the
guard, and it is set up where the guard is the only thing that can decide the
answer: two **nearly degenerate** Matérn terms (:data:`DEGENERATE_SCALES`)
fitted to :data:`SMOOTH_ONLY`, a spectrum whose deviation has exactly **one**
smooth component. The likelihood pins the two terms' *total* and says almost
nothing about how it is divided between them, so what divides it is the prior.

Under a flat prior on both amplitudes the fit spreads itself across both terms;
under :func:`~ampere.core.shrinkage_horseshoe` it does not. The statistic is
the smaller amplitude over the larger — 1 is "split evenly between two
components the truth does not have two of", 0 is "chose one" — and it is
reported both as a posterior median and as the posterior **mass** below a tenth,
because a mass is far less sensitive to how well an ensemble sampler explored a
hierarchical posterior than a quantile of it is.

Why a degenerate pair rather than a spurious term at some far-off length scale:
a term the data can positively *exclude* is switched off by the likelihood, not
by the prior, and measuring it would be measuring the likelihood. Both priors
switch such a term off and they agree to three digits when they do — which is
correct behaviour and no demonstration at all.

Run it as::

    python -m examples.m2_misspecification.many_lines
    python -m examples.m2_misspecification.many_lines --size 400 --shrinkage
    python -m examples.m2_misspecification.many_lines --figures /tmp/m2

Nothing here writes a file (ground rule 7): the tables go to stdout.
"""

from __future__ import annotations

import argparse
import dataclasses
from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.core import (
    Dataset,
    FittingProblem,
    GaussianFamily,
    Kernel,
    Likelihood,
    shrinkage_horseshoe,
    with_shrinkage,
)

from .generators import MANY_LINES, LineForest, SyntheticSpectrum, generate
from .study import (
    DATASET_LABEL,
    GP_AMPLITUDE_SCALE,
    MANY_LINES_EMCEE,
    MANY_LINES_KERNELS,
    PHYSICAL_NAMES,
    SEED,
    Diagnosis,
    EmceeBudget,
    Summary,
    _backend_module,
    build_kernel,
    build_likelihood,
    diagnose,
    model_for,
    run,
    summarise,
)

__all__ = [
    "ARMS",
    "DEGENERATE_LABELS",
    "DEGENERATE_SCALES",
    "HORSESHOE_GLOBAL_SCALE",
    "SHRINKAGE_EMCEE",
    "SMOOTH_ONLY",
    "SPARSE_FRACTION",
    "Comparison",
    "Shrinkage",
    "build_arm",
    "compare",
    "degenerate_pair",
    "horseshoe_kernel",
    "main",
    "shrinkage",
]

#: The four fits the comparison runs, in the order its table prints them. The
#: first is the standard likelihood; the other three are the flexible
#: likelihood with :data:`~examples.m2_misspecification.study.MANY_LINES_KERNELS`.
ARMS: tuple[str, ...] = ("standard", *MANY_LINES_KERNELS)

#: The global scale of :func:`shrinkage`'s horseshoe, Jy. It is the prior guess
#: at "how big is a noise component that is really there", and M2's own fitted
#: GP amplitudes are what it is read off: a few thousandths of a Jy in the
#: control scenario, a few hundredths where there is 7 % fringing to absorb.
#: Not the study's :data:`~examples.m2_misspecification.study.GP_AMPLITUDE_SCALE`
#: of 0.3, which is a prior on one component's amplitude and far too permissive
#: to be a *global* scale — a horseshoe whose global scale is ten times the
#: components it is shrinking is a horseshoe that shrinks nothing.
HORSESHOE_GLOBAL_SCALE = 0.05

#: The one-component truth :func:`shrinkage` fits: W5.8's scenario with the
#: line forest switched off, so the only deviation left is the smooth continuum
#: error. A ``Sum`` of a broad and a narrow term fitted to it has one component
#: to find and one that is spurious, which is exactly the situation the
#: sparsity prior exists for.
SMOOTH_ONLY: LineForest = dataclasses.replace(
    MANY_LINES,
    key="smooth_only",
    description="Smooth continuum error only (the forest switched off)",
    amplitude=0.0,
)

#: The two length-scale prior scales of :func:`degenerate_pair`, micron. Close
#: enough together that a 200-point spectrum cannot say which of the two
#: absorbed the smooth error — which is the situation a sparsity prior exists
#: for — and far enough apart that the two terms are not exchangeable, so the
#: posterior has one mode rather than two mirror images of one.
DEGENERATE_SCALES: tuple[float, float] = (0.010, 0.008)

#: The labels those two terms carry into the posterior.
DEGENERATE_LABELS: tuple[str, str] = ("first", "second")

#: "Chose one component" — the ratio below which the smaller amplitude is
#: counted as switched off, for the posterior-mass statistic.
SPARSE_FRACTION = 0.10

#: The shrinkage demonstration's budget. Longer than
#: :data:`~examples.m2_misspecification.study.MANY_LINES_EMCEE`, and it has to
#: be: the horseshoe's three extra levels give the posterior a funnel, which is
#: the geometry an ensemble sampler explores worst, and at a shorter budget the
#: measured factors moved by more than the effect being measured. About three
#: and a half minutes for the pair of fits.
SHRINKAGE_EMCEE = EmceeBudget(walkers=32, steps=2_000, burn_in=1_000)


def degenerate_pair(backend: str = "reference") -> Kernel:
    """Two nearly degenerate Matérn-3/2 terms, under a **flat** amplitude prior.

    The contrast fit, and the reason its amplitude prior is uniform rather than
    the study's half-normal: a half-normal is *itself* a shrinkage prior — the
    study says so, and it is why M2's control scenario works — so comparing a
    horseshoe against one would be comparing two shrinkage priors. A flat prior
    is the honest "no opinion", and it is what a reader who had not thought
    about sparsity would write.

    The ceiling is the study's own
    :data:`~examples.m2_misspecification.study.GP_AMPLITUDE_SCALE`: wide enough
    that the posterior never reaches it, so the prior is flat where the answer
    is.
    """
    module = _backend_module(backend)
    return module.Sum(
        *[
            module.Matern32(
                st.uniform(0.0, GP_AMPLITUDE_SCALE),
                st.halfnorm(scale=scale),
                amplitude_unit=u.Jy,
                length_scale_unit=u.micron,
            )
            for scale in DEGENERATE_SCALES
        ],
        labels=DEGENERATE_LABELS,
    )


def horseshoe_kernel(backend: str = "reference", *, tail: str = "regularised") -> Kernel:
    """:func:`degenerate_pair` with a regularised horseshoe over its amplitudes.

    Same two terms, same two length-scale priors, same flat ceiling underneath
    — :func:`~ampere.core.with_shrinkage` replaces the two amplitudes'
    declarations and adds the two scale levels, and changes nothing else. That
    is what makes the contrast a contrast.
    """
    return with_shrinkage(
        degenerate_pair(backend),
        shrinkage_horseshoe(
            tuple(f"{label}.amplitude" for label in DEGENERATE_LABELS),
            global_scale=HORSESHOE_GLOBAL_SCALE,
            tail=tail,
            unit=u.Jy,
        ),
    )


def _problem(data: SyntheticSpectrum, kernel: Kernel, *, backend: str, seed: int) -> FittingProblem:
    module = _backend_module(backend)
    likelihood = Likelihood(
        GaussianFamily(), module.GaussianProcessNoise(kernel, module.QuasisepGP())
    )
    dataset = Dataset(data.container(), likelihood=likelihood, label=DATASET_LABEL)
    return FittingProblem(model_for(backend, data.wavelength), [dataset], seed=seed)


def build_arm(
    data: SyntheticSpectrum,
    arm: str,
    *,
    backend: str = "reference",
    seed: int = SEED,
) -> FittingProblem:
    """One of :data:`ARMS`, as a fitting problem over *data*.

    ``"standard"`` is the ordinary chi-square likelihood; the other three are
    the flexible likelihood over the kernel of that name, composed entirely
    from *backend*'s own parts (W2.13).
    """
    if arm == "standard":
        likelihood = build_likelihood("standard", backend=backend)
        dataset = Dataset(data.container(), likelihood=likelihood, label=DATASET_LABEL)
        return FittingProblem(model_for(backend, data.wavelength), [dataset], seed=seed)
    if arm not in MANY_LINES_KERNELS:
        raise ValueError(f"unknown arm {arm!r}; the four are {list(ARMS)}.")
    return _problem(data, build_kernel(backend, arm), backend=backend, seed=seed)


@dataclasses.dataclass(frozen=True)
class Comparison:
    """Four fits of one spectrum, reduced to what the extension asserts on.

    Attributes
    ----------
    scenario, size
        Which spectrum was fitted, and at how many points.
    summaries
        ``{arm: {parameter name: Summary}}`` for each of :data:`ARMS`.
    diagnoses
        ``{arm: Diagnosis}`` — residual whiteness for the standard arm, GP
        localisation for the three flexible ones (``study.diagnose``'s scoping,
        ``diagnostics.md`` §3.1 and §4).
    band
        The band the line forest was confined to, so a reader of the table can
        see whether the localisation peak landed in it.
    entries
        ``{arm: entry}``, each the ``run``/``problem``/``data`` mapping the
        study passes around, carrying the derived ``residuals`` or
        ``gp_localisation`` group its diagnostic needed. Kept because
        :func:`~examples.m2_misspecification.figures.figure_many_lines` draws
        from the conditioned GP means in them, and re-deriving those from a
        second fit would make the figure and the table two experiments.
    data
        The spectrum every arm was fitted to.
    """

    scenario: str
    size: int
    summaries: dict[str, dict[str, Summary]]
    diagnoses: dict[str, Diagnosis]
    band: tuple[float, float]
    entries: dict[str, dict[str, Any]] = dataclasses.field(default_factory=dict, repr=False)
    data: SyntheticSpectrum | None = dataclasses.field(default=None, repr=False)

    @property
    def worst_bias(self) -> dict[str, float]:
        """The largest ``|median - truth|`` in posterior widths, per arm."""
        return {
            arm: max(
                float(summary.bias_in_widths or 0.0)
                for summary in per_arm.values()
                if summary.truth is not None
            )
            for arm, per_arm in self.summaries.items()
        }

    @property
    def coverage(self) -> dict[str, int]:
        """How many of the four physical parameters' 68 % intervals cover the truth."""
        return {
            arm: sum(
                1
                for summary in per_arm.values()
                if summary.truth is not None and summary.covers_truth
            )
            for arm, per_arm in self.summaries.items()
        }

    @property
    def localised_in_band(self) -> dict[str, bool]:
        """Whether each flexible arm's GP-localisation peak fell inside the band."""
        lower, upper = self.band
        return {
            arm: bool(
                diagnosis.localisation_peak is not None
                and lower <= diagnosis.localisation_peak <= upper
            )
            for arm, diagnosis in self.diagnoses.items()
            if diagnosis.localisation_peak is not None
        }

    def table(self) -> str:
        """The comparison as plain text, one row per physical parameter."""
        header = (
            f"M2 {self.scenario!r}, {self.size} points — bias in posterior widths "
            f"(68 % coverage in brackets)\n"
            f"{'parameter':<14}" + "".join(f"{arm:>16}" for arm in ARMS)
        )
        rows = [header, "-" * (14 + 16 * len(ARMS))]
        for name in PHYSICAL_NAMES:
            cells = []
            for arm in ARMS:
                summary = self.summaries[arm][name]
                offset = summary.bias_in_widths
                covered = "y" if summary.covers_truth else "n"
                cells.append(f"{offset:.2f} [{covered}]" if offset is not None else "-")
            rows.append(f"{name:<14}" + "".join(f"{cell:>16}" for cell in cells))
        rows.append("-" * (14 + 16 * len(ARMS)))
        worst = self.worst_bias
        covered = self.coverage
        rows.append(f"{'worst bias':<14}" + "".join(f"{worst[arm]:>16.2f}" for arm in ARMS))
        rows.append(f"{'covered / 4':<14}" + "".join(f"{covered[arm]:>16d}" for arm in ARMS))
        peaks = []
        for arm in ARMS:
            peak = self.diagnoses[arm].localisation_peak
            peaks.append("-" if peak is None else f"{peak:.5f}")
        rows.append(f"{'GP peak (um)':<14}" + "".join(f"{peak:>16}" for peak in peaks))
        rows.append(f"the line band is {self.band[0]:.3f} to {self.band[1]:.3f} um")
        return "\n".join(rows)


def compare(
    *,
    size: int = 200,
    backend: str = "reference",
    budget: EmceeBudget | None = None,
    seed: int = SEED,
    progress: bool = False,
    thin: int = 20,
) -> Comparison:
    """Fit the ``many_lines`` spectrum four ways and reduce them to a :class:`Comparison`.

    One spectrum, generated once and shared, so the four fits differ in their
    likelihood and in nothing else — ``study.run_study``'s discipline, for the
    same reason.
    """
    data = generate(MANY_LINES, size=size)
    chosen = budget or MANY_LINES_EMCEE
    summaries: dict[str, dict[str, Summary]] = {}
    diagnoses: dict[str, Diagnosis] = {}
    entries: dict[str, dict[str, Any]] = {}
    for arm in ARMS:
        problem = build_arm(data, arm, backend=backend, seed=seed)
        entry: dict[str, Any] = {
            "run": run(problem, chosen, progress=progress),
            "problem": problem,
            "data": data,
            "likelihood": "standard" if arm == "standard" else "flexible",
        }
        summaries[arm] = summarise(entry["run"], names=PHYSICAL_NAMES)
        diagnoses[arm] = diagnose(entry, thin=thin)
        entries[arm] = entry
    return Comparison(
        scenario=MANY_LINES.key,
        size=size,
        summaries=summaries,
        diagnoses=diagnoses,
        band=MANY_LINES.band,
        entries=entries,
        data=data,
    )


@dataclasses.dataclass(frozen=True)
class Shrinkage:
    """What a sparsity prior did to a redundant noise component.

    Attributes
    ----------
    size
        Points in the fitted spectrum.
    amplitudes
        ``{"flat" | "horseshoe": {label: Summary}}``, one entry per term of
        :func:`degenerate_pair`.
    ratios
        ``{"flat" | "horseshoe": array}`` — the smaller amplitude over the
        larger one, **per posterior draw**. Per draw rather than a ratio of
        summaries because that is what a posterior mass can be computed from,
        and the mass is the statistic that survives an ensemble sampler's
        difficulty with a hierarchical geometry.

    Notes
    -----
    The truth has one smooth component, so one of the two terms is redundant
    whichever it turns out to be. Which one it is is not part of the claim —
    the two are nearly degenerate, and a prior that chose the *named* one would
    be a prior that had been told the answer — so every statistic here is over
    ``min`` and ``max`` rather than over ``first`` and ``second``.
    """

    size: int
    amplitudes: dict[str, dict[str, Summary]]
    ratios: dict[str, Any] = dataclasses.field(default_factory=dict, repr=False)

    @property
    def ratio(self) -> dict[str, float]:
        """The posterior median of ``min(a) / max(a)``, per prior."""
        return {prior: float(np.median(values)) for prior, values in self.ratios.items()}

    @property
    def sparse_mass(self) -> dict[str, float]:
        """The posterior mass with ``min(a) / max(a)`` below :data:`SPARSE_FRACTION`."""
        return {
            prior: float(np.mean(np.asarray(values) < SPARSE_FRACTION))
            for prior, values in self.ratios.items()
        }

    @property
    def largest(self) -> dict[str, float]:
        """The median of the *larger* amplitude, per prior: the component the truth has."""
        return {
            prior: max(summary.median for summary in pair.values())
            for prior, pair in self.amplitudes.items()
        }

    def table(self) -> str:
        """The demonstration as plain text."""
        ratio = self.ratio
        mass = self.sparse_mass
        largest = self.largest
        rows = [
            f"Two nearly degenerate noise terms on a one-component truth, {self.size} points",
            (
                f"{'prior':<12}{'larger (Jy)':>14}{'median min/max':>16}"
                f"{f'P(min/max < {SPARSE_FRACTION})':>20}"
            ),
            "-" * 62,
        ]
        for prior in ("flat", "horseshoe"):
            rows.append(
                f"{prior:<12}{largest[prior]:>14.5f}{ratio[prior]:>16.3f}{mass[prior]:>20.3f}"
            )
        rows.append("-" * 62)
        rows.append(
            f"the horseshoe divides the redundant component by "
            f"{ratio['flat'] / ratio['horseshoe']:.1f} and multiplies the sparse mass by "
            f"{mass['horseshoe'] / mass['flat']:.1f}"
        )
        return "\n".join(rows)


def shrinkage(
    *,
    size: int = 200,
    backend: str = "reference",
    budget: EmceeBudget | None = None,
    seed: int = SEED,
    progress: bool = False,
    tail: str = "regularised",
) -> Shrinkage:
    """Fit two degenerate terms to a one-component truth, with and without the horseshoe.

    Both fits are over :data:`SMOOTH_ONLY` — W5.8's scenario with the line
    forest switched off — so there is one smooth component to find and two
    terms that can both find it. The two kernels are the same two Matérns with
    the same two length-scale priors and the same flat ceiling; only the prior
    on their amplitudes differs.
    """
    data = generate(SMOOTH_ONLY, size=size)
    chosen = budget or SHRINKAGE_EMCEE
    kernels = {
        "flat": degenerate_pair(backend),
        "horseshoe": horseshoe_kernel(backend, tail=tail),
    }
    names = {label: f"{DATASET_LABEL}.likelihood.{label}.amplitude" for label in DEGENERATE_LABELS}
    amplitudes: dict[str, dict[str, Summary]] = {}
    ratios: dict[str, Any] = {}
    for prior, kernel in kernels.items():
        problem = _problem(data, kernel, backend=backend, seed=seed)
        stored = run(problem, chosen, progress=progress)
        found = summarise(stored, names=list(names.values()))
        amplitudes[prior] = {label: found[name] for label, name in names.items()}
        drawn = np.vstack(
            [
                np.asarray(stored["posterior"][name].values, dtype=float).ravel()
                for name in names.values()
            ]
        )
        ratios[prior] = drawn.min(axis=0) / drawn.max(axis=0)
    return Shrinkage(size=size, amplitudes=amplitudes, ratios=ratios)


def main(argv: Any = None) -> int:
    """Print the comparison table, and the shrinkage table when asked for it."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--size", type=int, default=200)
    parser.add_argument("--backend", default="reference")
    parser.add_argument("--walkers", type=int, default=MANY_LINES_EMCEE.walkers)
    parser.add_argument("--steps", type=int, default=MANY_LINES_EMCEE.steps)
    parser.add_argument("--burn-in", type=int, default=MANY_LINES_EMCEE.burn_in)
    parser.add_argument("--progress", action="store_true")
    parser.add_argument(
        "--shrinkage",
        action="store_true",
        help="also run the horseshoe demonstration on a one-component truth",
    )
    parser.add_argument("--figures", default=None, help="write W5.8's figure into this directory")
    args = parser.parse_args(argv)
    budget = EmceeBudget(args.walkers, args.steps, args.burn_in)
    comparison = compare(
        size=args.size, backend=args.backend, budget=budget, progress=args.progress
    )
    print(comparison.table())
    if args.figures is not None:
        from . import figures

        written = figures.save_many_lines_figure(
            comparison.entries, comparison.data, args.figures, band=comparison.band
        )
        print(f"\nWrote {len(written)} figure to {args.figures}")
    if args.shrinkage:
        print()
        print(shrinkage(size=args.size, backend=args.backend, progress=args.progress).table())
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI
    raise SystemExit(main())
