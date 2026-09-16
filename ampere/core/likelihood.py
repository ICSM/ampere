"""Likelihood families and noise models: ``DEVELOPMENT_PLAN.md`` §4.4.

This is the contract in which ampere says *how well a prediction matches an
observation*. It has three separable pieces, and keeping them separable is the
whole point (``prior_art.md`` Tension 3: 3ML fuses instrument response and
likelihood into one opaque plugin, and buys extensibility at the cost of all
reuse; ampere factors them and pays for it with this spec):

- a :class:`LikelihoodFamily` — the sampling distribution of a single datum
  given a prediction and its noise. Registered through one method,
  ``log_prob(predicted, observed, noise)`` (``prior_art.md`` lesson 3M3:
  copy 3ML's minimalism, not its opacity), so a user adds a family without
  touching ampere;
- a :class:`NoiseModel` — what the noise *is*: independent per-sample
  uncertainties (:class:`IndependentNoise`) or a Gaussian process over the
  residuals (:class:`GaussianProcessNoise`), which is ampere's flexible,
  misspecification-robust likelihood;
- a :class:`GPSolver` — *how* the GP algebra is done. A swappable strategy:
  :class:`DenseGP` (exact, O(N³), the correctness anchor) and
  :class:`QuasisepGP` (exact, O(N) on ordered 1D data, over celerite2), plus
  the approximate strategies for the cases neither covers.

:class:`Likelihood` composes the three and is the object a ``Dataset`` (W1.7)
holds.

Two structural facts drive everything below.

**Float64, always.** Every array this module touches is cast to
:data:`~numpy.float64` (or ``complex128``). ``DEVELOPMENT_PLAN.md`` §7 and
``architecture.md`` §5 are explicit that GP Cholesky and quasiseparable solves
in float32 fail in ways that read as science bugs rather than numerical ones,
and that this is not a per-backend default to inherit.

**Analytic marginalisation is a property of the family/noise *combination*,
not of either alone.** The "add the GP covariance to the noise and integrate"
trick only closes for a Gaussian noise process. For Poisson, Student-t and the
rest, robustness needs a latent-GP formulation (``counts ~ Poisson(rate ·
exp(f))``, ``f ~ GP``) marginalised numerically — so a :class:`Likelihood`
declares its :class:`Marginalisation`, and refuses an engine that cannot
deliver it (``DEVELOPMENT_PLAN.md`` §4.4, "the Likelihood interface must
therefore declare whether it marginalises its noise process analytically or
introduces latent variables").

Nothing here imports torch, jax or any optional dependency
(``architecture.md`` §4 rule 1): it is numpy, scipy, ``astropy.units`` and
stdlib. celerite2 — a **base** dependency, not an extra (``architecture.md``
§2) — is imported lazily inside :class:`QuasisepGP`, so importing this module
still costs nothing beyond that list. The narrative spec is
``docs/design/contracts/likelihoods.md``, whose every example runs as a
doctest.
"""

from __future__ import annotations

import abc
import dataclasses
import enum
import functools
import math
from collections.abc import Mapping, Sequence
from typing import Any, ClassVar

import astropy.units as u
import numpy as np
import scipy.linalg
import scipy.special
import scipy.stats as st

from .exceptions import LikelihoodError
from .kernels import (
    DTYPE,
    NUMPY_OPS,
    SHO,
    ArrayOps,
    CeleriteRepresentation,
    Kernel,
    KernelSpec,
    Matern12,
    Matern32,
    Matern52,
    NumpyOps,
    Product,
    QuasiseparableTerm,
    RotationTerm,
    SpectralMixture,
    SquaredExponential,
    StationaryKernel,
    Sum,
    TermBuilder,
    _as_float64,
    _as_hyperparameter,
    _as_points,
    _check_finite,
    _positive,
    lookup_quasiseparable_term,
    quasiseparable_families,
    register_quasiseparable_term,
    registered_quasiseparable_terms,
    term_provenance_entries,
)
from .hsgp import (
    DEFAULT_BASIS_SIZE,
    DEFAULT_BOUNDARY_FACTOR,
    HilbertSpaceBasis,
    basis_matrix,
    basis_size,
    check_spectral_support,
    hilbert_basis,
    normalise_counts,
    spectral_values,
)
from .parameter import (
    Identity,
    Parameter,
    Parameterised,
    ParameterSet,
)
from .results_schema import FunctionSamples, Layout

__all__ = [
    "DTYPE",
    "NUMPY_OPS",
    "SHO",
    "ArrayOps",
    "CauchyFamily",
    "CeleriteRepresentation",
    "Censoring",
    "ChannelCoupling",
    "CholeskyCoupling",
    "ComplexGaussianFamily",
    "DenseGP",
    "GPConditional",
    "GPSolver",
    "GaussianFamily",
    "GaussianProcessNoise",
    "HilbertSpaceGP",
    "IndependentNoise",
    "InducingPointGP",
    "JointGaussianProcessNoise",
    "Kernel",
    "KernelSpec",
    "LatentDeclaration",
    "Likelihood",
    "LikelihoodFamily",
    "LimitKind",
    "Marginalisation",
    "Matern12",
    "Matern32",
    "Matern52",
    "NoiseModel",
    "NoiseParams",
    "NumpyOps",
    "PoissonFamily",
    "Product",
    "QuasisepGP",
    "QuasiseparableTerm",
    "RiceFamily",
    "RotationCoupling",
    "RotationTerm",
    "SpectralMixture",
    "SquaredExponential",
    "StationaryKernel",
    "StructuredGridGP",
    "StudentTFamily",
    "Sum",
    "TermBuilder",
    "VecchiaGP",
    "VonMisesFamily",
    "WindowedSparseGP",
    "family_named",
    "latent_parameter",
    "list_families",
    "lookup_quasiseparable_term",
    "quasiseparable_families",
    "register_family",
    "register_quasiseparable_term",
    "registered_quasiseparable_terms",
    "term_provenance_entries",
]

_LOG_2PI = math.log(2.0 * math.pi)


def _stacked_components(residual: np.ndarray) -> np.ndarray:
    """A complex residual as the ``(n, 2)`` real block the circular GP solves on.

    **W4.2.** The circular complex Gaussian's real and imaginary parts are
    independent real processes sharing one covariance ``K + diag(σ²)``, so the
    2N-dimensional real covariance is block diagonal with the *same* block
    twice. Writing the residual as two columns is what lets a solver exploit
    that: one factorisation, two solves, ``log|K + diag(σ²)|`` once and then
    doubled, rather than a dense 2N by 2N factorisation costing eight times as
    much and carrying a zero off-diagonal block it already knows is zero.

    A real residual passes through untouched, so every pre-W4.2 call is
    unchanged: a one-column right-hand side is the ``k = 1`` case of the same
    rule (:class:`GPSolver`).
    """
    values = np.asarray(residual)
    if values.dtype.kind != "c":
        return values
    return np.ascontiguousarray(np.column_stack([values.real, values.imag]), dtype=DTYPE)


def _components(residual: np.ndarray) -> int:
    """How many independent realisations a right-hand side carries: ``k``."""
    array = np.asarray(residual)
    return 1 if array.ndim == 1 else int(array.shape[1])


# ---------------------------------------------------------------------------
# Declarations: marginalisation and censoring
# ---------------------------------------------------------------------------


class Marginalisation(enum.Enum):
    """How a family/noise combination deals with its noise process.

    ``DEVELOPMENT_PLAN.md`` §4.4 requires this to be *declared*, because it
    decides which inference engines can run the problem at all.
    """

    #: The noise process is integrated out in closed form. Any engine works.
    ANALYTIC = "analytic"
    #: The noise process enters as latent variables that inference must sample
    #: or approximate. Effectively a modern-backend (HMC/VI) capability.
    LATENT = "latent"


class LimitKind(enum.IntEnum):
    """What a single observed sample asserts about the truth.

    An :class:`enum.IntEnum` so a per-sample classification is an ordinary
    integer array — which is exactly the aligned-auxiliary-array hook
    ``results_schema.md`` §8 left open for this contract (issue #11).
    """

    #: An ordinary measurement: the value is the datum.
    DETECTION = 0
    #: A non-detection: the truth is *below* the recorded value.
    UPPER_LIMIT = 1
    #: The truth is *above* the recorded value (a saturated pixel, say).
    LOWER_LIMIT = 2


def _readonly(array: np.ndarray) -> np.ndarray:
    array.flags.writeable = False
    return array


@dataclasses.dataclass(frozen=True)
class Censoring:
    """A per-sample statement of which observations are limits, not detections.

    Censoring is *not* masking. ``results_schema.md`` §8 draws the line and
    hands the other side to this contract: a masked sample carries exactly zero
    information, whereas a 3sigma non-detection says something quite definite about
    the source. Discarding it throws that away; recording it as a datum with a
    large error bar is simply the wrong likelihood.

    A ``Censoring`` is declared *alongside* an observed container rather than
    inside it, because it is a statement about how the observation constrains
    the model — a likelihood concern — and the same container may legitimately
    be analysed with and without it.

    Parameters
    ----------
    kinds
        Integer array of :class:`LimitKind` codes, aligned index-by-index with
        the observed container's ``values``.

    Examples
    --------
    >>> import numpy as np
    >>> import astropy.units as u
    >>> from ampere.core import PhotometricPoints
    >>> phot = PhotometricPoints(
    ...     ["WISE_W3", "WISE_W4"],
    ...     [12.1, 22.2] * u.um,
    ...     [0.4, 0.9] * u.Jy,
    ...     uncertainty=[0.05, 0.30] * u.Jy,
    ...     extra_coords={"limit_kind": np.array([0, 1])},
    ... )
    >>> censoring = Censoring.from_extra_coord(phot, "limit_kind")
    >>> censoring.n_censored, censoring.any_censored
    (1, True)
    """

    kinds: np.ndarray

    def __post_init__(self) -> None:
        raw = np.asarray(self.kinds)
        if raw.dtype.kind == "b" or raw.dtype.kind not in "iu":
            raise LikelihoodError(
                f"Censoring's kinds must be an integer array of LimitKind codes "
                f"({[k.name for k in LimitKind]}), got dtype {raw.dtype}. A boolean array is "
                f"ambiguous — say which kind of limit it means with "
                f"Censoring.upper_limits(flags) or Censoring.lower_limits(flags)."
            )
        allowed = {int(kind) for kind in LimitKind}
        unknown = sorted(set(np.unique(raw).tolist()) - allowed)
        if unknown:
            raise LikelihoodError(
                f"Censoring was given unknown limit codes {unknown}. The declared codes are "
                f"{ {kind.name: int(kind) for kind in LimitKind} }."
            )
        object.__setattr__(self, "kinds", _readonly(np.array(raw, dtype=np.int8, copy=True)))

    # -- construction --------------------------------------------------------

    @classmethod
    def from_extra_coord(cls, container: FunctionSamples, name: str = "limit_kind") -> Censoring:
        """Read a censoring declaration off a container's ``extra_coords``.

        This is the hook ``results_schema.md`` §8 reserved: per-sample aligned
        arrays are validated against the container's shape exactly as values
        and mask are, "which is precisely the machinery a per-point
        ``limit_kind`` array needs".
        """
        if name not in container.extra_coords:
            available = sorted(container.extra_coords)
            raise LikelihoodError(
                f"{type(container).__name__} has no extra coordinate {name!r} to read a censoring "
                f"declaration from; it carries {available}. Attach one at construction, e.g. "
                f"extra_coords={{{name!r}: np.array([0, 1, 0])}} with LimitKind codes."
            )
        return cls(np.asarray(container.extra_coords[name]).ravel())

    @classmethod
    def upper_limits(cls, flags: Any) -> Censoring:
        """Build a declaration from a boolean "is an upper limit" array."""
        return cls._from_flags(flags, LimitKind.UPPER_LIMIT)

    @classmethod
    def lower_limits(cls, flags: Any) -> Censoring:
        """Build a declaration from a boolean "is a lower limit" array."""
        return cls._from_flags(flags, LimitKind.LOWER_LIMIT)

    @classmethod
    def _from_flags(cls, flags: Any, kind: LimitKind) -> Censoring:
        booleans = np.asarray(flags)
        if booleans.dtype.kind != "b":
            raise LikelihoodError(
                f"Censoring.{kind.name.lower()}s() expects a boolean array of flags, got dtype "
                f"{booleans.dtype}."
            )
        codes = np.where(booleans, int(kind), int(LimitKind.DETECTION))
        return cls(codes.ravel())

    # -- introspection -------------------------------------------------------

    @property
    def n_samples(self) -> int:
        """Number of samples this declaration covers."""
        return int(self.kinds.size)

    @property
    def n_censored(self) -> int:
        """Number of samples that are limits rather than detections."""
        return int(np.count_nonzero(self.kinds != int(LimitKind.DETECTION)))

    @property
    def any_censored(self) -> bool:
        """Whether this declaration contains any limit at all."""
        return self.n_censored > 0

    def check_against(self, container: FunctionSamples) -> None:
        """Verify alignment with the container this declaration describes."""
        if self.n_samples != container.n_samples:
            raise LikelihoodError(
                f"Censoring covers {self.n_samples} samples but the "
                f"{type(container).__name__} it was given has {container.n_samples}. A censoring "
                f"declaration is aligned index-by-index with the observed container."
            )

    def __repr__(self) -> str:
        counts = {kind.name: int(np.count_nonzero(self.kinds == int(kind))) for kind in LimitKind}
        body = ", ".join(f"{name}={n}" for name, n in counts.items() if n)
        return f"Censoring({body})"


# ---------------------------------------------------------------------------
# GP solver strategies
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class GPConditional:
    """The GP conditioned on the residuals: mean and variance at each location.

    This is the object ``DEVELOPMENT_PLAN.md`` §4.8 and the diagnostics spec's
    family C consume — "the conditioned GP mean ... already localises where the
    model is deficient". The mean is **signed**, so it shows the direction of
    the local deficiency and not merely its size.

    **W4.2**: with a multi-realisation right-hand side (:class:`GPSolver`) the
    mean is ``(m, k)`` and the variance stays ``(m,)``, since a posterior
    variance does not depend on the data. :meth:`Likelihood.conditional`
    recombines the circular complex GP's two columns into a **complex** mean
    before returning, so a caller who conditioned on visibilities gets a signed
    deficiency in each component rather than two anonymous columns.
    """

    mean: np.ndarray
    variance: np.ndarray

    @property
    def standard_deviation(self) -> np.ndarray:
        """Pointwise 1sigma band on :attr:`mean`."""
        return np.sqrt(np.clip(self.variance, 0.0, None))


def _find_nested_product(
    kernel: Kernel, path: tuple[str, ...] = ()
) -> tuple[Product, tuple[str, ...]] | None:
    """The first :class:`Product` in *kernel*'s tree, and the label path to it.

    W5.2: a bare :class:`Product` is refused by name (:meth:`GPSolver.check_compatible`)
    because "products are not quasiseparable" is a sharper diagnosis than the
    generic ``not kernel.QUASISEPARABLE`` one -- but that check only fired when
    *kernel* itself was the ``Product``. A :class:`Sum` that merely *contains*
    one — ``Sum(Matern32(...), Product(...))`` — is not itself a ``Product``,
    so it fell through to the generic refusal, which names ``Sum`` rather than
    the term that is actually the problem. This walks the composite tree (via
    :attr:`Kernel.terms`, empty for a leaf) and returns the first ``Product``
    found, together with the dotted label path from the root, so the caller
    can name it directly. The empty path means *kernel* itself is the
    ``Product`` -- the case the original check already handled correctly.
    """
    if isinstance(kernel, Product):
        return kernel, path
    for label, child in kernel.terms:
        found = _find_nested_product(child, (*path, label))
        if found is not None:
            return found
    return None


class GPSolver(abc.ABC):
    """How the GP algebra is done. A strategy, chosen per problem.

    ``DEVELOPMENT_PLAN.md`` §4.4 makes the solver swappable behind the noise
    model precisely so the scaling story can change without the science code
    changing. Every strategy computes the same quantity — the Gaussian marginal
    log-likelihood of the residuals under ``K(θ) + diag(sigma²)`` — and
    :class:`DenseGP` is the definition of the right answer.

    A concrete strategy declares its own applicability through the class
    attributes below and enforces it in :meth:`check_compatible`, at
    composition time, loudly. It may branch on ``Axis.regular`` or
    ``Axis.log_regular`` for a fast path but must have a path when both are
    false (``results_schema.md`` §16).

    **The four capability flags reach here at W2.13** (ruled 2026-09-07,
    ``DEVELOPMENT_PLAN.md`` §2's "Realisation surface (W2.13)" row, fold-in 7;
    ``inference.md`` §10a), with the reference path's honest defaults. The
    consequence is deliberate and loud: a problem built from a native backend's
    models and steps but left with ``ampere.core.DenseGP`` — which factorises
    in scipy — now reports two backends and is refused at composition, because
    a nominally differentiable problem whose GP solve runs in numpy is a
    problem whose GP hyperparameters get no gradient at all. Pass the
    backend's own solver.

    **The right-hand side may carry several realisations (W4.2).** Every method
    below that takes a ``residual`` accepts either an ``(n,)`` array -- one
    realisation, which is every pre-W4.2 call -- or an ``(n, k)`` block of
    ``k`` realisations that are **independent of one another and share this one
    covariance**. The declared quantity is then the joint one:
    :meth:`log_marginal_likelihood` returns the sum of the ``k`` marginals, so
    the log-determinant enters ``k`` times and the factorisation happens once;
    :meth:`conditional_loo` returns one term per *sample*, summed over the
    ``k`` components, since the precision diagonal they share is what a
    leave-one-out conditional is built from; :meth:`condition` returns an
    ``(m, k)`` mean beside the single ``(m,)`` variance the components share;
    :meth:`latent_transform` maps ``(n, k)`` whitened draws to ``(n, k)``.

    That is not a convenience. It is what makes the **circular complex GP** of
    ``likelihoods.md`` §4 computable at the cost of one real solve rather than
    eight (see :class:`ComplexGaussianFamily`): the real 2N covariance of a
    circular complex Gaussian is ``K + diag(sigma^2)`` twice on the diagonal
    and zero off it, and a two-column right-hand side is that structure written
    arithmetically. A strategy that cannot take more than one column says so
    through :attr:`STACKED_RESIDUALS`, and
    :meth:`GaussianProcessNoise.check_compatible` refuses it by name on complex
    data rather than letting it silently score one component.
    """

    #: Neutral strategy name, for provenance and error messages.
    NAME: ClassVar[str] = ""
    #: Whether the strategy computes the exact marginal likelihood of the
    #: declared kernel, as opposed to an approximation of it.
    EXACT: ClassVar[bool] = True
    #: Whether it needs a single, ordered, one-dimensional coordinate axis.
    REQUIRES_ORDERED_1D: ClassVar[bool] = False
    #: Whether it needs the kernel to have a quasiseparable representation.
    REQUIRES_QUASISEPARABLE: ClassVar[bool] = False
    #: Whether an implementation exists in the reference (numpy) path.
    IMPLEMENTED: ClassVar[bool] = False
    #: Whether the strategy accepts an ``(n, k)`` right-hand side of ``k``
    #: realisations sharing one covariance (**W4.2**; see the class docstring).
    #: ``False`` by default, which is the honest answer for a strategy written
    #: before the rule existed and for every declared slot: a solver that
    #: flattened a two-column residual would score the two components against a
    #: single covariance of the wrong size, and a solver that took only the
    #: first column would silently drop the imaginary part of every visibility.
    #: :meth:`GaussianProcessNoise.check_compatible` reads it and refuses by
    #: name, so the flag is a declaration a user-written solver opts into once
    #: its own algebra handles the extra columns.
    STACKED_RESIDUALS: ClassVar[bool] = False

    #: Whether a gradient can be taken through this solver's linear algebra.
    DIFFERENTIABLE: ClassVar[bool] = False
    #: Whether it solves a batch of parameter vectors in one call.
    BATCHABLE: ClassVar[bool] = False
    #: Device its arrays live on. Never auto-detected (``architecture.md`` §5).
    DEVICE: ClassVar[str] = "cpu"
    #: Which rung of the capability ladder supplies it (W2.12's fourth flag).
    BACKEND: ClassVar[str] = "reference"

    def provenance_config(self) -> Mapping[str, Any]:
        """Backend-specific solver configuration, for the run's attrs only.

        **W2.13, fold-in 10** (ruled 2026-09-07), answering the question W2.4
        carried out of its slice 1: where does a torch solver record its dtype
        and device — and, later, a deliberate float32 opt-out — given that
        :meth:`Likelihood.to_spec` records a *dataclass* solver's fields and
        the cross-backend conformance rows compare component spec hashes?

        The answer is: not in the spec. A solver's ``jitter`` is a
        **declaration** — it changes the model, so it belongs in the spec and
        in the hash — while its dtype, its device and its precision policy are
        **configuration**: they change how the same declared model is computed,
        two backends legitimately differ on them, and folding them into the
        spec hash would break the one promise ``results.md`` §14 makes about
        two backends implementing one declaration. So they are class-level
        policy on the backend's solver (never dataclass fields), returned from
        here, and recorded by ``ampere.results.provenance.provenance_attrs``
        under ``ampere_solver_config`` — visible, never hashed.

        Returns
        -------
        Mapping[str, Any]
            JSON-normalisable values. Empty by default: the reference solvers
            have no configuration beyond what they declare.
        """
        return {}

    def latent_size(self, kernel: Kernel, n_samples: int) -> int:
        """How many whitened variables :meth:`latent_transform` takes (**W5.4**).

        ``inference.md`` §17.4 fixes the latent block at composition time, and
        until W5.4 fixed it at *one value per retained sample* — which was a
        property of the two exact solvers rather than of the contract. An
        approximate solver's whitening need not be square: a reduced-rank
        representation has a factor of shape ``(N, m)``, so ``m`` whitened
        variables produce ``N`` correlated ones and the latent block is ``m``
        (``docs/design/horizon_notes.md`` §2, question (b), settled at W5.4).

        The default is ``n_samples``, which is what every exact solver
        answers. An override must depend on the **declaration** alone —
        ``Likelihood.latent_declaration`` asks before any container is in
        hand — and must agree with what :meth:`latent_transform` accepts and
        with what ``LikelihoodFamily.sample`` draws, because those two are the
        halves ``simulate(observe=True)`` and the latent path have to share.
        """
        del kernel
        return int(n_samples)

    def check_compatible(self, kernel: Kernel, observed: FunctionSamples) -> None:
        """Composition-time check that this strategy can run this problem.

        **W4.5 widens two of the rules to the kernel's axis selection.** The
        single-unit rule is now applied by :meth:`Kernel.check_axes`, per leaf
        kernel, over the axes that leaf *selects* — a three-axis container is
        refused for a bare isotropic kernel exactly as it was, and accepted for
        ``Matern32(axes=("u", "v"))``. The ordered-1D rule counts the axes the
        kernel tree actually uses rather than the container's, so a
        quasiseparable solve over the spectral axis of a three-axis container
        is expressible; and a :class:`~ampere.core.kernels.Product` — anywhere
        in the kernel's tree, not only at the root, so a
        ``Sum(Matern32(...), Product(...))`` names the ``Product`` term
        rather than the enclosing ``Sum`` (W5.2) — is refused on a
        quasiseparable solver by name, before the generic ``QUASISEPARABLE``
        refusal, because "products are not quasiseparable" is a sharper
        diagnosis than "this kernel is not".
        """
        kind = type(observed).__name__
        # Declarative incompatibilities first: they are permanent facts about
        # the choice, whereas "not implemented yet" is temporary, and a user
        # who paired QuasisepGP with a non-quasiseparable kernel needs to hear
        # about the kernel rather than about Phase 2's schedule.
        if observed.LAYOUT is not Layout.POINTS:
            raise LikelihoodError(
                f"{self.NAME} was given a {kind}, whose layout is {observed.LAYOUT.value}. The "
                f"v1 GP solvers work on point-set containers (Spectrum, TimeSeries, "
                f"PhotometricPoints, VisibilitySet); gridded 2D+ data are the subject of the "
                f"SVGP / SKI / Vecchia strategy slots (DEVELOPMENT_PLAN.md §4.4, Phase 5)."
            )
        kernel.check_axes(observed, owner=self.NAME)
        if self.REQUIRES_QUASISEPARABLE:
            found = _find_nested_product(kernel)
            if found is not None and not found[1]:
                raise LikelihoodError(
                    f"{self.NAME} cannot lower a Product: a product of quasiseparable kernels is "
                    f"not quasiseparable. Where the factors act on different axes — which is "
                    f"what a Product is for — the result is not a function of one ordered "
                    f"coordinate at all, and where they act on the same one the semiseparable "
                    f"rank multiplies and is not recoverable from the factors' own "
                    f"representations. Use DenseGP, or replace the Product with a Sum, which is "
                    f"quasiseparable exactly when every term is."
                )
            if found is not None:
                _, path = found
                location = ".".join(path)
                raise LikelihoodError(
                    f"{self.NAME} cannot lower this {type(kernel).__name__}: its term "
                    f"{location!r} is a Product, and a product of quasiseparable kernels is not "
                    f"quasiseparable. Where the factors act on different axes — which is what a "
                    f"Product is for — the result is not a function of one ordered coordinate at "
                    f"all, and where they act on the same one the semiseparable rank multiplies "
                    f"and is not recoverable from the factors' own representations. Use DenseGP, "
                    f"or replace {location!r} with a Sum, which is quasiseparable exactly when "
                    f"every term is."
                )
        if self.REQUIRES_QUASISEPARABLE and not kernel.QUASISEPARABLE:
            raise LikelihoodError(
                f"{self.NAME} needs a kernel with an exact quasiseparable representation, but "
                f"{type(kernel).__name__} ({kernel.FAMILY}) has none. Use Matern32 — which is why "
                f"DEVELOPMENT_PLAN.md §2 made it the default — or switch to DenseGP."
            )
        if self.REQUIRES_ORDERED_1D:
            used = kernel.selected_axes([axis.name for axis in observed.axes])
            if len(used) != 1:
                raise LikelihoodError(
                    f"{self.NAME} needs one ordered coordinate axis, but a {kind} has "
                    f"{len(observed.axes)}: {[axis.name for axis in observed.axes]}."
                    + (
                        ""
                        if kernel.axes is None and not kernel.terms
                        else f" The kernel selects {list(used)}; select exactly one axis to "
                        f"reach the O(N) path."
                    )
                )
        if not self.IMPLEMENTED:
            raise LikelihoodError(self._unimplemented_message())

    def _unimplemented_message(self) -> str:
        return (
            f"{self.NAME} is a declared strategy slot with no implementation yet "
            f"(DEVELOPMENT_PLAN.md §4.4). Its interface is fixed so backends can fill it in; "
            f"until then use DenseGP, which computes the same quantity exactly."
        )

    @abc.abstractmethod
    def log_marginal_likelihood(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> float:
        """log N(residual; 0, K(values) + diag(variance))."""

    @abc.abstractmethod
    def condition(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
        at: np.ndarray | None = None,
    ) -> GPConditional:
        """The GP posterior at ``at`` (default: the data coordinates)."""

    @abc.abstractmethod
    def latent_transform(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        whitened: np.ndarray,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> np.ndarray:
        """Map whitened latent draws ``z`` to a GP draw ``f = L(θ) z``.

        This is the deterministic half of the latent-GP declaration (see
        :func:`latent_parameter`): the *prior* stays i.i.d. standard normal and
        all the correlation lives here, where the solver can impose it in
        whatever representation it uses.
        """

    def conditional_loo(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> np.ndarray:
        """Per-sample leave-one-out conditional log-density terms.

        The named GP decomposition of ``results.md`` §6 (ruled 2026-09-03,
        R2): term *i* is ``log N(y_i | mu_i^{-i}, sigma_i^{2,-i})`` — the
        density of sample *i* under the GP conditioned on every *other*
        retained sample — computable in closed form from the same factor
        :meth:`log_marginal_likelihood` uses. These terms are what
        ``arviz.loo``/``waic`` consume; they are a *different* decomposition
        from the factorised pointwise terms of independent noise and do not
        sum to the joint log-likelihood.

        The default refuses, naming the strategy — the same declared-slot
        discipline as an unimplemented solver. :class:`DenseGP` implements
        it; ``QuasisepGP`` owes an O(N) recursion in Phase 2.
        """
        raise LikelihoodError(
            f"{self.NAME} does not implement the leave-one-out conditional terms "
            f"(GPSolver.conditional_loo). DenseGP computes them exactly from its Cholesky; a "
            f"faster strategy must supply its own recursion before pointwise_log_prob can use "
            f"it."
        )

    def __repr__(self) -> str:
        return f"{type(self).__name__}()"


@dataclasses.dataclass(frozen=True)
class DenseGP(GPSolver):
    """Exact O(N³) dense Cholesky. The correctness anchor for every other solver.

    ``DEVELOPMENT_PLAN.md`` §4.4 calls this "current behaviour, any kernel,
    O(N³); correctness reference", and §4.6 makes "DenseGP↔QuasisepGP agreement
    on quasiseparable kernels" a conformance row. It has no performance goal
    whatsoever: it exists so that "what should the answer be?" has an
    implementation, not an argument.

    Parameters
    ----------
    jitter
        A standard deviation, in the data's own units, added in quadrature to
        the diagonal. Defaults to zero: a covariance that will not factorise is
        a fact about the model, and this contract does not hide it behind a
        silent numerical fudge. The error message tells you to set this when
        setting it is the right answer.
    """

    jitter: float = 0.0

    NAME: ClassVar[str] = "DenseGP"
    EXACT: ClassVar[bool] = True
    IMPLEMENTED: ClassVar[bool] = True
    STACKED_RESIDUALS: ClassVar[bool] = True

    def __post_init__(self) -> None:
        if not math.isfinite(self.jitter) or self.jitter < 0.0:
            raise LikelihoodError(f"DenseGP's jitter must be finite and >= 0, got {self.jitter!r}.")

    def _factor(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> tuple[np.ndarray, bool]:
        covariance = kernel.matrix(coordinates, coordinates, values)
        total = covariance + np.diag(variance + self.jitter**2)
        if not np.all(np.isfinite(total)):
            raise LikelihoodError(
                "the covariance matrix K + diag(sigma^2) contains non-finite entries, so it "
                "cannot be factorised. The usual cause is a kernel amplitude large enough that "
                "amplitude**2 overflows float64; constrain the amplitude prior to the data's "
                "own scale."
            )
        try:
            return scipy.linalg.cho_factor(total, lower=True)
        except scipy.linalg.LinAlgError as error:
            raise LikelihoodError(
                f"the covariance matrix K + diag(sigma^2) is not positive definite, so its "
                f"Cholesky factorisation failed ({error}). Usual causes: a zero or near-duplicate "
                f"observational uncertainty, coordinates that are closer together than float64 "
                f"can separate at this length-scale, or a kernel amplitude far above the data "
                f"scale. Pass DenseGP(jitter=...) — a standard deviation in the data's units — "
                f"if the matrix is merely ill-conditioned rather than wrong."
            ) from error

    def log_marginal_likelihood(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> float:
        """``log N(residual; 0, K + diag(variance))``, summed over the columns.

        **W4.2** generalises the right-hand side (see :class:`GPSolver`): with
        ``k`` columns this is the joint marginal of ``k`` independent
        realisations sharing one covariance, which is the circular complex
        Gaussian's 2N-dimensional density written without ever forming the
        2N by 2N matrix. One ``cho_factor``, ``k`` triangular solves, the
        log-determinant computed once and counted ``k`` times.
        """
        factor = self._factor(kernel, coordinates, variance, values)
        alpha = scipy.linalg.cho_solve(factor, residual)
        log_determinant = 2.0 * float(np.sum(np.log(np.abs(np.diag(factor[0])))))
        # ``np.sum(r * alpha)`` rather than ``r @ alpha``: the same number for a
        # vector, and the sum of the k quadratic forms for a block, whereas the
        # matrix product of two (n, k) blocks is not defined at all.
        quadratic = float(np.sum(residual * alpha))
        columns = _components(residual)
        return -0.5 * (quadratic + columns * log_determinant + residual.size * _LOG_2PI)

    def conditional_loo(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> np.ndarray:
        """The closed form from the same Cholesky (Sundararajan & Keerthi 2001).

        With ``A = (K + diag(variance + jitter^2))^{-1}``:
        ``sigma_i^{2,-i} = 1 / A_ii`` and ``mu_i^{-i} = y_i - [A r]_i / A_ii``,
        so ``log p_i = 0.5 log A_ii - [A r]_i^2 / (2 A_ii) - 0.5 log(2 pi)``.

        **W4.2**: with a ``(n, k)`` right-hand side there is still **one term
        per sample**, because a sample is what ``results.md`` §6's pointwise
        group is indexed by and a complex visibility is one sample with two
        components. The components share ``A``, hence share ``A_ii``, so term
        *i* is the sum of the ``k`` conditional densities at that sample and
        the normalisation and the log-precision are counted ``k`` times.
        """
        factor = self._factor(kernel, coordinates, variance, values)
        alpha = scipy.linalg.cho_solve(factor, residual)
        rows = np.shape(residual)[0]
        columns = _components(residual)
        precision_diagonal = np.diag(scipy.linalg.cho_solve(factor, np.eye(rows)))
        quadratic = np.sum(np.reshape(alpha, (rows, columns)) ** 2, axis=1)
        return np.asarray(
            columns * (0.5 * np.log(precision_diagonal) - 0.5 * _LOG_2PI)
            - quadratic / (2.0 * precision_diagonal),
            dtype=DTYPE,
        )

    def condition(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
        at: np.ndarray | None = None,
    ) -> GPConditional:
        points = _as_points(coordinates, "data coordinates")
        factor = self._factor(kernel, points, variance, values)
        if at is None:
            target = points
        else:
            target = _as_points(at, "conditioning grid", dimensions=points.shape[1])
        cross = kernel.matrix(target, points, values)
        # An (n, k) residual gives an (m, k) conditioned mean -- one column per
        # realisation -- beside the single (m,) variance the k realisations
        # share, since a posterior variance does not depend on the data at all.
        mean = cross @ scipy.linalg.cho_solve(factor, residual)
        solved = scipy.linalg.cho_solve(factor, cross.T)
        prior_variance = kernel.diagonal(target, values)
        posterior = prior_variance - np.einsum("ij,ji->i", cross, solved)
        return GPConditional(mean=mean, variance=posterior)

    def latent_transform(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        whitened: np.ndarray,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> np.ndarray:
        points = _as_points(coordinates, "data coordinates")
        covariance = kernel.matrix(points, points, values)
        if not np.all(np.isfinite(covariance)):
            raise LikelihoodError(
                "the kernel matrix K contains non-finite entries, so the whitening transform "
                "f = L z is undefined. The usual cause is a kernel amplitude large enough that "
                "amplitude**2 overflows float64."
            )
        scale = float(np.mean(np.diag(covariance))) or 1.0
        stabilised = covariance + np.eye(covariance.shape[0]) * (jitter * scale)
        try:
            lower = scipy.linalg.cholesky(stabilised, lower=True)
        except scipy.linalg.LinAlgError as error:
            raise LikelihoodError(
                f"the kernel matrix K is not positive definite, so the whitening transform "
                f"f = L z is undefined ({error}). Increase the jitter argument, or check the "
                f"kernel hyperparameters."
            ) from error
        return lower @ _as_float64(whitened, "whitened latent draws")


# ---------------------------------------------------------------------------
# celerite2 terms: the exact quasiseparable representations (W2.3, W4.5)
# ---------------------------------------------------------------------------


@functools.cache
def _ampere_term_type() -> Any:
    """The ``celerite2.terms.Term`` subclass that wraps an ampere kernel.

    **One wrapper, every family (W4.5).** Until W4.5 this module carried a
    hand-written ``celerite2`` term class per kernel — one, for Matérn-3/2,
    with its rank-2 algebra spelled out — and each differentiable backend
    carried a transcription of it. That does not scale to the five families
    W4.5 adds: five closed forms and five sets of generators would have become
    fifteen transcriptions of identical mathematics, with nothing checking that
    they stayed identical.

    So the mathematics moved to :mod:`ampere.core.kernels`, where each family's
    generators are a registered
    :class:`~ampere.core.kernels.CeleriteRepresentation` builder written once
    against :class:`~ampere.core.kernels.ArrayOps`, and this class became a
    *shim*: it asks the registry for the kernel's builder, calls it on the
    coordinates celerite2 hands it, and returns celerite2's ``(c, a, U, V)``.

    Built on first use rather than at import, because it has to subclass
    ``celerite2.terms.Term`` and ``ampere.core`` does not import celerite2 at
    module level (see :class:`QuasisepGP`). :func:`functools.cache` makes the
    class a singleton, so ``isinstance`` and celerite2's own caches behave.

    **Why ampere supplies its own term rather than using celerite2's.**
    celerite2's ``Matern32Term`` is an *approximation* and says so: the
    celerite basis ``e^{-c t}(a cos d t + b sin d t)`` has no ``t e^{-c t}``
    member, so its Matérn-3/2 is a limit in a parameter ``eps`` which at the
    default ``eps=0.01`` costs about 5e-3 in the log-likelihood — three orders
    of magnitude outside the conformance battery's ``cross_solver`` tolerance.
    The *solver* underneath needs no such approximation: it factorises any
    rank-J semiseparable matrix, and Matérn-1/2, -3/2 and -5/2, the SHO and the
    rotation pair all have exact representations in that form. See
    :class:`~ampere.core.kernels.CeleriteRepresentation` for the general shape
    and each ``*_representation`` function for its algebra.
    """
    # Imported here rather than at module level: celerite2 is a base
    # dependency (architecture.md §2), but ampere.core promises to import
    # nothing beyond numpy/scipy/astropy/stdlib, and a problem with no
    # quasiseparable GP in it should not pay ~20 ms to load a C extension it
    # never calls. Same idiom as the reference backend's pyphot import.
    import celerite2.terms

    class _AmpereTerm(celerite2.terms.Term):
        """An ampere :class:`Kernel` presented as a celerite2 term."""

        def __init__(self, kernel: Kernel, values: Mapping[str, Any]) -> None:
            self.kernel = kernel
            self.values = dict(values)
            self.builder = lookup_quasiseparable_term(kernel.FAMILY)

        def get_value(self, tau: Any) -> np.ndarray:
            separation = np.abs(np.atleast_1d(np.asarray(tau, dtype=DTYPE)))
            return np.asarray(self.kernel.value(separation, self.values), dtype=DTYPE)

        def get_celerite_matrices(
            self,
            x: Any,
            diag: Any,
            *,
            c: Any = None,
            a: Any = None,
            U: Any = None,
            V: Any = None,
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
            points = np.ascontiguousarray(np.atleast_1d(np.asarray(x, dtype=DTYPE)))
            diagonal = np.ascontiguousarray(np.atleast_1d(np.asarray(diag, dtype=DTYPE)))
            # x arrives sorted (QuasisepGP sorts before calling), which is what
            # the builders' midpoint centring assumes.
            representation = self.builder(self.kernel, self.values, points)
            return (
                np.ascontiguousarray(representation.decay, dtype=DTYPE),
                np.ascontiguousarray(diagonal + representation.marginal, dtype=DTYPE),
                np.ascontiguousarray(representation.left, dtype=DTYPE),
                np.ascontiguousarray(representation.right, dtype=DTYPE),
            )

    return _AmpereTerm


def _celerite_term(kernel: Kernel, values: Mapping[str, Any]) -> Any:
    """Build the exact celerite representation of a resolved kernel.

    The overflow guard is here rather than inside a builder because it applies
    to every family: ``k(0)`` is ``amplitude**2`` throughout ampere, so an
    amplitude above about 1e154 makes the marginal variance — and with it every
    generator — non-finite, and celerite2's factorisation returns quiet NaN
    where a Cholesky raises. A composite's amplitudes are qualified
    (``term0.amplitude``), so the last name segment is what is matched.
    """
    # The reference path computes in numpy whatever namespace the kernel was
    # declared in, which is what ``with_ops`` exists for (a torch kernel handed
    # to this solver must not build tensors celerite2's numpy driver cannot
    # read). Ordinary reference kernels are already numpy, and get themselves
    # back unchanged.
    kernel = kernel.with_ops(NUMPY_OPS)
    resolved = kernel.resolve(values)
    for name, value in resolved.items():
        if name.rsplit(".", 1)[-1] != "amplitude":
            continue
        amplitude = float(np.asarray(value, dtype=DTYPE))
        if not math.isfinite(amplitude * amplitude):
            raise LikelihoodError(
                f"the kernel's marginal variance amplitude**2 = {amplitude!r}**2 overflows "
                f"float64, so the quasiseparable representation cannot be built. Constrain the "
                f"amplitude prior to the data's own scale."
            )
    return _ampere_term_type()(kernel, values)


class _SolverSlot(GPSolver):
    """Base for the declared-but-unimplemented strategies."""

    IMPLEMENTED: ClassVar[bool] = False

    def log_marginal_likelihood(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> float:
        raise LikelihoodError(self._unimplemented_message())

    def condition(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
        at: np.ndarray | None = None,
    ) -> GPConditional:
        raise LikelihoodError(self._unimplemented_message())

    def latent_transform(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        whitened: np.ndarray,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> np.ndarray:
        raise LikelihoodError(self._unimplemented_message())


@dataclasses.dataclass(frozen=True)
class QuasisepGP(GPSolver):
    """Exact O(N) for ordered 1D data via a quasiseparable (celerite-class) solve.

    The strategy ``DEVELOPMENT_PLAN.md`` §4.4 names as the scaling answer, and
    the reason ``architecture.md`` §2 puts celerite2's numpy interface in the
    **base** install rather than behind an extra: ampere's base install --
    ``pip install .`` from a checkout (PyPI's ``ampere`` package is
    unrelated) -- with no extras, must be a scalable fitting environment, not
    one that still has the O(N³) problem. It is **exact**, not approximate: a
    Matérn-3/2 kernel has an exact rank-2 semiseparable representation (see
    :func:`~ampere.core.kernels.matern32_representation` for the algebra), so
    this recursion computes the same marginal likelihood :class:`DenseGP` does,
    in linear time. That equivalence is a conformance row (§4.6), which is why
    :class:`DenseGP` exists at all. Since **W4.5** the same holds for
    Matérn-1/2 and -5/2, the SHO, the rotation pair, and any :class:`Sum` of
    them: a sum of quasiseparable terms is quasiseparable, at a rank that is
    the sum of the terms'.

    A kernel reaches this path only if it declares ``QUASISEPARABLE`` **and**
    every family in its tree has a representation in the public registry
    (:func:`~ampere.core.kernels.register_quasiseparable_term`); anything else
    is refused by name at composition time — including a
    :class:`~ampere.core.kernels.Product`, which gets its own refusal because
    the reason is structural rather than a missing row. Coordinates need not
    arrive sorted — a Gaussian density
    is invariant under a simultaneous permutation of residuals, variances and
    coordinates, so this solver sorts internally and undoes the permutation on
    the way out.

    What is O(N) and what is not, stated plainly:

    * :meth:`log_marginal_likelihood` and :meth:`latent_transform` are O(N),
      which is what a sampler calls.
    * :meth:`condition` is O(N·M) for M evaluation points, because the
      cross-covariance block is dense by construction and each output needs
      every input. :class:`DenseGP` pays O(N³) for the same answer.
    * :meth:`conditional_loo` is **deferred** (recorded in
      ``DEVELOPMENT_PLAN.md`` §2, 2026-09-05) and refuses: the leave-one-out
      terms need the diagonal of ``(K + diag(sigma²))⁻¹``, and celerite2's
      public numpy API exposes no O(N) route to it. Use :class:`DenseGP` for
      ``pointwise_log_prob`` under a GP.

    Parameters
    ----------
    jitter
        A standard deviation, in the data's own units, added in quadrature to
        the diagonal — the same knob, meaning and default as
        :class:`DenseGP`'s. Zero by default: a covariance that will not
        factorise is a fact about the model, not something to hide.
    """

    jitter: float = 0.0

    NAME: ClassVar[str] = "QuasisepGP"
    EXACT: ClassVar[bool] = True
    REQUIRES_ORDERED_1D: ClassVar[bool] = True
    REQUIRES_QUASISEPARABLE: ClassVar[bool] = True
    IMPLEMENTED: ClassVar[bool] = True

    def __post_init__(self) -> None:
        if not math.isfinite(self.jitter) or self.jitter < 0.0:
            raise LikelihoodError(
                f"QuasisepGP's jitter must be finite and >= 0, got {self.jitter!r}."
            )

    def check_compatible(self, kernel: Kernel, observed: FunctionSamples) -> None:
        super().check_compatible(kernel, observed)
        # Every family in the tree needs a registered representation, not just
        # the root: a Sum lowers term by term (``sum_representation``), so one
        # unregistered term is enough to stop it, and it must be named here
        # rather than discovered inside the first factorisation.
        for leaf in kernel.leaves():
            lookup_quasiseparable_term(leaf.FAMILY, owner=self.NAME)

    # -- internals -----------------------------------------------------------

    def _axis(
        self, coordinates: np.ndarray, kernel: Kernel | None = None
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """The ``(n, 1)`` points, their bare axis, and the sorting permutation.

        ``kernel`` is passed since W4.5 so an axis-selecting kernel gets its own
        column out of a multi-axis container's coordinates. It is optional so a
        direct call with already-1-D coordinates still works.
        """
        points = _as_points(coordinates, "data coordinates")
        selected = points if kernel is None else kernel.select(points)
        if selected.shape[1] != 1:
            raise LikelihoodError(
                f"{self.NAME} needs one ordered coordinate per sample, but the coordinates have "
                f"{selected.shape[1]} per point. check_compatible refuses this at composition "
                f"time; a direct solver call reaches it here. Use DenseGP for 2D+ coordinates."
            )
        axis = np.ascontiguousarray(selected[:, 0])
        return points, axis, np.argsort(axis, kind="stable")

    def _factorise(
        self,
        kernel: Kernel,
        ordered_axis: np.ndarray,
        ordered_diagonal: np.ndarray,
        values: Mapping[str, Any],
        *,
        whitening: bool = False,
    ) -> Any:
        """A celerite2 ``GaussianProcess`` factorised on sorted coordinates."""
        # Lazy, for the reason _ampere_term_type() gives.
        import celerite2

        # celerite2.driver is the compiled extension: no stubs, so pyrefly
        # cannot resolve it. The name is public and documented.
        from celerite2.driver import LinAlgError  # type: ignore[missing-import]

        # celerite2's factorisation does not itself notice a diagonal that
        # cannot belong to a covariance: it returns NaN quietly, where
        # DenseGP's Cholesky raises. Check the precondition instead, so both
        # strategies refuse the same inputs with the same kind of message.
        if not np.all(np.isfinite(ordered_diagonal)):
            raise LikelihoodError(
                "the covariance matrix K + diag(sigma^2) contains non-finite entries, so it "
                "cannot be factorised. The usual cause is a kernel amplitude large enough that "
                "amplitude**2 overflows float64; constrain the amplitude prior to the data's "
                "own scale."
            )
        if np.any(ordered_diagonal < 0.0):
            raise LikelihoodError(
                "the diagonal handed to QuasisepGP contains negative entries, so K + "
                "diag(sigma^2) is not a covariance matrix at all. Variances are squares; this "
                "is a caller error rather than an ill-conditioned problem."
            )

        term = _celerite_term(kernel, values)
        gp = celerite2.GaussianProcess(term, mean=0.0)
        try:
            gp.compute(ordered_axis, diag=ordered_diagonal, check_sorted=False)
        except LinAlgError as error:
            if whitening:
                raise LikelihoodError(
                    f"the kernel matrix K is not positive definite in its quasiseparable "
                    f"representation, so the whitening transform f = L z is undefined "
                    f"({error}). Increase the jitter argument, or check the kernel "
                    f"hyperparameters."
                ) from error
            raise LikelihoodError(
                f"the covariance matrix K + diag(sigma^2) is not positive definite, so its "
                f"quasiseparable factorisation failed ({error}). Usual causes: a zero or "
                f"near-duplicate observational uncertainty, coordinates that are closer "
                f"together than float64 can separate at this length-scale, or a kernel "
                f"amplitude far above the data scale. Pass QuasisepGP(jitter=...) — a standard "
                f"deviation in the data's units — if the matrix is merely ill-conditioned "
                f"rather than wrong."
            ) from error
        return gp

    # -- the interface -------------------------------------------------------

    def log_marginal_likelihood(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> float:
        _, axis, order = self._axis(coordinates, kernel)
        residuals = _as_float64(residual, "residuals")
        diagonal = _as_float64(variance, "noise variances") + self.jitter**2
        gp = self._factorise(kernel, axis[order], diagonal[order], values)
        return float(gp.log_likelihood(residuals[order]))

    def condition(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
        at: np.ndarray | None = None,
    ) -> GPConditional:
        points, axis, order = self._axis(coordinates, kernel)
        residuals = _as_float64(residual, "residuals")
        diagonal = _as_float64(variance, "noise variances") + self.jitter**2
        gp = self._factorise(kernel, axis[order], diagonal[order], values)

        alpha = np.empty(order.size, dtype=DTYPE)
        alpha[order] = gp.apply_inverse(residuals[order])
        target = points if at is None else _as_points(at, "conditioning grid", dimensions=1)
        # The cross-covariance is dense whatever the solver: M outputs each
        # need all N inputs. Only the solve against it is O(N) per column.
        cross = kernel.matrix(target, points, values)
        solved = np.empty((order.size, cross.shape[0]), dtype=DTYPE)
        solved[order, :] = gp.apply_inverse(np.ascontiguousarray(cross.T[order, :]))
        prior_variance = kernel.diagonal(target, values)
        return GPConditional(
            mean=cross @ alpha,
            variance=prior_variance - np.einsum("ij,ji->i", cross, solved),
        )

    def latent_transform(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        whitened: np.ndarray,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> np.ndarray:
        points, axis, order = self._axis(coordinates, kernel)
        draws = _as_float64(whitened, "whitened latent draws")
        # The same stabilisation DenseGP applies: a jitter relative to the
        # kernel's own scale, so the two solvers factorise the same matrix.
        scale = float(np.mean(kernel.diagonal(points, values))) or 1.0
        gp = self._factorise(
            kernel,
            axis[order],
            np.full(order.size, jitter * scale, dtype=DTYPE),
            values,
            whitening=True,
        )
        transformed = np.empty(order.size, dtype=DTYPE)
        transformed[order] = gp.dot_tril(draws[order])
        return transformed

    def conditional_loo(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> np.ndarray:
        """Deferred at W2.3 (``DEVELOPMENT_PLAN.md`` §2, 2026-09-05).

        Every leave-one-out term needs ``A_ii`` for ``A = (K +
        diag(sigma²))⁻¹``, and celerite2's public numpy interface has no O(N)
        route to that diagonal — its own ``condition(...).variance`` forms the
        cross-covariance densely and costs O(N·M). Supplying one means
        reimplementing celerite2's internal factorisation convention, which is
        a coupling this contract declined to take on for a decomposition
        nothing yet stores by default.
        """
        raise LikelihoodError(
            f"{self.NAME} does not implement the leave-one-out conditional terms "
            f"(GPSolver.conditional_loo): the O(N) recursion for the diagonal of "
            f"(K + diag(sigma^2))^-1 was deferred at W2.3 and recorded in "
            f"DEVELOPMENT_PLAN.md §2. DenseGP computes them exactly from its Cholesky — use it "
            f"for pointwise_log_prob under a GP, or the per-dataset log_likelihood group."
        )


@dataclasses.dataclass(frozen=True)
class HilbertSpaceGP(GPSolver):
    r"""Approximate ``O(N m + m³)`` by a Hilbert-space (reduced-rank spectral) basis.

    **The first ``EXACT = False`` solver ampere implements** (W5.4), settling
    the Phase 5 bullet ``DEVELOPMENT_PLAN.md`` §5 opens with and the two
    contract questions ``docs/design/horizon_notes.md`` §2 left for the phase.
    The method is Solin & Särkkä (2020), in the practical form Riutort-Mayol
    et al. (2023) give: on a box that contains the data, the Dirichlet
    Laplacian's eigenfunctions diagonalise any stationary kernel through its
    spectral density, so

    .. math::
        K(\theta) \;\approx\; \Phi \operatorname{diag}\!\big(S_\theta\big)
        \Phi^{\mathsf T},

    with :math:`\Phi` a **fixed** ``(N, m)`` block of sines — it depends on the
    coordinates and the box, never on a fitted value — and :math:`S_\theta` the
    kernel's spectral density at the box's eigenvalues, which is where all the
    hyperparameter dependence lives. Writing :math:`\tilde\Phi = \Phi
    \operatorname{diag}(\sqrt{S_\theta})` and :math:`D = \operatorname{diag}
    (\sigma^2 + \text{jitter}^2)`, Woodbury gives every quantity this contract
    asks for from one ``m`` by ``m`` Cholesky of :math:`M = I_m + \tilde\Phi^{\mathsf
    T} D^{-1} \tilde\Phi`:

    * :math:`\log|K + D| = \log|D| + \log|M|`;
    * :math:`(K + D)^{-1} r = D^{-1} r - D^{-1}\tilde\Phi\, M^{-1}
      \tilde\Phi^{\mathsf T} D^{-1} r`;
    * :math:`\tilde\Phi^{\mathsf T} (K+D)^{-1} \tilde\Phi = I_m - M^{-1}`,
      which collapses the conditional variance to
      :math:`\|L_M^{-1}\tilde\Phi_*^{\mathsf T}\|^2` — non-negative by
      construction, which a subtracted quadratic form is not.

    **What it is for.** Three things, in order of how much they matter.
    It is the scaling answer for a kernel with *no* quasiseparable form —
    :class:`~ampere.core.kernels.SquaredExponential` above all, where the
    spectral density decays as a Gaussian and the approximation converges
    exponentially in ``m``. It is the scaling answer in **two and three axes**,
    where :class:`QuasisepGP` does not apply at all: the basis is a tensor
    product, so ``m = m₁·m₂(·m₃)`` and the method stays cheap exactly while the
    dimension stays low. And it is the natural **latent** representation under
    NUTS, because the whitened block it needs is ``m`` variables rather than
    ``N`` — see :meth:`latent_size`.

    **What it is not.** It is not exact, and the error does not go to zero in
    ``m`` alone: a finite box has a boundary, and the approximation is poor
    within about one length scale of it, so ``basis_size`` and
    ``boundary_factor`` must grow together. Riutort-Mayol et al. §3 give the
    rule of thumb — roughly ``c >= 1.2 max(length_scale)/S`` for a data
    half-extent ``S``, and ``m`` large enough that :math:`\pi m / (2 c S)`
    reaches several times the inverse length scale. ampere does not choose
    either for you, and will not: they *are* the
    approximation, and an approximation chosen silently is the failure mode
    ``DEVELOPMENT_PLAN.md`` §7 warns about. What ampere does instead is make
    them *visible* — they are dataclass fields, so ``Likelihood.to_spec``
    records them and two runs at different ``m`` have different spec hashes,
    and they are repeated in :meth:`provenance_config` so a stored run says in
    its attrs what approximation produced it.

    **The conformance battery holds it to a convergence claim, not a number**
    (horizon notes §2, question (a)): the marginal likelihood, the conditioned
    moments and the leave-one-out terms are compared against :class:`DenseGP`
    at a sequence of ``m``, inside an envelope that tightens as ``m`` grows and
    down to a floor set by the box. See ``tests/conformance/README.md`` §3.

    Parameters
    ----------
    basis_size
        Basis members **per axis**: an integer for a one-axis kernel, or one
        entry per selected axis. The total ``m`` is their product, and it is
        this declaration — not the sample count — that fixes the latent
        block's size, which is why it must be expressible without the data in
        hand (``inference.md`` §17.4, amended at W5.4).
    boundary_factor
        The box's half-width on each axis, as a multiple of the data's own
        half-extent. See :data:`~ampere.core.hsgp.DEFAULT_BOUNDARY_FACTOR`.
    jitter
        A standard deviation added in quadrature to the diagonal, the same
        knob, meaning and default as :class:`DenseGP`'s.

    Examples
    --------
    >>> solver = HilbertSpaceGP(basis_size=64, boundary_factor=2.0)
    >>> solver.basis_size
    (64,)
    >>> solver.latent_size(Matern32(0.3, 2.0), 500)
    64
    >>> HilbertSpaceGP(basis_size=(16, 16)).provenance_config()
    {'basis_size': [16, 16], 'boundary_factor': 2.0}
    """

    basis_size: int | tuple[int, ...] = DEFAULT_BASIS_SIZE
    boundary_factor: float = DEFAULT_BOUNDARY_FACTOR
    jitter: float = 0.0

    NAME: ClassVar[str] = "HilbertSpaceGP"
    #: The whole point. Every row that compares this solver against
    #: :class:`DenseGP` reads it and asks for a convergence tolerance.
    EXACT: ClassVar[bool] = False
    IMPLEMENTED: ClassVar[bool] = True
    #: The ``(n, k)`` right-hand side of the circular complex GP: the Woodbury
    #: solve takes extra columns exactly as the dense Cholesky does — one
    #: factorisation of ``M``, ``k`` solves, the log-determinant counted ``k``
    #: times.
    STACKED_RESIDUALS: ClassVar[bool] = True

    def __post_init__(self) -> None:
        object.__setattr__(self, "basis_size", normalise_counts(self.basis_size))
        if not math.isfinite(self.boundary_factor) or self.boundary_factor <= 0.0:
            raise LikelihoodError(
                f"HilbertSpaceGP's boundary_factor must be finite and > 0, got "
                f"{self.boundary_factor!r}. It is the box's half-width as a multiple of the "
                f"data's own half-extent, so it must be at least 1 for the box to contain the "
                f"data at all."
            )
        if not math.isfinite(self.jitter) or self.jitter < 0.0:
            raise LikelihoodError(
                f"HilbertSpaceGP's jitter must be finite and >= 0, got {self.jitter!r}."
            )

    @property
    def counts(self) -> tuple[int, ...]:
        """Basis members per axis, normalised to a tuple."""
        return normalise_counts(self.basis_size)

    # -- declarations --------------------------------------------------------

    def provenance_config(self) -> Mapping[str, Any]:
        """The approximation, for the run's attrs (``inference.md`` §10a, fold-in 10).

        Both numbers are *also* dataclass fields, and that is deliberate rather
        than a duplication. Fold-in 10's test is whether two backends may
        legitimately differ on a value: they may on a dtype or a device, which
        is why those stay out of the spec, and they may **not** on ``m`` or
        ``c`` — a torch and a jax run of one declared problem must build the
        same basis or they are not computing the same quantity. So the
        approximation is part of the declaration and enters
        ``ampere_spec_hash``, where it belongs: two runs at ``m = 8`` and
        ``m = 64`` are not the same model and must not hash alike. It is
        repeated here because ``ampere_solver_config`` is where a reader looks
        for "how was this computed", and an approximation is the first thing
        they should find there.
        """
        return {
            "basis_size": [int(count) for count in self.counts],
            "boundary_factor": float(self.boundary_factor),
        }

    def latent_size(self, kernel: Kernel, n_samples: int) -> int:
        """``m``, not ``N`` — the ruling W5.4 makes (horizon notes §2, question (b)).

        ``inference.md`` §17.4 fixes the latent block's size at composition,
        and said in passing that it was one value per retained sample. The
        second half was a property of the two exact solvers, not of the
        contract: the whitened block is whatever the solver's whitening takes,
        and this solver's takes ``m`` basis coefficients. ``N`` plays no part,
        which is the whole reason the method is worth having under NUTS —
        a 10⁴-sample spectrum fits with a few dozen latent dimensions instead
        of ten thousand.

        The size is read off the **declaration** alone, never off the data, so
        it is available where ``Likelihood.latent_declaration`` needs it and
        cannot drift from what :meth:`latent_transform` will accept.
        """
        del kernel, n_samples
        return basis_size(self.counts)

    def check_compatible(self, kernel: Kernel, observed: FunctionSamples) -> None:
        """The layout and axis rules, plus this solver's own two.

        The kernel must have a closed-form spectral density in every leaf
        (:func:`~ampere.core.hsgp.check_spectral_support` names the one that
        does not), and ``basis_size`` must have one entry per axis the kernel
        selects — the count that fixes the latent size, checked here against
        the container rather than discovered as a shape error inside a solve.
        """
        super().check_compatible(kernel, observed)
        selected = kernel.selected_axes([axis.name for axis in observed.axes])
        dimensions = len(selected)
        check_spectral_support(kernel, dimensions, owner=self.NAME)
        counts = self.counts
        if len(counts) != dimensions:
            raise LikelihoodError(
                f"{self.NAME} was declared with basis_size={self.basis_size!r} — "
                f"{len(counts)} axis count(s) — but the kernel selects {dimensions} axis/axes "
                f"{selected!r} of this {type(observed).__name__}. The basis is a tensor product "
                f"with one count per axis, and the total m is their product, so the two must "
                f"agree: pass basis_size={tuple([counts[0]] * dimensions)!r} for the same "
                f"resolution on each."
            )

    # -- internals -----------------------------------------------------------

    def _basis(self, kernel: Kernel, points: np.ndarray) -> HilbertSpaceBasis:
        return hilbert_basis(kernel.select(points), self.counts, self.boundary_factor)

    def _scaled_basis(
        self,
        kernel: Kernel,
        basis: HilbertSpaceBasis,
        points: np.ndarray,
        values: Mapping[str, Any],
    ) -> np.ndarray:
        r""":math:`\tilde\Phi = \Phi\operatorname{diag}(\sqrt{S_\theta})`, ``(n, m)``.

        Carried as one block rather than as ``Phi`` and ``S`` separately
        because every use wants the product, and because scaling the basis
        keeps a spectral density that has underflowed to zero out of a
        denominator: an unused basis member contributes a zero column and a
        unit diagonal to ``M``, which is exactly right.
        """
        density = np.asarray(spectral_values(kernel, basis, values), dtype=DTYPE)
        if not np.all(np.isfinite(density)) or np.any(density < 0.0):
            raise LikelihoodError(
                "the kernel's spectral density is not finite and non-negative at the basis "
                "frequencies, so the reduced-rank factorisation K ~ Phi diag(S) Phi^T is not a "
                "covariance. The usual cause is a kernel amplitude or length scale outside its "
                "prior's support; a spectral density is non-negative for every admissible "
                "hyperparameter (Bochner)."
            )
        matrix = np.asarray(basis_matrix(basis, kernel.select(points), NUMPY_OPS), dtype=DTYPE)
        return matrix * np.sqrt(density)[None, :]

    def _diagonal(self, variance: np.ndarray) -> np.ndarray:
        diagonal = _as_float64(variance, "noise variances") + self.jitter**2
        if not np.all(np.isfinite(diagonal)) or np.any(diagonal <= 0.0):
            raise LikelihoodError(
                f"{self.NAME} needs a strictly positive noise diagonal: the Woodbury identity it "
                f"solves by inverts diag(sigma^2 + jitter^2) directly, so a zero uncertainty is "
                f"a division by zero rather than an ill-conditioned matrix. Pass "
                f"{self.NAME}(jitter=...) — a standard deviation in the data's units — or use "
                f"DenseGP, whose Cholesky needs only K + diag(sigma^2) to be positive definite."
            )
        return diagonal

    def _factor(self, scaled: np.ndarray, diagonal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        r"""``(D⁻¹ tilde-Phi, cho_factor(M))`` — the one factorisation everything uses."""
        weighted = scaled / diagonal[:, None]
        capacitance = np.eye(scaled.shape[1], dtype=DTYPE) + scaled.T @ weighted
        try:
            factor = scipy.linalg.cho_factor(capacitance, lower=True)
        except scipy.linalg.LinAlgError as error:
            raise LikelihoodError(
                f"the reduced-rank capacitance matrix I + Phi^T D^-1 Phi is not positive "
                f"definite, so the Woodbury solve failed ({error}). It is positive definite for "
                f"every admissible hyperparameter, so this is a numerical rather than a "
                f"structural failure: reduce basis_size, or raise the jitter."
            ) from error
        return weighted, factor

    def _solve(
        self,
        scaled: np.ndarray,
        diagonal: np.ndarray,
        weighted: np.ndarray,
        factor: tuple[np.ndarray, bool],
        right: np.ndarray,
    ) -> np.ndarray:
        r"""``(K + D)⁻¹ right`` by Woodbury, for an ``(n,)`` or ``(n, k)`` block."""
        divisor = diagonal if right.ndim == 1 else diagonal[:, None]
        direct = right / divisor
        return direct - weighted @ scipy.linalg.cho_solve(factor, scaled.T @ direct)

    # -- the interface -------------------------------------------------------

    def log_marginal_likelihood(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> float:
        """The marginal likelihood **of the approximation**, summed over the columns.

        Not of the declared kernel: :attr:`EXACT` is ``False`` and this is what
        that means. It is the exact marginal likelihood of the rank-``m``
        Gaussian process ``Phi diag(S) Phi^T``, which converges to the declared
        one as the basis fills the box.
        """
        points = _as_points(coordinates, "data coordinates")
        residuals = _as_float64(residual, "residuals")
        basis = self._basis(kernel, points)
        scaled = self._scaled_basis(kernel, basis, points, values)
        diagonal = self._diagonal(variance)
        weighted, factor = self._factor(scaled, diagonal)
        alpha = self._solve(scaled, diagonal, weighted, factor, residuals)
        log_determinant = float(np.sum(np.log(diagonal))) + 2.0 * float(
            np.sum(np.log(np.abs(np.diag(factor[0]))))
        )
        quadratic = float(np.sum(residuals * alpha))
        columns = _components(residuals)
        return -0.5 * (quadratic + columns * log_determinant + residuals.size * _LOG_2PI)

    def conditional_loo(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> np.ndarray:
        r"""Sundararajan & Keerthi's identity, **exact in the approximation**.

        The item that landed this solver allowed it to be refused if the
        mathematics genuinely failed; it does not. Every leave-one-out term
        needs one number the dense path gets from its Cholesky — the diagonal
        of :math:`A^{-1}`, for :math:`A = K + D` — and Woodbury supplies it in
        closed form without ever forming :math:`A`:

        .. math::
            (A^{-1})_{ii} = \frac{1}{D_i}
                - \big\| L_M^{-1} \tilde\Phi_i^{\mathsf T} / D_i \big\|^2 ,

        at ``O(N m²)``. The identity is then the dense one, applied to the
        rank-``m`` covariance this solver actually scores, so the terms are as
        exact as its marginal likelihood is and converge with it. **W4.2**'s
        rule holds unchanged: one term per *sample*, the ``k`` components
        summed.
        """
        points = _as_points(coordinates, "data coordinates")
        residuals = _as_float64(residual, "residuals")
        basis = self._basis(kernel, points)
        scaled = self._scaled_basis(kernel, basis, points, values)
        diagonal = self._diagonal(variance)
        weighted, factor = self._factor(scaled, diagonal)
        alpha = self._solve(scaled, diagonal, weighted, factor, residuals)
        # ``weighted`` is D^-1 tilde-Phi, whose rows are the u_i above.
        triangular = scipy.linalg.solve_triangular(factor[0], weighted.T, lower=True)
        precision_diagonal = 1.0 / diagonal - np.sum(triangular * triangular, axis=0)
        rows = np.shape(residuals)[0]
        columns = _components(residuals)
        quadratic = np.sum(np.reshape(alpha, (rows, columns)) ** 2, axis=1)
        return np.asarray(
            columns * (0.5 * np.log(precision_diagonal) - 0.5 * _LOG_2PI)
            - quadratic / (2.0 * precision_diagonal),
            dtype=DTYPE,
        )

    def condition(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
        at: np.ndarray | None = None,
    ) -> GPConditional:
        r"""The posterior of the rank-``m`` process, at ``at`` or at the data.

        The prior variance is the **approximation's** — the row sums of
        :math:`\tilde\Phi^2` — rather than the declared kernel's ``k(0)``, and
        deliberately so: mixing the exact prior with the approximate posterior
        correction gives a variance that is not the variance of anything, and
        can go negative near the box's edge. Taken consistently, the posterior
        variance is :math:`\|L_M^{-1}\tilde\Phi_*^{\mathsf T}\|^2`, which is
        non-negative whatever ``m`` is.
        """
        points = _as_points(coordinates, "data coordinates")
        residuals = _as_float64(residual, "residuals")
        basis = self._basis(kernel, points)
        scaled = self._scaled_basis(kernel, basis, points, values)
        diagonal = self._diagonal(variance)
        weighted, factor = self._factor(scaled, diagonal)
        alpha = self._solve(scaled, diagonal, weighted, factor, residuals)
        if at is None:
            target = points
        else:
            target = _as_points(at, "conditioning grid", dimensions=points.shape[1])
        target_scaled = self._scaled_basis(kernel, basis, target, values)
        mean = target_scaled @ (scaled.T @ alpha)
        triangular = scipy.linalg.solve_triangular(factor[0], target_scaled.T, lower=True)
        posterior = np.sum(triangular * triangular, axis=0)
        return GPConditional(mean=mean, variance=np.asarray(posterior, dtype=DTYPE))

    def latent_transform(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        whitened: np.ndarray,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> np.ndarray:
        r"""``f = tilde-Phi z`` — the reduced-rank whitening, from ``m`` variables.

        The same factor the marginal-likelihood path scores against, which is
        what makes :func:`~ampere.core.dataset.FittingProblem.simulate` with
        ``observe=True`` and the latent-GP likelihood path agree under this
        solver: both go through here, and there is only one :math:`\tilde\Phi`.

        ``z`` has :meth:`latent_size` entries, not one per sample. The
        stabilising ``jitter`` argument the two exact solvers take is accepted
        and **unused**: :math:`\tilde\Phi\tilde\Phi^{\mathsf T}` is positive
        semi-definite by construction with no factorisation to stabilise,
        which is one of the quieter benefits of a reduced-rank representation.
        """
        del jitter
        points = _as_points(coordinates, "data coordinates")
        draws = _as_float64(whitened, "whitened latent draws")
        basis = self._basis(kernel, points)
        scaled = self._scaled_basis(kernel, basis, points, values)
        if np.shape(draws)[0] != scaled.shape[1]:
            raise LikelihoodError(
                f"{self.NAME} whitens {scaled.shape[1]} basis coefficient(s), but it was handed "
                f"{np.shape(draws)[0]} whitened value(s). This solver's latent block is the "
                f"basis size m, not the sample count — see GPSolver.latent_size."
            )
        return scaled @ draws


class WindowedSparseGP(_SolverSlot):
    """Approximate O(N) by tapering the kernel to zero beyond a cutoff radius.

    Starfish's actual strategy (``prior_art.md`` lesson S2): apply a Hann
    window to the kernel so it is *exactly* zero beyond ``r₀ ≈ 4l``, giving a
    genuinely sparse banded covariance that a sparse Cholesky factorises in
    roughly linear time. This is a qualitatively different family from
    :class:`QuasisepGP` — approximate but simple, needing no state-space
    machinery, and applicable to kernels with no quasiseparable form at all
    (including :class:`SquaredExponential`) and to more than one dimension.

    It is a named slot rather than the default because :class:`QuasisepGP` is
    exact for the kernel ampere actually recommends, and buys that exactness
    without the cutoff radius — one more approximation-controlling
    hyperparameter that nothing chooses for the user. See the likelihoods
    contract for the full argument.
    """

    NAME: ClassVar[str] = "WindowedSparseGP"
    EXACT: ClassVar[bool] = False


class InducingPointGP(_SolverSlot):
    """Sparse variational / inducing-point GP (SVGP): a Phase 5 slot for 2D+."""

    NAME: ClassVar[str] = "InducingPointGP"
    EXACT: ClassVar[bool] = False


class StructuredGridGP(_SolverSlot):
    """SKI / KISS-GP: structured-kernel interpolation. A Phase 5 slot for 2D+."""

    NAME: ClassVar[str] = "StructuredGridGP"
    EXACT: ClassVar[bool] = False


class VecchiaGP(_SolverSlot):
    """Vecchia / nearest-neighbour approximation. A Phase 5 slot for 2D+."""

    NAME: ClassVar[str] = "VecchiaGP"
    EXACT: ClassVar[bool] = False


# ---------------------------------------------------------------------------
# Noise models
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class NoiseParams:
    """What a :class:`NoiseModel` hands a family for one evaluation.

    Every array here covers the **retained** samples only: masking has already
    been applied by excision (see :meth:`Likelihood.log_prob`), so a family
    never sees a mask and cannot forget to honour one.
    """

    #: Per-sample standard deviation, float64. ``None`` when the family
    #: supplies its own dispersion (Poisson).
    sigma: np.ndarray | None
    #: Resolved values of every parameter the likelihood declares.
    values: Mapping[str, Any]
    #: Coordinates, ``(n, d)`` float64, when a correlated solve is needed.
    coordinates: np.ndarray | None = None
    kernel: Kernel | None = None
    solver: GPSolver | None = None
    #: Latent function values ``f``, when the combination is
    #: :attr:`Marginalisation.LATENT` and inference supplied them. Already
    #: **correlated**: the engine samples the whitened ``z`` a
    #: :func:`latent_parameter` declares and
    #: :meth:`GaussianProcessNoise.noise_params` applies
    #: :meth:`GPSolver.latent_transform` to it, so a family reads ``f`` and
    #: needs to know nothing about the kernel (W2.14).
    latent: np.ndarray | None = None
    #: :class:`LimitKind` codes for the retained samples, or ``None``.
    limits: np.ndarray | None = None
    #: Boolean inclusion indicator over the *full* containers (ruled
    #: 2026-09-02, W1.11 gap X-2), so a family carrying its own aligned
    #: per-sample data — a background spectrum, an instrumental template — can
    #: excise it the same way every array already in this record was excised.
    retain: np.ndarray | None = None

    @property
    def correlated(self) -> bool:
        """Whether a covariance beyond the diagonal is in play."""
        return self.kernel is not None

    @property
    def variance(self) -> np.ndarray:
        """``sigma**2``, or a loud failure if this noise model has no sigma."""
        if self.sigma is None:
            raise LikelihoodError(
                "the noise model supplied no per-sample sigma, because the observed container "
                "carries no uncertainties, so the diagonal of K + diag(sigma^2) is undefined. "
                "Attach uncertainties to the observed container, or give the noise model a "
                "'jitter' parameter. (A count-based family such as Poisson defines its own "
                "dispersion and never reaches this path.) "
                "Likelihood.check_alignment() catches this at composition time."
            )
        return self.sigma**2


class NoiseModel(Parameterised, abc.ABC):
    """What the noise is. A :class:`Parameterised`, so its knobs are parameters.

    ``parameters.md`` §13 instructs this contract that GP hyperparameters
    "are ordinary parameters on a ``Parameterised`` noise model, with ``Log``
    bijections", and that is what this class is: nothing about noise-model
    parameters is special, so tying, fixing, priors, plates, serialisation and
    W1.9's lowering all work on them unchanged.

    **A noise model receives the prediction as well as the observation**
    (ruled 2026-09-03, W1.11 gap X-1): :meth:`sigma` and :meth:`noise_params`
    take the *retained* predicted values as a keyword-only ``predicted``
    argument, passed at every call site. A noise whose magnitude depends on
    the model — a fractional model uncertainty, an analytically marginalised
    multiplicative calibration systematic, a model-variance weighting of
    counts — is a ``NoiseModel``, not a family: without the argument the only
    way to express one is to re-implement the sampling distribution, which
    welds noise to family, cannot be reused, and cannot reach the GP path.
    ``predicted`` defaults to ``None`` because a model that does not need it
    (:class:`IndependentNoise`, :class:`GaussianProcessNoise`) simply ignores
    it; the noise model sees the prediction but never the model's *parameters*
    beyond those the likelihood declares.

    **The four capability flags reach here at W2.13** (ruled 2026-09-07,
    ``DEVELOPMENT_PLAN.md`` §2's "Realisation surface (W2.13)" row, fold-in 7;
    ``inference.md`` §10a). Until then only models and instrument steps
    declared them, so a problem could report ``backend="jax"`` while its noise
    model computed in numpy — W2.5 recorded exactly that as a finding. The
    defaults are the reference path's honest answers, the same ones
    :class:`~ampere.core.transform.Model` takes: not differentiable, not
    batchable, on the CPU, supplied by the reference backend. A backend that
    ships a native noise model overrides all four, and
    :attr:`Likelihood.capability_parts` is what carries them onto the composed
    problem.
    """

    #: Whether this model induces correlations between samples.
    CORRELATED: ClassVar[bool] = False

    #: Whether a gradient can be taken through this noise model's evaluation.
    DIFFERENTIABLE: ClassVar[bool] = False
    #: Whether it evaluates a batch of parameter vectors in one call.
    BATCHABLE: ClassVar[bool] = False
    #: Device its arrays live on. Never auto-detected (``architecture.md`` §5).
    DEVICE: ClassVar[str] = "cpu"
    #: Which rung of the capability ladder supplies it (W2.12's fourth flag).
    BACKEND: ClassVar[str] = "reference"

    @abc.abstractmethod
    def sigma(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
        *,
        predicted: np.ndarray | None = None,
    ) -> np.ndarray | None:
        """Per-sample standard deviation on the retained samples.

        ``predicted`` is the retained predicted values — the identical,
        already-excised array the family's ``log_prob`` receives as its first
        argument (float64, or complex128 for a complex family; a noise model
        wanting an amplitude takes ``np.abs(predicted)`` itself). It is
        ``None`` only when no caller holds a prediction; every ampere call
        site passes it.
        """

    def check_compatible(self, family: LikelihoodFamily, observed: FunctionSamples) -> None:
        """Composition-time check. Subclasses extend; this checks uncertainties."""
        if family.REQUIRES_UNCERTAINTY and observed.uncertainty is None:
            raise LikelihoodError(
                f"the {family.NAME} family needs per-sample uncertainties, but the "
                f"{type(observed).__name__} it was given has none. Attach them at construction "
                f"(uncertainty=...), or use a family that models its own dispersion."
            )

    def noise_params(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
        *,
        predicted: np.ndarray | None = None,
        coordinates: np.ndarray | None = None,
        latent: np.ndarray | None = None,
        limits: np.ndarray | None = None,
    ) -> NoiseParams:
        """Assemble the per-evaluation noise description.

        The base implementation describes *uncorrelated* noise. A subclass that
        sets :attr:`CORRELATED` must override it, and is told so rather than
        being allowed to inherit a description that contradicts its own
        declaration — the resulting ``NoiseParams.correlated`` would be
        ``False`` and every family would take the uncorrelated branch.
        """
        if self.CORRELATED:
            raise LikelihoodError(
                f"{type(self).__name__} declares CORRELATED = True but inherits NoiseModel's "
                f"uncorrelated noise_params(), which reports no kernel and no solver. Every "
                f"family would then take its uncorrelated branch and the correlations would "
                f"silently do nothing. Override noise_params() to supply the kernel, solver and "
                f"coordinates (see GaussianProcessNoise)."
            )
        return NoiseParams(
            sigma=self.sigma(observed, retain, values, predicted=predicted),
            values=values,
            coordinates=None,
            kernel=None,
            solver=None,
            latent=latent,
            limits=limits,
            retain=retain,
        )

    def __repr__(self) -> str:
        declared = ", ".join(repr(self.parameters[name]) for name in self.parameters.names)
        return f"{type(self).__name__}({declared})"


def _observed_sigma(observed: FunctionSamples, retain: np.ndarray, owner: str) -> np.ndarray:
    if observed.uncertainty is None:
        raise LikelihoodError(
            f"{owner} needs the observed {type(observed).__name__}'s per-sample uncertainties, "
            f"but it has none. Attach them at construction (uncertainty=...)."
        )
    sigma = _as_float64(np.asarray(observed.uncertainty).ravel()[retain], "observed uncertainties")
    if np.any(sigma <= 0.0):
        raise LikelihoodError(
            f"{owner} was given zero or negative uncertainties on {int(np.sum(sigma <= 0.0))} "
            f"retained sample(s). A zero uncertainty is an infinitely precise measurement, which "
            f"no likelihood can normalise; mask the sample, or give it a real error bar."
        )
    return sigma


class IndependentNoise(NoiseModel):
    """Uncorrelated per-sample noise: the container's own uncertainties.

    The simple path, and it must stay simple —
    ``IndependentNoise()`` declares no parameters at all and does nothing but
    read ``observed.uncertainty``.

    Two optional knobs cover the cases a real fit needs. ``scale`` multiplies
    every uncertainty (the classic "the catalogue's error bars are
    underestimated by a factor" nuisance parameter); ``jitter`` adds a
    standard deviation in quadrature (an unmodelled noise floor). Both are
    ordinary parameters, so either may be fitted, fixed or tied.

    Examples
    --------
    >>> import numpy as np
    >>> import astropy.units as u
    >>> import scipy.stats as st
    >>> from ampere.core import Spectrum
    >>> noise = IndependentNoise()
    >>> noise.parameters.free_size
    0
    >>> inflated = IndependentNoise(scale=st.loguniform(0.5, 5.0))
    >>> inflated.parameters.free_names
    ('scale',)
    """

    CORRELATED: ClassVar[bool] = False

    def __init__(self, *, scale: Any = None, jitter: Any = None) -> None:
        if scale is not None:
            self.register_parameter(_as_hyperparameter("scale", scale, None))
        if jitter is not None:
            self.register_parameter(_as_hyperparameter("jitter", jitter, None))

    def sigma(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
        *,
        predicted: np.ndarray | None = None,
    ) -> np.ndarray | None:
        resolved = self.context({k: v for k, v in values.items() if k in self.parameters})
        if observed.uncertainty is None:
            if "jitter" not in resolved:
                return None
            floor = _positive(resolved["jitter"], "jitter", "IndependentNoise")
            return np.full(int(np.count_nonzero(retain)), floor, dtype=DTYPE)
        sigma = _observed_sigma(observed, retain, "IndependentNoise")
        if "scale" in resolved:
            sigma = sigma * _positive(resolved["scale"], "scale", "IndependentNoise")
        if "jitter" in resolved:
            floor = _positive(resolved["jitter"], "jitter", "IndependentNoise", allow_zero=True)
            sigma = np.sqrt(sigma**2 + floor**2)
        return sigma


class GaussianProcessNoise(NoiseModel):
    """Ampere's flexible likelihood: a GP over the residuals.

    The distinguishing feature of the package. The model's residuals are given
    a covariance ``K(θ) + diag(sigma²)``; structure the physical model cannot
    explain is absorbed by the GP instead of biasing the physical parameters
    (``prior_art.md`` §4: Starfish is the direct scientific ancestor of this
    idea, and ampere generalises and scales it).

    The kernel's hyperparameters become this object's parameters, flatly: the
    ``Parameter`` objects are the kernel's own, shared by name and not by
    identity (``prior_art.md`` Tension 1 — object-identity sharing is the
    gammapy footgun this contract avoids everywhere).

    Parameters
    ----------
    kernel
        The covariance function. :class:`Matern32` unless you have a reason.
    solver
        The strategy that does the algebra. :class:`DenseGP` by default,
        because it is the one that is right rather than the one that is fast.
    scale, jitter
        As :class:`IndependentNoise`: applied to the *diagonal* term before
        the kernel is added.

    Examples
    --------
    >>> import scipy.stats as st
    >>> noise = GaussianProcessNoise(Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0)))
    >>> noise.parameters.free_names
    ('amplitude', 'length_scale')
    >>> noise.solver
    DenseGP(jitter=0.0)
    >>> noise.CORRELATED
    True
    """

    CORRELATED: ClassVar[bool] = True

    def __init__(
        self,
        kernel: Kernel,
        solver: GPSolver | None = None,
        *,
        scale: Any = None,
        jitter: Any = None,
    ) -> None:
        if not isinstance(kernel, Kernel):
            raise LikelihoodError(
                f"GaussianProcessNoise needs a Kernel, got {type(kernel).__name__}. Kernels are "
                f"declared neutrally (family name plus Parameter hyperparameters) so that W1.9 "
                f"can lower them to celerite2 / tinygp / GPyTorch terms."
            )
        self._kernel = kernel
        self._solver = DenseGP() if solver is None else solver
        for parameter in kernel.parameters:
            self.register_parameter(parameter)
        if scale is not None:
            self.register_parameter(_as_hyperparameter("scale", scale, None))
        if jitter is not None:
            self.register_parameter(_as_hyperparameter("jitter", jitter, None))

    @property
    def kernel(self) -> Kernel:
        """The covariance function this noise model applies."""
        return self._kernel

    @property
    def solver(self) -> GPSolver:
        """The strategy that performs the GP algebra."""
        return self._solver

    def sigma(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
        *,
        predicted: np.ndarray | None = None,
    ) -> np.ndarray | None:
        resolved = self.context({k: v for k, v in values.items() if k in self.parameters})
        if observed.uncertainty is None:
            # A latent-GP combination (Poisson and friends) puts the GP on the
            # latent function, not on a residual, so there is no diagonal noise
            # term at all. Families that *do* need one say so through
            # REQUIRES_UNCERTAINTY and are refused at composition time.
            return None
        sigma = _observed_sigma(observed, retain, "GaussianProcessNoise")
        if "scale" in resolved:
            sigma = sigma * _positive(resolved["scale"], "scale", "GaussianProcessNoise")
        if "jitter" in resolved:
            floor = _positive(resolved["jitter"], "jitter", "GaussianProcessNoise", allow_zero=True)
            sigma = np.sqrt(sigma**2 + floor**2)
        return sigma

    def check_compatible(self, family: LikelihoodFamily, observed: FunctionSamples) -> None:
        super().check_compatible(family, observed)
        self._check_circular_solver(family, observed)
        self._solver.check_compatible(self._kernel, observed)
        self._check_hyperparameter_units(observed)

    def _check_circular_solver(self, family: LikelihoodFamily, observed: FunctionSamples) -> None:
        """Refuse, by name, a solver that cannot carry the circular GP's two columns.

        **W4.2.** A circular complex GP is two real processes sharing one
        covariance, so the solver is handed a two-column right-hand side
        (:class:`GPSolver`, :class:`ComplexGaussianFamily`). ``DenseGP``
        declares :attr:`GPSolver.STACKED_RESIDUALS`; the O(N) ``QuasisepGP``
        cannot, and the reason is structural rather than a missing afternoon's
        work: it needs one ordered one-dimensional coordinate
        (:attr:`GPSolver.REQUIRES_ORDERED_1D`), and the coordinates of a
        visibility measurement are a point in the ``(u, v)`` plane at a
        wavelength — there is no ordering of the plane under which a stationary
        kernel becomes a function of one coordinate, so the premise of the
        semiseparable recursion does not hold on this modality at all.

        Checked **before** the solver's own ``check_compatible`` so that this is
        the message a user sees, exactly as the ``Product`` refusal precedes the
        generic quasiseparable one: "select one axis to reach the O(N) path" is
        advice that cannot be taken here, and telling someone to do the
        impossible is worse than telling them nothing.
        """
        if not family.ALLOWS_COMPLEX or np.asarray(observed.values).dtype.kind != "c":
            return
        if self._solver.STACKED_RESIDUALS:
            return
        raise LikelihoodError(
            f"the {family.NAME} family on a complex {type(observed).__name__} is the circular "
            f"complex GP: the real and imaginary parts are two independent real processes "
            f"sharing one covariance, so the solve takes a two-column right-hand side, and "
            f"{self._solver.NAME} declares STACKED_RESIDUALS = False. "
            f"{self._solver.NAME} could not take them even in principle: it needs one ordered "
            f"one-dimensional coordinate, and a visibility lives at a point of the (u, v) plane "
            f"at a wavelength, which no ordering reduces to one coordinate. Use DenseGP, with "
            f'the kernel selecting the axes it acts on -- Matern32(axes=("u", "v")) for an '
            f'isotropic (u, v) kernel, or Product(Matern32(axes=("u", "v")), '
            f'Matern32(axes=("spectral_axis",))) for an error that is smooth in (u, v) and '
            f"sharp in wavelength."
        )

    def _check_hyperparameter_units(self, observed: FunctionSamples) -> None:
        """Delegated to the kernel tree since W4.5.

        The rule and its messages are unchanged; what changed is who knows the
        answer. A composite's hyperparameters are qualified
        (``term0.length_scale``), and a term that selects ``("spectral_axis",)``
        is measured in micron while its sibling on ``("u", "v")`` is
        dimensionless — so "the coordinate axis" is a per-leaf question, and
        :meth:`Kernel.check_units` is where each leaf answers it.
        """
        self._kernel.check_units(observed)

    def noise_params(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
        *,
        predicted: np.ndarray | None = None,
        coordinates: np.ndarray | None = None,
        latent: np.ndarray | None = None,
        limits: np.ndarray | None = None,
    ) -> NoiseParams:
        """The GP's own record, with the **whitening transform applied**.

        ``latent`` arrives whitened — ``latent_parameter`` declares ``z`` with
        an i.i.d. standard-normal prior and says in as many words that "the
        covariance enters through ``f = L(θ) z``, a deterministic transform
        owned by the ``GPSolver``" — and every family reads ``noise.latent``
        as ``f``. This is where the two meet, and it is the only place they
        can meet without a family knowing about kernels: the solver, the
        kernel, the retained coordinates and the resolved hyperparameters are
        all in hand here, and nothing downstream has them.

        Applying it anywhere else was tried and is wrong. Leaving it to the
        family means every latent-consuming family reimplements the transform;
        leaving it to ``Dataset.log_likelihood_of`` means the *dataset* forming
        the retained coordinates, which is ``Likelihood.log_prob``'s excision;
        leaving it out altogether is what this contract did until W2.14, and
        the measurable consequence was a log-likelihood **exactly flat in the
        kernel hyperparameters** — a latent fit that sampled the amplitude and
        the length scale against the prior alone and reported nothing wrong
        (``DEVELOPMENT_PLAN.md`` §2, 2026-09-08; §4.4 clarification).
        """
        bound = self.kernel_for(observed)
        return NoiseParams(
            sigma=self.sigma(observed, retain, values, predicted=predicted),
            values=values,
            coordinates=coordinates,
            kernel=bound,
            solver=self._solver,
            latent=self._realised_latent(bound, coordinates, latent, values),
            limits=limits,
            retain=retain,
        )

    def kernel_for(self, observed: FunctionSamples) -> Kernel:
        """This noise model's kernel, bound to *observed*'s axis order (W4.5).

        The one place a kernel's ``axes=("u", "v")`` becomes "columns 0 and 1
        of the stacked coordinate block". It has to be here and not in the
        solver: the solver is handed a bare ``(n, d)`` array and no longer
        knows what the columns are called, whereas this method has the
        container. Returns the kernel itself, allocating nothing, whenever
        nothing in the tree names axes — which is every problem written before
        W4.5.
        """
        return self._kernel.for_axes([axis.name for axis in observed.axes])

    def _realised_latent(
        self,
        kernel: Kernel,
        coordinates: np.ndarray | None,
        latent: np.ndarray | None,
        values: Mapping[str, Any],
    ) -> np.ndarray | None:
        """``f = L(θ) z`` on the retained coordinates, or ``None``.

        Both arrays are already the **retained** block: ``Likelihood.log_prob``
        excises before it calls this, and ``FittingProblem`` validation refuses
        a latent declaration whose size disagrees with the effective mask, so
        the block's size is an invariant rather than a hope. It is checked
        anyway, because the failure it would otherwise produce is a matrix-shape
        error from inside a solver.

        **W5.4** asks the solver how many whitened values it takes rather than
        assuming one per retained sample: an approximate solver's whitening is
        ``(N, m)`` rather than ``(N, N)``. :meth:`GPSolver.latent_size` answers
        ``n_samples`` for both exact solvers, so nothing about the two exact
        paths changes.
        """
        if latent is None:
            return None
        whitened = _as_float64(latent, "whitened latent GP values")
        if coordinates is None:
            raise LikelihoodError(
                "whitened latent GP values were supplied without the coordinates they live on, "
                "so the whitening transform f = L(theta) z cannot be applied. A correlated "
                "noise model is always handed coordinates by Likelihood.log_prob; a direct "
                "caller of noise_params must pass them too."
            )
        points = _as_points(coordinates, "data coordinates")
        expected = self._solver.latent_size(kernel, int(points.shape[0]))
        if whitened.shape != (expected,):
            per_sample = expected == points.shape[0]
            raise LikelihoodError(
                f"the whitened latent GP values have shape {whitened.shape} but "
                f"{self._solver.NAME} whitens {expected} value(s) for "
                f"{points.shape[0]} retained sample(s)."
                + (
                    " One latent value per retained sample."
                    if per_sample
                    else " This solver's latent block is its basis size, not the sample count "
                    "(GPSolver.latent_size)."
                )
            )
        return self._solver.latent_transform(kernel, points, whitened, values)


# ---------------------------------------------------------------------------
# Joint noise over a tuple of channels: the shared-grid intrinsic
# coregionalisation model (W5.9)
# ---------------------------------------------------------------------------


class ChannelCoupling(Parameterised, abc.ABC):
    """The ``T x T`` positive-definite matrix ``B`` of an intrinsic coregionalisation model.

    ``K = B ⊗ K_x`` is the covariance of ``T`` channels observed on one shared
    grid: ``K_x`` says how a channel correlates with itself along the grid, and
    ``B`` says how the channels correlate with each other. This class is ``B``,
    and it is a :class:`~ampere.core.parameter.Parameterised` for the reason
    every other knob in this contract is — tying, fixing, priors and W1.9's
    lowering then work on its parameters unchanged.

    Two things it must supply, and the second is why it is a class rather than
    a matrix. :meth:`matrix` is the declaration; :meth:`eigen` is the
    **eigendecomposition**, which is what makes the joint solve exact and
    O(N) rather than O((NT)³): rotating the ``T`` residual vectors by ``Qᵀ``
    decouples them into ``T`` scalar GPs sharing ``K_x``
    (:class:`JointGaussianProcessNoise`). A parameterisation that knows its own
    eigenvectors in closed form — :class:`RotationCoupling` does — hands them
    over without a numerical ``eigh`` at every evaluation, which matters
    because ``eigh``'s gradient is ill-conditioned exactly where two
    eigenvalues coincide, and "the two channels have equal variance" is a
    perfectly ordinary point of the posterior.

    ``xp`` is the array namespace the arithmetic happens in: :mod:`numpy` on
    the reference path, ``jax.numpy`` or :mod:`torch` when a backend lowers the
    same declaration for NUTS. Only ``cos``, ``sin``, ``exp``, ``stack`` and
    (for :class:`CholeskyCoupling`) ``linalg.eigh`` are used, which all three
    spell identically — so there is **one** implementation of each
    parameterisation and no per-backend transcription to drift.
    """

    #: Neutral name, for provenance and error messages.
    NAME: ClassVar[str] = ""

    def __init__(self, channels: int) -> None:
        count = int(channels)
        if count < 2:
            raise LikelihoodError(
                f"a channel coupling describes how two or more channels correlate, got "
                f"{count}. One channel is an ordinary GaussianProcessNoise."
            )
        self._channels = count

    @property
    def channels(self) -> int:
        """``T`` — how many channels this coupling spans."""
        return self._channels

    def resolved(self, values: Mapping[str, Any]) -> dict[str, Any]:
        """This coupling's own parameters, resolved out of a likelihood-wide mapping.

        The same filter-then-``context`` idiom :class:`IndependentNoise` and
        :class:`GaussianProcessNoise` use: the mapping handed down carries the
        kernel's hyperparameters too, and completing a declaration against
        names it does not own is exactly what the idiom exists to avoid.
        """
        return self.context({k: v for k, v in values.items() if k in self.parameters})

    @abc.abstractmethod
    def eigen(self, resolved: Mapping[str, Any], *, xp: Any = np) -> tuple[Any, Any]:
        """``(eigenvalues, Q)`` of ``B``, with ``B = Q diag(eigenvalues) Qᵀ``.

        The eigenvalues are strictly positive (``B`` is positive definite) and
        ``Q`` is orthogonal. Both are ``xp`` arrays of shape ``(T,)`` and
        ``(T, T)``; column ``s`` of ``Q`` is the eigenvector of eigenvalue
        ``s``.
        """

    def matrix(self, resolved: Mapping[str, Any], *, xp: Any = np) -> Any:
        """``B`` itself, reassembled from :meth:`eigen`.

        Used by the conformance battery — which materialises ``B ⊗ K_x`` and
        scores it densely as the definition of the right answer — and by
        anybody who wants to read the fitted cross-covariance off a posterior
        draw. Never on the hot path: the solve uses the eigendecomposition and
        never forms ``B``.
        """
        eigenvalues, rotation = self.eigen(resolved, xp=xp)
        return (rotation * eigenvalues) @ rotation.T

    def __repr__(self) -> str:
        declared = ", ".join(repr(self.parameters[name]) for name in self.parameters.names)
        return f"{type(self).__name__}({declared})"


class RotationCoupling(ChannelCoupling):
    """``T = 2``: a rotation angle and two log-variances — the physical parameterisation.

    ``B = Q(θ) diag(exp(v₀), exp(v₁)) Q(θ)ᵀ`` with ``Q(θ)`` the plane rotation.
    Three parameters, which is the full ``T(T+1)/2`` of a 2x2 positive-definite
    matrix — so nothing is given up — but expressed in the coordinates the
    *physics* is stated in rather than as three entries of a matrix that must
    then be checked for positive-definiteness. ``DEVELOPMENT_PLAN.md`` §5's
    Phase 5 bullet asks for exactly this: "``B`` parameterised physically (a
    rotation for Q/U leakage) rather than a free ``T(T+1)/2``".

    The two cases this ships for read directly off the parameters. For
    **astrometry** (W5.9's first customer) the angle is the position angle of a
    centroiding systematic's major axis on the sky and the two variances are
    its semi-axes: an error that is elongated along one direction and shared by
    the ``ra`` and ``dec`` channels is ``θ`` away from the axes and is
    invisible to two independent GPs, which can only inflate each axis
    separately. For **polarimetry** the angle is the Q/U leakage angle and the
    rotation is literally the one the instrument applies.

    ``θ`` is identified on ``[0, π)``: adding ``π`` is the same matrix, and
    adding ``π/2`` is the same matrix with the two variances exchanged — so a
    prior wider than ``[0, π)`` samples a label-switched copy of the same
    posterior. ``scipy.stats.uniform(0.0, np.pi)`` is the recommended prior;
    this class does not impose one, because a problem with a known instrumental
    angle should fix it.

    The **log**-variances, rather than variances with a positive support, are
    what a NUTS chain wants: the eigenvalues of a coupling matrix span orders
    of magnitude and a positivity constraint handled by a bijection on a log
    scale is the parameterisation ``lowering.md`` §3 recommends everywhere else
    for the same reason. They are variances in the observed values' **squared
    unit**, and they carry the joint model's overall scale — see
    :class:`JointGaussianProcessNoise` on why the kernel's own amplitude must
    not also be free.

    Parameters
    ----------
    angle
        ``θ``, radians. A frozen ``scipy.stats`` distribution to fit it, a
        number to hold it fixed.
    log_variance_0, log_variance_1
        The natural logarithms of ``B``'s two eigenvalues.

    Examples
    --------
    >>> import numpy as np, scipy.stats as st
    >>> coupling = RotationCoupling(
    ...     st.uniform(0.0, np.pi), st.norm(-7.0, 2.0), st.norm(-7.0, 2.0)
    ... )
    >>> coupling.parameters.free_names
    ('angle', 'log_variance_0', 'log_variance_1')
    >>> at = {"angle": 0.0, "log_variance_0": 0.0, "log_variance_1": np.log(4.0)}
    >>> bool(np.allclose(coupling.matrix(at), np.diag([1.0, 4.0])))
    True
    """

    NAME: ClassVar[str] = "rotation"

    def __init__(self, angle: Any, log_variance_0: Any, log_variance_1: Any) -> None:
        super().__init__(2)
        self.register_parameter(_as_hyperparameter("angle", angle, None))
        self.register_parameter(_as_hyperparameter("log_variance_0", log_variance_0, None))
        self.register_parameter(_as_hyperparameter("log_variance_1", log_variance_1, None))

    def eigen(self, resolved: Mapping[str, Any], *, xp: Any = np) -> tuple[Any, Any]:
        angle = resolved["angle"]
        cosine = xp.cos(angle)
        sine = xp.sin(angle)
        rotation = xp.stack([xp.stack([cosine, -sine]), xp.stack([sine, cosine])])
        variances = xp.stack([resolved["log_variance_0"], resolved["log_variance_1"]])
        return xp.exp(variances), rotation


class CholeskyCoupling(ChannelCoupling):
    """The general ``T``: ``B = L Lᵀ`` with ``L`` lower-triangular, positive diagonal.

    The fallback :class:`RotationCoupling` is the special case of. Any
    positive-definite ``B`` is ``L Lᵀ`` for exactly one such ``L``, so this
    parameterisation is complete, unconstrained (the ``T`` diagonal entries are
    fitted as logarithms and the ``T(T-1)/2`` strictly-lower ones are free on
    the whole line) and needs no positive-definiteness check at all — which is
    the whole reason to prefer it to fitting ``B``'s entries.

    What it gives up is the closed-form eigenvectors: :meth:`eigen` calls
    ``xp.linalg.eigh``, whose gradient is ill-conditioned where two eigenvalues
    coincide. That is not a problem for a ``B`` whose channels genuinely differ
    and it is a real one for a ``B`` near a multiple of the identity, so for
    ``T = 2`` prefer :class:`RotationCoupling`, which never forms an ``eigh``
    at all.

    Priors. ``off_diagonal`` entries are recommended a ``norm(0, s)`` whose
    scale is comparable with ``exp(log_diagonal / 2)``, which makes the implied
    prior on the channel *correlations* roughly flat rather than piled at
    ±1; ``log_diagonal`` entries take the same log-scale prior
    :class:`RotationCoupling`'s log-variances do. This class imposes neither.

    Parameters
    ----------
    channels
        ``T``.
    log_diagonal
        ``T`` priors or numbers: the logarithms of ``L``'s diagonal.
    off_diagonal
        ``T(T-1)/2`` priors or numbers, in row-major order over the strictly
        lower triangle (``(1,0), (2,0), (2,1), (3,0), ...``).

    Examples
    --------
    >>> import numpy as np
    >>> coupling = CholeskyCoupling(2, log_diagonal=[0.0, 0.0], off_diagonal=[0.5])
    >>> coupling.parameters.names
    ('log_diagonal_0', 'log_diagonal_1', 'off_diagonal_0')
    >>> resolved = coupling.resolved({})
    >>> np.allclose(coupling.matrix(resolved), np.array([[1.0, 0.5], [0.5, 1.25]]))
    True
    """

    NAME: ClassVar[str] = "cholesky"

    def __init__(
        self,
        channels: int,
        *,
        log_diagonal: Sequence[Any],
        off_diagonal: Sequence[Any] = (),
    ) -> None:
        super().__init__(channels)
        count = self.channels
        below = count * (count - 1) // 2
        if len(log_diagonal) != count:
            raise LikelihoodError(
                f"a CholeskyCoupling over {count} channels needs {count} log_diagonal "
                f"entries, got {len(log_diagonal)}."
            )
        if len(off_diagonal) != below:
            raise LikelihoodError(
                f"a CholeskyCoupling over {count} channels needs {below} off_diagonal "
                f"entries — the strictly lower triangle of L, row-major — got "
                f"{len(off_diagonal)}."
            )
        for index, declaration in enumerate(log_diagonal):
            self.register_parameter(_as_hyperparameter(f"log_diagonal_{index}", declaration, None))
        for index, declaration in enumerate(off_diagonal):
            self.register_parameter(_as_hyperparameter(f"off_diagonal_{index}", declaration, None))

    def eigen(self, resolved: Mapping[str, Any], *, xp: Any = np) -> tuple[Any, Any]:
        count = self.channels
        # A zero of whatever type the values are: a Python ``0.0`` in a
        # ``torch.stack`` is a TypeError and in a traced jax computation is a
        # constant of the wrong dtype, and this coupling has exactly one
        # implementation for all three namespaces.
        zero = 0.0 * resolved["log_diagonal_0"]
        rows = []
        below = 0
        for row in range(count):
            entries = []
            for column in range(count):
                if column < row:
                    entries.append(resolved[f"off_diagonal_{below + column}"])
                elif column == row:
                    entries.append(xp.exp(resolved[f"log_diagonal_{row}"]))
                else:
                    entries.append(zero)
            rows.append(xp.stack(entries))
            below += row
        lower = xp.stack(rows)
        eigenvalues, rotation = xp.linalg.eigh(lower @ lower.T)
        return eigenvalues, rotation


def _amplitude_names(kernel: Kernel) -> tuple[str, ...]:
    """Which of *kernel*'s parameters carry the covariance's overall scale.

    ``Kernel.VALUE_SCALED`` is the per-leaf declaration of which
    hyperparameters are measured in the observed values' unit — the amplitude,
    for every kernel that ships — and a composite qualifies its children's
    names (``term0.amplitude``), so the test is on the last segment.
    """
    scaled = {name for leaf in (*kernel.leaves(), kernel) for name in leaf.VALUE_SCALED}
    return tuple(name for name in kernel.parameters.names if name.rsplit(".", 1)[-1] in scaled)


def _stack_channels(arrays: Sequence[Any], what: str) -> np.ndarray:
    """``T`` equal-length vectors as one ``(n, T)`` block, float64."""
    columns = [_as_float64(np.asarray(array).ravel(), what) for array in arrays]
    sizes = {column.size for column in columns}
    if len(sizes) != 1:
        raise LikelihoodError(
            f"the {what} of a joint noise model's channels have different lengths "
            f"({sorted(sizes)}). The channels share one grid by declaration — that is what "
            f"makes B ⊗ K_x exact — so they have the same number of retained samples."
        )
    return np.ascontiguousarray(np.column_stack(columns))


def _channel_block(residuals: Any, channels: int) -> np.ndarray:
    """A sequence of ``T`` vectors, or an ``(n, T)`` block, as an ``(n, T)`` block."""
    array = np.asarray(residuals)
    if array.ndim == 2 and array.shape[1] == channels:
        return _as_float64(array, "residuals")
    return _stack_channels(list(residuals), "residuals")


class JointGaussianProcessNoise(NoiseModel):
    """One correlated process over ``T`` channels of one model on a shared grid.

    ``likelihoods.md`` §15's eleventh limitation, lifted for the case it is
    exact in. A noise model belongs to one :class:`Likelihood`, which belongs
    to one dataset, so a single correlated process over ``T``-vectors — the
    ``(RA, Dec)`` residuals of an astrometric solution perturbed together by a
    centroiding systematic, Stokes ``Q``/``U`` mixed by an instrumental
    leakage, a calibration error shared across bands — had no expression, and
    "two scalar GPs with tied hyperparameters" is a *different model*: it can
    inflate each channel but it cannot express the cross-covariance, which is
    the whole of the systematic.

    This class is that process, scoped to the case ``DEVELOPMENT_PLAN.md`` §5
    scopes it to: ``K = B ⊗ K_x``, the **shared-grid intrinsic**
    coregionalisation model, with ``B`` a ``TxT`` positive-definite
    :class:`ChannelCoupling` and ``K_x`` one ordinary :class:`Kernel` on the
    grid every channel shares.

    **Why it is still O(N).** Diagonalise ``B = Q Λ Qᵀ`` and rotate the
    residuals by ``Qᵀ``. Because ``(Qᵀ ⊗ I)(B ⊗ K_x)(Q ⊗ I) = Λ ⊗ K_x`` is
    block-diagonal, the ``T`` rotated residual vectors are **independent**, and
    each is scored by an ordinary scalar GP with covariance ``λ_s K_x +
    diag(sigma^2)`` — so the joint density is ``T`` calls to the bound
    :class:`GPSolver`, with ``QuasisepGP`` on an ordered one-dimensional grid
    giving ``T·O(N)`` rather than the ``O((NT)³)`` of the materialised matrix.
    Nothing is approximated: :meth:`log_prob` and a dense Cholesky of
    ``B ⊗ K_x + diag(sigma^2)`` agree to the solver tolerance, which is a
    conformance row.

    Each scalar solve is handed the **rescaled** problem ``u = r̃_s / √λ_s``
    against ``K_x + diag(sigma^2/λ_s)``, with ``-(N/2) log λ_s`` added back, rather
    than a kernel whose amplitude has been multiplied by ``√λ_s``. The two are
    the same number; the first works for *any* kernel — a ``Sum``, a kernel
    with no amplitude at all, a user's own — whereas the second would have to
    reach into a kernel tree and rescale it, which is neither general nor
    something a kernel's declaration promises.

    **The one restriction, and it is a restriction rather than an oversight.**
    The rotation only leaves the *diagonal* noise term diagonal when the
    channels share one per-sample variance: the rotated ``(s, s')`` block of
    ``diag(sigma_t²)`` is ``Σ_t Q_{ts} Q_{ts'} diag(sigma_t²)``, which is
    ``δ_{ss'} diag(sigma^2)`` when every ``sigma_t`` is the same vector and is a full
    coupling otherwise. So the channels must carry **equal uncertainties**
    (heteroscedastic along the grid as much as you like — it is the *channels*
    that must agree, not the samples), which is exactly what a shared-grid
    astrometric solution or a Stokes ``Q``/``U`` pair from one polarimeter
    produces. Unequal per-channel errors lose the Kronecker structure
    altogether and belong with the general LMC and mismatched grids on the
    dense/reduced-rank follow-on (``likelihoods.md`` §15). This is checked at
    composition, by name.

    **Where it lives.** Not in a :class:`Likelihood` — a likelihood scores one
    dataset and this scores ``T`` of them. It is declared on the
    :class:`~ampere.core.dataset.DatasetCollection` instead::

        DatasetCollection(
            {"ra": Dataset(...), "dec": Dataset(...)},
            joint={"astrom": JointGaussianProcessNoise(kernel, QuasisepGP(),
                                                       datasets=("ra", "dec"),
                                                       coupling=RotationCoupling(...))},
        )

    which is the first use of
    :meth:`~ampere.core.dataset.DatasetCollection.contributions` as something
    other than a sum: the group contributes **one** term, keyed by the group's
    own label, in place of its members' separate ones (``inference.md`` §4,
    the ``"joint"`` decomposition).

    Parameters
    ----------
    kernel
        ``K_x``, the covariance along the shared grid. Its hyperparameters
        become this noise model's own, as for
        :class:`GaussianProcessNoise`. Its overall **amplitude must not be
        free**: ``B ⊗ (a² K̃) = (a² B) ⊗ K̃``, so a free amplitude and a free
        ``B`` are exactly degenerate and the posterior has a ridge rather than
        a mode. Fix it — ``Matern32(1.0, length_scale_prior)`` — and let ``B``
        carry the scale, which is what ``B``'s log-variances are for. Refused
        at construction.
    solver
        The strategy each rotated scalar solve uses. ``DenseGP()`` by default,
        because that is the one that is right rather than the one that is fast;
        pass ``QuasisepGP()`` on an ordered one-dimensional grid — a
        ``TimeSeries``'s ``time`` axis — to get the O(N) path. The bound solver
        decides, once, for every rotated output: they all live on the same
        grid and see the same kernel, so a per-output choice could only differ
        by accident.
    datasets
        The labels of the datasets this process spans, **in channel order**:
        ``B``'s row ``t`` is ``datasets[t]``. Two or more, distinct, and the
        same count as the coupling's ``T``.
    coupling
        ``B``. :class:`RotationCoupling` for ``T = 2``,
        :class:`CholeskyCoupling` in general.
    scale, jitter
        As :class:`IndependentNoise`, applied to the channels' shared diagonal
        before ``B ⊗ K_x`` is added. One pair for the whole group, not one per
        channel: the exactness condition above is that the group has **one**
        diagonal, so a per-channel scale would be a per-channel diagonal and
        would take the group off the exact path.

    Examples
    --------
    >>> import numpy as np, scipy.stats as st
    >>> noise = JointGaussianProcessNoise(
    ...     Matern32(1.0, st.loguniform(1.0, 1e3), axes=("time",)),
    ...     QuasisepGP(),
    ...     datasets=("ra", "dec"),
    ...     coupling=RotationCoupling(
    ...         st.uniform(0.0, np.pi), st.norm(-7.0, 2.0), st.norm(-7.0, 2.0)
    ...     ),
    ... )
    >>> noise.datasets
    ('ra', 'dec')
    >>> noise.parameters.free_names
    ('length_scale', 'angle', 'log_variance_0', 'log_variance_1')
    >>> noise.JOINT, noise.CORRELATED
    (True, True)
    """

    CORRELATED: ClassVar[bool] = True
    #: This noise model scores several datasets at once. The flag is what
    #: :class:`Likelihood` refuses on and what
    #: :class:`~ampere.core.dataset.DatasetCollection` routes on; it is a
    #: declaration rather than an ``isinstance`` so that a user's own joint
    #: noise model composes without subclassing this one.
    JOINT: ClassVar[bool] = True

    def __init__(
        self,
        kernel: Kernel,
        solver: GPSolver | None = None,
        *,
        datasets: Sequence[str],
        coupling: ChannelCoupling,
        scale: Any = None,
        jitter: Any = None,
    ) -> None:
        if not isinstance(kernel, Kernel):
            raise LikelihoodError(
                f"JointGaussianProcessNoise needs a Kernel for K_x, got {type(kernel).__name__}."
            )
        if not isinstance(coupling, ChannelCoupling):
            raise LikelihoodError(
                f"JointGaussianProcessNoise needs a ChannelCoupling for B, got "
                f"{type(coupling).__name__}. Use RotationCoupling for two channels, "
                f"CholeskyCoupling in general."
            )
        labels = tuple(str(label) for label in datasets)
        if len(labels) < 2:
            raise LikelihoodError(
                f"a joint noise model spans two or more datasets, got {list(labels)}. One "
                f"channel is an ordinary GaussianProcessNoise."
            )
        if len(set(labels)) != len(labels):
            raise LikelihoodError(
                f"the datasets a joint noise model spans must be distinct, got {list(labels)}. "
                f"B's row t is datasets[t], so a repeated label is a channel correlated with "
                f"itself through two different rows."
            )
        if len(labels) != coupling.channels:
            raise LikelihoodError(
                f"the coupling B is {coupling.channels}x{coupling.channels} but the joint noise "
                f"model spans {len(labels)} datasets ({list(labels)})."
            )
        self._kernel = kernel
        self._solver = DenseGP() if solver is None else solver
        self._coupling = coupling
        self._datasets = labels
        for parameter in kernel.parameters:
            self.register_parameter(parameter)
        for parameter in coupling.parameters:
            if parameter.name in self.parameters:
                raise LikelihoodError(
                    f"the kernel and the coupling both declare a parameter named "
                    f"{parameter.name!r}. A noise model holds one flat namespace "
                    f"(parameters.md §12.4), so rename one of them."
                )
            self.register_parameter(parameter)
        if scale is not None:
            self.register_parameter(_as_hyperparameter("scale", scale, None))
        if jitter is not None:
            self.register_parameter(_as_hyperparameter("jitter", jitter, None))
        self._check_identifiable()

    # -- declarations --------------------------------------------------------

    @property
    def kernel(self) -> Kernel:
        """``K_x``: the covariance along the grid every channel shares."""
        return self._kernel

    @property
    def solver(self) -> GPSolver:
        """The strategy each rotated scalar solve goes through."""
        return self._solver

    @property
    def coupling(self) -> ChannelCoupling:
        """``B``: how the channels correlate with one another."""
        return self._coupling

    @property
    def datasets(self) -> tuple[str, ...]:
        """The dataset labels this process spans, in channel order."""
        return self._datasets

    @property
    def channels(self) -> int:
        """``T``."""
        return len(self._datasets)

    def _check_identifiable(self) -> None:
        """Refuse a free kernel amplitude beside a free ``B``.

        ``B ⊗ (a² K̃) = (a² B) ⊗ K̃`` exactly, so the overall amplitude of the
        kernel and the overall scale of the coupling are the same number
        written twice. Fitting both gives a posterior with a ridge along
        ``(a², B) -> (c a², B/c)``: a chain wanders along it, ``R̂`` never
        settles, and the marginal on either looks like its prior. One
        *fixed* amplitude anywhere in the tree pins the scale, which is why
        a ``Sum`` with one term's amplitude fixed and another's free is
        accepted — the relative amplitudes of a sum are identified.
        """
        amplitudes = _amplitude_names(self._kernel)
        if not amplitudes:
            return
        free = set(self._kernel.parameters.free_names)
        if not all(name in free for name in amplitudes):
            return
        if self._coupling.parameters.free_size == 0:
            return
        raise LikelihoodError(
            f"JointGaussianProcessNoise was given a kernel whose amplitude(s) "
            f"{list(amplitudes)} are all free, beside a coupling B with "
            f"{self._coupling.parameters.free_size} free parameter(s). Those are the same "
            f"degree of freedom twice: B ⊗ (a² K̃) = (a² B) ⊗ K̃, so the posterior has a ridge "
            f"along (a², B) -> (c a², B/c) rather than a mode, and both marginals come back "
            f"looking like their priors. Fix the kernel's amplitude — Matern32(1.0, "
            f"length_scale_prior) — and let B's log-variances carry the scale, which is what "
            f"they are for. (Fixing one term's amplitude in a Sum is enough: the relative "
            f"amplitudes of a sum are identified.)"
        )

    # -- the diagonal --------------------------------------------------------

    def sigma(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
        *,
        predicted: np.ndarray | None = None,
    ) -> np.ndarray | None:
        """The group's shared per-sample standard deviation, on one channel's container.

        Every channel returns the same vector by declaration — that is what
        :meth:`check_group` enforces — so which container it is asked of does
        not matter, and the group's own ``scale`` and ``jitter`` are applied
        here exactly as :class:`IndependentNoise` applies its own.
        """
        resolved = self.context({k: v for k, v in values.items() if k in self.parameters})
        if observed.uncertainty is None:
            if "jitter" not in resolved:
                return None
            floor = _positive(resolved["jitter"], "jitter", "JointGaussianProcessNoise")
            return np.full(int(np.count_nonzero(retain)), floor, dtype=DTYPE)
        sigma = _observed_sigma(observed, retain, "JointGaussianProcessNoise")
        if "scale" in resolved:
            sigma = sigma * _positive(resolved["scale"], "scale", "JointGaussianProcessNoise")
        if "jitter" in resolved:
            floor = _positive(
                resolved["jitter"], "jitter", "JointGaussianProcessNoise", allow_zero=True
            )
            sigma = np.sqrt(sigma**2 + floor**2)
        return sigma

    def noise_params(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
        *,
        predicted: np.ndarray | None = None,
        coordinates: np.ndarray | None = None,
        latent: np.ndarray | None = None,
        limits: np.ndarray | None = None,
    ) -> NoiseParams:
        """Refused: this noise model has no per-dataset record to give.

        A :class:`NoiseParams` describes the covariance of *one* container's
        samples, and this model's covariance couples ``T`` containers. Scoring
        one channel out of it would be scoring ``B_{tt} K_x``, which is the
        marginal of the joint process and is *not* the joint density — the
        cross-covariance is precisely what it drops, and the cross-covariance
        is the whole model.
        """
        raise LikelihoodError(
            f"JointGaussianProcessNoise scores {self.channels} channels "
            f"({list(self._datasets)}) jointly and has no single-dataset NoiseParams to give: "
            f"the cross-covariance between the channels is the model, and a per-dataset record "
            f"would silently drop it. Declare it on the DatasetCollection — "
            f"DatasetCollection({{...}}, joint={{'label': noise}}) — which routes the group "
            f"through log_prob() instead of through Likelihood.log_prob."
        )

    # -- composition-time checks ---------------------------------------------

    def check_compatible(self, family: LikelihoodFamily, observed: FunctionSamples) -> None:
        """Refused for the same reason :meth:`noise_params` is."""
        raise LikelihoodError(
            f"JointGaussianProcessNoise cannot be the noise model of a single Likelihood: it "
            f"spans {self.channels} datasets ({list(self._datasets)}). Pass it to the "
            f"DatasetCollection as joint={{'label': noise}}, and leave each member dataset's "
            f"own likelihood as Likelihood(GaussianFamily())."
        )

    def kernel_for(self, observed: FunctionSamples) -> Kernel:
        """``K_x`` bound to *observed*'s axis order — :meth:`GaussianProcessNoise.kernel_for`."""
        return self._kernel.for_axes([axis.name for axis in observed.axes])

    def check_group(
        self,
        observed: Mapping[str, FunctionSamples],
        likelihoods: Mapping[str, Likelihood],
        *,
        group: str,
    ) -> None:
        """Every composition-time rule this model has, checked once, loudly.

        Called by :class:`~ampere.core.dataset.DatasetCollection` at
        construction with the member datasets' observed containers and
        likelihoods, in this model's own channel order.
        """
        missing = [label for label in self._datasets if label not in observed]
        if missing:
            raise LikelihoodError(
                f"joint noise group {group!r} spans datasets {list(self._datasets)}, but "
                f"{missing} are not in this collection."
            )
        containers = [observed[label] for label in self._datasets]
        reference = containers[0]
        self._check_shared_grid(containers, group=group)
        for label in self._datasets:
            self._check_member_likelihood(label, likelihoods[label], group=group)
        self._check_shared_diagonal(containers, group=group)
        bound = self.kernel_for(reference)
        self._solver.check_compatible(bound, reference)
        bound.check_units(reference)

    def _check_shared_grid(self, containers: Sequence[FunctionSamples], *, group: str) -> None:
        """One grid, not ``T`` grids that happen to agree to nine decimal places.

        Exact equality, the same rule ``Likelihood.check_alignment`` applies
        between a predicted and an observed container and for the same reason
        (``transformations.md`` §10, W1.11 gap I-2): coordinates recomputed
        from first principles differ in their last bits, and a covariance built
        on one grid and applied to another is wrong in a way no tolerance
        describes. Mismatched grids are the general LMC's problem, not this
        model's (``likelihoods.md`` §15).
        """
        reference = containers[0]
        if reference.LAYOUT is not Layout.POINTS:
            raise LikelihoodError(
                f"joint noise group {group!r} is bound to a "
                f"{type(reference).__name__}, whose layout is {reference.LAYOUT.value}. The "
                f"shared-grid intrinsic model works on point-set containers, as every GP solver "
                f"in this contract does."
            )
        names = tuple(axis.name for axis in reference.axes)
        for label, container in zip(self._datasets[1:], containers[1:], strict=True):
            other = tuple(axis.name for axis in container.axes)
            if other != names:
                raise LikelihoodError(
                    f"joint noise group {group!r}: dataset {self._datasets[0]!r} has axes "
                    f"{list(names)} and dataset {label!r} has {list(other)}. B ⊗ K_x is a "
                    f"covariance over one grid shared by every channel."
                )
            if container.shape != reference.shape:
                raise LikelihoodError(
                    f"joint noise group {group!r}: dataset {self._datasets[0]!r} holds "
                    f"{reference.n_samples} sample(s) and dataset {label!r} holds "
                    f"{container.n_samples}. The channels share one grid."
                )
            for axis, reference_axis in zip(container.axes, reference.axes, strict=True):
                if not np.array_equal(
                    np.asarray(axis.values, dtype=DTYPE),
                    np.asarray(reference_axis.values, dtype=DTYPE),
                ):
                    raise LikelihoodError(
                        f"joint noise group {group!r}: the {axis.name!r} axis of dataset "
                        f"{label!r} is not identical to dataset {self._datasets[0]!r}'s. The "
                        f"shared-grid model needs one grid, compared exactly — build every "
                        f"channel's container from the same coordinate array rather than from "
                        f"two computations that agree to nine decimal places. Genuinely "
                        f"different grids are the general (non-Kronecker) LMC, recorded as a "
                        f"follow-on in likelihoods.md §15."
                    )
            if not _masks_agree(container, reference):
                raise LikelihoodError(
                    f"joint noise group {group!r}: dataset {label!r} and dataset "
                    f"{self._datasets[0]!r} mask different samples. A rotation mixes the "
                    f"channels sample by sample, so a sample masked in one channel and not in "
                    f"another has no rotated value at all; mask it in every channel, or in "
                    f"none."
                )

    def _check_member_likelihood(self, label: str, likelihood: Likelihood, *, group: str) -> None:
        """A member dataset supplies the family; the group supplies all the noise."""
        family = likelihood.family
        if not family.ANALYTIC_WITH_GP or not family.GP_ANALYTIC_IMPLEMENTED:
            raise LikelihoodError(
                f"joint noise group {group!r}: dataset {label!r} uses the {family.NAME!r} "
                f"family, whose GP marginalisation is not the analytic one this model needs. "
                f"B ⊗ K_x is marginalised in closed form under a Gaussian family and under no "
                f"other; a latent formulation over T coupled channels is not implemented."
            )
        if likelihood.censoring is not None:
            raise LikelihoodError(
                f"joint noise group {group!r}: dataset {label!r} carries a censoring "
                f"declaration. A limit is a statement about one sample's own sampling "
                f"distribution, and this model's samples are not independent across channels; "
                f"the composition is refused rather than silently scored as if they were."
            )
        noise = likelihood.noise
        if noise.CORRELATED:
            raise LikelihoodError(
                f"joint noise group {group!r}: dataset {label!r} also carries its own "
                f"correlated noise model ({type(noise).__name__}). The group *is* the "
                f"correlated process; a second one on the same residual would be counted twice. "
                f"Leave the member's likelihood as Likelihood(GaussianFamily())."
            )
        if len(noise.parameters.names) != 0:
            raise LikelihoodError(
                f"joint noise group {group!r}: dataset {label!r}'s own noise model "
                f"({type(noise).__name__}) declares parameters "
                f"{list(noise.parameters.names)}. The group owns the whole covariance, "
                f"diagonal included, so a per-channel scale or jitter would give the channels "
                f"different diagonals — and the rotation is exact only where they share one "
                f"(see this class's docstring). Declare scale=/jitter= on the joint noise "
                f"model instead, where they apply to the group."
            )

    def _check_shared_diagonal(self, containers: Sequence[FunctionSamples], *, group: str) -> None:
        """The channels' uncertainties must agree, because the rotation says so."""
        reference = containers[0]
        if reference.uncertainty is None:
            return
        first = np.asarray(reference.uncertainty, dtype=DTYPE).ravel()
        for label, container in zip(self._datasets[1:], containers[1:], strict=True):
            if container.uncertainty is None or not np.array_equal(
                np.asarray(container.uncertainty, dtype=DTYPE).ravel(), first
            ):
                raise LikelihoodError(
                    f"joint noise group {group!r}: dataset {label!r} and dataset "
                    f"{self._datasets[0]!r} carry different per-sample uncertainties. The "
                    f"rotation Q^T (x) I leaves diag(sigma^2) diagonal only when every channel "
                    f"has the same sigma vector; otherwise the rotated noise couples the outputs "
                    f"again and "
                    f"the T scalar solves are not the joint density. Heteroscedasticity *along* "
                    f"the grid is fine; it is the channels that must agree. Unequal per-channel "
                    f"errors need the dense or reduced-rank solver, recorded with the general "
                    f"LMC as a follow-on in likelihoods.md §15."
                )

    # -- evaluation ----------------------------------------------------------

    def eigen(self, values: Mapping[str, Any], *, xp: Any = np) -> tuple[Any, Any]:
        """``B``'s eigenvalues and eigenvectors at *values* — :meth:`ChannelCoupling.eigen`."""
        return self._coupling.eigen(self._coupling.resolved(values), xp=xp)

    def coupling_matrix(self, values: Mapping[str, Any], *, xp: Any = np) -> Any:
        """``B`` itself at *values*. Never on the hot path; see :meth:`ChannelCoupling.matrix`."""
        return self._coupling.matrix(self._coupling.resolved(values), xp=xp)

    def rotate(self, residuals: Any, values: Mapping[str, Any]) -> np.ndarray:
        """``R Q``: the ``T`` residual columns in ``B``'s eigenbasis.

        Column ``s`` of the result is ``Σ_t Q_{ts} r_t``, which is the ``s``-th
        block of ``(Qᵀ ⊗ I) r`` for a channel-major stacking — the rotated
        output the diagnostics and the pointwise group are indexed by.
        """
        block = _channel_block(residuals, self.channels)
        _, rotation = self.eigen(values)
        return np.ascontiguousarray(block @ np.asarray(rotation, dtype=DTYPE))

    def _rotated_problems(
        self, residuals: Any, variance: Any, values: Mapping[str, Any]
    ) -> tuple[list[tuple[np.ndarray, np.ndarray]], np.ndarray]:
        """The ``T`` decoupled scalar problems, rescaled to a unit eigenvalue.

        Returns ``[(u_s, w_s), ...]`` and the eigenvalues. Scoring ``u_s``
        against ``K_x + diag(w_s)`` and adding ``-(N/2) log λ_s`` is
        ``log N(r̃_s; 0, λ_s K_x + diag(sigma^2))`` — see the class docstring on why
        the rescaling happens here rather than inside the kernel.
        """
        eigenvalues = np.asarray(self.eigen(values)[0], dtype=DTYPE).ravel()
        if eigenvalues.size != self.channels:
            raise LikelihoodError(
                f"the coupling returned {eigenvalues.size} eigenvalue(s) for "
                f"{self.channels} channels."
            )
        if not np.all(np.isfinite(eigenvalues)) or np.any(eigenvalues <= 0.0):
            raise LikelihoodError(
                f"the coupling B has eigenvalues {eigenvalues.tolist()}, which are not all "
                f"finite and positive — so B is not a covariance at this parameter vector and "
                f"B ⊗ K_x is not a covariance either."
            )
        rotated = self.rotate(residuals, values)
        diagonal = _as_float64(np.asarray(variance).ravel(), "shared variance")
        return (
            [
                (rotated[:, index] / np.sqrt(eigenvalues[index]), diagonal / eigenvalues[index])
                for index in range(self.channels)
            ],
            eigenvalues,
        )

    def log_prob(
        self,
        residuals: Any,
        variance: Any,
        coordinates: np.ndarray,
        values: Mapping[str, Any],
        *,
        kernel: Kernel | None = None,
    ) -> float:
        """``log N(vec(R); 0, B ⊗ K_x + I_T ⊗ diag(variance))``.

        Parameters
        ----------
        residuals
            ``observed - predicted`` per channel, in this model's own channel
            order: a sequence of ``T`` vectors or one ``(n, T)`` block.
        variance
            The shared per-sample variance, ``(n,)``.
        coordinates
            The shared grid, ``(n, d)``.
        values
            The resolved parameter mapping — the kernel's hyperparameters and
            the coupling's together.
        kernel
            The axis-bound kernel, when the caller has already built it.
        """
        bound = self._kernel if kernel is None else kernel
        points = _as_points(coordinates, "data coordinates")
        problems, eigenvalues = self._rotated_problems(residuals, variance, values)
        resolved = bound.resolve(values)
        total = 0.0
        for index, (rotated, diagonal) in enumerate(problems):
            total += self._solver.log_marginal_likelihood(
                bound, points, rotated, diagonal, resolved
            )
            total -= 0.5 * rotated.size * float(np.log(eigenvalues[index]))
        return float(total)

    def pointwise_log_prob(
        self,
        residuals: Any,
        variance: Any,
        coordinates: np.ndarray,
        values: Mapping[str, Any],
        *,
        kernel: Kernel | None = None,
    ) -> np.ndarray:
        """Per-sample leave-one-out conditional terms, **per rotated output**.

        An ``(n, T)`` block: column ``s`` is
        :meth:`GPSolver.conditional_loo` for the ``s``-th rotated output.
        ``results.md`` §6's pointwise group is what consumes it, and the
        rotated outputs rather than the channels are what it is indexed by
        because the rotated outputs are the things that are independent: a
        leave-one-out conditional of channel ``ra`` alone would condition on
        ``dec``'s value at the same epoch without saying so.
        """
        bound = self._kernel if kernel is None else kernel
        points = _as_points(coordinates, "data coordinates")
        problems, eigenvalues = self._rotated_problems(residuals, variance, values)
        resolved = bound.resolve(values)
        columns = [
            self._solver.conditional_loo(bound, points, rotated, diagonal, resolved)
            - 0.5 * float(np.log(eigenvalues[index]))
            for index, (rotated, diagonal) in enumerate(problems)
        ]
        return np.ascontiguousarray(np.column_stack(columns))

    def sample(
        self,
        predicted: Sequence[Any],
        variance: Any,
        coordinates: np.ndarray,
        values: Mapping[str, Any],
        rng: np.random.Generator,
        *,
        kernel: Kernel | None = None,
    ) -> np.ndarray:
        """One **correlated** draw of all ``T`` channels: an ``(n, T)`` block.

        Drawn in the rotated basis and rotated back, which is the generative
        statement of the same factorisation :meth:`log_prob` scores with: each
        rotated output is an independent draw from
        ``N(0, λ_s K_x + diag(variance))`` — the GP part through
        :meth:`GPSolver.latent_transform`, the same whitening the latent
        declaration uses, and the diagonal part white — and ``R = R̃ Qᵀ``
        puts the correlation back. The solver's own jitter is folded in
        exactly as :meth:`GaussianFamily.sample` folds it in, and for the same
        reason: it is part of the covariance the density scores.
        """
        bound = self._kernel if kernel is None else kernel
        points = _as_points(coordinates, "data coordinates")
        eigenvalues, rotation = self.eigen(values)
        eigenvalues = np.asarray(eigenvalues, dtype=DTYPE).ravel()
        resolved = bound.resolve(values)
        diagonal = _as_float64(np.asarray(variance).ravel(), "shared variance")
        stabiliser = float(getattr(self._solver, "jitter", 0.0) or 0.0)
        size = int(points.shape[0])
        columns = []
        for index in range(self.channels):
            whitened = rng.standard_normal(self._solver.latent_size(bound, size))
            draw = np.sqrt(eigenvalues[index]) * np.asarray(
                self._solver.latent_transform(bound, points, whitened, resolved), dtype=DTYPE
            )
            white = diagonal + eigenvalues[index] * stabiliser**2
            columns.append(draw + np.sqrt(white) * rng.standard_normal(size))
        rotated = np.column_stack(columns)
        block = rotated @ np.asarray(rotation, dtype=DTYPE).T
        means = _stack_channels(predicted, "predicted values")
        return np.ascontiguousarray(means + block)

    # -- provenance ----------------------------------------------------------

    def to_spec(self) -> dict[str, Any]:
        """The declaration, as JSON-normalisable data — ``results.md`` §14's shape."""
        return {
            "noise": type(self).__name__,
            "joint": True,
            "datasets": list(self._datasets),
            "coupling": {
                "parameterisation": self._coupling.NAME,
                "channels": self._coupling.channels,
                "parameters": list(self._coupling.parameters.names),
            },
            "kernel": self._kernel.spec().to_dict(),
            "solver": self._solver.NAME,
            "parameters": list(self.parameters.names),
        }

    def __repr__(self) -> str:
        return (
            f"JointGaussianProcessNoise({type(self._kernel).__name__}, "
            f"{self._solver.NAME}, datasets={list(self._datasets)}, "
            f"coupling={type(self._coupling).__name__})"
        )


def _masks_agree(left: FunctionSamples, right: FunctionSamples) -> bool:
    """Whether two containers call the same samples valid."""
    return bool(
        np.array_equal(
            np.asarray(left.valid).ravel(),
            np.asarray(right.valid).ravel(),
        )
    )


# ---------------------------------------------------------------------------
# The latent-GP declaration
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class LatentDeclaration:
    """What a latent-GP formulation adds to a fitting problem's parameters.

    The latent values are declared **whitened**: ``z ~ Normal(0, 1)``
    independently, with the correlation supplied deterministically by
    :meth:`GPSolver.latent_transform` as ``f = L(θ) z``. See
    :func:`latent_parameter` for why this, and not a
    :class:`~ampere.core.parameter.HierarchicalPrior`, is the declaration.

    **Two sizes since W5.4**, because the two stopped being the same number
    when an approximate solver landed. :attr:`size` is the **retained sample
    count** this declaration was built for — what
    ``Dataset._declare_latent`` passes and what
    ``FittingProblem`` validation checks the effective mask against, so that a
    transformation which masks more than the data do is caught before a run
    starts. :attr:`whitened_size` is how many **whitened variables** the
    sampler actually carries, which is the solver's to say
    (:meth:`GPSolver.latent_size`): ``size`` for both exact solvers, and the
    basis size for :class:`HilbertSpaceGP`, whose whitening is ``(N, m)``
    rather than ``(N, N)``. The parameter itself always has
    :attr:`whitened_size` elements, so the sampler's dimension, ArviZ's
    coordinate and ``latent_transform``'s argument agree by construction.
    """

    parameter: Parameter
    size: int
    solver: GPSolver

    @property
    def whitened_size(self) -> int:
        """How many whitened variables the sampler carries: ``m``, not always ``N``."""
        shape = self.parameter.shape
        return 1 if not shape else int(shape[0])

    def as_parameter_set(self) -> ParameterSet:
        """The declaration as a set, ready for ``ParameterSet.merge``."""
        return ParameterSet([self.parameter])


def latent_parameter(name: str, size: int) -> Parameter:
    """Declare ``size`` whitened latent GP values as one array-valued parameter.

    ``parameters.md`` §13 asked this contract to answer whether an
    array-valued ``HierarchicalPrior`` scales as the latent-GP declaration for
    ``N ~ 10⁵``. The answer is that a ``HierarchicalPrior`` is not the right
    declaration at *any* N, for a reason that has nothing to do with size:
    ``parameters.md`` §12.1 states that an array-valued parameter's prior is
    i.i.d. across its elements, and a GP prior is precisely *not* i.i.d. — its
    whole content is the correlation between elements. Declaring ``f`` with a
    hierarchical Normal referencing the GP hyperparameters would silently
    describe white noise with a fitted variance.

    The declaration is therefore the standard non-centred (whitened)
    parameterisation, which needs no extension to ``parameters.md`` at all:

    * ``z`` is one array-valued parameter of shape ``(N,)`` with an i.i.d.
      standard-normal prior and an ``Identity`` bijection — exactly what
      §12.1 already supports, and exactly what a PPL wants for HMC geometry;
    * the covariance enters through ``f = L(θ) z``, a deterministic transform
      owned by the :class:`GPSolver`, where it can be a Cholesky factor
      (:class:`DenseGP`) or a state-space recursion (:class:`QuasisepGP`)
      without the declaration changing.

    On scale, ``N ~ 10⁵`` is fine as a *declaration* — one ``Parameter``, one
    shape tuple, one ``PriorSpec`` — and the flat-vector pack/unpack is a
    single reshape. Three things downstream are not fine, and the likelihoods
    contract records them as obligations rather than leaving them to be
    discovered: ``ParameterSet.free_labels()`` materialises one string per
    element (W1.8 must give ArviZ a dimension and a coordinate, not 10⁵ scalar
    names); the sampler dimension becomes ``N + k``, which no gradient-free
    engine can address (hence :meth:`Likelihood.check_engine`); and forming
    ``L`` densely is O(N³), so the latent path at that scale requires
    :class:`QuasisepGP`, not :class:`DenseGP`.
    """
    if not isinstance(size, (int, np.integer)) or int(size) < 1:
        raise LikelihoodError(
            f"a latent GP declaration needs a positive number of latent values, got {size!r}."
        )
    return Parameter(
        name,
        st.norm(0.0, 1.0),
        shape=(int(size),),
        bijection=Identity(),
        description=(
            "whitened latent GP values; the correlated draw is f = L(theta) @ z, applied by "
            "GPSolver.latent_transform"
        ),
    )


# ---------------------------------------------------------------------------
# Likelihood families and their registry
# ---------------------------------------------------------------------------

_FAMILIES: dict[str, type[LikelihoodFamily]] = {}


def register_family(cls: type[LikelihoodFamily]) -> type[LikelihoodFamily]:
    """Register a family under its :attr:`~LikelihoodFamily.NAME`.

    Usable as a decorator. ``DEVELOPMENT_PLAN.md`` §4.4 requires that "users
    can add their own without touching ampere", so this is public and the
    registry is an ordinary mutable mapping — a third-party family is a class
    with one method and a decorator.
    """
    name = cls.NAME
    if not name:
        raise LikelihoodError(
            f"{cls.__name__} must declare a non-empty NAME before it can be registered."
        )
    existing = _FAMILIES.get(name)
    if existing is not None and existing is not cls:
        raise LikelihoodError(
            f"a likelihood family named {name!r} is already registered "
            f"({existing.__module__}.{existing.__qualname__}). Family names are the key W1.8's "
            f"provenance records and W1.9's lowering table use; pick another."
        )
    _FAMILIES[name] = cls
    return cls


def family_named(name: str) -> type[LikelihoodFamily]:
    """The registered family class of that name."""
    try:
        return _FAMILIES[name]
    except KeyError:
        raise LikelihoodError(
            f"no likelihood family named {name!r} is registered; the registered families are "
            f"{list(list_families())}. Add your own with @register_family."
        ) from None


def list_families() -> tuple[str, ...]:
    """Every registered family name, in registration order."""
    return tuple(_FAMILIES)


class LikelihoodFamily(Parameterised, abc.ABC):
    """The sampling distribution of one datum, given a prediction and its noise.

    One method — :meth:`log_prob` — is the whole evaluation interface, with the
    signature ``DEVELOPMENT_PLAN.md`` §4.4 specifies. ``prior_art.md`` lesson
    3M3 is the reason it stays that small: 3ML's ``get_log_like()`` is the
    right *shape* of contract, and its opacity (the plugin also owning the
    instrument response) is the part not to copy.

    Everything else a family declares is a class attribute, so declaring costs
    nothing at evaluation time and composition can be checked before a sampler
    starts. A family may also be :class:`Parameterised` — Student-t's degrees
    of freedom are an ordinary fitted parameter.
    """

    #: Registry key. Also what appears in provenance and error messages.
    NAME: ClassVar[str] = ""
    #: Whether a GP covariance can be folded into this family's own noise
    #: process and marginalised in closed form. True for the Gaussian family
    #: and — ruled 2026-09-03 (the circular complex GP, ``likelihoods.md``
    #: §17 Q6) — for the complex Gaussian.
    ANALYTIC_WITH_GP: ClassVar[bool] = False
    #: Whether the closed form :attr:`ANALYTIC_WITH_GP` declares is actually
    #: implemented, as opposed to staged for a later phase. ``False`` makes
    #: :class:`Likelihood` refuse the composition with a message naming the
    #: phase that lands it — the same declared-but-staged discipline
    #: :attr:`IMPLEMENTED` applies to a whole family, applied to one
    #: combination (see :class:`ComplexGaussianFamily`).
    GP_ANALYTIC_IMPLEMENTED: ClassVar[bool] = True
    #: Whether :meth:`log_prob` actually implements the latent-conditional form
    #: — i.e. whether it reads ``noise.latent`` and refuses to proceed without
    #: it. **False by default, deliberately**: a family that declares
    #: :attr:`Marginalisation.LATENT` but ignores the latent values would
    #: silently return the *uncorrelated* likelihood, so every GP
    #: hyperparameter and every latent value an engine sampled would leave the
    #: log-probability untouched — a fit that runs, converges, and is wrong.
    #: :class:`Likelihood` therefore refuses a latent combination whose family
    #: has not opted in, exactly as it refuses an unimplemented family.
    CONSUMES_LATENT_GP: ClassVar[bool] = False
    #: Whether this family can consume a :class:`Censoring` declaration.
    SUPPORTS_CENSORING: ClassVar[bool] = False
    #: Whether it needs per-sample uncertainties on the observed container.
    REQUIRES_UNCERTAINTY: ClassVar[bool] = True
    #: Whether it accepts complex-valued containers (VisibilitySet).
    ALLOWS_COMPLEX: ClassVar[bool] = False
    #: Whether an implementation exists, as opposed to a declared slot.
    IMPLEMENTED: ClassVar[bool] = True

    @abc.abstractmethod
    def log_prob(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        noise: NoiseParams,
    ) -> float:
        """Total log-probability of ``observed`` given ``predicted`` and ``noise``.

        Both arrays cover the retained samples only, in float64 (or complex128
        for a complex family): masking has already been applied by excision.
        """

    def sample(
        self,
        predicted: np.ndarray,
        noise: NoiseParams,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """One draw of the retained observed values, given *predicted* and *noise*.

        The generative counterpart of :meth:`log_prob` (ruled 2026-09-02,
        ``inference.md`` §19 R3): a draw from the same distribution
        ``log_prob`` scores, over the retained samples, using the same
        ``NoiseParams``. ``simulate(observe=True)`` (W1.7) is the consumer.

        The default **refuses, specifically**: ``log_prob`` defines how a
        datum is scored, not how one is generated, and ampere will not guess a
        sampling distribution — a wrong guess would silently train an SBI
        posterior on the wrong forward model. A family for which the
        observation process is well defined overrides this with it; a user
        family may do the same, which is the supported route to
        ``observe=True`` draws for an exotic observation process.

        **Amended W3.14.** The refusal was the default for every family but
        the Gaussian one until Peter's ruling of 2026-09-09 (its use case
        arrived in W3.6: a counting experiment could not
        ``simulate(observe=True)``, so neither SBC nor SBI on count data
        worked without the user subclassing :class:`PoissonFamily`). It was
        written for families whose observation process is *genuinely
        ambiguous*, and three of the shipped families are not:
        :class:`PoissonFamily`, :class:`StudentTFamily` and
        :class:`ComplexGaussianFamily` each have exactly one generative form,
        and their own ``log_prob`` already fixes which one. Those three — and
        only those three, the ruling being a list rather than a principle
        applied by each implementer — override ``sample``. Every other family
        keeps the refusal below **word for word**, :class:`CauchyFamily`
        included, and a user who wants one of them to draw supplies the
        observation process by subclassing, exactly as before. The rule the
        refusal states is unchanged; what changed is the list of families the
        contract can honestly say it knows the answer for.
        """
        raise LikelihoodError(
            f"the {self.NAME or type(self).__name__} family does not implement sample(): its "
            f"log_prob defines how a datum is scored, not how one is generated, and ampere will "
            f"not guess a sampling distribution. Subclass {type(self).__name__} and override "
            f"sample(predicted, noise, rng) with the observation process, or use "
            f"simulate(observe=False) and draw observations from the predicted containers "
            f"yourself."
        )

    def check_observed(self, observed: FunctionSamples) -> None:
        """Composition-time precondition on the observed data. Default: none.

        The family's half of the obligation ``NoiseModel.check_compatible``
        already has (ruled 2026-09-02, W1.11 gap I-5): some preconditions are
        properties of the *sampling distribution* rather than of the noise — a
        circular family needs angles in radians, ``PoissonFamily`` needs
        integer counts, a Rice family needs non-negative amplitudes. Called by
        :meth:`Likelihood.check_alignment` on the **observed** container only:
        the unit check has already forced predicted and observed to agree on
        everything a container carries, and value-range properties genuinely
        differ between the two (a Poisson *rate* is not an integer). Checks
        that belong here are the ones otherwise forced into the hot loop or
        nowhere at all. Only retained samples should be held to a
        precondition — a masked sample carries zero information.
        """

    def marginalisation_with(
        self,
        noise: NoiseModel,
        censoring: Censoring | None = None,
    ) -> Marginalisation:
        """Whether this family plus that noise model closes analytically.

        The default rule implements ``DEVELOPMENT_PLAN.md`` §4.4 exactly: an
        uncorrelated noise model is always analytic; a correlated one is
        analytic only for a family whose own noise process is Gaussian; and
        censoring under a correlated noise model is latent regardless, because
        a censored multivariate Gaussian likelihood is an orthant probability
        with no closed form beyond a handful of dimensions.
        """
        if not noise.CORRELATED:
            return Marginalisation.ANALYTIC
        if censoring is not None and censoring.any_censored:
            return Marginalisation.LATENT
        return Marginalisation.ANALYTIC if self.ANALYTIC_WITH_GP else Marginalisation.LATENT

    def _unimplemented(self) -> LikelihoodError:
        return LikelihoodError(
            f"the {self.NAME} family is declared but not implemented "
            f"(DEVELOPMENT_PLAN.md §4.4 puts it in the interface design and stages the "
            f"implementation). Its declaration is live — list_families() reports it and "
            f"composition checks against it — but it cannot be evaluated yet."
        )

    def _gp_analytic_unimplemented(self) -> LikelihoodError:
        """The staged-combination refusal: declared analytic, not yet computable.

        **Generalised at W4.2.** Until then this message named the one case it
        had — the circular complex GP, staged for Phase 4 — and the wording was
        specific to it. Phase 4 implemented that case, so what is left is the
        *discipline* rather than an instance of it, and the message says what it
        is for: a family may declare :attr:`ANALYTIC_WITH_GP` before its closed
        form exists, and composition then refuses rather than falling back to a
        latent formulation the declaration does not describe.
        """
        return LikelihoodError(
            f"the {self.NAME} family with a correlated noise model declares "
            f"Marginalisation.ANALYTIC — its own noise process is Gaussian, so the GP covariance "
            f"folds into it and integrates out — but "
            f"{type(self).__name__}.GP_ANALYTIC_IMPLEMENTED is False, so the closed form is "
            f"declared rather than written. A declared-but-staged combination is refused here "
            f"instead of being quietly computed some other way, because the other ways are a "
            f"different model: use IndependentNoise with this family, or implement log_prob's "
            f"correlated branch and set GP_ANALYTIC_IMPLEMENTED = True."
        )

    def __repr__(self) -> str:
        declared = ", ".join(repr(self.parameters[name]) for name in self.parameters.names)
        return f"{type(self).__name__}({declared})"


def _wrap_to_pi(angles: np.ndarray) -> np.ndarray:
    """Wrap *angles* (radians) into ``(-pi, pi]``.

    The one line that separates a circular family from a catastrophic one: a
    phase residual is a point on a circle, and subtracting two angles that
    straddle the branch cut gives an error near ``2 pi`` where the truth is
    near zero. Written through ``angle(exp(i x))`` rather than through a
    modulo, because that is the expression that gets the boundary and the sign
    of ``-pi`` right without a special case.
    """
    return np.asarray(np.angle(np.exp(1j * np.asarray(angles, dtype=DTYPE))), dtype=DTYPE)


def _independent_sigma(noise: NoiseParams, family: str) -> np.ndarray:
    if noise.sigma is None:
        raise LikelihoodError(
            f"the {family} family needs per-sample uncertainties, but the noise model supplied "
            f"none — the observed container has no 'uncertainty'. Attach them at construction "
            f"(uncertainty=...), or give the noise model a 'jitter' parameter to stand in for "
            f"them. Likelihood.check_alignment() catches this at composition time."
        )
    return noise.sigma


def _location_scale_log_prob(
    distribution: Any,
    standardised: np.ndarray,
    sigma: np.ndarray,
    limits: np.ndarray | None,
    *shape_args: float,
) -> float:
    """Sum a location-scale family's log-density, honouring a censoring declaration.

    This is the whole of the censored-likelihood mathematics, written once for
    every family whose standardised residual has a scipy distribution: a
    detection contributes ``log f(z) - log sigma``; an upper limit contributes
    ``log F(z)``, the probability that the truth lies below the recorded value;
    a lower limit contributes ``log (1 - F(z))``. The Tobit construction, and
    the reason issue #11's "upper limits are fairly straightforward in other
    settings" is true for uncorrelated noise.
    """
    if limits is None:
        return float(np.sum(distribution.logpdf(standardised, *shape_args) - np.log(sigma)))
    total = 0.0
    detection = limits == int(LimitKind.DETECTION)
    if np.any(detection):
        kept = standardised[detection]
        total += float(np.sum(distribution.logpdf(kept, *shape_args) - np.log(sigma[detection])))
    upper = limits == int(LimitKind.UPPER_LIMIT)
    if np.any(upper):
        total += float(np.sum(distribution.logcdf(standardised[upper], *shape_args)))
    lower = limits == int(LimitKind.LOWER_LIMIT)
    if np.any(lower):
        total += float(np.sum(distribution.logsf(standardised[lower], *shape_args)))
    return total


@register_family
class GaussianFamily(LikelihoodFamily):
    """The Gaussian family: the one whose GP marginalisation closes.

    With :class:`IndependentNoise` this is the ordinary chi-square-plus-
    normalisation log-likelihood. With :class:`GaussianProcessNoise` it is the
    GP marginal likelihood ``log N(residual; 0, K(θ) + diag(sigma²))`` — the
    flexible likelihood, delegated to the :class:`GPSolver` strategy.

    It is the only family here with ``ANALYTIC_WITH_GP = True``, and that is
    not an accident of implementation effort: the "add the GP covariance to the
    noise and integrate" trick is a property of Gaussian conjugacy, and every
    other family needs the latent path.

    Examples
    --------
    >>> import numpy as np
    >>> import astropy.units as u
    >>> import scipy.stats as st
    >>> from ampere.core import Spectrum
    >>> data = Spectrum([1.0, 2.0, 3.0] * u.um, [1.0, 2.0, 3.0] * u.Jy,
    ...                 uncertainty=[0.1, 0.1, 0.1] * u.Jy)
    >>> model = data.with_values([1.05, 2.0, 2.9])
    >>> like = Likelihood(GaussianFamily(), IndependentNoise())
    >>> round(like.log_prob(model, data), 6)
    3.52594
    >>> like.marginalisation
    <Marginalisation.ANALYTIC: 'analytic'>
    """

    NAME: ClassVar[str] = "gaussian"
    ANALYTIC_WITH_GP: ClassVar[bool] = True
    SUPPORTS_CENSORING: ClassVar[bool] = True

    def log_prob(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        noise: NoiseParams,
    ) -> float:
        residual = observed - predicted
        if noise.correlated:
            if noise.limits is not None and np.any(noise.limits != int(LimitKind.DETECTION)):
                raise LikelihoodError(
                    "censored data under a correlated (GP) noise model do not marginalise "
                    "analytically: the likelihood is a multivariate-normal orthant probability, "
                    "which has no closed form beyond a few dimensions. This combination is "
                    "declared LATENT (see Likelihood.marginalisation) and needs a modern-backend "
                    "engine that can sample the truncated latent values."
                )
            assert noise.solver is not None and noise.kernel is not None  # narrowed by .correlated
            assert noise.coordinates is not None
            return noise.solver.log_marginal_likelihood(
                noise.kernel, noise.coordinates, residual, noise.variance, noise.values
            )
        sigma = _independent_sigma(noise, self.NAME)
        if noise.limits is not None:
            return _location_scale_log_prob(st.norm, residual / sigma, sigma, noise.limits)
        # The uncensored case is written out rather than delegated: it is the
        # hot loop of every ordinary fit ampere runs, and scipy's generic
        # machinery costs several times what the closed form does.
        return float(np.sum(-0.5 * ((residual / sigma) ** 2 + _LOG_2PI) - np.log(sigma)))

    def sample(
        self,
        predicted: np.ndarray,
        noise: NoiseParams,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """A draw from the same distribution :meth:`log_prob` scores.

        * uncorrelated noise — ``x = mu + sigma z``, with the noise model's own
          ``sigma``, so a fitted ``scale`` or ``jitter`` is already in it;
        * correlated (GP) noise — ``x = mu + L z1 + sigma z2``, where ``L``
          comes from :meth:`GPSolver.latent_transform`, the same whitening the
          latent declaration uses. That is a draw from
          ``N(mu, K + diag(sigma^2))`` **up to the numerical stabilisers**:
          the solver's own jitter is part of the covariance it scores, so it
          is folded into the draw here — omitting it would draw from a
          narrower distribution than the likelihood evaluates, and the error
          is not small at the jitter values the library's own error message
          tells a user to raise. ``latent_transform``'s own factorisation
          epsilon (a relative ``1e-10`` on the diagonal of ``K``, needed so a
          smooth kernel's near-singular matrix factorises at all) also
          inflates the drawn covariance, by an amount ten orders below the
          marginal variance; the scoring path does not carry it, and aligning
          the two exactly — drawing through the same factorisation the
          marginal likelihood forms — is recorded as a Phase 2 refinement.
        """
        realisation = np.asarray(predicted, dtype=DTYPE).copy()
        sigma = None if noise.sigma is None else np.asarray(noise.sigma, dtype=DTYPE)
        if noise.correlated:
            assert noise.solver is not None and noise.kernel is not None  # narrowed by .correlated
            assert noise.coordinates is not None
            # W5.4: how many whitened variables the draw needs is the
            # solver's to say -- ``N`` for both exact solvers, ``m`` for a
            # reduced-rank one -- so that this draw and the latent-GP
            # likelihood path go through one and the same whitening.
            whitened = rng.standard_normal(
                noise.solver.latent_size(noise.kernel, realisation.shape[0])
            )
            realisation = realisation + noise.solver.latent_transform(
                noise.kernel, noise.coordinates, whitened, noise.values
            )
            stabiliser = float(getattr(noise.solver, "jitter", 0.0) or 0.0)
            if stabiliser:
                floor = np.full(realisation.shape, stabiliser)
                sigma = floor if sigma is None else np.sqrt(sigma**2 + floor**2)
        if sigma is not None:
            realisation = realisation + sigma * rng.standard_normal(realisation.shape)
        return realisation


@register_family
class StudentTFamily(LikelihoodFamily):
    """Heavy-tailed robustness to outlying samples.

    ``DEVELOPMENT_PLAN.md`` §4.4 lists Student-t for outlier robustness. Note
    that this is robustness to *individual bad samples*, which is a different
    problem from the correlated-residual misspecification
    :class:`GaussianProcessNoise` addresses — the two compose, but only through
    the latent path, because a Student-t is a scale mixture of Gaussians and
    the mixture does not commute with a GP covariance in closed form.

    ``nu`` (degrees of freedom) is an ordinary parameter: fit it, fix it, or
    tie it across datasets.
    """

    NAME: ClassVar[str] = "student_t"
    ANALYTIC_WITH_GP: ClassVar[bool] = False
    SUPPORTS_CENSORING: ClassVar[bool] = True

    def __init__(self, nu: Any = 4.0) -> None:
        self.register_parameter(_as_hyperparameter("nu", nu, None))

    def log_prob(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        noise: NoiseParams,
    ) -> float:
        sigma = _independent_sigma(noise, self.NAME)
        resolved = self.context({k: v for k, v in noise.values.items() if k in self.parameters})
        nu = _positive(resolved["nu"], "nu", self.NAME)
        standardised = (observed - predicted) / sigma
        return _location_scale_log_prob(st.t, standardised, sigma, noise.limits, nu)

    def sample(
        self,
        predicted: np.ndarray,
        noise: NoiseParams,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """A draw from the same distribution :meth:`log_prob` scores (*W3.14*).

        Location-scale, with the family's own degrees of freedom and the noise
        model's ``sigma`` as the **scale** — the same two quantities
        :meth:`log_prob` standardises by, read the same way, so a fitted
        ``scale`` or ``jitter`` is in the draw exactly as it is in the
        density::

            x = mu + sigma * t_nu

        ``sigma`` is the scale, not the standard deviation: a Student-t's
        variance is ``sigma**2 * nu / (nu - 2)`` for ``nu > 2`` and is
        undefined below that. That is not an approximation to correct for
        here — it is what ``log_prob``'s ``logpdf(z, nu) - log(sigma)`` means,
        and a draw that matched the *variance* instead would be a draw from a
        different density from the one that will score it.

        A correlated (GP) noise model is refused **by name**: a Student-t is a
        scale mixture of Gaussians, that mixture does not commute with a GP
        covariance, and the marginal of the combination is not a Student-t at
        all. ``ampere.core.Likelihood`` already refuses the composition at
        construction (the family declares no ``CONSUMES_LATENT_GP``), so this
        guard is for a caller assembling :class:`NoiseParams` directly.
        """
        if noise.correlated:
            raise LikelihoodError(
                "the student_t family cannot draw an observation under a correlated (GP) noise "
                "model: a Student-t is a scale mixture of Gaussians, the mixture does not "
                "commute with a GP covariance, and the marginal of the two together is not a "
                "Student-t — so there is no location-scale draw to make. Compose this family "
                "with IndependentNoise, or subclass StudentTFamily and override "
                "sample(predicted, noise, rng) with the observation process you mean."
            )
        sigma = _independent_sigma(noise, self.NAME)
        resolved = self.context({k: v for k, v in noise.values.items() if k in self.parameters})
        nu = _positive(resolved["nu"], "nu", self.NAME)
        location = np.asarray(predicted, dtype=DTYPE)
        return location + sigma * rng.standard_t(nu, size=location.shape)


@register_family
class CauchyFamily(LikelihoodFamily):
    """Student-t with ``nu = 1``: maximal outlier tolerance, no finite variance.

    Registered separately because it takes no parameters and because a user
    reaching for "the Cauchy likelihood" should find it by that name.
    """

    NAME: ClassVar[str] = "cauchy"
    ANALYTIC_WITH_GP: ClassVar[bool] = False
    SUPPORTS_CENSORING: ClassVar[bool] = True

    def log_prob(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        noise: NoiseParams,
    ) -> float:
        sigma = _independent_sigma(noise, self.NAME)
        standardised = (observed - predicted) / sigma
        return _location_scale_log_prob(st.cauchy, standardised, sigma, noise.limits)


@register_family
class ComplexGaussianFamily(LikelihoodFamily):
    """Circular complex Gaussian: the standard interferometric noise model.

    ``results_schema.md`` §16 hands this contract the statement that a
    ``VisibilitySet``'s real-valued uncertainty "encodes the circular complex
    Gaussian only; a non-circular or amplitude/phase noise model must supply
    its own structure". This family is that model, and nothing more: the real
    and imaginary parts are independent ``Normal(0, sigma²)``, so

    ``log p = Σ [ -log(2π sigma²) - |y - μ|² / (2sigma²) ]``.

    A non-circular model (unequal real/imaginary variances, or a correlation
    between them) is a *different noise model* supplying a 2x2 structure per
    sample, and belongs to whoever needs it; the amplitude/phase formulations
    are :class:`RiceFamily` and :class:`VonMisesFamily`.

    :attr:`ANALYTIC_WITH_GP` is ``True`` (ruled 2026-09-03, ``likelihoods.md``
    §17 Q6), with the **circular complex GP** as the fixed meaning: one real
    kernel applied independently to the real and imaginary parts — equal
    component covariances, zero pseudo-covariance — which marginalises in
    closed form exactly as the real Gaussian does.

    **W4.2 implements it**, so :attr:`GP_ANALYTIC_IMPLEMENTED` is now ``True``
    and the composition that was refused with Phase 4 named is the flexible
    likelihood on visibilities. The closed form is the real 2N-dimensional
    Gaussian's, with the block structure used rather than materialised. Write
    ``S = K(θ) + diag(σ²)`` for the **per-component** covariance; the real
    ``2N`` covariance of ``(Re r, Im r)`` is ``diag(S, S)``, block diagonal
    with a zero off-diagonal block — that is exactly what circularity says —
    so

    ``log p = -½ [ rᵉ ᵀ S⁻¹ rᵉ + rⁱ ᵀ S⁻¹ rⁱ + 2 log|S| + 2N log 2π ]``

    with ``rᵉ`` and ``rⁱ`` the real and imaginary parts of ``observed -
    predicted``. One factorisation of the ``N by N`` matrix ``S``, two
    triangular solves against it, ``log|S|`` computed once and counted twice.
    Forming the ``2N by 2N`` matrix instead would cost eight times as much
    arithmetic to carry a zero block the model has already declared. The
    mechanism is :class:`GPSolver`'s ``(n, k)`` right-hand side, with ``k = 2``
    columns; :attr:`GPSolver.STACKED_RESIDUALS` is the declaration a solver
    makes that it can take them, and :meth:`GaussianProcessNoise.check_compatible`
    refuses one that cannot — the O(N) ``QuasisepGP`` — by name.

    **``sigma`` is the per-component standard deviation throughout**, as
    ``results_schema.md`` §16 says a :class:`~ampere.core.VisibilitySet`'s
    real-valued uncertainty is, and the kernel's ``amplitude`` is a
    per-component marginal standard deviation for the same reason: ``K``
    appears once per component in the expression above, so a correlated
    calibration error of RMS ``a`` in each of the real and the imaginary parts
    is the kernel with ``amplitude = a``. The total modulus variance
    ``E|r|²`` is ``2 (K_ii + σ²)``.
    """

    NAME: ClassVar[str] = "complex_gaussian"
    ALLOWS_COMPLEX: ClassVar[bool] = True
    ANALYTIC_WITH_GP: ClassVar[bool] = True

    def log_prob(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        noise: NoiseParams,
    ) -> float:
        if noise.correlated:
            return self._gp_log_prob(predicted, observed, noise)
        sigma = _independent_sigma(noise, self.NAME)
        residual = np.abs(observed - predicted)
        variance = sigma**2
        return float(np.sum(-(residual**2) / (2.0 * variance) - _LOG_2PI - np.log(variance)))

    def _gp_log_prob(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        noise: NoiseParams,
    ) -> float:
        """The circular complex GP marginal (**W4.2**). See the class docstring.

        Three lines of arithmetic, and the whole of the circular model is in the
        middle one: the complex residual becomes two real columns, and the
        solver is handed both against **one** covariance. The factorisation is
        shared because the solver sees one matrix and one two-column right-hand
        side, not because anything here caches a factor.
        """
        assert noise.solver is not None and noise.kernel is not None  # narrowed by .correlated
        assert noise.coordinates is not None
        residual = np.asarray(observed, dtype=np.complex128) - np.asarray(
            predicted, dtype=np.complex128
        )
        return noise.solver.log_marginal_likelihood(
            noise.kernel,
            noise.coordinates,
            _stacked_components(residual),
            noise.variance,
            noise.values,
        )

    def sample(
        self,
        predicted: np.ndarray,
        noise: NoiseParams,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """A draw from the same distribution :meth:`log_prob` scores (*W3.14*).

        Circular symmetry is the whole content of the family, so it is the
        whole content of the draw: the real and the imaginary parts each get
        an independent ``Normal(0, sigma**2)``, with no correlation between
        them and no difference in their variances::

            x = mu + sigma * (z_re + i z_im)

        **``sigma`` is the per-component standard deviation, not the total**,
        and that convention is read off :meth:`log_prob` rather than assumed.
        Its density is ``-|y - mu|**2 / (2 sigma**2) - log(2 pi) - log(sigma**2)``,
        which is exactly the joint density of two independent
        ``Normal(0, sigma**2)`` components — ``(2 pi sigma**2)**-1
        exp(-(r_re**2 + r_im**2) / (2 sigma**2))``. The total variance
        ``E|x - mu|**2`` is therefore ``2 sigma**2``. Halving the scale here
        to make ``sigma`` the total would draw from a distribution the
        library's own density does not score, and the error would look like a
        factor nobody could see in an amplitude plot.

        **Under a GP (W4.2)** the same circular model governs the draw, which
        means it is *not* one complex GP realisation but two real ones::

            x = mu + (L z1 + sigma w1) + i (L z1' + sigma w1')

        with ``L`` from :meth:`GPSolver.latent_transform` and the two
        whitened blocks drawn independently. The two components share the
        matrix ``L`` — that is the "equal component covariances" half of
        circularity — and share nothing else — that is the "zero
        pseudo-covariance" half. Drawing one realisation and using it for both
        components would give a draw perfectly correlated between real and
        imaginary parts, whose modulus statistics the density above does not
        score; the two-column whitened block is what keeps the draw and the
        density the same distribution. The solver's own ``jitter`` is folded
        into the diagonal exactly as :meth:`GaussianFamily.sample` folds it,
        and for the same reason.
        """
        sigma = _independent_sigma(noise, self.NAME)
        mean = np.asarray(predicted, dtype=np.complex128)
        if not noise.correlated:
            real = rng.standard_normal(mean.shape)
            imaginary = rng.standard_normal(mean.shape)
            return mean + sigma * (real + 1j * imaginary)
        assert noise.solver is not None and noise.kernel is not None  # narrowed by .correlated
        assert noise.coordinates is not None
        # W5.4: ``latent_size`` rows, two columns -- the circular pair shares
        # one covariance and so one whitening, whatever its rank.
        whitened = rng.standard_normal((noise.solver.latent_size(noise.kernel, mean.shape[0]), 2))
        correlated = noise.solver.latent_transform(
            noise.kernel, noise.coordinates, whitened, noise.values
        )
        stabiliser = float(getattr(noise.solver, "jitter", 0.0) or 0.0)
        scale = sigma if not stabiliser else np.sqrt(sigma**2 + stabiliser**2)
        independent = rng.standard_normal((mean.shape[0], 2))
        components = np.asarray(correlated, dtype=DTYPE) + scale[:, None] * independent
        return mean + components[:, 0] + 1j * components[:, 1]


@register_family
class PoissonFamily(LikelihoodFamily):
    """Counting statistics — and the family that exercises the latent path.

    With :class:`IndependentNoise` this is the plain Poisson log-pmf: the
    prediction is the expected count, the observation is an integer, and there
    is no sigma anywhere (:attr:`REQUIRES_UNCERTAINTY` is ``False``).

    With :class:`GaussianProcessNoise` it is the case
    ``DEVELOPMENT_PLAN.md`` §4.4 singles out. There is no covariance matrix to
    add: robustness needs a latent-GP formulation, ``counts ~ Poisson(rate ·
    exp(f))`` with ``f ~ GP``, marginalised numerically. So the combination
    declares :attr:`Marginalisation.LATENT`, :meth:`Likelihood.check_engine`
    refuses a gradient-free engine, and ``log_prob`` requires the latent values
    to be supplied — it is the *conditional* likelihood given ``f``, which is
    what an HMC/VI backend evaluates per gradient step.

    Examples
    --------
    >>> import numpy as np
    >>> import astropy.units as u
    >>> import scipy.stats as st
    >>> from ampere.core import Spectrum
    >>> counts = Spectrum([1.0, 2.0, 3.0] * u.um, [4.0, 7.0, 2.0])
    >>> rate = counts.with_values([3.5, 6.0, 2.5])
    >>> analytic = Likelihood(PoissonFamily(), IndependentNoise())
    >>> analytic.marginalisation
    <Marginalisation.ANALYTIC: 'analytic'>
    >>> round(analytic.log_prob(rate, counts), 6)
    -5.010413
    >>> flexible = Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0)))
    >>> flexible.marginalisation
    <Marginalisation.LATENT: 'latent'>
    >>> flexible.check_engine(differentiable=False, engine="emcee")
    Traceback (most recent call last):
        ...
    ampere.core.exceptions.LikelihoodError: emcee cannot run this likelihood: the poisson family
    with a GaussianProcessNoise noise model introduces N latent values ...
    """

    NAME: ClassVar[str] = "poisson"
    ANALYTIC_WITH_GP: ClassVar[bool] = False
    CONSUMES_LATENT_GP: ClassVar[bool] = True
    REQUIRES_UNCERTAINTY: ClassVar[bool] = False

    def check_observed(self, observed: FunctionSamples) -> None:
        """Integer counts, checked once at composition (ruled 2026-09-02, X-3).

        This ran in ``log_prob`` before the ``check_observed`` hook existed —
        an O(N) scan of the *data* on every evaluation, which §13's
        compile-once/evaluate-many split says is the wrong place. The
        ``rate > 0`` guard stays in the hot loop, because it is a property of
        the prediction and genuinely varies per draw.
        """
        counts = np.asarray(observed.values).ravel()
        kept = counts[np.asarray(observed.valid).ravel()]
        if np.any(kept < 0.0) or not np.all(kept == np.round(kept)):
            raise LikelihoodError(
                "the poisson family needs non-negative integer counts, but the observed values "
                "are not integral. Counts are counts; if the data are rates, multiply by the "
                "exposure in the instrument chain (W1.5) rather than here."
            )

    def log_prob(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        noise: NoiseParams,
    ) -> float:
        counts = np.asarray(observed)
        rate = np.asarray(predicted, dtype=DTYPE)
        if noise.correlated:
            if noise.latent is None:
                raise LikelihoodError(
                    "a Poisson likelihood with a GaussianProcessNoise model is a latent-variable "
                    "model: it is only defined given the latent function f, as "
                    "counts ~ Poisson(rate * exp(f)). Supply latent=... (an inference engine "
                    "that samples f does), or use IndependentNoise for the analytic case."
                )

            latent = _as_float64(noise.latent, "latent GP values")
            if latent.shape != rate.shape:
                raise LikelihoodError(
                    f"the latent GP values have shape {latent.shape} but there are {rate.shape} "
                    f"retained samples. One latent value per retained sample."
                )
            rate = rate * np.exp(latent)
        if np.any(rate <= 0.0):
            raise LikelihoodError(
                "the poisson family needs a strictly positive expected count; the model "
                "predicted a value <= 0."
            )
        return float(np.sum(st.poisson.logpmf(counts, rate)))

    def sample(
        self,
        predicted: np.ndarray,
        noise: NoiseParams,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """A draw from the same distribution :meth:`log_prob` scores (*W3.14*).

        ``rng.poisson(rate)``, returned as float because that is the dtype the
        containers hold — the values are integral, and
        :meth:`check_observed` will hold them to that when they are scored.

        **What ``predicted`` means here is what it means in**
        :meth:`log_prob`: the expected count *before* the latent is applied.

        * :class:`IndependentNoise` — the rate is ``predicted`` itself, and
          the draw is ``counts ~ Poisson(predicted)``. Mean and variance are
          both ``predicted``, which is what the conformance rows assert.
        * :class:`GaussianProcessNoise` — the family's model is
          ``counts ~ Poisson(predicted * exp(f))``, so the rate is
          ``predicted * exp(f)``. ``f`` is the latent the caller supplied, if
          one was: that is the **conditional** draw, at the same ``f``
          ``log_prob`` scores at, which is what ``simulate(observe=True)``
          needs for the whitened ``z`` in θ and the observation to describe
          one model rather than two. With no latent supplied, a fresh GP
          realisation is drawn through :meth:`GPSolver.latent_transform` —
          the same whitening the latent declaration uses, exactly as
          :meth:`GaussianFamily.sample`'s GP branch draws its own. The
          marginal of *that* draw is over-dispersed relative to a Poisson
          (mean ``predicted * exp(K_ii / 2)``), and correctly so: it is a
          log-normal mixture of Poissons, which is what the model says.

        A non-positive rate is refused **by name** rather than left to numpy,
        whose message for it names neither the family nor the prediction. The
        guard is ``rate <= 0``, the same as ``log_prob``'s (ruled 2026-09-10 on
        W3.14's finding): a rate of exactly zero is a well-defined draw — a
        point mass at zero counts — but the density declines to score it, and
        a draw the fitting likelihood cannot score is exactly what
        ``simulate(observe=True)`` must never hand an SBC or SBI consumer.
        """
        rate = np.asarray(predicted, dtype=DTYPE)
        if noise.correlated:
            latent = noise.latent
            if latent is None:
                assert noise.solver is not None and noise.kernel is not None  # .correlated
                assert noise.coordinates is not None
                whitened = rng.standard_normal(
                    noise.solver.latent_size(noise.kernel, rate.shape[0])
                )
                latent = noise.solver.latent_transform(
                    noise.kernel, noise.coordinates, whitened, noise.values
                )
            values = _as_float64(latent, "latent GP values")
            if values.shape != rate.shape:
                raise LikelihoodError(
                    f"the latent GP values have shape {values.shape} but there are {rate.shape} "
                    f"retained samples. One latent value per retained sample."
                )
            rate = rate * np.exp(values)
        if not np.all(np.isfinite(rate)) or np.any(rate <= 0.0):
            raise LikelihoodError(
                "the poisson family needs a finite, positive expected count to draw from; "
                "the model predicted a value that is zero, negative or not finite. A count is "
                "a count: constrain the prediction to the positive half-line (a Log bijection "
                "on the norm, or a positive-support prior) rather than clipping the rate here, "
                "which would draw from a distribution log_prob does not score -- the guard is "
                "log_prob's own, so every draw this returns is one the fit can score."
            )
        return np.asarray(rng.poisson(rate), dtype=DTYPE)


@register_family
class RiceFamily(LikelihoodFamily):
    """Rician amplitude noise — polarised intensity, debiased visibility amplitudes.

    **Declared, not implemented, with no scheduled phase** (amended **W5.2**:
    the previous wording promised "the implementation is Phase 4's", and
    Phase 4 closed 2026-09-15 without implementing it). The interface is
    fixed (ruled 2026-09-03, ``likelihoods.md`` §17 Q3) because fixing it cost
    nothing and a future implementation should not have to relitigate the
    shape: the *model* predicts the underlying complex value — which is what
    an interferometric model actually produces — and an ``Amplitude`` step in
    the instrument chain takes the modulus, so this family receives real,
    non-negative amplitudes as its prediction. It consumes the *same*
    per-sample sigma as the Gaussian family (it is the amplitude of a
    circular complex Gaussian), so a ``VisibilitySet`` needs no extra
    structure to support it. Fixing the interface is not the same as
    scheduling the work; ask for it, with the use case that needs it, if you
    want it put on a phase's plan.
    """

    NAME: ClassVar[str] = "rice"
    IMPLEMENTED: ClassVar[bool] = False

    def log_prob(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        noise: NoiseParams,
    ) -> float:
        raise self._unimplemented()

    def _unimplemented(self) -> LikelihoodError:
        """The sharper refusal **W5.2** gives this one declared-but-scheduleless family.

        Overrides the base :meth:`LikelihoodFamily._unimplemented`, whose
        wording ("stages the implementation") reads as a promise of a phase —
        true of the *mechanism* a declared-but-unimplemented family uses, but
        not of this family's own situation, which has no phase attached.
        :meth:`Likelihood.__init__`'s composition-time refusal calls this same
        method, so both entry points (composing a ``Likelihood`` and calling
        ``log_prob`` on a bare instance) say the same true thing.
        """
        return LikelihoodError(
            f"the {self.NAME} family is declared but not implemented, and unlike a family "
            f"staged for a specific phase, it has no scheduled one. The interface is fixed "
            f"(likelihoods.md §17 Q3) so a future implementation has a shape to target, but "
            f"nobody has asked for it with a use case yet — that is what would put it on a "
            f"phase's plan. It is visible in list_families() so the target set is on record."
        )


@register_family
class VonMisesFamily(LikelihoodFamily):
    """Wrapped/von Mises phase noise — closure phases, position angles.

    **Implemented at W4.1**, with the interferometric modality that needed it;
    the interface is the one fixed at the freeze. The observed and predicted
    values are angles in **radians** and the residual is *wrapped*, not
    subtracted; the concentration is ``kappa = 1/sigma**2`` **per sample**
    from the container's own uncertainties (ruled 2026-09-03,
    ``likelihoods.md`` §17 Q4) — exact in the small-sigma limit, which is
    where closure-phase practice lives. A fitted global ``kappa`` that ignores
    the per-sample uncertainties is a different model, and a user family if
    anyone wants it; the case it is really meant to serve ("the pipeline
    underestimates its closure-phase errors") is already
    ``IndependentNoise(scale=...)``, which gives ``kappa = 1/(s sigma)**2``
    with ``s`` an ordinary fitted parameter.

    The density is the normalised one::

        log p = kappa cos(delta) - kappa - log(2 pi) - log(i0e(kappa))

    with ``delta`` the residual wrapped into ``(-pi, pi]``. The normalisation
    is not decoration and is not a constant that cancels: an *unnormalised*
    wrapped Gaussian has mass over the circle that depends on sigma, sigma
    varies per triangle in every real dataset, and at sigma = pi it assigns
    *less* density to a perfect match than the uniform distribution does —
    which is impossible for a density on the circle (``interferometry.md``
    §5). :func:`scipy.special.i0e` is the exponentially scaled ``I0``, which
    is what keeps the normalisation finite at the large ``kappa`` a
    well-measured closure phase produces.

    Why a wrapped family at all, rather than a Gaussian on the wrapped
    residual: on phases, an ordinary :class:`GaussianFamily` computes
    ``observed - predicted`` unwrapped, so a 2° error straddling the branch
    cut is charged as a 358° one — a 5 000-nat penalty on a triangle that fits
    perfectly, silently, and no sampler recovers from it.
    """

    NAME: ClassVar[str] = "von_mises"
    IMPLEMENTED: ClassVar[bool] = True

    def check_observed(self, observed: FunctionSamples) -> None:
        """The container must hold angles in radians.

        Gap I-5's own example (``interferometry.md`` §5): a von Mises family
        handed degrees is silently accepted by every other check — the units
        agree with each other, the shapes agree, the kinds agree — and returns
        a number that is nonsense by a factor of ``(180/pi)**2`` in the
        concentration. ``None`` is accepted as "bare radians", because the
        arithmetic here is unit-free and a container built from plain arrays
        is an ordinary thing to have; any other unit is refused by name, with
        the one-call fix.
        """
        unit = observed.unit
        if unit is None or unit == u.rad:
            return
        raise LikelihoodError(
            f"the {self.NAME} family scores angles in radians, but the observed "
            f"{type(observed).__name__} is in {unit}. Its residual is wrapped into (-pi, pi] and "
            f"its concentration is 1/sigma**2, so a container in {unit} is not merely rescaled — "
            f"it is a different distribution. Convert once at composition time with "
            f".to_unit(u.rad); units never reach the hot loop (DEVELOPMENT_PLAN.md §7)."
        )

    def log_prob(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        noise: NoiseParams,
    ) -> float:
        if noise.correlated:
            raise LikelihoodError(
                f"the {self.NAME} family with a GaussianProcessNoise model is a latent-variable "
                f"model — a GP added to a *wrapped* observable does not marginalise in closed "
                f"form — and this family does not implement the latent-conditional log_prob "
                f"(CONSUMES_LATENT_GP is False). A latent GP on closure phases is a real and "
                f"useful thing, reachable on the native backends under NUTS or VI the way "
                f"PoissonFamily's is, and it is Phase 5's (W5.1). Use IndependentNoise here."
            )
        sigma = _independent_sigma(noise, self.NAME)
        kappa = 1.0 / sigma**2
        delta = _wrap_to_pi(np.asarray(observed, dtype=DTYPE) - np.asarray(predicted, dtype=DTYPE))
        return float(
            np.sum(kappa * (np.cos(delta) - 1.0) - _LOG_2PI - np.log(scipy.special.i0e(kappa)))
        )

    def sample(
        self,
        predicted: np.ndarray,
        noise: NoiseParams,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """A draw from the same distribution :meth:`log_prob` scores (*W4.1*).

        ``rng.vonmises(mu, kappa)`` per sample, with ``mu`` the predicted angle
        and ``kappa = 1/sigma**2`` read exactly as the density reads it — a
        fitted ``scale`` or ``jitter`` is in the draw because it is in the
        ``sigma`` the density uses. The result is in ``(-pi, pi]``, numpy's own
        range for the distribution and this family's own convention for an
        angle; a predicted angle outside that range draws around its wrapped
        position, which is the same statement the wrapped residual makes.

        This is the family the sampling-form principle (``likelihoods.md`` §3,
        ruled by Peter 2026-09-10) brings in with its likelihood rather than
        after it: a data type that can be fitted but not simulated breaks SBC
        and SBI for exactly the users who brought the data type, and closure
        phases are that data type here. A correlated noise model is refused by
        the same message the density refuses it with — there is no marginal to
        draw from either.
        """
        if noise.correlated:
            raise LikelihoodError(
                f"the {self.NAME} family cannot draw under a correlated noise model: a GP added "
                f"to a wrapped observable is a latent-variable model whose latent-conditional "
                f"form this family does not implement (Phase 5, W5.1). Use IndependentNoise."
            )
        sigma = _independent_sigma(noise, self.NAME)
        mean = np.asarray(predicted, dtype=DTYPE)
        return np.asarray(rng.vonmises(mean, 1.0 / sigma**2), dtype=DTYPE)


# ---------------------------------------------------------------------------
# Likelihood: the composition
# ---------------------------------------------------------------------------


class Likelihood(Parameterised):
    """A family, a noise model and (optionally) a censoring declaration.

    This is the object a ``Dataset`` holds (W1.7) and the thing that turns
    "predicted vs observed containers" into a number. It owns three
    responsibilities the pieces below it deliberately do not:

    1. **Masking.** The effective mask is the union of the predicted and
       observed containers' masks — a sample is used only if *both* are valid —
       read through ``weights()`` and applied by **excision**. See
       :meth:`log_prob`.
    2. **Composition-time checking.** :meth:`check_alignment` is the O(N) pass
       that belongs once, not per evaluation.
    3. **The marginalisation declaration** and the engine capability check
       built on it.

    Its parameters are the family's and the noise model's, in one flat
    :class:`~ampere.core.parameter.ParameterSet` — not a nested merge, because
    ``parameters.md`` §12.4 records that ``merge`` is not associative and W1.7
    must merge every dataset in one call. A name collision between the family
    and the noise model raises here rather than being silently qualified.

    Examples
    --------
    >>> import numpy as np
    >>> import astropy.units as u
    >>> import scipy.stats as st
    >>> from ampere.core import Spectrum
    >>> flexible = Likelihood(
    ...     GaussianFamily(),
    ...     GaussianProcessNoise(Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0))),
    ... )
    >>> flexible.parameters.free_names
    ('amplitude', 'length_scale')
    >>> flexible.marginalisation
    <Marginalisation.ANALYTIC: 'analytic'>
    """

    def __init__(
        self,
        family: LikelihoodFamily,
        noise: NoiseModel | None = None,
        *,
        censoring: Censoring | None = None,
    ) -> None:
        if not isinstance(family, LikelihoodFamily):
            raise LikelihoodError(
                f"Likelihood needs a LikelihoodFamily, got {type(family).__name__}. Register a "
                f"custom family with @register_family; it needs one method, log_prob."
            )
        if not family.IMPLEMENTED:
            # W5.2: delegated to the family's own _unimplemented() rather than
            # a message hardcoded here, so a family with no scheduled phase
            # (RiceFamily) can say so plainly instead of composition
            # promising a schedule that does not exist for it.
            raise family._unimplemented()
        noise = IndependentNoise() if noise is None else noise
        if not isinstance(noise, NoiseModel):
            raise LikelihoodError(f"Likelihood needs a NoiseModel, got {type(noise).__name__}.")
        # W5.9. A joint noise model's covariance couples several datasets, and a
        # Likelihood scores one; accepting it here would silently score the
        # marginal B_tt K_x of one channel, which drops the cross-covariance
        # that is the whole model. Declared rather than isinstance-checked, so a
        # user's own joint noise model is refused here too.
        if getattr(noise, "JOINT", False):
            raise LikelihoodError(
                f"{type(noise).__name__} declares JOINT = True: it is one correlated process "
                f"over several datasets' channels, and a Likelihood scores one dataset. Declare "
                f"it on the DatasetCollection instead — DatasetCollection({{...}}, "
                f"joint={{'label': noise}}) — and leave each member dataset's own likelihood as "
                f"Likelihood(GaussianFamily()). Scoring one channel through a Likelihood would "
                f"drop the cross-covariance between the channels, which is the model."
            )
        if censoring is not None and not family.SUPPORTS_CENSORING:
            raise LikelihoodError(
                f"the {family.NAME} family does not consume a censoring declaration, but one was "
                f"given. Families that do: "
                f"{sorted(n for n, c in _FAMILIES.items() if c.SUPPORTS_CENSORING)}."
            )
        self._family = family
        self._noise = noise
        self._censoring = censoring
        # Latency that follows from the family and noise model alone is a fact
        # about the pair, so it is refused here. Latency that only censoring
        # induces depends on whether a limit survives the mask, so it is
        # refused in check_alignment, which has the data.
        if family.marginalisation_with(noise) is Marginalisation.LATENT:
            self._refuse_unconsumed_latent(family, noise, censored=False)
        # A combination whose closed-form GP marginalisation is declared but
        # staged (the circular complex GP, Phase 4) is refused the same way an
        # unimplemented family is: at composition, with the schedule named.
        if (
            noise.CORRELATED
            and family.marginalisation_with(noise) is Marginalisation.ANALYTIC
            and not family.GP_ANALYTIC_IMPLEMENTED
        ):
            raise family._gp_analytic_unimplemented()
        for parameter in noise.parameters:
            self.register_parameter(parameter)
        for parameter in family.parameters:
            if parameter.name in self.parameters:
                raise LikelihoodError(
                    f"the {family.NAME} family and the {type(noise).__name__} noise model both "
                    f"declare a parameter named {parameter.name!r}. A Likelihood holds one flat "
                    f"namespace (parameters.md §12.4: merge is not associative, so this contract "
                    f"does not nest one); rename one of them."
                )
            self.register_parameter(parameter)

    @staticmethod
    def _refuse_unconsumed_latent(
        family: LikelihoodFamily, noise: NoiseModel, *, censored: bool
    ) -> None:
        """Refuse a latent combination whose family does not implement one.

        The failure this prevents is the worst kind available to a likelihood
        contract: a family that declares :attr:`Marginalisation.LATENT` but
        whose ``log_prob`` never reads ``noise.latent`` returns the
        *uncorrelated* likelihood, so every GP hyperparameter and every latent
        value an engine samples leaves the log-probability untouched. Their
        posteriors come back as their priors and the physical parameters are
        biased exactly as they would have been under the rigid likelihood — a
        fit that runs, converges, and is wrong, which is what this whole
        contract exists to make impossible.
        """
        if family.CONSUMES_LATENT_GP:
            return
        because = (
            "with a censoring declaration on correlated noise"
            if censored
            else f"with a {type(noise).__name__} noise model"
        )
        raise LikelihoodError(
            f"the {family.NAME} family {because} is a latent-variable model (see "
            f"Marginalisation.LATENT), but {type(family).__name__} does not implement the "
            f"latent-conditional log_prob — its CONSUMES_LATENT_GP is False. Evaluating it would "
            f"return the *uncorrelated* likelihood, leaving the GP hyperparameters and every "
            f"latent value with no effect on the log-probability at all: a fit that runs, "
            f"converges, and is wrong. Use IndependentNoise with this family, or the gaussian "
            f"family, whose GP marginalises in closed form. "
            f"family.marginalisation_with(noise, censoring) reports the declaration without "
            f"composing anything."
        )

    # -- declarations --------------------------------------------------------

    @property
    def family(self) -> LikelihoodFamily:
        """The sampling distribution."""
        return self._family

    @property
    def noise(self) -> NoiseModel:
        """The noise model."""
        return self._noise

    @property
    def censoring(self) -> Censoring | None:
        """The censoring declaration, if any."""
        return self._censoring

    @property
    def capability_parts(self) -> tuple[object, ...]:
        """The objects whose capability declarations this likelihood depends on.

        **W2.13, fold-in 7** (ruled 2026-09-07; ``inference.md`` §10a). The
        noise model always, and the GP solver as well when a GP is declared —
        the two pieces a log-likelihood's arithmetic actually goes through
        beyond the family. :attr:`Dataset.capability_parts` appends these to
        the instrument steps, so
        :func:`~ampere.core.dataset.declared_capabilities` sees them and
        ``FittingProblem.differentiable``/``.backend`` become true statements
        about the whole evaluation rather than about its first half.

        The **family is deliberately absent**, and that is the ruling rather
        than an oversight: a :class:`LikelihoodFamily` is a declaration of a
        sampling distribution, evaluated by whichever path the problem is on
        (the reference path in numpy, a backend's realisation natively from
        the same closed form), so it has no backend of its own to declare. A
        flag on it would either be a fiction or force every user family to
        pick a backend before it could be composed.

        **The kernel joined at W3.8** (ruled by Peter 2026-09-08 on W2.4
        slice 3's carried finding), reversing what this docstring used to say.
        The argument for leaving it out was that the kernel is consumed *by*
        the solver, which is the piece that decides whether the covariance is
        built in numpy or natively — and that argument does not survive
        contact with the code: a solver builds the covariance by *calling*
        ``Kernel.matrix``, so a numpy kernel inside a native solver returns a
        numpy array, the solver converts it, and the conversion detaches the
        graph. The composed problem then declared ``differentiable=True``
        while its kernel amplitude and length scale got no gradient at all —
        the same silent degradation fold-in 7 removed for solvers and noise
        models, one level down. It is here now, so the flags cover the whole
        of what one evaluation passes through and a foreign kernel is a
        backend disagreement like any other.

        One consequence recorded deliberately, as W2.1 recorded the same one
        for ``describe``: this is a class attribute of a
        :class:`~ampere.core.parameter.Parameterised`, so ``capability_parts``
        joins the names a family's or noise model's parameter may not shadow.
        """
        if isinstance(self._noise, GaussianProcessNoise):
            return (self._noise, self._noise.solver, self._noise.kernel)
        return (self._noise,)

    @property
    def marginalisation(self) -> Marginalisation:
        """Whether this combination integrates its noise process in closed form.

        Answered from the *declaration* alone, before any data are seen, and so
        deliberately conservative about censoring: a limit that a container's
        mask happens to exclude still counts here, because this property does
        not know which container it will meet. Use
        :meth:`marginalisation_for` once the observed container is in hand —
        which is the answer W1.7 should act on.
        """
        return self._family.marginalisation_with(self._noise, self._censoring)

    def marginalisation_for(self, observed: FunctionSamples) -> Marginalisation:
        """The marginalisation this combination actually needs for these data.

        Identical to :attr:`marginalisation` except that censoring is read
        through the mask, so a limit on a masked sample does not force the
        latent path. §9's rule is that masking beats censoring — a masked
        sample contributes nothing whatever its limit kind — and that has to
        hold for the declaration as well as for the arithmetic, or a problem
        that is analytic in fact gets refused a gradient-free engine.
        """
        if self._censoring is None:
            return self._family.marginalisation_with(self._noise, None)
        self._censoring.check_against(observed)
        retained = Censoring(np.asarray(self._censoring.kinds)[np.asarray(observed.valid).ravel()])
        return self._family.marginalisation_with(self._noise, retained)

    def latent_declaration(self, size: int, name: str = "latent") -> LatentDeclaration:
        """The latent parameters this combination adds, for ``size`` samples.

        Raises when the combination is analytic: an analytic likelihood has no
        latent values, and quietly returning an empty declaration would let a
        caller build a sampler dimension that does nothing.
        """
        if self.marginalisation is Marginalisation.ANALYTIC:
            raise LikelihoodError(
                f"the {self._family.NAME} family with a {type(self._noise).__name__} noise model "
                f"marginalises analytically, so it declares no latent values. "
                f"Likelihood.marginalisation tells you which path you are on."
            )
        if not isinstance(self._noise, GaussianProcessNoise):
            raise LikelihoodError(
                f"a latent declaration needs a GaussianProcessNoise model to supply the "
                f"whitening transform, but this Likelihood has a {type(self._noise).__name__}."
            )
        # W5.4: *size* is the retained-sample count the caller has, and it
        # stays the declaration's ``size`` because that is what the mask
        # invariance is checked against. How many *whitened variables* it
        # implies is the solver's to say — the same number for both exact
        # solvers, and the basis size for a reduced-rank one — and that is
        # what the parameter is shaped by.
        return LatentDeclaration(
            parameter=latent_parameter(
                name, self._noise.solver.latent_size(self._noise.kernel, int(size))
            ),
            size=int(size),
            solver=self._noise.solver,
        )

    def check_engine(
        self,
        *,
        differentiable: bool,
        engine: str = "this engine",
        observed: FunctionSamples | None = None,
    ) -> None:
        """Refuse an engine that cannot deliver this likelihood's marginalisation.

        ``DEVELOPMENT_PLAN.md`` §4.4: "gradient-free samplers cannot
        realistically handle hundreds of latent values, so non-Gaussian +
        flexible-GP robustness is effectively a modern-backend capability". A
        silent downgrade here would be a fit that runs, converges to something,
        and is wrong.

        Pass ``observed`` when the data are to hand — W1.7's ``Dataset`` always
        has them — so that a censored sample the mask excludes is not counted
        against a gradient-free engine.
        """
        required = self.marginalisation if observed is None else self.marginalisation_for(observed)
        if required is Marginalisation.ANALYTIC or differentiable:
            return
        raise LikelihoodError(
            f"{engine} cannot run this likelihood: the {self._family.NAME} family with a "
            f"{type(self._noise).__name__} noise model introduces N latent values that inference "
            f"must marginalise numerically, one per retained sample. That needs a "
            f"gradient-based engine (HMC/NUTS or VI on the torch/jax backends). Options: use "
            f"IndependentNoise with this engine, or a Gaussian family, whose GP marginalises in "
            f"closed form."
        )

    def to_spec(self) -> dict[str, Any]:
        """The declarative description of this likelihood, as plain data.

        Ruled 2026-09-03 (``results.md`` §15 R7, confirmed by W1.13's
        consolidated serialisation review) — the promotion of
        ``ampere.results.describe_likelihood``'s assembly onto the object
        that knows itself: the family name and class, the noise-model class,
        the marginalisation declaration, the parameters' spec
        (:meth:`ParameterSet.to_spec`), and — for a GP — the kernel spec and
        the solver's declaration (name, class, exactness, dataclass
        configuration: ``DenseGP``'s ``jitter`` changes the number the same
        θ scores, so it is part of the identity). Everything is JSON-able.

        **A spec describes the declaration; per-sample and bulk content is
        provenance's business.** Censoring appears as its counts only, and
        buffers not at all: the code positions of 10⁵ limits and the bytes of
        an opacity table are content, fingerprinted by
        ``ampere.results.provenance`` — which composes this mapping and adds
        the hashes (``describe_likelihood``). Two likelihoods differing only
        in kernel family or solver configuration — invisible to
        ``ParameterSet.to_spec()`` alone — are distinguishable here, which is
        the correctness argument R7 was granted on.
        """
        described: dict[str, Any] = {
            "family": self._family.NAME or type(self._family).__name__,
            "family_class": type(self._family).__name__,
            "noise": type(self._noise).__name__,
            "marginalisation": self.marginalisation.value,
            "parameters": self.parameters.to_spec(),
        }
        if isinstance(self._noise, GaussianProcessNoise):
            described["kernel"] = self._noise.kernel.spec().to_dict()
            solver = self._noise.solver
            solver_described: dict[str, Any] = {
                "name": solver.NAME or type(solver).__name__,
                "class": type(solver).__name__,
                "exact": bool(solver.EXACT),
            }
            if dataclasses.is_dataclass(solver) and not isinstance(solver, type):
                solver_described["config"] = {
                    field.name: getattr(solver, field.name) for field in dataclasses.fields(solver)
                }
            described["solver"] = solver_described
        if self._censoring is not None:
            described["censoring"] = {
                "n_samples": int(self._censoring.n_samples),
                "n_censored": int(self._censoring.n_censored),
            }
        return described

    # -- composition-time checking -------------------------------------------

    def check_alignment(self, predicted: FunctionSamples, observed: FunctionSamples) -> None:
        """Verify predicted and observed describe the same samples. O(N), once.

        This is deliberately not called from :meth:`log_prob`: coordinate
        comparison is a composition-time obligation, matching
        ``results_schema.md`` §10's compile-once/evaluate-many split. W1.7's
        ``Dataset`` calls it when the problem is assembled.
        """
        self._check_kinds(predicted, observed)
        if predicted.shape != observed.shape:
            raise LikelihoodError(
                f"the predicted {type(predicted).__name__} has shape {predicted.shape} but the "
                f"observed one has {observed.shape}. A likelihood compares samples index by "
                f"index; a resampling step belongs in the instrument chain (W1.5)."
            )
        if predicted.unit != observed.unit:
            raise LikelihoodError(
                f"the predicted {type(predicted).__name__} is in {predicted.unit} but the "
                f"observed one is in {observed.unit}. Convert once at composition time with "
                f".to_unit(...); units never reach the hot loop (DEVELOPMENT_PLAN.md §7)."
            )
        # Ruled 2026-09-02 (W1.11 gap I-1): kind equality does not imply
        # comparability — VisibilitySet is legal with real or complex values,
        # so a complex prediction could otherwise be fitted against real
        # amplitudes with no warning anywhere.
        if predicted.values.dtype.kind != observed.values.dtype.kind:
            raise LikelihoodError(
                f"the predicted {type(predicted).__name__} holds "
                f"{'complex' if predicted.values.dtype.kind == 'c' else 'real'} values but the "
                f"observed one holds "
                f"{'complex' if observed.values.dtype.kind == 'c' else 'real'} values. A "
                f"likelihood compares like with like: if you mean to fit amplitudes, take the "
                f"modulus in the instrument chain (W1.5) so both sides are real."
            )
        for left, right in zip(predicted.axes, observed.axes, strict=True):
            if left != right:
                raise LikelihoodError(
                    f"the predicted and observed {left.name!r} axes differ. The likelihood "
                    f"assumes they index the same samples; make the instrument chain (W1.5) "
                    f"produce the observed container's own coordinates."
                )
        if self._censoring is not None:
            self._censoring.check_against(observed)
            if not self._family.SUPPORTS_CENSORING:
                raise LikelihoodError(
                    f"the {self._family.NAME} family cannot consume a censoring declaration."
                )
            # Censoring on correlated noise turns the likelihood into a
            # multivariate-normal orthant probability, which is latent. Whether
            # it actually does depends on whether a limit survives the mask, so
            # this is the earliest point the question can honestly be asked.
            if self.marginalisation_for(observed) is Marginalisation.LATENT:
                self._refuse_unconsumed_latent(self._family, self._noise, censored=True)
        self._noise.check_compatible(self._family, observed)
        self._family.check_observed(observed)
        if observed.values.dtype.kind == "c" and not self._family.ALLOWS_COMPLEX:
            raise LikelihoodError(
                f"the {self._family.NAME} family holds real values, but the observed "
                f"{type(observed).__name__} is complex. Use the complex_gaussian family, or an "
                f"explicit real projection (amplitude, phase, real, imag)."
            )

    def _check_kinds(self, predicted: FunctionSamples, observed: FunctionSamples) -> None:
        related = isinstance(predicted, type(observed)) or isinstance(observed, type(predicted))
        # A subclass may legitimately be compared with its base, but only if it
        # kept the same axis signature; otherwise the axis loop below would be
        # comparing different things under the same name.
        if not related or len(predicted.axes) != len(observed.axes):
            raise LikelihoodError(
                f"the predicted container is a {type(predicted).__name__} but the observed one "
                f"is a {type(observed).__name__}. A likelihood compares like with like; channel "
                f"binding (ModelResult.require) is where a kind mismatch should have been caught."
            )

    # -- evaluation ----------------------------------------------------------

    def log_prob(
        self,
        predicted: FunctionSamples,
        observed: FunctionSamples,
        values: Mapping[str, Any] | Sequence[float] | np.ndarray | None = None,
        *,
        latent: np.ndarray | None = None,
    ) -> float:
        """log p(observed | predicted, parameters).

        **The mask convention, stated once.** Masks arrive through
        ``weights()`` — ``results_schema.md`` §16 asked this contract to pick
        one of ``weights()`` and ``masked_uncertainty()`` and say which, and it
        is ``weights()``. The effective weight is the product of the predicted
        and observed containers' weights, so a sample counts only if both call
        it valid, and because those weights are exactly 0 or 1 the product is
        an inclusion indicator. Retained samples are then **excised**: they are
        the only rows and columns that enter the covariance at all.

        Excision, not the infinite-variance limit, is the rule — for *every*
        noise model, not only the correlated one. ``masked_uncertainty()`` and
        a zero weight are equivalent statements about a **chi-square term**,
        which is exactly what ``results_schema.md`` §7 claims for them; they
        are not equivalent statements about a normalised log-density, and a
        likelihood is a normalised log-density. Each Gaussian term carries a
        ``-log sigma`` alongside its chi-square, so
        ``-1/2 log(2 pi) - log sigma - 1/2 (r/sigma)**2`` tends to ``-inf``
        rather than to zero as ``sigma`` grows. The GP case is the same
        statement with a log-determinant in place of the ``-log sigma``: the
        limit approaches the excised value minus ``1/2 log(2 pi sigma**2)``,
        and so diverges. Since §7 defines a masked sample as carrying *exactly*
        zero information, only excision delivers that, and a 0/1 weight
        multiplying a whole term is excision written arithmetically.

        A fully masked pair returns ``0.0``: no data, no information, no
        contribution — which is the correct limit of the same rule, not a
        special case.

        Masked coordinates are still available as evaluation locations for
        :meth:`conditional`, so a GP-localisation diagnostic can show the
        conditioned mean where the data were excluded.
        """
        if latent is not None and self.marginalisation is Marginalisation.ANALYTIC:
            raise LikelihoodError(
                "latent values were supplied, but this likelihood marginalises analytically and "
                "has none. Check Likelihood.marginalisation before allocating them."
            )
        self._check_shapes(predicted, observed)
        resolved = self.context(values)

        weights = np.asarray(observed.weights()).ravel() * np.asarray(predicted.weights()).ravel()
        retain = weights > 0.0
        if not np.any(retain):
            return 0.0

        observed_values = self._retained(observed, retain, "observed")
        predicted_values = self._retained(predicted, retain, "predicted")
        limits = self._retained_limits(retain)
        coordinates = self._coordinates(observed, retain) if self._noise.CORRELATED else None

        noise = self._noise.noise_params(
            observed,
            retain,
            resolved,
            predicted=predicted_values,
            coordinates=coordinates,
            latent=latent,
            limits=limits,
        )
        return float(self._family.log_prob(predicted_values, observed_values, noise))

    def pointwise_log_prob(
        self,
        predicted: FunctionSamples,
        observed: FunctionSamples,
        values: Mapping[str, Any] | Sequence[float] | np.ndarray | None = None,
    ) -> np.ndarray:
        """Per-retained-sample log-likelihood terms — the §4.4 addition R2 granted.

        Ruled 2026-09-03 (``results.md`` §15 R2): granted at the freeze, **not
        stored by default** — a run's ``log_likelihood`` group stays per
        dataset, and these terms are written only by an explicit call. Two
        decompositions, each with the name ``results.md`` §6 reserves for it:

        * **independent noise** (``"factorised"``): the family's own
          ``log_prob`` evaluated pointwise — exact, and the terms sum to
          :meth:`log_prob`. Works for every family, censoring included, and
          for a user family carrying its own aligned data: each single-sample
          call receives a ``retain`` selecting exactly that sample from the
          full containers, so X-2's excision contract holds per term.
        * **a GP** (``"conditional_loo"``): the leave-one-out conditional
          terms ``log N(y_i | mu_i^{-i}, sigma_i^{2,-i})``, from the same
          Cholesky the marginal likelihood forms
          (:meth:`GPSolver.conditional_loo`). These are what ``arviz.loo``
          consumes; they are a *different* decomposition and deliberately do
          **not** sum to the joint value, which has no per-observation
          factorisation under a GP.

        A latent combination is refused: its per-observation terms are
        conditional on latent values that belong to inference, not to this
        method.
        """
        if self.marginalisation is Marginalisation.LATENT:
            raise LikelihoodError(
                f"the {self._family.NAME} family with a {type(self._noise).__name__} noise model "
                f"is a latent-variable likelihood, whose per-observation terms are conditional "
                f"on latent values inference owns; pointwise_log_prob has no unconditional "
                f"answer to give. Use the per-dataset log_likelihood group instead."
            )
        self._check_shapes(predicted, observed)
        resolved = self.context(values)
        weights = np.asarray(observed.weights()).ravel() * np.asarray(predicted.weights()).ravel()
        retain = weights > 0.0
        if not np.any(retain):
            return np.zeros(0, dtype=DTYPE)

        observed_values = self._retained(observed, retain, "observed")
        predicted_values = self._retained(predicted, retain, "predicted")
        limits = self._retained_limits(retain)

        if self._noise.CORRELATED:
            # The Gaussian family, and -- since W4.2 -- the complex Gaussian
            # one. A complex residual becomes the circular GP's two real
            # columns, and `conditional_loo` returns one term per *sample*
            # summing the two components, which is what `results.md` §6's
            # pointwise group is indexed by.
            coordinates = self._coordinates(observed, retain)
            noise = self._noise.noise_params(
                observed,
                retain,
                resolved,
                predicted=predicted_values,
                coordinates=coordinates,
                limits=limits,
            )
            assert noise.solver is not None and noise.kernel is not None
            assert noise.coordinates is not None
            return noise.solver.conditional_loo(
                noise.kernel,
                noise.coordinates,
                _stacked_components(observed_values - predicted_values),
                noise.variance,
                noise.values,
            )

        indices = np.flatnonzero(retain)
        terms = np.empty(indices.size, dtype=DTYPE)
        for position, index in enumerate(indices):
            single = np.zeros(retain.size, dtype=bool)
            single[index] = True
            noise = self._noise.noise_params(
                observed,
                single,
                resolved,
                predicted=predicted_values[position : position + 1],
                limits=None if limits is None else limits[position : position + 1],
            )
            terms[position] = self._family.log_prob(
                predicted_values[position : position + 1],
                observed_values[position : position + 1],
                noise,
            )
        return terms

    def conditional(
        self,
        predicted: FunctionSamples,
        observed: FunctionSamples,
        values: Mapping[str, Any] | Sequence[float] | np.ndarray | None = None,
        *,
        at: np.ndarray | None = None,
    ) -> GPConditional:
        """The conditioned GP mean and variance — the localisation diagnostic.

        ``DEVELOPMENT_PLAN.md`` §4.8's third diagnostic family: "the
        conditioned GP mean ... already localises where the model is
        deficient". Defaults to the *full* coordinate set including masked
        samples, because "what would the GP have said here?" is exactly the
        question a user asks about a region they excluded.
        """
        if not isinstance(self._noise, GaussianProcessNoise):
            raise LikelihoodError(
                f"a conditioned GP mean needs a GaussianProcessNoise model; this Likelihood has "
                f"a {type(self._noise).__name__}, which induces no correlations to condition."
            )
        self._check_shapes(predicted, observed)
        resolved = self.context(values)
        weights = np.asarray(observed.weights()).ravel() * np.asarray(predicted.weights()).ravel()
        retain = weights > 0.0
        if not np.any(retain):
            raise LikelihoodError(
                "every sample is masked, so there is nothing to condition the GP on."
            )
        coordinates = self._coordinates(observed, retain)
        predicted_values = self._retained(predicted, retain, "predicted")
        residual = self._retained(observed, retain, "observed") - predicted_values
        # Since W4.2 a complex residual conditions as the circular GP's two
        # columns, and the conditioned mean comes back complex: one signed
        # deficiency in the real part and one in the imaginary part is what
        # localises a calibration error on a visibility, where the modulus of
        # the mean would hide which way the error went.
        complex_residual = np.asarray(residual).dtype.kind == "c"
        # The prediction is passed here too (ruled 2026-09-03, X-1), so W1.12's
        # diagnostics see the same effective sigma the fit used.
        sigma = self._noise.sigma(observed, retain, resolved, predicted=predicted_values)
        if sigma is None:
            raise LikelihoodError(
                f"conditioning the GP on residuals needs the observed "
                f"{type(observed).__name__}'s per-sample uncertainties, but it has none. A "
                f"latent-GP combination localises through its latent posterior instead, which "
                f"is inference's output rather than this method's."
            )
        # Default to the *full* axis, masked samples included; otherwise hand
        # `at` to the solver untouched, so one place (``_as_points``) decides
        # what a bare 1-D array of coordinates means.
        target = (
            self._coordinates(observed, np.ones(observed.n_samples, dtype=bool))
            if at is None
            else at
        )
        conditioned = self._noise.solver.condition(
            self._noise.kernel_for(observed),
            coordinates,
            _stacked_components(residual),
            sigma**2,
            resolved,
            at=target,
        )
        if not complex_residual:
            return conditioned
        columns = np.asarray(conditioned.mean, dtype=DTYPE)
        return GPConditional(mean=columns[:, 0] + 1j * columns[:, 1], variance=conditioned.variance)

    # -- internals -----------------------------------------------------------

    def _check_shapes(self, predicted: FunctionSamples, observed: FunctionSamples) -> None:
        if predicted.shape != observed.shape:
            raise LikelihoodError(
                f"the predicted {type(predicted).__name__} has shape {predicted.shape} but the "
                f"observed one has {observed.shape}. Call check_alignment() at composition time "
                f"to catch this before a sampler starts."
            )

    def _retained(self, container: FunctionSamples, retain: np.ndarray, role: str) -> np.ndarray:
        flat = np.asarray(container.values).ravel()[retain]
        if flat.dtype.kind == "c":
            if not self._family.ALLOWS_COMPLEX:
                raise LikelihoodError(
                    f"the {self._family.NAME} family holds real values, but the {role} container "
                    f"is complex. Use the complex_gaussian family."
                )
            complex_values = np.ascontiguousarray(flat, dtype=np.complex128)
            return _check_finite(complex_values, f"{role} values")
        return _as_float64(flat, f"{role} values")

    def _retained_limits(self, retain: np.ndarray) -> np.ndarray | None:
        if self._censoring is None:
            return None
        if self._censoring.n_samples != retain.size:
            raise LikelihoodError(
                f"the censoring declaration covers {self._censoring.n_samples} samples but the "
                f"containers hold {retain.size}."
            )
        return np.asarray(self._censoring.kinds)[retain]

    def _coordinates(self, observed: FunctionSamples, retain: np.ndarray) -> np.ndarray:
        if observed.LAYOUT is not Layout.POINTS:
            raise LikelihoodError(
                f"a correlated noise model needs point-set coordinates, but a "
                f"{type(observed).__name__} has a {observed.LAYOUT.value} layout. Gridded 2D+ "
                f"data are the SVGP / SKI / Vecchia strategy slots (Phase 5)."
            )
        stacked = np.column_stack([np.asarray(axis.values, dtype=DTYPE) for axis in observed.axes])
        return np.ascontiguousarray(stacked[retain])

    def __repr__(self) -> str:
        censored = "" if self._censoring is None else f", {self._censoring!r}"
        return (
            f"Likelihood({type(self._family).__name__}, {type(self._noise).__name__}"
            f"{censored}, marginalisation={self.marginalisation.value})"
        )
