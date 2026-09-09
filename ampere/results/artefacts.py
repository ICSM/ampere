"""Trained-artefact caching keyed on the problem's own hashes.

``DEVELOPMENT_PLAN.md`` §7's trap, named with evidence ("the recent SBI
caching bugs on master"): an SBI posterior, embedding net or emulator reused
against a problem it was not trained on is silently wrong. This module is
the fix — a content-addressed store whose key is built entirely from hashes
:mod:`ampere.results.provenance` already knows how to compute, so a stale
artefact cannot be served without every one of those hashes agreeing first.

Placement, decided here (W3.5)
-------------------------------
The item asks for this decided against ``results.md`` §8's placement logic;
that section turned out to be "The plotting surface", not a placement table —
the closest thing this repository has to one is
``docs/design/contracts/diagnostics.md`` §6's "Module placement summary"
(Family / Module / New extra? / New heavy dependency?), and that is the logic
applied here, with the mismatch recorded for whoever reads this next.

Applying it: this module needs no new extra and no new heavy dependency —
:mod:`pickle` and :mod:`json` are the standard library, and neither ``sbi``
nor ``torch`` is imported anywhere below (an artefact is opaque bytes to this
module; only the caller that built it, and the caller that later unpickles
it, needs those packages installed). What it *does* need, heavily, is
:mod:`ampere.results.provenance` — :func:`~ampere.results.provenance.
spec_hashes`, :func:`~ampere.results.provenance.hash_of`,
:func:`~ampere.results.provenance.hash_container`,
:func:`~ampere.results.provenance.model_hash` (W3.12; the composition below
used to be inline here) and :func:`~ampere.results.provenance.package_versions`
are the whole of how a key is built — and the training-set format this cache
is meant to sit beside (``ampere.results.training``, W2.8) already lives here,
keyed by the same
``ampere_spec_hash``. ``ampere.inference`` already imports from
``ampere.results`` for exactly this reason (``_sbi.py`` writes training sets
through it), so nothing about backend-neutrality is at stake either way — the
deciding weight is that this module *is* provenance machinery, reusing
almost none of ``ampere.inference.engine``'s machinery (no ``Engine``
subclass, no sampler, no problem-evaluation cache) and all of
``ampere.results.provenance``'s. Hence ``ampere.results.artefacts``, not
``ampere.inference._cache``.

The key, and one deliberate departure from the item's literal text
--------------------------------------------------------------------
The item's text names the key ingredients as "the problem's joint spec hash
(``spec_hashes``' ``"spec"`` entry), the observed-data hash, the encoding
layout hash, the method, the estimator architecture, the budget, the seed,
and the sbi/torch versions" (the dispatcher's brief adds "the rounds", also
included below). Read literally, ``spec_hashes(problem)["spec"]`` is
``hash_of(problem.parameters.to_spec())`` — the *parameter declaration*
alone. ``provenance.py``'s own module docstring says why that is not enough
for this exact purpose, in so many words: "the same hash is the cache key the
plan's §7 trap list asks for [...] That is why :func:`problem_fingerprint`
covers the likelihood family, the solver and the kernel as well as the
parameters — two runs differing only in Matérn-3/2 versus squared-exponential,
or ``DenseGP`` versus ``QuasisepGP``, have *identical* parameter specs, and a
hash that could not tell them apart would happily serve a stale emulator or
SBI posterior." A kernel or solver swap that leaves every parameter's name and
prior unchanged is exactly a "model changed" case the item's own Accept line
requires to be a miss, and ``spec_hashes["spec"]`` alone would not catch it.

So :class:`ArtefactKey` keeps ``prior_hash`` exactly as named (literally
``spec_hashes(problem)["spec"]``) and adds ``model_hash`` beside it, kept
separate on purpose so a miss can name *which* of "prior", "model" or "data"
moved rather than one blended hash that can only say "something did".

**W3.12**: ``model_hash`` was originally built here, inline, as a composition
of :func:`~ampere.results.provenance.model_fingerprint` and
:func:`~ampere.results.provenance.dataset_fingerprint` with their
``"parameters"``/``"observed"`` entries stripped out. That composition is
now :func:`~ampere.results.provenance.model_hash` itself, promoted into
:mod:`ampere.results.provenance` so that
:func:`~ampere.results.training.append_training_set` can check the same
digest a training set was written against, not a second copy of the same
recipe. This module no longer computes a fingerprint of its own anywhere —
``prior_hash``, ``model_hash`` and ``data_hash`` are each one call into
:mod:`ampere.results.provenance`.

Refusals
--------
:func:`artefact_key` takes every ingredient as a required, keyword-only
argument — there is no default that would let a caller build a key with one
left out, which is what "a partial key is never accepted" means here. Nothing
in :class:`ArtefactStore` takes a ``force=`` that skips comparing the stored
sidecar's ingredients to the requested key's: :meth:`ArtefactStore.get`
either finds an artefact whose sidecar's ingredients match the key's, byte
for byte, or it returns ``None`` and the caller trains. A stored artefact
that cannot be trusted — a hand-edited or truncated sidecar, a truncated or
otherwise unpicklable ``.pkl`` — is exactly the same ``None``, with an
:class:`ArtefactCacheWarning` naming what went wrong, never a raised
exception: "a hit that fails to load is a miss, never an error" (the item
text), because the fallback (train again) is always safe and a raised
exception in the middle of what looks like a cache-hit fast path is the one
outcome a caller cannot design around.
"""

from __future__ import annotations

import dataclasses
import json
import pickle
import warnings
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any, TypeVar

from ampere.core.dataset import FittingProblem
from ampere.core.exceptions import ResultsError

from .provenance import (
    hash_container,
    hash_of,
    model_hash,
    package_versions,
    spec_hashes,
)

__all__ = [
    "ARTEFACT_CACHE_SCHEMA_VERSION",
    "ArtefactCacheWarning",
    "ArtefactKey",
    "ArtefactStore",
    "artefact_key",
]

_T = TypeVar("_T")

#: Bumped whenever the *shape* of the sidecar this module writes changes —
#: never when an ingredient's own hash recipe changes (that already moves
#: every digest it touches, which is the whole mechanism). Recorded in the
#: sidecar and not hashed, exactly as ``ampere.results.provenance``'s own
#: schema version is recorded and not hashed.
ARTEFACT_CACHE_SCHEMA_VERSION = 1

#: The two packages the item text names by name. A caller training an
#: emulator rather than an SBI posterior may pass a different tuple to
#: :func:`artefact_key`; the default is what every ``SBIEngine`` run needs.
_DEFAULT_PACKAGES: tuple[str, ...] = ("sbi", "torch")


class ArtefactCacheWarning(UserWarning):
    """A stored artefact could not be trusted, so this call trained instead.

    Raised — never an exception — when a store entry addressed by a key's own
    digest turns out not to back that key after all: a sidecar whose recorded
    ingredients disagree with the key just built (hand-edited, or written by
    an older, incompatible version of this module), or a ``.pkl`` that fails
    to unpickle (truncated, corrupted, or written by an interpreter that
    cannot read this one's pickle protocol). Plan §7's rule for a stale cache
    is "invalidate automatically", not "except this once", so every one of
    these is a silent miss with this warning attached, and the caller —
    typically :meth:`ArtefactStore.train_or_load` — proceeds exactly as if
    nothing had ever been stored under this key.
    """


@dataclasses.dataclass(frozen=True, slots=True)
class ArtefactKey:
    """Every ingredient a trained-artefact cache key hashes, named apart.

    Built by :func:`artefact_key`, never by hand — the module docstring
    above records why each field exists and why ``prior_hash`` and
    ``model_hash`` are two fields rather than one. Kept apart deliberately:
    :attr:`ArtefactStore.diff` compares two keys' :meth:`ingredients`
    field by field, and a caller reading "``model_hash`` changed" learns
    something a single combined digest could not have told them.
    """

    #: ``ampere.results.provenance.spec_hashes(problem)["spec"]`` — the
    #: problem's own merged parameter/prior declaration, in the order
    #: ``lowering.md`` §9.2 makes load-bearing. Named literally as the item
    #: text's "the problem's joint spec hash".
    prior_hash: str
    #: The model/likelihood identity the prior hash does not cover — see the
    #: module docstring's "one deliberate departure" section.
    model_hash: str
    #: The same recipe :func:`~ampere.results.provenance.provenance_attrs`
    #: already uses for ``ampere_data_hash`` — reused, not reinvented.
    data_hash: str
    #: The encoding layout name, or (once W3.3 lands) its hash. ``"flat"``
    #: until then, per :data:`ampere.inference.SUMMARY_LAYOUT`.
    layout: str
    #: ``"npe"``, ``"nle"`` or ``"nre"``.
    method: str
    #: The estimator's architecture name (``"maf"``, ``"resnet"``, …), or a
    #: user-supplied builder's ``__name__``.
    architecture: str
    #: Simulations per round.
    budget: int
    #: Rounds of simulate-and-train.
    rounds: int
    #: ``problem.seed``, or the literal string ``"entropy"`` for an unseeded
    #: problem. An unseeded problem is still cacheable — the same declaration
    #: still maps to the same slot — but the "identical posterior" guarantee
    #: the item's Accept line makes is for a seeded one specifically.
    seed: int | str
    #: ``{"python": ..., "sbi": ..., "torch": ...}`` from
    #: :func:`~ampere.results.provenance.package_versions` — an artefact
    #: trained under one library version and reused under another is exactly
    #: the silent staleness plan §7 warns about.
    versions: Mapping[str, str]
    #: **W3.4's carried gap, closed at W3.12.** TMNRE's marginal order — ``1``
    #: for 1-D marginals only, ``2`` to add the pairs — changes what is
    #: trained, so two TMNRE runs differing only here must not share a cache
    #: key. ``None`` for every other method, which has no marginals to
    #: choose; last field, so every pre-W3.4 :meth:`ingredients` digest for a
    #: non-TMNRE key is unchanged (see below).
    marginals: int | None = None
    #: **W3.4's carried gap, closed at W3.12.** The truncation fraction ε that
    #: shapes each round's restricted-prior box — a different ε gives a
    #: different box and so a different trained estimator. ``None`` for every
    #: other method. Last field, for the same reason ``marginals`` is.
    truncation_epsilon: float | None = None

    def ingredients(self) -> dict[str, Any]:
        """Every field as a JSON-safe mapping — what the sidecar records.

        ``marginals`` and ``truncation_epsilon`` are written **only when
        set** (i.e. for a TMNRE key): a non-TMNRE key's ingredients — and so
        its :meth:`digest` — are byte-for-byte what they were before W3.4's
        two fields existed, since the two keys simply do not appear in the
        mapping rather than appearing as ``null``.
        """
        ingredients: dict[str, Any] = {
            "prior_hash": self.prior_hash,
            "model_hash": self.model_hash,
            "data_hash": self.data_hash,
            "layout": self.layout,
            "method": self.method,
            "architecture": self.architecture,
            "budget": int(self.budget),
            "rounds": int(self.rounds),
            "seed": self.seed,
            "versions": dict(self.versions),
        }
        if self.marginals is not None:
            ingredients["marginals"] = int(self.marginals)
        if self.truncation_epsilon is not None:
            ingredients["truncation_epsilon"] = float(self.truncation_epsilon)
        return ingredients

    def digest(self) -> str:
        """The cache key proper — one hash of every ingredient together.

        This is what the artefact and its sidecar are filed under
        (:meth:`ArtefactStore._paths`), so two keys differing in *any* field
        never collide on disk, and there is nothing short of "every field
        agrees" that reads as a hit.
        """
        return hash_of(self.ingredients())


def artefact_key(
    problem: FittingProblem,
    *,
    layout: str,
    method: str,
    architecture: str,
    budget: int,
    rounds: int,
    packages: Sequence[str] = _DEFAULT_PACKAGES,
    marginals: int | None = None,
    truncation_epsilon: float | None = None,
) -> ArtefactKey:
    """Build the one complete :class:`ArtefactKey` for *problem* and a run's settings.

    Every argument beyond *problem* is required and keyword-only: "a partial
    key is never accepted" (the item text) is enforced here, at the one place
    a key is built, so a caller cannot construct one that silently omits an
    ingredient and compares equal to a run that differed in it. The two
    exceptions are *marginals* and *truncation_epsilon*, which are optional
    because they mean nothing outside TMNRE — see their own parameters below.

    Parameters
    ----------
    problem
        The composed problem the artefact would be trained against. Its
        prior, model and data hashes are read from
        :mod:`ampere.results.provenance`, never recomputed by hand.
    layout
        The encoding layout name or hash the estimator's summary vector
        used — ``"flat"`` until W3.3 lands the coordinate-value-mask
        encoding.
    method
        ``"npe"``, ``"nle"`` or ``"nre"``.
    architecture
        The estimator's architecture name, or a builder's ``__name__``.
    budget
        Simulations per round; must be at least 1.
    rounds
        Rounds of simulate-and-train; must be at least 1.
    packages
        Which installed package versions join the key — the default is the
        two the item names, ``("sbi", "torch")``; an emulator cache might
        pass a different tuple.
    marginals
        **W3.4's carried gap, closed here.** TMNRE's marginal order, ``1`` or
        ``2``; ``None`` (the default) for every other method, which has no
        marginals to choose. Written into :meth:`ArtefactKey.ingredients`
        only when given, so a non-TMNRE key's digest is unaffected by this
        parameter existing at all.
    truncation_epsilon
        **W3.4's carried gap, closed here.** TMNRE's truncation fraction,
        strictly between 0 and 1; ``None`` for every other method. A
        different ε gives a different restricted-prior box and so a
        different trained estimator, which must not share a cache key with
        another ε's.

    Raises
    ------
    ampere.core.exceptions.ResultsError
        If *layout*, *method* or *architecture* is empty, if *budget* or
        *rounds* is less than 1, if *marginals* is given and not ``1`` or
        ``2``, or if *truncation_epsilon* is given and not strictly between 0
        and 1 — the ingredients a caller could otherwise leave meaninglessly
        blank or nonsensical.
    """
    layout = str(layout)
    method = str(method)
    architecture = str(architecture)
    if not layout:
        raise ResultsError(
            "artefact_key needs a non-empty layout= (the encoding it was built for)."
        )
    if not method:
        raise ResultsError("artefact_key needs a non-empty method= ('npe', 'nle' or 'nre').")
    if not architecture:
        raise ResultsError(
            "artefact_key needs a non-empty architecture= (the estimator's own name)."
        )
    if int(budget) < 1:
        raise ResultsError(f"artefact_key needs a budget of at least 1, got {budget!r}.")
    if int(rounds) < 1:
        raise ResultsError(f"artefact_key needs at least one round, got {rounds!r}.")
    if marginals is not None and int(marginals) not in (1, 2):
        raise ResultsError(
            f"artefact_key's marginals= is TMNRE's marginal order: 1 (1-D marginals only) or 2 "
            f"(the 2-D pairs too), or None for a method with no marginals to choose. Got "
            f"{marginals!r}."
        )
    if truncation_epsilon is not None and not (0 < float(truncation_epsilon) < 1):
        raise ResultsError(
            f"artefact_key's truncation_epsilon= is a fraction of each 1-D marginal's own "
            f"integral, strictly between 0 and 1, or None for a method with no truncation. Got "
            f"{truncation_epsilon!r}."
        )

    data_hashes = {
        label: hash_container(problem.datasets[label].observed) for label in problem.datasets
    }
    return ArtefactKey(
        prior_hash=spec_hashes(problem)["spec"],
        model_hash=model_hash(problem),
        data_hash=hash_of(data_hashes),
        layout=layout,
        method=method,
        architecture=architecture,
        budget=int(budget),
        rounds=int(rounds),
        seed=problem.seed if problem.seed is not None else "entropy",
        versions=package_versions(extra=tuple(packages)),
        marginals=None if marginals is None else int(marginals),
        truncation_epsilon=None if truncation_epsilon is None else float(truncation_epsilon),
    )


class ArtefactStore:
    """A directory of trained artefacts, content-addressed by :class:`ArtefactKey`.

    ``ArtefactStore(root)`` — *root* is any directory outside the repository
    (a tests-only cache uses ``tmp_path``; a real one, a user-chosen path
    such as ``~/.cache/ampere`` — this module has no opinion). Every entry is
    two files, both named by :meth:`ArtefactKey.digest`: ``<digest>.pkl``,
    the pickled artefact, and ``<digest>.json``, a sidecar recording the
    key's own :meth:`~ArtefactKey.ingredients` so a later reader — or a miss
    — can see exactly what produced it without unpickling anything.

    Nothing here is SBI-specific. sbi 0.27's own trained posterior objects —
    verified directly against the installed package — pickle and unpickle
    with plain :mod:`pickle` and round-trip ``potential()`` and ``sample()``
    exactly under a fixed seed; this store never imports ``sbi`` or ``torch``
    to make that true, because pickling is generic and the artefact's own
    classes carry what they need to reconstruct themselves, exactly as any
    other :mod:`pickle` payload does. A future emulator cache (design horizon
    (c), Phase 5 or later) uses the same store unchanged.
    """

    def __init__(self, root: str | Path) -> None:
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def _paths(self, key: ArtefactKey) -> tuple[Path, Path]:
        stem = key.digest()
        return self.root / f"{stem}.pkl", self.root / f"{stem}.json"

    def get(self, key: ArtefactKey) -> Any | None:
        """The artefact stored under *key*, or ``None`` on any kind of miss.

        A miss is: no file at this digest (never trained, or a different key
        entirely); a sidecar that cannot be parsed or whose
        :meth:`~ArtefactKey.ingredients` disagree with *key*'s own (a
        corrupted or hand-edited sidecar — trusting it anyway is exactly the
        stale-cache trap this module exists to close, so it never happens
        even though the digest matched); or a ``.pkl`` that fails to
        unpickle. The last two warn with :class:`ArtefactCacheWarning` and
        are never raised.
        """
        artefact_path, sidecar_path = self._paths(key)
        if not artefact_path.exists() or not sidecar_path.exists():
            return None
        try:
            sidecar = json.loads(sidecar_path.read_text())
        except (OSError, json.JSONDecodeError) as error:
            warnings.warn(
                f"the cache sidecar at {sidecar_path} could not be read ({error!r}); treating "
                f"artefact {key.digest()} as a miss.",
                ArtefactCacheWarning,
                stacklevel=2,
            )
            return None
        if sidecar.get("ingredients") != key.ingredients():
            warnings.warn(
                f"the cache sidecar at {sidecar_path} does not match the key it is filed "
                f"under (hand-edited, or written by an incompatible version of this module); "
                f"treating artefact {key.digest()} as a miss rather than trusting it.",
                ArtefactCacheWarning,
                stacklevel=2,
            )
            return None
        try:
            with artefact_path.open("rb") as handle:
                return pickle.load(handle)
        except Exception as error:  # any unpickling failure is a miss, not a crash
            warnings.warn(
                f"the cached artefact at {artefact_path} could not be loaded ({error!r}); "
                f"treating artefact {key.digest()} as a miss and training again.",
                ArtefactCacheWarning,
                stacklevel=2,
            )
            return None

    def put(
        self, key: ArtefactKey, artefact: Any, *, attrs: Mapping[str, Any] | None = None
    ) -> None:
        """Store *artefact* under *key*, plus a JSON sidecar of the key's ingredients.

        *attrs* are additional, purely informational fields folded into the
        sidecar (a training-loss trace, a wall-clock time) — never part of
        the key, never compared on a later :meth:`get`. Written to a
        temporary file and renamed into place, so a crash mid-write leaves
        the previous entry (if any) intact rather than a half-written one
        that would later fail to unpickle and be reported, correctly but
        confusingly, as a corrupted cache.
        """
        artefact_path, sidecar_path = self._paths(key)
        self.root.mkdir(parents=True, exist_ok=True)
        tmp_artefact = artefact_path.with_suffix(artefact_path.suffix + ".tmp")
        with tmp_artefact.open("wb") as handle:
            pickle.dump(artefact, handle)
        tmp_artefact.replace(artefact_path)

        sidecar: dict[str, Any] = {
            "schema_version": ARTEFACT_CACHE_SCHEMA_VERSION,
            "digest": key.digest(),
            "ingredients": key.ingredients(),
        }
        if attrs:
            sidecar["attrs"] = dict(attrs)
        payload = json.dumps(sidecar, indent=2, sort_keys=True)
        tmp_sidecar = sidecar_path.with_suffix(sidecar_path.suffix + ".tmp")
        tmp_sidecar.write_text(payload)
        tmp_sidecar.replace(sidecar_path)

        # The one entry point outside content-addressing: recorded so a miss
        # can be explained (`diff`) even though its own digest never matches
        # a *different* key's path.
        latest_tmp = (self.root / "latest.json").with_suffix(".json.tmp")
        latest_tmp.write_text(payload)
        latest_tmp.replace(self.root / "latest.json")

    def diff(self, key: ArtefactKey) -> dict[str, tuple[Any, Any]]:
        """Which ingredient(s) of *key* differ from this store's last :meth:`put`.

        Compares *key*'s :meth:`~ArtefactKey.ingredients` against
        ``latest.json`` — the sidecar of the most recent :meth:`put` this
        store instance (or a previous one pointed at the same *root*) made —
        field by field. This is "a miss can say which ingredient changed"
        (the item text): call it after a :meth:`get` miss to name what moved.

        Empty means either nothing has ever been stored here, or every
        ingredient agrees with the last write — in which case a miss is
        something :meth:`get` already reported with an
        :class:`ArtefactCacheWarning` (a corrupted sidecar or artefact), not
        a changed ingredient.
        """
        latest_path = self.root / "latest.json"
        if not latest_path.exists():
            return {}
        try:
            latest = json.loads(latest_path.read_text())
        except (OSError, json.JSONDecodeError):
            return {}
        old = latest.get("ingredients", {})
        new = key.ingredients()
        names = sorted(set(old) | set(new))
        return {
            name: (old.get(name), new.get(name)) for name in names if old.get(name) != new.get(name)
        }

    def train_or_load(
        self,
        key: ArtefactKey,
        train: Callable[[], _T],
        *,
        attrs: Mapping[str, Any] | None = None,
    ) -> tuple[_T, bool]:
        """The seam ``SBIEngine(cache=store)`` calls around training.

        A hit: :meth:`get` returns the stored artefact, *train* is never
        called, and this returns ``(artefact, True)``. A miss — key never
        seen, or a corrupted entry :meth:`get` already warned about — calls
        *train*, stores what it returns via :meth:`put`, and returns
        ``(artefact, False)``. There is no ``force=`` parameter: the only way
        to make this retrain is for *key* to actually not match a trustworthy
        stored entry.
        """
        cached = self.get(key)
        if cached is not None:
            return cached, True
        artefact = train()
        self.put(key, artefact, attrs=attrs)
        return artefact, False
