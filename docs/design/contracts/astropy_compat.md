# Ampere v2 — Astropy Interop Contract (W4.6)

Status: **bound by W4.6, 2026-09-11**. A **post-freeze §4 addition** —
`DEVELOPMENT_PLAN.md` §4.7, the last §4 contract to be implemented, made
concrete. Recorded in the plan's decision log the same day. Every sentence
below describes `ampere/core/astropy_compat.py` as it exists; the acceptance
rows are `tests/core/test_astropy_compat.py`,
`tests/core/test_astropy_engines.py` and
`tests/core/test_astropy_backend_hook.py`, and every `>>>` example here is
executed by the last of the three (`TestThisPage`).

Where this document and `DEVELOPMENT_PLAN.md` disagree, the plan wins.

Builds on `parameters.md` (`Parameter`'s three states, `Tie`, priors),
`transformations.md` §4–§6 (`Model`, channels, `compile_for`, the four
capability flags), `results_schema.md` (the container kinds), `inference.md`
§10 and §18 (which engines a capability declaration reaches), and the
2026-09-01 ruling in `DEVELOPMENT_PLAN.md` §2, *"Curated astropy→native
translation: opt-in only, never silent"*.

## 1. What this document is for

`astropy.modeling` is where most of astronomy's analytic models already live,
and a great many people arrive at ampere with one in hand. §4.7's promise is
that such a model "just works": the parameters, their bounds, their fixed flags
and their ties are *already declared* on the astropy side, so ampere translates
that declaration rather than asking for it again.

This contract fixes what that translation is, and — equally — what it is not.
It is **not** a promise of gradients, and it is **not** a licence to substitute
a native lookalike. Both of those are capability claims, and a silently wrong
capability claim is the failure class this architecture exists to prevent.

## 2. The surface

One function and one class, both exported from `ampere.core`:

```
from_astropy(model, *, kind=None, priors=None, channel="default",
             grid=None, output_unit=None, equivalencies=()) -> AdaptedAstropyModel
```

`AdaptedAstropyModel` is an `ampere.core.Model` with one channel. Two helpers
support the backend half (§5): `astropy_components(model)` and
`translation_refusal(backend, model, table)`.

| member | meaning |
|---|---|
| `model` | any `astropy.modeling.Model`, compound models included. **Deep-copied**, never mutated: the adapter writes parameter values into its own copy on every evaluation. |
| `kind` | the `FunctionSamples` subclass the channel emits (§4). Declared by the caller; inferred only where the grid's units settle it. |
| `priors` | `astropy parameter name -> prior`, overriding the translation of §3. |
| `channel` | the emitted channel's name. |
| `grid` | the coordinates to evaluate on, or `None` to take them from `compile_for` (§6). |
| `output_unit`, `equivalencies` | the emitted container's unit, and how to reach it (§4). |

```python
>>> import astropy.units as u, numpy as np, scipy.stats as st
>>> from astropy.modeling.models import PowerLaw1D
>>> from ampere.core import from_astropy
>>> grid = np.array([1.0, 2.0, 4.0]) * u.micron
>>> model = from_astropy(
...     PowerLaw1D(amplitude=2.0 * u.Jy, x_0=1.0 * u.micron, alpha=1.0),
...     grid=grid,
...     priors={"amplitude": st.loguniform(0.1, 10.0), "x_0": 1.0, "alpha": st.norm(1.0, 0.5)},
...     output_unit=u.Jy,
... )
>>> model.kind.__name__, model.parameters.free_names
('Spectrum', ('amplitude', 'alpha'))
>>> model(amplitude=2.0, alpha=1.0)["default"].values.tolist()
[2.0, 1.0, 0.5]

```

## 3. Parameter translation (binding)

Each astropy `Parameter` becomes exactly one of three things, tried in this
order:

| astropy declaration | ampere |
|---|---|
| named in `priors=` | whatever the caller passed — a frozen `scipy.stats` distribution (fitted), a number (held fixed), or a ready-made `Parameter` under the same name. Wins outright, including over `fixed=True`. |
| `tied` (a callable) | **not a parameter**: an `AstropyTie` (§3.2). |
| `fixed = True` | a frozen `Parameter` at its value. |
| finite `bounds` | `Parameter(name, scipy.stats.uniform(lo, hi - lo))`. |
| free, unbounded | **refused, by name.** |

Three properties of this table are the contract:

**3.1 astropy's names are kept exactly.** A compound model declares
`temperature_0` and `temperature_1`, and so does the posterior. Renaming them
would break the tie callables, which read the astropy model *by attribute*.

**3.2 `tied` is a derived quantity, not an ampere `Tie`.** The words collide
and the concepts do not. §4.1's `Tie` collapses several parameter **sites**
into one free parameter: an equality, one sampler dimension, one prior.
astropy's `tied` is an arbitrary Python function of the whole model, so the
quantity it determines is in none of `Parameter`'s three states — not free (no
prior, no dimension), not fixed (its value moves with the fit), not deferred
(no tie group supplies it). Declaring it a `Parameter` would therefore be a lie
whichever state one picked. It is recorded as an `AstropyTie` instead,
recomputed from the astropy model on every evaluation in `param_names` order —
which is the order `astropy.modeling.fitting` applies ties in — and reported in
`describe()` so it reaches provenance. A `priors=` entry for a tied parameter
is refused: nothing would ever consult it.

```python
>>> from astropy.modeling.models import Gaussian1D, Const1D
>>> def twice_the_amplitude(m):
...     return float(m.amplitude_0.value) * 2.0
>>> compound = Gaussian1D(3.0, 6.0, 1.5) + Const1D(0.3)
>>> compound.mean_0.tied = twice_the_amplitude
>>> for name in ("amplitude_0", "stddev_0", "amplitude_1"):
...     getattr(compound, name).bounds = (0.01, 12.0)
>>> from ampere.core import Spectrum
>>> wrapped = from_astropy(compound, grid=grid, kind=Spectrum, output_unit=u.Jy)
>>> wrapped.parameters.free_names          # mean_0 is derived, not sampled
('amplitude_0', 'stddev_0', 'amplitude_1')
>>> sorted(wrapped.ties)
['mean_0']
>>> _ = wrapped(amplitude_0=4.0, stddev_0=1.5, amplitude_1=0.3)
>>> float(wrapped.astropy_model.mean_0.value)     # the tie, applied
8.0
>>> float(compound.mean_0.value)                 # the caller's model, untouched
6.0

```

A Python callable is a black box, so a tie changes nothing about the capability
declaration of §7: an adapted model carrying one is exactly as differentiable
as one without, which is to say not at all.

**3.3 A free, unbounded parameter is refused rather than given a default.**
This is the table's sharpest edge and it is deliberate. `BlackBody.temperature`
ships with `bounds = (0, None)` — half-open, and therefore *not* a uniform
prior over anything. An improper default here would be a fit silently becoming
a different fit. The refusal names the parameter and the three remedies (a
`priors=` entry, finite `bounds` on the astropy model, or `fixed = True`).

A parameter's unit is astropy's own, and `bounds` are always bare numbers in
it — which is exactly ampere's rule that a prior is declared numerically in the
parameter's declared unit, so no conversion is needed or attempted.

## 4. Units (binding)

Both conversions happen **once, at configuration time**, never per evaluation
(`DEVELOPMENT_PLAN.md` §7's units trap).

**4.1 The input grid.** A model that declares `input_units` is handed a
`Quantity`, so astropy performs its own conversion through its
`input_units_equivalencies`. This is not politeness: `BlackBody.input_units` is
**Hz**, so handing it bare micron numbers evaluates a blackbody between 1 and
30 Hz and returns a plausible-looking array. A unitless grid for a model with
declared input units is therefore refused, and so is a unitless grid for any
kind whose axis requires a unit. A model that declares **no** input units is
handed bare numbers **in the grid's own unit** — astropy's unitless models
produce nonsense units from a `Quantity` input rather than refusing, so the
adapter does not give them one.

**4.2 The output.** The astropy model's output unit is probed once; the factor
taking it to `output_unit` is computed once against the grid, checked to be
linear in the value (so that hoisting it is exact — a magnitude or decibel
equivalency is refused by name rather than silently applied per draw), and then
multiplied in. Exactly two equivalencies are ever in force: the caller's own,
and `astropy.units.spectral_density` against the spectral axis where the kind
has one, which is the exact axis-determined conversion between flux-density
conventions rather than an assumption about the source. `output_unit=None`
keeps whatever astropy returned; a model with no units at all has its numbers
*declared* to be in `output_unit` rather than converted into it, and that is the
one place the adapter takes the caller's word.

**4.3 No solid angle is ever invented.** `astropy.modeling`'s `BlackBody` emits
a **surface brightness** — everything per steradian. A surface brightness
reaches a flux density only through a solid angle, and the adapter refuses to
supply one:

```python
>>> from astropy.modeling.models import BlackBody
>>> from ampere.core.exceptions import CompositionError
>>> try:
...     from_astropy(
...         BlackBody(temperature=3000.0 * u.K, scale=1.0 * u.Jy / u.sr),
...         grid=grid, output_unit=u.Jy,
...         priors={"temperature": st.uniform(100.0, 9900.0), "scale": st.loguniform(0.1, 10.0)},
...     )
... except CompositionError as exc:
...     print(exc)
from_astropy() cannot express this model's output, which astropy returns in Jy / sr, in
the requested output_unit=Jy. No equivalency in force makes the conversion, and ampere
will not invent one. A surface brightness (anything per steradian — which is everything
astropy's BlackBody emits) becomes a flux density only through a solid angle: state it.
The usual answer is equivalencies=[astropy.units.dimensionless_angles()], which declares
a solid angle of exactly one steradian and is the convention
ampere.backends.reference.BlackBody's dimensionless 'scale' carries; put the real solid
angle in the astropy model's own scale parameter.

```

The convention, stated: `scale=1.0*u.Jy/u.sr` puts astropy's `B_nu` in Jy/sr,
and `equivalencies=[u.dimensionless_angles()]` declares a solid angle of
**exactly one steradian**. That is the same convention
`ampere.backends.reference.BlackBody` carries in its dimensionless `scale` —
the factor that absorbs the solid angle and the distance dilution together — so
the two are the same physics written twice, and
`TestItAgreesWithTheReferenceModels` compares their `log_prob` at the
conformance suite's `exact` tolerance class (1e-12 relative) rather than
approximately. The real solid angle belongs in the astropy model's own `scale`
parameter, where it is fitted or fixed like any other quantity.

## 5. Opt-in translation (binding)

The 2026-09-01 ruling, implemented. `ampere.core.from_astropy` **never**
substitutes anything: it evaluates the caller's actual astropy model, always.
A *native* equivalent is reachable only through a backend-scoped hook —
`ampere.backends.torch.from_astropy`, `ampere.backends.jax.from_astropy` —
which consults a curated table, decomposes a compound model into its leaves,
and **raises for any leaf it has no row for** rather than falling back to the
black box.

Both directions of silence are forbidden, and for symmetric reasons. Falling
back silently hands somebody who asked for a differentiable model one that is
not, discovered several composition steps later when `NUTSEngine` refuses.
Substituting silently hands somebody a model they did not write: a curated
`BlackBody` is not guaranteed numerically identical to astropy's.

The refusal is `CapabilityError` — nothing about the caller's model is
malformed, the black-box route consumes it happily, and one specific path
cannot serve it, which is what that class is for; it is also a
`NotImplementedError`, which is what "this row has not landed yet" means.
(`LoweringError` was the other candidate and was rejected: its constructor is
prior-family-shaped, rendering "prior family 'BlackBody'".)

**W4.6 landed the hook and an empty table. W4.7 fills it**: `BlackBody`,
`PowerLaw1D`, `BrokenPowerLaw1D`, `Polynomial1D`, `Gaussian1D`, `Const1D`, and
their compound sums, products, differences and ratios (astropy's `+ - * /`;
its other two composition operators, `|` and `&`, have no elementwise meaning
as one channel's flux and are refused by name, same as a leaf outside the
table). The physics is written once, backend-neutral, in
`ampere.core.astropy_translations` — a two-method namespace
(`exp`, `where`) plus each backend's own tested `planck_jy`, in the same
spirit `kernels.py`'s `ArrayOps` fixed for the GP kernels (W4.5) — and both
`ampere.backends.torch.astropy` and `ampere.backends.jax.astropy` build their
`TRANSLATIONS` table from that one dict, so "one table serves both backends".
Parameters, priors, frozen-ness, the channel and the kind all come from
`ampere.core.astropy_compat.translate_astropy_parameters` — the same
translation §3's table describes — so a native and a black-box fit of the
same astropy model cannot disagree on what a bound or a fixed value means.

A compound model carrying an astropy `tied=` parameter is refused on *both*
native backends: a tie is an arbitrary Python callable evaluated on plain
floats, so there is no gradient through it, and no way to evaluate it at all
while jax is tracing. The black-box route (§3.2) is unaffected — it still
applies a tie exactly as astropy defines it.

```python
>>> from ampere.core import astropy_components, translation_refusal
>>> [type(part).__name__ for part in astropy_components(compound)]
['Gaussian1D', 'Const1D']
>>> print(translation_refusal("torch", compound, {}))
ampere.backends.torch.from_astropy() has no native translation for Const1D, Gaussian1D.
Translation is opt-in and never silent (DEVELOPMENT_PLAN.md §2, ruled 2026-09-01), so
this refuses rather than falling back to the black-box adapter: you asked for a
differentiable model and would have been given one that is not. Use
ampere.core.from_astropy(), which evaluates your actual astropy model on the reference
path — gradient-free engines and SBI, never NUTS or VI — or extend the table, which
currently holds: nothing yet (W4.7 adds the curated rows).

```

The empty `{}` above is deliberate — it is `translation_refusal`'s own table
argument, not a backend's, so this page's doctest needs neither extra
installed to run. `ampere.backends.torch.TRANSLATIONS` and
`ampere.backends.jax.TRANSLATIONS`, in an environment with the extra, hold
the six curated classes above; `tests/core/test_astropy_backend_hook.py`
holds both to refusing a model outside them (`Sersic1D`, which the table
does not curate), by name, exactly as this page's example does.

## 6. Kinds and the grid (binding)

The adapter builds three of `results_schema.md`'s kinds, which are the ones an
`astropy.modeling` model maps onto without the adapter inventing structure:

| kind | astropy `n_inputs` |
|---|---|
| `Spectrum` | 1 |
| `TimeSeries` | 1 |
| `Image` | 2 |

`PhotometricPoints` needs a filter list the model does not have; `Cube` and
`VisibilitySet` arrive with their own modality rather than through a generic
adapter. `n_outputs` must be 1: one output is one channel, and the adapter will
not choose kinds or names for several.

**The kind is the caller's declaration.** It is inferred only where the grid's
own units settle it unambiguously — a length, frequency or energy axis is a
`Spectrum`, a time axis a `TimeSeries`, two axes an `Image` — and **never**
from the model's class name. Everything else is refused by name, with the list
of what can be built:

```python
>>> try:
...     from_astropy(Gaussian1D(1.0, 2.0, 0.5))
... except CompositionError as exc:
...     print(exc)
from_astropy() cannot infer the ModelResult kind for Gaussian1D: no grid was given, so
there are no axis units to read, and the negotiated grid is not known until
compile_for(). Declare it — from_astropy(model, kind=Spectrum, ...) — choosing from
Spectrum (1 axis), TimeSeries (1 axis), Image (2 axes). The kind is never guessed from
the model's class name: §4.7 makes it the caller's declaration, and the wrong kind is a
fit that runs and is silently wrong.

```

**The grid** is either the caller's (`grid=`) or the negotiated one.
`compile_for` follows the same template pattern the reference backend's
spectral models use (`transformations.md` §14): one container built once per
channel and refilled with `with_values` on every evaluation, so successive
results share their axes by identity and the unit factor is computed once. A
channel whose requirements name *some* of the kind's axes but not all of them,
with nothing in hand for the rest, is `CompositionError` — the loud option
ruled 2026-09-03 for a model that engages with a requirement it cannot honour.
Each axis is registered as a **buffer** under its axis name, so it reaches
provenance; an astropy parameter colliding with an axis name is refused.

## 7. Capabilities (binding)

```
DIFFERENTIABLE = False    BATCHABLE = False    DEVICE = "cpu"    BACKEND = "reference"
```

All four declared rather than inherited, as every piece of the reference path
declares them (W2.12). The consequences, which the docstring and
`docs/source/astropy.rst` (W4.8) both state:

- **reachable**: `EmceeEngine`, `DynestyEngine`, `ZeusEngine` and `SBIEngine`.
  SBI is the engine a wrapped external model exists for, and it is the one this
  adapter most obviously unlocks.
- **not reachable, ever**: `NUTSEngine` and `VIEngine`. There is no gradient
  through a Python callable, so `ampere.core.realise` refuses the problem by
  name. That is the contract, not a gap to be closed later — closing it is what
  §5's opt-in hook is for, and it closes it by using a *different model*.

`BATCHABLE = False` deserves its precision. astropy models **are** vectorised
over their input grid — that is why one evaluation covers a whole spectrum —
but nothing in `astropy.modeling` takes a stack of *parameter* vectors in one
call, and `BATCHABLE` is a claim about θ (`transformations.md` §4).
`simulate_many` therefore parallelises an adapted model by process rather than
by vectorisation, which is what a black box wants; the adapter pickles for
exactly that reason, and `tests/core/test_astropy_engines.py` runs an
`SBIEngine` budget through a `ProcessExecutor` to prove it.

`describe()` reports the astropy class, its leaf components, its parameter
names, its ties, the kind, the channel and the unit convention — everything
that changes what the model computes and is neither a parameter nor a buffer,
so that two fits of two different astropy models cannot share a cache key
(`results.md` §13.13's limitation 13).

## 8. Limitations (guidance)

1. **Astropy model sets** (`n_models > 1`) are refused: a channel is one
   model's output. Wrap each member separately.
2. **Multi-output models** are refused, as §6 says.
3. **Array-valued astropy parameters** translate (the ampere `Parameter` takes
   the shape), but the prior applies i.i.d. to each element, which is
   `parameters.md`'s general rule and not always what a polynomial's
   coefficients want.
4. **A tie reading another tied parameter** sees that one's value from this
   same evaluation only if it comes earlier in `param_names`. This is astropy's
   own ordering rule, inherited rather than improved on; ampere does not
   topologically sort astropy's ties.
5. **The output unit may not depend on parameter values.** It is probed once
   and checked cheaply (by identity of the unit, not by conversion) on every
   evaluation; a model that changes output unit with θ is refused when it does.
6. **The curated table is empty** until W4.7.
