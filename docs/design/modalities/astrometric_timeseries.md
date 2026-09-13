# Modality sketch (d) — astrometric time series

Status: **design sketch, dispositioned at the freeze (W1.13)** — every gap and requirement below carries its status line. Part of W1.11. Not a contract; a
composition worked example against the merged W1.3–W1.6 code. Code cited
here is `ampere.core` as merged at this repository's `master`.

Astrometric time series (positional wobble from an unseen companion,
parallax + proper motion fitting, timing residuals) is one of the two
modalities `DEVELOPMENT_PLAN.md` Phase 4 names as following the
interferometry template once that lands. This sketch checks the template
fits *before* that — specifically, whether `TimeSeries` plus the
already-merged `Likelihood`/`NoiseModel` contract can express the one thing
astrometry is chosen for in `DEVELOPMENT_PLAN.md` §5's target-scope list:
correlated timing residuals as a genuine time-domain flexible likelihood,
not just an independent-noise fit with a fancy container.

*(Amended W4.8: landed at W4.9, by the interferometry template of
`interferometry.rst` rather than as a fresh sketch-to-code translation —
`ampere.backends.reference.astrometry.ReflexOrbit` and `.EpochSample` are
this sketch's `ReflexOrbit`/`EpochSampling` under the shipped naming
convention, essentially unchanged in shape. The joint 2-vector GP gap below
is exactly as deferred: W4.9 ships two independent `QuasisepGP`s, one per
coordinate, not the coregionalised pair this section asks for. See
`docs/source/astrometry.rst` for what W4.9 found that this sketch could not
have anticipated — chiefly, that a periodic model needs its own discussion
of prior width, which is not a composition question at all.)*

## 1. The model: a reflex orbit, on two channels

`results_schema.md` §15.2 rules that a vector-valued observable (Stokes
I/Q/U/V, a model's flux and its optical depth) is **separate channels, not
extra columns** — "one value array per container." A 2D sky position is the
same shape of decision: RA and Dec are two `TimeSeries` channels, not one
`TimeSeries` with a `(N, 2)`-shaped value array (which `FunctionSamples`
does not support in any case — every kind here is single-valued).

```python
class ReflexOrbit(Model):
    """Linear proper motion + a periodic reflex wobble, per coordinate."""
    def __init__(self):
        self.register_parameter(Parameter("pmra", st.norm(0.0, 5.0), unit=u.mas / u.yr, value=1.0))
        self.register_parameter(Parameter("pmdec", st.norm(0.0, 5.0), unit=u.mas / u.yr, value=-0.5))
        self.register_parameter(Parameter("period", st.loguniform(50.0, 2000.0), unit=u.day, value=400.0))
        self.register_parameter(Parameter("phase", st.uniform(0.0, 2 * np.pi), value=0.5))
        self.register_parameter(Parameter("amp_ra", st.uniform(0.0, 2.0), unit=u.mas, value=0.5))
        self.register_parameter(Parameter("amp_dec", st.uniform(0.0, 2.0), unit=u.mas, value=0.3))
        self._templates = {
            "ra": TimeSeries(np.array([0.0, 200.0, 400.0, 600.0]) * u.day, np.zeros(4) * u.mas),
            "dec": TimeSeries(np.array([0.0, 200.0, 400.0, 600.0]) * u.day, np.zeros(4) * u.mas),
        }

    def compile_for(self, requirements):
        for channel, asked in requirements.items():
            if channel in self._templates and "time" in asked:
                grid = asked["time"].coordinates()
                self._templates[channel] = TimeSeries(grid, np.zeros(grid.size) * u.mas)
        return self

    def evaluate(self, **values):
        t = self._templates["ra"].time.values
        phase = 2 * np.pi * t / values["period"] + values["phase"]
        return ModelResult({
            "ra": self._templates["ra"].with_values(
                values["pmra"] * t / 365.25 + values["amp_ra"] * np.sin(phase)
            ),
            "dec": self._templates["dec"].with_values(
                values["pmdec"] * t / 365.25 + values["amp_dec"] * np.cos(phase)
            ),
        })
```

`period` and `phase` are shared, at the language level, by construction —
both channels read the same `values["period"]`/`values["phase"]` inside one
`evaluate()` call. That is the "one model, several channels" pattern
(sketches (a), (b)); the physical tying here needs no `shared_as` or `Tie`
at all, because there is only one `ParameterSet` in the picture.

## 2. The instrument: epoch sampling, and nothing else

`transformations.md` §10's standard-library table names exactly this row:
"Epoch sampling | `TimeSeries` → `TimeSeries` | none | `points` at the
observed epochs." It needs no parameters — the observation times are
buffers-by-construction (they are the requirement itself, not the
transformation's payload):

```python
class EpochSampling(Transformation):
    """Pin the model onto the observed epochs; no free parameters."""
    ACCEPTS = (TimeSeries,)
    def __init__(self, epochs, **kwargs):
        super().__init__(**kwargs)
        self._epochs = np.asarray(epochs, dtype=float)
    def requirements(self):
        return (AxisRequirement("time", unit=u.day, points=self._epochs),)
    def apply(self, samples, values):
        return samples  # already on the requested epochs once compiled

epochs = np.array([0.0, 200.0, 400.0, 600.0])
ra_instrument = Instrument([EpochSampling(epochs, label="epochs_ra")], channel="ra", label="astrom_ra")
dec_instrument = Instrument([EpochSampling(epochs, label="epochs_dec")], channel="dec", label="astrom_dec")
```

```pycon
>>> asked = negotiate([ra_instrument, dec_instrument])
>>> compiled = model.compile_for(asked)
>>> result = compiled(pmra=1.0, pmdec=-0.5, period=400.0, phase=0.5, amp_ra=0.5, amp_dec=0.3)
<ModelResult ra: TimeSeries[4], dec: TimeSeries[4]>
>>> ra_pred = ra_instrument(result)
>>> np.round(ra_pred.values, 3).tolist()
[0.24, 0.308, 1.335, 1.403]
```

`AxisRequirement`'s two independent halves (`transformations.md` §7) read
naturally here: astrometric epochs are irregular by nature — telescope
scheduling, weather, target visibility windows — so `points=` (exact
coordinates) is the correct requirement, not `intervals` + a cadence. This
is the same requirement shape sketch (a)'s gap 1 is about; the astrometric
case does not hit that gap, because each channel has exactly one
requester.

## 3. The likelihood: time-domain flexible noise, per coordinate

This is the composition's real target: a `TimeSeries` residual process with
its own correlation structure — unmodelled short-period orbital terms,
DCR-like systematics, guide-star jitter — is precisely the misspecification
case `DEVELOPMENT_PLAN.md` §1 names ampere for, now in the time domain
rather than the wavelength domain. `TimeSeries` is `Layout.POINTS`
(`results_schema.md` §4), so nothing in `likelihoods.md` needs to change at
all: `GaussianProcessNoise`'s `Layout.POINTS`-only restriction
(`likelihoods.md` §15.4) is not a restriction here, it is exactly the right
shape.

```pycon
>>> ra_like = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.05, 100.0)))
>>> dec_like = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.05, 100.0)))
>>> ra_like.check_alignment(ra_pred, ra_obs)
>>> round(ra_like.log_prob(ra_pred, ra_obs), 4)
6.9282
```

The composition — one model, two `TimeSeries` channels sharing orbital
parameters by construction, two epoch-sampling instrument chains, two
time-domain flexible likelihoods — closes without needing anything new.
Stress-testing it further, in the direction real astrometric pipelines
actually need, is where the gap below appears.

## Interface gaps

### Gap — no cross-channel (vector-valued) correlated noise model

*Dispositioned at the freeze (W1.13): **deferred to Phase 5**, where the
approximate/advanced GP strategies live — the `JointGP` slot proposed
below rides with that work, and the freeze precludes nothing (it is an
additive solver strategy plus a noise model spanning datasets, both
extension points the contracts already reserve). The limitation itself is
now recorded in `likelihoods.md` §15, as the amendment asked.*

**What is missing.** Real astrometric residuals are not two independent 1D
noise processes. A guide-star jitter, an uncorrected chromatic term, or an
unmodelled short-period companion all perturb **both** coordinates at the
same epoch together — the physically natural noise model is a single
correlated process over 2-vectors (position residuals), not two separate
scalar processes that happen to share hyperparameters. The same shape
recurs in `results_schema.md` §15.2's own examples: Stokes Q/U with
correlated instrumental leakage, or simultaneous multi-band photometry with
a shared, correlated calibration systematic.

`results_schema.md` §15.2 requires this to be represented as two (or more)
*separate channels* — "one value array per container" — which is correct
at the container level (nothing about a `FunctionSamples` needs to change).
But it pushes the correlation entirely onto the likelihood side, and there
is currently no likelihood-side answer: `Likelihood.log_prob`
(`likelihoods.md` §3–4) takes exactly one `(predicted, observed)` container
pair, and `GaussianProcessNoise._coordinates` (`likelihoods.md` §7) reads
exactly one container's own axes to build its covariance. There is no
`NoiseModel` — not even a declared-but-unimplemented slot, unlike `Rice`
or `VonMisesFamily` — that spans two named channels/datasets and
correlates their residuals. Tying kernel *hyperparameters* across channels
(sketch (b)'s échelle example) gives two GPs with the same amplitude and
length-scale; it does not give one GP over the pair, because the two
processes would still be evaluated, Cholesky-factorised and conditioned
independently. A genuinely joint 2-vector GP (say, one isotropic kernel in
angular-separation-on-sky space, or a intrinsic-coregionalisation-model
kernel with a fitted RA/Dec cross-correlation) cannot be expressed at all
today.

**Why this matters for the freeze.** `DEVELOPMENT_PLAN.md` §4.4's own
future-strategy list (SVGP, SKI, Vecchia for 2D+) is about *scaling* an
already-single-channel GP to gridded or high-dimensional coordinates
(sketch (e)'s territory); it is not about a GP that spans more than one
named channel. This is a genuinely different axis of extension that
neither `likelihoods.md` nor `DEVELOPMENT_PLAN.md` §4.4 currently names.

**Proposed amendment.** Record this explicitly rather than leaving it
implicit:
1. Add a row to `likelihoods.md` §15 (deliberate limitations): "No
   cross-channel/cross-dataset correlated noise. A `NoiseModel` is scoped to
   one container; a physically joint noise process over two or more named
   channels (vector astrometry, Stokes parameters with correlated leakage,
   multi-band photometry with a shared calibration systematic) is not
   expressible. The extension point is a `NoiseModel` bound to a tuple of
   channels rather than one, with a coregionalisation-style kernel — Phase 5
   territory, alongside the other multi-dimensional GP strategies."
2. Name it as a slot in `likelihoods.md` §7's solver-strategy table
   (alongside `InducingPointGP`, `StructuredGridGP`, `VecchiaGP`) so a
   future reader sees it was considered rather than rediscovers the gap:
   e.g. `JointGP` — not exact, applies to 2+ channels, Phase 5.

### Requirements on W1.7

*Dispositioned at the freeze (W1.13): **discharged by W1.7** (merged
2026-09-02) — one model evaluation per draw, distributed to every dataset
(pinned by `tests/core/test_dataset.py`'s evaluation-count test); the
extension point for the vector-GP gap is deferred with the gap itself
(Phase 5, above).*

- **`DatasetCollection` must support several `Dataset`s sharing one `Model`
  evaluation.** RA and Dec are two channels of a *single* `ModelResult`; a
  naive per-`Dataset` architecture that re-evaluates the model once per
  dataset would risk parameter drift between the two evaluations (if, say,
  a stochastic buffer or a not-fully-deterministic `compile_for` were
  involved) and doubles physics evaluation for no reason. W1.7 should state
  explicitly that a `DatasetCollection` evaluates the shared `Model` once
  per parameter draw and hands each `Dataset` its own channel via
  `ModelResult.require(channel, kind)`, rather than leaving "one model,
  many datasets" to arise by accident from independently-written `Dataset`
  objects that each hold their own model instance.
- **`DatasetCollection` should provide the extension point for gap 1
  above**, when it lands: a `NoiseModel`/`Likelihood` scoped to a *group* of
  datasets rather than one is naturally a `DatasetCollection`-level concept
  (it needs to see every dataset bound to the group's channels before it
  can build a joint covariance), not something W1.6's per-dataset
  `Likelihood` can grow into on its own.
