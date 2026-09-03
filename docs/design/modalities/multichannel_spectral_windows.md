# Modality sketch (b) — multi-channel low-res SED + CO-windows, and échelle orders

Status: **design sketch, dispositioned at the freeze (W1.13)** — every gap and requirement below carries its status line. Part of W1.11. Not a contract; a
composition worked example against the merged W1.4–W1.6 code. Code cited
here is `ampere.core` as merged at this repository's `master`.

This sketch builds on the case that is *already* the acceptance criterion
for W1.4 and the headline worked example for W1.5: a dusty envelope's
broad, coarse SED continuum and its narrow, finely-sampled CO-line windows,
in `results_schema.md` §11 and `transformations.md` §8. That composition is
already implemented and doctested — re-deriving it here would add nothing.
What this sketch adds is the two things the acceptance example does not
cover: **per-channel likelihoods** (§11/§8 stop at the predicted
container) and the **overlapping-échelle-orders** ruling from
`results_schema.md` §17 Q2, which is the same "two channels of the same
kind" pattern in a different instrument and is worth stress-testing in its
own right.

## 1. Recap: the acceptance composition, extended with likelihoods

The model and instruments are exactly `results_schema.md` §11 /
`transformations.md` §8 (see those documents for the full negotiation
worked example — coverage, density union, the preserved gap between the CO
windows). What is new here is closing the loop to a likelihood per channel:

```python
class DustyEnvelope(Model):
    LINES = (866.96, 1300.40)
    def __init__(self):
        self.register_parameter(Parameter("temperature", st.uniform(50.0, 450.0), unit=u.K, value=300.0))
        self.register_parameter(Parameter("line_flux", st.uniform(0.0, 5.0), value=1.0))
        self._templates = {
            "sed_lowres": Spectrum(np.geomspace(1.0, 200.0, 6) * u.micron, np.zeros(6) * u.Jy),
            "co_windows": Spectrum(
                np.concatenate([np.linspace(866.90, 867.00, 3), np.linspace(1300.30, 1300.50, 3)]) * u.micron,
                np.zeros(6) * u.Jy,
            ),
        }
    def evaluate(self, **values):
        out = {}
        for channel, template in self._templates.items():
            lam = template.spectral_axis.values
            flux = values["temperature"] * lam ** -1.8
            for centre in self.LINES:
                flux = flux + values["line_flux"] * np.exp(-0.5 * ((lam - centre) / 0.02) ** 2)
            out[channel] = template.with_values(flux)
        return ModelResult(out)
```

```pycon
>>> model = DustyEnvelope()
>>> result = model(temperature=300.0, line_flux=2.0)
<ModelResult sed_lowres: Spectrum[6], co_windows: Spectrum[6]>
```

**Deliberately different noise models per channel**, chosen for physical
reasons rather than uniformity: the SED continuum has a long correlation
length in log-wavelength (flux calibration and continuum-placement
systematics vary slowly across a decade of wavelength), while the CO
windows — narrow, line-dominated, well-modelled locally — get a much
shorter length-scale tied to the line width itself.

```pycon
>>> sed_like = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(1.0, 30.0)))
>>> co_like = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.05, 0.02)))
>>> sed_like.check_alignment(result["sed_lowres"], sed_obs)
>>> round(sed_like.log_prob(result["sed_lowres"], sed_obs), 4)
-14.3505
>>> co_like.check_alignment(result["co_windows"], co_obs)
>>> round(co_like.log_prob(result["co_windows"], co_obs), 4)
11.8955
```

This is the point `results_schema.md` §14 makes about named channels: a
kind-only schema could not distinguish `sed_lowres` from `co_windows`
(both `Spectrum`), and here that same distinction is exactly what lets each
channel carry an independent, physically appropriate noise-model
declaration. Nothing new was needed — `Likelihood` is already scoped to one
`(predicted, observed)` container pair, so "two channels, two likelihoods"
is the default, not a special case.

## 2. The échelle-order case: two channels of the same kind, one instrument

`results_schema.md` §17 Q2 rules: "overlapping échelle orders are two
channels, which matches how they arrive on the data side (two separate
observed items); no merge helper is owed." That ruling is stated but not
worked through end to end anywhere in the merged specs. Doing so here
stress-tests it against a case the SED/CO-windows example does not cover:
two channels that come from *one physical spectrograph*, so their noise
processes are not merely similar, they are **the same detector's own
systematics**, and a fit that treats them as independent throws away real
information.

```python
class StellarSpectrum(Model):
    LINE = 656.3  # nm
    def __init__(self):
        self.register_parameter(Parameter("teff_index", st.uniform(0.5, 1.5), value=1.0))
        self._templates = {
            "order_43": Spectrum(np.linspace(655.5, 657.0, 4) * u.nm, np.zeros(4) * u.Jy),
            "order_44": Spectrum(np.linspace(656.8, 658.2, 4) * u.nm, np.zeros(4) * u.Jy),
        }
    def evaluate(self, **values):
        out = {}
        for channel, template in self._templates.items():
            lam = template.spectral_axis.values
            flux = np.full_like(lam, values["teff_index"]) - 0.3 * np.exp(-0.5 * ((lam - self.LINE) / 0.3) ** 2)
            out[channel] = template.with_values(flux)
        return ModelResult(out)
```

The two orders overlap in wavelength (655.5–657.0 nm and 656.8–658.2 nm)
and are refused as one `Spectrum` — `Spectrum` requires strictly increasing
coordinates and duplicates are rejected outright (`results_schema.md` §5) —
so §17 Q2's ruling is not optional here, it is the *only* way to represent
this data at all. Each order gets its own blaze-function correction (a
per-order sensitivity curve — a buffer, not a parameter, per
`architecture.md` §6):

```python
class BlazeCorrection(Transformation):
    ACCEPTS = (Spectrum,)
    def __init__(self, blaze, **kwargs):
        super().__init__(**kwargs)
        self.register_buffer("blaze", blaze)
    def apply(self, samples, values):
        return samples.with_values(samples.values * self.context(values)["blaze"])

order43 = Instrument([BlazeCorrection(np.array([0.9, 1.0, 1.0, 0.9]), label="blaze43")],
                      channel="order_43", label="order43")
order44 = Instrument([BlazeCorrection(np.array([0.85, 1.0, 1.0, 0.85]), label="blaze44")],
                      channel="order_44", label="order44")
```

**The interesting part: one spectrograph, one correlated-noise process,
tied across both channels.** Two `GaussianProcessNoise` instances sharing a
*prior* would still be two independent GPs at fit time — the wrong model,
since the whole point is that both orders see the same detector. The fix
is `likelihoods.md` §6's own instruction, applied across channels rather
than within one: kernel hyperparameters are ordinary `Parameter`s, and
`parameters.md` §8's `shared_as` ties them by declaration.

```pycon
>>> amplitude = Parameter("amplitude", st.loguniform(1e-3, 1e-1), shared_as="echelle_amplitude")
>>> length_scale = Parameter("length_scale", st.loguniform(0.05, 2.0), shared_as="echelle_length_scale")
>>> like_43 = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(amplitude, length_scale)))
>>> like_44 = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(amplitude, length_scale)))
```

Note the constraint this ran into and satisfied immediately:
`likelihoods.md` §6 fixes a kernel hyperparameter's *local* name
(`"amplitude"`, `"length_scale"` — they are part of `KernelSpec`, which
W1.9 lowers to term keywords), so the tie label must be carried separately,
via `shared_as`, rather than by renaming the parameter itself. Passing
`Parameter("echelle_amplitude", ..., shared_as=...)` is refused
(`LikelihoodError`, "rename it with `.rename('amplitude')`"); passing
`Parameter("amplitude", ..., shared_as="echelle_amplitude")` is exactly
right, and this is a real, easy-to-hit distinction worth this sketch
recording as a lesson even though it is not itself a gap — the contract's
error message already names the fix.

```pycon
>>> mapping = ParameterSet.merge({
...     "star": star.parameters, "order43_instr": order43.parameters, "order44_instr": order44.parameters,
...     "order43_like": like_43.parameters, "order44_like": like_44.parameters,
... })
>>> mapping.merged.names
('star.teff_index', 'echelle_amplitude', 'echelle_length_scale')
>>> mapping.tied_names
('echelle_amplitude', 'echelle_length_scale')
>>> mapping.merged.free_size
3
```

Three free parameters, not five: `teff_index` plus one shared amplitude and
one shared length-scale, exactly the "one sampler dimension with several
binding sites" `parameters.md` §8 promises for tying. The blaze buffers
contribute nothing to the sampler dimension at all (they are buffers, not
parameters), which is the correct outcome — the calibration curve is known
from the instrument, not fitted.

This composition — two channels of the same kind from one physical
instrument, tied noise-model hyperparameters across them via the *same*
mechanism `parameters.md` uses for tying physical parameters — closes
cleanly. Unlike sketch (a), pushing on it (different channel names per
order, rather than one shared channel) did not surface the label-collision
failure mode found there: each order is its own channel with its own
explicitly-labelled instrument, so there is no shared-channel default-label
collision to hit. The distinction is worth naming precisely because it
shows the échelle case and the "two catalogues, one SED channel" case in
sketch (a) look superficially similar (two things reading "the same
output") but are structurally different: échelle orders are two *channels*
composed by one model, while sketch (a)'s photometric catalogues are two
*instruments* reading one *channel* — and it is specifically the latter
shape that gap 1/2 in sketch (a) attack.

## Interface gaps

None found. This sketch genuinely stress-tested the composition in two
directions — adding per-channel likelihoods with physically distinct noise
models to the already-accepted SED/CO-windows example, and cross-channel
tying of a noise model's own hyperparameters (not just a model's physical
parameters) across two channels of the same kind from one instrument — and
both closed using only mechanisms `parameters.md` §8 and `likelihoods.md`
§6 already specify. The near-miss (a kernel hyperparameter's tie label vs.
its local name) is caught by the existing `LikelihoodError` message, not a
gap.

### Requirements on W1.7

*Dispositioned at the freeze (W1.13): **discharged by W1.7** (merged
2026-09-02) — cross-dataset ties are declared at `FittingProblem`
composition and collapse in the single top-level merge, exercised by
`tests/core/test_dataset.py`'s tied-GP-hyperparameter cases.*

- **A `DatasetCollection` joining `order_43` and `order_44` must perform the
  cross-dataset kernel-hyperparameter tie shown above in its single,
  top-level `ParameterSet.merge` call** — not as two independent
  per-dataset merges later combined, since `parameters.md` §12.4 (merge is
  not associative) rules that out. This is the same requirement sketch (a)
  names for instrument/likelihood components generally; this sketch adds
  that it must hold for noise-model components too, which is not obviously
  the same code path in every plausible W1.7 implementation (an
  implementation that merges "physical model parameters" and "per-dataset
  likelihood parameters" through separate mechanisms would break this
  case).
- **No new requirement beyond what sketch (a) and (d) already name** for
  negotiation ordering or `check_alignment`/`check_engine` calls — this
  composition does not exercise anything beyond a second, independent
  instance of the same pattern.
