# Modality sketch (a) — joint spectrum + photometry (the v1 slice)

Status: **design sketch, dispositioned at the freeze (W1.13)** — every gap and requirement below carries its status line. Part of W1.11. Not a contract; a
composition worked example against the merged W1.3–W1.6 code, written to
find interface gaps before the Phase 1 spec freeze. Code cited here is
`ampere.core` as merged at this repository's `master`; nothing in this
sketch is itself implemented or tested by the conformance suite.

This is the case `DEVELOPMENT_PLAN.md` §5 calls the v1 vertical slice: one
physical model, observed both spectroscopically and photometrically. It is
also the least exotic of the six sketches, so it is where composition
*should* close cleanly — and mostly does. The two gaps below were found by
deliberately pushing on the one place a real SED fit always has more than
one instrument reading the same output: several photometric catalogues.

## 1. The model: one SED, two channels

Per `results_schema.md` §4.2, one channel per *distinguishable sampling* of
the SED, not per physical quantity: a fine grid for the spectrograph, a
coarse grid wide enough for synthetic photometry. Both are `Spectrum` —
échelle-order-style "two channels of the same kind" (sketch (b)), one level
simpler because there is no gap between them to preserve.

```python
class GreyBodySED(Model):
    """One SED, exposed on two channels: a fine grid for the spectrograph,
    a coarse grid for synthetic photometry."""
    def __init__(self):
        self.register_buffer("beta", 1.8)
        self.register_parameter(
            Parameter("temperature", st.uniform(100.0, 900.0), unit=u.K, value=300.0)
        )
        self._templates = {
            "spectrum": Spectrum(np.geomspace(2.0, 20.0, 5) * u.micron, np.zeros(5) * u.Jy),
            "sed": Spectrum(np.geomspace(1.0, 30.0, 6) * u.micron, np.zeros(6) * u.Jy),
        }

    def compile_for(self, requirements):
        for channel, asked in requirements.items():
            if channel in self._templates and "spectral_axis" in asked:
                grid = asked["spectral_axis"].coordinates()
                self._templates[channel] = Spectrum(grid, np.zeros(grid.size) * u.Jy)
        return self

    def evaluate(self, **values):
        out = {}
        for channel, template in self._templates.items():
            lam = template.spectral_axis.values
            flux = values["temperature"] * lam ** -self.buffers["beta"].value
            out[channel] = template.with_values(flux)
        return ModelResult(out)
```

`Model`, `Parameterised.register_buffer`/`register_parameter`, `context`,
`compile_for` and the bare-`self` default are all `transformations.md` §3,
§9 as merged; nothing here is new vocabulary.

## 2. Two instruments

**Spectrograph chain**: kind-preserving, so it needs only the calibration
scale legacy hard-codes as `scaleFac` (`transformations.md` §3).

```python
class CalibrationScale(Transformation):
    ACCEPTS = (Spectrum,)
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.register_parameter(Parameter("scale", st.lognorm(0.05), value=1.0))
    def apply(self, samples, values):
        return samples.with_values(samples.values * self.context(values)["scale"])

spex = Instrument([CalibrationScale(label="spec_cal")], channel="spectrum", label="spex")
```

**Synthetic-photometry chain — contract shape only.** W0.9 (pyphot ≥2
forward-migration) has not landed, so the real filter-curve loading stays
Phase 2. What is sketched here is exactly `transformations.md` §4's own
worked example: filter pivots and a response matrix as **buffers** (never
parameters — `architecture.md` §6's question answers itself), `PRODUCES =
PhotometricPoints`.

```python
class SyntheticPhotometry(Transformation):
    """Contract shape only: real filter curves arrive via pyphot post-W0.9."""
    ACCEPTS = (Spectrum,)
    PRODUCES = PhotometricPoints
    def __init__(self, filters, pivots, response, tabulated_wavelength, **kwargs):
        super().__init__(**kwargs)
        self.register_buffer("pivots", pivots, unit=u.micron)
        self.register_buffer("response", response)
        self._filters = tuple(filters)
        self._tabulated_wavelength = np.asarray(tabulated_wavelength, dtype=float)

    def requirements(self):
        # points=, not intervals+max_step: the response matrix's columns
        # are positionally tied to one specific tabulation (see §4 below).
        return (AxisRequirement("spectral_axis", unit=u.micron,
                                 points=self._tabulated_wavelength),)

    def apply(self, samples, values):
        ctx = self.context(values)
        weights = ctx["response"] / ctx["response"].sum(axis=1, keepdims=True)
        return PhotometricPoints(
            self._filters, ctx["pivots"] * u.micron, weights @ samples.values * samples.unit,
            mask=propagate_mask(samples, weights),
        )

wise = Instrument(
    [CalibrationScale(label="phot_cal"),
     SyntheticPhotometry(["blue", "red"], np.array([1.6, 20.0]), response, filter_grid)],
    channel="sed", label="wise",
)
```

Calibration must come **before** the kind change (`ACCEPTS = (Spectrum,)`
on both steps) — `Instrument` checks this at construction
(`transformations.md` §5) and refuses the reverse order outright.

## 3. Negotiation, once

```pycon
>>> asked = negotiate([wise, spex])
>>> sorted(asked)
['sed', 'spectrum']
>>> compiled = model.compile_for(asked)
>>> result = compiled(temperature=300.0)
<ModelResult spectrum: Spectrum[5], sed: Spectrum[6]>
```

The `sed` channel comes back at exactly 6 points — `filter_grid`'s own
length — because the requirement was `points=`, not `intervals` +
`max_step`. That distinction is worth stating as a design lesson even
though it is not itself a gap: a step whose `apply` reads a
construction-time buffer **positionally** (a response matrix, an RMF/ARF)
must publish `points=` naming the exact tabulation that buffer was built
against. A density-only requirement (`intervals` + `max_step`) only
guarantees *coverage*, and `negotiate`'s union is free to satisfy it with
any grid dense enough — which need not have the same length or spacing as
the buffer's own columns. `transformations.md` §10's standard-library table
already says this for RMF/ARF ("points at the matrix's own tabulated
energies"); this sketch confirms the same rule applies to synthetic
photometry and is worth stating there too (see gap 1 for what happens when
it is followed correctly but the channel is shared).

## 4. Evaluating and comparing

```pycon
>>> phot_pred = wise(result, {"phot_cal.scale": 1.0})
>>> phot_pred.filters.tolist(), np.round(phot_pred.values, 3).tolist()
(['blue', 'red'], [194.089, 3.505])
>>> spec_pred = spex(result, {"spec_cal.scale": 1.0})
>>> spec_pred.n_samples
5
```

Per-dataset likelihoods: a flexible GP on the spectrum (the misspecification
robustness this package exists for — a real spectrograph has a wavelength
solution and continuum-placement systematics no photometric point does),
plain independent noise on the photometry (three or four broadband points
carry no exploitable correlation structure at this resolution):

```pycon
>>> flexible = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(2.0, 3.0)))
>>> flexible.check_alignment(spec_pred, spec_obs)
>>> round(flexible.log_prob(spec_pred, spec_obs), 4)
-13.0582
>>> simple = Likelihood(GaussianFamily(), IndependentNoise())
>>> simple.check_alignment(phot_pred, phot_obs)
>>> round(simple.log_prob(phot_pred, phot_obs), 4)
-5.3583
```

And the joint parameter space a `DatasetCollection` (W1.7) needs is one
`ParameterSet.merge` call away, exactly as `parameters.md` §13 and
`likelihoods.md` §16 describe it:

```pycon
>>> mapping = ParameterSet.merge({
...     "model": model.parameters, "spex": spex.parameters, "wise": wise.parameters,
...     "spectrum_like": flexible.parameters, "phot_like": simple.parameters,
... })
>>> mapping.merged.names
('model.temperature', 'spex.spec_cal.scale', 'wise.phot_cal.scale',
 'spectrum_like.amplitude', 'spectrum_like.length_scale')
```

(`phot_like` contributes no names: `IndependentNoise()` with neither
`scale` nor `jitter` declared has zero free parameters, which is correct,
not a bug — the simple case stays simple, per `likelihoods.md` §5.)

This composition — one model, two named channels, two instrument chains,
two per-dataset likelihoods, one merged parameter space — closes without
needing anything not already in W1.3–W1.6. The two gaps below appear only
once the sketch is pushed past the textbook case, to what every real SED
fit actually does: put more than one photometric catalogue on the same
channel.

## Interface gaps

### Gap 1 — a shared channel's negotiated grid can silently break a step's own buffer alignment

*Dispositioned at the freeze (W1.13) and **ruled by Peter, 2026-09-03**:
the lookup is the approved approach — `Axis.locate` (matching within
`COORDINATE_RTOL`, per the W1.11 review's correction) lands with W2.1,
whose standard resampling/photometry steps are its first consumer. It is
additive, so the freeze precludes nothing; the landing PR carries the
decision-log entry ground rule 9 requires.*

***Closed by W2.1.** `Axis.locate` is specified in `results_schema.md` §5 and
implemented in `ampere/core/results_schema.py`; `COORDINATE_RTOL` moved down to
that module (re-exported from `transform.py`) so the lookup and the union's
collapsing cannot drift apart. `transformations.md` §10 now records the pattern
as required for any step whose buffer is tied to specific `points=`, and the
reference backend's `SyntheticPhotometry` is the first consumer. That step
also carries a required `detector=` convention per filter (photon-counting or
energy), ruled by Peter 2026-09-05 and recorded in `DEVELOPMENT_PLAN.md` §2:
the reference backend is the conformance oracle, so the weighting it computes
is the one every other backend must reproduce. The scenario
below — a second instrument unioning its own points onto a bound channel — is a
live test in `tests/backends/test_reference.py`. Only the resampling half of the
original wording turned out not to need the lookup: `Resample` publishes a
density requirement (`intervals` + `max_step`), builds its weights from whatever
grid arrives, and so has no positional buffer to keep aligned.*

**What breaks.** Real SED fits bind several photometric instruments (2MASS,
WISE, IRAS, ...) to the *same* model channel. Give a second instrument its
own, independently tabulated `points=` requirement on that channel:

```pycon
>>> twomass = Instrument(
...     [SyntheticPhotometry(["J", "H"], np.array([1.25, 1.65]), twomass_response, twomass_grid)],
...     channel="sed", label="twomass",
... )
>>> asked2 = negotiate([wise, twomass, spex])
>>> asked2["sed"]["spectral_axis"]
<AxisRequirement 'spectral_axis' micron 10 point(s)>
>>> compiled2 = model.compile_for(asked2)
>>> result2 = compiled2(temperature=300.0)
>>> result2["sed"].n_samples
10
>>> wise(result2, {"phot_cal.scale": 1.0})
Traceback (most recent call last):
    ...
ValueError: matmul: Input operand 1 has a mismatch in its core dimension 0, ... (size 10 is different from 6)
```

`negotiate`'s union of `points=` requirements is, correctly per
`transformations.md` §7, the *superset* of both instruments' exact
coordinates (6 + 4 distinct values here). The model (correctly, per
limitation 13.3: "`negotiate` does not build containers... the model
builds the container") compiles one 10-point channel satisfying both. But
`wise`'s `SyntheticPhotometry.apply` reads `samples.values` **positionally**
against its own 6-column `response` buffer, assuming the channel it is
handed is exactly its own tabulation. Once a second instrument's points are
unioned in, that assumption is false, and the failure is a bare `matmul`
shape error with no reference to negotiation, channels, or the second
instrument at all — exactly the "confusing shape error inside a likelihood"
`results_schema.md` §16 designed composition-time kind checking to prevent,
now recurring one layer down, inside a step's own buffer.

**Root cause.** There is no supported way for a step to recover, from the
(possibly larger, possibly reordered) compiled container it receives at
`apply()` time, the indices of the subset of coordinates matching the
`points=` it originally published. `Axis` has no lookup method; a step
built around a positional buffer has no choice but to assume the compiled
grid *is* its own tabulation, which is true only when it is the sole
requester of that channel.

**Proposed amendment.** Add a lookup to `Axis`
(`results_schema.md` §2's object table):

```
Axis.locate(values) -> np.ndarray[intp]
```

Index lookup matching within `COORDINATE_RTOL`, not exactly: `union`'s
`_dedupe` collapses coordinates that coincide to within that relative
tolerance (`transformations.md` §7), so when two instruments publish
near-coincident but non-identical points only one representative survives,
and a step's own published value may differ from the surviving coordinate
by up to `COORDINATE_RTOL`. `locate` must therefore match within the same
tolerance, raising a `SchemaError` naming any value with no match inside
it (which would indicate a negotiation defect, not a usage error). A
buffer-matrix step then does:

```python
idx = samples.axis("spectral_axis").locate(self._tabulated_wavelength)
weights @ samples.values[idx]
```

instead of assuming positional alignment. `transformations.md` §10's
standard-library table (RMF/ARF, synthetic photometry) should document this
as the required pattern for any step whose buffer is tied to specific
`points=`, since both entries in that table share exactly this shape.

### Gap 2 — two instruments on one channel can default to the same component label, and the collision is silent

*Dispositioned at the freeze (W1.13): **landed.** W1.7's
`DatasetCollection` refuses duplicate dataset labels with the message this
gap asked for, and the 2026-09-03 ruling on `transformations.md` §15 Q4
went further: when more than one instrument reads a channel, distinct
instrument labels are required, checked at problem composition
(`inference.md` §8).*

**What breaks.** `Instrument.label` defaults to the channel name
(`transformations.md` §5). Two instruments bound to the same channel that
both omit `label=` therefore default to the *same* label:

```pycon
>>> wise_v2 = Instrument(
...     [SyntheticPhotometry(["blue", "red"], np.array([1.6, 20.0]), response, filter_grid)],
...     channel="sed",
... )
>>> twomass_v2 = Instrument(
...     [SyntheticPhotometry(["J", "H"], np.array([1.25, 1.65]), twomass_response, twomass_grid)],
...     channel="sed",
... )
>>> wise_v2.label, twomass_v2.label
('sed', 'sed')
```

`transformations.md` §15 Q4 already names this as an open question
("nothing currently checks that at the `negotiate` step ... Should it?").
This sketch supplies the concrete failure mode the question was missing: if
a `DatasetCollection` (W1.7) builds its component mapping the natural way —
keyed by `Instrument.label`, since §14 names that label as "the component
label to merge its `parameters` under" —

```pycon
>>> components = {}
>>> for instr in (wise_v2, twomass_v2):
...     components[instr.label] = instr.parameters
>>> list(components)
['sed']
```

the second instrument's parameters silently overwrite the first's *before*
`ParameterSet.merge` is ever called — no `CompositionError`, no
`TyingError`, just one instrument's calibration parameter vanishing from
the fit. `merge`'s own qualified-name collision detection (`parameters.md`
§8) cannot catch this, because the collision happens one step earlier, in
however W1.7 chooses to build the dict `merge` is handed.

**Proposed amendment.** This is squarely a requirement on W1.7 rather than
a change to W1.5 or W1.3 (see below), but the shape of the fix is worth
recording here since it was found here: `DatasetCollection` construction
must collect every instrument's (and likelihood's) label across all
datasets it is about to join and raise, before calling
`ParameterSet.merge`, if two are equal on a channel that has more than one
instrument bound to it. Relying on `merge`'s own collision check is not
sufficient, because the loss happens in the dict-building step that
precedes it.

### Requirements on W1.7

*Dispositioned at the freeze (W1.13): **all four discharged by W1.7**
(merged 2026-09-02) — `inference.md` §8 is the record: one `negotiate` and
one `compile_for` per model at `FittingProblem` construction;
`check_alignment` and `check_engine(observed=…)` called; one merge with
label uniqueness refused at `DatasetCollection`; and the channel/label
collision made loud (Gap 2 above).*

- **Must call `negotiate()` and `Model.compile_for()` before the first
  evaluation.** `transformations.md` §14 already flags this: "this contract
  defines them but names no caller." This sketch's composition only
  produces the negotiated 6/5-point grids because something called
  `negotiate` and `compile_for` explicitly; a `Dataset`/`DatasetCollection`
  that skips this silently falls back to whatever fixed grid the model
  declared by default, which may not match the observed data at all (a
  shape mismatch `Likelihood.check_alignment` would then catch, but later
  and less informatively than negotiation would have).
- **Must call `Likelihood.check_alignment(predicted_template, observed)`
  once at `Dataset` construction**, and `Likelihood.check_engine(...,
  observed=...)` when an engine is selected, per `likelihoods.md` §16 —
  passing `observed` is not optional (§9's masking-beats-censoring rule
  depends on it).
- **Must build its joint parameter mapping in one `ParameterSet.merge`
  call** across the model and every instrument/likelihood component
  (`parameters.md` §12.4: merge is not associative), and — per gap 2 above —
  **must validate label uniqueness itself** before that call, since neither
  `merge` nor `negotiate` checks it.
- **Must reject (or auto-disambiguate, loudly) instruments that collide on
  both channel and label** when several catalogues bind one channel — gap 2
  is not hypothetical for the v1 slice; it is the default outcome of the
  most common real fit shape (one SED channel, several photometric
  catalogues) unless every instrument is given an explicit label.
