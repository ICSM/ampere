# Modality sketch (e) — IFU cube

Status: **design sketch, dispositioned at the freeze (W1.13)** — every gap and requirement below carries its status line. Part of W1.11. Not a contract; a
composition worked example against the merged W1.4–W1.6 code. Code cited
here is `ampere.core` as merged at this repository's `master`.

`Cube` (`results_schema.md` §12) is the modality where the "no regular-grid
assumption anywhere" stance (`results_schema.md` §4.2) meets its most
concrete test: an IFU cube is exactly the case where the schema's `GRID`
layout earns its keep (memory), and exactly the case where the flexible
likelihood — the package's whole reason for existing — is architecturally
hardest to reach, because `GaussianProcessNoise` is declared
`Layout.POINTS`-only. This sketch works that tension through concretely
rather than leaving it as a citation of `likelihoods.md` §15.4.

## 1. The model and instrument: spatial PSF on a Cube

```python
class DustyDiskCube(Model):
    def __init__(self, x, y, wavelength):
        self.register_buffer("x", x, unit=u.arcsec)
        self.register_buffer("y", y, unit=u.arcsec)
        self.register_buffer("wavelength", wavelength, unit=u.micron)
        self.register_parameter(Parameter("t0", st.uniform(50.0, 500.0), unit=u.K, value=200.0))

    def evaluate(self, **values):
        ctx = self.context(values)
        r = np.hypot(ctx["x"][:, None], ctx["y"][None, :])
        temperature = ctx["t0"] / (1.0 + r)               # toy radial cooling, shape (nx, ny)
        lam = ctx["wavelength"]
        flux = temperature[:, :, None] * lam[None, None, :] ** -1.8
        cube = Cube(ctx["x"] * u.arcsec, ctx["y"] * u.arcsec, ctx["wavelength"] * u.micron, flux * u.Jy)
        return ModelResult(cube)                          # bare container -> DEFAULT_CHANNEL
```

`Cube`'s axis order is `(x, y, spectral_axis)` (`results_schema.md` §12),
so `values.shape == (nx, ny, n_wave)` and per-spaxel spectra are
C-contiguous — the rationale `results_schema.md` §17 Q4 states for the
choice. This sketch is input to that open question, so it is worth
recording that this model's natural construction (`temperature[:, :, None]
* lam[None, None, :]`) fits that layout exactly: nothing in writing this
toy model wanted the FITS `(spectral, x, y)` order instead. That is one
data point in favour of the existing decision, not a refutation of it — a
model dominated by *spatial* operations (the PSF step below) would want the
opposite, which is the trade-off §17 Q4 already names; this sketch does not
find a reason to reopen it.

The instrument is a spatial PSF applied slice-by-slice — kind-preserving
(`Cube` → `Cube`), so it needs no `PRODUCES` override, but it is the first
sketch to hit a `Layout.GRID` transformation's mask obligation:

```python
class SpatialPSF(Transformation):
    """A spatial blur applied slice-by-slice; Cube is Layout.GRID, so mask
    propagation has no 1-D propagate_mask() helper to reuse."""
    ACCEPTS = (Cube,)
    def __init__(self, kernel, **kwargs):
        super().__init__(**kwargs)
        self.register_buffer("kernel", kernel)

    def apply(self, samples, values):
        from scipy.signal import convolve2d
        kernel = self.context(values)["kernel"]
        values_arr = samples.values
        out = np.empty_like(values_arr)
        for k in range(values_arr.shape[-1]):
            out[:, :, k] = convolve2d(values_arr[:, :, k], kernel, mode="same")
        new_mask = None
        if samples.mask is not None:
            bad = samples.mask.any(axis=-1)
            grown = convolve2d(bad.astype(float), (kernel != 0).astype(float), mode="same") > 0
            new_mask = np.repeat(grown[:, :, None], values_arr.shape[-1], axis=-1)
        return Cube(samples.x.quantity(), samples.y.quantity(), samples.spectral_axis.quantity(),
                    out * samples.unit, mask=new_mask)
```

```pycon
>>> instrument = Instrument([SpatialPSF(kernel)], label="ifu")
>>> result = model(t0=200.0)
<ModelResult default: Cube[27]>
>>> predicted = instrument(result)
>>> predicted.shape
(3, 3, 3)
```

The `mask=None` case costs nothing (`transformations.md` §6's "free in the
common case" claim holds for `Layout.GRID` too, when there is no mask to
propagate at all). The masked case is where limitation 13.5 bites — see
gap 2.

## 2. The likelihood: `IndependentNoise` works; the flexible likelihood does not, by design

```pycon
>>> whole_cube_like = Likelihood(GaussianFamily(), IndependentNoise())
>>> whole_cube_like.check_alignment(predicted, observed)
>>> round(whole_cube_like.log_prob(predicted, observed), 4)
-6.1336
```

```pycon
>>> Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.1, 2.0))).log_prob(predicted, observed)
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: a correlated noise model needs point-set coordinates, but a
Cube has a grid layout. Gridded 2D+ data are the SVGP / SKI / Vecchia strategy slots (Phase 5).
```

This is not a bug and not a gap: `likelihoods.md` §15.4 states the
restriction, `DEVELOPMENT_PLAN.md`'s Phase 5 bullet ("approximate GP
strategies for images/IFU ... behind the `NoiseModel` interface") already
names the intended future home, and the error message is exactly the loud,
correctly-attributed refusal `architecture.md` §7's "fast path may be
taken, never required" rule calls for. The question this sketch actually
needed to answer is: **what does misspecification-robust IFU fitting look
like *today*, before Phase 5's gridded GP strategies exist?**

## 3. The workaround that exists today: explode into per-spaxel channels

The only route to a flexible likelihood on IFU data with the merged
contracts is to decompose the cube into one `Spectrum` channel per spaxel —
each `Layout.POINTS`, each eligible for `GaussianProcessNoise` — tying the
kernel hyperparameters across all of them exactly as sketch (b) ties them
across échelle orders, because it is genuinely the same detector's
systematics at every spaxel:

```pycon
>>> amplitude = Parameter("amplitude", st.loguniform(1e-2, 1.0), shared_as="ifu_amp")
>>> length_scale = Parameter("length_scale", st.loguniform(0.5, 10.0), shared_as="ifu_len")
>>> spaxel_likes = {}
>>> for i in range(predicted.shape[0]):
...     for j in range(predicted.shape[1]):
...         spec_pred = Spectrum(predicted.spectral_axis.quantity(), predicted.values[i, j, :] * predicted.unit)
...         spec_obs = Spectrum(observed.spectral_axis.quantity(), observed.values[i, j, :] * observed.unit,
...                              uncertainty=observed.uncertainty[i, j, :] * observed.unit)
...         like = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(amplitude, length_scale)))
...         like.check_alignment(spec_pred, spec_obs)
...         spaxel_likes[f"spaxel_{i}_{j}"] = like
>>> len(spaxel_likes)
9
>>> mapping = ParameterSet.merge({"model": model.parameters, **{k: v.parameters for k, v in spaxel_likes.items()}})
>>> mapping.merged.free_size, mapping.tied_names
(3, ('ifu_amp', 'ifu_len'))
```

This **works** — 3 × 3 spaxels correctly collapse to 3 free parameters
(`t0` plus one shared amplitude and one shared length-scale) — and it is
mechanically identical to sketch (b)'s échelle tying. It is also the source
of the two gaps below: a 3×3 toy cube already needed 9 separate `Likelihood`
Python objects (and, in a real fit, 9 separate `Spectrum` channels on the
model side, 9 `Dataset`s, and 9 instrument bindings); a realistic IFU with
tens of thousands of spaxels needs tens of thousands of each, all
hand-constructed the same way.

## Interface gaps

### Gap 1 — no "plate of datasets": per-spaxel decomposition does not scale ergonomically

*Dispositioned at the freeze (W1.13): **deferred to Phase 5** — the
hierarchical/population implementation phase, which this construct
serves for IFU and population fitting alike. What Phase 1 landed is the
substrate: `Binding.index` routes one element of a plate's array-valued
parameter to one component (H-2, ruled 2026-09-02), so per-member wiring
is expressible today at the cost of writing it out; the ergonomic
plate-of-datasets construct is additive on top and the freeze precludes
it nowhere. Its per-spaxel diagnostics half rides with it.*

**What is missing.** `parameters.md` §9's `Plate` solves exactly this
problem *for parameters*: "N members sharing hyperparameters" expands into
one array-valued parameter rather than requiring N hand-written
`Parameter` objects. There is no analogous construct for **datasets**. The
per-spaxel workaround above requires the caller to construct N_spaxel
`Spectrum` channels (or slice them from a `Cube` by hand, as done here), N
`Dataset`-shaped bindings, and N `Likelihood` objects, tying each one's
kernel hyperparameters by passing around two shared `Parameter` objects —
all before `DatasetCollection` (W1.7) exists to do any of this
automatically. This is mechanically sound (the 9-spaxel example above
proves it) but does not scale as *written code*: nothing amortises the
construction, iteration, or bookkeeping across the plate.

**Why it is a real gap and not just "W1.7 will have a loop".** A loop that
builds N `Likelihood` objects is fine computationally (each is a small,
independent GP solve — O(N_λ³) per spaxel, not O(N_λ · N_spaxel)³, so there
is no hidden scaling catastrophe hiding in the workaround itself). The gap
is that nothing in the contracts *names* "many structurally identical
datasets sharing tied noise-model hyperparameters" as a first-class
pattern, so every backend and every user re-derives the same
`{Parameter(shared_as=...): ...}` boilerplate `Plate` was introduced to
avoid on the parameter side. It also has no natural home for reporting
(W1.8's ArviZ coordinates): `Plate`'s array-valued parameter already gives
`free_labels()` a size-N_spaxel problem for the *plate's own* hyperparameter
count (small — two here), but N_spaxel separate `Likelihood` objects give
W1.8 no equivalent structure to hang a "one point per spaxel" diagnostic
plot (§4.8's GP-localisation family) from.

**Proposed amendment.** Not a change to W1.3–W1.6 — the mechanism they
provide is sufficient, as demonstrated above. This is a design requirement
on W1.7 (see below): `DatasetCollection` should support declaring a
templated dataset (instrument chain + likelihood family + noise-model
kernel) applied across a named plate of channels/coordinates (e.g. every
spaxel of a `Cube`), constructing the per-member `Dataset`s and tying their
noise-model hyperparameters automatically, analogous to how `Plate` expands
into per-member `Parameter`s. `results_schema.md` §17's already-open
"nested result channels" question (routed to W1.7, qualified flat names
like `"spaxel.03.07"`) is the same underlying need on the *channel-naming*
side of this problem and should be resolved together with it, not
independently.

### Gap 2 — mask propagation on `Layout.GRID` has no supported helper

*Dispositioned at the freeze (W1.13): **deferred to Phase 5**, beside the
gridded-container strategies that are its consumers (`transformations.md`
§13.5 already names the limitation; the proposed `propagate_mask_grid` is
an additive helper the freeze precludes nowhere). PSF convolution's
standard-library slot (`transformations.md` §10) should land it — or an
equivalent — when the image steps are implemented.*

***Closed at W5.5*, as that disposition anticipated and by the route it
named: `ampere.core.propagate_mask_grid` landed beside `propagate_mask`,
and `PSFConvolution` — the §10 image slot, filled in the same item — is its
first caller. Two details of the shipped helper differ from the sketch's
proposal below and are worth reading beside it. It takes the kernel's
**half-support in pixels**, per axis, rather than the kernel itself: the
ANY rule cares only which inputs reach an output, and a support is what a
step already knows in order to publish its padding, so passing the array
would have invited two answers to the question of what counts as "reached"
(non-zero? above a threshold?). And it dilates **separably**, axis by axis,
which is the same answer as a rectangular-support convolution at `O(k N)`
rather than `O(N²)` — so the "grow the bad-spaxel mask by the kernel's
support" logic this section predicted every author would reinvent now has
one shared, tested primitive, with the boundary behaviour fixed by
construction rather than by each caller's choice of `mode=`. What the
helper deliberately does **not** cover, and what `transformations.md`
§13.5's amended text now says: a grid step whose influence is not a local
neighbourhood — an arbitrary warp, a non-separable resampling — still has
to flatten its own mask.*

**What is missing.** `transformations.md` limitation 13.5 already states
this precisely: "`propagate_mask`'s influence matrix is `(n_out, n_in)`,
which fits `Layout.POINTS`. A `Layout.GRID` container's mask must be
flattened by the transformation itself." `SpatialPSF.apply` above has to
hand-roll a `scipy.signal.convolve2d`-based mask-growth rule to satisfy the
"a step must not drop a mask" obligation `Transformation.__call__` enforces
(`transformations.md` §6) — there is no ready-made 2D (or general N-D)
equivalent of `propagate_mask`'s ANY-rule matrix multiply for a `Cube` or
`Image` step. This sketch confirms the limitation is not merely
theoretical: any real IFU spatial-response step (a PSF, a distortion
resampling, a bad-spaxel dilation) hits it immediately, and every author of
such a step will independently reinvent the same "grow the bad-spaxel mask
by the kernel's support" logic, with no shared, tested primitive and no
guarantee two independent implementations agree on the ANY rule's exact
edge behaviour (e.g. `mode="same"` boundary handling).

**Proposed amendment.** Add an N-D counterpart to `propagate_mask`
(`ampere.core.transform`), e.g.:

```
propagate_mask_grid(mask, kernel, axis=None) -> np.ndarray | None
```

implementing the same ANY rule as `propagate_mask` (an output element is
masked if any input element it depends on, per the kernel's support, is
masked) via binary dilation (`scipy.ndimage.binary_dilation` or equivalent)
rather than a dense `(n_out, n_in)` matrix, which would be wasteful for a
spatial kernel. `axis=` would let a step apply the rule only across the
spatial axes and broadcast across the spectral one, exactly as
`SpatialPSF` needed to. This keeps `results_schema.md` §16's obligation
("a transformation must propagate masks") satisfiable for `Layout.GRID`
steps with the same ready-made, tested confidence `Layout.POINTS` steps
already have, rather than leaving it as a documented but unsupported
requirement.

### Requirements on W1.7

*Dispositioned at the freeze (W1.13): items 1 and 2 defer with Gap 1
(Phase 5); item 3 asked for nothing beyond it. W1.7 as merged supports the
written-out per-spaxel composition today, which is this sketch's own
finding.*

- **`DatasetCollection` needs the "plate of datasets" construct from gap 1**
  before IFU fitting (or any modality needing many structurally identical
  datasets — population fitting, sketch (f), shares this need) is anything
  but hand-rolled boilerplate. This should be designed alongside, not after,
  the "nested result channels" qualified-naming question `results_schema.md`
  §17 already routes to W1.7.
- **`Dataset` construction for a plate member must still call
  `Likelihood.check_alignment` and negotiate/compile_for per member** (the
  same obligations sketches (a), (b) and (d) name) — a batching construct
  must not skip these checks for the sake of amortising them; batching the
  *construction code*, not the *validation*, is the goal.
- **No requirement beyond gap 1** on capability flags or failure signalling
  specific to this modality — a per-spaxel non-positive-definite covariance
  (`likelihoods.md` §16's obligation to convert that exception to `-inf`
  with a recorded reason) is exactly as relevant here as in any other
  flexible-likelihood dataset, and W1.7's general failure-signalling
  requirement (already named by W1.6, not repeated here) covers it without
  anything IFU-specific.
