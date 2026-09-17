# Ampere v2 — Encoding Contract (W3.3)

Status: **bound by W3.3, 2026-09-09** (drafted the same day by Fable, at
Peter's request, while W3.1 slice 2's gates ran). A **post-freeze §4
addition** — `DEVELOPMENT_PLAN.md` §5's Phase 3 second bullet ("a canonical
coordinate–value–mask tensor encoding of containers for embedding
networks") made concrete. Recorded in the plan's decision log the same day;
**W3.3 landed the code** (`ampere/core/encoding.py`, the wrappers and
`layout=` in `ampere/inference/_sbi.py`) and every sentence below now
describes what exists. Sentences the code proved wrong were amended in
place and are marked ***Amended W3.3***, with the reason; the document is
frozen at the W3.3 merge, except §7, which W3.11 amended in place
(marked ***Amended W3.11***) for two wrapper defaults ruled by Peter on
2026-09-09 — no contract change. §3–§7 are binding; §8–§10 are guidance.

Where this document and `DEVELOPMENT_PLAN.md` disagree, the plan wins.

Builds on `results_schema.md` (containers: `axes`, `values`, `unit`,
`uncertainty`, `mask` with `True` = excluded, `extra_coords`; complex values
only where `ALLOW_COMPLEX`), `inference.md` §13 (`simulate`, `Simulation`,
the batched form and `ContainerBatch`), the W3.2 decision-log row (the
`"flat"` summary, `_summary_of`), and `docs/design/horizon_notes.md` §4–5
and its follow-up (why σ and the coordinate must be features; why the
packing is the contract and the network is not).

## 1. What this document is for

An embedding network in the SBI layer receives **one tensor** per
observation — that is the whole of what the `sbi` package's interface
allows — so everything a network could need to condition on must travel
inside it: where each sample sits (coordinates), what was measured
(values), how well (uncertainties), whether it counts (masks), which
dataset it came from, and, later, under what instrument settings
(context). This contract fixes **the packing**: the rows, the column
groups, their order, their standardisation, and the frozen description
(the *layout*) whose hash says whether two tensors mean the same thing.
It deliberately fixes **nothing about the network**: every embedding —
sbi's set and transformer nets today, neural-process or operator encoders
later — is a module over the unpacked view (§7).

Two consequences are the point. A set-based network over this packing is
invariant to how many samples there are and where they sit, which is what
amortisation across differently-sampled datasets needs. And because
uncertainties are columns, a network can condition on the error bars of
the observation at hand, which is what amortisation across noise
realisations needs — provided the simulator varied them (the observation
context, reserved in `WORK_ITEMS.md`'s Phase 3 section; not this item's).

## 2. Vocabulary

- **Row**: one sample of one dataset's observed container — one
  wavelength of a spectrum, one band of a photometric set, one pixel of
  an image, one baseline of a visibility set. Rows are the unit set
  embeddings are permutation-invariant over.
- **Column group**: a named, contiguous block of columns with one meaning
  (§3). Groups are fixed in order; their widths are fixed by the layout.
- **Layout**: the frozen, hashable description of the packing for one
  fitting problem (§5). Two tensors with the same layout hash are
  comparable; a network trained on one layout refuses another by name.
  ***Amended W3.3***: the class is spelled **`EncodingLayout`**, not
  `Layout`. `ampere.core.Layout` already names the container schema's
  `POINTS`/`GRID` enum and is frozen at `spec-v1.0`, so a second `Layout`
  in the same namespace would be a collision on the one word that has to
  stay unambiguous. `EncodingError` likewise lives in `encoding.py` rather
  than `core/exceptions.py` — this is a post-freeze addition, and its
  vocabulary arrives with it — and subclasses `ContractError`, so a caller
  catching that catches this.
- **Encoded**: the `(rows, columns)` float64 array for one observation, or
  `(count, rows, columns)` for a batch, plus its layout. ***Amended
  W3.3***: it is **always** `(count, rows, columns)`, with `count = 1` for
  one observation, so the observation the posterior is conditioned on and
  the rows a network trains on have the same rank and the same code path —
  which is the reason W3.2's `_summary_of` was already two-dimensional.
  `Encoded.matrix` is the `(count, rows * columns)` view, and is what a
  `"flat"` layout hands `sbi`.
- **Kind of layout**: `"set"` (rows as above; the general case) or
  `"flat"` (W3.2's fixed-size summary — one row, the masked-out values of
  every dataset concatenated in `datasets` order; a degenerate layout
  under the same object so its hash covers it too).

## 3. The packing: column groups, in order

For `kind="set"`, every row carries these groups, in this order. Widths
in brackets; `A` is the layout's coordinate count (§5).

1. **`dataset`** [1] — the dataset's index in `datasets` order (an integer
   stored as float). The network is expected to embed it; the layout maps
   index to label so the tensor is interpretable.
2. **`coordinate`** [A] — the sample's coordinate on each of its
   container's axes, in the container kind's axis order, **standardised
   per dataset and per axis** to `[-1, 1]` over that dataset's observed
   range (§4). A dataset with fewer axes than `A` fills the rest with 0;
   the layout records each dataset's axis count.
3. **`coordinate_features`** [A · 2B] — `sin` and `cos` of the standardised
   coordinate at `B` fixed frequencies `π·2^k`, `k = 0…B−1` (NeRF-style
   Fourier features); `B` is a layout parameter, default 4, and `B = 0`
   removes the group. This is how a network with no positional embedding
   sees *where* a sample is, at more than one scale. Ordered axis-major,
   then band, then `(sin, cos)`. ***Amended W3.3***: an axis a dataset does
   not have contributes **zeros** here rather than the `sin(0) = 0`,
   `cos(0) = 1` its zero-filled coordinate column would give. An axis that
   does not exist must not look like an axis sitting at the centre of its
   range.
4. **`axis_identity`** [A] — ***added W5.11***. One small integer code per
   coordinate column, saying what that column physically *is* for the
   dataset whose rows it sits on: the axis's physical type, read from its
   unit. A dataset with fewer axes than `A` leaves the rest at `0`
   (`absent`), which is a different code from `1` (`unknown`: an axis that
   exists and whose physical type the table cannot name). **The code table
   is contract** — the numbers below are frozen into every layout's hash
   and every trained network's input, so a code may be *added* as the next
   free integer but never renumbered or reused:

   | code | name | what it is |
   |---:|---|---|
   | 0 | `absent` | no axis in this column for this dataset (padding) |
   | 1 | `unknown` | an axis with no unit, or a physical type not below |
   | 2 | `dimensionless` | a dimensionless coordinate |
   | 3 | `length` | wavelength, a spatial coordinate in metres, a baseline |
   | 4 | `frequency` | a spectral axis in Hz |
   | 5 | `energy` | a spectral axis in keV |
   | 6 | `time` | a light curve's time axis |
   | 7 | `angle` | an image's `x`/`y` on the sky |
   | 8 | `spatial frequency` | interferometric `u`, `v` (wavelengths or rad⁻¹) |
   | 9 | `wavenumber` | a spectral axis in cm⁻¹ |
   | 10 | `speed` | a velocity axis |
   | 11 | `temperature` | |
   | 12 | `mass` | |

   The code comes **from the unit**, through
   `astropy.units.get_physical_type`, because that is the one statement of
   what a coordinate is that every container already carries: a `Spectrum`
   in `um` and one in `Hz` hold different quantities in the same column and
   the codes say so. The single exception is `spatial frequency`, which
   astropy cannot name — a radian is dimensionless, so `u` in wavelengths
   is `dimensionless` and `u` in `rad**-1` is `unknown` — and which a
   container kind therefore declares by putting `rad**-1` in its
   `AxisSpec.equivalent_units`; an axis of such a kind whose unit is
   dimensionless-equivalent is coded 8. A `u` given in metres is a baseline
   **length** and is coded 3, which is the honest answer.
5. **`value`** [1 or 2] — the whitened value `y/σ` (real and imaginary
   columns when the container is complex, else one column). Where the
   dataset has no uncertainties, the asinh-scaled value (§4) stands in and
   the `has_sigma` feature (group 9) says so.
6. **`value_asinh`** [1 or 2] — `asinh(y / s)` with `s` the dataset's
   value scale (§4): the magnitude of the measurement on a scale that is
   linear near zero and logarithmic far from it, so a network can tell a
   bright line from a faint one without the raw unit.
7. **`log_sigma`** [1] — `log(σ / s)`; 0 where there is no σ.
8. **`mask`** [1] — 1 for a row that counts, 0 for a row the dataset's
   `effective_mask` excludes **and** for every padded row. There is
   exactly one mask convention in the contract; §7 says how wrappers
   convert it.
9. **`set_features`** [3] — per-dataset quantities broadcast onto every
   row of that dataset: `log(N_valid)` (valid rows in this dataset),
   `has_sigma` (1/0), and `is_complex` (1/0). They exist because a masked
   mean (§7) deliberately removes the row count from the pooled
   embedding, and the count is information.
10. **`context`** [0 today] — the reserved observation-context group.
    Width 0 until the context item lands; its presence at position 10 is
    what makes adding it a layout-hash change and nothing else.
    ***W5.10*** amortises over the observation context **without widening
    it**: the context prior varies the σ-pattern of each simulated draw,
    and the per-row `log_sigma` and whitened-value columns are what the
    network then sees it through. The group stays reserved for a context
    that is *not* already a column.

The values in groups 5–7 are computed from the container the row belongs
to, with masked samples' values **left as the container holds them** and
then zeroed by the mask at unpack time (§7) — the contract does not invent
values for masked samples, mirroring `inference.md` §13's rule for
`simulate`.

***Amended W3.3***, one clause: a masked sample's encoded value is left as
the container holds it *unless the result is not finite*, in which case
the column is **0**. Masking a sample is precisely what a user does about
a NaN, so a masked sample legitimately holds one; and a non-finite number
in the stored tensor poisons the training loss, the netCDF round trip and
`sbi`'s own input validation even where the mask says the row is absent —
zeroing it downstream is too late. Finite values on masked rows are kept,
so the round trip is exact there too. The mirror-image rule is a refusal:
a **retained** sample whose value is not finite, or whose σ is not
strictly positive and finite, is refused by name (§6) rather than
sanitised, because there the number is the data.

Rows are ordered dataset by dataset in `datasets` order, and within a
dataset in the container's storage order (flattened C-order for multi-axis
kinds). Row order carries **no meaning** — a set embedding must be
invariant to it, and the conformance row for the wrapper asserts it — but
it is fixed so that an encoding is reproducible bit for bit.

**Padding.** The layout fixes a row cap `R` (§5). An observation with
fewer rows is padded to `R` with all-zero rows and `mask = 0`; one with
more is **refused by name** (the remedy is a larger cap in a new layout,
or masking). Padding is what lets a batch be one rectangular tensor.

***Amended W3.3***, on where the padding sits and what "fewer rows" can
mean in v1. The padded rows are at the **tail**, after every dataset's
block, and each dataset's row count is itself a field of the layout, so a
container whose sample count differs from the layout's is refused at
encode (§6) rather than encoded into a shifted block. That is what keeps
`per_dataset`'s row bounds (§7) exact: a per-dataset capacity, which is
what a genuinely differently-sampled observation would need, is §9.1's
reserved extension and a field rather than a version bump. So `R` above
the observation's own row count is legal and costs only zero rows today;
it becomes useful when that extension lands.

## 4. Standardisation, and why it lives here

`sbi` standardises its input tensor column-wise from training statistics
by default. On this packing that would z-score the mask, the dataset
index and the coordinate columns, and padded rows would enter every
column's mean. So the engine passes `z_score_x="none"` for any `"set"`
layout and **standardisation is the encoding's**, with every statistic a
field of the layout, computed **once from the observed data at layout
construction** and never from the simulations:

- coordinate range per dataset per axis: `(min, max)` of the observed
  axis values, mapped to `[-1, 1]`;
- value scale `s` per dataset: the median of `|y|` over valid observed
  samples (falling back to the median σ, then 1.0), so `value_asinh` and
  `log_sigma` are dimensionless and O(1);
- σ enters only as `y/σ` and `log(σ/s)`.

Computing from the observation rather than the training set is
deliberate: it makes the layout a function of the *problem*, fixed before
any simulation is drawn, identical at training and at inference, and
independent of the budget. The cost is that a simulation whose values
wander far outside the observed scale is encoded at the tails of `asinh`
— acceptable, and visible, rather than a silent renormalisation.

## 5. The layout

A frozen, plain-data record (`EncodingLayout`) with, at least:

- `kind` (`"set"` | `"flat"`), `version` (`ENCODING_VERSION`, **2 since
  W5.11**);
- `datasets`: for each label in order — `kind` (container class name),
  `axes` (names, count), `axis_codes` (***added W5.11***: one §3 group 4
  code per axis, in the same order), `rows` (valid sample count of the
  observation), `is_complex`, `has_sigma`, the coordinate ranges and the
  value scale;
- `coordinates` (`A`, the maximum axis count), `fourier_bands` (`B`),
  `row_cap` (`R`; default the total observed row count, so the default
  layout fits the observation exactly and a differently-sampled
  observation needs a wider cap chosen on purpose);
- `columns`: the resolved column groups with their offsets and widths
  (derivable, stored for interpretability);
- `hash`: `hash_of(to_dict())` (`ampere.results.provenance.hash_of`),
  over everything above.

***Amended W3.3***, three corrections to that list.

1. **`rows` is the total sample count, not the valid one**, and `valid` is
   a second field beside it. §3 requires a masked sample to be a row with
   `mask = 0` — that is the whole reason the mask is a column — so the
   valid count cannot also be the row count. `row_cap` defaults to the sum
   of `rows`.
2. **The mask itself is in the layout**, as each dataset's `excluded` flat
   indices. It has to be. `Dataset.effective_mask` is resolved lazily, at
   the first prediction, so a driver that encodes the observation before
   simulating and the draws afterwards would read `None` once and a mask
   the next time — two different packings inside one run, under one hash.
   Freezing it at layout construction (from `effective_mask` where it is
   resolved, from the observed container's own mask where it is not) makes
   the guarantee structural. A differently *masked* observation is
   therefore a different layout, and `compare` says so.
3. **`hash_of` is imported inside the method that calls it.**
   `ampere.core` imports no other ampere namespace at module scope, and
   `ampere.results` imports `ampere.core`. The recipe is unchanged — the
   same canonical JSON and the same personalised BLAKE2b digest every
   other ampere fingerprint uses — so a layout hash sits beside a spec hash
   in a run's attrs and means the same kind of thing.

The layout is recorded in a run's attrs (`ampere_encoding_layout`, the
dict; `ampere_encoding_hash`) and in a training set's attrs, is one
ingredient of W3.5's artefact key, and is what `SBIEngine(layout=…)`
selects — `"flat"`, `"set"`, or a `Layout` instance for a widened cap or
a different `B`. **A network trained under one hash refuses an
observation encoded under another, by name, saying which field
differs.**

## 6. Refusals

By name, with the remedy: more rows than `row_cap`; a layout whose
dataset labels, kinds, axis counts or complex-ness differ from the
problem's (a *different* observation of the same shape passes — that is
amortisation; a different *shape* does not); a non-numeric extra
coordinate (ignored with the layout recording that it was); `"flat"` on a
complex container (W3.2's refusal, kept until a use case says otherwise).

***Amended W3.3.*** The landed list, in full, each naming the remedy:

- a `row_cap` below the observation's own row count (raise the cap, or
  mask what should not be fitted) — checked at layout construction, so a
  cap that could never work is refused before any encoding happens;
- a layout whose dataset labels, container kinds, axis names, row counts,
  grid shapes, complex-ness, `has_sigma` or **mask** differ from the
  problem's, reported field by field with the layout's own hash in the
  message;
- a container whose sample count differs from the layout's for that
  dataset (see §3's padding amendment);
- a batch whose datasets disagree on how many draws they carry;
- a dataset the layout knows and the observations lack — a simulation that
  produced no observation is not an empty one;
- a retained sample whose value is not finite, or whose σ is not strictly
  positive and finite (mask it, or declare the dataset without
  uncertainties and let the `asinh` value stand in);
- a problem with no datasets, and one where every sample of every dataset
  is masked out — the encoding would be empty, and there would be nothing
  to condition on;
- `"flat"` on a complex container, **at layout construction** rather than
  at the first encode: complex-ness is a fact about the problem, so it can
  be refused before any observation is packed.

One entry of the drafted list did **not** land: a non-numeric extra
coordinate is ignored, and the layout records *nothing* about it. Recording
which extra coordinates were ignored would put a container's free-form
per-sample labels into the layout hash, so renaming a filter would become a
different layout and therefore a retrained network. §9.5 already says extra
coordinates are not features in v1; silence about them is the consistent
reading, and the container is where they remain visible.

## 7. `unpack`, wrappers, and the one mask convention

`unpack(x, layout) -> Unpacked`: named views of the groups (`dataset`,
`coordinate`, `coordinate_features`, `value`, `value_asinh`, `log_sigma`,
`mask`, `set_features`, `context`), `features` (every non-mask group
concatenated, the tensor a per-row network consumes), `valid` (the
boolean mask), and **`per_dataset`**: the same views sliced by dataset in
`datasets` order (the layout knows each dataset's row offset and count),
with, for a grid kind (`Image`, `Cube`), a `grid_shape` so the slice can
be reshaped back to its axes exactly — rows are stored in C-order for
precisely this reason. `per_dataset` is what lets a hierarchical
embedding (§9.1) be a wrapper rather than a packing change. It is the **only** way a network is meant to read the
tensor, and it is the backend-neutral half; the torch-side wrappers in
`ampere.inference` adapt to what each embedding wants:

- sbi's `PermutationInvariantEmbedding` treats an all-NaN row as absent
  and pools with a masked mean: the wrapper sets padded and masked rows to
  NaN in `features` and passes the rest; no mask column reaches the net.
- sbi's `TransformerEmbedding` takes an `attention_mask`: the wrapper
  passes `valid`, constructs the net with `is_causal=False`,
  `pos_emb="none"` and explicit dropout, and returns the first element of
  its tuple output.
- A user module gets `Unpacked` if it declares it wants one (a class
  attribute), and the raw tensor otherwise.

The cast from float64 to the network's float32 happens once, in the
wrapper. Pooling is the **masked mean**, never the sum, so the pooled
embedding does not scale with row count; the count is `set_features`'.

***Amended W3.3***, from reading what sbi 0.27 actually does. Both nets
work, both are used, and both needed more of the wrapper than the draft
assumed.

- `PermutationInvariantEmbedding` pools with a **sum** by default
  (`aggregation_fn="sum"`), so it is constructed with `"mean"`. Its own
  valid-row count is computed as
  `isnan(x).sum(dim=1).reshape(-1)[:num_batch]`, which reads the *first*
  batch element's count for every element as soon as `x` has more than one
  feature column. Under this contract that is exactly right, because the
  mask is a property of the **layout** and is therefore identical across a
  batch — but it is right by construction rather than by luck, and a
  future packing with a per-draw mask must stop using that aggregation.
- `TransformerEmbedding.forward` **discards `attention_mask` unless
  `is_causal` is true** (`else: attention_mask = None`), and its
  aggregation reads the **last token** rather than pooling. So the mask
  cannot be delegated on the non-causal path the draft (rightly) requires:
  the wrapper zeroes masked and padded rows' *tokens* after its
  projection, which makes the output independent of what a masked row
  holds whatever the net does with the mask. Its `feature_space_dim` is
  the transformer's **model** dimension, not an input width, so the
  wrapper owns a linear projection from the packing's feature width into
  it. `forward` returns a bare tensor in 0.27 despite its `Tuple`
  annotation. (The last-token read, and passing `attention_mask` to
  `forward` at all, are superseded below — *Amended W3.11*.)
- **A user module gets the raw tensor**, and `unpack(x, layout)` is public
  for it to call. The class-attribute protocol was dropped: it adds a
  second way to write an embedding, and a second thing every future
  wrapper has to remember, in exchange for one line the module can write
  itself.
- `unpack` also accepts a **flattened** `(..., rows · columns)` tensor and
  reshapes it, because sbi flattens `x` on some of its paths and no
  wrapper should have to know which one it is on.

***Amended W3.11***, from W3.3's two open questions (ruled by Peter
2026-09-09; both are defaults, not findings — a later embedding study
tests how the right choice depends on the data, model and problem
structure, and writes up the result as guidance).

- **Default output width.** `"flat"`'s CNN/FC embeddings keep the legacy
  default, `2 · free_size`. For `"set"`/`"transformer"` the default is
  now `max(2 · free_size, 32)`: the pooled vector is the whole of the
  conditioning vector's capacity rather than an intermediate width, and
  both shipped nets end in a ReLU, so a narrow, randomly-initialised one
  can emit all zeros before a single gradient step. A user's own
  `output_dim=` always overrides either default.
- **The transformer's readout.** `TransformerEmbedding.forward` reads the
  **last token** as its summary (`self.aggregator(hidden_states[:, -1,
  :])`). Under full attention with no positional embedding that token is
  a function of every row, but *which* row happens to sit last is a
  property of row order, not of the observation — two encodings of one
  set that differ only in row order train and condition on different
  summaries. The wrapper now builds a **masked-mean readout head**
  instead: permutation-invariant, and using every retained token rather
  than one. (A learned CLS token is the alternative — more parameters for
  the same information — and was not chosen.)

  There is no hook for this in sbi 0.27's public interface — `forward`
  returns a bare tensor, not the hidden states it pooled — so the wrapper
  never calls `net.forward` for the transformer at all, and instead calls
  the net's own body modules directly, in the order and with the
  arguments `forward` itself uses before its last-token slice:
  `net.preprocess` (identity, on the non-ViT path this wrapper always
  builds), each block in `net.layers` (called as `block(hidden_states,
  attention_mask=None, position_ids=None)`, taking element 0 of its
  tuple), `net.norm` (the final normalisation), and `net.aggregator` (the
  linear layer to `final_emb_dimension`, applied to the pooled vector
  rather than to the last token); `net.is_causal` is checked and required
  `False`, since the masked mean assumes attention was never causally
  restricted. Because these are read by name rather than through a
  public API, `ampere/inference/_sbi.py`'s `_transformer_masked_mean`
  documents each one and raises `AttributeError` naming whichever is
  missing or wrong; `tests/inference/test_sbi.py`'s
  `test_the_readout_depends_on_named_sbi_attributes` fails loudly first
  if a future sbi renames or restructures any of them.

## 8. Stability

`ENCODING_VERSION` is bumped when a column group is added, removed,
reordered or redefined; the reserved `context` group exists so that the
first real context is a width change under the current version, not a
version bump. ***W5.11*** is the first such bump: adding the
`axis_identity` group took the version from 1 to 2, and **every layout
hash written under version 1 moved once**. That is what the version is
for — a stored layout, training set or trained network from before the
bump is refused *by name*, saying that the axis identity is the
difference, rather than reported as a hash that does not match
(`EncodingLayout.from_dict`, `append_training_set`). The code table
itself is frozen the same way: a code may be added as the next free
integer, never renumbered. The packing is independent of `CONTAINER_SCHEMA_VERSION`: it reads
containers through their public attributes, and a container schema change
that keeps `axes`/`values`/`uncertainty`/`mask` leaves every encoding and
every hash unchanged.

***W5.10 moved no hash at all***, and it is worth saying why, because it is
the clearest statement of what this packing is for. Amortising over the
observation context means simulating each draw at a σ pattern drawn from a
context prior (`inference.md` §13). A drawn observation carries its own
uncertainties, so the *values* of the `log_sigma` and `value` columns vary
across the budget while the *layout* — every field of §5, and therefore the
hash — is unchanged: it is a function of the observed problem, and the
problem did not change. A network trained on that budget reads the context
through columns it was already reading. The reserved `context` group stays
width 0, for a context that is genuinely not already a column.

## 9. Deliberate limitations of this draft

1. **One tensor for the whole collection, and one pooling rule.** A long
   spectrum beside a few photometric points shares one row set, and
   full attention is quadratic in rows; an image (a 256×256 `Image` is
   65 536 rows) makes a set embedding with full attention infeasible and
   even a DeepSets pass memory-heavy at SBI batch sizes. **The packing is
   not the obstacle — the single pooling stage is**, and the extension
   path is fixed here so it is not re-derived (Peter's question,
   2026-09-09). (a) *Per-kind first-stage encoders* over `per_dataset`
   slices: a set or attention encoder for irregular 1D kinds, a CNN or
   ConvCNP over the reshaped grid for `Image`/`Cube` (the right inductive
   bias for a grid anyway, and linear in pixels), an identity for a
   handful of photometric points — each emitting a small number `K` of
   *tokens* per dataset rather than one vector, so nothing is pooled away
   before datasets meet. (b) *A second stage over the union of tokens*
   with dataset-identity embeddings — a set transformer, or a
   Perceiver-style latent bottleneck (`M` latent tokens cross-attending
   to all tokens, cost `O(N·M)` rather than `O(N²)`) when the token count
   is still large. Cross-dataset synthesis — a line strength read against
   a photometric colour, an image's morphology against a spectrum's
   slope — happens in stage (b), which sees every dataset's tokens; it is
   not lost by splitting stage (a) per dataset, because stage (a) does
   not pool to a scalar. Both stages are wrappers over `unpack`; the only
   layout addition they need is a per-dataset row cap, which is a field,
   not a version bump. Until then, `row_cap` and chunking are the
   controls, and an `Image` observation under a flat set embedding is
   *allowed* but expensive — the refusal is on rows, not on kind.
2. **Coordinate features are per axis.** Multi-axis kinds get independent
   Fourier features per axis, not mixed 2D features; sufficient for a
   first encoder, and a wrapper may add its own.
3. **No context yet** (group 10 is width 0). W5.10 amortises over the
   observation context through the σ columns instead; see §3 group 10.
4. **Standardisation from the observation only.** A layout built for one
   observation and reused for a very differently-scaled one is honest but
   not optimal; the remedy is a layout built from the new observation,
   which by construction has a different hash and therefore a retrained
   network — which is the correct outcome.
5. **Extra coordinates are not features** in v1 (they may be non-numeric
   and are rarely what a network should condition on); a numeric extra
   coordinate a user wants seen is a later layout option.
6. **The `coordinate` group aligns axes by position, not by name.** *(Added
   W4.8, from W4.3's carried finding.)* §3 group 2 packs "the sample's
   coordinate on each of its container's axes, in the container kind's axis
   order" — column `k` is whichever axis a kind declares in position `k`,
   with no check that two datasets' column `k` mean the same physical
   quantity. That is harmless while every kind in a collection has axes in
   a broadly comparable order (a `Spectrum`'s one spectral axis, an
   `Image`'s `x`/`y`), but W4.3's first complex, multi-axis customer —
   `VisibilitySet`'s `(u, v, spectral_axis)` and `ClosurePhases`' five —
   makes the assumption visible rather than incidental: nothing stops a
   layout mixing a kind whose column 0 is a spatial coordinate with one
   whose column 0 is `u`, and a network sees two unrelated quantities in
   the same input slot. Whether the packing should carry axis *identity*
   (a physical-type or axis-name feature per column, not just position) is
   a genuine design question this draft leaves open rather than answers —
   a Phase 5 question, not a defect in what shipped. ***Closed by W5.11***,
   which answers it with a physical-type code per column: §3 group 4 is the
   packing's own statement of what each coordinate column holds, the code
   table there is contract, and `DatasetLayout.axis_codes` puts it in the
   layout and therefore in the hash. What remains open is deliberately
   narrower — the codes name a *physical type*, not a role, so two `length`
   axes of one collection (a wavelength and a baseline) still share a code,
   and a kind that wants them distinguished needs a new code rather than a
   new mechanism.

## 10. Handoffs

- **W3.3 landed** `ampere/core/encoding.py` (`EncodingLayout`, `encode`,
  `encode_observations`, `decode`, `unpack`, `EncodingError`), the two
  wrappers and `layout=` in `ampere/inference/_sbi.py`, and bound this
  text; the decision-log row is amended with what changed. `_summary_of`
  survives as a thin adapter over `encode`, so the `"flat"` packing has one
  implementation rather than two. `decode` inverts the packing through the
  `value_asinh` column rather than the whitened one, because that column is
  defined whether or not a dataset has uncertainties; it reconstructs
  coordinates, values, uncertainties and the mask, and
  `Decoded.as_container(template)` borrows the axes, unit, extra
  coordinates and metadata the encoding does not carry.
- **W3.5** hashes the layout into the artefact key.
- **W3.6** encodes the calibration batches with the run's own layout, so
  a coverage test is a test of the network as trained.
- **W5.11 landed** the axis identity: `AXIS_TYPE_CODES` and
  `axis_type_code` in `ampere/core/encoding.py`, `DatasetLayout.axis_codes`
  in the layout and its hash, the `axis_identity` column group (§3 group
  4), `ENCODING_VERSION` 2, and the two by-name refusals of a pre-W5.11
  record (`EncodingLayout.from_dict`,
  `ampere.results.training.append_training_set`), which share one message
  through `axis_identity_complaint`.
- **The context item** (reserved) fills group 10 and adds the context prior
  to `simulate_many`. ***W5.10*** added the prior without filling the
  group; see §3 group 10.
- **Later encoders** (horizon notes, follow-up §4/5) are wrappers over
  `unpack`; none needs this document to change.
