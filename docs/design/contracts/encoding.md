# Ampere v2 — Encoding Contract (draft, W3.3)

Status: **draft, 2026-09-09** (Fable, at Peter's request, while W3.1
slice 2's gates ran). A **post-freeze §4 addition** — `DEVELOPMENT_PLAN.md`
§5's Phase 3 second bullet ("a canonical coordinate–value–mask tensor
encoding of containers for embedding networks") made concrete. Recorded in
the plan's decision log the same day; **W3.3 lands the code and binds this
text**, amending it in place where the code proves it wrong, and the
document is frozen at the W3.3 merge. Until then the binding parts are
§3–§7; §8–§10 are guidance.

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
- **Encoded**: the `(rows, columns)` float64 array for one observation, or
  `(count, rows, columns)` for a batch, plus its layout.
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
   sees *where* a sample is, at more than one scale.
4. **`value`** [1 or 2] — the whitened value `y/σ` (real and imaginary
   columns when the container is complex, else one column). Where the
   dataset has no uncertainties, the asinh-scaled value (§4) stands in and
   the `has_sigma` feature (group 8) says so.
5. **`value_asinh`** [1 or 2] — `asinh(y / s)` with `s` the dataset's
   value scale (§4): the magnitude of the measurement on a scale that is
   linear near zero and logarithmic far from it, so a network can tell a
   bright line from a faint one without the raw unit.
6. **`log_sigma`** [1] — `log(σ / s)`; 0 where there is no σ.
7. **`mask`** [1] — 1 for a row that counts, 0 for a row the dataset's
   `effective_mask` excludes **and** for every padded row. There is
   exactly one mask convention in the contract; §7 says how wrappers
   convert it.
8. **`set_features`** [3] — per-dataset quantities broadcast onto every
   row of that dataset: `log(N_valid)` (valid rows in this dataset),
   `has_sigma` (1/0), and `is_complex` (1/0). They exist because a masked
   mean (§7) deliberately removes the row count from the pooled
   embedding, and the count is information.
9. **`context`** [0 today] — the reserved observation-context group.
   Width 0 until the context item lands; its presence at position 9 is
   what makes adding it a layout-hash change and nothing else.

The values in groups 4–6 are computed from the container the row belongs
to, with masked samples' values **left as the container holds them** and
then zeroed by the mask at unpack time (§7) — the contract does not invent
values for masked samples, mirroring `inference.md` §13's rule for
`simulate`.

Rows are ordered dataset by dataset in `datasets` order, and within a
dataset in the container's storage order (flattened C-order for multi-axis
kinds). Row order carries **no meaning** — a set embedding must be
invariant to it, and the conformance row for the wrapper asserts it — but
it is fixed so that an encoding is reproducible bit for bit.

**Padding.** The layout fixes a row cap `R` (§5). An observation with
fewer rows is padded to `R` with all-zero rows and `mask = 0`; one with
more is **refused by name** (the remedy is a larger cap in a new layout,
or masking). Padding is what lets a batch be one rectangular tensor.

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

A frozen, plain-data record (`Layout`) with, at least:

- `kind` (`"set"` | `"flat"`), `version` (`ENCODING_VERSION`, 1);
- `datasets`: for each label in order — `kind` (container class name),
  `axes` (names, count), `rows` (valid sample count of the observation),
  `is_complex`, `has_sigma`, the coordinate ranges and the value scale;
- `coordinates` (`A`, the maximum axis count), `fourier_bands` (`B`),
  `row_cap` (`R`; default the total observed row count, so the default
  layout fits the observation exactly and a differently-sampled
  observation needs a wider cap chosen on purpose);
- `columns`: the resolved column groups with their offsets and widths
  (derivable, stored for interpretability);
- `hash`: `hash_of(to_dict())` (`ampere.results.provenance.hash_of`),
  over everything above.

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

## 7. `unpack`, wrappers, and the one mask convention

`unpack(x, layout) -> Unpacked`: named views of the groups (`dataset`,
`coordinate`, `coordinate_features`, `value`, `value_asinh`, `log_sigma`,
`mask`, `set_features`, `context`), `features` (every non-mask group
concatenated, the tensor a per-row network consumes) and `valid` (the
boolean mask). It is the **only** way a network is meant to read the
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

## 8. Stability

`ENCODING_VERSION` is bumped when a column group is added, removed,
reordered or redefined; the reserved `context` group exists so that the
first real context is a width change under version 1, not a version
bump. The packing is independent of `CONTAINER_SCHEMA_VERSION`: it reads
containers through their public attributes, and a container schema change
that keeps `axes`/`values`/`uncertainty`/`mask` leaves every encoding and
every hash unchanged.

## 9. Deliberate limitations of this draft

1. **One tensor for the whole collection.** A long spectrum beside a few
   photometric points shares one row set; attention is quadratic in rows.
   The extension is per-dataset sub-encoders over per-dataset row caps,
   concatenated (a "multi-set" wrapper) — a later item, not a packing
   change, since each sub-tensor is this packing restricted to one
   dataset.
2. **Coordinate features are per axis.** Multi-axis kinds get independent
   Fourier features per axis, not mixed 2D features; sufficient for a
   first encoder, and a wrapper may add its own.
3. **No context yet** (group 9 is width 0).
4. **Standardisation from the observation only.** A layout built for one
   observation and reused for a very differently-scaled one is honest but
   not optimal; the remedy is a layout built from the new observation,
   which by construction has a different hash and therefore a retrained
   network — which is the correct outcome.
5. **Extra coordinates are not features** in v1 (they may be non-numeric
   and are rarely what a network should condition on); a numeric extra
   coordinate a user wants seen is a later layout option.

## 10. Handoffs

- **W3.3** lands `ampere/core/encoding.py` (`Layout`, `encode`, `decode`,
  `unpack`), the wrappers in `ampere/inference/_sbi.py`, `layout=` on
  `SBIEngine`, and binds this text; the decision-log row is amended with
  what changed.
- **W3.5** hashes the layout into the artefact key.
- **W3.6** encodes the calibration batches with the run's own layout, so
  a coverage test is a test of the network as trained.
- **The context item** (reserved) fills group 9 and adds the context prior
  to `simulate_many`.
- **Later encoders** (horizon notes, follow-up §4/5) are wrappers over
  `unpack`; none needs this document to change.
