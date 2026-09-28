<!-- Research memo for Phase 6 W6.13 (4) and (5), D12. Written by an Opus research agent on the orchestrator's brief 2026-09-28 and reviewed by the orchestrator (Fable), who spot-checked the load-bearing claims against PyPI's JSON API and the GitHub API the same day: astrostarfish 0.4.2 (2022-06-07, requires_python >=3.7,<3.10); hyperion 0.9.11 (2024-04-22, >=3.9); miepython 3.3.0 (2026-07-28, >=3.10); Starfish PR #161 merged 2026-07-22 by pscicluna. The conda-forge Python-version claims for hyperion are the agent's and were not re-checked. The recommendations are put to Peter; the twins follow what he accepts. -->

# External codes for the PHOENIX-star and carbon-star v2 twins: are there better-maintained alternatives?

Research memo, 2026-09-28. All dates and versions below were read from PyPI's JSON API, the GitHub API, conda-forge (anaconda.org API) or Zenodo on 2026-09-28, unless marked *unverified*.

## §1 The question

Two legacy examples depend on external codes that are brittle on Python 3.12–3.14:

- `examples/examples_paper/phoenixstar.py` uses **Starfish** (`download_PHOENIX_models`, `PHOENIXGridInterfaceNoAlpha`, `HDF5Creator`, `Emulator.from_grid/train/load_flux`) to fetch PHOENIX-ACES spectra for Teff 5000–8000 K, log g 4–5, [Fe/H] 0–0.5, train a PCA+GP emulator and evaluate spectra at arbitrary parameters. The model then rescales to a luminosity, applies CCM89 extinction (`extinction`) and produces Gaia/2MASS/WISE photometry (`pyphot`) and a Gaia-RVS-like spectrum (R≈11 000, 0.842–0.872 µm).
- `examples/cstar_model_test_sbi_v2.py` uses **Hyperion** through `ampere.models.Hyperion.HyperionCStarRTModel`: `AnalyticalYSOModel` with a tabulated stellar spectrum, one power-law envelope (ρ∝r⁻²; mass, r_in, r_out, r_0 free), a 1-D spherical-polar grid (`set_spherical_polar_grid_auto(251,1,1)`), raytracing, modified random walk, a peeled SED, MPI with `nproc=70`. Dust is built on every call by writing a parameter file and shelling out to the **`bhmie` executable** (a separate Fortran program, `hyperion-rt/bhmie`), then `BHDust` + `set_lte_emissivities`. The two free abundances (amorphous carbon, SiC) change the dust, so bhmie runs per simulation.

Are there better-maintained alternatives the v2 twins could switch to?

## §2 Starfish

### Incumbent state

- **PyPI trap**: `pip install Starfish` installs an *unrelated* package (spacetx image-based transcriptomics; 0.4.0, 2026-01-27; deps `slicedimage`, `regional`…). The astronomy package is **`astrostarfish`** on PyPI: latest **0.4.2, 2022-06-07**, `requires_python >=3.7,<3.10`, `numpy<2`, `astropy<6`, `nptyping==1.*`. **It cannot be pip-installed on 3.12–3.14 from PyPI.**
- **GitHub master is fixed**: [Starfish-develop/Starfish PR #161](https://github.com/Starfish-develop/Starfish/pull/161) ("Support Python 3.11–3.14 and NumPy 2", authored by pscicluna, merged 2026-07-22) moved to PEP 621, `requires-python >=3.11`, classifiers 3.11–3.14, floors instead of caps, dropped `nptyping`. CI on master green on 2026-07-22. **No release has been cut since**, so the route is `pip install git+https://github.com/Starfish-develop/Starfish@<sha>`. Before July 2026 the last code commit was 2023-03-17; 30 open issues. Licence BSD-3 (repo) / "BSD-4-Clause" (pyproject — inconsistent).
- **Data source**: the Göttingen PHOENIX server is up (a GET of `lte05000-4.50-0.0…HiRes.fits` returned 200, 6.3 MB, today; HEAD requests return 500, so naive link checkers will report it down). The example's range is roughly 26 Teff × 3 log g × 2 [Fe/H] ≈ 156 files ≈ 1 GB (*estimate from grid spacing, not measured*).
- **Emulator file size**: a Starfish emulator stores eigenspectra over the full processed wavelength grid; for PHOENIX HiRes that is likely several MB or more unless the wavelength range is truncated (*not measured*). Fitting under the 1 MB ruling needs a truncated/binned grid anyway.

### Candidates

| Name | Covers the need? | Last release | Python | Licence | Weight | Verdict |
|---|---|---|---|---|---|---|
| **Starfish** (git master) | Yes, fully (download, grid, PCA+GP emulator) | PyPI 0.4.2, 2022-06-07 (py<3.10); master 2026-07-22 supports 3.11–3.14 | 3.11–3.14 from git only | BSD | Pure Python; ~1 GB PHOENIX download | Viable as a one-off training tool, pinned to a commit |
| [gollum](https://github.com/BrownDwarf/gollum) | PHOENIX reader/dashboard; docs: "we do not *yet* support interpolation"; needs local files | 0.4.3, 2025-08-17 | `>=3.8`; CI modernised 2026-02 | MIT | Pure Python, pulls bokeh | No (no interpolation) |
| [pystellibs](https://github.com/mfouesneau/pystellibs) | Interpolates bundled libraries (BaSeL, Kurucz 2004, ATLAS9-Munari hi-res, BT-Settl low-res, TLUSTY, Rauch). No PHOENIX-ACES; Munari hi-res is optical only, so no single library covers RVS at R≈11 000 *and* WISE | GitHub v1.0.0, 2025-10-01; **not on PyPI** | `>=3.9` | MIT | ~235 MB of library FITS in the repo | No for a PHOENIX twin |
| [expecto](https://github.com/bmorris3/expecto) | Fetches one PHOENIX spectrum at a grid point as a `specutils` spectrum; no interpolation | 0.1.4, 2025-02-21 | `>=3.5` | BSD-3 | Pure Python | Useful as the *download* step only |
| [speclib](https://github.com/brackham/speclib) (Rackham) | PHOENIX and others, `Spectrum.from_grid(..., interpolate=True)`, synthetic photometry | GitHub v0.1.0b12, 2026-08-21; very active (2026-09-09). PyPI `speclib` is an unrelated GPL package | `>=3.11,<3.14` (**no 3.14**); `specutils<2` | MIT | Pure Python; downloads grids | Closest drop-in, but beta, git-only, excludes 3.14, not differentiable |
| [synphot](https://github.com/spacetelescope/synphot_refactor)/[stsynphot](https://github.com/spacetelescope/stsynphot_refactor) | `grid_to_spec('phoenix', T, [M/H], logg)` interpolates STScI's PHOENIX (Allard) catalogue at **~2 Å resolution** (R≈4000 at RVS) | synphot 1.7.0, 2026-03-07; stsynphot 1.5.1, 2026-03-09 | `>=3.10` | BSD-3 | Pure Python; TRDS data tarballs | Fine for photometry, too coarse for the R≈11 000 spectrum |
| pysynphot | Predecessor of synphot | 2.0.0, 2021-09-07 | `>=3.6` | BSD | Compiled | No (superseded) |
| [blase](https://github.com/gully/blase) | PyTorch semi-empirical clone of individual PHOENIX spectra; not a grid emulator | GitHub v0.3, 2022-06-23; last push 2025-05 (paper revisions). PyPI `blase` is an unrelated blazar package | — | MIT | torch | No |
| The Payne / MINESweeper | Payne: neural emulator trained on APOGEE H-band ([repo](https://github.com/tingyuansen/The_Payne), no licence, no PyPI, last push 2025-10). MINESweeper: repo not found; PyPI `minesweeper` is an unrelated 2016 package | — | — | — | — | No; the *idea* (small NN emulator) is what the own-emulator route implements |
| [smart](https://github.com/chihchunhsu/smart) | Instrument-specific PHOENIX-ACES/BT-Settl interpolation for RV work | GitHub only, push 2026-05-28 | — | MIT | Precomputed grids | No (instrument-specific) |
| **ampere's own emulator** (horizon (c)) | PCA+GP (or PCA+MLP) trained in the example on a downloaded PHOENIX subset, restricted to what the data need | n/a | whatever ampere supports | ampere's | Training needs the ~1 GB download once; the shipped file does not | **Recommended** |

### Recommendation (confidence: medium-high)

**Switch the twin to ampere's own emulator, shipped pre-trained.** Rationale: (i) no maintained package offers "PHOENIX-ACES, interpolated, R≈11 000 plus 0.3–5 µm, on 3.12–3.14" — speclib is closest but beta, git-only and excludes 3.14; (ii) the pre-trained file under 1 MB is only reachable by restricting the representation (e.g. a binned SED of a few hundred points for photometry plus the RVS window, ~400 pixels, a handful of PCA components: tens of kB), which Starfish's emulator format does not do by default; (iii) an emulator written against the core contracts is differentiable on the torch and jax backends, which Starfish's (numpy/sklearn) is not. Training is a separate script: download with a few lines of `astropy.utils.data.download_file` (or `expecto`), PCA with numpy, GP or small MLP per weight. Starfish (git master, pinned) remains available as a cross-check and is what the legacy example keeps.

## §3 Hyperion

### Incumbent state

- **PyPI**: `hyperion` **0.9.11, 2024-04-22**, `requires_python >=3.9`, abi3 wheels (`cp39-abi3`, linux x86_64 and macOS x86_64), so `pip install hyperion` succeeds on 3.13/3.14 — but the wheel contains the **Python front end only**. The Fortran binaries (`hyperion_sph` etc.) need a Fortran compiler, HDF5 and MPI; the classic failure is issue #224 ("`hyperion_car: not found`").
- **conda-forge** ships both: `hyperion` 0.9.11 (builds py39–py313, **no py314**, linux-64/osx-64/osx-arm64, last upload 2025-05-09) and **`hyperion-fortran`** 0.9.11 (with `mpich`, `hdf5`; last upload 2026-03-25; no Windows). So in pixi's Python 3.13 environments Hyperion is installable today without a compiler; for 3.14, pip's abi3 wheel plus conda-forge `hyperion-fortran` should work (*untested*).
- **Maintenance**: very active — ~40 commits in July–August 2026 by T. Robitaille (MPI, spherical-grid, peel-off and `PowerLawEnvelope` mass/ρ₀ fixes; last push 2026-08-04), CI matrix 3.11–3.13. No release since 2024, so the fixes are unreleased. BSD-2. 64 open issues.
- **The real brittleness is `bhmie`**: [hyperion-rt/bhmie](https://github.com/hyperion-rt/bhmie) (last push 2022-11-01) is a Fortran program on neither PyPI nor conda-forge; it must be compiled by hand and on `PATH`. Hyperion itself does not need it: `IsotropicDust(nu, albedo, chi)` and `HenyeyGreensteinDust(nu, albedo, chi, g, p_lin_max)` accept arrays (verified in `hyperion/dust/dust_type.py`), so opacities can be computed in Python with [miepython](https://pypi.org/project/miepython/) (3.3.0, 2026-07-28, `>=3.10`, MIT, pure Python + numba) from the example's `.optc` files.
- **Cost**: Monte Carlo noise and minutes per model at `nproc=70`; for SBI the cost is the simulation budget, not the dependency.

### Candidates

| Name | Covers the need? | Last release / activity | Python | Licence | Weight | Verdict |
|---|---|---|---|---|---|---|
| **Hyperion** | Yes (3-D MC, used in 1-D) | 0.9.11, 2024-04-22; active 2026-08 | front end abi3 ≥3.9; conda-forge ≤3.13 | BSD-2 | Fortran + MPI + HDF5 (conda-forge binaries) | Keep, with bhmie replaced |
| [DUSTY](https://github.com/ivezic/dusty) v4 (legacy `ampere/models/Dusty.py` wraps it) | Yes, and better suited: exact 1-D spherical solution, no MC noise, seconds per model, user n,k files, r⁻² or wind density | No releases; last commit 2021-01-12 | n/a (Fortran exe) | BSD-3 | gfortran build by hand; not packaged anywhere | Physics fits best, packaging worst. Wrappers: [pyDusty](https://github.com/mgomezAstro/pyDusty) (push 2026-08-24, no licence), [DustyPY](https://github.com/gtomass/DustyPY) (2026-02, MIT); both git-only |
| RADMC-3D ([radmc3d-2.0](https://github.com/dullemond/radmc3d-2.0)) | Yes (1-D spherical supported, user opacities) | No releases; push 2026-09-03 | radmc3dPy not on PyPI | no licence file detected by GitHub | Fortran build by hand | No gain over Hyperion |
| [MCFOST](https://github.com/cpinte/mcfost) | Disc-oriented MC code; spherical envelopes possible | v4.1.14, 2026-08-23 | `pymcfost` 0.1.1 (2023) is post-processing only | custom (NOASSERTION) | Fortran binaries | No |
| [SKIRT 9](https://github.com/SKIRT/SKIRT9) | Yes (general MC, 1-D geometries, custom dust) | active, push 2026-09-21; no GitHub releases | PTS toolkit not on PyPI | AGPL-3.0 | C++ build with CMake | No (heavier, AGPL) |
| [sedfitter](https://github.com/astrofrog/sedfitter) | Fits pre-computed Robitaille YSO grids; no AGB/carbon-star grid | 1.5, 2026-08-02 | `>=3.11` | BSD-2 | Pure Python | No for this physics; its grid-fitting pattern is relevant |
| [2-DUST](https://github.com/sundarjhu/2-DUST) | Axisymmetric dust RT (the code behind GRAMS) | push 2021-05-20 | n/a | GPL-3.0 | Fortran | No |
| MoDust / "More of DUSTY" (Groenewegen) | DUSTY-based fitting driver | *unverified — no public repository found* | — | — | Fortran | No |
| Pyradex | **Not applicable**: non-LTE molecular *line* RT (RADEX), no dust continuum; PyPI 0.2.2, last upload 2016 | — | — | — | — | No |
| Pure-Python 1-D dust RT | None maintained found in searches | — | — | — | — | — |
| **GRAMS carbon grid** via [DESK](https://github.com/s-goldman/Dusty-Evolved-Star-Kit) ([Zenodo 10.5281/zenodo.14448621](https://doi.org/10.5281/zenodo.14448621), CC-BY-4.0) | 12 244 2-Dust models of LMC carbon stars (Srinivasan et al. 2011): amorphous carbon + **fixed** 10 % SiC, KMH sizes, COMARCS photospheres, r_in ∈ {1.5,3,4.5,7,12} R★, 26 τ(11.3 µm) from 0.001–4, 0.2–200 µm. `grams-carbon_models.hdf5` 25.5 MB + outputs 0.7 MB | DESK 1.9.1, 2025-11-21 (only needed for its reader, not required) | h5py + numpy | CC-BY-4.0 (grid), BSD (DESK) | 26 MB download, no compiled code | Strong lightweight alternative; loses free composition |
| Nanni et al. (2019) dust-growth carbon grids (same Zenodo record) | Self-consistent wind models, LMC/SMC | as above | h5py | CC-BY-4.0 | 700–750 MB each | Too heavy for an example |

### Recommendation (confidence: medium)

**Hybrid, keeping Hyperion as the physics**: no candidate is both better maintained and packaged for 3.12–3.14 while supporting free dust composition. Hyperion is actively developed and conda-forge supplies compiled binaries for Python ≤3.13, which is ampere's pixi default; the unpackaged piece is `bhmie`, which the twin can drop by computing opacities with `miepython` and passing arrays to `IsotropicDust`/`HenyeyGreensteinDust`. Put Hyperion in an optional pixi feature (conda-forge `hyperion` + `hyperion-fortran`), and let the SBI engine's cached training set make the expensive simulation a one-off. **If Peter prefers zero compiled dependencies**, the twin should instead interpolate the GRAMS carbon grid. The data (an LMC OGLE carbon AGB star with an IRS spectrum) is exactly GRAMS's target population, but the SiC fraction is fixed at 10 % and r_out/r_0/stellar mass are not free, so the parameter set changes and the SiC-abundance inference the legacy example demonstrates is lost. DUSTY is the physically cleanest 1-D code, but it has no releases since 2021 and no packaging, so it is not "better maintained".

## §4 What each v2 twin would look like

**PHOENIX twin.** A training script (run once, not in CI) downloads the PHOENIX-ACES subset from Göttingen, bins it to a coarse 0.3–5 µm SED grid plus the RVS window at R≈11 000, fits a PCA with a GP (or small MLP) per component weight over (Teff, log g, [Fe/H]), and writes a sub-MB `.npz`. The example loads that file into an emulator `Model` written against the core contracts (numpy reference, torch and jax twins, so NUTS works), applies luminosity scaling and CCM89 extinction (`extinction` 0.4.9, 2026-08-09, wheels through 3.15), and fits photometry plus the RVS spectrum with `SBIEngine` and a gradient-based sampler for comparison. Starfish is not imported; a comment points to the legacy example for the Starfish route.

**Carbon-star twin (recommended hybrid).** A v2 `Model` wraps Hyperion's `AnalyticalYSOModel` with the same shell parameters, building dust in Python (miepython over the `.optc` n,k tables and power-law size distributions, mixed by the abundance parameters) instead of calling `bhmie`. It runs in an optional pixi environment with conda-forge Hyperion and fits with `SBIEngine`, caching the simulation bank so reruns and calibration reuse it. The GRAMS-grid alternative would be a numpy/jax interpolator over (Teff, r_in, τ₁₁.₃) with luminosity scaling, fetching the 26 MB file from Zenodo at first use; it would be seconds end to end and differentiable, but it would not fit abundances.

## §5 Open questions for Peter

1. Carbon star: keep Hyperion physics (optional pixi env, bhmie replaced by miepython), or switch to the GRAMS grid and accept a fixed 10 % SiC and a changed parameter set? Or offer both, with GRAMS as the fast default?
2. Should Starfish cut a release (0.5.0) from master, given that PR #161 is yours and PyPI `astrostarfish` still requires Python <3.10? That would fix the legacy example's install story independently of the twin.
3. PHOENIX twin emulator: is PCA+GP acceptable, or should it be PCA+MLP (cheaper, easier to make differentiable on every backend)? And is fetching the ~1 GB grid for retraining acceptable as a documented manual step?
4. Is a 26 MB runtime download from Zenodo (CC-BY-4.0, needs citation of Srinivasan et al. 2011 and DESK) acceptable for an example, with caching outside the repo?
5. Hyperion's many 2026 fixes (including `PowerLawEnvelope` mass/ρ₀) are unreleased; conda-forge still ships 0.9.11 from 2024. Pin 0.9.11, or ask upstream about a 0.9.12/1.0 release?

## Side findings (out of scope, not fixed)

- `ampere/models/Dusty.py`: `dusty_run_dusty()` is defined without `self` and calls `subprocess.call("dusty dusty_model.inp")` without `shell=True` (a single string argument without the shell fails on POSIX). The legacy DUSTY wrapper probably does not run as written.
- `HyperionCStarRTModel.__init__` references an undefined name `temperature` in its `elif` branch.
- `phoenixstar.py` uses `np.trapz`, deprecated since NumPy 2.0.
