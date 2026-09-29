# Changelog

## Unreleased

### Fixed
- **`Dist_Stats.JSD` computed the wrong quantity.** Since the first release it evaluated
  `-1/2 sum[P log(2Q/(P+Q)) + Q log(2P/(P+Q))]` -- P and Q swapped inside the logarithms --
  and returned its square root. That equals the square root of half the symmetric (Jeffreys)
  KL divergence *minus* the Jensen-Shannon divergence: unbounded and not a Jensen-Shannon
  quantity (a 1415-atom melting frame gave 1.04 against a true JS distance of 0.09). It now
  returns the Jensen-Shannon **distance** in base 2, bounded in [0, 1], matching
  `scipy.spatial.distance.jensenshannon(P, Q, base=2)`. **JSD values from earlier versions
  are not comparable with new ones**; they can be recomputed from the stored per-frame PDFs.
- **`Dist_Stats.JSD` and `Dist_Stats.Kullback` normalise their inputs.** PDFs are stored as
  densities (they integrate to 1 over r, so their values sum to 1/dr, ~33 on the default
  grid); divergences are defined between probability vectors. Both now normalise to unit sum,
  which also makes them independent of the grid spacing. KL values change by that factor.

## 1.3.0.dev0 — command line, sparse adjacency, parallel frames (2026-09-16)

### Added
- **`sapphire` command line interface** (`Sapphire.cli`, installed as a console script).
  `sapphire run` analyses a trajectory, `sapphire expand` converts sparse adjacency matrices
  back to dense text, `sapphire quantities` lists what can be asked for. Built on the existing
  `Sapphire.api.Config`, so a run is reproducible from one TOML file.
  - Prerequisites are resolved automatically: asking for `cna_sigs` alone used to log a
    `KeyError` for `FullCut`, then another for `Adj`, write nothing and exit **0**. It now
    pulls in `pdf` and `adj`, and a run that produced nothing exits non-zero.
  - `--no-quote` / `SAPPHIRE_NO_QUOTE` skips the log header's quote lookup, which reaches the
    network from inside `Process.__init__`.
  - `-j/--jobs N` analyses frames across N processes.
  - `--adj-format {npz,text}`, `--dry-run`, `--strict`, `--frames START:END:STEP`.
- **Parallel frame analysis** (`Sapphire.parallel`). Workers analyse contiguous chunks in
  private directories; results merge in frame order, and `analyse()` then runs once over the
  merged output exactly as after a serial run. Output is byte-identical to serial, verified
  for mono- and bimetallic runs at 2, 3, 4 and 8 jobs.

### Changed
- **Adjacency is written sparse (`scipy.sparse` npz) by default**, ~165x smaller than dense
  text — roughly 160 GB to 1 GB over a 20 000-frame run at N = 2000. `Reader` loads either
  encoding transparently and still returns dense arrays; `--adj-format text` and
  `sapphire expand` recover the old format byte for byte. See `docs/FILE_CONTRACT.md`.
- **CNA signature tables are rectangular.** The masterkey grows as new signatures appear, so
  rows written early were narrower than rows written late and `np.loadtxt` could not read the
  file. It is now laid out against the final masterkey with zero padding; no counts change.
- `Process` reads only the frames it will analyse rather than the whole trajectory.

### Fixed
- `Utilities/ExtendXYZ.py` paired `Traj[i]` with quantity row `i`, so extended-xyz output used
  the **wrong atoms** for every row whenever `Step > 1`. Fixed by the sliced read above.
- `Process.analyse` sized the concertedness series for one more entry than the loop ever
  assigned, so element 0 was written to file as a hard `0.0` indistinguishable from a
  measurement. A three-frame run produced *only* that fabricated value. Concertedness spans a
  lag of two and is defined from the fourth analysed frame; the series is now `n-3` long.
  `Graphing.Plot_Funcs.h_c` plots collectivity alone when concertedness is undefined.
- `Utilities/Initial.py` — a `wikiquote` failure propagated out of `Process.__init__`, so a
  compute node with no outbound route could not construct a run. Now bounded by a timeout and
  never raised. The version in the log header was hardcoded `1.0.0`; it reads the package.

### Performance
Frame cost fell from 0.328 s to 0.055 s serially and 0.015 s across 8 processes (N = 1415),
with every cutoff and all derived quantities bit-identical throughout.
- `CNA/FrameSignature` densified the sparse adjacency and scanned two full N-length columns
  per bonded pair, making the routine **quadratic in N** for physics that is linear in it
  (per-atom cost grew 49.7 → 270.6 µs between N = 147 and N = 6525). It now keeps an adjacency
  list; cost is flat in N.
- The longest-chain search behind a signature's `t` is cached on the bond graph. A 923-atom
  icosahedron presents 11 distinct graphs across 9804 calls — a 99.9% hit rate.
- `Post_Process/Adjacent` called `np.savetxt` once per atom (1415 calls per frame); matrices
  are now written in one pass.
- `Post_Process/Kernels._gaussian_sum` writes the Gaussian out instead of calling
  `scipy.stats.norm.pdf` (~2.3x). Values agree to ~8e-15 relative; every cutoff is unchanged
  bit for bit.
- `Post_Process/AtomicEnvironment` densified the adjacency for `ele_nn` and `lae`; both
  operations work on sparse directly (3.2 GB per frame at N = 20 000, for nothing).

### Packaging
- `packages.find` defaulted to namespace discovery, so every directory under `main/Sapphire`
  became a package: 1.1.0 and 1.2.0 shipped ~30 bogus importables such as
  `Sapphire.Tutorials.08_Ensemble_Averaging.ensemble.seed0` and the unimportable
  `Sapphire.Sapphire-logos`. `namespaces = false` plus `exclude-package-data` for tutorial run
  outputs; the wheel drops from 12.5 MB to 8.1 MB and declares 10 packages.

## 1.2.0 — released 2026-09-15
Zenodo DOI in citation metadata and README; documentation moved out of internal planning
notes; local build artefacts untracked. Release tagging is now guarded: the workflow refuses
to build unless the tag, `pyproject.toml`, `CITATION.cff` and `Sapphire.__version__` agree and
the version is not a development one, publishing is idempotent and verified against PyPI, and
the GitHub release waits on a successful upload.

## 1.1.0 — released 2026-08-31
Everything below (the 2026 restoration, Phases 1–8) constitutes release 1.1.0.

## 1.1.0.dev0 — restoration (2026-08-28)

### Repository
- Removed ~400 tracked build artefacts (`build/`, `dist/`, `.eggs/`, `*.egg-info`, `__pycache__`, `.ipynb_checkpoints`), a copied `.git/` internals dump (`HEAD`, `config`, `description`, `hooks/`, `info/`), duplicate `README`/`setup.py` under `main/`, and an Excel lock file.
- Removed the never-initialised `Raffy` submodule.
- New `legacy/` for quarantined code, `examples/` for driver templates, `tests/`, `docs/`.

### Packaging
- `setup.py` → `pyproject.toml`; `requires-python >= 3.10`; version `1.1.0.dev0`.
- Core dependencies reduced to `numpy scipy ase pandas networkx`; everything else is an extra (`plot`, `changepoint`, `ml`, `cna`, `light`, `quote`, `notebooks`, `dev`, `all`). Rationale for the dropped packages is summarised in `docs/ML_POTENTIALS.md`.

### Importability (behaviour-preserving)
- `Light/Epsilon_DFT.py` — Ni Johnson & Christy table was pasted as two tab-separated columns and did not parse (SyntaxError). Rebuilt as three 150-point arrays; values unchanged.
- 24 Python-2-style implicit imports rewritten as `Sapphire.*` absolute imports.
- `CNA/Model.py` — `sklearn.externals.joblib` → `joblib`.
- Glob-generated `__all__` replaced with explicit lists in every `__init__.py`.
- `IO/input.py`, `Graphing/Plot_input.py` were scripts executing on import → `examples/`.
- `Emerald/`, `CNA/main.py`, `CNA/FrameSignature_Ovito.py`, `IO/ExtendXYZ.py`, `IO/OutWrite.py`, `IO/*.sh` → `legacy/`.

### Bug fixes (may change behaviour — each was a crash or a silent no-op before)
- `Utilities/ExtendXYZ.py` — a second `class Extend` (an Au/Pt one-off script) shadowed the real one, so `Process` could never write extended-xyz output. Script moved to `legacy/Utilities/ExtendXYZ_AuPt_script.py`. `Process.py` now imports `Sapphire.Utilities.ExtendXYZ`.
- `Utilities/Pattern_Clean.py` — `Process(Pattern_Input=None)` (the documented default) raised `TypeError`; defaults are now pre-filled. `FROM_MEMORY` cleaner wrote to `System` instead of `Pattern_Input`.
- `Post_Process/Stats.py` — Jensen–Shannon `calculate()` referenced bare `P`, `Q` (NameError). Now uses `self.P`, `self.Q`; loop vectorised, same formula.
- `Post_Process/Mass_Activity.py` — `mass_NP` was undefined; now `len(traj[j]) * mass_cu`. **Please sanity-check this physics choice.**
- `Graphing/Reader.py` — `Get_Heights_Ovito` used `Metadata` instead of its `CNAs` argument; `Quant.lower is 'masterkey'` (always False) → `Quant.lower() == 'masterkey'`.
- `Graphing/Plot_Funcs.py` — `autolabel` now takes `ax`; `ax2Ticks`/`tick_function` typos → `ax3Ticks`/`self.tick_function`; `is` on string literals → `==`.
- `Graphing/Plotter.py` — Lindemann/CoM helper block (no `self`, six undefined globals) → `legacy/Graphing/Lindemann_CoM_script.py`.
- `CNA/Model.py` — `np.maximum(np.max(X, axis=0))` (TypeError) → `np.max(X, axis=0)`; six bare names → `self.*`; ndarray writes wrapped in `str()`.
- 50 unused imports removed (ruff F401).

### Tests
- `tests/` with import sweep, `DistFuncs` numerics, and a `Process` smoke run on `Tutorials/CNA/Au561.xyz`.

### numpy 2 / end-to-end fixes found by the smoke test
- `np.trapz` (removed in numpy 2.0) → `np.trapezoid` with fallback, in `Kernels`, `Plot_Funcs`, `Read_Plot`. Before this, **the PDF failed, so no cutoff was found, so adjacency/CN/CNA all silently failed** for every user on numpy ≥ 2.
- `CNA/FrameSignature.py` — row-iteration over the scipy sparse adjacency matrix raised under numpy 2; the matrix is densified once and neighbour lookups are vectorised (same results, much faster).
- `Process.Initialising` — `All_Times`/`Band` were only set in the bimetallic branch, so **every monometallic run crashed**; missing `Quantities` groups (`Homo`/`Hetero`) now default to `{}`.
- `Process.__init__` — `self.metadata = {}` so `analyse()`/`write_meta()` no longer raise; note the in-memory metadata design was superseded by `Time_Dependent/` file output in 1.0 and those two methods need a Reader-based redesign (open item).

### Follow-up (2026-08-28, after first push)
- CI was failing at the `ruff` step (24 residual findings, exit code 1) before `pytest` ran. All cleared: unused result bindings in `Process.calculate`, unused `as e`, dead `freq`/`fig`/`ax`/`Plot`/`f` locals, and `for self.x in …` in `IO/Output.py`.
- `Mass_Activity` — nanoparticle mass is now Σ over atoms of the species' relative atomic mass (`ase Atoms.get_masses()`), converted to mg via `AMU_MG`. Replaces the single-species `mass_cu` constant.

## Phase 6 — file-backed metadata, strict mode, logging (2026-08-28)
- **`Sapphire.IO.Reader`** — reads a run directory (`Time_Dependent/`, `CNA/`, `Adjacency/`, `Exec/`) back into arrays keyed by the metadata names (`nn`, `pdf`, `cna_sigs`, `hocomAu`, `adj`, …). Tolerant of absent quantities. Any external code that writes `<frame> <values…>` under `Time_Dependent/` is readable the same way.
- `Process.load_metadata()` / `Process.reader()`; `analyse()`, `write_meta()` (now `Metadata.pkl`) and the extended-xyz writer run on the Reader output instead of an empty dict.
- `Process(strict=True)` re-raises instead of log-and-continue; 17 copy-pasted `except` blocks collapsed into `Process._report()`, which also emits a `logging` warning.
- `Sapphire.Utilities.log` — 31 `print()` calls now go through the `Sapphire` logger.
- `Stats.Dist_Stats.Kullback/JSD` returned the divergence *object*, never a value; now return the number. `JSD_Dist` no longer mutates its inputs.
- `Utilities/ExtendXYZ.py` — `Names.pop(name)` (TypeError) fixed; writes a valid multi-frame xyz (`N`, comment, atoms per frame) rather than a trailing-count layout.
- Tests: `test_reader.py` (parsing, shapes, adjacency ↔ NN consistency, analyse/JSD, write_meta, strict, extend_xyz), `test_graphing.py`.

## Phase 7.1 — bimetallic paths, Graphing on Reader, exemplar figures (2026-08-28)

### Bugs found by the first-ever bimetallic run (AuPt sample) and fixed
- `DistFuncs.CoM_Dist` ignored its `CoM` argument (`self.CoM = Positions`), so **every `CoMDist` ever written was a per-atom 3-vector, not a distance**; `get_CoM` returned a scalar (`np.average` without `axis=0`). Homo CoM distances now use the sub-species centre.
- `Process`: Homo- and Hetero-RDF calls passed positions positionally into `System`; the Homo-RDF block was indented inside the `except` of the Full RDF (so it only ran when Full RDF failed); the LAE/mixing block was a stub with wrong keyword names under a copy-pasted "Gyration" error message.
- `Post_Process/AtomicEnvironment.py` rewritten as a coherent module: `Mix` (mixing parameter + homo/hetero bond counts), `LAE` (hetero-neighbour histograms per species), `Ele_NN` (per-atom neighbour counts by species). Wired into `Process` and verified bond-for-bond against the full adjacency.
- `analyse()`: called `Adjacent.R/Collectivity/Concertedness`, which live in `Stats.Mobility` (AttributeError swallowed → collectivity always 0); looped `range(1, (End-Start)/Step)` instead of over loaded frames; treated `'collect': None` as "not requested"; matched `'pdf' in key` so divergences were computed on the `*space` grids too. Collectivity/concertedness and every statistic are now written to `Time_Dependent/` (`Collectivity`, `Concertedness`, `Stats/<Stat><quantity>`).
- `Stats.Mobility.R` accepted only sparse input and summed the *signed* neighbour difference (lose one + gain one = "unchanged"); now any change counts.
- `IO.Reader`: per-frame matrix sets (`Adjacency/File<n>`, `HeAdjFile<n>`, `HomoAdjAuFile<n>`) collapse to one 3-D array per key; `Time_Dependent/Stats/*` pass through; `masterkey` stays as strings.

### Graphing
- `Graphing.Read_Meta` rewritten as an adapter over `IO.Reader` producing the legacy layout `Plot_Funcs` expects (`(space, heights)` per frame, `R_Cut`, `Cut<X>`, normalised `cna_sigs` with right-padded ragged rows, KDE `CoMDist`/`MidCoMDist<X>` over `CoMSpace`, `h`/`c`, `SimTime`/`Temp` from `frame_dt`/`temperature`). Multi-run averaging = mean/std across `iter_dir`.
- `Plot_Funcs`: `inspect.getargspec` (removed in 3.11) → `getfullargspec`; `sys.exit` → exceptions; output paths via `os.path.join`; `agcn_heat` ticks matched to data; `plot_stats` time axis from `SimTime`; `cna_traj` tuple labels; `h_c` aligned to frame pairs; missing quantities skipped with a log line instead of crashing.
- `examples/make_figures.py` regenerates the 26-figure exemplar set into `assets/` (gitignored; provenance in `assets/LOG.md`).

### Tests
- Shared 3-frame AuPt fixture (`tests/conftest.py`), `tests/test_bimetallic.py` (shapes, bond bookkeeping, mixing parameter, LAE, collectivity, melting raises JSD), Graphing figure-rendering tests (12 figures).

## Phase 7.3–7.7 (2026-08-28)
- **`Post_Process/Mass_Activity.py`** rewritten to Rossi, Asara & Baletto, ChemPhysChem 2019 (Eqs. 5–7): volcano `A_l=exp(3.14α−23.40)`, `A_r=exp(−4.96α+42.18)`, sites with GCN > 6, `MA = j_flat/ρ_sites · Σ A / M_NP`; the 4.107 A/mg prefactor is reproduced from 2 mA cm⁻² and 1.503·10¹⁵ cm⁻² and the mass is species-weighted (alloys). The old file could not run (`beta` function used as a scalar). **Note:** the printed branch coefficients intersect at GCN 8.10, not the stated 8.33; we switch branches where they meet (continuous) — `apex=8.33` reproduces the literal reading.
- **`Post_Process/Morphology.py`** — Emerald's surface/core peeling, shell thickness, faceting ratio, surface area, volume, solid angle; vectorised; tests reproduce the Mackay shells [252, 162, 92, 42, 12, 1] of Au561.
- **`CNA/Classify.py`** — fingerprint → bulk-pattern features → SVC (classes fcc/Ih/Dh/amorphous), filename labelling, save/load; test trains on ASE Ih/Dh/Oct clusters and predicts held-out sizes.
- **`Utilities/errors.py`** — single `report()` honouring strict mode; `Process._report`, `Adjacent`, `Kernels`, `DistFuncs`, `FrameSignature` route through it (silent `except: pass/return None` no longer silent).
- `Adjacency_Matrix`: `Type` defaults to `'Full'` and no write without a `System` (standalone use, as in the tutorials, previously produced an empty list). `FrameSignature`: no write without `System`; `Exec/Masterkey` is overwritten per frame rather than appended (it accumulated duplicates).
- Docs: `docs/FILE_CONTRACT.md` (run-directory format for external producers), `examples/from_lammps.py`, `docs/ML_POTENTIALS.md` (brief + proposal), `legacy/README.md` marks Emerald and CNA/main superseded.

## Phase 7.2 — tutorials as a graded course (2026-08-29)
- Eight new executable notebooks under `Tutorials/0N_*/Tutorial.ipynb` (build & morphology → PDDF/cutoff → adjacency/aGCN/mass activity → CNA signatures/patterns/classifier → `Process` on a bimetallic melt → divergences/collectivity/change-points → shape → ensemble averaging). Each runs offline on bundled data; all outputs verified physically (Mackay shells, five Ih signatures with exact counts, JSD change-points at 5 and 10 ns). Previous notebooks (one real, five byte-identical clones, one stub, the Emerald script) moved to `legacy/Tutorials/`. `Data/Au561.xyz` is the bundled reference structure.
- `.github/workflows/notebooks.yml` executes every tutorial on push/PR and monthly.
- `Process(overwrite=True)` (default) clears previous result files from the Sapphire-owned output folders before a run — results were appended, so re-running in the same directory doubled every series.
- `FrameSignature`: `cna_sigs` without `cna_patterns` raised `AttributeError: Fingerprint` (silently, per frame) — attribute now always defined.
- `analyse()` no longer re-ingests its own earlier statistics (`JSDpdf`) as distributions.
- Cleaners (`System_Clean`, `Pattern_Clean`) log at debug level instead of `warnings.warn` (they already write to `Sapphire_Info.txt`); `Stats.KB_Dist` masks zero bins instead of emitting divide warnings; `data.synthetic` uses `FixCom` per ASE ≥ 3.28.

## Decisions (2026-08-29)
- **Mass_Activity volcano apex:** branches switch at their intersection, GCN = 8.096 (Rob, 2026-08-29). Literal 8.33 available via `apex=8.33`.
- **ML potentials:** scrap `tensorflow`/`mir-flare`/`ray`; adopt MACE foundation models as an optional extra (`pip install -e .[mlpot]`), see `Potentials/MLCalculator.py` and Tutorial 09.

## Phase 8 — legacy corners, performance, API (2026-08-29)

### Performance (Task 4) — 2-frame bimetallic benchmark 17.7 s → 1.7 s, outputs bit-identical except the aGCN fix below
- `DistFuncs.Euc_Dist`/`Hetero` → `scipy.spatial.distance.pdist/cdist` (13 M Python calls to `distance()` removed); `RDF.calculate` → `bincount`; `Kernels` Gaussian/Epanechnikov → chunked broadcasting, Uniform → sorted search; `Adjacent.calculate_adj` → `squareform`; `FrameSignature.T` → small DFS (networkx only for the rare cyclic bond graphs). Verified against a pickled baseline of 64 quantities.
- **aGCN was wrong.** The original `agcn_generator` walked `scipy.sparse.find` output assuming column-grouped rows; `find` sorts row-major, so each atom's own CN was summed: aGCN = CN²/12 (up to 16.3 in a liquid; a (111) terrace gave 6.75 instead of 7.5). Now aGCN_i = Σ_{j∈N(i)} CN_j / 12 exactly (tests: terrace 7.5, Ih vertex 4.33, bulk 12). Surface-area and surface-atom counts, which derive from aGCN, change accordingly. Whether the 2020 scipy behaved identically could not be verified — treat older Sapphire aGCN output with caution.
- CNA signature columns are now consistent across frames: `Process` carries a running masterkey so each frame's row is a prefix-compatible extension of the previous one; `Exec/Masterkey` holds the final key.

### Legacy corners (Task 3)
- `IO/Output.py`, `Graphing/Read_Plot.py`, `Graphing/Plotter.py` → `legacy/` (unused since the file-output refactor).
- `Potentials/GuptaPotential.py` — RGL energies (numpy) + `GuptaCalculator` (ASE, numerical forces); tests reproduce the Cleri–Rosato cohesive energies of Au/Ag/Cu/Ni/Pt/Pd within 6 %.
- `Light/Epsilon_DFT`: Au table has 150 n but 149 k values (source omission) — classes now truncate to the common range with a warning; tests check plasmonic ε′ < 0 for Ag/Au/Cu and the Ag interband edge. `ClassSpec` (pyGDM2) remains untested without the optional dependency.
- `Process` boilerplate docstrings replaced; documentation coverage is now tracked by the API reference rather than by a spreadsheet.

### API (Task 5)
- `Sapphire.api.Config` (dataclass, validated against `Utilities.Supported`, TOML round-trip) and `api.run(trajectory, out_dir, ...)` → `Reader`; bimetallic Homo/Hetero defaults inferred from the species; every run writes `sapphire_config.toml` for reproducibility.

### 1.1.0 packaging (2026-08-31)
- PyPI distribution name **`sapphire-nano`** (`sapphire` was taken in 2018); the import stays `import Sapphire`.
- Wheel/sdist verified: bundled tutorial samples and potentials included (`*.xyz.gz` added to package data); a **bare** `pip install sapphire-nano` imports and runs (`ruptures` import in `Stats` made lazy — it is the `changepoint` extra).
- `CITATION.cff` (preferred citation: Jones et al., Faraday Discuss. 2023, 242, 326); `release.yml` builds on tag `v*`, publishes to PyPI via trusted publishing and attaches the dist to a GitHub release.
- ASCII-logo strings made raw (SyntaxWarning on import gone).
