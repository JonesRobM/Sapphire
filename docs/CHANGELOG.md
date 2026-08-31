# Changelog

## 1.1.0 — released 2026-08-31
Everything below (the 2026 restoration, Phases 1–8) constitutes release 1.1.0.

## 1.1.0.dev0 — restoration (2026-08-28)

### Repository
- Removed ~400 tracked build artefacts (`build/`, `dist/`, `.eggs/`, `*.egg-info`, `__pycache__`, `.ipynb_checkpoints`), a copied `.git/` internals dump (`HEAD`, `config`, `description`, `hooks/`, `info/`), duplicate `README`/`setup.py` under `main/`, and an Excel lock file.
- Removed the never-initialised `Raffy` submodule.
- New `legacy/` for quarantined code, `examples/` for driver templates, `tests/`, `docs/`.

### Packaging
- `setup.py` → `pyproject.toml`; `requires-python >= 3.10`; version `1.1.0.dev0`.
- Core dependencies reduced to `numpy scipy ase pandas networkx`; everything else is an extra (`plot`, `changepoint`, `ml`, `cna`, `light`, `quote`, `notebooks`, `dev`, `all`). See `docs/DEPENDENCY_NOTES.md`.

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
- `Process` boilerplate (bluepy) docstrings replaced; `docs/DOCSTRING_AUDIT.md` supersedes `Logs/DocStrings.xlsx` (which can be dropped).

### API (Task 5)
- `Sapphire.api.Config` (dataclass, validated against `Utilities.Supported`, TOML round-trip) and `api.run(trajectory, out_dir, ...)` → `Reader`; bimetallic Homo/Hetero defaults inferred from the species; every run writes `sapphire_config.toml` for reproducibility.

### 1.1.0 packaging (2026-08-31)
- PyPI distribution name **`sapphire-nano`** (`sapphire` was taken in 2018); the import stays `import Sapphire`.
- Wheel/sdist verified: bundled tutorial samples and potentials included (`*.xyz.gz` added to package data); a **bare** `pip install sapphire-nano` imports and runs (`ruptures` import in `Stats` made lazy — it is the `changepoint` extra).
- `CITATION.cff` (preferred citation: Jones et al., Faraday Discuss. 2023, 242, 326); `release.yml` builds on tag `v*`, publishes to PyPI via trusted publishing and attaches the dist to a GitHub release.
- ASCII-logo strings made raw (SyntaxWarning on import gone).
