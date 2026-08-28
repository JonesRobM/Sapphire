# Changelog

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
