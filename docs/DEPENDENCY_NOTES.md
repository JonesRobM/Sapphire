# Dependency notes

Recorded during the 2026 restoration so that original intent is not lost.

## Core (always installed)
`numpy`, `scipy`, `ase`, `pandas`, `networkx` — the only packages imported by
`Sapphire.Process` and `Sapphire.Post_Process`.

## Optional extras
| extra | packages | used by |
|---|---|---|
| `plot` | matplotlib, seaborn | `Sapphire.Graphing` |
| `changepoint` | ruptures | `Post_Process/Stats.py` |
| `ml` | scikit-learn, joblib | `CNA/Model.py` (SVM pattern classifier) |
| `cna` | ovito | `legacy/CNA/*` only — kept as an extra for resurrecting them |
| `light` | pygdm2 | `Sapphire.Light.ClassSpec` |
| `quote` | wikiquote | `Utilities/Initial.py` start-up quote (network call, guarded) |

## Removed from `install_requires` (were listed but never imported)
| package | original intent (best recollection) | status |
|---|---|---|
| `tensorflow` | ML classification of CNA patterns beyond the SVM | never wired in |
| `mir-flare` | FLARE Gaussian-process force fields for on-the-fly MD | never wired in |
| `ray` | parallel frame processing in `Process.calculate` | never wired in; `multiprocessing` is used instead |
| `numba` | JIT for `DistFuncs` / `Adjacent` hot loops | never wired in; a future performance phase could revisit |
| `sklearn` | typo-alias of `scikit-learn` | blocked on PyPI since 2023 |
| `networkx` (×3) | duplicate entries | collapsed to one |
