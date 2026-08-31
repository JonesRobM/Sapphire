# Docstring audit (from `Logs/DocStrings.xlsx`, 2021)

The spreadsheet was a 2021 progress tracker: one sheet per sub-package with columns
*Sub-Module · Description · Status · Functions · Written · % complete*. Only the `Post_Process`
and `CNA` sheets carried real entries; the others were unfilled copies of the CNA template.
Recorded here so the binary can be dropped from the repository.

| Sub-module (2021) | Description recorded | 2021 status | 2026 status |
|---|---|---|---|
| Post_Process/Adjacent | adjacency matrix in sparse form; adjacency-dependent quantities (aGCN) | not started | documented; vectorised 2026-08-29 |
| Post_Process/Kernels | kernel functions for KDE of distributions | in progress | documented; vectorised 2026-08-29 |
| Post_Process/AtomicEnvironment | semi-local (hetero) environment properties | "complete" (was a stub) | rewritten 2026-08-28 |
| Post_Process/DistFuncs | distance-based building blocks | in progress | documented; vectorised 2026-08-29 |
| CNA/FramePattern | parse CNAP data for the SVC | not started | superseded by `CNA/Classify.py` |
| CNA/FrameSignature | per-frame CNA signatures | in progress | documented; `T()` rewritten without networkx |
| CNA/main | Armand's pattern-dictionary pipeline | complete | `legacy/`; superseded by `CNA/Classify.py` |
| CNA/Model | SVC model | in progress | kept; `Classify.Classifier` is the supported path |
| CNA/PatternWriter | write pattern movies | in progress | unchanged |

Coverage today is tracked by the API reference (`mkdocs build` warns on undocumented public
objects) rather than by hand.
