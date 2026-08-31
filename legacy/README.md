# legacy/

Code quarantined during the 2026 modernisation and kept for reference.
Nothing here is packaged or imported by `Sapphire`. Kept because the scientific
intent may be worth resurrecting once remembered.

| Path | Why it is here |
|---|---|
| `Emerald/Constants.py`, `Emerald/main.py` | Morphology prototype (CN/GCN, surface–core peeling, thickness, volume, faceting). **Superseded 2026-08-28 by `Sapphire/Post_Process/Morphology.py`** (vectorised, tested against the Mackay shells of Au561). Kept for the record. |
| `CNA/main.py` | Armand's batch CNA-pattern → SVC pipeline (OVITO-based, interactive prompts). **Superseded 2026-08-28 by `Sapphire/CNA/Classify.py`** built on Sapphire's own `FrameSignature`. |
| `CNA/FrameSignature_Ovito.py` | ovito-only variant superseded by `Sapphire/CNA/FrameSignature.py` (ASE-based). |
| `IO/ExtendXYZ.py` | Byte-identical duplicate of `Sapphire/Utilities/ExtendXYZ.py`. |
| `IO/*.sh` | HPC job-submission templates for a specific cluster. |
| `IO/Output.py` | Pre-1.0 metadata writer (22 broad excepts); not imported anywhere since the file-output refactor. Moved 2026-08-29. |
| `Graphing/Read_Plot.py`, `Graphing/Plotter.py` | PhD-era plotting scripts for specific Au/Pt runs (hard-coded species and directory layouts); superseded by `Graphing.Read_Meta` + `Plot_Funcs` on `IO.Reader`. Moved 2026-08-29. |
| `Tutorials/*` | The 2022 notebooks (see `Tutorials/ReadMe.md` for the current course). |
