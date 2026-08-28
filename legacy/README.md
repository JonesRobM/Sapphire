# legacy/

Code quarantined during the 2026 restoration (see `docs/RESTORATION_PLAN.md`).
Nothing here is packaged or imported by `Sapphire`. Kept because the scientific
intent may be worth resurrecting once remembered.

| Path | Why it is here |
|---|---|
| `Emerald/Constants.py`, `Emerald/main.py` | A script (runs an analysis of `Au887_MDh453_rlx.xyz` at import time) with no imports; 232 undefined names. Not a module. |
| `CNA/main.py` | ovito driver importing `cna_signatures`, `cna_patterns`, `xyz_io` — modules that never existed in the repo. |
| `CNA/FrameSignature_Ovito.py` | ovito-only variant superseded by `Sapphire/CNA/FrameSignature.py` (ASE-based). |
| `IO/ExtendXYZ.py` | Byte-identical duplicate of `Sapphire/Utilities/ExtendXYZ.py`. |
| `IO/*.sh` | HPC job-submission templates for a specific cluster. |
