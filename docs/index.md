# Sapphire

**Sapphire** turns molecular-dynamics trajectories of metallic nanoparticles into structure
descriptors and statistics: pair-distance distributions and cutoffs, adjacency and coordination,
atop generalised coordination numbers and the ORR mass-activity model built on them,
common-neighbour-analysis signatures and patterns, chemical ordering in alloys, shape and
morphology, and the divergences and change-points that locate transitions along a run.

* Runs on any trajectory ASE can read (xyz, ASE `.traj`, LAMMPS dumps, …).
* Writes plain-text per-frame results you can read back with `Sapphire.IO.Reader` or any tool
  (see the [file contract](FILE_CONTRACT.md)).
* Nine executable [tutorials](tutorials/01_Build_and_Inspect.ipynb) build from a single cluster to
  ensembles of runs.

```python
from Sapphire.api import run
reader = run("movie.xyz", "out/", quantities=["pdf", "adj", "nn", "agcn", "cna_sigs"], frames=(0, 100, 5))
agcn = reader.load("agcn")          # (frames, atoms)
```

Developed in the Baletto group (King's College London) for the study of nanoalloys; restored and
modernised in 2026. GPL-3.0.
