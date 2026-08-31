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

Developed in the [Baletto group](http://balettogroup.org) (King's College London); restored and
modernised in 2026. GPL-3.0.

**Citing Sapphire:** R. M. Jones, K. Rossi, C. Zeni, M. Vanzan, I. Vasiljevic, A. Santana-Bonilla
and F. Baletto, *Structural characterisation of nanoalloys for (photo)catalytic applications with
the Sapphire library*, Faraday Discuss., 2023, **242**, 326–352,
[doi:10.1039/D2FD00097K](https://doi.org/10.1039/D2FD00097K). Software archive:
[doi:10.5281/zenodo.22211283](https://doi.org/10.5281/zenodo.22211283).
