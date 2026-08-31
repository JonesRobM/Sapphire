# SAPPHIRE

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22211283.svg)](https://doi.org/10.5281/zenodo.22211283)
[![PyPI](https://img.shields.io/pypi/v/sapphire-nano)](https://pypi.org/project/sapphire-nano/)
[![CI](https://github.com/JonesRobM/Sapphire/actions/workflows/ci.yml/badge.svg)](https://github.com/JonesRobM/Sapphire/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-jonesrobm.github.io%2FSapphire-blue)](https://jonesrobm.github.io/Sapphire/)
[![License: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0-lightgrey)](LICENSE)

![Sapphire-logos_black](https://user-images.githubusercontent.com/52043020/154812758-c0184aa5-c2d0-4b5c-b23f-5955221d1d79.png)

**Sapphire** is a post-processing environment for the structural characterisation of metallic
nanoparticles and nanoalloys from molecular-dynamics trajectories. It turns frames into
physics: pair-distance distributions and per-frame nearest-neighbour cutoffs, adjacency and
coordination, atop generalised coordination numbers (aGCN) and the GCN-based oxygen-reduction
mass-activity model, common-neighbour-analysis (CNA) signatures and patterns with an SVM
structure classifier, chemical ordering in alloys (mixing parameter, LAE, per-species
neighbour counts), shape and shell-by-shell morphology, and the distribution divergences and
change-point statistics that locate melting and other transitions along a run.

It reads anything [ASE](https://wiki.fysik.dtu.dk/ase/) can read and writes plain per-frame
text files that any tool can consume. Documentation, rendered tutorials and the full API
reference live at **https://jonesrobm.github.io/Sapphire/**.

## Citing Sapphire

If Sapphire contributes to your work, please cite the method paper:

> R. M. Jones, K. Rossi, C. Zeni, M. Vanzan, I. Vasiljevic, A. Santana-Bonilla and F. Baletto,
> *Structural characterisation of nanoalloys for (photo)catalytic applications with the Sapphire
> library*, **Faraday Discussions**, 2023, **242**, 326–352.
> [doi:10.1039/D2FD00097K](https://doi.org/10.1039/D2FD00097K)

```bibtex
@article{Jones2023Sapphire,
  author  = {Jones, Robert M. and Rossi, Kevin and Zeni, Claudio and Vanzan, Mirko and
             Vasiljevic, Igor and Santana-Bonilla, Alejandro and Baletto, Francesca},
  title   = {Structural characterisation of nanoalloys for (photo)catalytic applications
             with the Sapphire library},
  journal = {Faraday Discussions},
  year    = {2023},
  volume  = {242},
  pages   = {326--352},
  doi     = {10.1039/D2FD00097K},
}
```

To cite the software itself (a specific archived version), use the Zenodo DOI:
[10.5281/zenodo.22211283](https://doi.org/10.5281/zenodo.22211283) resolves to the latest
release; v1.1.0 is [10.5281/zenodo.22211284](https://doi.org/10.5281/zenodo.22211284).
GitHub's "Cite this repository" button (from `CITATION.cff`) gives both formats.

The GCN-based mass-activity model implemented in `Post_Process.Mass_Activity` follows
Rossi, Asara & Baletto, *ChemPhysChem* 2019, **20**, 3037
([doi:10.1002/cphc.201900564](https://doi.org/10.1002/cphc.201900564)), building on Rück,
Bandarenka, Calle-Vallejo & Gagliardi, *J. Phys. Chem. Lett.* 2018, **9**, 4463.

## Installation

Sapphire supports Python 3.10+.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[plot,changepoint]"      # core + plotting + change-point analysis
```

Extras: `ml` (CNA structure classifier), `mlpot` (MACE foundation-model potentials; pulls in
PyTorch), `light` (pyGDM2 optics), `quote`, `notebooks`, `docs`, `dev`, `all`. Verify with
`python -c "import Sapphire; print(Sapphire.__version__)"` and, with the `dev` extra, `pytest`.

## Quick start

```python
from Sapphire.api import run
from Sapphire.Tutorials import data

xyz = data.sample("AuPt", "work/")            # bundled Au80Pt20 melting trajectory (70 frames)
r = run(xyz, "work/out/", quantities=["pdf", "adj", "nn", "agcn", "cna_sigs"],
        frames=(0, 70, 7), statistics={"JSD": ["pdf"]})
r.load("agcn")                                 # (frames, atoms) atop generalised coordination numbers
```

Results are plain per-frame text files (see the
[file contract](https://jonesrobm.github.io/Sapphire/FILE_CONTRACT/)) readable with
`Sapphire.IO.Reader` or any other tool. The classic two-dictionary interface to
`Sapphire.Process` is unchanged (`examples/run_analysis.py`); `examples/from_lammps.py`
ingests a LAMMPS dump.

## Tutorials

Nine executable notebooks in `main/Sapphire/Tutorials/`, run in CI and rendered on the
[documentation site](https://jonesrobm.github.io/Sapphire/):

| # | Topic |
|---|---|
| 01 | Build a cluster; CN/GCN; surface–core peeling (`Morphology`) |
| 02 | Pair-distance KDE, RDF, deriving the cutoff |
| 03 | Adjacency, aGCN and the ORR mass-activity volcano |
| 04 | CNA signatures, patterns, structure classifier |
| 05 | Bimetallic trajectory with `Process`, `Reader`, `Graphing` |
| 06 | Divergences, collectivity, change-point detection |
| 07 | Shape: inertia, radii of gyration, radial density |
| 08 | Ensemble averaging over runs |
| 09 | MD with a MACE foundation-model potential |

The bundled samples are down-sampled from four 14 ns bimetallic MD data sets published as a
[GitHub release](https://github.com/JonesRobM/Sapphire/releases/tag/data-v1); fetch the full
trajectories with `Sapphire.Tutorials.data.fetch(...)`.

## Repository layout

| Path | Contents |
|---|---|
| `main/Sapphire/` | the package (`api`, `Process`, `Post_Process`, `CNA`, `IO`, `Graphing`, `Potentials`, `Light`, `Utilities`, `Tutorials`) |
| `examples/` | driver-script templates |
| `tests/` | pytest suite (import sweep, numerics, end-to-end runs on bundled data) |
| `docs/` | mkdocs site: guides, file contract, API reference, changelog (`mkdocs serve`) |
| `legacy/` | earlier code kept for reference; not installed |

## Authors and acknowledgements

Sapphire was written in the [Baletto group](http://balettogroup.org) at King's College London.

* **Robert M. Jones** — lead author and maintainer ([Robert.M.Jones@kcl.ac.uk](mailto:Robert.M.Jones@kcl.ac.uk))
* **Francesca Baletto** — principal investigator
* **Claudio Zeni**, **Kevin Rossi**, **Mirko Vanzan**, **Igor Vasiljevic**,
  **Alejandro Santana-Bonilla** — co-authors of the Sapphire paper
* **Matteo Tibberi**, **Armand Aquier** — contributors to early versions of the CNA-pattern
  pipeline and morphology tools (their originals are preserved under `legacy/`)

The library builds on [ASE](https://wiki.fysik.dtu.dk/ase/), numpy/scipy, and optionally
[ruptures](https://centre-borelli.github.io/ruptures-docs/), scikit-learn,
[pyGDM2](https://wiechapeter.gitlab.io/pyGDM2-doc/) and [MACE](https://mace-docs.readthedocs.io/).

Sapphire is research software in active development: if something looks wrong, please open an
issue with whatever you can share — trajectories, logs, or just the surprise.

## Licence

GPL-3.0 (see `LICENSE`).
