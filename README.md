# SAPPHIRE
Post - Processing Software
![Sapphire-logos_black](https://user-images.githubusercontent.com/52043020/154812758-c0184aa5-c2d0-4b5c-b23f-5955221d1d79.png)

Sapphire is a multi-faceted platform engineered to facilitate the design, characterisation, and classification of metallic nanoparticles in addition to computing their dynamics and energetics. Integrating these two aspects of nanoparticle simulation is the comprehensive set of analysis techniques offered by Sapphire which utilises pre-exisiting, highly-efficient, pythonic libraries such as ase (atomic simulation environment.

Please, also be aware that Sapphire is still in a perpetual state of development. Should you recover any strange results, contact one of the authors of this document with any information you are able to provide about the problem and a hot-fix will be implemented if necessary. Otherwise, feedback on ease of use is always appreciated.

## Installation

Sapphire supports Python 3.10+.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[plot,changepoint]"      # core + plotting + change-point analysis
# other extras: ml (SVM pattern classifier), light (pyGDM2 optics), quote, notebooks, dev, all
```

Verify:

```bash
python -c "import Sapphire; print(Sapphire.__version__)"
pytest            # needs the dev extra
```

## Quick start

```python
from Sapphire.api import run
from Sapphire.Tutorials import data

xyz = data.sample("AuPt", "work/")            # bundled Au80Pt20 melting trajectory (70 frames)
r = run(xyz, "work/out/", quantities=["pdf", "adj", "nn", "agcn", "cna_sigs"],
        frames=(0, 70, 7), statistics={"JSD": ["pdf"]})
r.load("agcn")                                 # (frames, atoms) atop generalised coordination numbers
```

Results are plain per-frame text files (`docs/FILE_CONTRACT.md`) readable with `Sapphire.IO.Reader`
or any other tool. The classic two-dictionary interface to `Sapphire.Process` is unchanged
(`examples/run_analysis.py`).

## Tutorials

Nine executable notebooks in `main/Sapphire/Tutorials/`, run in CI, rendered on the
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

## Repository layout

| Path | Contents |
|---|---|
| `main/Sapphire/` | the package (`Process`, `Post_Process`, `CNA`, `IO`, `Graphing`, `Potentials`, `Light`, `Utilities`, `Tutorials`) |
| `examples/` | driver-script templates |
| `tests/` | pytest suite (import sweep, numerics, end-to-end smoke run on Au561) |
| `docs/` | mkdocs site: guides, file contract, API reference, changelog (`mkdocs serve`) |
| `legacy/` | quarantined code kept for reference; not installed |

Please check out our group page to learn more about the exciting work we do at the nanoscale! 

http://balettogroup.org

Authors:

1. Robert M. Jones[^1].

2. Matteo Tibberi

3. Armand Aquier

4. Claudio Zeni

5. Francesca Baletto

[^1]: Robert.M.Jones@kcl.ac.uk
