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
from Sapphire import Process

System = {"base_dir": "./run/", "movie_file_name": "movie.xyz", "extend_xyz": None,
          "Homo": None, "Hetero": False,
          "Start": 0, "End": 100, "Step": 1, "Skip": 1, "UniformPDF": False, "Band": 0.05}
Quantities = {"Full": {"pdf": None, "adj": None, "nn": None, "agcn": None, "cna_sigs": None}}

Process.Process(System=System, Quantities=Quantities)
# results: ./run/Time_Dependent/<Quantity>  (one line per frame)
```

Complete templates live in `examples/`; worked notebooks in `main/Sapphire/Tutorials/`.
Tutorial MD data is fetched on demand — see `main/Sapphire/Tutorials/Data/ReadMe.md`.

## Repository layout

| Path | Contents |
|---|---|
| `main/Sapphire/` | the package (`Process`, `Post_Process`, `CNA`, `IO`, `Graphing`, `Potentials`, `Light`, `Utilities`, `Tutorials`) |
| `examples/` | driver-script templates |
| `tests/` | pytest suite (import sweep, numerics, end-to-end smoke run on Au561) |
| `docs/` | restoration plan, changelog, dependency notes, history-rewrite runbook |
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
