# Getting started

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[plot,changepoint]"
```

Extras: `ml` (CNA structure classifier), `mlpot` (MACE foundation-model potentials, pulls in
PyTorch), `light` (pyGDM2 optics), `quote`, `notebooks`, `dev`, `docs`, `all`.

## First analysis

```python
from Sapphire.api import run
from Sapphire.Tutorials import data

xyz = data.sample("AuPt", "work/")                       # bundled 70-frame Au80Pt20 melt
r = run(xyz, "work/out/", quantities=["pdf", "adj", "nn", "agcn", "cna_sigs", "mix"],
        frames=(0, 70, 7), statistics={"JSD": ["pdf"]})
print(sorted(r.available()))
print(r.load("mix"))                                     # mixing parameter per frame
```

The same run expressed as the classic two dictionaries is in `examples/run_analysis.py`;
`examples/from_lammps.py` ingests a LAMMPS dump; `examples/make_figures.py` renders the standard
figure set.

## Where results go

```
out/
  Time_Dependent/   one line per frame per quantity      (NN, AGCN, PDF, RCut, Mixing, ...)
  Time_Dependent/Stats/   statistics from analyse()      (JSDpdf, Collectivity, ...)
  Adjacency/        one dense matrix per frame
  CNA/              Signatures, Patterns
  Exec/             Masterkey
  sapphire_config.toml, Sapphire_Info.txt, Sapphire_Errors.log
```

`Sapphire.IO.Reader.Reader("out/").load("nn")` gives a `(frames, atoms)` array; `load_all()` a dict.
