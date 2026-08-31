# Roadmap: machine-learned potentials

Sapphire analyses trajectories; it does not train force fields. This page records which
machine-learned potentials Sapphire interoperates with, and why the dependencies listed in
early versions were dropped.

**Current position (1.1.0):** MACE foundation models are supported through the optional
`mlpot` extra and `Sapphire.Potentials.MLCalculator`; see Tutorial 09. FLARE, TensorFlow and
Ray are not used.

## What the 2019–2021 dependencies were for

Early `setup.py` files listed `tensorflow`, `mir-flare`, `ray` and `numba`, none of which were
ever imported. Their intended roles were:

* **mir-flare** (FLARE, Vandermause/Kozinsky group) — sparse Gaussian-process force fields
  trained *on the fly* during MD, with Bayesian uncertainty to trigger DFT calls. The plan was
  to run Sapphire's structural descriptors (CNA, aGCN, PDDF) on FLARE trajectories and/or use
  them as inputs to inference.
* **tensorflow** — a neural-network alternative to the SVM in `CNA/Model.py` for classifying
  structures from CNA-pattern fingerprints, and possibly an NN potential.
* **ray** — parallelising `Process.calculate` across frames/nodes.
* **numba** — JIT for the O(N²) loops in `DistFuncs`, `Adjacent`, `FrameSignature`.

## Landscape (early 2026)
| Need | 2019 answer | Current answer | Why |
|---|---|---|---|
| ML interatomic potential for MD of metallic NPs | FLARE (GP, on-the-fly) | **MACE** (equivariant MPNN, `mace-torch`), plus **foundation models** `MACE-MP-0` / `MACE-MPA` trained on Materials Project — usable *zero-shot* for Au/Pt/Ag/Cu/Ni/Pd clusters with an ASE calculator; fine-tune on a few hundred DFT frames for accuracy. Alternatives: **NequIP/Allegro**, **ACE/pacemaker** (fast, linear), **SevenNet**, **Orb**. | Foundation models remove the "no training data" barrier entirely; all expose `ase.Calculator`, so `Sapphire.Tutorials.data.synthetic()` can swap EMT → MACE with one line. |
| Uncertainty-driven active learning | FLARE's GP variance | Ensembles (MACE committee), `ase` + `mace` active-learning scripts, or FLARE++ (still maintained) | If on-the-fly DFT is really wanted, FLARE++ is alive; otherwise fine-tuning a foundation model is cheaper. |
| Structure classification from CNA patterns | SVM (sklearn) → planned TF NN | Keep **scikit-learn** (SVC / gradient boosting); a NN adds nothing at this data scale. `CNA/Classify.py` now does this end-to-end. | Hundreds of labelled structures, ~10 features: classical ML is correct. |
| Descriptors for ML on structures | hand CNA/aGCN | **SOAP / ACE descriptors** via `dscribe` or `ase`-compatible libs feed directly into sklearn; Sapphire's CNA/aGCN remain excellent *interpretable* features alongside them. | |
| Parallel frames | ray | `multiprocessing` (already used) or `joblib`; ray only if you go multi-node. | |
| Hot loops | numba | vectorise with numpy/scipy (`cdist`, sparse) — done in `Morphology`, `FrameSignature`; numba remains a legitimate later step for `RDF`/`Kernels`. | |

## Decisions taken
1. **Dropped:** `tensorflow`, `mir-flare`, `ray` from the roadmap. Record here; they are already out
   of `pyproject.toml`.
2. **Adopted** (optional `mlpot` extra): `mace-torch` (+ `torch`) — heavy, hence optional.
   Provide `Sapphire.Potentials.MLCalculator("mace-mp-0")` returning an ASE calculator, and make
   `data.synthetic(calculator=...)` accept it. A tutorial then shows: build cluster → short MD
   with MACE-MP → Sapphire analysis → compare with the Gupta result. That is a genuinely modern
   demonstration and costs ~1 day.
3. `numba` remains a documented performance option, not a dependency.
4. Revisit FLARE++ only if a project needs on-the-fly DFT.

Sources to check when you decide: MACE (Batatia et al., NeurIPS 2022; MACE-MP-0, Batatia et al.
2023, arXiv:2401.00096), FLARE (Vandermause et al., npj Comput. Mater. 2020), ACE (Drautz 2019;
Lysogorskiy et al. 2021), NequIP (Batzner et al., Nat. Commun. 2022).
