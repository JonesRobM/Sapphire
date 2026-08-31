# Sapphire tutorials

A graded course; each notebook runs offline on bundled data in under a minute or two.

| # | Notebook | What you learn |
|---|---|---|
| 01 | Build_and_Inspect | icosahedra, CN/GCN, surface–core peeling, faceting, volume (`Morphology`) |
| 02 | PDDF_and_Cutoff | pair-distance KDE, RDF, deriving R_cut and why it is per-frame |
| 03 | Adjacency_CN_aGCN | adjacency matrix, aGCN, the ORR volcano and mass activity |
| 04 | CNA_Signatures_and_Patterns | (r,s,t) signatures, patterns, drawing them, SVC structure classifier |
| 05 | Trajectory_Analysis_with_Process | `Process` on a bimetallic melt, `Reader`, chemical ordering, `Graphing` |
| 06 | Statistics_and_Changepoints | JSD / KL / KS along a run, collectivity, PELT change-points → T_melt |
| 07 | Shape_MoI_and_Radii | inertia ratio, radii of gyration per species, radial density |
| 08 | Ensemble_Averaging | several runs → mean ± std with `Read_Meta`, temperature averaging |
| 09 | ML_Potentials | MD with a MACE foundation model (falls back to EMT), analysed unchanged |

Notebooks are executed in CI (`.github/workflows/notebooks.yml`); run one yourself with
`jupyter nbconvert --to notebook --execute Tutorial.ipynb`. The previous notebooks live in
`legacy/Tutorials/` for reference.
