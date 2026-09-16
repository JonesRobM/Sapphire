# The Sapphire run-directory contract

Sapphire's analysis and plotting read **files**, not Python objects. Any tool that writes files in
this layout is a first-class data source for `Sapphire.IO.Reader`, `Process.analyse`, and
`Graphing` — whether the numbers came from Sapphire, LAMMPS, OVITO, or a notebook.

```
<base_dir>/
  Time_Dependent/<Quantity>[<Species>]        one line per frame:  <frame> <v1> <v2> ...
  Time_Dependent/Stats/<Statistic><quantity>  one line per frame:  <frame> <value>
  Time_Dependent/HeAdjFile<frame>[.npz]       N_A×N_B adjacency between the two species
  Adjacency/File<frame>[.npz]                 N×N 0/1 adjacency for the whole cluster
  Adjacency/HomoAdj<Species>File<frame>[.npz] N_X×N_X adjacency within species X
  CNA/Signatures                              <frame> <count per masterkey entry ...>
  CNA/Patterns                                <frame> <per-atom pattern tuple repr ...>
  Exec/Masterkey                              "000 100 200 211 ... 555 666" (r s t triples)
  Sapphire_Info.txt / Sapphire_Errors.log     human-readable run log / caught exceptions
```

## Line format
* First token: integer frame index (the trajectory frame, not a running counter).
* Remaining tokens, one of:
  * scalars — `12 12 11 …` or `3.1959…`;
  * vectors — `[x y z]` (numpy repr, any whitespace inside the brackets);
  * CNA patterns — `((12, (5, 5, 5)),)` or `((2, (5, 5, 5)), (10, (4, 2, 2)))`;
  * bare words — masterkey labels.
* Rows are **rectangular**. Where a per-frame vocabulary grows (CNA signatures), the table is
  laid out against the final `Exec/Masterkey` at the end of the run and short rows are padded
  with zeros, so `np.loadtxt` reads it directly. Runs written before 1.3.0 may be ragged;
  consumers that right-pad defensively stay correct for both.
* Species-resolved files carry the symbol as a suffix (`HomoPDFAu`, `HomoCoMDistPt`). The
  Reader key is the table name plus suffix (`hopdfAu`).

## Per-frame matrices

A matrix file is identified by its **name**, not its directory: anything ending `File<frame>`
(optionally `.npz`). The full and per-species matrices live under `Adjacency/`, the hetero one
under `Time_Dependent/`.

Two interchangeable encodings:

| | file | contents |
|---|---|---|
| sparse (**default**) | `File7.npz` | `scipy.sparse.save_npz` of a CSR 0/1 matrix |
| dense text | `File7` | one row per atom, `0`/`1` space separated, `%d` |

Adjacency is ~12 non-zeros per row, so dense text spends O(N²) bytes on O(N) information —
8 MB per frame at N = 2000, about 160 GB over a 20 000-frame run, against roughly 1 GB sparse.
`Sapphire.IO.Reader` loads either transparently and always hands back a dense
`numpy` array indexed `[frame, i, j]`, so a consumer using the Reader need not care which is
on disk.

To choose the encoding, or to convert:

```bash
sapphire run movie.xyz -o run/ --adj-format text   # write dense text instead
sapphire expand run/ -o dense/                     # sparse -> dense text, layout preserved
sapphire expand run/ --frames 0:1:1 --stdout       # one matrix to stdout
```

## Table of names
`IO/OutputInfoFull.py`, `OutputInfoHomo.py`, `OutputInfoHetero.py`, `OutputInfoExec.py` map
metadata keys (`pdf`, `rcut`, `hocomAu`, …) to `Dir` + `File`. Add a dictionary there to teach
the Reader a new quantity; anything under `Time_Dependent/Stats/` needs no entry.

## Minimal external producer
```python
import numpy as np
with open("run/Time_Dependent/MyQuantity", "w") as f:
    for frame, values in enumerate(my_per_frame_arrays):
        f.write(f"{frame} " + " ".join(map(str, values)) + "\n")
```
Then `Reader("run/").load("myquantity")` — after adding `myquantity = {'Dir': 'Time_Dependent/',
'File': 'MyQuantity', ...}` to `OutputInfoFull.py`, or write it under `Time_Dependent/Stats/`
to skip that step.
