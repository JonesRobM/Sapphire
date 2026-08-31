# The Sapphire run-directory contract

Sapphire's analysis and plotting read **files**, not Python objects. Any tool that writes files in
this layout is a first-class data source for `Sapphire.IO.Reader`, `Process.analyse`, and
`Graphing` — whether the numbers came from Sapphire, LAMMPS, OVITO, or a notebook.

```
<base_dir>/
  Time_Dependent/<Quantity>[<Species>]        one line per frame:  <frame> <v1> <v2> ...
  Time_Dependent/<Quantity>File<frame>        one dense matrix per frame (rows = atoms)
  Time_Dependent/Stats/<Statistic><quantity>  one line per frame:  <frame> <value>
  Adjacency/File<frame>                       dense N×N 0/1 adjacency for the whole cluster
  Adjacency/HomoAdj<Species>File<frame>       N_X×N_X adjacency within species X
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
* Rows may be **ragged** across frames when a per-frame vocabulary grows (CNA signatures);
  consumers right-pad with zeros.
* Species-resolved files carry the symbol as a suffix (`HomoPDFAu`, `HomoCoMDistPt`). The
  Reader key is the table name plus suffix (`hopdfAu`).

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
