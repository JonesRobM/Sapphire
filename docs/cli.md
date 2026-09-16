# Command line

Installing the package puts a `sapphire` command on your path. It drives the same
`Sapphire.api.Config` the Python interface uses, so anything you can do from a notebook you can
do from a shell script, and a run is reproducible from one TOML file.

Two things shape its behaviour, both from running unattended on a batch node: **a run that
could not compute what you asked for exits non-zero**, and **nothing optional is allowed to
fail a run**.

## sapphire run

```bash
sapphire run movie.xyz -o out/ -q pdf,adj,nn,agcn,cna_sigs --frames 0:1000:10
sapphire run --config out/sapphire_config.toml          # repeat an earlier run
sapphire run movie.xyz --dry-run                        # resolve the plan, compute nothing
```

| Option | |
|---|---|
| `-o, --out-dir DIR` | where results go (default `sapphire_run/`) |
| `-q, --quantities A,B,C` | see `sapphire quantities` |
| `--homo A,B` / `--hetero A,B` | per-species and between-species quantities |
| `--species Au,Pt` | default: inferred from the first frame |
| `--frames START:END:STEP` | any field may be empty — `::10`, `5:`, `:100:` |
| `-j, --jobs N` | analyse frames across N processes |
| `--adj-format {npz,text}` | adjacency encoding on disk (default `npz`) |
| `--band A` | KDE bandwidth in ångström for the cutoff (default 0.05) |
| `--no-quote` | skip the quote lookup in the log header |
| `--strict` | stop at the first failing quantity |
| `--no-overwrite` | keep results from a previous run |
| `--dry-run` | print the resolved plan and exit |

### Prerequisites are resolved for you

The cutoff that decides whether two atoms are adjacent is read off the pair-distance
distribution, so everything downstream of the adjacency matrix needs `pdf` in the same run.
Ask for `cna_sigs` on its own and the missing pieces are added, and said so:

```console
$ sapphire run movie.xyz -q cna_sigs --dry-run
sapphire: also computing adj, pdf (required by what you asked for)
quantities : pdf, adj, cna_sigs
```

`sapphire quantities` lists every quantity with its prerequisites.

### Exit codes

| | |
|---|---|
| `0` | everything asked for reached disk |
| `1` | the run finished but a quantity produced no output — see `Sapphire_Errors.log` |
| `2` | bad arguments, missing trajectory, unknown quantity |
| `130` | interrupted |

### Running offline

`Process` decorates its log header with a quote fetched over the network. On a node with no
outbound route that is a needless delay, so it is bounded by a short timeout, never raised, and
can be turned off:

```bash
sapphire run movie.xyz --no-quote          # or export SAPPHIRE_NO_QUOTE=1
```

`SAPPHIRE_QUOTE_TIMEOUT` sets the timeout in seconds (default 3).

## Analysing frames in parallel

```bash
sapphire run movie.xyz -o out/ -j 8
```

Frames are independent, and every cross-frame quantity — collectivity, concertedness, JSD,
autocorrelation, change points — is computed afterwards from the finished per-frame results.
Workers take contiguous chunks, and the merge restores frame order, so **output is identical to
a serial run**; only the wall clock changes.

Speedup improves with trajectory length, because each worker's start-up is amortised over more
frames. On 350 frames of the bundled 1415-atom sample: 39.0 s at `-j 1`, 12.7 s at `-j 4`,
8.6 s at `-j 8`.

## sapphire expand

Adjacency is stored sparse, which is about 165x smaller — roughly 1 GB rather than 160 GB over
a 20 000-frame run at 2000 atoms. When you want to read a matrix by eye or feed one to a tool
that expects text, expand it:

```bash
sapphire expand run/ -o dense/           # every matrix, run layout preserved
sapphire expand run/ --frames 2:6:2      # a subset, written next to the npz files
sapphire expand run/ --frames 0:1:1 --stdout | less
```

Output is byte-identical to what `--adj-format text` would have written. Full, per-species and
hetero matrices are all covered, wherever they live in the run directory.

`Sapphire.IO.Reader` reads both encodings and always hands back dense arrays, so this is only
needed for tools that bypass the Reader.
