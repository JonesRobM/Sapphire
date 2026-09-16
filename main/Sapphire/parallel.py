"""Analyse a trajectory's frames across several processes.

Every per-frame quantity is independent -- ``Process.result_cache`` is per-frame scratch --
and every cross-frame quantity (collectivity, concertedness, JSD, autocorrelation,
changepoints) already lives in :meth:`Process.analyse`, which reads finished results back
off disk. So the frame loop parallelises without disturbing the science.

Two things stop it being a bare ``Pool.map``:

* **Shared output files.** ``calculate`` appends a line per frame to files like
  ``Time_Dependent/AGCN``. Concurrent appends from several processes interleave and land
  in arbitrary order. Each worker therefore writes into a private directory and the
  results are merged afterwards, in frame order.

* **The CNA masterkey.** ``Process`` carries signatures forward between frames
  (``self.Masterkey = tuple(cna.Sigs.keys())``), so a frame's signature columns depend on
  the frames before it. Workers cannot share that. Each starts from the same bulk key and
  extends it independently; because the key only ever *appends*, a row of width ``w``
  corresponds to the first ``w`` entries of that worker's final key, which is enough to
  recover every (signature, count) pair and re-emit them against a single merged key.

The merged signature table is rectangular. A serial run's is ragged whenever the key
grows mid-run, so the two agree exactly when no new signature appears -- the usual case,
since the bulk key already covers the twenty standard signatures -- and the parallel one
is the better-formed of the two when they differ.
"""
from __future__ import annotations

import os
import pathlib
import re
import shutil
from concurrent.futures import ProcessPoolExecutor

# A per-frame matrix names its frame, and may be sparse (.npz) or dense text. These are
# identified by NAME rather than by directory: the full and homo matrices live under
# Adjacency/, but the hetero ones (HeAdjFile7.npz) are written into Time_Dependent/
# alongside ordinary per-frame text. Merging one of those as text corrupts it.
_PER_FRAME = re.compile(r".*File\d+(?:\.npz)?$")
_LOG_FILES = ('Sapphire_Info.txt', 'Sapphire_Errors.log')


def plan_chunks(start, end, step, jobs):
    """Split ``range(start, end, step)`` into at most ``jobs`` contiguous (start, end, step).

    Contiguous rather than strided so each worker reads one slice of the trajectory
    instead of seeking across the whole file.
    """
    frames = list(range(start, end, step))
    if not frames:
        return []
    jobs = max(1, min(jobs, len(frames)))
    per, extra = divmod(len(frames), jobs)
    chunks, lo = [], 0
    for k in range(jobs):
        n = per + (1 if k < extra else 0)
        block = frames[lo:lo + n]
        lo += n
        # end is exclusive and must not swallow the next chunk's first frame
        chunks.append((block[0], block[-1] + 1, step))
    return chunks


def _run_chunk(payload):
    """Analyse one chunk in its own directory. Runs in a worker process."""
    system, quantities, pattern_input, strict = payload
    from Sapphire import Process
    Process.Process(System=system, Quantities=quantities, Pattern_Input=pattern_input,
                    strict=strict, overwrite=True)
    return system['base_dir']


def _frame_of(line):
    """Leading frame number of an output line, or None if it does not start with one."""
    head = line.split(' ', 1)[0]
    try:
        return int(head)
    except ValueError:
        return None


def _merge_frame_indexed(sources, target):
    """Concatenate per-frame lines from every worker and restore frame order."""
    numbered, unnumbered = [], []
    for src in sources:
        try:
            text = src.read_text()
        except UnicodeDecodeError as exc:
            # Something binary reached the text merge -- a per-frame matrix whose name the
            # _PER_FRAME pattern did not recognise. Concatenating it would corrupt it, so
            # say which file rather than surfacing a bare decode error.
            raise ValueError(
                f"{src.name} is not text and cannot be merged line by line; "
                f"if it is a per-frame matrix, _PER_FRAME must match its name"
            ) from exc
        for line in text.splitlines():
            if not line.strip():
                continue
            f = _frame_of(line)
            (numbered if f is not None else unnumbered).append((f, line))
    numbered.sort(key=lambda pair: pair[0])
    target.parent.mkdir(parents=True, exist_ok=True)
    body = [line for _, line in unnumbered] + [line for _, line in numbered]
    target.write_text(''.join(line + '\n' for line in body))


def _read_masterkey(path):
    return path.read_text().split() if path.is_file() else []


def rewrite_signatures(worker_dirs, out_dir):
    """Re-emit signature counts from one or more runs against a single masterkey.

    ``Process`` carries the CNA masterkey forward between frames, so a frame's row is only
    as wide as the key was when that frame was written: the table is ragged whenever a new
    signature turns up mid-run, and ``np.loadtxt`` cannot read it. Since the key only ever
    *appends*, a row of width w is labelled by the first w signatures of that run's final
    key, which is enough to recover every (signature, count) pair and lay them out
    rectangularly against the union key.

    Used for both a parallel merge (several source directories) and a serial run's
    finalisation (one directory, rewritten in place). Every source is read before anything
    is written, so in-place is safe.
    """
    worker_dirs = [pathlib.Path(w) for w in worker_dirs]
    out_dir = pathlib.Path(out_dir)
    rows, merged_key = {}, []
    for wd in worker_dirs:
        key = _read_masterkey(wd / 'Exec' / 'Masterkey')
        for name in key:
            if name not in merged_key:
                merged_key.append(name)
        sig_file = wd / 'CNA' / 'Signatures'
        if not sig_file.is_file():
            continue
        for line in sig_file.read_text().splitlines():
            if not line.strip():
                continue
            parts = line.split()
            frame, counts = int(parts[0]), parts[1:]
            rows[frame] = dict(zip(key[:len(counts)], counts))

    if not rows and not merged_key:
        return
    (out_dir / 'Exec').mkdir(parents=True, exist_ok=True)
    # Reproduce the serial file exactly. When every worker agreed on the key -- the usual
    # case, since the bulk key already covers the standard signatures -- copy one verbatim
    # rather than guessing at the writer's separators. Only a key that actually grew needs
    # rebuilding, and then we follow the observed convention of a trailing separator.
    sources = [wd / 'Exec' / 'Masterkey' for wd in worker_dirs]
    present = [p for p in sources if p.is_file()]
    texts = {p.read_text() for p in present}
    if len(texts) == 1:
        (out_dir / 'Exec' / 'Masterkey').write_bytes(present[0].read_bytes())
    else:
        (out_dir / 'Exec' / 'Masterkey').write_text(''.join(k + ' ' for k in merged_key))
    if rows:
        (out_dir / 'CNA').mkdir(parents=True, exist_ok=True)
        with open(out_dir / 'CNA' / 'Signatures', 'w') as f:
            for frame in sorted(rows):
                counts = rows[frame]
                f.write(str(frame) + ' ' + ' '.join(counts.get(k, '0') for k in merged_key) + '\n')


def merge(worker_dirs, out_dir, skip_names=()):
    """Combine worker output directories into ``out_dir``.

    ``skip_names`` are files in a worker directory that are inputs rather than results
    (the trajectory each worker was pointed at).
    """
    worker_dirs = [pathlib.Path(w) for w in worker_dirs]
    out_dir = pathlib.Path(out_dir)
    skip = set(skip_names) | set(_LOG_FILES)

    # Per-frame matrices carry their frame in the name, so they only need moving --
    # whatever directory they are in, and whether they are npz or text.
    for wd in worker_dirs:
        for p in sorted(wd.rglob('*')):
            if p.is_file() and _PER_FRAME.fullmatch(p.name):
                target = out_dir / p.relative_to(wd)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(p), str(target))

    # Signatures need the masterkey treatment before the generic text pass.
    rewrite_signatures(worker_dirs, out_dir)

    handled = {pathlib.PurePath('CNA/Signatures'), pathlib.PurePath('Exec/Masterkey')}
    remaining = set()
    for wd in worker_dirs:
        for p in wd.rglob('*'):
            if not p.is_file() or p.is_symlink():
                continue
            rel = p.relative_to(wd)
            if rel.name in skip or rel in handled or _PER_FRAME.fullmatch(rel.name):
                continue
            remaining.add(rel)
    for rel in sorted(remaining):
        sources = [wd / rel for wd in worker_dirs if (wd / rel).is_file()]
        if sources:
            _merge_frame_indexed(sources, out_dir / rel)

    for name in _LOG_FILES:
        chunks = [(wd / name).read_text() for wd in worker_dirs if (wd / name).is_file()]
        if chunks:
            with open(out_dir / name, 'a') as f:
                for i, text in enumerate(chunks):
                    f.write(f"\n--- worker {i} ---\n{text}")


def run_parallel(system, quantities, jobs, pattern_input=None, strict=False):
    """Run ``system``'s frames across ``jobs`` processes, then merge into its base_dir.

    Returns the number of chunks actually used (fewer than ``jobs`` for short runs).
    """
    base = pathlib.Path(system['base_dir'])
    chunks = plan_chunks(system['Start'], system['End'], system['Step'], jobs)
    if len(chunks) <= 1:
        from Sapphire import Process
        Process.Process(System=system, Quantities=quantities, Pattern_Input=pattern_input,
                        strict=strict, overwrite=True)
        return 1

    work_root = base / '.sapphire_workers'
    if work_root.exists():
        shutil.rmtree(work_root)
    work_root.mkdir(parents=True)

    payloads, worker_dirs = [], []
    for k, (lo, hi, step) in enumerate(chunks):
        wd = work_root / f'w{k}'
        wd.mkdir()
        # Process builds its input path as base_dir + movie_file_name, so the trajectory
        # has to be reachable from inside each worker directory. A symlink avoids copying
        # a multi-gigabyte file once per worker.
        src = base / system['movie_file_name']
        link = wd / system['movie_file_name']
        try:
            link.symlink_to(src.resolve())
        except OSError:                                  # filesystems without symlinks
            shutil.copy(src, link)
        sub = dict(system)
        sub['base_dir'] = str(wd) + os.sep
        sub['Start'], sub['End'], sub['Step'] = lo, hi, step
        payloads.append((sub, quantities, pattern_input, strict))
        worker_dirs.append(wd)

    with ProcessPoolExecutor(max_workers=len(chunks)) as pool:
        list(pool.map(_run_chunk, payloads))

    merge(worker_dirs, base, skip_names=(system['movie_file_name'],))
    shutil.rmtree(work_root, ignore_errors=True)
    return len(chunks)
