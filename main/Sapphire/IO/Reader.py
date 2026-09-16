"""Read Sapphire run output back into memory.

Since the 1.0 refactor every calculator appends one line per frame to
``<base_dir>/<Dir>/<File>`` (see ``IO/OutputInfo*.py`` for the table). The line format is::

    <frame> <token> <token> ...

where a token is a scalar (``12``, ``3.19``), a bracketed vector (``[x y z]``), a CNA
pattern tuple (``((2, (5, 5, 5)), (10, (4, 2, 2)))``) or a bare word (masterkey entries).
Adjacency matrices are the exception: one dense ``N x N`` file per frame,
``Adjacency/File<frame>``.

:class:`Reader` maps those files back to the metadata keys the rest of Sapphire uses
(``'nn'``, ``'pdf'``, ``'cna_sigs'``, ``'hocomAu'``, ...) so that ``Process.analyse``,
``write_meta`` and the extended-xyz writer can run on the files rather than on an
in-memory dictionary. It is deliberately tolerant: a quantity that was not requested
simply is not present in :meth:`available`.

The same directory layout is the natural drop-in point for quantities produced by other
codes — write a ``<frame> <values...>`` file under ``Time_Dependent/`` and it is readable.
"""
from __future__ import annotations

import ast
import os
import pathlib
import re
from inspect import getmembers

import numpy as np

from Sapphire.IO import OutputInfoExec, OutputInfoFull, OutputInfoHetero, OutputInfoHomo

_VEC = re.compile(r"\[([^\]]*)\]")
# optional species prefix + per-frame index, either dense text (File7) or sparse (File7.npz)
_FRAME = re.compile(r"([A-Z][a-z]*)?File(\d+)(?:\.npz)?")


def _frame_index(path) -> int:
    """Frame number from a per-frame matrix filename, whichever format it is in."""
    stem = path.name[:-4] if path.name.endswith(".npz") else path.name
    return int(stem.split("File")[-1])


def _load_matrix(path) -> np.ndarray:
    """One frame's matrix as a dense array, from sparse npz or dense text.

    Adjacency is written sparse by default; callers downstream (tests, tutorials,
    Plot_Funcs) expect dense arrays indexed [frame, i, j], so densify on read.
    """
    if path.name.endswith(".npz"):
        import scipy.sparse as spa
        return spa.load_npz(path).toarray().astype(np.int32)
    return np.loadtxt(path, dtype=np.int32, ndmin=2)


def _table() -> dict[str, tuple[str, str]]:
    """key -> (Dir, File) for every quantity Sapphire knows how to write."""
    out: dict[str, tuple[str, str]] = {}
    for mod in (OutputInfoFull, OutputInfoHomo, OutputInfoHetero, OutputInfoExec):
        for name, val in getmembers(mod):
            if isinstance(val, dict) and "Dir" in val and "File" in val and not val.get("Exec"):
                out[name] = (val["Dir"], val["File"])
    return out


def _split_top_level(s: str) -> list[str]:
    """Split a whitespace-separated sequence of parenthesised tuple reprs at depth 0."""
    parts, depth, cur = [], 0, []
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch.isspace() and depth == 0:
            if cur:
                parts.append("".join(cur)); cur = []
            continue
        cur.append(ch)
    if cur:
        parts.append("".join(cur))
    return parts


def parse_line(line: str):
    """Return ``(frame, payload)`` for one output line."""
    head, _, rest = line.strip().partition(" ")
    frame = int(head)
    rest = rest.strip()
    if not rest:
        return frame, np.array([])
    if rest.startswith("["):
        return frame, np.array([np.fromstring(m, sep=" ") for m in _VEC.findall(rest)])
    if rest.startswith("("):
        return frame, [ast.literal_eval(p) for p in _split_top_level(rest)]
    try:
        return frame, np.fromstring(rest, sep=" ")
    except ValueError:  # pragma: no cover - defensive
        return frame, rest.split()


class Reader:
    """Load a Sapphire run directory.

    >>> r = Reader("run/")
    >>> r.available()            # {'nn': PosixPath('run/Time_Dependent/NN'), ...}
    >>> nn = r.load("nn")        # shape (frames, atoms)
    >>> meta = r.load_all()      # dict ready for Process.analyse / ExtendXYZ
    """

    def __init__(self, base_dir: str | os.PathLike):
        self.base = pathlib.Path(base_dir)
        self._table = _table()
        self._available = None

    def refresh(self):
        """Forget the cached directory listing; call after a run adds files."""
        self._available = None
        return self

    # ------------------------------------------------------------------ discovery
    def available(self) -> dict[str, pathlib.Path | list[pathlib.Path]]:
        """Map metadata key -> file (or list of per-frame files for matrix quantities).

        Species-suffixed homo files (``HomoCoMAu``) become ``hocomAu``; per-frame matrix
        sets (``Adjacency/File0``, ``Time_Dependent/HeAdjFile7``, ``Adjacency/HomoAdjPtFile0``)
        collapse to one key (``adj``, ``headj``, ``hoadjPt``) holding the ordered file list.
        """
        # A run directory holds one matrix file per frame, so this walk is O(frames):
        # ~6 s at 20 000 frames. load() and frames() each call it, and load_all() calls it
        # once per key, so the scan is cached and invalidated explicitly via refresh().
        if self._available is not None:
            return self._available
        found: dict = {}
        claimed: set = set()
        # Longest File names first so 'HomoCoMDist' is not claimed by 'HomoCoM'.
        for key, (d, f) in sorted(self._table.items(), key=lambda kv: -len(kv[1][1])):
            directory = self.base / d
            if not directory.is_dir():
                continue
            for p in sorted(directory.iterdir()):
                if not p.is_file() or p in claimed or not p.name.startswith(f):
                    continue
                suffix = p.name[len(f):]
                m = _FRAME.fullmatch(suffix)
                if m:  # per-frame matrix file
                    k = key + (m.group(1) or "")
                    found.setdefault(k, []).append(p); claimed.add(p)
                elif p.stat().st_size > 0 and (suffix == "" or d.startswith("Time_Dependent")):
                    found[key + suffix] = p; claimed.add(p)
        stats = self.base / "Time_Dependent" / "Stats"
        if stats.is_dir():  # analyse() output: one series per file, key = file name (e.g. JSDpdf)
            for p in sorted(stats.iterdir()):
                if p.is_file() and p.stat().st_size > 0 and p not in claimed:
                    found[p.name] = p
        adj = self.base / "Adjacency"
        if adj.is_dir():
            files = sorted(adj.glob("File*"), key=_frame_index)
            if files:
                found["adj"] = files
        for v in found.values():
            if isinstance(v, list):
                v.sort(key=_frame_index)
        self._available = found
        return found

    def frames(self, key: str) -> np.ndarray:
        path = self.available()[key]
        if isinstance(path, list):
            return np.array([_frame_index(p) for p in path])
        return np.array([int(line.split(" ", 1)[0]) for line in path.read_text().splitlines() if line.strip()])

    # ------------------------------------------------------------------ loading
    def load(self, key: str):
        """Return the quantity as an array indexed ``[frame, ...]`` (a list for CNA patterns)."""
        if key == "masterkey":
            return self.masterkey()  # CNA signature labels, e.g. "421" -> keep as strings
        available = self.available()
        path = available.get(key)
        if path is None:
            raise KeyError(f"{key!r} not present in {self.base}; have {sorted(available)}")
        if isinstance(path, list):  # per-frame matrices
            return np.array([_load_matrix(p) for p in path])
        payloads = [parse_line(l)[1] for l in path.read_text().splitlines() if l.strip()]
        if payloads and isinstance(payloads[0], list):
            return payloads
        try:
            arr = np.array(payloads)
            return arr[:, 0] if arr.ndim == 2 and arr.shape[1] == 1 else arr
        except ValueError:  # ragged (e.g. NAtoms changes between frames)
            return payloads

    def load_all(self) -> dict:
        return {k: self.load(k) for k in self.available()}

    def masterkey(self) -> list[str]:
        p = self.base / "Exec" / "Masterkey"
        return p.read_text().split() if p.exists() else []
