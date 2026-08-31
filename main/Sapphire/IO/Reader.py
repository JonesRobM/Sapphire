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
_FRAME = re.compile(r"([A-Z][a-z]*)?File(\d+)")  # optional species prefix + per-frame index


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

    # ------------------------------------------------------------------ discovery
    def available(self) -> dict[str, pathlib.Path | list[pathlib.Path]]:
        """Map metadata key -> file (or list of per-frame files for matrix quantities).

        Species-suffixed homo files (``HomoCoMAu``) become ``hocomAu``; per-frame matrix
        sets (``Adjacency/File0``, ``Time_Dependent/HeAdjFile7``, ``Adjacency/HomoAdjPtFile0``)
        collapse to one key (``adj``, ``headj``, ``hoadjPt``) holding the ordered file list.
        """
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
            files = sorted(adj.glob("File*"), key=lambda p: int(p.name[4:]))
            if files:
                found["adj"] = files
        for v in found.values():
            if isinstance(v, list):
                v.sort(key=lambda p: int(p.name.split("File")[-1]))
        return found

    def frames(self, key: str) -> np.ndarray:
        path = self.available()[key]
        if isinstance(path, list):
            return np.array([int(p.name.split("File")[-1]) for p in path])
        return np.array([int(line.split(" ", 1)[0]) for line in path.read_text().splitlines() if line.strip()])

    # ------------------------------------------------------------------ loading
    def load(self, key: str):
        """Return the quantity as an array indexed ``[frame, ...]`` (a list for CNA patterns)."""
        if key == "masterkey":
            return self.masterkey()  # CNA signature labels, e.g. "421" -> keep as strings
        path = self.available().get(key)
        if path is None:
            raise KeyError(f"{key!r} not present in {self.base}; have {sorted(self.available())}")
        if isinstance(path, list):  # per-frame matrices
            return np.array([np.loadtxt(p, dtype=np.int32, ndmin=2) for p in path])
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
