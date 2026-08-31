"""Structure classification from CNA patterns (successor to legacy/CNA/main.py + Model.py).

Pipeline (Armand's design, on Sapphire's own CNA rather than OVITO):

1. :func:`fingerprint` — per-atom CNA patterns of one structure via the PDDF-derived cutoff,
   `Adjacency_Matrix` and `FrameSignature.CNA`.
2. :func:`pattern_counts` — a structure's feature vector: how many atoms carry each pattern in a
   *masterkey* of patterns. The bulk patterns that separate the classical motifs are

   * ``((12,(4,2,1)),)``                 fcc bulk
   * ``((6,(4,2,1)),(6,(4,2,2)))``        hcp / twin plane
   * ``((10,(4,2,2)),(2,(5,5,5)))``       five-fold axis (decahedral / icosahedral spine)
   * ``((12,(5,5,5)),)``                  icosahedral centre
3. :class:`Classifier` — an SVC on those counts normalised per atom; classes 0 fcc, 1 Ih, 2 Dh,
   3 amorphous (same encoding as the original ``Model.py``). Save/load with joblib.

Label reference structures by filename (`Ih`, `Dh`, `To/Co/Oc`, `Am`) or explicitly.
"""
from __future__ import annotations

import os
import re

import numpy as np

from Sapphire.CNA.FrameSignature import CNA
from Sapphire.Post_Process import Adjacent, DistFuncs, Kernels

CLASSES = {0: "fcc", 1: "icosahedral", 2: "decahedral", 3: "amorphous"}
BULK_PATTERNS = (
    ((12, (4, 2, 1)),),
    ((6, (4, 2, 1)), (6, (4, 2, 2))),
    ((10, (4, 2, 2)), (2, (5, 5, 5))),
    ((12, (5, 5, 5)),),
)
_LABEL_RULES = ((r"Ih|Ico", 1), (r"Dh|Deca", 2), (r"Am", 3), (r"To|Co|Oc|fcc|FCC", 0))


def label_from_name(name):
    for pat, lab in _LABEL_RULES:
        if re.search(pat, os.path.basename(name)):
            return lab
    raise ValueError(f"cannot infer a structure class from {name!r}")


def _canon(pattern):
    """Order-independent form of a pattern tuple ((count, sig), ...)."""
    return tuple(sorted(pattern, key=lambda cs: (-cs[0], cs[1])))


def cutoff(positions, band=0.05):
    """First minimum of the Gaussian-KDE pair-distance distribution."""
    dist = DistFuncs.Euc_Dist(positions)
    return Kernels.Gauss(dist, band).ReturnRCut()


def fingerprint(atoms_or_positions, r_cut=None):
    """Per-atom CNA patterns for one structure. Returns a list of canonical pattern tuples."""
    pos = getattr(atoms_or_positions, "positions", atoms_or_positions)
    pos = np.asarray(pos, dtype=float)
    dist = DistFuncs.Euc_Dist(pos)
    if r_cut is None:
        r_cut = Kernels.Gauss(dist, 0.05).ReturnRCut()
    adj = Adjacent.Adjacency_Matrix(Positions=pos, Distances=dist, R_Cut=r_cut).ReturnAdj()
    cna = CNA(System=None, Adj=adj, Fingerprint=True)
    cna.calculate()
    return [_canon(p) for p in cna.Fingerprint]


def pattern_counts(patterns, masterkey):
    """Count atoms per masterkey pattern -> 1-D array aligned with masterkey."""
    idx = {_canon(p): i for i, p in enumerate(masterkey)}
    out = np.zeros(len(masterkey))
    for p in patterns:
        i = idx.get(_canon(p))
        if i is not None:
            out[i] += 1
    return out


def features(patterns, masterkey=BULK_PATTERNS):
    """Normalised feature vector: fraction of atoms carrying each masterkey pattern."""
    c = pattern_counts(patterns, masterkey)
    return c / max(len(patterns), 1)


class Classifier:
    def __init__(self, masterkey=BULK_PATTERNS, **svc_kwargs):
        from sklearn.svm import SVC
        self.masterkey = tuple(masterkey)
        self.clf = SVC(**{"kernel": "rbf", "C": 10.0, **svc_kwargs})
        self.fitted = False

    def fit(self, structures, labels, r_cut=None):
        X = np.array([features(fingerprint(s, r_cut), self.masterkey) for s in structures])
        self.clf.fit(X, np.asarray(labels))
        self.fitted = True
        return self

    def predict(self, structures, r_cut=None):
        X = np.array([features(fingerprint(s, r_cut), self.masterkey) for s in structures])
        return self.clf.predict(X)

    def predict_one(self, structure, r_cut=None):
        return CLASSES[int(self.predict([structure], r_cut)[0])]

    def save(self, path):
        import joblib
        joblib.dump({"masterkey": self.masterkey, "clf": self.clf}, path)

    @classmethod
    def load(cls, path):
        import joblib
        d = joblib.load(path)
        obj = cls(masterkey=d["masterkey"])
        obj.clf, obj.fitted = d["clf"], True
        return obj


def train_from_directory(xyz_dir, r_cut=None, **svc_kwargs):
    """Fit a classifier on every *.xyz in a directory, labelling by filename."""
    from ase.io import read
    files = sorted(f for f in os.listdir(xyz_dir) if f.endswith(".xyz"))
    structures = [read(os.path.join(xyz_dir, f)) for f in files]
    labels = [label_from_name(f) for f in files]
    return Classifier(**svc_kwargs).fit(structures, labels, r_cut), dict(zip(files, labels))
