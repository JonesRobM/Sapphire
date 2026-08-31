"""Chemical-ordering descriptors for bimetallic clusters.

All three calculators take the adjacency matrices produced by
:class:`Sapphire.Post_Process.Adjacent.Adjacency_Matrix` for one frame:

* ``HoAdj[X]`` — (N_X, N_X) homo adjacency of species X,
* ``HeAdj``    — (N_A, N_B) hetero adjacency between species A (rows) and B (columns),
* ``Adj``      — (N, N) full adjacency, ``Elements`` the matching symbol list.

Outputs go to ``<base_dir>/Time_Dependent/`` one line per frame, following the
``IO/OutputInfo*`` tables, and are read back by :class:`Sapphire.IO.Reader.Reader`.

Definitions
-----------
Mixing parameter (Baletto & Ferrando style)::

    mix = (N_AA + N_BB - N_AB) / (N_AA + N_BB + N_AB)

+1: fully segregated (no hetero bonds), -1: perfectly alternating.

LAE (local atomic environment): for each species, the histogram over 0..12 of how
many *hetero* neighbours each atom has. Two rows per frame (one per species).

Ele_NN: for every atom (full ordering), the number of neighbours of species X.
"""
import os

import numpy as np

from Sapphire.Utilities.log import get_logger

log = get_logger('Sapphire.Post_Process.AtomicEnvironment')


def _dense(m):
    return np.asarray(m.todense() if hasattr(m, "todense") else m)


def _append(system, attributes, frame, values, suffix=""):
    """Append ``frame v1 v2 ...`` to the file described by an OutputInfo entry."""
    directory = system['base_dir'] + attributes['Dir']
    os.makedirs(directory, exist_ok=True)
    path = directory + attributes['File'] + suffix
    with open(path, 'a') as f:
        f.write(str(frame) + ' ' + ' '.join(str(v) for v in np.atleast_1d(values)) + '\n')


class Mix:
    """Mixing parameter plus homo / hetero bond counts for one frame."""

    def __init__(self, System=None, Frame=None, HoAdj=None, HeAdj=None,
                 HomoBonds=False, HeteroBonds=False):
        self.System, self.Frame = System, Frame
        self.HoAdj = {k: _dense(v) for k, v in (HoAdj or {}).items()}
        self.HeAdj = _dense(HeAdj)
        self.HomoBonds, self.HeteroBonds = HomoBonds, HeteroBonds
        self.calculate()
        if System is not None:
            self.write()

    def calculate(self):
        self.HoBonds = {k: int(a.sum() // 2) for k, a in self.HoAdj.items()}  # symmetric -> /2
        self.HeBonds = int(self.HeAdj.sum())
        n_ho = sum(self.HoBonds.values())
        total = n_ho + self.HeBonds
        self.Mix_Param = (n_ho - self.HeBonds) / total if total else np.nan
        return self.Mix_Param

    def write(self):
        from Sapphire.IO import OutputInfoHetero as He
        _append(self.System, He.mix, self.Frame, self.Mix_Param)
        if self.HeteroBonds:
            _append(self.System, He.hetero_bonds, self.Frame, self.HeBonds)
        if self.HomoBonds:
            from Sapphire.IO import OutputInfoHomo as Ho
            for specie, n in self.HoBonds.items():
                _append(self.System, Ho.homo_bonds, self.Frame, n, suffix=specie)


class LAE:
    """Histogram (0..12) of hetero-neighbour counts, per species."""

    def __init__(self, System=None, Frame=None, HeAdj=None, Species=None, MaxCN=12):
        self.System, self.Frame, self.Species = System, Frame, list(Species)
        self.HeAdj = _dense(HeAdj)
        self.MaxCN = MaxCN
        self.calculate()
        if System is not None:
            self.write()

    def calculate(self):
        bins = np.arange(self.MaxCN + 2)
        self.Hist = {
            self.Species[0]: np.histogram(self.HeAdj.sum(axis=1), bins=bins)[0],  # A atoms' B-neighbours
            self.Species[1]: np.histogram(self.HeAdj.sum(axis=0), bins=bins)[0],  # B atoms' A-neighbours
        }
        return self.Hist

    def write(self):
        from Sapphire.IO import OutputInfoHetero as He
        for specie, h in self.Hist.items():
            _append(self.System, He.lae, self.Frame, h, suffix=specie)


class Ele_NN:
    """For every atom, how many neighbours of each species it has."""

    def __init__(self, System=None, Frame=None, Adj=None, Elements=None, Species=None):
        self.System, self.Frame = System, Frame
        self.Adj = _dense(Adj)
        self.Elements = np.asarray(Elements)
        self.Species = list(Species)
        self.calculate()
        if System is not None:
            self.write()

    def calculate(self):
        self.EleNN = {x: self.Adj[:, self.Elements == x].sum(axis=1).astype(int) for x in self.Species}
        return self.EleNN

    def write(self):
        from Sapphire.IO import OutputInfoHetero as He
        for specie, counts in self.EleNN.items():
            _append(self.System, He.ele_nn, self.Frame, counts, suffix=specie)
