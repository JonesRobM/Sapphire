"""Gupta (RGL / second-moment tight-binding) potential.

    E_i = Σ_j A_ij exp[−p_ij (r_ij/r0_ij − 1)]  −  sqrt( Σ_j ξ_ij² exp[−2 q_ij (r_ij/r0_ij − 1)] )

Parameters ``[p, q, A (eV), ξ (eV), r0 (Å)]`` per species / pair come from
:mod:`Sapphire.Potentials.GuptaParameters` (Cleri & Rosato 1993 and Baletto-group alloy fits).
This is a plain numpy implementation for analysis (energies, per-atom energies); for MD use an
ASE calculator (see :class:`GuptaCalculator`).
"""
from __future__ import annotations

import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from scipy.spatial.distance import cdist

from Sapphire.Potentials import GuptaParameters as GP


def parameter_set(species):
    """Return the parameter dict for one or two species, e.g. ('Au',) or ('Au', 'Pt')."""
    species = tuple(sorted(set(species)))
    if len(species) == 1:
        return getattr(GP, f"{species[0]}_parameters")
    for a, b in (species, species[::-1]):
        name = f"{a}{b}_parameters"
        if hasattr(GP, name):
            return getattr(GP, name)
    raise KeyError(f"no Gupta parameters for {species}")


def _pair_params(params, s1, s2):
    if s1 == s2:
        return params[s1]
    for key in ((s1, s2), (s2, s1)):
        if key in params and len(params[key]) == 5:
            return params[key]
    raise KeyError(f"no cross parameters for {s1}-{s2}")


def per_atom_energies(positions, symbols, params=None, r_cut=None):
    """Per-atom Gupta energies (eV) for a finite cluster; ``r_cut`` defaults to 3·max(r0)."""
    pos = np.asarray(positions, dtype=float)
    symbols = list(symbols)
    params = params or parameter_set(symbols)
    kinds = sorted(set(symbols))
    P = {(a, b): _pair_params(params, a, b) for a in kinds for b in kinds}
    r0max = max(v[4] for v in P.values())
    r_cut = r_cut or 3.0 * r0max
    d = cdist(pos, pos)
    np.fill_diagonal(d, np.inf)
    sym = np.asarray(symbols)
    rep = np.zeros(len(pos)); band = np.zeros(len(pos))
    for a in kinds:
        ia = sym == a
        for b in kinds:
            ib = sym == b
            p, q, A, xi, r0 = P[(a, b)]
            r = d[np.ix_(ia, ib)]
            mask = r < r_cut
            x = np.where(mask, r / r0 - 1.0, 0.0)
            rep[ia] += (A * np.exp(-p * x) * mask).sum(1)
            band[ia] += (xi**2 * np.exp(-2 * q * x) * mask).sum(1)
    return rep - np.sqrt(band)


def total_energy(positions, symbols, **kw):
    return float(per_atom_energies(positions, symbols, **kw).sum())


class GuptaCalculator(Calculator):
    """ASE calculator (energies + numerical forces) for quick relaxations of small clusters."""

    implemented_properties = ["energy", "energies", "forces"]

    def __init__(self, params=None, r_cut=None, **kwargs):
        super().__init__(**kwargs)
        self.params, self.r_cut = params, r_cut

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        e = per_atom_energies(atoms.positions, atoms.get_chemical_symbols(), self.params, self.r_cut)
        self.results["energies"] = e
        self.results["energy"] = float(e.sum())
        if "forces" in properties:
            self.results["forces"] = self.calculate_numerical_forces(atoms, d=1e-4)
