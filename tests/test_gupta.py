"""Gupta potential against the numbers it was fitted to (Cleri & Rosato 1993 cohesive energies)."""
import numpy as np
import pytest
from ase.cluster import Octahedron
from ase.build import bulk

from Sapphire.Potentials import GuptaPotential as G

COHESIVE = {"Au": -3.78, "Ag": -2.96, "Cu": -3.54, "Ni": -4.44, "Pt": -5.86, "Pd": -3.94}   # eV/atom, experiment


@pytest.mark.parametrize("element", sorted(COHESIVE))
def test_bulk_cohesive_energy(element):
    """Energy of an interior atom of a large fcc chunk must be the bulk cohesive energy."""
    a = bulk(element, "fcc").cell.lengths()[0] * np.sqrt(2)     # conventional lattice constant
    chunk = Octahedron(element, length=11, cutoff=3, latticeconstant=a)
    e = G.per_atom_energies(chunk.positions, chunk.get_chemical_symbols())
    centre = np.argmin(np.linalg.norm(chunk.positions - chunk.positions.mean(0), axis=1))
    assert e[centre] == pytest.approx(COHESIVE[element], rel=0.06)


def test_alloy_parameters_and_calculator():
    from ase.cluster import Icosahedron
    c = Icosahedron("Au", noshells=3, latticeconstant=4.08)
    syms = ["Pt" if i % 5 == 0 else "Au" for i in range(len(c))]
    c.set_chemical_symbols(syms)
    params = G.parameter_set(("Au", "Pt"))
    assert ("Au", "Pt") in params or ("Pt", "Au") in params
    c.calc = G.GuptaCalculator()
    e = c.get_potential_energy()
    assert -6 * len(c) < e < -2 * len(c)                        # a few eV per atom, bound
    assert c.get_forces().shape == (len(c), 3)
