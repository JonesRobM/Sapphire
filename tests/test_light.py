"""Tabulated dielectric functions: shapes and the physics everyone knows (noble-metal plasmonics)."""
import numpy as np
import pytest

from Sapphire.Light import Epsilon_DFT as D


@pytest.mark.parametrize("metal", ["Ag", "Au", "Cu", "Ni", "Pt", "Pd"])
def test_tables_consistent(metal):
    c = getattr(D, metal)()
    assert len(c.wl) == len(c.n_real) == len(c.n_imag) >= 49
    assert np.all(np.diff(c.wl) > 0)                      # wavelengths ascending (nm)
    eps = c.epsilon(500.0)
    assert np.isfinite(eps) and eps.imag > 0             # absorbing


def test_noble_metals_are_plasmonic_in_the_visible():
    for metal in ("Ag", "Au", "Cu"):
        eps = getattr(D, metal)().epsilon(600.0)
        assert eps.real < -5                              # Drude-like negative permittivity


def test_silver_interband_edge():
    ag = D.Ag()
    assert ag.epsilon(320.0).real > ag.epsilon(500.0).real   # ε' rises towards the ~3.9 eV interband edge
