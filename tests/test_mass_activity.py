import numpy as np
import pytest

from Sapphire.Post_Process import Mass_Activity as MA


def test_volcano_shape():
    a = MA.site_activity([6.5, 7.5, MA.GCN_APEX, 9.0, 11.0])
    assert a[1] == pytest.approx(np.exp(3.14 * 7.5 - 23.40))       # Pt(111)-like site ~ 1.16
    assert a[2] == a.max()                                          # apex where the branches meet
    assert a[3] < a[2] and a[4] < a[3]                              # right branch falls
    assert MA.site_activity([5.9, 6.0])[1] == 0.0                   # threshold: GCN must exceed 6


def test_volcano_is_continuous_at_switch():
    a = MA.branch_intersection()
    assert 8.0 < a < 8.2                                            # printed coefficients meet at ~8.10, not 8.33
    left = MA.site_activity([a - 1e-9]); right = MA.site_activity([a + 1e-9])
    assert left[0] == pytest.approx(right[0], rel=1e-6)


def test_prefactor_reproduces_paper_constant():
    """Eq. 7: 4.107 A/mg for pure Pt."""
    n = 561
    pref = MA.prefactor_A_per_mg(MA.masses_from_symbols(["Pt"] * n))
    assert pref * n == pytest.approx(4.107, rel=2e-3)


def test_alloy_mass_is_species_weighted():
    pure = MA.prefactor_A_per_mg(MA.masses_from_symbols(["Pt"] * 10))
    mixed = MA.prefactor_A_per_mg(MA.masses_from_symbols(["Pt"] * 5 + ["Ni"] * 5))
    assert mixed > pure                                             # Ni is lighter -> more A per mg


def test_mass_activity_series_shape_and_threshold():
    agcn = np.array([[12.0] * 8 + [7.5] * 2, [12.0] * 10])
    ma = MA.mass_activity_series(agcn, ["Pt"] * 10)
    assert ma.shape == (2,) and ma[0] > 0 and ma[1] < 1e-6         # all-bulk frame: right branch ~0
