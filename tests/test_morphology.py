"""Morphology on the 5-shell Mackay icosahedron Au561 (shells of 10k²+2 atoms + centre)."""
import numpy as np
import pytest
from ase.io import read

from Sapphire.Post_Process import Morphology as M
from tests.conftest import AU561

R_CUT = 3.2      # between Au 1st (2.88 Å) and 2nd (4.08 Å) shells
R_AU = 1.44


@pytest.fixture(scope="module")
def au561():
    return read(AU561).positions


def test_coordination_and_gcn(au561):
    a = M.adjacency(au561, R_CUT)
    cn = M.coordination(a)
    assert cn.max() == 12 and (cn == 12).sum() == 309
    gcn = M.generalised_coordination(a)
    assert gcn.max() == pytest.approx(12.0) and gcn.min() > 0


def test_peeling_reproduces_mackay_shells(au561):
    pops = [len(s) for s in M.peel(au561, R_CUT)]
    assert pops == [252, 162, 92, 42, 12, 1]


def test_faceting_ratio_and_area(au561):
    a = M.adjacency(au561, R_CUT)
    assert M.faceting_ratio(a) == pytest.approx(100 * 252 / 561)
    area = M.surface_area(a, R_AU)
    sphere = 4 * np.pi * (M.radii(au561).max() + R_AU) ** 2
    assert 0.3 * sphere < area < 1.5 * sphere


def test_thickness_and_volume(au561):
    lo, hi = M.shell_thickness(au561, R_CUT)
    assert 2.0 < lo < 3.5 and 2.0 < hi < 3.5          # one Mackay shell ≈ 2.5–2.9 Å radially
    a = M.adjacency(au561, R_CUT)
    v = M.volume(au561, a, R_AU)
    v_sphere = 4 / 3 * np.pi * (M.radii(au561).max() + R_AU) ** 3
    assert 0.4 * v_sphere < v < 1.2 * v_sphere


def test_describe_keys(au561):
    d = M.describe(au561, R_CUT, R_AU)
    assert d["n_atoms"] == 561 and d["shell_populations"][0] == 252
