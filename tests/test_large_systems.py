"""The large-system routes (KD-tree pairs, sparse adjacency) must give the same answers as the
all-pairs code. Lowering ``LARGE_N`` sends a small random cluster down both routes."""
import numpy as np
import pytest
from scipy.spatial.distance import cdist, pdist

from Sapphire.Post_Process import Adjacent, AtomicEnvironment, DistFuncs

RNG = np.random.default_rng(7)
POS = RNG.uniform(-8, 8, size=(400, 3))
ELEMENTS = np.where(RNG.random(400) < 0.6, "Au", "Pd").tolist()


@pytest.fixture
def large(monkeypatch):
    monkeypatch.setattr(DistFuncs, "LARGE_N", 10)


def test_pairs_within_matches_pdist():
    i, j, d = DistFuncs.pairs_within(POS, 3.0)
    full = pdist(POS)
    np.testing.assert_allclose(np.sort(d), np.sort(full[full < 3.0]))
    assert np.all(i < j)


def test_pairs_within_cross_matches_cdist():
    a, b = POS[:150], POS[150:]
    _, _, d = DistFuncs.pairs_within(a, 3.0, b)
    full = cdist(a, b).ravel()
    np.testing.assert_allclose(np.sort(d), np.sort(full[full < 3.0]))


def test_max_pair_distance_via_hull(large):
    assert DistFuncs.max_pair_distance(POS) == pytest.approx(pdist(POS).max(), rel=1e-12)


def test_euc_dist_r_max_keeps_short_range(large):
    full = pdist(POS)
    np.testing.assert_allclose(np.sort(DistFuncs.Euc_Dist(POS, r_max=6.2)), np.sort(full[full < 6.2]))


@pytest.mark.parametrize("kind, metals", [("Full", ["Au", "Pd"]), ("Homo", ["Au"]), ("Hetero", ["Au", "Pd"])])
def test_adjacency_same_both_routes(monkeypatch, kind, metals):
    def adj():
        return Adjacent.Adjacency_Matrix(Positions=POS, Distances=None, R_Cut=2.9, Type=kind,
                                         Metals=metals, Elements=ELEMENTS).ReturnAdj()
    small = adj()
    monkeypatch.setattr(DistFuncs, "LARGE_N", 10)
    big = adj()
    assert big.shape == small.shape
    assert (big != small).nnz == 0


def test_mix_accepts_sparse():
    ho = {x: Adjacent.Adjacency_Matrix(Positions=POS, Distances=None, R_Cut=2.9, Type="Homo",
                                       Metals=[x], Elements=ELEMENTS).ReturnAdj() for x in ("Au", "Pd")}
    he = Adjacent.Adjacency_Matrix(Positions=POS, Distances=None, R_Cut=2.9, Type="Hetero",
                                   Metals=["Au", "Pd"], Elements=ELEMENTS).ReturnAdj()
    sparse = AtomicEnvironment.Mix(HoAdj=ho, HeAdj=he).Mix_Param
    dense = AtomicEnvironment.Mix(HoAdj={k: v.toarray() for k, v in ho.items()}, HeAdj=he.toarray()).Mix_Param
    assert sparse == pytest.approx(dense)
