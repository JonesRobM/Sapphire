import numpy as np
import pytest

from Sapphire.Post_Process import DistFuncs

TETRAHEDRON = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float)


def test_distance():
    assert DistFuncs.distance(np.zeros(3), np.array([3.0, 4.0, 0.0])) == pytest.approx(5.0)


def test_get_com_is_centroid():
    np.testing.assert_allclose(DistFuncs.get_CoM(TETRAHEDRON), np.zeros(3))
    np.testing.assert_allclose(DistFuncs.get_CoM(TETRAHEDRON + 2.0), np.full(3, 2.0))


def test_euc_dist_regular_tetrahedron():
    d = np.asarray(DistFuncs.Euc_Dist(TETRAHEDRON))
    assert d.shape == (6,)
    np.testing.assert_allclose(d, np.sqrt(8.0))


def test_euc_dist_homo_subset():
    elements = ["Au", "Au", "Pt", "Pt"]
    d = np.asarray(DistFuncs.Euc_Dist(TETRAHEDRON, homo=True, specie="Au", elements=elements))
    assert d.shape == (1,)


def test_comdist_from_origin():
    d = np.asarray(DistFuncs.CoMDist(TETRAHEDRON, CoM=np.zeros(3)))
    np.testing.assert_allclose(d, np.sqrt(3.0))
