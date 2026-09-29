"""Divergences between per-frame distributions (Post_Process.Stats.Dist_Stats).

The JSD used to be computed with P and Q swapped inside the logarithms, on unnormalised
densities: it returned sqrt(Jeffreys/2 - JSD), which is unbounded (a real 1415-atom melting
frame gave 1.04 against a true Jensen-Shannon distance of 0.09). These pin the definitions.
"""
import numpy as np
import pytest
from scipy.spatial.distance import jensenshannon
from scipy.stats import entropy

from Sapphire.Post_Process.Stats import Dist_Stats


@pytest.fixture
def pair():
    rng = np.random.default_rng(7)
    return rng.random(200), rng.random(200)


def test_jsd_matches_scipy_base2(pair):
    P, Q = pair
    assert Dist_Stats.JSD(P, Q) == pytest.approx(jensenshannon(P, Q, base=2), rel=1e-12)


def test_jsd_bounds():
    P = np.array([1.0, 1.0, 0.0, 0.0])
    Q = np.array([0.0, 0.0, 1.0, 1.0])
    assert Dist_Stats.JSD(P, P) == pytest.approx(0.0, abs=1e-12)
    assert Dist_Stats.JSD(P, Q) == pytest.approx(1.0)          # disjoint supports


def test_jsd_is_symmetric(pair):
    P, Q = pair
    assert Dist_Stats.JSD(P, Q) == pytest.approx(Dist_Stats.JSD(Q, P))


def test_jsd_does_not_depend_on_density_scaling(pair):
    # a density on a grid (sums to 1/dr) and the same shape as a histogram give one answer
    P, Q = pair
    assert Dist_Stats.JSD(33.3 * P, 33.3 * Q) == pytest.approx(Dist_Stats.JSD(P, Q))


def test_jsd_is_not_the_old_swapped_formula(pair):
    P, Q = pair
    eps = 1e-6
    old = np.sqrt(-0.5 * np.sum((P + eps) * np.log(2 * (Q + eps) / (P + Q + 2 * eps))
                                + (Q + eps) * np.log(2 * (P + eps) / (P + Q + 2 * eps))))
    assert Dist_Stats.JSD(P, Q) != pytest.approx(old, rel=1e-3)


def test_kl_matches_scipy_on_normalised_inputs(pair):
    P, Q = pair
    assert Dist_Stats.Kullback(P, Q) == pytest.approx(entropy(P, Q), rel=1e-12)   # entropy normalises
    assert Dist_Stats.Kullback(10 * P, 10 * Q) == pytest.approx(Dist_Stats.Kullback(P, Q))


def test_empty_distribution_is_an_error():
    with pytest.raises(ValueError):
        Dist_Stats.JSD(np.zeros(5), np.ones(5))
