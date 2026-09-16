"""The hand-written kernels must agree with the textbook definitions.

``_gaussian_sum`` used to call ``scipy.stats.norm.pdf``; it now writes the Gaussian out,
which is ~2.3x faster because it skips argument validation and frozen-distribution
machinery a fixed scalar bandwidth does not need. scipy stays here as an independent
oracle, so the replacement is pinned to something that was not derived from it.
"""
import numpy as np
import pytest
from scipy.stats import norm

from Sapphire.Post_Process.Kernels import _CHUNK, _epanechnikov_sum, _gaussian_sum

# Reordering a large summation and factoring the normalisation out moves the last couple
# of bits. Measured worst case on real pair-distance data is ~8e-15 relative.
TOL = 1e-12


def reference_gaussian(data, space, band):
    """sum_i N(space; data_i, band), straight from scipy."""
    return norm.pdf(space[None, :], np.asarray(data)[:, None], band).sum(axis=0)


@pytest.mark.parametrize("n_data, n_space, band", [
    (1, 50, 0.05),
    (10, 100, 0.05),
    (1000, 200, 0.05),
    (1000, 200, 0.5),          # a wide kernel
    (1000, 200, 0.01),         # a narrow one
])
def test_gaussian_matches_scipy(n_data, n_space, band):
    rng = np.random.default_rng(0)
    data = rng.random(n_data) * 10
    space = np.linspace(0, 10, n_space)
    np.testing.assert_allclose(_gaussian_sum(data, space, band),
                               reference_gaussian(data, space, band), rtol=TOL, atol=0)


def test_gaussian_across_the_chunk_boundary():
    """The sum is accumulated in blocks; the seam must not drop or double-count anything."""
    rng = np.random.default_rng(1)
    space = np.linspace(0, 10, 64)
    for n in (_CHUNK - 1, _CHUNK, _CHUNK + 1, 2 * _CHUNK + 7):
        data = rng.random(n) * 10
        np.testing.assert_allclose(_gaussian_sum(data, space, 0.05),
                                   reference_gaussian(data, space, 0.05), rtol=TOL, atol=0)


def test_gaussian_is_normalised():
    """One data point integrates to one, so the prefactor is right."""
    space = np.linspace(-5, 5, 20001)
    density = _gaussian_sum(np.array([0.0]), space, 0.5)
    assert np.trapezoid(density, space) == pytest.approx(1.0, rel=1e-6)


def test_gaussian_peaks_at_the_datum():
    space = np.linspace(0, 10, 1001)
    density = _gaussian_sum(np.array([4.0]), space, 0.1)
    assert space[np.argmax(density)] == pytest.approx(4.0, abs=0.01)


def test_gaussian_scales_with_the_number_of_points():
    """It is a sum, not a mean: ten copies give ten times the density."""
    space = np.linspace(0, 10, 501)
    one = _gaussian_sum(np.array([5.0]), space, 0.2)
    ten = _gaussian_sum(np.full(10, 5.0), space, 0.2)
    np.testing.assert_allclose(ten, 10 * one, rtol=TOL, atol=0)


def test_empty_data_gives_zero_density():
    space = np.linspace(0, 10, 32)
    np.testing.assert_array_equal(_gaussian_sum(np.array([]), space, 0.05), np.zeros_like(space))


def test_epanechnikov_is_compactly_supported():
    """It must vanish beyond one bandwidth, unlike the Gaussian."""
    space = np.linspace(0, 10, 1001)
    density = _epanechnikov_sum(np.array([5.0]), space, 0.5)
    assert density[np.abs(space - 5.0) > 0.5].max() == 0.0
    assert density[np.abs(space - 5.0) < 0.4].min() > 0.0
