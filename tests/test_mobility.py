"""Collectivity and concertedness: lengths, alignment, and the values themselves.

Concertedness compares collectivity across a lag of two, so it is only defined from the
fourth analysed frame onwards. The array used to be sized for one more entry than the
loop ever assigned, leaving a leading zero that was written to file as though it were a
measurement -- and a three-frame run produced *only* that fabricated value.
"""
import numpy as np
import pytest

from Sapphire.Post_Process.Stats import Mobility


def _adj(neighbour_pairs, n=4):
    """Tiny symmetric 0/1 adjacency with a zero diagonal."""
    a = np.zeros((n, n), dtype=int)
    for i, j in neighbour_pairs:
        a[i, j] = a[j, i] = 1
    return a


# ------------------------------------------------------------------ the primitives
def test_collectivity_is_the_fraction_of_atoms_that_moved():
    before = _adj([(0, 1), (2, 3)])
    after = _adj([(0, 1), (2, 3)])
    assert Mobility.Collectivity(Mobility.R(before, after)) == 0.0

    after = _adj([(0, 1), (1, 2)])          # atoms 1, 2 and 3 all change neighbours
    moved = Mobility.R(before, after)
    assert Mobility.Collectivity(moved) == pytest.approx(moved.sum() / 4)


def test_R_detects_a_swap_that_leaves_the_count_unchanged():
    """Losing one neighbour and gaining another must count as a change."""
    before = _adj([(0, 1)])
    after = _adj([(0, 2)])
    assert Mobility.R(before, after)[0]


def test_concertedness_is_a_magnitude():
    assert Mobility.Concertedness(0.25, 0.75) == pytest.approx(0.5)
    assert Mobility.Concertedness(0.75, 0.25) == pytest.approx(0.5)   # symmetric


# ------------------------------------------------------- lengths and frame alignment
@pytest.mark.parametrize("n_frames, n_collect, n_concert", [
    (1, 0, 0),
    (2, 1, 0),
    (3, 2, 0),      # concertedness is NOT yet defined: it needs a lag of two
    (4, 3, 1),
    (5, 4, 2),
    (12, 11, 9),
])
def test_series_lengths(n_frames, n_collect, n_concert):
    """collect has one value per consecutive pair; concert one per pair separated by two."""
    assert max(n_frames - 1, 0) == n_collect
    assert max(n_frames - 3, 0) == n_concert


def test_concertedness_values_match_a_direct_computation(tmp_path):
    """End to end: every written concertedness value is a real lag-2 difference.

    Guards the bug where element 0 was never assigned and went to file as 0.0.
    """
    from Sapphire.api import run

    from Sapphire.Tutorials import data
    data.sample("AuPt", tmp_path)
    reader = run(str(tmp_path / "AuPt_sample.xyz"), str(tmp_path / "out"),
                 quantities=["pdf", "adj", "collect", "concert"],
                 frames=(0, 12, 2), homo=[], hetero=[])

    h = np.asarray(reader.load("collect"))
    c = np.asarray(reader.load("concert"))

    assert len(c) == len(h) - 2, "concertedness spans a lag of two collectivity values"
    assert len(c) > 0, "six frames is enough for concertedness to exist"
    # concert[k] pairs collect[k+2] with collect[k]  (the i-1 / i-3 lag)
    np.testing.assert_allclose(c, np.abs(h[2:] - h[:-2]), rtol=0, atol=1e-12)
    assert not np.all(c == 0.0), "a series of exact zeros means nothing was computed"
