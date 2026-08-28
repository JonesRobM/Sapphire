"""Graphing must at least consume real Reader output without error."""
import numpy as np

from Sapphire.Graphing import Reader as GReader


def test_get_heights_from_signature_dicts():
    masterkey = [(4, 2, 1), (4, 2, 2), (5, 5, 5)]
    frames = [{(4, 2, 1): 6, (5, 5, 5): 1}, {(4, 2, 2): 3}]
    h = GReader.Get_Heights_Ovito(frames, masterkey, Norm=False)
    assert h.shape == (2, 3)
    np.testing.assert_array_equal(h[0], [6, 0, 1])
    np.testing.assert_array_equal(h[1], [0, 3, 0])
