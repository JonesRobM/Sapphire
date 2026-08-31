"""Homo / Hetero paths on the bundled AuPt sample (2 frames: 300 K and ~1000 K)."""
import numpy as np
import pytest



@pytest.fixture(scope="module")
def run(aupt_run):
    return aupt_run


def test_no_swallowed_exceptions(run):
    work, _ = run
    log_file = work / "Sapphire_Errors.log"
    log = log_file.read_text() if log_file.exists() else ""     # absent == nothing was ever reported
    assert "Exception raised" not in log and "Error raised" not in log, log


def test_species_and_shapes(run):
    _, proc = run
    m = proc.metadata
    assert proc.Species == ["Au", "Pt"] and proc.NAtoms[0] == 1415
    assert m["adj"].shape == (3, 1415, 1415)
    assert m["hoadjAu"].shape == (3, 1149, 1149) and m["hoadjPt"].shape == (3, 266, 266)
    assert m["headj"].shape == (3, 1149, 266)
    assert m["comdist"].shape == (3, 1415) and m["hocomdistAu"].shape == (3, 1149)


def test_bond_bookkeeping_is_consistent(run):
    """Full adjacency must equal homo blocks + hetero block, bond for bond."""
    _, proc = run
    m = proc.metadata
    for f in range(3):
        n_full = m["adj"][f].sum() // 2
        n_homo = m["hoadjAu"][f].sum() // 2 + m["hoadjPt"][f].sum() // 2
        n_he = m["headj"][f].sum()
        assert n_full == n_homo + n_he
        assert m["hetero_bonds"][f] == n_he
        assert m["homo_bondsAu"][f] + m["homo_bondsPt"][f] == n_homo
        expected_mix = (n_homo - n_he) / (n_homo + n_he)
        assert m["mix"][f] == pytest.approx(expected_mix)
        assert -1 <= m["mix"][f] <= 1


def test_lae_and_ele_nn(run):
    _, proc = run
    m = proc.metadata
    assert m["laeAu"].shape == (3, 13) and m["laeAu"][0].sum() == 1149
    assert m["laePt"][0].sum() == 266
    assert m["ele_nnAu"].shape == (3, 1415)
    np.testing.assert_array_equal(m["ele_nnAu"][0] + m["ele_nnPt"][0], m["nn"][0])


def test_com_distances_are_distances(run):
    _, proc = run
    d = proc.metadata["comdist"][0]
    assert d.ndim == 1 and (d >= 0).all() and 5 < d.max() < 25   # ~1400-atom cluster radius


def test_melting_increases_divergence(run):
    _, proc = run
    jsd = proc.metadata["JSDpdf"]
    assert jsd[0] == pytest.approx(0.0, abs=1e-9) and jsd[-1] > 0.5 and jsd[1] > 0


def test_collectivity_and_statistics_keys(run):
    _, proc = run
    m = proc.metadata
    # frames 0 -> 34 -> 68 of a melting run: neighbourhoods must change
    assert "collect" in m and m["collect"][0] > 0.0 and m["collect"][1] > 0.0
    jsd_keys = {k for k in m if k.startswith("JSD")}
    assert "JSDpdf" in jsd_keys and "JSDhopdfAu" in jsd_keys
    assert not any("space" in k for k in jsd_keys)   # never the *space grids
