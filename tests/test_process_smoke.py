"""End-to-end: run `Sapphire.Process` on the bundled 561-atom Au icosahedron.

Since the 1.0 refactor, results are written to ``<base_dir>/Time_Dependent/<Quantity>``
(one line per frame, leading integer = frame index), so that is what we check.
"""
import numpy as np
import pytest

from Sapphire import Process

FULL = {
    "euc": None, "rdf": None, "pdf": None, "com": None, "comdist": None,
    "adj": None, "nn": None, "agcn": {"Write_Movie": False},
    "gyration": None, "stat_radius": None, "surf_area": None, "surf_atoms": None,
    "moi": None, "cna_sigs": None, "cna_patterns": None,
}


def _read(path):
    rows = [line.split() for line in path.read_text().splitlines() if line.strip()]
    return rows


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    import shutil
    from tests.conftest import AU561
    work = tmp_path_factory.mktemp("au561")
    shutil.copy(AU561, work / "Au561.xyz")
    System = {
        "base_dir": str(work) + "/", "movie_file_name": "Au561.xyz", "extend_xyz": None,
        "Homo": None, "Hetero": False,
        "Start": 0, "End": 1, "Step": 1, "Skip": 1, "UniformPDF": False, "Band": 0.05,
    }
    proc = Process.Process(System=System, Quantities={"Full": FULL})
    proc.analyse({})
    return work, proc


def test_species_and_size(run):
    _, proc = run
    assert proc.NAtoms == [561]
    assert proc.Species == ["Au"]


def test_no_swallowed_exceptions(run):
    """Sapphire logs exceptions instead of raising; the log must contain none."""
    work, _ = run
    log = (work / "Sapphire_Errors.log").read_text() if (work / "Sapphire_Errors.log").exists() else ""
    assert "Exception raised" not in log, log


@pytest.mark.parametrize("name", ["Time_Dependent/NN", "Time_Dependent/AGCN", "Time_Dependent/RCut",
                                  "Time_Dependent/PDF", "Time_Dependent/RDF", "Time_Dependent/CoM",
                                  "Time_Dependent/CoMDist", "Time_Dependent/Gyration",
                                  "Time_Dependent/Stat_Radius", "Time_Dependent/Surf_Area",
                                  "Time_Dependent/Surf_Atoms", "Time_Dependent/MoI",
                                  "CNA/Signatures", "CNA/Patterns", "Exec/Masterkey"])
def test_output_written(run, name):
    work, _ = run
    f = work / name
    assert f.exists() and f.stat().st_size > 0, f"{name} not written"


def test_coordination_of_perfect_icosahedron(run):
    work, _ = run
    row = _read(work / "Time_Dependent" / "NN")[0]
    nn = np.asarray(row[1:], dtype=float)
    assert nn.shape == (561,)
    assert nn.max() == 12                     # bulk atoms
    assert nn.min() >= 6                      # vertex atoms of an Ih
    assert (nn == 12).sum() == 309            # 5-shell Ih: 309 interior, 252 surface atoms


def test_first_shell_cutoff_is_physical(run):
    work, _ = run
    rcut = float(_read(work / "Time_Dependent" / "RCut")[0][1])
    assert 2.9 < rcut < 3.6                   # between Au 1st (2.88 Å) and 2nd (4.08 Å) shells


def test_gyration_radius(run):
    work, _ = run
    rg = float(_read(work / "Time_Dependent" / "Gyration")[0][1])
    assert 9.5 < rg < 11.5                    # Au561 Ih, R ≈ 13.4 Å  ->  Rg ≈ sqrt(3/5) R


def test_adjacency_matrix_saved(run):
    work, _ = run
    assert any(p.stat().st_size > 0 for p in (work / "Adjacency").iterdir())


def test_cna_signatures_of_icosahedron(run):
    """Au561 is a perfect 5-shell Mackay icosahedron: 12 five-fold axes x 10 (5,5,5) bonds = 120."""
    work, _ = run
    masterkey = (work / "Exec" / "Masterkey").read_text().split()
    sigs = np.asarray(_read(work / "CNA" / "Signatures")[0][1:], dtype=float)
    assert len(sigs) == len(masterkey)
    assert sigs[masterkey.index("555")] == 120
    assert sigs[masterkey.index("421")] > 0 and sigs[masterkey.index("422")] > 0


def test_agcn_is_generalised_coordination(run):
    """aGCN_i = sum over neighbours j of CN_j / 12 (Calle-Vallejo): terrace 7.5, Ih vertex 4.33, bulk 12."""
    work, _ = run
    agcn = np.asarray(_read(work / "Time_Dependent" / "AGCN")[0][1:], dtype=float)
    nn = np.asarray(_read(work / "Time_Dependent" / "NN")[0][1:], dtype=float)
    assert agcn.max() <= 12.0 + 1e-9
    interior = agcn[nn == 12]
    assert interior.max() == 12.0 and interior.min() >= 9.5         # deep bulk 12; sub-surface 11.25; under edges/vertices 9.83
    assert np.allclose(agcn[nn == 6], 52 / 12, atol=1e-3)          # 12 vertices: one CN-12 below, five CN-8 edges
    facet = agcn[nn == 9]
    # (111) facet atoms: 3 CN-12 below + 6 in-plane neighbours of CN 8/9 -> 7.167 or 7.333 on a 5-shell Ih
    # (a facet large enough for all-CN-9 neighbours would give the textbook terrace value 7.5)
    assert len(facet) and facet.min() >= 7.1 and facet.max() <= 7.5
