"""IO.Reader reads a run directory back; Process.analyse / write_meta / strict work on it."""
import pickle
import shutil

import numpy as np
import pytest

from Sapphire import Process
from Sapphire.IO.Reader import Reader, parse_line
from tests.conftest import AU561

SYSTEM = {"movie_file_name": "Au561.xyz", "extend_xyz": None, "Homo": None, "Hetero": False,
          "Start": 0, "End": 1, "Step": 1, "Skip": 1, "UniformPDF": False, "Band": 0.05}
FULL = {"pdf": None, "rdf": None, "com": None, "comdist": None, "adj": None, "nn": None,
        "agcn": None, "gyration": None, "cna_sigs": None, "cna_patterns": None}


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    work = tmp_path_factory.mktemp("reader")
    shutil.copy(AU561, work / "Au561.xyz")
    proc = Process.Process(System={"base_dir": str(work) + "/", **SYSTEM}, Quantities={"Full": FULL})
    return work, proc


def test_parse_line_scalar_vector_pattern():
    f, v = parse_line("3 12 12 11")
    assert f == 3 and v.tolist() == [12, 12, 11]
    f, v = parse_line("0 [0. 0. 0.] [ 2.4  0. -1.5]")
    assert v.shape == (2, 3) and v[1, 2] == -1.5
    f, v = parse_line("0 ((12, (5, 5, 5)),) ((2, (5, 5, 5)), (10, (4, 2, 2)))")
    assert v[0] == ((12, (5, 5, 5)),) and v[1][1] == (10, (4, 2, 2))


def test_available_maps_files_to_keys(run):
    work, _ = run
    keys = Reader(work).available()
    assert {"nn", "pdf", "pdfspace", "com", "comdist", "gyration", "cna_sigs", "pattern_indices", "adj", "rcut"} <= set(keys)


def test_load_shapes(run):
    work, _ = run
    r = Reader(work)
    assert r.load("nn").shape == (1, 561)
    assert r.load("com").shape == (1, 3)          # one centre-of-mass vector per frame
    assert r.load("gyration").shape == (1,)
    assert r.load("adj").shape == (1, 561, 561)
    assert r.load("cna_sigs").shape == (1, len(r.masterkey()))
    pats = r.load("pattern_indices")
    assert len(pats) == 1 and len(pats[0]) == 561 and pats[0][0] == ((12, (5, 5, 5)),)


def test_adjacency_is_symmetric_and_matches_nn(run):
    work, _ = run
    r = Reader(work)
    adj = r.load("adj")[0]
    assert (adj == adj.T).all()
    np.testing.assert_array_equal(adj.sum(axis=1), r.load("nn")[0])


def test_missing_key_is_informative(run):
    work, _ = run
    with pytest.raises(KeyError, match="hopdf"):
        Reader(work).load("hopdf")


def test_analyse_statistics_and_write_meta(run):
    work, proc = run
    proc.analyse({"JSD": ["pdf"], "Kullback": ["rdf"]})
    assert "JSDpdf" in proc.metadata and "Kullbackrdf" in proc.metadata
    assert proc.metadata["JSDpdf"].shape == (1,)
    assert proc.metadata["JSDpdf"][0] == pytest.approx(0.0, abs=1e-9)   # frame 0 vs itself
    path = proc.write_meta()
    with open(path, "rb") as f:
        meta = pickle.load(f)
    assert meta["nn"].shape == (1, 561)
    log_file = work / "Sapphire_Errors.log"
    log = log_file.read_text() if log_file.exists() else ""
    assert "Exception raised" not in log and "Error raised" not in log, log


def test_strict_reraises(tmp_path, monkeypatch):
    shutil.copy(AU561, tmp_path / "Au561.xyz")
    from Sapphire.Post_Process import Kernels

    def boom(*a, **k):
        raise RuntimeError("forced failure")
    monkeypatch.setattr(Kernels, "Gauss", boom)
    with pytest.raises(RuntimeError, match="forced failure"):
        Process.Process(System={"base_dir": str(tmp_path) + "/", **SYSTEM},
                        Quantities={"Full": {"pdf": None}}, strict=True)
    # non-strict: same failure is logged, not raised
    shutil.copy(AU561, tmp_path / "Au561.xyz")
    Process.Process(System={"base_dir": str(tmp_path) + "/", **SYSTEM}, Quantities={"Full": {"pdf": None}})
    assert "forced failure" in (tmp_path / "Sapphire_Errors.log").read_text()


def test_extend_xyz_written(tmp_path):
    shutil.copy(AU561, tmp_path / "Au561.xyz")
    system = {"base_dir": str(tmp_path) + "/", **SYSTEM, "extend_xyz": ["nn", "agcn"]}
    proc = Process.Process(System=system, Quantities={"Full": {"pdf": None, "adj": None, "nn": None, "agcn": None}})
    proc.analyse({})
    out = tmp_path / "Extended.xyz"
    assert out.exists() and out.stat().st_size > 0
    header = out.read_text().splitlines()[:3]
    assert header[0].strip() == "561"
    assert len(header[2].split()) == 4 + 2   # symbol, x, y, z, nn, agcn
