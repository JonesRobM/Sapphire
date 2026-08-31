import pathlib
import shutil

import pytest

REPO = pathlib.Path(__file__).resolve().parents[1]
AU561 = REPO / "main" / "Sapphire" / "Tutorials" / "Data" / "Au561.xyz"


@pytest.fixture
def au561_workdir(tmp_path):
    """A scratch directory containing the 561-atom Au icosahedron, laid out the way
    `Sapphire.Process` expects (`base_dir` with trailing separator + movie file)."""
    shutil.copy(AU561, tmp_path / "Au561.xyz")
    return tmp_path


AUPT_QUANTITIES = {
    "Full": {"pdf": None, "rdf": None, "com": None, "comdist": None, "adj": None, "nn": None, "agcn": None,
             "cna_sigs": None, "cna_patterns": None, "gyration": None, "collect": None, "concert": None},
    "Homo": {"hopdf": None, "hordf": None, "hocom": None, "hoadj": None, "hocomdist": None, "homobonds": None},
    "Hetero": {"hepdf": None, "herdf": None, "headj": None, "mix": None, "lae": None, "ele_nn": None,
               "heterobonds": None},
}


@pytest.fixture(scope="session")
def aupt_run(tmp_path_factory):
    """A 3-frame bimetallic run on the bundled AuPt sample (300 K -> ~1000 K). Shared by the
    bimetallic and graphing tests; ~10 s."""
    from Sapphire import Process
    from Sapphire.Tutorials import data
    work = tmp_path_factory.mktemp("aupt")
    data.sample("AuPt", work)
    System = {"base_dir": str(work) + "/", "movie_file_name": "AuPt_sample.xyz", "extend_xyz": None,
              "Homo": ["Au", "Pt"], "Hetero": True,
              "Start": 0, "End": 70, "Step": 34, "Skip": 1, "UniformPDF": False, "Band": 0.05}
    proc = Process.Process(System=System, Quantities=AUPT_QUANTITIES)
    proc.analyse({"JSD": ["pdf"]})
    return work, proc
