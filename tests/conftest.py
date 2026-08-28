import pathlib
import shutil

import pytest

REPO = pathlib.Path(__file__).resolve().parents[1]
AU561 = REPO / "main" / "Sapphire" / "Tutorials" / "CNA" / "Au561.xyz"


@pytest.fixture
def au561_workdir(tmp_path):
    """A scratch directory containing the 561-atom Au icosahedron, laid out the way
    `Sapphire.Process` expects (`base_dir` with trailing separator + movie file)."""
    shutil.copy(AU561, tmp_path / "Au561.xyz")
    return tmp_path
