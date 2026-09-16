"""Adjacency is stored sparse by default; dense text stays available on request.

The matrix has ~12 nonzeros per row, so dense text spends O(N^2) bytes on O(N)
information -- 8 MB per frame at N=2000, 160 GB over a 20k-frame run. What must not
change is the data: whichever format is on disk, the Reader has to hand back the same
dense arrays, and `sapphire expand` has to reproduce the original text byte for byte.
"""
import numpy as np
import pytest
import scipy.sparse as spa

from Sapphire.IO.Reader import _frame_index, _load_matrix
from Sapphire.Post_Process.Adjacent import write_adjacency
from Sapphire.cli import main


@pytest.fixture
def matrix():
    rng = np.random.default_rng(0)
    a = (rng.random((40, 40)) < 0.2).astype(int)
    a = np.triu(a, 1)
    a = a + a.T                       # symmetric, zero diagonal, as Adjacent builds it
    return a


# ----------------------------------------------------------------- the two formats
def test_sparse_is_the_default_and_round_trips(tmp_path, matrix):
    written = write_adjacency(str(tmp_path / "File0"), spa.csr_matrix(matrix))
    assert written.endswith(".npz")
    np.testing.assert_array_equal(_load_matrix(tmp_path / "File0.npz"), matrix)


def test_dense_text_on_request_round_trips(tmp_path, matrix):
    written = write_adjacency(str(tmp_path / "File0"), spa.csr_matrix(matrix), dense_text=True)
    assert not written.endswith(".npz")
    np.testing.assert_array_equal(_load_matrix(tmp_path / "File0"), matrix)


def test_both_formats_agree(tmp_path, matrix):
    write_adjacency(str(tmp_path / "File0"), spa.csr_matrix(matrix))
    write_adjacency(str(tmp_path / "File1"), spa.csr_matrix(matrix), dense_text=True)
    np.testing.assert_array_equal(_load_matrix(tmp_path / "File0.npz"),
                                  _load_matrix(tmp_path / "File1"))


def test_sparse_is_dramatically_smaller(tmp_path, matrix):
    write_adjacency(str(tmp_path / "File0"), spa.csr_matrix(matrix))
    write_adjacency(str(tmp_path / "File1"), spa.csr_matrix(matrix), dense_text=True)
    assert (tmp_path / "File0.npz").stat().st_size < (tmp_path / "File1").stat().st_size


# ---------------------------------------------------------------- filename parsing
@pytest.mark.parametrize("name, expected", [
    ("File0", 0), ("File7", 7), ("File12", 12),
    ("File0.npz", 0), ("File7.npz", 7), ("File12.npz", 12),
    ("HomoAdjPtFile3", 3), ("HomoAdjPtFile3.npz", 3),
])
def test_frame_index_handles_both_formats(tmp_path, name, expected):
    assert _frame_index(tmp_path / name) == expected


# ------------------------------------------------------------------ sapphire expand
def test_expand_reproduces_the_dense_text_exactly(tmp_path, matrix, capsys):
    adj = tmp_path / "Adjacency"
    adj.mkdir()
    for i in range(3):
        write_adjacency(str(adj / f"File{i}"), spa.csr_matrix(matrix))
        write_adjacency(str(tmp_path / f"reference{i}"), spa.csr_matrix(matrix), dense_text=True)

    out = tmp_path / "dense"
    assert main(["expand", str(tmp_path), "-o", str(out)]) == 0
    # -o preserves the run's layout, so Adjacency/ and Time_Dependent/ matrices
    # (the hetero ones live in the latter) stay distinguishable.
    for i in range(3):
        assert (out / "Adjacency" / f"File{i}").read_bytes() == (tmp_path / f"reference{i}").read_bytes()


def test_expand_respects_a_frame_range(tmp_path, matrix):
    adj = tmp_path / "Adjacency"
    adj.mkdir()
    for i in range(6):
        write_adjacency(str(adj / f"File{i}"), spa.csr_matrix(matrix))
    out = tmp_path / "dense"
    assert main(["expand", str(tmp_path), "-o", str(out), "--frames", "2:6:2"]) == 0
    assert sorted(p.name for p in (out / "Adjacency").iterdir()) == ["File2", "File4"]


def test_expand_without_a_run_directory_exits_2(tmp_path, capsys):
    assert main(["expand", str(tmp_path)]) == 2
    assert "no adjacency matrices" in capsys.readouterr().err


def test_expand_says_so_when_already_dense(tmp_path, matrix, capsys):
    adj = tmp_path / "Adjacency"
    adj.mkdir()
    write_adjacency(str(adj / "File0"), spa.csr_matrix(matrix), dense_text=True)
    assert main(["expand", str(tmp_path)]) == 0
    assert "already holds dense text" in capsys.readouterr().err


# ------------------------------------------------------------------------ the config
def test_bad_adj_format_is_rejected(tmp_path):
    from Sapphire.api import Config
    traj = tmp_path / "t.xyz"
    traj.write_text("1\n\nAu 0.0 0.0 0.0\n")
    with pytest.raises(ValueError, match="adj_format"):
        Config(trajectory=str(traj), adj_format="parquet").validate()


def test_adj_format_reaches_the_legacy_system_dict():
    from Sapphire.api import Config
    for fmt in ("npz", "text"):
        system, _ = Config(trajectory="x", adj_format=fmt).to_legacy(n_frames=3, species=["Au"])
        assert system["adj_format"] == fmt


def test_expand_finds_hetero_matrices_outside_the_adjacency_directory(tmp_path, matrix):
    """Full and homo matrices live under Adjacency/, hetero ones under Time_Dependent/.

    A scan restricted to Adjacency/ silently skips a third of the matrices.
    """
    (tmp_path / "Adjacency").mkdir()
    (tmp_path / "Time_Dependent").mkdir()
    write_adjacency(str(tmp_path / "Adjacency" / "File0"), spa.csr_matrix(matrix))
    write_adjacency(str(tmp_path / "Adjacency" / "HomoAdjAuFile0"), spa.csr_matrix(matrix))
    write_adjacency(str(tmp_path / "Time_Dependent" / "HeAdjFile0"), spa.csr_matrix(matrix))

    out = tmp_path / "dense"
    assert main(["expand", str(tmp_path), "-o", str(out)]) == 0
    assert (out / "Adjacency" / "File0").is_file()
    assert (out / "Adjacency" / "HomoAdjAuFile0").is_file()
    assert (out / "Time_Dependent" / "HeAdjFile0").is_file()
