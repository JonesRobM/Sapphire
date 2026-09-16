"""Analysing frames across processes must not change the answer.

Per-frame work is independent and the cross-frame quantities live in ``Process.analyse``,
which reads finished results off disk -- so the only risks are mechanical: interleaved
writes to the shared per-frame output files, and the CNA masterkey, which a serial run
carries forward from frame to frame.
"""
import numpy as np
import pytest

from Sapphire.parallel import plan_chunks


# ------------------------------------------------------------------ chunk planning
@pytest.mark.parametrize("start, end, step, jobs, expected_frames", [
    (0, 12, 1, 4, [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]]),
    (0, 10, 1, 3, [[0, 1, 2, 3], [4, 5, 6], [7, 8, 9]]),       # remainder to the early chunks
    (0, 5, 1, 1, [[0, 1, 2, 3, 4]]),
    (0, 20, 5, 2, [[0, 5], [10, 15]]),                          # stride preserved
    (3, 9, 2, 2, [[3, 5], [7]]),                                # non-zero start
])
def test_chunks_cover_every_frame_exactly_once(start, end, step, jobs, expected_frames):
    chunks = plan_chunks(start, end, step, jobs)
    got = [list(range(lo, hi, st)) for lo, hi, st in chunks]
    assert got == expected_frames
    flat = [f for block in got for f in block]
    assert flat == list(range(start, end, step))
    assert len(flat) == len(set(flat)), "a frame must not be analysed twice"


def test_more_jobs_than_frames_is_capped():
    chunks = plan_chunks(0, 3, 1, 16)
    assert len(chunks) == 3
    assert [f for lo, hi, st in chunks for f in range(lo, hi, st)] == [0, 1, 2]


def test_no_frames_gives_no_chunks():
    assert plan_chunks(0, 0, 1, 4) == []


# --------------------------------------------------------------- serial == parallel
@pytest.fixture(scope="module")
def serial_and_parallel(tmp_path_factory):
    """The same six frames, analysed serially and across three processes."""
    from Sapphire.Tutorials import data
    from Sapphire.api import run

    work = tmp_path_factory.mktemp("par")
    data.sample("AuPt", work)
    traj = str(work / "AuPt_sample.xyz")
    quantities = ["pdf", "adj", "nn", "agcn", "cna_sigs", "collect", "concert"]
    common = dict(quantities=quantities, frames=(0, 6, 1), homo=[], hetero=[])
    ser = run(traj, str(work / "serial"), jobs=1, **common)
    par = run(traj, str(work / "parallel"), jobs=3, **common)
    return ser, par


@pytest.mark.parametrize("key", ["nn", "agcn", "pdf", "rcut", "collect", "concert"])
def test_series_match(serial_and_parallel, key):
    ser, par = serial_and_parallel
    np.testing.assert_allclose(np.asarray(ser.load(key), dtype=float),
                               np.asarray(par.load(key), dtype=float))


def test_adjacency_matrices_match(serial_and_parallel):
    ser, par = serial_and_parallel
    np.testing.assert_array_equal(ser.load("adj"), par.load("adj"))
    np.testing.assert_array_equal(ser.frames("adj"), par.frames("adj"))


def test_frames_are_in_order_and_complete(serial_and_parallel):
    """Workers finish out of order; the merge has to restore frame order."""
    _, par = serial_and_parallel
    assert par.frames("adj").tolist() == [0, 1, 2, 3, 4, 5]
    assert np.asarray(par.load("nn")).shape[0] == 6


def test_cna_signatures_carry_the_same_counts(serial_and_parallel):
    """A serial run's signature table is ragged when the masterkey grows mid-run; the
    merged one is rectangular. The counts per signature must still agree."""
    ser, par = serial_and_parallel
    skey, pkey = ser.masterkey(), par.masterkey()
    s, p = np.asarray(ser.load("cna_sigs")), np.asarray(par.load("cna_sigs"))
    assert set(skey) <= set(pkey)
    for frame in range(s.shape[0]):
        s_counts = dict(zip(skey, np.asarray(s[frame]).ravel()))
        p_counts = dict(zip(pkey, np.asarray(p[frame]).ravel()))
        for sig, count in s_counts.items():
            assert p_counts[sig] == count, f"frame {frame}, signature {sig}"


def test_worker_scratch_is_cleaned_up(serial_and_parallel, tmp_path_factory):
    _, par = serial_and_parallel
    assert not (par.base / ".sapphire_workers").exists()


# ------------------------------------------------------------------------- validation
def test_jobs_below_one_is_rejected(tmp_path):
    from Sapphire.api import Config
    traj = tmp_path / "t.xyz"
    traj.write_text("1\n\nAu 0.0 0.0 0.0\n")
    with pytest.raises(ValueError, match="jobs"):
        Config(trajectory=str(traj), jobs=0).validate()


# ------------------------------------------------------- bimetallic (the real case)
@pytest.fixture(scope="module")
def bimetallic_serial_and_parallel(tmp_path_factory):
    """A two-species run, serial and across three processes.

    Worth its own fixture: the homo and hetero quantities write per-frame matrices that
    a monometallic run never produces, and the hetero ones (HeAdjFile7.npz) land in
    Time_Dependent/ rather than Adjacency/. A merge that classified per-frame files by
    directory treated those as text and corrupted them.
    """
    from Sapphire.Tutorials import data
    from Sapphire.api import run

    work = tmp_path_factory.mktemp("bipar")
    data.sample("AuPt", work)
    traj = str(work / "AuPt_sample.xyz")
    common = dict(
        quantities=["pdf", "adj", "nn", "agcn", "cna_sigs"],
        homo=["hopdf", "hoadj", "hocomdist", "homobonds"],
        hetero=["hepdf", "headj", "mix", "lae", "ele_nn", "heterobonds"],
        species=["Au", "Pt"], frames=(0, 6, 1),
    )
    return (run(traj, str(work / "s"), jobs=1, **common),
            run(traj, str(work / "p"), jobs=3, **common))


@pytest.mark.parametrize("key", ["adj", "hoadjAu", "hoadjPt", "headj"])
def test_bimetallic_matrices_survive_the_merge(bimetallic_serial_and_parallel, key):
    """headj lives in Time_Dependent/, so it is the one a directory-based merge broke."""
    ser, par = bimetallic_serial_and_parallel
    np.testing.assert_array_equal(ser.load(key), par.load(key))


@pytest.mark.parametrize("key", ["nn", "agcn", "mix", "ele_nnAu", "homo_bondsAu", "hetero_bonds"])
def test_bimetallic_series_match(bimetallic_serial_and_parallel, key):
    ser, par = bimetallic_serial_and_parallel
    np.testing.assert_allclose(np.asarray(ser.load(key), dtype=float),
                               np.asarray(par.load(key), dtype=float))


def test_bimetallic_produces_the_same_keys(bimetallic_serial_and_parallel):
    ser, par = bimetallic_serial_and_parallel
    assert set(ser.available()) == set(par.available())
