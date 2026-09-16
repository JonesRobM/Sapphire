"""Tests for the ``sapphire`` command line front end."""
import argparse
import os

import pytest

from Sapphire.api import REQUIRES, resolve_dependencies
from Sapphire.cli import main, parse_frames


# ---------------------------------------------------------------- frame slices
@pytest.mark.parametrize("spec, expected", [
    (None, (0, None, 1)),
    ("0:100:10", (0, 100, 10)),
    ("5", (5, None, 1)),
    ("5:", (5, None, 1)),
    ("5:10", (5, 10, 1)),
    ("::2", (0, None, 2)),
    (":100:", (0, 100, 1)),
])
def test_parse_frames_accepts(spec, expected):
    assert parse_frames(spec) == expected


@pytest.mark.parametrize("spec", ["10:5:1", "0:10:0", "0:10:-1", "a:b", "1:2:3:4"])
def test_parse_frames_rejects(spec):
    with pytest.raises(argparse.ArgumentTypeError):
        parse_frames(spec)


# ---------------------------------------------------------- dependency closure
def test_cna_pulls_in_its_prerequisites():
    """`adj` needs the cutoff from `pdf`; without it Process writes nothing."""
    resolved, added = resolve_dependencies(["cna_sigs"])
    assert resolved == ["pdf", "adj", "cna_sigs"]
    assert set(added) == {"pdf", "adj"}


def test_dependencies_precede_their_dependents():
    resolved, _ = resolve_dependencies(["cna_patterns"])
    for q, deps in REQUIRES.items():
        if q in resolved:
            for dep in deps:
                assert resolved.index(dep) < resolved.index(q), f"{dep} must precede {q}"


def test_already_complete_request_is_untouched():
    resolved, added = resolve_dependencies(["pdf", "adj", "agcn"])
    assert resolved == ["pdf", "adj", "agcn"]
    assert added == []


def test_independent_quantities_need_nothing():
    resolved, added = resolve_dependencies(["com", "gyration"])
    assert resolved == ["com", "gyration"]
    assert added == []


# ------------------------------------------------------------------ exit codes
def test_missing_trajectory_exits_2(capsys):
    assert main(["run", "definitely_not_here.xyz", "-q", "pdf"]) == 2
    assert "no such trajectory" in capsys.readouterr().err


def test_unknown_quantity_exits_2(tmp_path, capsys):
    traj = tmp_path / "t.xyz"
    traj.write_text("1\n\nAu 0.0 0.0 0.0\n")
    assert main(["run", str(traj), "-q", "not_a_quantity"]) == 2
    assert "unknown quantities" in capsys.readouterr().err


def test_jobs_reaches_the_plan(tmp_path, capsys):
    """--jobs is parsed and carried into the run rather than silently ignored."""
    traj = tmp_path / "t.xyz"
    traj.write_text("1\n\nAu 0.0 0.0 0.0\n")
    assert main(["run", str(traj), "-q", "com", "--jobs", "4", "--dry-run"]) == 0


def test_jobs_below_one_is_refused(tmp_path, capsys):
    traj = tmp_path / "t.xyz"
    traj.write_text("1\n\nAu 0.0 0.0 0.0\n")
    assert main(["run", str(traj), "-q", "com", "--jobs", "0"]) == 2
    assert "jobs" in capsys.readouterr().err


def test_dry_run_reports_plan_without_computing(tmp_path, capsys):
    traj = tmp_path / "t.xyz"
    traj.write_text("1\n\nAu 0.0 0.0 0.0\n")
    out = tmp_path / "out"
    assert main(["run", str(traj), "-o", str(out), "-q", "cna_sigs", "--dry-run"]) == 0
    printed = capsys.readouterr().out
    assert "pdf, adj, cna_sigs" in printed
    assert not out.exists(), "--dry-run must not create the output directory"


def test_quantities_listing(capsys):
    assert main(["quantities"]) == 0
    out = capsys.readouterr().out
    assert "cna_sigs" in out and "requires pdf, adj" in out


# ----------------------------------------------------------------- quote opt-out
def test_no_quote_sets_the_env_var(tmp_path, monkeypatch):
    monkeypatch.delenv("SAPPHIRE_NO_QUOTE", raising=False)
    traj = tmp_path / "t.xyz"
    traj.write_text("1\n\nAu 0.0 0.0 0.0\n")
    main(["run", str(traj), "-q", "com", "--no-quote", "--dry-run"])
    assert os.environ.get("SAPPHIRE_NO_QUOTE") == "1"


def test_quote_is_skipped_without_network(monkeypatch):
    """The quote must never reach the network when opted out, and never raise."""
    from Sapphire.Utilities import Initial
    monkeypatch.setenv("SAPPHIRE_NO_QUOTE", "1")
    assert Initial.Info("")._quote_() == Initial.Info._NO_QUOTE


def test_quote_survives_a_broken_lookup(monkeypatch):
    """A network failure decorating the log header must not fail the analysis."""
    import sys
    import types
    from Sapphire.Utilities import Initial

    monkeypatch.delenv("SAPPHIRE_NO_QUOTE", raising=False)
    broken = types.ModuleType("wikiquote")

    def explode(*_a, **_k):
        raise OSError("no route to host")

    broken.quotes = explode
    broken.random_titles = explode
    monkeypatch.setitem(sys.modules, "wikiquote", broken)
    assert Initial.Info("")._quote_() == Initial.Info._NO_QUOTE


def test_info_reports_the_package_version():
    from Sapphire import __version__
    from Sapphire.Utilities import Initial
    assert __version__ in Initial.Info("")._version_()
