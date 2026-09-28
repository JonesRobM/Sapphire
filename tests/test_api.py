import numpy as np
import pytest

from Sapphire import api
from tests.conftest import AU561


def test_config_validation(tmp_path):
    with pytest.raises(ValueError):
        api.Config(str(AU561), quantities=["nope"]).validate()
    with pytest.raises(FileNotFoundError):
        api.Config("missing.xyz").validate()


def test_config_toml_roundtrip(tmp_path):
    cfg = api.Config(str(AU561), out_dir=str(tmp_path), frames=(0, 10, 2), quantities=["pdf", "nn"], statistics={"JSD": ["pdf"]})
    path = cfg.to_toml(tmp_path / "c.toml")
    back = api.Config.from_toml(path)
    assert back.frames == (0, 10, 2) and back.quantities == ["pdf", "nn"] and back.statistics == {"JSD": ["pdf"]}


def test_run_monometallic(tmp_path):
    r = api.run(AU561, tmp_path / "out", quantities=["pdf", "adj", "nn", "agcn"], frames=(0, 1, 1))
    assert r.load("nn").shape == (1, 561) and (tmp_path / "out" / "sapphire_config.toml").exists()
    assert r.load("agcn").max() <= 12


def test_run_bimetallic_defaults(tmp_path):
    from Sapphire.Tutorials import data
    xyz = data.sample("AuPt", tmp_path)
    r = api.run(xyz, tmp_path / "out", quantities=["pdf", "adj", "nn"], frames=(0, 70, 69), statistics={"JSD": ["pdf"]})
    keys = r.available()
    assert {"hopdfAu", "hopdfPt", "hepdf", "mix", "laeAu", "JSDpdf"} <= set(keys)
    assert np.asarray(r.load("mix")).shape == (2,)


def test_statistics_do_not_load_unused_matrices(tmp_path, monkeypatch):
    # analyse() used to read back every result, and per-frame matrices come back dense:
    # frames * N^2 ints, ~11 GB for a 1400-frame, 1415-atom movie.
    from Sapphire.Tutorials import data
    from Sapphire.IO.Reader import Reader
    xyz = data.sample("AuPt", tmp_path)
    loaded = []
    real = Reader.load
    monkeypatch.setattr(Reader, "load", lambda self, key: loaded.append(key) or real(self, key))
    r = api.run(xyz, tmp_path / "out", quantities=["pdf", "adj", "nn"], frames=(0, 4, 1),
                homo=[], hetero=[], statistics={"JSD": ["pdf"]})
    assert "adj" not in loaded and "pdf" in loaded
    assert len(r.load("JSDpdf")) == 4


def test_every_frame_logs_its_timing_from_a_nonzero_start(tmp_path):
    # The per-frame timing indexed atom counts by frame number rather than by position in
    # the analysed slice, overrunning the list for any run (or parallel chunk) not starting at 0.
    from Sapphire.Tutorials import data
    xyz = data.sample("AuPt", tmp_path)
    api.run(xyz, tmp_path / "out", quantities=["pdf"], frames=(2, 6, 1), homo=[], hetero=[])
    info = (tmp_path / "out" / "Sapphire_Info.txt").read_text()
    assert info.count("for each atom") == 4
