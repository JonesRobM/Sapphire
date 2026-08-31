import pytest

from Sapphire.Potentials import MLCalculator as ML
from Sapphire.Tutorials import data


def test_unknown_model_is_rejected():
    with pytest.raises(ValueError):
        ML.ml_calculator("not-a-model")


def test_missing_extra_gives_install_hint():
    if ML.available():
        pytest.skip("mace installed")
    with pytest.raises(ImportError, match="mlpot"):
        ML.ml_calculator("mace-mp-0")


def test_fallback_and_synthetic_accepts_calculator():
    calc, name = ML.calculator_or_emt()
    frames = data.synthetic(n_frames=3, noshells=2, calculator=calc)
    assert len(frames) == 3 and len(frames[0]) == 13
