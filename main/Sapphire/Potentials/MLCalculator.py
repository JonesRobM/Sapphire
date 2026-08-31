"""Machine-learned interatomic potentials as ASE calculators.

Sapphire is a post-processing package; it does not train potentials. What it offers is a
one-liner to obtain a *foundation-model* calculator so that a cluster can be relaxed or run
in MD before analysis — replacing the toy EMT potential in ``Tutorials.data.synthetic``
and the classical Gupta parameters in ``Potentials.Gupta`` when accuracy matters.

Install with ``pip install -e ".[mlpot]"`` (MACE + torch, ~2 GB with CUDA wheels).

>>> from Sapphire.Potentials.MLCalculator import ml_calculator
>>> calc = ml_calculator("mace-mp-0")       # downloads weights on first use (~100 MB)
>>> atoms.calc = calc
"""
from __future__ import annotations

KNOWN = {
    # name -> (loader, kwargs)  All return an ase.calculators.calculator.Calculator.
    "mace-mp-0": ("mace.calculators", "mace_mp", {"model": "medium", "default_dtype": "float32"}),
    "mace-mp-0-small": ("mace.calculators", "mace_mp", {"model": "small", "default_dtype": "float32"}),
    "mace-mpa-0": ("mace.calculators", "mace_mp", {"model": "medium-mpa-0", "default_dtype": "float32"}),
    "mace-off": ("mace.calculators", "mace_off", {"model": "medium", "default_dtype": "float32"}),
}


def available() -> bool:
    try:
        import mace  # noqa: F401
        return True
    except ImportError:
        return False


def ml_calculator(name: str = "mace-mp-0", device: str = "cpu", **overrides):
    """Return an ASE calculator for a named foundation model.

    Raises ImportError with install instructions if the ``mlpot`` extra is missing.
    """
    if name not in KNOWN:
        raise ValueError(f"unknown model {name!r}; choose from {sorted(KNOWN)}")
    module, func, kwargs = KNOWN[name]
    try:
        import importlib
        loader = getattr(importlib.import_module(module), func)
    except ImportError as e:  # pragma: no cover - depends on optional extra
        raise ImportError("machine-learned potentials need the 'mlpot' extra: pip install -e '.[mlpot]'") from e
    return loader(device=device, **{**kwargs, **overrides})


def calculator_or_emt(name: str = "mace-mp-0", **kw):
    """The named ML calculator if installed, else ASE's EMT (with a note). Handy in tutorials/CI."""
    if available():
        return ml_calculator(name, **kw), name
    from ase.calculators.emt import EMT
    return EMT(), "EMT (install the 'mlpot' extra for %s)" % name
