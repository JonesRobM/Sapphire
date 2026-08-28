"""Optical response: dielectric functions and pyGDM2-based spectra (requires the `light` extra)."""

from . import Epsilon_DFT as dft
from . import Epsilon_ExpClass as exp

__all__ = ['ClassSpec', 'Epsilon_DFT', 'Epsilon_ExpClass', 'dft', 'exp']

try:  # pyGDM2 is optional
    from .ClassSpec import Spectrum  # noqa: F401
    __all__.append('Spectrum')
except ImportError:  # pragma: no cover
    pass
