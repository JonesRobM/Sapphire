"""Access to the tutorial MD trajectories.

The full data set (four bimetallic systems, two 14 ns runs each, ~232 MB
compressed) is too large to live in git. Three tiers are offered:

* :func:`sample`   — a down-sampled trajectory shipped with the package
  (every 20th frame of one run per system, ~1 MB each). Works offline.
* :func:`fetch`    — download the full tarball for a system from the GitHub
  release assets into ``~/.cache/sapphire`` and return the extracted directory.
* :func:`synthetic`— build a short, seeded Langevin trajectory of a small
  cluster with ASE's EMT calculator. No download, deterministic, ideal for tests.
"""
from __future__ import annotations

import gzip
import hashlib
import os
import pathlib
import shutil
import tarfile
import urllib.request

SYSTEMS = ("AgNi", "AuPd", "AuPt", "CuNi")

RELEASE_URL = "https://github.com/JonesRobM/Sapphire/releases/download/data-v1/{system}.tar.gz"

# SHA-256 of the original tarballs (computed 2026-08-28). Fill in after upload if they change.
SHA256 = {
    "AgNi": "cea6dedbf1b77699f44a1d38994b0a3775b36c5ff741848adba51bc6d256b88d",
    "AuPd": "16562723b869d21f1435746a0264ee21c616a26386b484165936b1d3b75173c3",
    "AuPt": "ae497a9cc4016c477e10a5ea0a90ee3b0fdd4765ea192fd1a6b8615f5cf970a9",
    "CuNi": "91ff5663bc37e48c9b817201b6acb6a0c00e95381732a3245e9ef415fe5b2068",
}

_HERE = pathlib.Path(__file__).resolve().parent
_CACHE = pathlib.Path(os.environ.get("SAPPHIRE_CACHE", pathlib.Path.home() / ".cache" / "sapphire"))


def _check(system: str) -> None:
    if system not in SYSTEMS:
        raise ValueError(f"unknown system {system!r}; choose from {SYSTEMS}")


def sample(system: str, out_dir: str | os.PathLike | None = None) -> pathlib.Path:
    """Return the path to the bundled down-sampled ``.xyz`` for *system*.

    If *out_dir* is given the (gzip-compressed) sample is decompressed there and
    that path is returned, ready for ``System['movie_file_name']``.
    """
    _check(system)
    src = _HERE / "Data" / "samples" / f"{system}_sample.xyz.gz"
    if not src.exists():
        raise FileNotFoundError(f"sample for {system} missing at {src}")
    if out_dir is None:
        return src
    out = pathlib.Path(out_dir) / f"{system}_sample.xyz"
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(src, "rb") as f_in, open(out, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    return out


def fetch(system: str, cache: str | os.PathLike | None = None, force: bool = False) -> pathlib.Path:
    """Download and extract the full trajectory set for *system*; returns the directory."""
    _check(system)
    cache = pathlib.Path(cache) if cache else _CACHE
    cache.mkdir(parents=True, exist_ok=True)
    tgz = cache / f"{system}.tar.gz"
    target = cache / f"{system}_Traj"
    if target.exists() and not force:
        return target
    url = RELEASE_URL.format(system=system)
    print(f"Sapphire: downloading {url} -> {tgz}")
    urllib.request.urlretrieve(url, tgz)
    expected = SHA256.get(system)
    if expected:
        got = hashlib.sha256(tgz.read_bytes()).hexdigest()
        if got != expected:
            tgz.unlink()
            raise OSError(f"checksum mismatch for {system}: {got} != {expected}")
    with tarfile.open(tgz) as tar:
        tar.extractall(cache, filter="data")
    return target


def synthetic(n_frames: int = 50, element: str = "Cu", noshells: int = 2,
              temperature_K: float = 600.0, seed: int = 0, out: str | os.PathLike | None = None):
    """Short seeded Langevin MD of an icosahedral cluster with ASE's EMT potential.

    Returns the list of :class:`ase.Atoms`; also writes an ``.xyz`` if *out* is given.
    Requires only ``ase`` (EMT supports Al, Cu, Ag, Au, Ni, Pd, Pt).
    """
    import numpy as np
    from ase import units
    from ase.calculators.emt import EMT
    from ase.cluster import Icosahedron
    from ase.io import write
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

    rng = np.random.RandomState(seed)
    atoms = Icosahedron(element, noshells=noshells)
    atoms.calc = EMT()
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature_K, rng=rng)
    dyn = Langevin(atoms, 5 * units.fs, temperature_K=temperature_K, friction=0.02, rng=rng)
    frames = []
    for _ in range(n_frames):
        dyn.run(10)
        frames.append(atoms.copy())
    if out is not None:
        write(out, frames)
    return frames
