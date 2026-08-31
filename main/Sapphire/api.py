"""A compact front door to Sapphire.

The historical interface is two nested dictionaries (``System`` / ``Quantities``) handed to
:class:`Sapphire.Process.Process`. That remains fully supported; this module wraps it in a typed
:class:`Config` and a one-call :func:`run` that returns a :class:`~Sapphire.IO.Reader.Reader`.

>>> from Sapphire.api import run
>>> reader = run("movie.xyz", "out/", quantities=["pdf", "adj", "nn", "agcn", "cna_sigs"], frames=(0, 100, 5))
>>> reader.load("agcn").shape
(20, 561)

Configs can be saved/loaded as TOML (``Config.to_toml`` / ``Config.from_toml``) so an analysis
is reproducible from one small text file.
"""
from __future__ import annotations

import dataclasses
import os
import pathlib
import sys
from dataclasses import dataclass, field

from Sapphire.Utilities.Supported import Supported

_SUP = Supported()
FULL_KEYS, HOMO_KEYS, HETERO_KEYS = set(_SUP.Full()), set(_SUP.Homo()), set(_SUP.Hetero())
FULL_KEYS |= {"collect", "concert", "euc", "pos"}
STAT_KEYS = {"JSD", "Kullback", "PStat"}

DEFAULT_QUANTITIES = ["pdf", "rdf", "adj", "nn", "agcn", "cna_sigs", "cna_patterns", "com", "comdist", "gyration"]
DEFAULT_HOMO = ["hopdf", "hoadj", "hocomdist", "homobonds"]
DEFAULT_HETERO = ["hepdf", "headj", "mix", "lae", "ele_nn", "heterobonds"]


@dataclass
class Config:
    """Everything needed for one analysis run."""

    trajectory: str
    out_dir: str = "sapphire_run/"
    frames: tuple = (0, None, 1)                  # (start, end, step); end=None -> all
    quantities: list = field(default_factory=lambda: list(DEFAULT_QUANTITIES))
    homo: list | None = None                       # per-species quantities (needs species)
    hetero: list | None = None                     # between-species quantities (bimetallic only)
    species: list | None = None                    # e.g. ["Au", "Pt"]; inferred from frame 0 if None
    band: float = 0.05                             # KDE bandwidth (Å) for the PDDF cutoff
    statistics: dict = field(default_factory=dict) # e.g. {"JSD": ["pdf"]}
    extend_xyz: list | None = None                 # per-atom quantities to write as extra xyz columns
    strict: bool = False
    overwrite: bool = True

    # ------------------------------------------------------------------ validation
    def validate(self):
        bad = [q for q in self.quantities if q not in FULL_KEYS]
        bad += [q for q in (self.homo or []) if q not in HOMO_KEYS]
        bad += [q for q in (self.hetero or []) if q not in HETERO_KEYS]
        bad_stats = [s for s in self.statistics if s not in STAT_KEYS]
        if bad or bad_stats:
            raise ValueError(f"unknown quantities {bad} / statistics {bad_stats}; see Sapphire.Utilities.Supported")
        if not os.path.isfile(self.trajectory):
            raise FileNotFoundError(self.trajectory)
        return self

    # ------------------------------------------------------------------ conversion
    def to_legacy(self, n_frames=None, species=None):
        """The (System, Quantities) dictionaries `Process` consumes."""
        species = species or self.species
        start, end, step = self.frames
        bimetallic = species is not None and len(species) > 1
        homo = self.homo if self.homo is not None else (DEFAULT_HOMO if bimetallic else None)
        hetero = self.hetero if self.hetero is not None else (DEFAULT_HETERO if bimetallic else None)
        system = {
            "base_dir": os.path.join(os.path.abspath(self.out_dir), ""),
            "movie_file_name": os.path.basename(self.trajectory),
            "extend_xyz": self.extend_xyz,
            "Homo": list(species) if (bimetallic and homo) else None,
            "Hetero": bool(bimetallic and hetero),
            "Start": start, "End": end if end is not None else n_frames, "Step": step, "Skip": 1,
            "UniformPDF": False, "Band": self.band,
        }
        quantities = {"Full": {q: None for q in self.quantities}}
        if system["Homo"]:
            quantities["Homo"] = {q: None for q in homo}
        if system["Hetero"]:
            quantities["Hetero"] = {q: None for q in hetero}
        return system, quantities

    def to_dict(self):
        return dataclasses.asdict(self)

    def to_toml(self, path):
        import json
        lines, tables = [], []
        for k, v in self.to_dict().items():
            if v is None:
                continue
            if isinstance(v, (list, tuple)):
                lines.append(f"{k} = {json.dumps([x for x in v if x is not None])}")
            elif isinstance(v, dict):
                tables.append(f"\n[{k}]")   # tables must follow every top-level key
                tables += [f"{kk} = {json.dumps(vv)}" for kk, vv in v.items()]
            elif isinstance(v, bool):
                lines.append(f"{k} = {'true' if v else 'false'}")
            elif isinstance(v, str):
                lines.append(f"{k} = {json.dumps(v)}")
            else:
                lines.append(f"{k} = {v}")
        pathlib.Path(path).write_text("\n".join(lines + tables) + "\n")
        return path

    @classmethod
    def from_toml(cls, path):
        if sys.version_info >= (3, 11):
            import tomllib
        else:  # pragma: no cover
            import tomli as tomllib
        d = tomllib.loads(pathlib.Path(path).read_text())
        if "frames" in d:
            f = list(d["frames"]) + [None] * (3 - len(d["frames"]))
            d["frames"] = tuple(f[:3])
        return cls(**d)


def run(trajectory, out_dir="sapphire_run/", **kwargs):
    """Analyse ``trajectory`` into ``out_dir`` and return a Reader over the results.

    Keyword arguments are :class:`Config` fields. The trajectory is copied/converted next to the
    output if it is not already an xyz/traj file ASE can index.
    """
    cfg = Config(trajectory=str(trajectory), out_dir=str(out_dir), **kwargs)
    return run_config(cfg)


def run_config(cfg: Config):
    import shutil
    from ase.io import read
    from Sapphire import Process

    cfg.validate()
    frames = read(cfg.trajectory, index=":")
    species = cfg.species or sorted(set(frames[0].get_chemical_symbols()))
    base = os.path.join(os.path.abspath(cfg.out_dir), "")
    os.makedirs(base, exist_ok=True)
    target = os.path.join(base, os.path.basename(cfg.trajectory))
    if os.path.abspath(cfg.trajectory) != os.path.abspath(target):
        shutil.copy(cfg.trajectory, target)
    system, quantities = cfg.to_legacy(n_frames=len(frames), species=species)
    proc = Process.Process(System=system, Quantities=quantities, strict=cfg.strict, overwrite=cfg.overwrite)
    if cfg.statistics or cfg.extend_xyz or "collect" in cfg.quantities:
        proc.analyse(cfg.statistics or {})
    cfg.to_toml(os.path.join(base, "sapphire_config.toml"))
    return proc.reader()
