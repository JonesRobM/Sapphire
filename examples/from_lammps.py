"""Analyse a LAMMPS dump (or any ASE-readable trajectory) with Sapphire.

    python examples/from_lammps.py dump.lammpstrj --species Au Pt --frames 0 1000 10

ASE reads the dump; if it lacks chemical symbols (LAMMPS writes numeric types) map them with
--types 1=Au 2=Pt. The trajectory is written as extxyz next to the run output and `Process`
is driven exactly as for Sapphire-native data; results land in the standard run directory
(see docs/FILE_CONTRACT.md).
"""
import argparse
import os

from ase.io import read, write

from Sapphire import Process


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("trajectory")
    ap.add_argument("--format", default=None, help="ASE format name, e.g. lammps-dump-text (default: guess)")
    ap.add_argument("--types", nargs="*", default=[], help="LAMMPS type -> symbol, e.g. 1=Au 2=Pt")
    ap.add_argument("--species", nargs="*", default=None, help="species for Homo/Hetero analysis")
    ap.add_argument("--frames", nargs=3, type=int, default=[0, None, 1], metavar=("START", "END", "STEP"))
    ap.add_argument("--out", default="run_from_lammps/")
    a = ap.parse_args()

    frames = read(a.trajectory, index=":", format=a.format)
    if a.types:
        mapping = {int(k): v for k, v in (t.split("=") for t in a.types)}
        for atoms in frames:
            atoms.set_chemical_symbols([mapping[int(n)] for n in atoms.get_atomic_numbers()])
    os.makedirs(a.out, exist_ok=True)
    movie = os.path.join(a.out, "movie.xyz")
    write(movie, frames)

    start, end, step = a.frames
    species = a.species or sorted(set(frames[0].get_chemical_symbols()))
    bimetallic = len(species) > 1
    System = {"base_dir": os.path.join(a.out, ""), "movie_file_name": "movie.xyz", "extend_xyz": None,
              "Homo": species if bimetallic else None, "Hetero": bimetallic,
              "Start": start, "End": end or len(frames), "Step": step, "Skip": 1, "UniformPDF": False, "Band": 0.05}
    Quantities = {"Full": {"pdf": None, "rdf": None, "adj": None, "nn": None, "agcn": None, "cna_sigs": None,
                           "cna_patterns": None, "com": None, "comdist": None, "gyration": None}}
    if bimetallic:
        Quantities["Homo"] = {"hopdf": None, "hoadj": None, "hocomdist": None}
        Quantities["Hetero"] = {"hepdf": None, "headj": None, "mix": None, "lae": None}
    proc = Process.Process(System=System, Quantities=Quantities)
    proc.analyse({"JSD": ["pdf"]})
    print("done ->", a.out, "; keys:", sorted(proc.reader().available()))


if __name__ == "__main__":
    main()
