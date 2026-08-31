"""Regenerate the exemplar figure set in assets/ from the bundled AuPt sample.

    python examples/make_figures.py [--frames 10] [--out assets/]

Runs `Sapphire.Process` (Full + Homo + Hetero quantities) on the down-sampled AuPt
melting trajectory, then renders the standard `Graphing.Plot_Funcs` set via the
`IO.Reader`-backed `Graphing.Read_Meta`. Provenance is recorded in assets/LOG.md.
"""
import argparse
import os
import pathlib
import tempfile
import time

import numpy as np

from Sapphire import Process
from Sapphire.Graphing import Plot_Funcs
from Sapphire.Graphing import Reader as GReader
from Sapphire.Tutorials import data

QUANTITIES = {
    "Full": {"pdf": None, "rdf": None, "com": None, "comdist": None, "adj": None, "nn": None, "agcn": None,
             "cna_sigs": None, "cna_patterns": None, "gyration": None, "stat_radius": None, "surf_area": None,
             "surf_atoms": None, "moi": None, "collect": None, "concert": None},
    "Homo": {"hopdf": None, "hordf": None, "hocom": None, "hoadj": None, "hocomdist": None, "homobonds": None},
    "Hetero": {"hepdf": None, "herdf": None, "headj": None, "mix": None, "lae": None, "ele_nn": None,
               "heterobonds": None},
}
FIGURES = {
    "agcn_heat": {"Name": "agcn_heat.png"},
    "prdf_plot": {"Names": ["pdf", "rdf"], "Frames": [0, 5, 9], "He": True, "Ho": ["Au", "Pt"]},
    "plot_stats": {"Stats": ["JSD", "Kullback"], "Quants": ["pdf", "rdf"], "Species": ["Au", "Pt"]},
    "com_plot_bi": {"Dists": ["CoMDist", "MidCoMDist"], "Species": ["Au", "Pt"], "Frames": [0, 9]},
    "cna_plot": {"Frames": [0, 9]},
    "agcn_histo": {"Frames": [0, 9]},
    "com_full_plot": {"Frames": [0, 9]},
    "cum_com": {"Frames": [0, 5, 9]},
    "cna_traj": {"Sigs": [(4, 2, 1), (4, 2, 2), (5, 5, 5), (3, 1, 1), (2, 0, 0)]},
    "h_c": {},
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", type=int, default=10, help="frames to analyse out of the 70 in the sample")
    ap.add_argument("--out", default="assets/", help="figure directory")
    ap.add_argument("--workdir", default=None, help="where the Process run writes (default: temp dir)")
    a = ap.parse_args()
    out = str(pathlib.Path(a.out).resolve()) + "/"
    work = a.workdir or tempfile.mkdtemp(prefix="sapphire_fig_")
    data.sample("AuPt", work)
    step = max(1, 70 // a.frames)
    t0 = time.time()
    System = {"base_dir": os.path.join(work, ""), "movie_file_name": "AuPt_sample.xyz", "extend_xyz": None,
              "Homo": ["Au", "Pt"], "Hetero": True,
              "Start": 0, "End": 70, "Step": step, "Skip": 1, "UniformPDF": False, "Band": 0.05}
    proc = Process.Process(System=System, Quantities=QUANTITIES)
    proc.analyse({"JSD": ["pdf", "rdf"], "Kullback": ["pdf"]})
    print(f"Process: {time.time() - t0:.0f}s, {len(proc.All_Times)} frames -> {work}")
    # sample = every 20th MD frame = 200 ps; melting ramp +50 K / ns from 300 K
    t_ps = np.asarray(proc.All_Times) * 200.0
    PS = {"base_dir": System["base_dir"], "plot_dir": out, "frame_dt": 200.0, "temperature": 300.0 + 50.0 * t_ps / 1000.0}
    meta, err = GReader.Read_Meta(PS).Average()
    Plot_Funcs.Plot_Funcs(meta, Quantities=FIGURES, System=PS, Errors=err).Make_Plots()
    print("figures:", len([f for f in os.listdir(out) if f.endswith(".png")]), "->", out)


if __name__ == "__main__":
    main()
