"""Oxygen-reduction mass activity from the atop generalised coordination number (aGCN).

Model: Rossi, Asara & Baletto, *ChemPhysChem* **20**, 3037 (2019), doi:10.1002/cphc.201900564,
building on Rück, Bandarenka, Calle-Vallejo & Gagliardi, *J. Phys. Chem. Lett.* **9**, 4463 (2018).

* Eq. 5  GCN(j) = Σ_{i ∈ neighbours of j} CN(i) / CN_max,  CN_max = 12 (fcc atop site).
  Sapphire computes exactly this per atom as ``agcn`` (``Adjacent.agcn_generator``).
* Eq. 6  relative activity of a site with GCN = α (volcano, apex at α = 8.33, unity ≈ Pt(111) at 7.5)::

      A_l(α) = exp( 3.14 α − 23.40)   α ≤ 8.33
      A_r(α) = exp(−4.96 α + 42.18)   α > 8.33

  Only sites with GCN above a threshold count (paper: > 6; a stricter variant uses > 7.5).
* j̃_NP = Σ_sites A(α_i)  — current relative to one Pt(111) site.
* Eq. 7  MA_NP = j_flat / ρ_sites · j̃_NP / M_NP.  With j_flat = 2 mA cm⁻² and
  ρ_sites = 1.503·10¹⁵ cm⁻² this is 4.107 A mg⁻¹ · j̃_NP / N_NP for pure Pt. Here the mass is
  Σ_atoms m_i (species-weighted, from ASE), so the same expression covers alloys.

The volcano parameters are for Pt; for other metals pass your own ``branches``.
"""
from __future__ import annotations

import numpy as np
from ase.data import atomic_masses, atomic_numbers

J_FLAT_A_PER_CM2 = 2.0e-3          # Pt(111) reference current density, A cm^-2
SITE_DENSITY_PER_CM2 = 1.503e15    # atop sites on Pt(111), cm^-2
AMU_MG = 1.66053906660e-21         # 1 u in mg
PT_BRANCHES = ((3.14, -23.40), (-4.96, 42.18))   # (slope, intercept) of ln A: left, right


def branch_intersection(branches=PT_BRANCHES):
    """GCN at which the two printed branches are equal (continuous volcano)."""
    (ml, cl), (mr, cr) = branches
    return (cr - cl) / (ml - mr)


# DECISION (R. M. Jones, 2026-08-29): the paper states the apex is at GCN = 8.33, but the printed
# coefficients intersect at 8.096 (A = 6.9); at 8.33 they give 15.7 (left) vs 2.4 (right). Sapphire
# switches branches where they meet, so A(α) is continuous. Pass apex=8.33 for the literal reading.
# Recorded in docs/CHANGELOG.md.
GCN_APEX = branch_intersection()


def site_activity(gcn, branches=PT_BRANCHES, apex=None, gcn_min=6.0):
    """Relative ORR activity of each site (Eq. 6). Sites with GCN <= gcn_min contribute 0."""
    g = np.asarray(gcn, dtype=float)
    (ml, cl), (mr, cr) = branches
    if apex is None:
        apex = branch_intersection(branches)
    a = np.where(g <= apex, np.exp(ml * g + cl), np.exp(mr * g + cr))
    return np.where(g > gcn_min, a, 0.0)


def relative_current(gcn, **kw):
    """j̃_NP: summed relative activity of one frame's sites."""
    return float(site_activity(gcn, **kw).sum())


def prefactor_A_per_mg(masses_u):
    """j_flat / ρ_sites / M_NP  in A mg⁻¹ for a particle whose atoms have masses ``masses_u`` (u).
    For N Pt atoms this equals 4.107 / N."""
    m_mg = float(np.sum(masses_u)) * AMU_MG
    return J_FLAT_A_PER_CM2 / SITE_DENSITY_PER_CM2 / m_mg


def masses_from_symbols(symbols):
    return np.array([atomic_masses[atomic_numbers[s]] for s in symbols])


def mass_activity(gcn, symbols=None, masses_u=None, **kw):
    """Mass activity (A mg⁻¹) of one frame from per-atom aGCN and composition (Eq. 7)."""
    if masses_u is None:
        if symbols is None:
            raise ValueError("give either symbols or masses_u")
        masses_u = masses_from_symbols(symbols)
    return prefactor_A_per_mg(masses_u) * relative_current(gcn, **kw)


def mass_activity_series(agcn_frames, symbols, **kw):
    """Per-frame mass activity for a trajectory: ``agcn_frames`` is (frames, atoms)."""
    masses = masses_from_symbols(symbols)
    return np.array([prefactor_A_per_mg(masses) * relative_current(f, **kw) for f in agcn_frames])


def gcn_histogram(agcn_frames, bins=np.arange(0, 12.25, 0.25)):
    """Site-count histogram per frame on a fixed GCN grid -> (frames, bins). For heat maps."""
    return np.array([np.histogram(f, bins=bins)[0] for f in agcn_frames]), bins


def from_run(base_dir, gcn_min=6.0, write=True):
    """Compute the mass-activity series from a Sapphire run directory (needs ``agcn``).

    Reads ``Time_Dependent/AGCN`` and the trajectory's symbols via the run's ``Exec``/movie,
    writes ``Time_Dependent/Stats/MassActivity`` and returns (frames, MA in A/mg).
    """
    import os
    from ase.io import read
    from Sapphire.IO.Reader import Reader
    r = Reader(base_dir)
    agcn = np.asarray(r.load("agcn"), dtype=float)
    frames = r.frames("agcn")
    movie = next((f for f in os.listdir(base_dir) if f.endswith((".xyz", ".traj"))), None)
    if movie is None:
        raise FileNotFoundError("no trajectory file next to the run output to read symbols from")
    symbols = read(os.path.join(base_dir, movie), index=0).get_chemical_symbols()
    ma = mass_activity_series(agcn, symbols, gcn_min=gcn_min)
    if write:
        os.makedirs(os.path.join(base_dir, "Time_Dependent", "Stats"), exist_ok=True)
        with open(os.path.join(base_dir, "Time_Dependent", "Stats", "MassActivity"), "w") as f:
            for fr, v in zip(frames, ma):
                f.write(f"{fr} {v}\n")
    return frames, ma
