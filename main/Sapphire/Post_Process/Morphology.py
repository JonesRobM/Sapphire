"""Shell-by-shell morphology of a finite cluster (successor to the `Emerald` prototype).

Everything operates on one frame: an ``(N, 3)`` position array plus either a nearest-neighbour
cutoff or a ready-made adjacency matrix. All routines are vectorised with numpy.

Concepts (Emerald's definitions, kept):

* ``cn``          coordination number = neighbours within ``r_cut``.
* ``gcn``         generalised coordination number = Σ_neighbours cn / 12 (fcc atop convention).
* surface atom    cn ≤ ``surface_cn`` (default 10); everything else is core.
* peeling         remove the surface, recompute on what is left, repeat → concentric shells.
* shell thickness difference in radial extent between the outer shell and the next one.
* faceting ratio  fraction of atoms that are surface atoms (%).
* surface area    Σ_surface 4π r_at² (1 − gcn/12)   (Emerald) — cf. ACS Catal. 10, 3911 (2020).
* volume          sphere fitted inside the surface plus one cone per surface atom (Emerald's estimate).
* solid angle     per surface atom, ρ̄² π / r², with ρ̄ the mean neighbour distance.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial.distance import cdist


def pair_distances(positions):
    return cdist(positions, positions)


def adjacency(positions, r_cut):
    """Boolean (N, N) adjacency for |r_ij| <= r_cut, zero diagonal."""
    d = pair_distances(positions)
    a = d <= r_cut
    np.fill_diagonal(a, False)
    return a


def coordination(adj):
    return np.asarray(adj).sum(axis=1)


def generalised_coordination(adj, cn_max=12):
    a = np.asarray(adj, dtype=float)
    return a @ coordination(adj) / cn_max


def centre(positions):
    return np.asarray(positions, dtype=float).mean(axis=0)


def radii(positions, origin=None):
    p = np.asarray(positions, dtype=float)
    return np.linalg.norm(p - (centre(p) if origin is None else origin), axis=1)


def surface_mask(adj, surface_cn=10):
    return coordination(adj) <= surface_cn


def peel(positions, r_cut, surface_cn=10, max_shells=50):
    """Successively strip surface atoms.

    Returns a list of index arrays (into the original positions), outermost shell first; the
    last entry is whatever is left when nothing qualifies as surface any more (the core).
    For a perfect n-shell Mackay icosahedron this reproduces the 10k²+2 shell populations.
    """
    pos = np.asarray(positions, dtype=float)
    remaining = np.arange(len(pos))
    shells = []
    for _ in range(max_shells):
        if len(remaining) == 0:
            break
        a = adjacency(pos[remaining], r_cut)
        surf = surface_mask(a, surface_cn)
        if not surf.any() or surf.all():
            shells.append(remaining)
            return shells
        shells.append(remaining[surf])
        remaining = remaining[~surf]
    if len(remaining):
        shells.append(remaining)
    return shells


def shell_thickness(positions, r_cut, surface_cn=10):
    """(Δr_min, Δr_max): radial gaps between the outer shell and the shell beneath it."""
    shells = peel(positions, r_cut, surface_cn)
    if len(shells) < 2:
        return 0.0, 0.0
    r = radii(positions)
    r1, r2 = r[shells[0]], r[shells[1]]
    return float(r1.min() - r2.min()), float(r1.max() - r2.max())


def faceting_ratio(adj, surface_cn=10):
    """Percentage of atoms on the surface."""
    m = surface_mask(adj, surface_cn)
    return 100.0 * m.sum() / len(m)


def surface_area(adj, r_atom, surface_cn=10, cn_max=12):
    """Σ_surface 4π r_atom² (1 − gcn/12)  (Å² if r_atom in Å)."""
    gcn = generalised_coordination(adj, cn_max)
    m = surface_mask(adj, surface_cn)
    return float(np.sum(4 * np.pi * r_atom**2 * (1 - gcn[m] / cn_max)))


def solid_angle(positions, adj, surface_cn=10):
    """Emerald's per-surface-atom solid angle: π ρ̄² / r², ρ̄ = mean neighbour distance."""
    pos = np.asarray(positions, dtype=float)
    d = pair_distances(pos)
    m = surface_mask(adj, surface_cn)
    r = radii(pos)
    out = np.zeros(len(pos))
    for i in np.flatnonzero(m):
        nb = np.flatnonzero(adj[i])
        if len(nb) > 2 and r[i] > 0:
            rho = d[i, nb].sum() / (len(nb) - 2)   # Emerald divides by m-2
            out[i] = np.pi * rho**2 / r[i] ** 2
    return out


def volume(positions, adj, r_atom, surface_cn=10):
    """Emerald's volume estimate: inscribed sphere + a cone per surface atom + half-atom caps."""
    pos = np.asarray(positions, dtype=float)
    m = surface_mask(adj, surface_cn)
    r = radii(pos)[m]
    r_min = r.min() - r_atom
    v_core = 4 / 3 * np.pi * r_min**3
    v_caps = 2 / 3 * np.pi * r_atom**3 * len(r)
    h = r - r_min
    inner = np.sqrt(np.maximum((np.sqrt(r**2 + r_atom**2) - r_min) ** 2 - h**2, 0.0))
    r2 = r_atom - inner
    v_cones = np.sum(np.pi / 3 * (r_atom**2 + r_atom * r2 + r2**2) * h)
    return float(v_core + v_caps + v_cones)


def describe(positions, r_cut, r_atom, surface_cn=10):
    """One-call summary dictionary for a frame."""
    a = adjacency(positions, r_cut)
    shells = peel(positions, r_cut, surface_cn)
    return {
        "n_atoms": len(positions),
        "cn": coordination(a),
        "gcn": generalised_coordination(a),
        "surface_mask": surface_mask(a, surface_cn),
        "shell_populations": [len(s) for s in shells],
        "faceting_ratio_pct": faceting_ratio(a, surface_cn),
        "surface_area": surface_area(a, r_atom, surface_cn),
        "shell_thickness": shell_thickness(positions, r_cut, surface_cn),
        "volume": volume(positions, a, r_atom, surface_cn),
    }
