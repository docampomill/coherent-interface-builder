"""Hexagonal-lattice geometry: how the metal is matched onto the MXene and how
each layer is filled with atoms.

In-plane vectors use the 120° hexagonal basis a1 = a(1, 0), a2 = a(-1/2, √3/2).
A supercell vector [i, j] means i·a1 + j·a2. The second supercell vector is the
first one rotated by 120°, so every supercell stays hexagonal.
"""
import cmath
import math
from dataclasses import dataclass

import numpy as np

SQ3 = math.sqrt(3)
ROT120 = complex(-0.5, SQ3 / 2)
SITE_FRACTIONS = {"A": (0.0, 0.0), "B": (1 / 3, 2 / 3), "C": (2 / 3, 1 / 3)}


def to_xy(z):
    return [z.real, z.imag]


def hex_vector(i, j, a=1.0):
    return a * (i + j * ROT120)


def n_cells(i, j):
    return i * i - i * j + j * j


def wrap_twist(theta):
    """Map an angle to (-30°, 30°], the range that is distinct for two hexagonal lattices."""
    theta = theta % 60.0
    return theta - 60.0 if theta > 30.0 + 1e-9 else theta


@dataclass
class Match:
    supercell: np.ndarray     # rows A1, A2 in Å, MXene frame (repeats included)
    mxene_basis: np.ndarray   # rows a1, a2
    metal_basis: np.ndarray   # rows b1, b2, already strained and rotated
    strain: float             # metal in-plane strain, fraction
    twist: float              # metal rotation relative to the MXene, degrees
    n_mx: int                 # MXene primitive cells per supercell (before repeats)
    n_metal: int              # metal primitive cells per supercell (before repeats)

    @property
    def coherent(self):
        return abs(self.twist) < 1e-6 and self.n_mx == self.n_metal


def match(p):
    """Supercell and strained metal lattice for resolved inputs."""
    i, j = p.mxene_supercell
    k, l = p.metal_supercell
    s1 = hex_vector(i, j, p.mxene_a)
    f = s1 / hex_vector(k, l)        # complex scale-and-rotate that maps the metal vector onto s1
    n1, n2 = p.repeat
    return Match(
        supercell=np.array([to_xy(s1 * n1), to_xy(s1 * ROT120 * n2)]),
        mxene_basis=np.array([to_xy(hex_vector(1, 0, p.mxene_a)), to_xy(hex_vector(0, 1, p.mxene_a))]),
        metal_basis=np.array([to_xy(f), to_xy(f * ROT120)]),
        strain=abs(f) / (p.metal_a / math.sqrt(2)) - 1,
        twist=wrap_twist(math.degrees(cmath.phase(f))),
        n_mx=n_cells(i, j),
        n_metal=n_cells(k, l),
    )


def search(a_mx, a_metal_bulk, atoms_per_mx, atoms_per_metal, max_strain=0.03, max_atoms=400,
           twist=None, twist_tol=1.0, limit=15):
    """Coincidence cells of two hexagonal lattices, metal strained onto the MXene.

    Returns the smallest cells first, keeping only those that lower |strain|
    (a Pareto front of size vs strain). `twist` restricts the metal rotation.
    """
    a_m = a_metal_bulk / math.sqrt(2)
    n_max = max_atoms // max(atoms_per_mx, 1)
    i_max = int(math.sqrt(4 * n_max / 3)) + 1
    vectors = [(i, j) for i in range(1, i_max + 1) for j in range(i) if n_cells(i, j) <= n_max]
    found = []
    for i, j in vectors:
        n_mx = n_cells(i, j)
        length = a_mx * math.sqrt(n_mx)
        for k, l in vectors:
            n_m = n_cells(k, l)
            atoms = n_mx * atoms_per_mx + n_m * atoms_per_metal
            if atoms > max_atoms or math.gcd(math.gcd(i, j), math.gcd(k, l)) > 1:
                continue
            strain = length / (a_m * math.sqrt(n_m)) - 1
            if abs(strain) > max_strain:
                continue
            theta = wrap_twist(math.degrees(cmath.phase(hex_vector(i, j)) - cmath.phase(hex_vector(k, l))))
            if twist is not None and abs(wrap_twist(theta - twist)) > twist_tol:
                continue
            found.append(dict(mxene=[i, j], metal=[k, l], strain=round(100 * strain, 3),
                              twist=round(theta, 3), n_atoms=atoms, n_mx=n_mx, n_metal=n_m))
    found.sort(key=lambda c: (c["n_atoms"], abs(c["strain"])))
    front, seen, best = [], set(), math.inf
    for c in found:
        key = (c["n_mx"], c["n_metal"], c["twist"])
        if key in seen:
            continue
        seen.add(key)
        if abs(c["strain"]) < best - 1e-6:
            front.append(c)
            best = abs(c["strain"])
    return front[:limit]


def fill(basis, offset, supercell, tol=1e-6):
    """Every lattice point n·b1 + m·b2 + offset inside the supercell.

    Returns fractional coordinates in [0, 1), sorted, and the expected count.
    """
    corners = np.array([[0, 0], [1, 0], [0, 1], [1, 1]]) @ supercell
    nm = (corners - offset) @ np.linalg.inv(basis)
    lo = np.floor(nm.min(axis=0)).astype(int) - 1
    hi = np.ceil(nm.max(axis=0)).astype(int) + 1
    n, m = np.meshgrid(np.arange(lo[0], hi[0] + 1), np.arange(lo[1], hi[1] + 1), indexing="ij")
    points = np.stack([n.ravel(), m.ravel()], axis=1) @ basis + offset
    frac = points @ np.linalg.inv(supercell)
    frac -= np.floor(frac + tol)
    frac[frac < 0] = 0.0
    frac = np.unique(np.round(frac, 12), axis=0)
    expected = round(abs(np.linalg.det(supercell) / np.linalg.det(basis)))
    return frac, expected
