"""Reference data and the default geometry rules.

Lengths in Å. These values only set the *starting* geometry; relax before
reading anything physical from it. The rules match Yi's builder so that
structures from both tools agree for the same inputs.
"""
import math

from ase.data import covalent_radii, atomic_numbers

TRANSITION_METALS = [
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
]
X_ELEMENTS = ["C", "N"]                      # Greek layers
LATIN_ELEMENTS = TRANSITION_METALS + ["O"]   # metal or termination layers

# In-plane lattice constant a of M(n+1)X(n) MXenes.
# Ti, V, Cr, Zr, Nb, Mo, Hf, Ta carbides and Ti nitride: values from Yi's presets.
# The rest: rock-salt MX lattice constant / sqrt(2).
MXENE_A = {
    ("Sc", "C"): 3.30, ("Sc", "N"): 3.18,
    ("Ti", "C"): 3.07, ("Ti", "N"): 2.99,
    ("V", "C"): 2.90, ("V", "N"): 2.93,
    ("Cr", "C"): 2.83, ("Cr", "N"): 2.93,
    ("Y", "C"): 3.52, ("Y", "N"): 3.45,
    ("Zr", "C"): 3.27, ("Zr", "N"): 3.24,
    ("Nb", "C"): 3.13, ("Nb", "N"): 3.10,
    ("Mo", "C"): 2.99, ("Mo", "N"): 2.94,
    ("Hf", "C"): 3.27, ("Hf", "N"): 3.20,
    ("Ta", "C"): 3.10, ("Ta", "N"): 3.07,
    ("W", "C"): 3.00, ("W", "N"): 2.92,
}

# fcc metals: bulk lattice constant (Å) and Poisson ratio.
FCC_METALS = {
    "Au": (4.078, 0.44), "Ag": (4.086, 0.37), "Cu": (3.615, 0.34), "Pt": (3.924, 0.38),
    "Pd": (3.891, 0.39), "Ni": (3.524, 0.31), "Al": (4.050, 0.35),
}

# Spacing rules
D_MX_OVER_A = 0.385      # M–X layer spacing / a, close to relaxed Ti3C2
D_SAME_OVER_A = 0.816    # spacing of two adjacent layers of the same kind / a, sqrt(2/3)
BOND_M_O = 2.00          # M–O bond length used to place O layers
MIN_HEIGHT = 0.60        # floor for bond-derived heights
METAL_GAP = 2.30         # MXene–metal vertical gap
ONTOP_EXTRA_GAP = 0.30   # added when the first metal layer sits on top of the outer MXene atoms

PRESETS = {
    "Ti3C2O2": ["O", "Ti", "C", "Ti", "C", "Ti", "O"],
    "Ti2CO2": ["O", "Ti", "C", "Ti", "O"],
    "Ti3CNO2": ["O", "Ti", "C", "Ti", "N", "Ti", "O"],
    "Ti4N3O2": ["O", "Ti", "N", "Ti", "N", "Ti", "N", "Ti", "O"],
    "V2CO2": ["O", "V", "C", "V", "O"],
    "Nb2CO2": ["O", "Nb", "C", "Nb", "O"],
    "Mo2CO2": ["O", "Mo", "C", "Mo", "O"],
    "Mo2TiC2O2": ["O", "Mo", "C", "Ti", "C", "Mo", "O"],
    "Ti3C2 (bare)": ["Ti", "C", "Ti", "C", "Ti"],
}


def is_greek(element):
    return element in X_ELEMENTS


def mxene_a(element, x):
    """Default a for one M–X pair; estimated from covalent radii if not tabulated."""
    if (element, x) in MXENE_A:
        return MXENE_A[(element, x)], False
    bond = 0.92 * (covalent_radii[atomic_numbers[element]] + covalent_radii[atomic_numbers[x]])
    return round(math.sqrt(2) * bond, 3), True


def default_mxene_a(layers):
    """Mean of the M–X values over all metal and X layers (Vegard-like). Returns (a, estimated)."""
    metals = [e for e in layers if e in TRANSITION_METALS]
    xs = [e for e in layers if is_greek(e)] or ["C"]
    pairs = [mxene_a(m, x) for m in metals for x in xs]
    a = sum(p[0] for p in pairs) / len(pairs)
    return round(a, 4), any(p[1] for p in pairs)


def default_dz(lower, upper, a):
    """Vertical spacing between two adjacent MXene layers (element, site)."""
    (el1, site1), (el2, site2) = lower, upper
    if "O" in (el1, el2) and (el1 in TRANSITION_METALS or el2 in TRANSITION_METALS):
        lateral = 0.0 if site1 == site2 else a / math.sqrt(3)
        return math.sqrt(max(BOND_M_O**2 - lateral**2, MIN_HEIGHT**2))
    if is_greek(el1) != is_greek(el2):
        return D_MX_OVER_A * a
    return D_SAME_OVER_A * a


def default_metal_dz(metal_a, nu, strain):
    """fcc(111) layer spacing with a Poisson response to the in-plane strain."""
    d0 = metal_a / math.sqrt(3)
    d = d0 * (1 - 2 * nu / (1 - nu) * strain)
    return min(1.3 * d0, max(0.7 * d0, d))
