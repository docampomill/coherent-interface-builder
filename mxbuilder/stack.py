"""Stacking sites and the resolution of "auto" inputs into explicit numbers.

Sites are the three positions of a close-packed plane in the hexagonal cell:
A = (0, 0), B = (1/3, 2/3), C = (2/3, 1/3) in fractions of a1, a2.
Greek letters (α, β, γ) are the same sites for C/N layers.
"""
from dataclasses import replace

from . import lattice
from .elements import (FCC_METALS, METAL_GAP, ONTOP_EXTRA_GAP, default_dz, default_metal_dz,
                       default_mxene_a, is_greek)
from .schema import SITES, validate

GREEK = {"A": "α", "B": "β", "C": "γ"}


def label(element, site):
    return GREEK[site] if is_greek(element) else site


def next_site(site, step=1):
    return SITES[(SITES.index(site) + step) % 3]


def third_site(a, b):
    """The site used by neither neighbour; with one neighbour, the next site in the cycle."""
    if b is None or a == b:
        return next_site(a)
    return next(s for s in SITES if s not in (a, b))


def mxene_sites(layers, stacking):
    """Default site of each MXene layer, bottom → top.

    The core (everything between the outer O layers) follows the stacking:
    ABC cycles A → B → C through every layer (AβCαB…), ABA keeps metals on A
    and C/N on B (AβAβA…). Each outer O sits on the site free of its two
    nearest layers, which for ABC is the fcc hollow over the inner metal.
    """
    lo = next(i for i, el in enumerate(layers) if el != "O")
    hi = max(i for i, el in enumerate(layers) if el != "O")
    sites = [None] * len(layers)
    for k in range(lo, hi + 1):
        if stacking == "ABC":
            sites[k] = SITES[(k - lo) % 3]
        else:
            sites[k] = "B" if is_greek(layers[k]) else "A"
    for i in range(lo - 1, -1, -1):
        sites[i] = third_site(sites[i + 1], sites[i + 2] if i + 2 < len(layers) else None)
    for i in range(hi + 1, len(layers)):
        sites[i] = third_site(sites[i - 1], sites[i - 2] if i >= 2 else None)
    return sites


def metal_start_sites(p):
    """First-layer sites of the top and bottom metal slabs.

    Auto continues the stacking from the outer two MXene layers. The bottom
    slab mirrors the top one, so a symmetric MXene gives an inversion-symmetric sandwich.
    """
    s, n = p.mxene_sites, len(p.mxene_sites)
    top_auto = third_site(s[-1], s[-2] if n > 1 else None)
    bottom_auto = third_site(s[0], s[1] if n > 1 else None)
    if p.metal_site is None:
        return top_auto, bottom_auto
    shift = SITES.index(p.metal_site) - SITES.index(top_auto)
    return p.metal_site, next_site(bottom_auto, -shift)


def default_fixed_layers(p):
    """Outer two metal layers of each slab (none in periodic cells)."""
    if p.geometry == "periodic" or p.metal_layers == 0:
        return []
    n_fix = min(2, p.metal_layers - 1)
    n_mx, n_m = len(p.mxene_layers), p.metal_layers
    if p.geometry == "sandwich":
        total = n_mx + 2 * n_m
        return list(range(n_fix)) + list(range(total - n_fix, total))
    total = n_mx + n_m
    return list(range(total - n_fix, total))


def resolve(p):
    """Fill every null input. Returns (explicit Inputs, warnings)."""
    errors = validate(p)
    if errors:
        raise ValueError("; ".join(errors))
    warnings = []

    # MXene: a null list, or a null entry in a list, means "use the rule"
    auto = mxene_sites(p.mxene_layers, p.stacking)
    p = replace(p, mxene_sites=[given or rule for given, rule in zip(p.mxene_sites or auto, auto)])
    if p.mxene_a is None:
        a, estimated = default_mxene_a(p.mxene_layers)
        if estimated:
            warnings.append("No tabulated a for this composition; a was estimated from covalent radii. Set it or use MLIP lattice constants.")
        p = replace(p, mxene_a=a)
    pairs = list(zip(p.mxene_layers, p.mxene_sites))
    given = p.mxene_dz or [None] * len(pairs)
    dz = [0.0] + [given[i] if given[i] is not None else round(default_dz(pairs[i - 1], pairs[i], p.mxene_a), 6)
                  for i in range(1, len(pairs))]
    p = replace(p, mxene_dz=dz)

    if p.metal_a is None:
        p = replace(p, metal_a=FCC_METALS[p.metal][0])
    if p.metal_site is None:
        p = replace(p, metal_site=metal_start_sites(p)[0])
    if p.gap is None:
        on_top = lattice.match(p).coherent and p.rbt == [0, 0] and p.metal_site == p.mxene_sites[-1]
        p = replace(p, gap=METAL_GAP + (ONTOP_EXTRA_GAP if on_top else 0.0))
    if p.metal_dz is None:
        strain = lattice.match(p).strain
        p = replace(p, metal_dz=round(default_metal_dz(p.metal_a, FCC_METALS[p.metal][1], strain), 6))
    if p.fixed_layers is None:
        p = replace(p, fixed_layers=default_fixed_layers(p))
    return p, warnings
