"""Inputs → atoms, built layer by layer from the bottom.

Every atom carries its layer index as its ASE tag, so constraints and analysis
refer to layers, never to atom indices.
"""
from collections import Counter
from dataclasses import dataclass, field

import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms

from . import lattice
from .elements import TRANSITION_METALS, is_greek
from .stack import label, metal_start_sites, next_site, resolve

MAX_ATOMS = 3000


@dataclass
class Layer:
    element: str
    site: str
    part: str                       # "mxene" or "metal"
    z: float
    shift: np.ndarray = field(default_factory=lambda: np.zeros(2))   # in-plane offset (Å), the metal RBT
    index: int = 0
    n_atoms: int = 0


@dataclass
class Interface:
    atoms: Atoms
    layers: list
    inputs: object                  # resolved Inputs
    match: lattice.Match
    warnings: list

    @property
    def formula(self):
        p = self.inputs
        mx = mxene_formula(p.mxene_layers)
        n = self.match.n_mx * p.repeat[0] * p.repeat[1]
        cells = f" ×{n}" if n > 1 else ""
        if p.metal_layers == 0:
            return mx + cells
        metal = f"{p.metal}{p.metal_layers}"
        return (f"{metal}/{mx}/{metal}" if p.geometry == "sandwich" else f"{metal}/{mx}") + cells

    def summary(self):
        p, m = self.inputs, self.match
        a, b, c, _, _, gamma = self.atoms.cell.cellpar()
        n_if = 0 if p.metal_layers == 0 else (1 if p.geometry == "slab" else 2)
        return dict(formula=self.formula, n_atoms=len(self.atoms),
                    cell=dict(a=round(float(a), 4), b=round(float(b), 4), c=round(float(c), 4), gamma=round(float(gamma), 3)),
                    area=round(float(np.linalg.norm(np.cross(*self.atoms.cell[:2]))), 4),
                    strain=round(100 * m.strain, 3) if p.metal_layers else None,
                    twist=round(m.twist, 3) if p.metal_layers else None,
                    n_interfaces=n_if)


def mxene_formula(elements):
    counts = Counter(elements)
    order = [e for e in dict.fromkeys(elements) if e in TRANSITION_METALS]
    order += [e for e in dict.fromkeys(elements) if is_greek(e)] + (["O"] if "O" in counts else [])
    return "".join(e + (str(counts[e]) if counts[e] > 1 else "") for e in order)


def stack_layers(p, m):
    """All layers bottom → top with their site and height (before vacuum is added)."""
    layers, z = [], 0.0
    for element, site, dz in zip(p.mxene_layers, p.mxene_sites, p.mxene_dz):
        z += dz
        layers.append(Layer(element, site, "mxene", z))
    if p.metal_layers:
        top_site, bottom_site = metal_start_sites(p)
        rbt = np.array(p.rbt, dtype=float) @ m.mxene_basis
        z_top, z_bottom = layers[-1].z, layers[0].z
        layers += [Layer(p.metal, next_site(top_site, k), "metal", z_top + p.gap + k * p.metal_dz, rbt)
                   for k in range(p.metal_layers)]
        if p.geometry == "sandwich":
            bottom = [Layer(p.metal, next_site(bottom_site, -k), "metal", z_bottom - p.gap - k * p.metal_dz, -rbt)
                      for k in range(p.metal_layers)]
            layers = bottom[::-1] + layers
    for i, layer in enumerate(layers):
        layer.index = i
    return layers


def build(inputs):
    p, warnings = resolve(inputs)
    m = lattice.match(p)
    layers = stack_layers(p, m)

    zs = [layer.z for layer in layers]
    padding = p.gap if p.geometry == "periodic" else p.vacuum
    c = max(zs) - min(zs) + padding
    for layer in layers:
        layer.z += padding / 2 - min(zs)

    a1, a2 = (np.linalg.norm(v) for v in m.supercell)
    cell = np.array([[a1, 0, 0], [-a2 / 2, a2 * lattice.SQ3 / 2, 0], [0, 0, c]])

    rows = []   # (species rank, layer, fy, fx, element, z); sorting on the first four fixes the atom order
    species = list(dict.fromkeys(layer.element for layer in layers))
    for layer in layers:
        basis = m.mxene_basis if layer.part == "mxene" else m.metal_basis
        offset = np.array(lattice.SITE_FRACTIONS[layer.site]) @ basis + layer.shift
        frac, expected = lattice.fill(basis, offset, m.supercell)
        if len(frac) != expected:
            raise RuntimeError(f"Layer {layer.index} ({layer.element}) gave {len(frac)} atoms, expected {expected}")
        layer.n_atoms = len(frac)
        rows += [(species.index(layer.element), layer.index, fy, fx, layer.element, layer.z) for fx, fy in frac]
    rows.sort(key=lambda r: r[:4])
    if len(rows) > MAX_ATOMS:
        raise ValueError(f"{len(rows)} atoms is above the limit of {MAX_ATOMS}")

    frac = np.array([[r[3], r[2], r[5] / c] for r in rows])
    atoms = Atoms([r[4] for r in rows], scaled_positions=frac, cell=cell, pbc=True, tags=[r[1] for r in rows])
    atoms.positions += 0.0      # turns -0.0 into 0.0 in the files

    bad = [i for i in p.fixed_layers if not 0 <= i < len(layers)]
    if bad:
        raise ValueError(f"fixed_layers {bad} are outside 0–{len(layers) - 1}")
    fixed = [i for i, t in enumerate(atoms.get_tags()) if t in set(p.fixed_layers)]
    if fixed:
        atoms.set_constraint(FixAtoms(indices=fixed))

    warnings += geometry_warnings(p, m)
    return Interface(atoms, layers, p, m, warnings)


def geometry_warnings(p, m):
    notes = []
    if p.metal_layers and abs(m.strain) > 0.05:
        notes.append(f"Metal in-plane strain is {100 * m.strain:.1f}%. A coincidence cell can lower it.")
    if p.metal_layers and p.gap < 1.2:
        notes.append(f"A gap of {p.gap:.2f} Å is very short and may give large initial forces.")
    if p.geometry != "periodic" and p.vacuum < 10:
        notes.append("Vacuum below 10 Å lets periodic images interact along c.")
    if p.metal_layers and p.geometry == "slab":
        notes.append("A one-sided slab is polar along c; use a dipole correction in DFT, or the sandwich geometry.")
    return notes


def layer_label(layer):
    return f"{layer.element}({label(layer.element, layer.site)})"
