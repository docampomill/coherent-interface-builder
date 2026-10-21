"""mxene_au_box.py with its variables as inputs. MXenes and metals are the rows of the spreadsheet.

    python builder.py [inputs.json]      writes POSCAR, POSCAR_<mxene>, POSCAR_<metal> and outputs.json
"""
import hashlib
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import openpyxl
from pymatgen.core import Element, Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

SPREADSHEET = Path(__file__).parent / "data" / "MXene_FCC_L12_strain_matrix.xlsx"

DEFAULTS = dict(
    mxene="Ti3C2",       # MXenes sheet, Formula
    stacking="ABC",      # ABC = AβC (octahedral), ABA = AβA (trigonal prismatic)
    vacuum=20.0,         # Å, below the MXene only
    metal="Au",          # Metals sheet, Symbol
    metal_layers=4,      # metal (001) atomic planes
    gap=2.3,             # Å, top MXene layer to first metal layer
    shift=[0.0, 0.0],    # Å, in-plane shift of the metal along x and y
    tol=0.02,            # largest |strain| allowed on the metal along x and along y
    max_replicas=10,     # largest number of MXene replicas tried per direction
)

SITES = {"A": (0, 0), "B": (0, 1 / 3), "C": (1 / 2, 1 / 6)}    # hexagonal A, B, C sites in the a × √3a cell
D_MX = 0.385                                                  # M–X layer spacing / a
BOND_MO = 2.0                                                 # Å, M–O bond that sets the O layer height
SPACE_GROUP = {"FCC": ("Fm-3m", 225), "L12": ("Pm-3m", 221)}


def mxene_layers(formula):
    """Bottom → top elements of an ordered MXene (Ti3C2, Mo2TiC2, Ti2VCNO2); None for solid solutions and Tx."""
    if not re.fullmatch(r"([A-Z][a-z]?\d*)+", formula):
        return None
    tokens = [(e, int(n or 1)) for e, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula)]
    if not all(Element.is_valid_symbol(e) for e, _ in tokens):
        return None
    metals = [(e, n) for e, n in tokens if e not in ("C", "N", "O")]
    xs = [e for e, n in tokens if e in ("C", "N") for _ in range(n)]
    outer = next((e for e, n in metals if n == 2 and len(metals) > 1), metals[0][0])    # outer metal: the one with count 2
    inner = [e for e, n in metals for _ in range(n) if e != outer] if len(metals) > 1 else [outer] * (len(xs) - 1)
    ms = [outer] + inner + [outer]
    core = [ms[0]] + [e for pair in zip(xs, ms[1:]) for e in pair]
    o = ["O"] if dict(tokens).get("O") else []
    return o + core + o


def read_spreadsheet(path=SPREADSHEET):
    book = openpyxl.load_workbook(path, read_only=True, data_only=True)
    mxenes, metals = {}, {}
    for f, family, termination, a, source, *_ in book["MXenes"].iter_rows(min_row=3, values_only=True):
        if f and isinstance(a, (int, float)) and mxene_layers(str(f)):
            mxenes[f] = dict(formula=f, family=family, termination=termination, a=float(a), source=source)
    for symbol, name, structure, a, _, source, *_ in book["Metals"].iter_rows(min_row=3, values_only=True):
        if symbol and isinstance(a, (int, float)) and re.fullmatch(r"[A-Z][a-z]?(3[A-Z][a-z]?)?", symbol):
            metals[symbol] = dict(symbol=symbol, name=name, structure=structure, a=float(a), source=source)
    book.close()
    return mxenes, metals


MXENES, METALS = read_spreadsheet()
SHA256 = hashlib.sha256(SPREADSHEET.read_bytes()).hexdigest()


def build(inputs):
    """inputs: dict, missing keys take DEFAULTS. Returns (interface Structure, outputs dict)."""
    p = {**DEFAULTS, **inputs}
    if p["mxene"] not in MXENES or p["metal"] not in METALS:
        raise ValueError(f"{p['mxene']} or {p['metal']} is not a row of {SPREADSHEET.name}")

    # ---------------- MXene: conventional (orthogonal) cell, vacuum only below ----------------
    row = MXENES[p["mxene"]]
    a, vac = row["a"], p["vacuum"]
    layers = mxene_layers(row["formula"])

    # AβC: sites cycle A → B → C with the central layer on A; AβA: metals on A, C/N on B.
    # O sits on the site free of its two nearest layers.
    core = [i for i, e in enumerate(layers) if e != "O"]
    centre = (core[0] + core[-1]) // 2
    sites = [None if e == "O" else "ABC"[(i - centre) % 3] if p["stacking"] == "ABC" else "B" if e in ("C", "N") else "A"
             for i, e in enumerate(layers)]
    if layers[0] == "O":
        sites[0] = ({"A", "B", "C"} - {sites[1], sites[2]}).pop()
    if layers[-1] == "O":
        sites[-1] = ({"A", "B", "C"} - {sites[-2], sites[-3]}).pop()

    heights = [0.0]                                    # above the bottom layer, Å
    for lower, upper in zip(layers, layers[1:]):
        heights.append(heights[-1] + (math.sqrt(BOND_MO**2 - a**2 / 3) if "O" in (lower, upper) else D_MX * a))

    c = heights[-1] + vac                              # top layer sits exactly at z = c
    lattice = Lattice.orthorhombic(a, np.sqrt(3) * a, c)
    species, coords, layer = [], [], []
    for i, (s, site, h) in enumerate(zip(layers, sites, heights)):
        x, y = SITES[site]
        for dx, dy in [(0, 0), (0.5, 0.5)]:            # the rectangular cell holds two copies of each layer
            species.append(s)
            coords.append([(x + dx) % 1, (y + dy) % 1, (vac + h) / c])
            layer.append(i)
    mxene = Structure(lattice, species, coords, to_unit_cell=False, site_properties={"layer": layer})

    # ---------------- Metal: conventional cubic cell from its space group ----------------
    metal_row = METALS[p["metal"]]
    a_m = metal_row["a"]
    symbol, number = SPACE_GROUP[metal_row["structure"]]
    if metal_row["structure"] == "FCC":                # metal at 4a, 4 atoms
        metal = Structure.from_spacegroup(symbol, Lattice.cubic(a_m), [metal_row["symbol"]], [[0, 0, 0]])
    else:                                              # L12 A3B: B at 1a, A at 3c
        A, B = re.fullmatch(r"([A-Z][a-z]?)3([A-Z][a-z]?)", metal_row["symbol"]).groups()
        metal = Structure.from_spacegroup(symbol, Lattice.cubic(a_m), [B, A], [[0, 0, 0], [0, 0.5, 0.5]])
    assert SpacegroupAnalyzer(metal, symprec=1e-3).get_space_group_number() == number

    # ---------------- In-plane replicas giving (approximately) the same box ----------------
    replicas = []                                      # per direction: (MXene replicas, metal replicas, strain on metal)
    for L_mx, L_m in zip(mxene.lattice.abc[:2], metal.lattice.abc[:2]):
        for m in range(1, p["max_replicas"] + 1):
            n = max(1, round(m * L_mx / L_m))
            strain = m * L_mx / (n * L_m) - 1
            if abs(strain) <= p["tol"]:
                replicas.append((m, n, strain))
                break
        else:
            raise ValueError(f"No match within {p['tol']:.1%} up to {p['max_replicas']} replicas; raise the tolerance or max replicas")
    (mx_x, m_x, eps_x), (mx_y, m_y, eps_y) = replicas

    # ---------------- Location: metal (001) on MXene (0001), metal [100] along x ----------------
    gap, shift, n_layers = p["gap"], p["shift"], p["metal_layers"]
    metal_super = metal.copy()
    metal_super.make_supercell([m_x, m_y, math.ceil(n_layers / 2)])       # two (001) planes per cubic cell
    metal_super.remove_sites([k for k, s in enumerate(metal_super) if s.coords[2] > (n_layers - 1) * a_m / 2 + 1e-3])

    box_x, box_y = mx_x * mxene.lattice.a, mx_y * mxene.lattice.b
    z_top = c + gap + metal_super.cart_coords[:, 2].max()                  # vacuum only below the MXene
    interface_lattice = Lattice.orthorhombic(box_x, box_y, z_top)

    species, coords, layer, part = [], [], [], []
    for i in range(mx_x):                              # MXene replicas, copied by hand: make_supercell
        for j in range(mx_y):                          # would wrap the top layer from z = c to z = 0
            for s in mxene:
                species.append(s.species_string)
                coords.append(s.coords + [i * mxene.lattice.a, j * mxene.lattice.b, 0])
                layer.append(s.properties["layer"])
                part.append("mxene")

    for s in metal_super:                              # same fractional x, y in the matched box = strained metal
        species.append(s.species_string)
        x = (s.frac_coords[0] + shift[0] / box_x) % 1 * box_x
        y = (s.frac_coords[1] + shift[1] / box_y) % 1 * box_y
        coords.append([x, y, c + gap + s.coords[2]])
        layer.append(len(layers) + round(s.coords[2] / (a_m / 2)))
        part.append("metal")

    interface = Structure(interface_lattice, species, coords, coords_are_cartesian=True, to_unit_cell=False,
                          site_properties={"layer": layer, "part": part})

    # ---------------- Outputs ----------------
    z = interface.cart_coords[:, 2]
    stack = []
    for i in range(len(layers) + n_layers):
        on = [k for k, L in enumerate(layer) if L == i]
        names = [species[k] for k in on]
        stack.append(dict(index=i, part=part[on[0]], site=sites[i] if i < len(layers) else None,
                          species={e: names.count(e) for e in dict.fromkeys(names)}, z=round(float(z[on].mean()), 4)))

    r = lambda v, n=4: round(float(v), n)
    outputs = dict(
        inputs=p,
        mxene=dict(**row, stacking=p["stacking"], layers=layers, sites=sites, cell=[r(a), r(np.sqrt(3) * a), r(c)],
                   thickness=r(heights[-1])),
        metal=dict(**metal_row, space_group=f"{symbol} ({number})", layers=n_layers,
                   strained_cell=[r(box_x / m_x), r(box_y / m_y), r(a_m)]),
        match=dict(x=dict(mxene=mx_x, metal=m_x, length=r(box_x), strain=r(eps_x, 6)),
                   y=dict(mxene=mx_y, metal=m_y, length=r(box_y), strain=r(eps_y, 6))),
        interface=dict(formula=interface.composition.formula, n_atoms=len(interface), box=[r(box_x), r(box_y), r(z_top)],
                       layers=stack),
        spreadsheet=dict(file=SPREADSHEET.name, sha256=SHA256),
        files={"POSCAR": interface.get_sorted_structure().to(fmt="poscar"),
               f"POSCAR_{p['mxene']}": mxene.to(fmt="poscar"),
               f"POSCAR_{p['metal']}": metal.to(fmt="poscar")},
    )
    return interface, outputs


if __name__ == "__main__":
    inputs = json.loads(Path(sys.argv[1]).read_text()) if len(sys.argv) > 1 else {}
    _, outputs = build(inputs)
    for name, text in outputs["files"].items():
        Path(name).write_text(text)
    Path("outputs.json").write_text(json.dumps(outputs, indent=2) + "\n")
    print(outputs["interface"]["formula"], outputs["match"])
