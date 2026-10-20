"""The workflow's inputs and outputs: one place, used by the web app, the CLI,
the papermill notebook and the Sim2L.

`None` means "auto": the value comes from the rules in elements.py and is
written out explicitly in resolved.json.
"""
import hashlib
import json
from dataclasses import asdict, dataclass, field, fields

from . import __version__
from .elements import FCC_METALS, LATIN_ELEMENTS, X_ELEMENTS

SCHEMA_VERSION = 1
SITES = ("A", "B", "C")
HEADS = ("omat_pbe", "oc20_usemppbe", "matpes_r2scan")


def spec(group, units="", choices=None, low=None, high=None, doc="", **kw):
    return field(metadata=dict(group=group, units=units, choices=choices, min=low, max=high, doc=doc), **kw)


@dataclass
class Inputs:
    # MXene, bottom → top
    mxene_layers: list[str] = spec("MXene", doc="Element of each layer, bottom → top (TM or O; C or N)",
                              default_factory=lambda: ["O", "Ti", "C", "Ti", "C", "Ti", "O"])
    stacking: str = spec("MXene", choices=("ABC", "ABA"), doc="AβC (octahedral) or AβA (trigonal prismatic)", default="ABC")
    mxene_sites: list[str | None] | None = spec("MXene", choices=SITES, doc="Site of each layer; null (or a null entry) = from stacking", default=None)
    mxene_dz: list[float | None] | None = spec("MXene", units="Å", doc="Spacing of each layer to the one below (first entry 0); null (or a null entry) = rules", default=None)
    mxene_a: float | None = spec("MXene", units="Å", low=2.0, high=4.5, doc="In-plane lattice constant; null = table", default=None)

    # Metal
    metal: str = spec("Metal", choices=tuple(FCC_METALS), doc="fcc metal, (111) planes parallel to the MXene", default="Au")
    metal_layers: int = spec("Metal", low=0, high=20, doc="Number of metal (111) layers; 0 = MXene only", default=5)
    metal_a: float | None = spec("Metal", units="Å", low=3.0, high=5.0, doc="Bulk fcc lattice constant; null = table", default=None)
    metal_dz: float | None = spec("Metal", units="Å", low=1.0, high=4.0, doc="Metal layer spacing; null = d111 with Poisson response", default=None)
    gap: float | None = spec("Metal", units="Å", low=0.8, high=8.0, doc="Vertical gap between MXene and metal; null = rule", default=None)
    metal_site: str | None = spec("Metal", choices=SITES, doc="Site of the first metal layer; null = continue the stacking", default=None)
    rbt: list[float] = spec("Metal", doc="Rigid-body translation of the metal, fractions of the MXene a1, a2",
                     default_factory=lambda: [0.0, 0.0])

    # Cell
    geometry: str = spec("Cell", choices=("slab", "sandwich", "periodic"),
                         doc="Metal on top; on both faces; or no vacuum (two interfaces through PBC)", default="slab")
    mxene_supercell: list[int] = spec("Cell", doc="MXene supercell vector [i, j] in the 120° hexagonal basis",
                                 default_factory=lambda: [1, 0])
    metal_supercell: list[int] = spec("Cell", doc="Metal supercell vector [k, l] matched onto it",
                                 default_factory=lambda: [1, 0])
    repeat: list[int] = spec("Cell", low=1, high=10, doc="Extra repeats along a1, a2", default_factory=lambda: [1, 1])
    vacuum: float = spec("Cell", units="Å", low=0.0, high=60.0, doc="Vacuum along c (slab, sandwich)", default=15.0)

    # Relaxation
    relax: bool = spec("Relax", doc="Relax with the MLIP; false = build only", default=True)
    model: str = spec("Relax", choices=("mh-1",), doc="MACE foundation model", default="mh-1")
    head: str = spec("Relax", choices=HEADS, doc="MACE-MH-1 output head", default="omat_pbe")
    dispersion: bool = spec("Relax", doc="Add D3(BJ) dispersion", default=False)
    fmax: float = spec("Relax", units="eV/Å", low=0.001, high=1.0, doc="Force convergence", default=0.05)
    max_steps: int = spec("Relax", low=1, high=10000, doc="Maximum optimizer steps", default=500)
    optimizer: str = spec("Relax", choices=("FIRE", "BFGS"), doc="ASE optimizer", default="FIRE")
    fixed_layers: list[int] | None = spec("Relax", doc="Layer indices (bottom = 0) held fixed; null = outer 2 metal layers", default=None)
    relax_cell: bool = spec("Relax", doc="Also relax a, b and γ (c stays fixed)", default=False)
    mlip_lattice: bool = spec("Relax", doc="Take auto lattice constants from the MLIP instead of the table", default=False)
    work_of_adhesion: bool = spec("Relax", doc="Rigid-separation work of adhesion", default=True)

    schema_version: int = SCHEMA_VERSION

    @classmethod
    def from_dict(cls, data):
        known = {f.name for f in fields(cls)}
        unknown = set(data) - known
        if unknown:
            raise ValueError(f"Unknown inputs: {', '.join(sorted(unknown))}")
        return cls(**data)

    def to_dict(self):
        return asdict(self)

    def canonical_json(self):
        return json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":"))

    def run_id(self):
        """Same inputs and same code version → same id (and the same run folder)."""
        text = self.canonical_json() + __version__
        return hashlib.sha256(text.encode()).hexdigest()[:12]


def validate(p):
    """Return a list of problems; empty when the inputs can be built."""
    errors = []
    for el in p.mxene_layers:
        if el not in LATIN_ELEMENTS and el not in X_ELEMENTS:
            errors.append(f"{el} is not an allowed MXene layer element")
    if not any(el not in X_ELEMENTS and el != "O" for el in p.mxene_layers):
        errors.append("The MXene needs at least one transition-metal layer")
    n = len(p.mxene_layers)
    for name in ("mxene_sites", "mxene_dz"):
        value = getattr(p, name)
        if value is not None and len(value) != n:
            errors.append(f"{name} has {len(value)} entries for {n} layers")
    if p.mxene_sites and any(s not in SITES for s in p.mxene_sites if s is not None):
        errors.append("Sites must be A, B or C")
    for f in fields(p):
        meta, value = f.metadata, getattr(p, f.name)
        if not meta or value is None:
            continue
        if meta["choices"] and not isinstance(value, list) and value not in meta["choices"]:
            errors.append(f"{f.name} must be one of {', '.join(meta['choices'])}")
        values = value if isinstance(value, list) else [value]
        for v in values:
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                if meta["min"] is not None and v < meta["min"] or meta["max"] is not None and v > meta["max"]:
                    errors.append(f"{f.name} = {v} is outside [{meta['min']}, {meta['max']}]")
    i, j = p.mxene_supercell
    k, l = p.metal_supercell
    if i * i - i * j + j * j < 1 or k * k - k * l + l * l < 1:
        errors.append("Supercell vectors must be non-zero")
    return errors


def describe_inputs():
    """Inputs as rows: name, group, type, default, units, allowed, description."""
    rows = []
    for f in fields(Inputs):
        if not f.metadata:
            continue
        default = f.default_factory() if callable(f.default_factory) else f.default
        m = f.metadata
        allowed = ", ".join(m["choices"]) if m["choices"] else (
            f"{m['min']} – {m['max']}" if m["min"] is not None else "")
        rows.append(dict(name=f.name, group=m["group"], type=f.type.__name__ if isinstance(f.type, type) else str(f.type),
                         default=default, units=m["units"], allowed=allowed, description=m["doc"]))
    return rows


# name: (type, units, description)
OUTPUTS = {
    "run_id": ("str", "", "Hash of inputs.json and the code version; also the run folder name"),
    "formula": ("str", "", "Interface label, e.g. Au5/Ti3C2O2"),
    "n_atoms": ("int", "", "Atoms in the cell"),
    "cell": ("dict", "Å, °", "a, b, c, gamma"),
    "area": ("float", "Å²", "Interface area of the cell"),
    "strain": ("float", "%", "In-plane strain of the metal (the MXene keeps its a)"),
    "twist": ("float", "°", "Rotation of the metal lattice relative to the MXene"),
    "n_interfaces": ("int", "", "1 for slab, 2 for sandwich and periodic"),
    "mlip_lattice": ("dict", "Å", "MLIP-derived mxene_a and metal_a (when mlip_lattice is on)"),
    "layers": ("list", "Å", "Per layer: element, site, part, atoms, z and spacing before (and after) relaxation"),
    "gaps_initial": ("dict", "Å", "MXene–metal gaps before relaxation"),
    "gaps_relaxed": ("dict", "Å", "MXene–metal gaps after relaxation"),
    "energy_initial": ("float", "eV", "MLIP energy of the starting structure"),
    "energy_final": ("float", "eV", "MLIP energy after relaxation"),
    "converged": ("bool", "", "fmax reached within max_steps"),
    "n_steps": ("int", "", "Optimizer steps taken"),
    "fmax_final": ("float", "eV/Å", "Largest residual force"),
    "work_of_adhesion": ("float", "J/m²", "(E_metal + E_MXene − E_interface) / (n_interfaces · area)"),
    "energy_metal": ("float", "eV", "Metal slab alone, relaxed-interface geometry"),
    "energy_mxene": ("float", "eV", "MXene alone, relaxed-interface geometry"),
    "trace": ("list", "eV, eV/Å", "Energy and fmax per step"),
    "files": ("dict", "", "Structure files: initial.* and relaxed.* (POSCAR, extxyz, cif, lammps-data)"),
    "warnings": ("list", "", "Notes about the starting geometry"),
}
