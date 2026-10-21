"""MACE relaxation of the interface built by builder.py.

    python relax.py [inputs.json]      writes POSCAR, POSCAR_relaxed, POSCAR_<mxene>, POSCAR_<metal> and outputs.json
"""
import functools
import json
import sys
from importlib import metadata
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.mixing import SumCalculator
from ase.constraints import FixAtoms
from ase.filters import FrechetCellFilter
from ase.optimize import FIRE
from pymatgen.core import Structure

import builder

DEFAULTS = dict(
    model="mh-1",              # MACE foundation model
    head="omat_pbe",           # output head of the MH models
    dispersion=True,           # add D3(BJ)
    fmax=0.01,                 # eV/Å
    max_steps=500,             # FIRE steps
    fixed_layers=2,            # outermost metal layers held fixed
    relax_cell=True,           # relax a, b and γ; c stays fixed
    work_of_adhesion=False,    # from the rigidly separated MXene and metal
)
HEADS = {"omat_pbe": "PBE(+U)", "oc20_usemppbe": "PBE slabs", "matpes_r2scan": "r²SCAN"}
IN_PLANE = [True, True, False, False, False, True]     # FrechetCellFilter mask: a, b, γ relax; c fixed
EV_PER_A2_TO_J_PER_M2 = 16.021766


@functools.lru_cache(maxsize=2)
def calculator(model, head, dispersion):
    """MACE in float64. D3(BJ) is added here, not through mace_mp, which would pass `head` to D3 too."""
    from mace.calculators import mace_mp
    calc = mace_mp(model=model, device="cpu", default_dtype="float64", **({"head": head} if model.startswith("mh") else {}))
    if not dispersion:
        return calc
    import torch
    from ase.units import Bohr
    from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator
    d3 = TorchDFTD3Calculator(device="cpu", damping="bj", xc="r2scan" if head == "matpes_r2scan" else "pbe",
                              dtype=torch.float64, cutoff=40.0 * Bohr)
    return SumCalculator([calc, d3])


def version(package):
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return None


def relax(inputs, progress=None, should_stop=None):
    """Build, then relax. Returns (relaxed Structure, outputs of builder.build plus "relaxation" and POSCAR_relaxed)."""
    p = {**DEFAULTS, **inputs}
    structure, outputs = builder.build(inputs)
    outputs["inputs"].update({k: p[k] for k in DEFAULTS})

    layers = np.array(structure.site_properties["layer"])
    n_mx = len(outputs["mxene"]["layers"])
    n_all = n_mx + outputs["metal"]["layers"]
    atoms = Atoms([s.species_string for s in structure], positions=structure.cart_coords,
                  cell=structure.lattice.matrix, pbc=True)
    fixed = np.flatnonzero(layers >= max(n_mx, n_all - p["fixed_layers"]))      # top metal layers only
    atoms.set_constraint(FixAtoms(indices=fixed))
    atoms.calc = calculator(p["model"], p["head"], p["dispersion"])

    target = FrechetCellFilter(atoms, mask=IN_PLANE) if p["relax_cell"] else atoms
    opt = FIRE(target, logfile=None)
    trace, converged, stopped = [], False, False
    for converged in opt.irun(fmax=p["fmax"], steps=p["max_steps"]):
        forces = target.get_forces()
        trace.append(dict(step=opt.nsteps, energy=float(atoms.get_potential_energy()),
                          fmax=float(np.sqrt((forces**2).sum(axis=1).max()))))
        if progress:
            progress(trace[-1])
        if should_stop and should_stop():
            stopped = True
            break

    z = atoms.positions[:, 2]
    a, b, _, _, _, gamma = atoms.cell.cellpar()
    result = dict(optimizer="FIRE", fixed_atoms=len(fixed), converged=bool(converged) and not stopped, stopped=stopped,
                  n_steps=opt.nsteps, energy_initial=trace[0]["energy"], energy_final=float(atoms.get_potential_energy()),
                  fmax_final=trace[-1]["fmax"], cell=dict(a=round(a, 4), b=round(b, 4), gamma=round(gamma, 3)),
                  gap=round(float(z[layers == n_mx].mean() - z[layers == n_mx - 1].mean()), 4),
                  versions={pkg: version(pkg) for pkg in ("mace-torch", "torch", "torch-dftd", "ase")}, trace=trace)

    if p["work_of_adhesion"] and not stopped:
        metal = layers >= n_mx
        energies = []
        for mask in (metal, ~metal):
            part = atoms[mask]
            part.set_constraint()
            part.cell[2, 2] += 20.0                    # extra vacuum keeps each part isolated
            part.calc = atoms.calc
            energies.append(float(part.get_potential_energy()))
        area = float(np.linalg.norm(np.cross(atoms.cell[0], atoms.cell[1])))
        w = (sum(energies) - result["energy_final"]) / area
        result.update(energy_metal=energies[0], energy_mxene=energies[1], area=round(area, 4),
                      work_of_adhesion=round(w * EV_PER_A2_TO_J_PER_M2, 5))

    relaxed = Structure(atoms.cell[:], atoms.get_chemical_symbols(), atoms.positions, coords_are_cartesian=True,
                        to_unit_cell=False, site_properties=structure.site_properties)
    outputs["relaxation"] = result
    outputs["comparison"] = compare(structure, relaxed)
    outputs["files"]["POSCAR_relaxed"] = relaxed.get_sorted_structure().to(fmt="poscar")
    return relaxed, outputs


def compare(initial, relaxed, top=15):
    """Built vs relaxed structure (same atom order): displacements, layer heights and spacings, gap and cell."""
    layers = np.array(initial.site_properties["layer"])
    parts = np.array(initial.site_properties["part"])
    frac = relaxed.lattice.get_fractional_coords(relaxed.cart_coords - initial.cart_coords)
    frac[:, :2] -= np.round(frac[:, :2])               # minimum image in-plane
    d = relaxed.lattice.get_cartesian_coords(frac)
    dr = np.linalg.norm(d, axis=1)
    z0, z1 = initial.cart_coords[:, 2], relaxed.cart_coords[:, 2]
    r = lambda v, n=4: round(float(v), n)

    rows, zm = [], []                                  # zm: mean z per layer, built and relaxed
    for i in range(layers.max() + 1):
        on = layers == i
        zm.append((z0[on].mean(), z1[on].mean()))
        rows.append(dict(index=i, part=str(parts[on][0]), z_initial=r(zm[i][0]), z_relaxed=r(zm[i][1]),
                         dz=r(zm[i][1] - zm[i][0]),
                         spacing_initial=r(zm[i][0] - zm[i - 1][0]) if i else None,
                         spacing_relaxed=r(zm[i][1] - zm[i - 1][1]) if i else None,
                         rumpling=r(np.ptp(z1[on])),                          # z range within the relaxed layer
                         slide=r(np.linalg.norm(d[on, :2].mean(axis=0))),     # mean in-plane displacement
                         mean_displacement=r(dr[on].mean())))

    tagged = relaxed.copy(site_properties={"atom": list(range(len(relaxed)))})
    poscar = np.empty(len(relaxed), dtype=int)                                # atom number in POSCAR_relaxed, from 1
    poscar[tagged.get_sorted_structure().site_properties["atom"]] = np.arange(1, len(relaxed) + 1)
    largest = [dict(atom=int(poscar[k]), species=relaxed[k].species_string, layer=int(layers[k]), part=str(parts[k]),
                    displacement=[r(v) for v in d[k]], norm=r(dr[k])) for k in np.argsort(-dr, kind="stable")[:top]]

    n_mx = int(layers[parts == "mxene"].max()) + 1    # first metal layer
    cell = lambda s: dict(a=r(s.lattice.a), b=r(s.lattice.b), gamma=r(s.lattice.gamma, 3))
    return dict(
        max_displacement=largest[0],
        rms_displacement=r(np.sqrt((dr**2).mean())),
        max_displacement_part={p: r(dr[parts == p].max()) for p in ("mxene", "metal")},
        gap=dict(initial=r(zm[n_mx][0] - zm[n_mx - 1][0]), relaxed=r(zm[n_mx][1] - zm[n_mx - 1][1])),
        mxene_thickness=dict(initial=r(zm[n_mx - 1][0] - zm[0][0]), relaxed=r(zm[n_mx - 1][1] - zm[0][1])),
        cell=dict(initial=cell(initial), relaxed=cell(relaxed)),
        layers=rows,
        largest=largest,
    )


if __name__ == "__main__":
    inputs = json.loads(Path(sys.argv[1]).read_text()) if len(sys.argv) > 1 else {}
    _, outputs = relax(inputs, progress=lambda s: print(f"{s['step']:5d}  {s['energy']:.6f} eV  {s['fmax']:.4f} eV/Å"))
    for name, text in outputs["files"].items():
        Path(name).write_text(text)
    Path("outputs.json").write_text(json.dumps(outputs, indent=2) + "\n")
    c = outputs["comparison"]
    print(f"max |Δr| {c['max_displacement']['norm']:.3f} Å ({c['max_displacement']['species']} #{c['max_displacement']['atom']}), "
          f"RMS {c['rms_displacement']:.3f} Å, gap {c['gap']['initial']:.3f} → {c['gap']['relaxed']:.3f} Å")
