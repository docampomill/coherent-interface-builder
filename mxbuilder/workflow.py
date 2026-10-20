"""The workflow: inputs → files and outputs.json.

The web app, the CLI and the papermill notebook all call `run` (and `scan`),
so every route executes the same code. Each run lives in runs/<run_id>/:

    inputs.json       the inputs as given (null = auto)
    resolved.json     every value actually used; feeding it back gives the same structure
    initial.*         POSCAR (.vasp), .extxyz, .cif, .lmp
    relaxed.*         same formats, after relaxation
    relax.log, relax.traj
    outputs.json      see schema.OUTPUTS
    environment.json  package versions, model file and checksum
"""
import hashlib
import json
import platform
from dataclasses import replace
from importlib import metadata
from pathlib import Path

import numpy as np

from . import __version__, analysis, export, relax
from .build import build
from .stack import resolve

PACKAGES = ("ase", "numpy", "torch", "mace-torch", "e3nn", "torch-dftd", "flask")


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + "\n")


def environment(model=None):
    env = dict(mxbuilder=__version__, python=platform.python_version(), platform=platform.platform())
    for pkg in PACKAGES:
        try:
            env[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            env[pkg] = None
    if model:
        env["model"] = relax.model_record(model)
    return env


def quiet(event):
    pass


def run(inputs, runs_dir="runs", progress=quiet, should_stop=None, force=False):
    """Build, and relax if inputs.relax. Returns the outputs dict (cached by run id)."""
    rid = inputs.run_id()
    folder = Path(runs_dir) / rid
    done = folder / "outputs.json"
    if done.exists() and not force:
        cached = json.loads(done.read_text())
        if not cached.get("stopped"):
            return cached
    folder.mkdir(parents=True, exist_ok=True)
    write_json(folder / "inputs.json", inputs.to_dict())
    out = dict(run_id=rid)

    p, calc = inputs, None
    if inputs.relax:
        progress(dict(stage="loading", message=f"Loading MACE {inputs.model} ({inputs.head})"))
        calc = relax.calculator(inputs.model, inputs.head, inputs.dispersion)
        if inputs.mlip_lattice:
            progress(dict(stage="lattice", message="Relaxing bulk metal and free MXene for lattice constants"))
            a_mx, a_metal = relax.mlip_lattice(inputs, calc)
            out["mlip_lattice"] = dict(mxene_a=a_mx, metal_a=a_metal)
            p = replace(inputs, mxene_a=inputs.mxene_a or a_mx, metal_a=inputs.metal_a or a_metal)

    iface = build(p)
    periodic = iface.inputs.geometry == "periodic"
    write_json(folder / "resolved.json", iface.inputs.to_dict())
    initial = iface.atoms.copy()
    files = dict(initial=export.write_all(initial, folder, "initial"))
    out.update(iface.summary(), warnings=iface.warnings,
               layers=analysis.layer_rows(initial, iface.layers),
               gaps_initial=analysis.gaps(initial, iface.layers, periodic))

    if inputs.relax:
        atoms = iface.atoms
        progress(dict(stage="relaxing", message="Relaxing"))
        result = relax.relax(atoms, calc, inputs.fmax, inputs.max_steps, inputs.optimizer, inputs.relax_cell,
                             logfile=str(folder / "relax.log"), trajectory=str(folder / "relax.traj"),
                             on_step=lambda s: progress(dict(stage="relaxing", **s)), should_stop=should_stop)
        files["relaxed"] = export.write_all(atoms, folder, "relaxed")
        trace = result["trace"]
        out.update(energy_initial=trace[0]["energy"], energy_final=float(atoms.get_potential_energy()),
                   converged=result["converged"], stopped=result["stopped"], n_steps=result["n_steps"],
                   fmax_final=trace[-1]["fmax"], trace=trace,
                   layers=analysis.compare(initial, atoms, iface.layers),
                   gaps_relaxed=analysis.gaps(atoms, iface.layers, periodic))
        if inputs.relax_cell:
            a, b, _, _, _, gamma = atoms.cell.cellpar()
            out["cell_relaxed"] = dict(a=round(float(a), 4), b=round(float(b), 4), gamma=round(float(gamma), 3))
        if inputs.work_of_adhesion and p.metal_layers and not result["stopped"]:
            progress(dict(stage="adhesion", message="Computing the work of adhesion"))
            metal = np.array([iface.layers[t].part == "metal" for t in atoms.get_tags()])
            e_metal, e_mxene = relax.separated_energies(atoms, calc, metal)
            area = float(np.linalg.norm(np.cross(*atoms.cell[:2])))
            w = (e_metal + e_mxene - out["energy_final"]) / (out["n_interfaces"] * area)
            out.update(energy_metal=e_metal, energy_mxene=e_mxene,
                       work_of_adhesion=round(w * relax.EV_PER_A2_TO_J_PER_M2, 5))

    out["files"] = files
    write_json(folder / "environment.json", environment(inputs.model if inputs.relax else None))
    write_json(done, out)
    progress(dict(stage="done", message="Done"))
    return out


def scan_id(inputs, grid, relax_z):
    text = inputs.canonical_json() + f"|scan|{grid}|{relax_z}|" + __version__
    return "scan-" + hashlib.sha256(text.encode()).hexdigest()[:12]


def scan(inputs, runs_dir="runs", grid=6, relax_z=False, progress=quiet, should_stop=None):
    """Energy vs rigid-body translation of the metal over one MXene primitive cell.

    The gap stays at its resolved value. relax_z lets atoms move along z only.
    """
    sid = scan_id(inputs, grid, relax_z)
    folder = Path(runs_dir) / sid
    done = folder / "scan.json"
    if done.exists():
        return json.loads(done.read_text())
    folder.mkdir(parents=True, exist_ok=True)
    write_json(folder / "inputs.json", inputs.to_dict())

    progress(dict(stage="loading", message=f"Loading MACE {inputs.model} ({inputs.head})"))
    calc = relax.calculator(inputs.model, inputs.head, inputs.dispersion)
    base, _ = resolve(inputs)
    points = []
    for i in range(grid):
        for j in range(grid):
            if should_stop and should_stop():
                return dict(scan_id=sid, stopped=True)
            shift = [round(i / grid, 6), round(j / grid, 6)]
            iface = build(replace(base, rbt=shift))
            energy = relax.single_point(iface.atoms, calc, relax_z)
            gap = analysis.gaps(iface.atoms, iface.layers).get("top")
            points.append(dict(rbt=shift, energy=energy, gap=gap))
            progress(dict(stage="scanning", message=f"Shift {len(points)} of {grid * grid}",
                          done=len(points), total=grid * grid))

    summary = iface.summary()
    e_min = min(pt["energy"] for pt in points)
    for pt in points:
        pt["delta"] = round(1000 * (pt["energy"] - e_min) / (summary["n_interfaces"] * summary["area"]), 4)
    best = min(points, key=lambda pt: pt["energy"])
    out = dict(scan_id=sid, grid=grid, relax_z=relax_z, unit="meV/Å²", best=best["rbt"], points=points)
    write_json(done, out)
    progress(dict(stage="done", message="Done"))
    return out
