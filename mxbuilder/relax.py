"""MACE-MH-1 calculator, relaxation and the energies derived from it."""
import functools
import hashlib
import math
from dataclasses import replace
from pathlib import Path

import numpy as np
from ase.build import bulk
from ase.constraints import FixAtoms, FixCartesian
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS, FIRE

EV_PER_A2_TO_J_PER_M2 = 16.021766
IN_PLANE = [True, True, False, False, False, True]    # a, b, gamma relax; c fixed


@functools.lru_cache(maxsize=4)
def calculator(model="mh-1", head="omat_pbe", dispersion=False, device="cpu"):
    """MACE in float64; D3(BJ) is added here rather than through mace_mp, which would pass `head` to D3 too."""
    from mace.calculators import mace_mp
    calc = mace_mp(model=model, head=head, device=device, default_dtype="float64")
    if not dispersion:
        return calc
    import torch
    from ase.calculators.mixing import SumCalculator
    from ase.units import Bohr
    from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator
    xc = "r2scan" if head == "matpes_r2scan" else "pbe"
    d3 = TorchDFTD3Calculator(device=device, damping="bj", xc=xc, dtype=torch.float64, cutoff=40.0 * Bohr)
    return SumCalculator([calc, d3])


@functools.lru_cache(maxsize=4)
def model_record(model):
    """Where the model file lives and its SHA-256, for environment.json."""
    try:
        from mace.calculators.foundations_models import download_mace_mp_checkpoint
        path = Path(download_mace_mp_checkpoint(model))
        return dict(name=model, file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    except Exception as exc:    # the record is informative only
        return dict(name=model, error=str(exc))


def max_force(forces):
    return float(np.sqrt((forces ** 2).sum(axis=1).max()))


def relax(atoms, calc, fmax=0.05, max_steps=500, optimizer="FIRE", relax_cell=False,
          logfile=None, trajectory=None, on_step=None, should_stop=None):
    """Relax in place. Returns converged, stopped, n_steps and the energy/fmax trace."""
    atoms.calc = calc
    target = FrechetCellFilter(atoms, mask=IN_PLANE) if relax_cell else atoms
    opt = {"FIRE": FIRE, "BFGS": BFGS}[optimizer](target, logfile=logfile, trajectory=trajectory)
    trace, converged, stopped = [], False, False
    for converged in opt.irun(fmax=fmax, steps=max_steps):
        trace.append(dict(step=opt.nsteps, energy=float(atoms.get_potential_energy()),
                          fmax=max_force(target.get_forces())))
        if on_step:
            on_step(trace[-1])
        if should_stop and should_stop():
            stopped = True
            break
    return dict(converged=bool(converged) and not stopped, stopped=stopped, n_steps=opt.nsteps, trace=trace)


def separated_energies(atoms, calc, metal_mask):
    """Energies of the metal and the MXene alone, frozen at the interface geometry."""
    energies = []
    for mask in (metal_mask, ~metal_mask):
        part = atoms[mask]
        part.set_constraint()
        part.cell[2, 2] += 20.0          # extra vacuum keeps the part isolated, also for periodic cells
        part.calc = calc
        energies.append(float(part.get_potential_energy()))
    return energies


def single_point(atoms, calc, relax_z=False):
    """Energy at fixed geometry, or after relaxing z only (x and y held)."""
    atoms.calc = calc
    if relax_z:
        everyone = range(len(atoms))
        atoms.set_constraint(atoms.constraints + [FixCartesian(everyone, mask=(True, True, False))])
        FIRE(atoms, logfile=None).run(fmax=0.05, steps=200)
    return float(atoms.get_potential_energy())


def mlip_lattice(p, calc):
    """MLIP equilibrium a of the free MXene (1×1, in-plane cell relaxed) and of the bulk fcc metal."""
    from .build import build
    from .elements import FCC_METALS

    metal = bulk(p.metal, "fcc", a=p.metal_a or FCC_METALS[p.metal][0])
    metal.calc = calc
    FIRE(FrechetCellFilter(metal, hydrostatic_strain=True), logfile=None).run(fmax=1e-3, steps=300)
    a_metal = metal.cell.lengths()[0] * math.sqrt(2)

    free = replace(p, metal_layers=0, geometry="slab", mxene_supercell=[1, 0], metal_supercell=[1, 0],
                   repeat=[1, 1], rbt=[0.0, 0.0], fixed_layers=[])
    mx = build(free).atoms
    mx.calc = calc
    FIRE(FrechetCellFilter(mx, mask=IN_PLANE), logfile=None).run(fmax=1e-3, steps=500)
    return round(float(mx.cell.lengths()[0]), 5), round(float(a_metal), 5)


def fixed_atoms(atoms):
    return [int(i) for c in atoms.constraints if isinstance(c, FixAtoms) for i in c.index]
