"""Layer-resolved numbers: heights, spacings, gaps and what changed on relaxation.

Heights are layer means (atoms grouped by tag). Gaps are differences of layer
means across an MXene/metal boundary.
"""
import numpy as np
from ase.constraints import FixAtoms
from ase.geometry import find_mic

from .stack import label


def layer_z(atoms, n_layers):
    z, tags = atoms.positions[:, 2], atoms.get_tags()
    return np.array([z[tags == i].mean() for i in range(n_layers)])


def layer_rows(atoms, layers):
    """One dict per layer, bottom → top."""
    z, tags = atoms.positions[:, 2], atoms.get_tags()
    zmean = layer_z(atoms, len(layers))
    fixed = {int(i) for c in atoms.constraints if isinstance(c, FixAtoms) for i in c.index}
    rows = []
    for i, layer in enumerate(layers):
        mine = np.flatnonzero(tags == i)
        rows.append(dict(index=i, element=layer.element, site=layer.site, label=label(layer.element, layer.site),
                         part=layer.part, n_atoms=len(mine), fixed=bool(fixed & set(mine.tolist())),
                         z=round(float(zmean[i]), 4),
                         spacing=round(float(zmean[i] - zmean[i - 1]), 4) if i else None,
                         rumpling=round(float(np.ptp(z[mine])), 4)))
    return rows


def gaps(atoms, layers, periodic=False):
    """MXene–metal gaps: 'top' (metal above), 'bottom' (metal below), and 'periodic' across c."""
    zmean = layer_z(atoms, len(layers))
    out = {}
    for i in range(1, len(layers)):
        lower, upper = layers[i - 1].part, layers[i].part
        if lower != upper:
            out["top" if upper == "metal" else "bottom"] = round(float(zmean[i] - zmean[i - 1]), 4)
    if periodic and layers[0].part != layers[-1].part:
        across = atoms.cell[2, 2] + zmean[0] - zmean[-1]
        out["periodic"] = round(float(across), 4)
    return out


def compare(initial, relaxed, layers):
    """Add relaxed heights and per-layer displacements to the layer rows."""
    rows = layer_rows(initial, layers)
    z_after = layer_z(relaxed, len(layers))
    disp, _ = find_mic(relaxed.positions - initial.positions, relaxed.cell, pbc=[True, True, False])
    tags = initial.get_tags()
    for row in rows:
        i = row["index"]
        d = disp[tags == i]
        row["z_relaxed"] = round(float(z_after[i]), 4)
        row["spacing_relaxed"] = round(float(z_after[i] - z_after[i - 1]), 4) if i else None
        row["slide"] = round(float(np.linalg.norm(d[:, :2].mean(axis=0))), 4)
        row["mean_displacement"] = round(float(np.linalg.norm(d, axis=1).mean()), 4)
    return rows
