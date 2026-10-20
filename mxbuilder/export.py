"""Structure files. Fixed layers become selective dynamics (POSCAR) and
move_mask (extxyz); layer indices travel as tags in extxyz."""
import io
from pathlib import Path

from ase.io import write

FORMATS = {               # key: (extension, ASE format)
    "poscar": ("vasp", "vasp"),
    "extxyz": ("extxyz", "extxyz"),
    "cif": ("cif", "cif"),
    "lammps": ("lmp", "lammps-data"),
}


def to_text(atoms, fmt):
    buf = io.StringIO()
    if fmt == "poscar":
        write(buf, atoms, format="vasp", direct=True)
    elif fmt == "lammps":
        species = list(dict.fromkeys(atoms.get_chemical_symbols()))
        write(buf, atoms, format="lammps-data", specorder=species, masses=True, atom_style="atomic")
    elif fmt == "cif":
        raw = io.BytesIO()
        write(raw, atoms, format="cif")
        return raw.getvalue().decode()
    else:
        write(buf, atoms, format="extxyz")
    return buf.getvalue()


def write_all(atoms, folder, stem):
    """Write every format as <stem>.<ext>; returns {format: filename}."""
    names = {}
    for fmt, (ext, _) in FORMATS.items():
        name = f"{stem}.{ext}"
        Path(folder, name).write_text(to_text(atoms, fmt))
        names[fmt] = name
    return names
