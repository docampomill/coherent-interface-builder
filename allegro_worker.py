"""Allegro (nequip) energies, forces and stresses for relax.py, in the Python that has nequip.

    /opt/allegro/bin/python allegro_worker.py MODEL.nequip.zip

nequip needs e3nn 0.6 and MACE's model loads only with e3nn 0.4.4, so Allegro runs here, in its own environment, and
relax.AllegroCalculator talks to it. One JSON object per line each way. The first reply says the model is loaded:
    {"ready": true, "versions": {"nequip": ..., "nequip-allegro": ..., "e3nn": ..., "torch": ...}}
then, for each request {"numbers": [...], "positions": [[x, y, z], ...], "cell": [[...], ...], "pbc": [...]}:
    {"energy": E, "forces": [[...], ...], "stress": [6 Voigt components, eV/Å³]}  or  {"error": "..."}
Anything nequip or torch prints goes to stderr, so stdout carries only the replies.
"""
import json
import os
import sys
import traceback


def main():
    reply = os.fdopen(os.dup(sys.stdout.fileno()), "w", buffering=1)   # the protocol's stdout
    os.dup2(sys.stderr.fileno(), sys.stdout.fileno())                  # everything else to stderr
    sys.stdout = sys.stderr

    try:
        import numpy as np
        from ase import Atoms
        try:                                                           # nequip >= 0.17
            from nequip.integrations.ase import NequIPCalculator
        except ImportError:                                            # nequip 0.16
            from nequip.ase import NequIPCalculator
        model = sys.argv[1]
        if model.endswith(".zip"):                                     # the package, as the handoff loads it
            calc = NequIPCalculator._from_saved_model(model, device="cpu", chemical_species_to_atom_type_map=True)
        else:                                                          # a TorchScript compile
            calc = NequIPCalculator.from_compiled_model(model, device="cpu", chemical_species_to_atom_type_map=True)
        from importlib import metadata
        versions = {}
        for package in ("nequip", "nequip-allegro", "e3nn", "torch"):
            try:
                versions[package] = metadata.version(package)
            except metadata.PackageNotFoundError:
                versions[package] = None
    except Exception:
        reply.write(json.dumps({"error": traceback.format_exc()}) + "\n")
        return
    reply.write(json.dumps({"ready": True, "versions": versions}) + "\n")

    for line in sys.stdin:
        try:
            r = json.loads(line)
            atoms = Atoms(numbers=r["numbers"], positions=r["positions"], cell=r["cell"], pbc=r["pbc"])
            atoms.calc = calc
            out = dict(energy=float(atoms.get_potential_energy()), forces=atoms.get_forces().tolist())
            if all(r["pbc"]):
                out["stress"] = np.asarray(atoms.get_stress(voigt=True)).tolist()
        except Exception:
            out = {"error": traceback.format_exc()}
        reply.write(json.dumps(out) + "\n")


if __name__ == "__main__":
    main()
