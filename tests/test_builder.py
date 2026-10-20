import json

import numpy as np
import pytest
from ase.calculators.lj import LennardJones

from mxbuilder import export, lattice, relax, workflow
from mxbuilder.build import build, layer_label
from mxbuilder.schema import Inputs
from mxbuilder.stack import resolve


def sequence(iface):
    return " ".join(layer_label(layer) for layer in iface.layers)


def test_matches_yi_builder():
    """Same spacings, vacuum and stacking (mirror labels) as Yi's Au5/Ti3C2O2 run of 2026-10-16."""
    iface = build(Inputs(mxene_a=3.07))
    z = np.array([layer.z for layer in iface.layers])
    yi = [0.926481, 1.18195, 1.18195, 1.18195, 1.18195, 0.926481, 2.3, 2.115247, 2.115247, 2.115247, 2.115247]
    assert np.allclose(np.diff(z), yi, atol=1e-5)
    assert iface.atoms.cell[2, 2] == pytest.approx(32.34175, abs=1e-4)
    # Yi's file reads A C B A C B A | C B A C B; ours is the same with A and C swapped (a mirror)
    swap = str.maketrans("AC", "CA")
    ours = "".join(layer.site for layer in iface.layers)
    assert ours.translate(swap) == "ACBACBACBACB"


def test_sites_from_stacking():
    assert sequence(build(Inputs(metal_layers=0))) == "O(C) Ti(A) C(β) Ti(C) C(α) Ti(B) O(C)"
    assert sequence(build(Inputs(metal_layers=0, stacking="ABA"))) == "O(C) Ti(A) C(β) Ti(A) C(β) Ti(A) O(C)"
    assert sequence(build(Inputs(metal_layers=0, mxene_layers=["O", "V", "C", "V", "O"]))) == "O(C) V(A) C(β) V(C) O(A)"


def test_per_layer_overrides():
    p = Inputs(mxene_sites=[None, None, None, None, None, "C", None], mxene_dz=[0, None, 1.3, None, None, None, None])
    r, _ = resolve(p)
    assert r.mxene_sites[5] == "C" and r.mxene_dz[2] == 1.3 and r.mxene_dz[1] == pytest.approx(0.926481, abs=1e-6)


def test_deterministic_and_reproducible():
    a, b = build(Inputs()), build(Inputs())
    assert export.to_text(a.atoms, "poscar") == export.to_text(b.atoms, "poscar")
    assert Inputs().run_id() == Inputs().run_id() != Inputs(gap=2.4).run_id()
    again = build(Inputs.from_dict(json.loads(json.dumps(a.inputs.to_dict()))))   # resolved.json fed back
    assert export.to_text(again.atoms, "extxyz") == export.to_text(a.atoms, "extxyz")


def test_fixed_layers_become_selective_dynamics():
    iface = build(Inputs())
    assert iface.inputs.fixed_layers == [10, 11]
    flags = [line.split()[3] for line in export.to_text(iface.atoms, "poscar").splitlines()[9:]]
    assert flags.count("F") == 2


def test_sandwich_is_inversion_symmetric():
    iface = build(Inputs(geometry="sandwich"))
    atoms = iface.atoms
    middle = atoms.positions[atoms.get_tags() == 8][0]          # middle Ti of Ti3C2O2
    frac = atoms.cell.scaled_positions(2 * middle - atoms.positions)
    frac[:, :2] %= 1.0
    original = np.round(atoms.get_scaled_positions(wrap=True), 6) % 1.0
    inverted = np.round(frac, 6) % 1.0
    key = lambda arr: sorted(map(tuple, arr))
    assert key(original) == key(inverted)
    assert iface.inputs.fixed_layers == [0, 1, 15, 16]


def test_periodic_gap_across_the_boundary():
    from mxbuilder import analysis
    iface = build(Inputs(geometry="periodic", gap=2.5))
    gaps = analysis.gaps(iface.atoms, iface.layers, periodic=True)
    assert gaps == {"top": 2.5, "periodic": 2.5}


def test_coincidence_cell():
    found = lattice.search(3.07, 4.078, 7, 5, max_strain=0.02, max_atoms=400)
    c = found[0]
    iface = build(Inputs(mxene_supercell=c["mxene"], metal_supercell=c["metal"]))
    s = iface.summary()
    assert s["n_atoms"] == c["n_atoms"]
    assert s["strain"] == pytest.approx(c["strain"], abs=1e-3)
    assert s["twist"] == pytest.approx(c["twist"], abs=1e-3)
    assert all(layer.n_atoms in (c["n_mx"], c["n_metal"]) for layer in iface.layers)


def test_validation():
    with pytest.raises(ValueError):
        build(Inputs(mxene_layers=["O", "Xx", "O"]))
    with pytest.raises(ValueError):
        build(Inputs(mxene_sites=["A"]))
    with pytest.raises(ValueError):
        Inputs.from_dict({"not_an_input": 1})


@pytest.fixture
def lj(monkeypatch):
    monkeypatch.setattr(relax, "calculator", lambda *a, **k: LennardJones(sigma=2.2, epsilon=0.02, rc=6.0))
    monkeypatch.setattr(relax, "model_record", lambda model: dict(name=model, note="test stand-in"))


def test_workflow_run_and_cache(tmp_path, lj):
    p = Inputs(max_steps=5)
    out = workflow.run(p, tmp_path)
    folder = tmp_path / p.run_id()
    for name in ["inputs.json", "resolved.json", "outputs.json", "environment.json", "initial.vasp", "relaxed.lmp", "relax.traj"]:
        assert (folder / name).exists(), name
    assert out["n_steps"] <= 5 and "work_of_adhesion" in out and len(out["layers"]) == 12
    assert workflow.run(p, tmp_path) == json.loads((folder / "outputs.json").read_text())


def test_scan_and_relax_z_keep_xy(tmp_path, lj):
    out = workflow.scan(Inputs(), tmp_path, grid=2, relax_z=True)
    assert len(out["points"]) == 4 and min(pt["delta"] for pt in out["points"]) == 0
    atoms = build(Inputs()).atoms
    xy = atoms.positions[:, :2].copy()
    relax.single_point(atoms, LennardJones(sigma=2.2, epsilon=0.02, rc=6.0), relax_z=True)
    assert np.allclose(atoms.positions[:, :2], xy)
