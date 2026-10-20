"""Command line: the same functions the web app calls.

    python -m mxbuilder inputs                  inputs table (markdown)
    python -m mxbuilder outputs                 outputs table (markdown)
    python -m mxbuilder build inputs.json -o out    structure files only, no MLIP
    python -m mxbuilder run inputs.json         build + relax, files in runs/<run_id>/
    python -m mxbuilder scan inputs.json --grid 6 [--relax-z]
    python -m mxbuilder match inputs.json --max-strain 3 --max-atoms 400 [--twist 0]
"""
import argparse
import json
from pathlib import Path

from . import export, lattice, workflow
from .build import build
from .schema import OUTPUTS, Inputs, describe_inputs
from .stack import resolve


def load(path):
    return Inputs.from_dict(json.loads(Path(path).read_text())) if path else Inputs()


def markdown(rows, columns):
    lines = ["| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
    cell = lambda v: "" if v is None or v == "" else str(v).replace("|", "\\|")
    lines += ["| " + " | ".join(cell(r[c]) for c in columns) + " |" for r in rows]
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="mxbuilder", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("inputs")
    sub.add_parser("outputs")
    b = sub.add_parser("build")
    b.add_argument("inputs", nargs="?")
    b.add_argument("-o", "--out", default=".")
    r = sub.add_parser("run")
    r.add_argument("inputs", nargs="?")
    r.add_argument("--runs", default="runs")
    r.add_argument("--force", action="store_true")
    s = sub.add_parser("scan")
    s.add_argument("inputs", nargs="?")
    s.add_argument("--runs", default="runs")
    s.add_argument("--grid", type=int, default=6)
    s.add_argument("--relax-z", action="store_true")
    m = sub.add_parser("match")
    m.add_argument("inputs", nargs="?")
    m.add_argument("--max-strain", type=float, default=3.0, help="percent")
    m.add_argument("--max-atoms", type=int, default=400)
    m.add_argument("--twist", type=float, default=None, help="degrees; omit for any")
    m.add_argument("--twist-tol", type=float, default=1.0)
    args = ap.parse_args(argv)

    if args.cmd == "inputs":
        print(markdown(describe_inputs(), ["name", "group", "type", "default", "units", "allowed", "description"]))
    elif args.cmd == "outputs":
        rows = [dict(name=k, type=t, units=u, description=d) for k, (t, u, d) in OUTPUTS.items()]
        print(markdown(rows, ["name", "type", "units", "description"]))
    elif args.cmd == "build":
        iface = build(load(args.inputs))
        Path(args.out).mkdir(parents=True, exist_ok=True)
        names = export.write_all(iface.atoms, args.out, "initial")
        print(iface.formula, f"{len(iface.atoms)} atoms ->", ", ".join(names.values()))
        for note in iface.warnings:
            print("note:", note)
    elif args.cmd == "run":
        progress = lambda e: print(f"  step {e['step']:4d}  E = {e['energy']:.5f} eV  fmax = {e['fmax']:.4f}") \
            if "step" in e else print(e["message"])
        out = workflow.run(load(args.inputs), args.runs, progress=progress, force=args.force)
        keys = ["run_id", "formula", "n_atoms", "strain", "energy_final", "converged", "gaps_relaxed", "work_of_adhesion"]
        print(json.dumps({k: out.get(k) for k in keys}, indent=2))
    elif args.cmd == "scan":
        out = workflow.scan(load(args.inputs), args.runs, args.grid, args.relax_z,
                            progress=lambda e: print(e["message"]))
        print("lowest energy at rbt =", out["best"])
    elif args.cmd == "match":
        p, _ = resolve(load(args.inputs))
        per_metal = p.metal_layers * (2 if p.geometry == "sandwich" else 1)
        found = lattice.search(p.mxene_a, p.metal_a, len(p.mxene_layers), per_metal,
                               args.max_strain / 100, args.max_atoms, args.twist, args.twist_tol)
        print(markdown(found, ["mxene", "metal", "strain", "twist", "n_atoms"]))
