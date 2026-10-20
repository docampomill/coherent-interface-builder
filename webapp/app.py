"""Flask front end. The science lives in mxbuilder; this file only moves JSON and files.

    flask --app webapp.app run            (development)
    gunicorn -w 1 --threads 4 -b 0.0.0.0:8000 webapp.app:app
"""
import io
import json
import os
import re
import zipfile
from pathlib import Path

import numpy as np
from ase.io import read
from flask import Flask, abort, jsonify, render_template, request, send_file, send_from_directory

from mxbuilder import __version__, analysis, export, lattice
from mxbuilder.build import build
from mxbuilder.elements import FCC_METALS, PRESETS, TRANSITION_METALS, X_ELEMENTS
from mxbuilder.schema import HEADS, Inputs
from mxbuilder.stack import resolve

from .jobs import Busy, Jobs

RUNS = Path(os.environ.get("MXB_RUNS", "runs")).resolve()
RUN_ID = re.compile(r"(scan-)?[0-9a-f]{12}")

app = Flask(__name__)
jobs = Jobs(RUNS)


def posted_inputs():
    return Inputs.from_dict(request.get_json()["inputs"])


def structure(atoms):
    return dict(symbols=atoms.get_chemical_symbols(), positions=np.round(atoms.positions, 5).tolist(),
                cell=np.round(atoms.cell.array, 5).tolist(), tags=atoms.get_tags().tolist())


def run_folder(run_id):
    if not RUN_ID.fullmatch(run_id) or not (RUNS / run_id).is_dir():
        abort(404)
    return RUNS / run_id


@app.errorhandler(ValueError)
def bad_inputs(err):
    return jsonify(error=str(err)), 400


@app.errorhandler(Busy)
def busy(err):
    return jsonify(error=str(err)), 409


@app.get("/")
def index():
    options = dict(presets=PRESETS, metals=TRANSITION_METALS, greek=X_ELEMENTS, fcc=list(FCC_METALS),
                   heads=HEADS, defaults=Inputs().to_dict())
    return render_template("index.html", options=options, version=__version__)


@app.post("/api/build")
def api_build():
    inputs = posted_inputs()
    iface = build(inputs)
    periodic = iface.inputs.geometry == "periodic"
    return jsonify(run_id=inputs.run_id(), resolved=iface.inputs.to_dict(), summary=iface.summary(),
                   warnings=iface.warnings, layers=analysis.layer_rows(iface.atoms, iface.layers),
                   gaps=analysis.gaps(iface.atoms, iface.layers, periodic), structure=structure(iface.atoms))


@app.post("/api/match")
def api_match():
    body = request.get_json()
    p, _ = resolve(Inputs.from_dict(body["inputs"]))
    per_metal = p.metal_layers * (2 if p.geometry == "sandwich" else 1)
    found = lattice.search(p.mxene_a, p.metal_a, len(p.mxene_layers), per_metal,
                           float(body["max_strain"]) / 100, int(body["max_atoms"]),
                           body.get("twist"), float(body.get("twist_tol", 1.0)))
    return jsonify(candidates=found)


@app.post("/api/export/<fmt>")
def api_export(fmt):
    if fmt not in export.FORMATS:
        abort(404)
    iface = build(posted_inputs())
    data = io.BytesIO(export.to_text(iface.atoms, fmt).encode())
    return send_file(data, as_attachment=True, download_name=f"initial.{export.FORMATS[fmt][0]}", mimetype="text/plain")


@app.post("/api/run")
def api_run():
    return jsonify(job=jobs.start("run", posted_inputs()))


@app.post("/api/scan")
def api_scan():
    body = request.get_json()
    grid = min(max(int(body.get("grid", 6)), 2), 12)
    return jsonify(job=jobs.start("scan", posted_inputs(), grid=grid, relax_z=bool(body.get("relax_z"))))


@app.get("/api/jobs/<job_id>")
def api_job(job_id):
    job = jobs.get(job_id)
    since = int(request.args.get("since", 0))
    if job is None:                       # finished in an earlier session: read it from disk
        folder = run_folder(job_id)
        done = folder / ("scan.json" if job_id.startswith("scan-") else "outputs.json")
        if not done.exists():
            abort(404)
        reply = dict(id=job_id, kind="scan" if job_id.startswith("scan-") else "run", state="done",
                     events=[], n_events=0, result=json.loads(done.read_text()))
    else:
        reply = dict(id=job.id, kind=job.kind, state=job.state, events=job.events[since:],
                     n_events=len(job.events), error=job.error, result=job.result)
    relaxed = RUNS / job_id / "relaxed.extxyz"
    if reply["state"] in ("done", "stopped") and relaxed.exists():
        reply["structure"] = structure(read(relaxed))
    return jsonify(reply)


@app.post("/api/jobs/<job_id>/stop")
def api_stop(job_id):
    job = jobs.get(job_id) or abort(404)
    job.stop.set()
    return jsonify(ok=True)


@app.get("/api/runs/<run_id>/files/<name>")
def api_file(run_id, name):
    return send_from_directory(run_folder(run_id), name, as_attachment=True)


@app.get("/api/runs/<run_id>/archive")
def api_archive(run_id):
    folder = run_folder(run_id)
    data = io.BytesIO()
    with zipfile.ZipFile(data, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(folder.iterdir()):
            zf.write(path, f"{run_id}/{path.name}")
    data.seek(0)
    return send_file(data, as_attachment=True, download_name=f"{run_id}.zip", mimetype="application/zip")
