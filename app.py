"""Flask front end for builder.py and relax.py.

    flask --app app run --port 8000
    gunicorn -w 1 --threads 4 -b 0.0.0.0:8000 app:app      (one worker: the relaxation job lives in it)
"""
import threading
from collections import defaultdict

import numpy as np
from flask import Flask, jsonify, render_template, request

import builder
import planes
import relax

app = Flask(__name__)
job = dict(id=0, state="idle", events=[], result=None, error=None, stop=threading.Event())    # one relaxation at a time
lock = threading.Lock()


def structure_json(structure):
    return dict(symbols=[s.species_string for s in structure], positions=np.round(structure.cart_coords, 5).tolist(),
                cell=np.round(structure.lattice.matrix, 5).tolist(),
                layers=structure.site_properties["layer"], parts=structure.site_properties["part"])


@app.errorhandler(ValueError)
def bad_inputs(err):
    return jsonify(error=str(err)), 400


@app.get("/")
def index():
    mxenes, metals = defaultdict(list), defaultdict(list)
    for row in builder.MXENES.values():
        mxenes[row["family"]].append(row)
    for row in builder.METALS.values():
        metals[row["structure"]].append(row)
    options = dict(mxenes=mxenes, metals=metals, heads=relax.HEADS, defaults={**builder.DEFAULTS, **relax.DEFAULTS},
                   spreadsheet=dict(file=builder.SPREADSHEET.name, sha256=builder.SHA256))
    return render_template("index.html", options=options)


@app.post("/api/build")
def api_build():
    structure, outputs = builder.build(request.get_json()["inputs"])
    return jsonify(outputs=outputs, structure=structure_json(structure))


@app.post("/api/planes")
def api_planes():
    body = request.get_json()
    structure, outputs = builder.build(body["inputs"])
    cell = outputs["mxene"]["matrix"] if body["part"] == "mxene" else outputs["metal"]["matrix"]
    result = planes.atomic_planes(structure, body["part"], cell, planes.miller(body["hkl"]))
    result["polygons"] = planes.plane_polygons(structure, result, [int(v) for v in body.get("repeat", [1, 1])])
    return jsonify(result)


def run_relaxation(j, inputs):
    try:
        structure, outputs = relax.relax(inputs, progress=j["events"].append, should_stop=j["stop"].is_set)
        initial, _ = builder.build(inputs)             # same atom order as the relaxed structure
        j["result"] = dict(outputs=outputs, structure=structure_json(structure), initial=structure_json(initial))
        j["state"] = "stopped" if outputs["relaxation"]["stopped"] else "done"
    except Exception as exc:
        j["error"], j["state"] = f"{type(exc).__name__}: {exc}", "error"


@app.post("/api/relax")
def api_relax():
    global job
    with lock:
        if job["state"] == "running":
            return jsonify(error="A relaxation is running. Stop it first."), 409
        job = dict(id=job["id"] + 1, state="running", events=[], result=None, error=None, stop=threading.Event())
    threading.Thread(target=run_relaxation, args=(job, request.get_json()["inputs"]), daemon=True).start()
    return jsonify(id=job["id"])


@app.get("/api/relax")
def api_relax_status():
    since = int(request.args.get("since", 0))
    return jsonify(id=job["id"], state=job["state"], events=job["events"][since:], n_events=len(job["events"]),
                   error=job["error"], result=job["result"])


@app.post("/api/relax/stop")
def api_relax_stop():
    job["stop"].set()
    return jsonify(ok=True)
