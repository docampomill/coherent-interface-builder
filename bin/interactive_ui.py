"""The form of interactive.ipynb: the web app's inputs and results as ipywidgets. Every build runs the
coherentmxmetal Sim2L with simtool.Run.

    from interactive_ui import app; app()
"""
from IPython.display import HTML, display

loading = display(HTML("⏳ Loading the interface builder…"), display_id=True)   # while the imports below run

import base64
import contextlib
import html
import json
import threading
import time
import uuid
import warnings
from pathlib import Path

import ipywidgets as w
import matplotlib.pyplot as plt
from ase.data import atomic_numbers
from ase.data.colors import jmol_colors
from ase.io import read
from ase.visualize.plot import plot_atoms
from io import BytesIO, StringIO
from simtool import searchForSimTool, getSimToolInputs, Run
from simtool.experiment import get_experiment

loc = searchForSimTool("coherentmxmetal")
schema = getSimToolInputs(loc)
GREEK = {"A": "α", "B": "β", "C": "γ"}
SPIN_CSS = "@keyframes cib-spin { to { transform: rotate(360deg) } } .cib-spin { display: inline-block; width: 1em; height: 1em; border: 2px solid #ccc; border-top-color: #1f4e79; border-radius: 50%; animation: cib-spin 0.8s linear infinite; vertical-align: -0.15em; margin-right: 0.4em }"
SPIN = "<span class='cib-spin'></span>"                  # loading icon, styled by SPIN_CSS
WIDE, HALF = w.Layout(width="100%"), w.Layout(width="48%")
STYLE = {"description_width": "initial"}


def value(name):
    return getattr(schema, name).value


def choice(name, label, options=None, layout=WIDE):
    return w.Dropdown(options=options or getattr(schema, name).options, value=value(name), description=label,
                      style=STYLE, layout=layout)


def toggle(name, label, options):
    return w.ToggleButtons(options=options, value=value(name), description=label, style=STYLE)


def number(name, label, step=0.1, layout=HALF):
    kind = w.IntText if isinstance(value(name), int) else w.FloatText
    return kind(value=value(name), description=label, step=step, style=STYLE, layout=layout)


def check(name, label):
    return w.Checkbox(value=value(name), description=label, indent=False)


# ---------- inputs, in the sections of the web app ----------
f = dict(
    mxene=choice("mxene", "MXene"),
    stacking=toggle("stacking", "Stacking", [("AβC · octahedral", "ABC"), ("AβA · trigonal prismatic", "ABA")]),
    mxene_rotation=toggle("mxene_rotation", "Orientation about z", [("0° · C right of A", 0), ("180° · B ↔ C", 180)]),
    metal=choice("metal", "Metal", layout=HALF),
    metal_layers=number("metal_layers", "Layers, (111) planes", 1),
    metal_site=choice("metal_site", "First layer site", [("auto · free of the top two MXene layers", "auto"),
                      ("A · over the bottom MXene layer", "A"), ("B", "B"), ("C", "C")], HALF),
    metal_stacking=choice("metal_stacking", "Stacking", [("A → B → C · continues an MXene at 0°", "ABC"),
                          ("A → C → B · continues an MXene at 180°", "ACB")], HALF),
    gap=number("gap", "Gap (Å)"),
    shift_x=number("shift_x", "Shift x (Å)"),
    shift_y=number("shift_y", "Shift y (Å)"),
    tol=number("tol", "Strain tolerance", 0.005),
    max_replicas=number("max_replicas", "Max metal replicas", 1),
    vacuum=number("vacuum", "Vacuum below (Å)", 0.5),
    model=w.Text(value=value("model"), description="MACE model", style=STYLE, layout=HALF),
    head=choice("head", "Head", [("omat_pbe · PBE(+U)", "omat_pbe"), ("oc20_usemppbe · PBE slabs", "oc20_usemppbe"),
                ("matpes_r2scan · r²SCAN", "matpes_r2scan")], HALF),
    dispersion=check("dispersion", "Add D3(BJ) dispersion (needs torch-dftd)"),
    fmax=number("fmax", "fmax (eV/Å)", 0.01),
    max_steps=number("max_steps", "Max steps", 10),
    fixed_layers=number("fixed_layers", "Fix outermost metal layers", 1),
    relax_cell=check("relax_cell", "Relax in-plane cell (c stays fixed)"),
    work_of_adhesion=check("work_of_adhesion", "Report work of adhesion"),
)
auto = w.Checkbox(value=True, description="Rebuild on every change", indent=False)
build_button = w.Button(description="Build", button_style="primary", icon="cubes")
relax_button = w.Button(description="Relax", button_style="warning", icon="play")
status = w.HTML()
relax_status = w.HTML()
progress = w.IntProgress(min=0, max=1, layout=w.Layout(width="100%", visibility="hidden"))


def section(step, title, *children, hint=None):
    head = w.HTML(f"<h3 style='margin:0.6em 0 0.2em'><span style='display:inline-block;width:1.5em;height:1.5em;"
                  f"border-radius:50%;background:#1f4e79;color:white;text-align:center;margin-right:0.4em'>{step}</span>"
                  f"{title}</h3>")
    rows = [head, *children] + ([w.HTML(f"<small style='color:#666'>{hint}</small>")] if hint else [])
    return w.VBox(rows)


row = lambda *c: w.HBox(c, layout=w.Layout(justify_content="space-between", flex_flow="row wrap"))
form = w.VBox([
    section(1, "MXene", f["mxene"], f["stacking"], f["mxene_rotation"]),
    section(2, "Metal", row(f["metal"], f["metal_layers"]), row(f["metal_site"], f["metal_stacking"])),
    section(3, "Location", row(f["gap"], f["vacuum"]), row(f["shift_x"], f["shift_y"]), row(f["tol"], f["max_replicas"]),
            hint="Metal (111) on MXene (0001), metal [1-10] along x. The gap is from the top MXene layer to the first "
                 "metal (111) plane; the shift moves the metal from its site. The metal is the substrate and keeps its lattice; the MXene takes the in-plane strain."),
    w.HBox([build_button, auto]),
    section(4, "Relaxation", row(f["model"], f["head"]), f["dispersion"], row(f["fmax"], f["max_steps"]),
            row(f["fixed_layers"]), f["relax_cell"], f["work_of_adhesion"], relax_button, relax_status, progress,
            hint="MACE runs inside the Sim2L; a relaxation can take minutes. Runs with the same inputs come from the cache."),
], layout=w.Layout(flex="0 1 440px", min_width="360px", padding="0 1em 0 0"))

# ---------- results ----------
readout = w.HTML()
view = w.ToggleButtons(options=[("Front (x–z)", "-90x"), ("Side (y–z)", "-90y,-90z"), ("Top (x–y)", ""), ("All three", "all")],
                       value="all", style={"button_width": "auto"})
repeat = w.ToggleButtons(options=[("1×1", 1), ("2×2", 2), ("3×3", 3)], value=1, style={"button_width": "auto"})
show = w.ToggleButtons(options=[("Built", "built"), ("Relaxed", "relaxed"), ("Side by side", "split")], value="built",
                       style={"button_width": "auto"})
show_cell = w.Checkbox(value=True, description="Cell", indent=False, layout=w.Layout(width="auto"))
viewer = w.Image(format="png", layout=w.Layout(width="100%", height="auto"))    # scales to the width available
tables, files = w.HTML(), w.HTML()
results = w.VBox([readout, w.HBox([view, repeat]), w.HBox([show, show_cell]), viewer, tables, files],
                 layout=w.Layout(flex="1 1 560px", min_width="0"))
state = dict(outputs=None, built=None, relaxed=None)


def collect(relax):
    inputs = getSimToolInputs(loc)
    for name, widget in f.items():
        getattr(inputs, name).value = widget.value
    inputs.run_relaxation.value = relax
    return inputs


def sub(formula):
    return "".join(f"<sub>{c}</sub>" if c.isdigit() else c for c in formula)


def dot(element):
    rgb = "".join(f"{int(255 * c):02x}" for c in jmol_colors[atomic_numbers[element]])
    return f"<span style='display:inline-block;width:0.8em;height:0.8em;border-radius:50%;background:#{rgb};" \
           f"border:1px solid #888;margin-right:0.3em'></span>"


def table(title, head, rows, hint=""):
    th = "".join(f"<th style='text-align:{'left' if i == 0 else 'right'}'>{h}</th>" for i, h in enumerate(head))
    body = "".join("<tr>" + "".join(f"<td style='text-align:{'left' if i == 0 else 'right'}'>{c}</td>"
                                    for i, c in enumerate(r)) + "</tr>" for r in rows)
    hint = f"<small style='color:#666'>{hint}</small>" if hint else ""
    return f"<h4 style='margin:1em 0 0.3em'>{title}</h4>{hint}<table><thead><tr>{th}</tr></thead><tbody>{body}</tbody></table>"


fmt = lambda v, d=3: "–" if v is None else f"{v:.{d}f}"
pct = lambda v: f"{100 * v:+.2f} %".replace("-", "−")
chip = lambda text: f"<span style='background:#e8eef5;border-radius:1em;padding:0.1em 0.7em;margin-right:0.3em'>{text}</span>"


def render():
    o = state["outputs"]
    mx, mt, it, m = o["mxene"], o["metal"], o["interface"], o["match"]
    label = lambda e, s: GREEK[s] if e in ("C", "N") else s
    seq = " ".join(f"{e}({label(e, s)})" for e, s in zip(mx["layers"], mx["sites"])) + " &nbsp;|&nbsp; " + \
          " ".join(f"{mt['symbol']}({s})" for s in mt["sites"])
    bx, by, bz = it["box"]
    readout.value = (f"<h2 style='margin:0'>{sub(mt['symbol'])}(111) / {sub(mx['formula'])}(0001)</h2>"
                     f"<p style='font-family:monospace;margin:0.3em 0'>{seq}</p><p>"
                     + "".join(chip(c) for c in [it["formula"], f"{it['n_atoms']} atoms", f"{bx:.2f} × {by:.2f} × {bz:.2f} Å",
                                                 f"εx {pct(m['x']['strain'])}", f"εy {pct(m['y']['strain'])}"])
                     + f"</p><small style='color:#666'>MXene a = {mx['a']} Å ({html.escape(mx['source'])}); "
                       f"{mt['name']} a = {mt['a']} Å ({html.escape(mt['source'])})</small>")

    parts = [table("Box", ["", "MXene replicas", "Metal replicas", "Box (Å)", "MXene strain"],
                   [[axis, f"{m[axis]['mxene']} × {mx['cell'][i]:.4f}", f"{m[axis]['metal']} × {mt['cell'][i]:.4f}",
                     f"{m[axis]['length']:.4f}", pct(m[axis]["strain"])] for i, axis in enumerate("xy")])]
    L = it["layers"]
    parts.append(table("Layers", ["#", "Layer", "Atoms", "z (Å)", "Spacing (Å)"], [
        [layer["index"], "".join(dot(e) for e in layer["species"]) + f"{list(layer['species'])[0]} "
         f"({label(list(layer['species'])[0], layer['site'])})",
         " ".join(f"{e}{n}" for e, n in layer["species"].items()), fmt(layer["z"]),
         fmt(layer["z"] - L[layer["index"] - 1]["z"]) if layer["index"] else "–"] for layer in reversed(L)]))
    r, c = o.get("relaxation"), o.get("comparison")
    if r:
        parts.append(table("Relaxation", ["", ""], [
            ["Converged", f"{r['converged']} after {r['n_steps']} FIRE steps, fmax {r['fmax_final']:.4f} eV/Å"],
            ["Energy (eV)", f"{r['energy_initial']:.4f} → {r['energy_final']:.4f}"],
            ["Cell a, b (Å), γ (°)", f"{r['cell']['a']}, {r['cell']['b']}, {r['cell']['gamma']}"],
            ["Gap (Å)", f"{c['gap']['initial']} → {c['gap']['relaxed']}"],
            ["MXene thickness (Å)", f"{c['mxene_thickness']['initial']} → {c['mxene_thickness']['relaxed']}"],
            ["Max |Δr| (Å)", f"{c['max_displacement']['norm']} ({c['max_displacement']['species']} "
                             f"#{c['max_displacement']['atom']}), RMS {c['rms_displacement']}"],
            ["Work of adhesion (J/m²)", fmt(r.get("work_of_adhesion"))],
        ]))
        parts.append(table("Built vs relaxed, layer by layer",
                           ["#", "Part", "Δz", "Spacing built", "Spacing relaxed", "Rumpling", "Slide", "⟨|Δr|⟩"],
                           [[x["index"], x["part"], fmt(x["dz"]), fmt(x["spacing_initial"]), fmt(x["spacing_relaxed"]),
                             fmt(x["rumpling"]), fmt(x["slide"]), fmt(x["mean_displacement"])] for x in reversed(c["layers"])],
                           "Lengths in Å. Rumpling is the z range of the relaxed layer, slide its mean in-plane displacement."))
    tables.value = "".join(parts)

    links = []
    for name, text in state["files"].items():
        data = base64.b64encode(text.encode()).decode()
        links.append(f"<a download='{name}' href='data:text/plain;base64,{data}' style='margin-right:1em'>{name}</a>")
    files.value = "<h4 style='margin:1em 0 0.3em'>Files</h4>" + "".join(links) + \
                  "<br><small style='color:#666'>Also written next to this notebook.</small>"


def draw(*_):
    if state["built"] is None:
        return
    both = show.value == "split" and state["relaxed"] is not None
    shown = [("built", state["built"]), ("relaxed", state["relaxed"])] if both else \
            [(show.value, state[show.value] if state[show.value] is not None else state["built"])]
    views = [("front (x–z)", "-90x"), ("side (y–z)", "-90y,-90z"), ("top (x–y)", "")] if view.value == "all" else \
            [(dict((v, k) for k, v in view.options)[view.value].lower(), view.value)]
    fig, axes = plt.subplots(len(shown), len(views), figsize=(3.4 * len(views) + 1, 4.5 * len(shown)), squeeze=False)
    for i, (which, atoms) in enumerate(shown):
        for j, (title, rotation) in enumerate(views):
            plot_atoms(atoms.repeat((repeat.value, repeat.value, 1)), axes[i, j], rotation=rotation, radii=0.5,
                       show_unit_cell=2 if show_cell.value else 0)
            axes[i, j].set_title(f"{which}, {title}" if len(shown) > 1 else title)
            axes[i, j].set_axis_off()
    fig.tight_layout()
    png = BytesIO()
    fig.savefig(png, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    viewer.value = png.getvalue()


def watch(path, done, start):
    """Shows the last FIRE step the Sim2L wrote to progress.jsonl, like the web app, until done is set."""
    while not done.wait(0.5):
        lines = path.read_text().splitlines() if path.exists() else []
        try:
            s = json.loads(lines[-1])
        except (IndexError, ValueError):                       # nothing yet, or a line half written
            continue
        progress.value = min(s["step"], progress.max)
        relax_status.value = (f"Step {s['step']} · E = {s['energy']:.4f} eV · fmax = {s['fmax']:.4f} eV/Å"
                              f" · {time.time() - start:.0f} s")


def run(relax):
    build_button.disabled = relax_button.disabled = True
    status.value = "" if relax else f"{SPIN}<i>Building…</i>"
    name, done = uuid.uuid4().hex, threading.Event()       # the run folder is RUNS/<name>, so it can be watched;
                                                           # absolute, as papermill changes directory into it
    watcher = threading.Thread(target=watch, args=(Path(get_experiment(), name, "progress.jsonl").resolve(), done, time.time()),
                               daemon=True)
    if relax:
        relax_status.value = f"{SPIN}<i>Loading MACE {html.escape(f['model'].value)}…</i>"
        progress.max, progress.value, progress.layout.visibility = f["max_steps"].value, 0, "visible"
        watcher.start()
    try:
        with contextlib.redirect_stdout(StringIO()), contextlib.redirect_stderr(StringIO()), warnings.catch_warnings(), \
                w.Output():                                    # never shown: papermill's progress bar
            warnings.simplefilter("ignore")                    # simtool's and papermill's run log
            r = Run(loc, collect(relax), runName=name)
        done.set()
        if watcher.is_alive():
            watcher.join()                                     # its last step line must not overwrite the summary
        poscar_relaxed = r.read("POSCAR_relaxed")
        state["outputs"] = r.read("summary")
        mx, mt = state["outputs"]["mxene"]["formula"], state["outputs"]["metal"]["symbol"]
        state["files"] = {name: text for name, text in [("POSCAR", r.read("POSCAR")), (f"POSCAR_{mx}", r.read("POSCAR_mxene")),
                          (f"POSCAR_{mt}", r.read("POSCAR_metal")), ("POSCAR_relaxed", poscar_relaxed),
                          ("outputs.json", json.dumps(state["outputs"], indent=2) + "\n")] if text}
        for name, text in state["files"].items():
            Path(name).write_text(text)
        state["built"] = read(StringIO(r.read("POSCAR")), format="vasp")
        state["relaxed"] = read(StringIO(poscar_relaxed), format="vasp") if poscar_relaxed else None
        if relax and state["relaxed"] is not None:
            show.value = "relaxed"
        elif state["relaxed"] is None:
            show.value = "built"
        status.value = ""
        if relax:
            rx, i = state["outputs"]["relaxation"], state["outputs"]["inputs"]
            verdict = "stopped" if rx["stopped"] else "converged" if rx["converged"] else "not converged"
            relax_status.value = (f"{i['model']} ({i['head']}){' + D3(BJ)' if i['dispersion'] else ''}: {verdict} after "
                                  f"{rx['n_steps']} steps, E = {rx['energy_final']:.4f} eV, fmax = {rx['fmax_final']:.4f} eV/Å"
                                  + (" (from the cache)" if r.cached else ""))
        render()
        draw()
    except Exception as error:
        lines = [line for line in str(error).splitlines() if line.strip()] or [type(error).__name__]
        status.value = f"<span style='color:#b00020'>{html.escape(lines[-1])}</span>"   # last line of the Sim2L traceback
        if relax:
            relax_status.value = ""
    finally:
        done.set()
        progress.layout.visibility = "hidden"
        build_button.disabled = relax_button.disabled = False


build_button.on_click(lambda _: run(False))
relax_button.on_click(lambda _: run(True))
for name in ("mxene", "stacking", "mxene_rotation", "metal", "metal_layers", "metal_site", "metal_stacking",
             "gap", "shift_x", "shift_y", "tol", "max_replicas", "vacuum"):
    f[name].observe(lambda _: auto.value and run(False), names="value")
for widget in (view, repeat, show, show_cell):
    widget.observe(draw, names="value")

def app():
    """Builds the interface with the current inputs, then shows the form in place of the loading message."""
    global loading
    display(HTML("<style>.container { width: 98% !important; }"                  # classic Notebook
                 " :root, body { --jp-notebook-max-width: 100% !important; }"     # Notebook 7
                 f" {SPIN_CSS}</style>"))
    message = HTML(f"<p>{SPIN}Building the first interface…</p>")
    if loading:
        loading.update(message)
    else:
        loading = display(message, display_id=True)
    run(False)
    loading.update(HTML(""))
    loading = None                                         # a second app() shows its own message
    display(w.VBox([status, w.HBox([form, results], layout=w.Layout(flex_flow="row wrap"))]))
