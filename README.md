# Coherent MXene/metal interface builder

Builds an MXene/metal interface in one orthogonal box, relaxes it with a MACE machine-learned potential and compares the built and relaxed structures. It runs as a web app with a 3D viewer, or from the command line.

## What it does

1. **MXene.** Ordered MXene (e.g. Ti₃C₂O₂, Mo₂TiC₂O₂) in its conventional a × √3a × c cell, with AβC (octahedral) or AβA (trigonal prismatic) stacking, optionally turned 180° about z. Vacuum is only below the MXene, so its top layer is at z = c.
2. **Metal.** Conventional cubic cell built from its space group: FCC `Fm-3m` (4a), or L1₂ `Pm-3m` (B at 1a, A at 3c). It is then turned so that [111] is along z: a stack of close-packed (111) layers, d₁₁₁ = a_cubic/√3 apart, in the same a × √3a form as the MXene with a → a₁₁₁ (a₁₁₁ × √3a₁₁₁ for FCC, 2a₁₁₁ × 2√3a₁₁₁ for L1₂, whose B atoms form a 2 × 2 superlattice in each layer).
3. **Supercell.** The metal is the substrate and the MXene is added on it. Along x and along y, the tool finds the smallest whole number of metal replicas that matches a whole number of MXene replicas within the strain tolerance. The metal keeps its lattice constant and the MXene takes the strain (see [Supercell matching](#supercell-matching)).
4. **Placement.** Metal (111) on MXene (0001), metal [1-10] along x, layer by layer: the first metal (111) plane sits at a set gap above the top MXene layer, on site A, B or C of the MXene grid, and the metal layers then follow A → B → C (or its twin). An optional in-plane shift moves the metal off the site.
5. **Relaxation.** MACE (default MACE-MH-1, `omat_pbe` head) in float64, optional D3(BJ) dispersion and FIRE. The outermost metal layers are held fixed, and the in-plane cell can relax too. Optionally reports the work of adhesion.
6. **Built vs relaxed.** Per-atom displacements, layer heights and spacings, the gap, the MXene thickness and the cell, before and after relaxation (see below).
7. **Atomic planes.** Shows the (hkl) planes of the strained MXene cell or of the metal cell in the viewer.

Lattice constants come from `data/MXene_FCC_L12_strain_matrix.xlsx` (41 ordered MXenes, 23 metals). Solid solutions, Ti₃C₂Tₓ and Ni₃(Al,Ti) are left out.

Fixed geometry rules for the starting structure:
- M–X layer spacing is 0.385 a.
- The O height comes from a 2.0 Å M–O bond.
- Sites are counted from the bottom MXene layer, which is on A. AβC cycles A → B → C upward; AβA puts the metals on A and C/N on B.
- O sits on the site that neither of its two nearest layers uses.
- `mxene_rotation` = 180 turns the MXene 180° about z. A stays and B ↔ C swap, so AβC becomes AγB (A C B from the bottom) and the C layer moves from the right of A to its left.
- With `metal_site` = `auto`, the first metal layer takes the same rule: the site free of the top two MXene layers, so the close-packed stacking carries on into the metal.

## Supercell matching

The MXene basal plane and the metal (111) plane are both hexagonal nets, written as rectangular cells: a × √3a for the MXene and a₁₁₁ × √3a₁₁₁ for the metal (`a_111 (A)` in the spreadsheet; twice that for L1₂). The two in-plane directions are matched separately. For each direction, with cell lengths L_MX (a along x, √3a along y) and L_M (a₁₁₁ or √3a₁₁₁):

1. Try n = 1, 2, … up to `max_replicas` metal replicas.
2. For each n, take the whole number of MXene replicas closest to the same length: m = round(n · L_M / L_MX), at least 1.
3. The MXene strain is ε = n · L_M / (m · L_MX) − 1.
4. Keep the first n with |ε| ≤ `tol`. This gives the smallest supercell that fits the tolerance.

The box is then n_x · a₁₁₁ by n_y · √3a₁₁₁ (twice that for L1₂). The metal keeps its lattice constant, and the MXene is strained in-plane by ε_x and ε_y to fill the box exactly. Its layer heights are those of the unstrained MXene. Both cells have the same √3 aspect ratio, so x and y give the same m, n and ε: the MXene takes an equi-biaxial strain. At m = n = 1 it is ε = −δ/(1 + δ), where δ = (a_MXene − a₁₁₁)/a₁₁₁ is the value in the spreadsheet's `Strain` sheet (the metal strain if the metal were strained instead). If no n up to `max_replicas` works, the build stops with an error; raise `tol` or `max_replicas`.

For example:
- V₂C (a = 2.881 Å) on Au (a₁₁₁ = 2.884 Å): 1 × 1, ε = +0.09 %.
- Ti₃C₂ (a = 3.104 Å) on Au: 1 × 1 would be −7.1 %. The first coincidence within 2 % is 11 MXene cells on 12 Au cells (ε = +1.35 %), so it needs `max_replicas` ≥ 12; 13 : 14 gives +0.05 %.

The *Box* panel in the web app and `outputs.json → match` report m, n, the box length and ε for each direction.

## Run

With Podman (the image includes CPU PyTorch and the MACE-MH-1 weights):

```bash
podman build -t coherent-interface-builder:2.1.0 .
podman run --rm -p 8000:8000 coherent-interface-builder:2.1.0          # open http://localhost:8000
podman run --rm -p 8000:8000 -v "$PWD/data:/app/data:Z" coherent-interface-builder:2.1.0   # use your own copy of the spreadsheet
```

Without the container (Python 3.11 or newer):

```bash
pip install torch --index-url https://download.pytorch.org/whl/cpu   # or the CUDA build
pip install -r requirements.txt
flask --app app run --port 8000
```

From the command line, for example on a cluster:

```bash
python builder.py inputs.json     # build only: POSCAR, POSCAR_<mxene>, POSCAR_<metal>, outputs.json
python relax.py inputs.json       # build and relax: adds POSCAR_relaxed and the comparison to outputs.json
```

`inputs.json` holds only the inputs you want to change; everything else takes its default. Without a file, both scripts use the defaults. The web app can save `inputs.json`. It can also load `inputs.json`, or an `outputs.json` to restore that run's inputs.

## nanoHUB SimTool

`simtool/coherentmxmetal.ipynb` is the same workflow as a nanoHUB [SimTool](https://simtool.readthedocs.io/) (Sim2L): a notebook that declares its inputs and outputs in YAML, checks the inputs, and saves the outputs to nanoHUB's results database. Runs with the same inputs come from the cache.

- **Inputs:** the inputs above, with `shift` split into `shift_x` and `shift_y`, plus `run_relaxation` (off by default: build only). Lengths take units, e.g. `gap = "0.25 nm"`.
- **Outputs:** `formula`, `n_atoms`, `box`, `strain_x`, `strain_y`, `match`, the POSCARs as text, and after a relaxation `energy`, `gap_relaxed` and `work_of_adhesion` (empty otherwise). `summary` holds everything in `outputs.json`.
- `simtool/builder.py` and `simtool/relax.py` are links to the top-level files. The notebook lists them and `sim2l.py` as `EXTRA_FILES`, so each run gets them, and `builder.py` finds the spreadsheet through the link.
- `landing_page.ipynb` is the page nanoHUB opens: abstract, links and references.
- `bin/interactive.ipynb` is a form like the web app (ipywidgets): menus for the MXene, metal and placement, a rebuild after each change, a Relax button with live progress (step, energy and fmax, read from the `progress.jsonl` the Sim2L writes in its run folder at every FIRE step), the structure views, tables and POSCAR downloads. Every build runs the Sim2L.
- Both notebooks call `simtool.Run(..., venue="local")` on nanoHUB, which runs the Sim2L with `submit --local` in the user's run folder. For a published Sim2L simtool's default goes through nanoHUB's shared-cache helpers (`ionhelper`), which failed in this tool (`No module named 'simtool'` / `'hubzero'`). Local runs use the user's own cache, and let the interactive builder show relaxation progress. Without `submit` (off nanoHUB), simtool's default is used.
- `bin/workflow.ipynb` sets the inputs, runs the Sim2L with `simtool.Run`, and shows the match, a side and top view and the POSCARs. It runs in `anaconda-7`'s `python3` kernel, which has `simtool` and `papermill`.
- `middleware/invoke` (nanoHUB's Sim2L template) opens the landing page in a nanoHUB Jupyter session with `anaconda-7`. The Sim2L notebook itself runs in that module's `torchsim` kernel (Python 3.12, torch 2.8, mace-torch 0.3.14, pymatgen 2025.10, ase 3.27). That kernel has neither `simtool` nor `yamlmagic`, so `simtool/sim2l.py` stands in for them there (it uses them when they are installed); the caller's `simtool` reads the outputs as usual. mace-torch 0.3.14 does not know the name `mh-1`, so `relax.py` then loads MACE-MH-1 from its download URL.
- The MXene and metal option lists in the notebook are the spreadsheet rows. Update them if you add rows.

Besides `requirements.txt` (without Flask and gunicorn), the SimTool needs `simtool`, `papermill` and `yamlmagic`. To test it locally, from the repository root:

```python
from simtool import searchForSimTool, getSimToolInputs, Run
loc = searchForSimTool("coherentmxmetal")
inputs = getSimToolInputs(loc)
inputs.gap.value = "0.25 nm"
r = Run(loc, inputs)            # writes RUNS/<id>/
r.read("POSCAR"), r.read("strain_x")
```

## Inputs (`inputs.json`)

| Name | Default | Meaning |
|---|---|---|
| `mxene` | `V2C` | Formula, a row of the MXenes sheet |
| `stacking` | `ABC` | `ABC` = AβC (octahedral), `ABA` = AβA (trigonal prismatic) |
| `mxene_rotation` | 0 | degrees about z, `0` or `180` (B ↔ C) |
| `vacuum` | 20.0 | Å, below the MXene only |
| `metal` | `Au` | Symbol, a row of the Metals sheet |
| `metal_layers` | 4 | metal (111) atomic planes |
| `metal_site` | `auto` | site of the first metal layer on the MXene grid: `A` (in-plane position of the bottom MXene layer), `B`, `C`, or `auto` (free of the top two MXene layers) |
| `metal_stacking` | `ABC` | `ABC` climbs A → B → C (continues an MXene at 0°), `ACB` climbs A → C → B (continues one at 180°) |
| `gap` | 2.3 | Å, top MXene layer to first metal (111) layer |
| `shift` | [0, 0] | Å, in-plane shift of the metal along x, y, from its site |
| `tol` | 0.02 | largest \|strain\| on the MXene along x and y |
| `max_replicas` | 10 | largest number of metal replicas tried per direction |
| `model` | `mh-1` | MACE foundation model |
| `head` | `omat_pbe` | output head of the MH models: `omat_pbe` (PBE(+U)), `oc20_usemppbe` (PBE slabs), `matpes_r2scan` (r²SCAN) |
| `dispersion` | false | add D3(BJ) |
| `fmax` | 0.01 | eV/Å, force convergence |
| `max_steps` | 500 | FIRE steps |
| `fixed_layers` | 2 | outermost metal layers held fixed |
| `relax_cell` | true | relax a, b and γ; c stays fixed |
| `work_of_adhesion` | false | from the rigidly separated MXene and metal |

## Outputs (`outputs.json`)

- `inputs`: every input used, defaults filled in.
- `mxene`: spreadsheet row, layers, sites, the cell before and after strain, the strained cell vectors and the thickness.
- `metal`: spreadsheet row (with a₁₁₁), space group, d₁₁₁, orientation, layer sites, the (111) cell and the cubic cell vectors.
- `match`: replicas, box length and MXene strain along x and y.
- `interface`: formula, atom count, box and per-layer table.
- `spreadsheet`: file name and SHA-256.
- `files`: `POSCAR` (interface), `POSCAR_<mxene>`, `POSCAR_<metal>`.

After a relaxation, outputs.json also has:
- `relaxation`: optimizer, fixed atoms, converged, steps, initial and final energy, final fmax, cell a, b and γ, MXene–metal gap, package versions and the energy/fmax per step. Also the work of adhesion in J/m², when asked.
- `comparison`: built vs relaxed, below.
- `files.POSCAR_relaxed`.

## Built vs relaxed

After a relaxation, the viewer has four modes:

- **Built**: the starting structure.
- **Relaxed**: the relaxed structure.
- **Overlay**: the relaxed atoms over the faded built ones, with displacement arrows (scale ×1 to ×20, or off).
- **Side by side**: built on the left, relaxed on the right, with one shared camera.

Hovering over an atom in the overlay shows its displacement.

The *Built vs relaxed* panel and `outputs.json → comparison` give:

| Key | Content |
|---|---|
| `max_displacement` | largest \|Δr\|, with the atom number in `POSCAR_relaxed` (from 1), element and layer |
| `rms_displacement` | RMS \|Δr\| over all atoms |
| `max_displacement_part` | largest \|Δr\| in the MXene and in the metal |
| `gap` | MXene–metal gap (mean z of the first metal layer minus mean z of the top MXene layer), built and relaxed |
| `mxene_thickness` | bottom to top MXene layer, built and relaxed |
| `cell` | a, b and γ, built and relaxed |
| `layers` | per layer: mean z, Δz, spacing to the layer below, rumpling (z range in the relaxed layer), slide (mean in-plane displacement) and mean \|Δr\| |
| `largest` | the 15 atoms that moved most, with Δx, Δy, Δz and \|Δr\| |

Displacements are Cartesian, using the minimum image in-plane. When the in-plane cell relaxes, the change in cell size adds to them. All lengths are in Å.

## Atomic planes

(hkl) refer to the MXene's strained a × √3a × c cell (turned with the MXene, so at 180° its a and b point along −x and −y) or to the metal's cubic cell, whose [111] is along z, so metal (111) planes are the horizontal layers. The MXene c includes the vacuum, so MXene planes with l ≠ 0 tilt with it.

## Files

| File | Content |
|---|---|
| `builder.py` | builds the interface from the inputs |
| `relax.py` | MACE relaxation, work of adhesion and the built-vs-relaxed comparison |
| `planes.py` | atomic planes for the viewer |
| `app.py` | Flask app |
| `templates/`, `static/` | page, styles, script and bundled 3Dmol.js |
| `data/` | lattice-constant spreadsheet |
| `simtool/` | nanoHUB SimTool notebook (`coherentmxmetal.ipynb`) |
| `landing_page.ipynb` | nanoHUB landing page |
| `bin/interactive.ipynb` | form on nanoHUB, like the web app, that runs the Sim2L |
| `bin/interactive_ui.py` | the form's code (ipywidgets); the notebook cell only calls `app()` |
| `bin/workflow.ipynb` | runs the Sim2L on nanoHUB |
| `middleware/invoke` | nanoHUB launcher |
| `src/Makefile`, `doc/`, `examples/` | nanoHUB tool layout (nothing to compile) |
