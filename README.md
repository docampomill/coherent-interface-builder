# MXene/metal interface builder

Builds MXene/metal interfaces layer by layer from the bottom and relaxes them with MACE-MH-1.
The science is in the Python package `mxbuilder`. A small Flask app, a command line and a
Podman container sit on top of it. The inputs and outputs are listed below so the same
workflow can become a papermill notebook and a nanoHUB Sim2L.

## Run it

```bash
podman build -t mxbuilder:0.1.0 .
mkdir -p runs
podman run --rm -p 8000:8000 -v "$PWD/runs:/app/runs:Z" mxbuilder:0.1.0
# open http://localhost:8000
```

Command line and tests in the same image:

```bash
podman run --rm -v "$PWD/runs:/app/runs:Z" mxbuilder:0.1.0 python -m mxbuilder run examples/au5_ti3c2o2.json
podman run --rm mxbuilder:0.1.0 pytest -q tests
```

Without the container (Python 3.11 or newer):

```bash
pip install torch --index-url https://download.pytorch.org/whl/cpu   # or the CUDA build
pip install -r requirements.txt
flask --app webapp.app run --port 8000
```

On a Mac, Podman runs a Linux VM, so relaxations run on the CPU. MACE cannot use Apple MPS in float64 anyway.

## Layout

```
mxbuilder/        the science; no Flask code
  schema.py       inputs and outputs: names, defaults, units, limits
  elements.py     reference data and the default geometry rules
  stack.py        A/B/C sites, and turning "auto" inputs into numbers
  lattice.py      coherent and coincidence supercells, filling a layer with atoms
  build.py        inputs → ase.Atoms, layer by layer from the bottom
  analysis.py     layer heights, spacings, gaps, displacements
  export.py       POSCAR, extxyz, CIF, LAMMPS data
  relax.py        MACE-MH-1, relaxation, work of adhesion, MLIP lattice constants
  workflow.py     run() and scan(): inputs → run folder and outputs.json
  cli.py          python -m mxbuilder …
webapp/           Flask routes, background jobs, page (HTML, CSS, JS, bundled 3Dmol.js)
tests/            pytest; no model needed (the MLIP is replaced by a stand-in)
examples/         input files
```

## Conventions

- **Cell.** Hexagonal, a₁ = a(1, 0, 0), a₂ = a(−½, √3/2, 0), γ = 120°, c normal to the interface.
- **Sites.** A = (0, 0), B = (⅓, ⅔), C = (⅔, ⅓) in fractions of a₁, a₂. Latin letters mark metal and O layers, Greek letters (α, β, γ) mark C/N layers.
- **Layers** are numbered from 0 at the bottom. Each atom's ASE tag is its layer index. Fixed atoms, the layer table and all analysis refer to layer indices, never atom indices.
- **Stacking.** AβC cycles A → B → C through every layer: O(C) Ti(A) C(β) Ti(C) C(α) Ti(B) O(C). AβA keeps the metal layers on A and C/N on B. An outer O sits on the site that neither of its two neighbouring layers uses. For AβC that is the fcc hollow over the inner metal.
- **Metal.** fcc(111). By default its first layer continues the stacking. `metal_site` overrides that, and `rbt` shifts the whole slab rigidly in fractions of the MXene cell. In a sandwich the bottom slab mirrors the top one, so a symmetric MXene gives an inversion-symmetric cell.
- **Matching.** The MXene keeps its a and the metal takes all the strain. Supercell `[i, j]` means i·a₁ + j·a₂ and is matched to `[k, l]` of the metal. Twist is reported in (−30°, 30°].
- **Gap.** The difference of layer-mean z across an MXene/metal boundary.
- **Atom order in files.** Species in order of first appearance from the bottom, then layer, then y, x. The same inputs give byte-identical files.
- **Units.** Å, eV, eV/Å, J/m², and meV/Å² for the shift scan.
- **Starting-geometry rules.** These are the same as Yi's builder and only set the starting point:
  - d(M–X) = 0.385·a
  - O height from a 2.00 Å M–O bond to the hollow site
  - gap 2.30 Å, plus 0.30 Å when on top
  - metal spacing d₁₁₁ with a Poisson response to the strain
  - outer two metal layers fixed

## Reproducibility

Each run lives in `runs/<run id>/`:

| File | Content |
|---|---|
| `inputs.json` | Inputs as given; `null` means auto |
| `resolved.json` | Every value actually used. Building from it gives the same structure even if a default changes later |
| `initial.{vasp,extxyz,cif,lmp}` | Structure before relaxation |
| `relaxed.{vasp,extxyz,cif,lmp}` | Structure after relaxation |
| `relax.log`, `relax.traj` | ASE optimizer log and trajectory |
| `outputs.json` | The outputs below |
| `environment.json` | Python and package versions, MACE model file and its SHA-256 |

- The run id is the SHA-256 of the canonical `inputs.json` plus the code version. A repeated request returns the stored outputs instead of rerunning, the same way a Sim2L caches.
- The image pins its packages, installs CPU PyTorch, downloads the MACE-MH-1 weights at build time and saves `pip freeze` as `/app/pip-freeze.txt`. Tag the image and keep it (`podman save`) with the data it produced.
- MACE runs in float64.

## Inputs

`python -m mxbuilder inputs` prints this table from `mxbuilder/schema.py`.

| name | group | type | default | units | allowed | description |
|---|---|---|---|---|---|---|
| mxene_layers | MXene | list[str] | ['O', 'Ti', 'C', 'Ti', 'C', 'Ti', 'O'] |  |  | Element of each layer, bottom → top (TM or O; C or N) |
| stacking | MXene | str | ABC |  | ABC, ABA | AβC (octahedral) or AβA (trigonal prismatic) |
| mxene_sites | MXene | list[str \| None] \| None |  |  | A, B, C | Site of each layer; null (or a null entry) = from stacking |
| mxene_dz | MXene | list[float \| None] \| None |  | Å |  | Spacing of each layer to the one below (first entry 0); null (or a null entry) = rules |
| mxene_a | MXene | float \| None |  | Å | 2.0 – 4.5 | In-plane lattice constant; null = table |
| metal | Metal | str | Au |  | Au, Ag, Cu, Pt, Pd, Ni, Al | fcc metal, (111) planes parallel to the MXene |
| metal_layers | Metal | int | 5 |  | 0 – 20 | Number of metal (111) layers; 0 = MXene only |
| metal_a | Metal | float \| None |  | Å | 3.0 – 5.0 | Bulk fcc lattice constant; null = table |
| metal_dz | Metal | float \| None |  | Å | 1.0 – 4.0 | Metal layer spacing; null = d111 with Poisson response |
| gap | Metal | float \| None |  | Å | 0.8 – 8.0 | Vertical gap between MXene and metal; null = rule |
| metal_site | Metal | str \| None |  |  | A, B, C | Site of the first metal layer; null = continue the stacking |
| rbt | Metal | list[float] | [0.0, 0.0] |  |  | Rigid-body translation of the metal, fractions of the MXene a1, a2 |
| geometry | Cell | str | slab |  | slab, sandwich, periodic | Metal on top; on both faces; or no vacuum (two interfaces through PBC) |
| mxene_supercell | Cell | list[int] | [1, 0] |  |  | MXene supercell vector [i, j] in the 120° hexagonal basis |
| metal_supercell | Cell | list[int] | [1, 0] |  |  | Metal supercell vector [k, l] matched onto it |
| repeat | Cell | list[int] | [1, 1] |  | 1 – 10 | Extra repeats along a1, a2 |
| vacuum | Cell | float | 15.0 | Å | 0.0 – 60.0 | Vacuum along c (slab, sandwich) |
| relax | Relax | bool | True |  |  | Relax with the MLIP; false = build only |
| model | Relax | str | mh-1 |  | mh-1 | MACE foundation model |
| head | Relax | str | omat_pbe |  | omat_pbe, oc20_usemppbe, matpes_r2scan | MACE-MH-1 output head |
| dispersion | Relax | bool | False |  |  | Add D3(BJ) dispersion |
| fmax | Relax | float | 0.05 | eV/Å | 0.001 – 1.0 | Force convergence |
| max_steps | Relax | int | 500 |  | 1 – 10000 | Maximum optimizer steps |
| optimizer | Relax | str | FIRE |  | FIRE, BFGS | ASE optimizer |
| fixed_layers | Relax | list[int] \| None |  |  |  | Layer indices (bottom = 0) held fixed; null = outer 2 metal layers |
| relax_cell | Relax | bool | False |  |  | Also relax a, b and γ (c stays fixed) |
| mlip_lattice | Relax | bool | False |  |  | Take auto lattice constants from the MLIP instead of the table |
| work_of_adhesion | Relax | bool | True |  |  | Rigid-separation work of adhesion |

The shift scan (`workflow.scan`) takes the same inputs plus `grid` (shifts per cell edge) and `relax_z` (relax along z only). It writes `runs/scan-<id>/scan.json` with ΔE per interface area for each shift.

## Outputs

`python -m mxbuilder outputs` prints this table.

| name | type | units | description |
|---|---|---|---|
| run_id | str |  | Hash of inputs.json and the code version; also the run folder name |
| formula | str |  | Interface label, e.g. Au5/Ti3C2O2 |
| n_atoms | int |  | Atoms in the cell |
| cell | dict | Å, ° | a, b, c, gamma |
| area | float | Å² | Interface area of the cell |
| strain | float | % | In-plane strain of the metal (the MXene keeps its a) |
| twist | float | ° | Rotation of the metal lattice relative to the MXene |
| n_interfaces | int |  | 1 for slab, 2 for sandwich and periodic |
| mlip_lattice | dict | Å | MLIP-derived mxene_a and metal_a (when mlip_lattice is on) |
| layers | list | Å | Per layer: element, site, part, atoms, z and spacing before (and after) relaxation |
| gaps_initial | dict | Å | MXene–metal gaps before relaxation |
| gaps_relaxed | dict | Å | MXene–metal gaps after relaxation |
| energy_initial | float | eV | MLIP energy of the starting structure |
| energy_final | float | eV | MLIP energy after relaxation |
| converged | bool |  | fmax reached within max_steps |
| n_steps | int |  | Optimizer steps taken |
| fmax_final | float | eV/Å | Largest residual force |
| work_of_adhesion | float | J/m² | (E_metal + E_MXene − E_interface) / (n_interfaces · area) |
| energy_metal | float | eV | Metal slab alone, relaxed-interface geometry |
| energy_mxene | float | eV | MXene alone, relaxed-interface geometry |
| trace | list | eV, eV/Å | Energy and fmax per step |
| files | dict |  | Structure files: initial.* and relaxed.* (POSCAR, extxyz, cif, lammps-data) |
| warnings | list |  | Notes about the starting geometry |

## Toward papermill and the Sim2L

Every input is a scalar or a short list, so a notebook's parameters cell is just the inputs table:

```python
# parameters
mxene_layers = ["O", "Ti", "C", "Ti", "C", "Ti", "O"]
metal = "Au"
metal_layers = 5
dispersion = True
```

```python
from mxbuilder import workflow
from mxbuilder.schema import Inputs

out = workflow.run(Inputs(mxene_layers=mxene_layers, metal=metal, metal_layers=metal_layers,
                          dispersion=dispersion), "runs")
```

Each row of the inputs table maps to one Sim2L input (number, integer, choice, boolean or list) with the limits shown. Each row of the outputs table maps to one Sim2L output. One thing to settle for the Sim2L is how to expose "auto" (`null`): either as an empty value or as an explicit "use default" switch.

## Compared with Yi's builder

- **Same starting geometry.** For Au₅/Ti₃C₂O₂ at a = 3.07 Å, MACE-MH-1 (omat_pbe) gives −77.857175 eV here and −77.857175 eV in Yi's run. With D3 the energies are −82.995888 and −82.995888 eV. The site letters are mirrored: A↔C.
- **Inputs per layer.** The input is the list of layers (element, site and spacing per layer), not only layer counts.
- **Runs keyed by inputs.** Runs are stored by input hash with `resolved.json` and `environment.json`.
- **New geometries.** Sandwich and periodic cells, a rigid-shift scan, and coincidence cells with an optional twist target.
- **New output.** LAMMPS data.
- **Not yet.** Metal planes other than (111), hcp/bcc metals, mixed-occupancy layers, and OH/F/Cl terminations. Element lists and rules live in `elements.py` and `stack.py`.
