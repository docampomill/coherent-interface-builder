# MXene/metal interface builder

Builds an MXene/metal interface in one orthogonal box, relaxes it with a MACE machine-learned potential and compares the built and relaxed structures. It runs as a web app with a 3D viewer, or from the command line.

## What it does

1. **MXene.** Ordered MXene (e.g. Ti₃C₂O₂, Mo₂TiC₂O₂) in its conventional a × √3a × c cell, with AβC (octahedral) or AβA (trigonal prismatic) stacking. Vacuum is only below the MXene, so its top layer is at z = c.
2. **Metal.** Conventional cubic cell built from its space group: FCC `Fm-3m` (4a), or L1₂ `Pm-3m` (B at 1a, A at 3c).
3. **Box.** Along x and along y, the tool finds the smallest number of MXene replicas that matches a whole number of metal replicas within the strain tolerance. The metal takes the strain.
4. **Placement.** Metal (001) on MXene (0001), metal [100] along x, at a set gap above the top MXene layer, with an optional in-plane shift.
5. **Relaxation.** MACE (default MACE-MH-1, `omat_pbe` head) in float64, optional D3(BJ) dispersion and FIRE. The outermost metal layers are held fixed, and the in-plane cell can relax too. Optionally reports the work of adhesion.
6. **Built vs relaxed.** Per-atom displacements, layer heights and spacings, the gap, the MXene thickness and the cell, before and after relaxation (see below).
7. **Atomic planes.** Shows the (hkl) planes of the MXene cell or of the strained metal cell in the viewer.

Lattice constants come from `data/MXene_FCC_L12_strain_matrix.xlsx` (41 ordered MXenes, 23 metals). Solid solutions, Ti₃C₂Tₓ and Ni₃(Al,Ti) are left out.

Fixed geometry rules for the starting structure:
- M–X layer spacing is 0.385 a.
- The O height comes from a 2.0 Å M–O bond.
- AβC puts the central layer on A. AβA puts the metals on A and C/N on B.
- O sits on the site that neither of its two nearest layers uses.

## Run

With Podman (the image includes CPU PyTorch and the MACE-MH-1 weights):

```bash
podman build -t mxbuilder:0.8.0 .
podman run --rm -p 8000:8000 mxbuilder:0.8.0          # open http://localhost:8000
podman run --rm -p 8000:8000 -v "$PWD/data:/app/data:Z" mxbuilder:0.8.0   # use your own copy of the spreadsheet
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

## Inputs (`inputs.json`)

| Name | Default | Meaning |
|---|---|---|
| `mxene` | `Ti3C2` | Formula, a row of the MXenes sheet |
| `stacking` | `ABC` | `ABC` = AβC (octahedral), `ABA` = AβA (trigonal prismatic) |
| `vacuum` | 20.0 | Å, below the MXene only |
| `metal` | `Au` | Symbol, a row of the Metals sheet |
| `metal_layers` | 4 | metal (001) atomic planes |
| `gap` | 2.3 | Å, top MXene layer to first metal layer |
| `shift` | [0, 0] | Å, in-plane shift of the metal along x, y |
| `tol` | 0.02 | largest \|strain\| on the metal along x and y |
| `max_replicas` | 10 | largest number of MXene replicas tried per direction |
| `model` | `mh-1` | MACE foundation model |
| `head` | `omat_pbe` | output head of the MH models: `omat_pbe` (PBE(+U)), `oc20_usemppbe` (PBE slabs), `matpes_r2scan` (r²SCAN) |
| `dispersion` | true | add D3(BJ) |
| `fmax` | 0.01 | eV/Å, force convergence |
| `max_steps` | 500 | FIRE steps |
| `fixed_layers` | 2 | outermost metal layers held fixed |
| `relax_cell` | true | relax a, b and γ; c stays fixed |
| `work_of_adhesion` | false | from the rigidly separated MXene and metal |

## Outputs (`outputs.json`)

- `inputs`: every input used, defaults filled in.
- `mxene`: spreadsheet row, layers, sites, cell and thickness.
- `metal`: spreadsheet row, space group and strained cell.
- `match`: replicas, box length and metal strain along x and y.
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

(hkl) refer to the MXene a × √3a × c cell or to the metal's strained cubic cell. The MXene c includes the vacuum, so MXene planes with l ≠ 0 tilt with it.

## Files

| File | Content |
|---|---|
| `builder.py` | builds the interface from the inputs |
| `relax.py` | MACE relaxation, work of adhesion and the built-vs-relaxed comparison |
| `planes.py` | atomic planes for the viewer |
| `app.py` | Flask app |
| `templates/`, `static/` | page, styles, script and bundled 3Dmol.js |
| `data/` | lattice-constant spreadsheet |
