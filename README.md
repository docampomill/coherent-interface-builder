# MXene/metal interface builder

Builds an MXene/metal interface in one orthogonal box, relaxes it with a MACE machine-learned potential and compares the built and relaxed structures. It runs as a web app with a 3D viewer, or from the command line.

## What it does

1. **MXene.** Ordered MXene (e.g. Ti₃C₂O₂, Mo₂TiC₂O₂) in its conventional a × √3a × c cell, with AβC (octahedral) or AβA (trigonal prismatic) stacking. Vacuum is only below the MXene, so its top layer is at z = c.
2. **Metal.** Conventional cubic cell built from its space group: FCC `Fm-3m` (4a), or L1₂ `Pm-3m` (B at 1a, A at 3c). It is then turned so that [111] is along z: a stack of close-packed (111) layers, d₁₁₁ = a_cubic/√3 apart, in the same a × √3a form as the MXene with a → a₁₁₁ (a₁₁₁ × √3a₁₁₁ for FCC, 2a₁₁₁ × 2√3a₁₁₁ for L1₂, whose B atoms form a 2 × 2 superlattice in each layer).
3. **Supercell.** Along x and along y, the tool finds the smallest whole number of MXene replicas that matches a whole number of metal replicas within the strain tolerance. The metal takes the strain (see [Supercell matching](#supercell-matching)).
4. **Placement.** Metal (111) on MXene (0001), metal [1-10] along x, layer by layer: the first metal (111) plane sits at a set gap above the top MXene layer, on site A, B or C of the MXene grid, and the metal layers then follow A → B → C (or its twin). An optional in-plane shift moves the metal off the site.
5. **Relaxation.** MACE (default MACE-MH-1, `omat_pbe` head) in float64, optional D3(BJ) dispersion and FIRE. The outermost metal layers are held fixed, and the in-plane cell can relax too. Optionally reports the work of adhesion.
6. **Built vs relaxed.** Per-atom displacements, layer heights and spacings, the gap, the MXene thickness and the cell, before and after relaxation (see below).
7. **Atomic planes.** Shows the (hkl) planes of the MXene cell or of the strained metal cell in the viewer.

Lattice constants come from `data/MXene_FCC_L12_strain_matrix.xlsx` (41 ordered MXenes, 23 metals). Solid solutions, Ti₃C₂Tₓ and Ni₃(Al,Ti) are left out.

Fixed geometry rules for the starting structure:
- M–X layer spacing is 0.385 a.
- The O height comes from a 2.0 Å M–O bond.
- AβC puts the central layer on A. AβA puts the metals on A and C/N on B.
- O sits on the site that neither of its two nearest layers uses.
- With `metal_site` = `auto`, the first metal layer takes the same rule: the site free of the top two MXene layers, so the close-packed stacking carries on into the metal.

## Supercell matching

The MXene basal plane and the metal (111) plane are both hexagonal nets, written as rectangular cells: a × √3a for the MXene and a₁₁₁ × √3a₁₁₁ for the metal (`a_111 (A)` in the spreadsheet; twice that for L1₂). The two in-plane directions are matched separately. For each direction, with cell lengths L_MX (a along x, √3a along y) and L_M (a₁₁₁ or √3a₁₁₁):

1. Try m = 1, 2, … up to `max_replicas` MXene replicas.
2. For each m, take the whole number of metal replicas closest to the same length: n = round(m · L_MX / L_M), at least 1.
3. The metal strain is ε = m · L_MX / (n · L_M) − 1.
4. Keep the first m with |ε| ≤ `tol`. This gives the smallest supercell that fits the tolerance.

The box is then m_x · a by m_y · √3a. The MXene keeps its lattice constant, and the metal is strained by ε_x and ε_y to fill the box exactly. Both cells have the same √3 aspect ratio, so x and y give the same m, n and ε: the metal (111) net takes an equi-biaxial strain, and at m = n = 1 it is the value in the spreadsheet's `Strain` sheet. If no m up to `max_replicas` works, the build stops with an error; raise `tol` or `max_replicas`.

For example:
- Au (a₁₁₁ = 2.884 Å) on V₂C (a = 2.881 Å): 1 × 1, ε = −0.09 %.
- Au on Ti₃C₂ (a = 3.104 Å): 1 × 1 would be +7.6 %. The first coincidence within 2 % is 11 MXene cells on 12 Au cells (ε = −1.3 %), so it needs `max_replicas` ≥ 11; 13 : 14 gives −0.05 %.

The *Box* panel in the web app and `outputs.json → match` report m, n, the box length and ε for each direction.

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
| `mxene` | `V2C` | Formula, a row of the MXenes sheet |
| `stacking` | `ABC` | `ABC` = AβC (octahedral), `ABA` = AβA (trigonal prismatic) |
| `vacuum` | 20.0 | Å, below the MXene only |
| `metal` | `Au` | Symbol, a row of the Metals sheet |
| `metal_layers` | 4 | metal (111) atomic planes |
| `metal_site` | `auto` | site of the first metal layer on the MXene grid: `A`, `B`, `C`, or `auto` (free of the top two MXene layers) |
| `metal_stacking` | `ABC` | `ABC` continues the MXene's A → B → C sense upward, `ACB` is the twin |
| `gap` | 2.3 | Å, top MXene layer to first metal (111) layer |
| `shift` | [0, 0] | Å, in-plane shift of the metal along x, y, from its site |
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
- `metal`: spreadsheet row (with a₁₁₁), space group, d₁₁₁, orientation, layer sites, the (111) cell before and after strain and the strained cubic cell vectors.
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

(hkl) refer to the MXene a × √3a × c cell or to the metal's strained cubic cell, whose [111] is along z, so metal (111) planes are the horizontal layers. The MXene c includes the vacuum, so MXene planes with l ≠ 0 tilt with it.

## Files

| File | Content |
|---|---|
| `builder.py` | builds the interface from the inputs |
| `relax.py` | MACE relaxation, work of adhesion and the built-vs-relaxed comparison |
| `planes.py` | atomic planes for the viewer |
| `app.py` | Flask app |
| `templates/`, `static/` | page, styles, script and bundled 3Dmol.js |
| `data/` | lattice-constant spreadsheet |
