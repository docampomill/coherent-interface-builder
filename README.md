# MXene/metal interface builder

`mxene_au_box.py` as a web app. The MXene and the metal are each built from their own crystal system with pymatgen, then matched in one orthogonal box:

1. **MXene**: conventional a × √3a × c cell, two copies of each layer, vacuum only below (top layer at z = c).
2. **Metal**: conventional cubic cell from its space group, FCC `Fm-3m` (4a) or L1₂ `Pm-3m` (B at 1a, A at 3c).
3. **Box**: along x and along y, the smallest m MXene replicas ≈ n metal replicas with |strain| ≤ tolerance. The metal takes the strain.
4. **Location**: metal (001) on MXene (0001), metal [100] along x, at a gap above the top layer, shifted in-plane.
5. **Relaxation** (`relax.py`): MACE (default MACE-MH-1, omat_pbe head) in float64, optional D3(BJ), FIRE, top metal layers fixed, optional in-plane cell relaxation.

Lattice constants come from `data/MXene_FCC_L12_strain_matrix.xlsx`: 41 ordered MXenes and 23 metals. Solid solutions, Ti3C2Tx and Ni3(Al,Ti) are left out.

## Run

```bash
podman build -t mxbuilder:0.8.0 .
podman run --rm -p 8000:8000 mxbuilder:0.8.0          # http://localhost:8000
python builder.py inputs.json                          # build only, without the web app
python relax.py inputs.json                            # build and relax, e.g. on a cluster
```

## Inputs (`inputs.json`)

| Name | Default | Meaning |
|---|---|---|
| `mxene` | `Ti3C2` | MXenes sheet, Formula |
| `stacking` | `ABC` | `ABC` = AβC (octahedral), `ABA` = AβA (trigonal prismatic) |
| `vacuum` | 20.0 | Å, below the MXene only |
| `metal` | `Au` | Metals sheet, Symbol |
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

Fixed rules: M–X layer spacing 0.385 a; O height from a 2.0 Å M–O bond; AβC puts the central layer on A, AβA puts metals on A and C/N on B; O sits on the site free of its two nearest layers.

## Outputs (`outputs.json`)

`inputs`, `mxene` (row, layers, sites, cell, thickness), `metal` (row, space group, strained cell), `match` (replicas, box length and metal strain along x and y), `interface` (formula, atoms, box, layers), `spreadsheet` (file, SHA-256) and `files`: `POSCAR` (interface), `POSCAR_<mxene>`, `POSCAR_<metal>`.

After a relaxation, `relaxation` (optimizer, fixed atoms, converged, steps, initial and final energy, final fmax, cell a, b, γ, MXene–metal gap, package versions, energy/fmax per step, and the work of adhesion in J/m² when asked) and `files.POSCAR_relaxed`.

Loading `outputs.json` in the app restores its inputs.

## Atomic planes

(hkl) refer to the MXene a × √3a × c cell or to the metal's strained cubic cell. The MXene c includes the vacuum, so MXene planes with l ≠ 0 tilt with it.
