"""Atomic planes (hkl) of the MXene or of the metal inside the interface.

Indices refer to that crystal's conventional orthogonal cell: a × √3a × c for the MXene
(its c includes the vacuum, so planes with l ≠ 0 tilt with it) and the strained cubic cell
for the metal. An atomic plane is a set of atoms with the same projection on the plane normal;
planes of one family repeat every d_hkl.
"""
import itertools
import math

import numpy as np
from pymatgen.core import Lattice


def miller(indices):
    indices = [int(v) for v in indices]
    if len(indices) != 3 or not any(indices):
        raise ValueError("Give three indices, not all zero")
    g = math.gcd(*indices)
    return [v // g for v in indices]


def atomic_planes(structure, part, cell, hkl, tol=0.05):
    """Distinct atomic planes of one part ("mxene" or "metal") within one period d."""
    g = np.array(hkl, dtype=float) @ Lattice.orthorhombic(*cell).reciprocal_lattice_crystallographic.matrix
    d = 1.0 / np.linalg.norm(g)
    n = g * d
    proj = (structure.cart_coords @ n) % d
    planes = []
    for k in sorted((k for k, s in enumerate(structure) if s.properties["part"] == part), key=lambda k: proj[k]):
        p = proj[k]
        match = next((pl for pl in planes if min(abs(p - pl["position"]), d - abs(p - pl["position"])) < tol), None)
        if match is None:
            match = dict(position=float(p), atoms=[])
            planes.append(match)
        match["atoms"].append(k)
    for pl in planes:
        species = [structure[k].species_string for k in pl["atoms"]]
        pl["species"] = {e: species.count(e) for e in dict.fromkeys(species)}
        pl["layers"] = sorted({int(structure[k].properties["layer"]) for k in pl["atoms"]})
        pl["position"] = round(pl["position"], 4)
    return dict(hkl=hkl, part=part, normal=(np.round(n, 6) + 0.0).tolist(), d=round(float(d), 5), planes=planes)


def polygon(n, offset, origin, edges):
    """Corners of the plane n·r = offset inside the box origin + [0,1]³·edges, ordered around the outline."""
    corners = [origin + np.array(f) @ edges for f in itertools.product((0, 1), repeat=3)]
    points = []
    for a, b in itertools.combinations(range(8), 2):
        if bin(a ^ b).count("1") != 1:
            continue
        sa, sb = corners[a] @ n - offset, corners[b] @ n - offset
        if abs(sa) < 1e-9:
            points.append(corners[a])
        elif sa * sb < 0:
            points.append(corners[a] + (corners[b] - corners[a]) * sa / (sa - sb))
    unique = []
    for q in points:
        if all(np.linalg.norm(q - u) > 1e-6 for u in unique):
            unique.append(q)
    if len(unique) < 3:
        return None
    centre = np.mean(unique, axis=0)
    u = unique[0] - centre
    u /= np.linalg.norm(u)
    w = np.cross(n, u)
    unique.sort(key=lambda q: math.atan2((q - centre) @ w, (q - centre) @ u))
    return [np.round(q, 4).tolist() for q in unique]


def plane_polygons(structure, result, repeat=(1, 1), pad=1.0, limit=150):
    """Polygons for every plane of the family crossing the part's region: repeat[0] × repeat[1] cells
    in-plane and the part's z range plus pad."""
    n, d = np.array(result["normal"]), result["d"]
    z = np.array([s.coords[2] for s in structure if s.properties["part"] == result["part"]])
    m = structure.lattice.matrix
    origin = np.array([0.0, 0.0, z.min() - pad])
    edges = np.array([m[0] * repeat[0], m[1] * repeat[1], [0.0, 0.0, np.ptp(z) + 2 * pad]])
    corners = [origin + np.array(f) @ edges for f in itertools.product((0, 1), repeat=3)]
    lo, hi = min(c @ n for c in corners), max(c @ n for c in corners)
    shapes = []
    for index, plane in enumerate(result["planes"]):
        for k in range(math.floor((lo - plane["position"]) / d), math.ceil((hi - plane["position"]) / d) + 1):
            outline = polygon(n, plane["position"] + k * d, origin, edges)
            if outline:
                shapes.append(dict(plane=index, vertices=outline))
            if len(shapes) >= limit:
                return shapes
    return shapes
