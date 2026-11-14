#!/usr/bin/env python3
"""Fetch the Gaskell shape model of (433) Eros from the NASA Planetary Data System.

    python3 examples/python/fetch_eros_shape.py --outdir data/asteroids/eros_gaskell

The shape model was derived from images of the NEAR Shoemaker multispectral imager by
R. Gaskell. The Small Bodies Node of the Planetary Data System serves it at four resolutions
(Q = 64, 128, 256 and 512). The script downloads the smallest, Q = 64 with 49,152 triangular
plates (4.2 MB), and the planetary constants file of the model (a few hundred bytes), and
checks the size and the SHA-256 checksum of both. The table lists the vertices of the model in
the order of its implicitly connected quadrilateral grid, which makes a coarser model by
keeping every second vertex, as the documentation of the model describes. ``--q 32`` (the
default) does that and gives 12,288 plates, a mesh that PTM3D can spin up for 500 rotations in
about 10 GB of memory. ``--q 64`` keeps all 49,152 plates. The script writes a Wavefront OBJ
file in kilometers, with the quadrilaterals cut along their shorter diagonals, and checks the
mesh: it must be closed, with every edge shared by two plates, and the plates must be wound
consistently and point outward. It writes ``shape.obj`` and ``eros_alex.tpc`` to --outdir and prints the pole, the prime meridian
and the period of the model in the form that ``generic_body.spin`` takes.

The shape file is in the body-fixed frame of the constants file. Its pole is
RA 11.363, Dec 17.232 degrees and its prime meridian is W = 326.08 + 1639.38928 d degrees, with
d the days from J2000. The constants file measures W from the ascending node of the equator on
the celestial equator. PTM3D measures it from the J2000 ecliptic x axis projected onto the
equator, or from the ecliptic y axis when the pole is within 26 degrees of the x axis, as it is
for Eros. ``prime_meridian_deg`` in the configuration is therefore W at J2000 plus the angle
between those axes (-19.943 degrees for Eros). The script prints it.

Data: Gaskell, R., Gaskell Eros Shape Model V1.1, NASA Planetary Data System,
doi:10.26033/d0gq-9427 (urn:nasa:pds:gaskell.ast-eros.shape-model). Public domain (CC0).
Run it from the repository root. ``data/asteroids`` is not tracked by git.
"""
import argparse
import hashlib
import math
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np

BASE = 'https://sbnarchive.psi.edu/pds4/non_mission/gaskell.ast-eros.shape-model_V1_1'
FILES = {   # name: (path on the server, size in bytes, SHA-256)
    'ver64q.tab': ('data/vertex/ver64q.tab', 4246671,
                   'de97008ff2c865039a51aa29ec13441ac3c87bd5cb27879091cb488c74fcb6dd'),
    'eros_alex.tpc': ('document/eros_alex.tpc.txt', 1353, None),
}


def download(name, outdir):
    path, size, digest = FILES[name]
    target = Path(outdir) / name
    if not target.is_file() or target.stat().st_size != size:
        request = urllib.request.Request(f'{BASE}/{path}', headers={'User-Agent': 'ptm3d-fetch'})
        with urllib.request.urlopen(request, timeout=120) as response:
            target.write_bytes(response.read())
    data = target.read_bytes()
    if len(data) != size:
        sys.exit(f'{target}: {len(data)} bytes, expected {size}. Delete it and run again.')
    found = hashlib.sha256(data).hexdigest()
    if digest and found != digest:
        sys.exit(f'{target}: SHA-256 {found} differs from the recorded {digest}.')
    if not digest:
        print(f'{name}: {size} bytes, SHA-256 {found}')
    return target


def read_table(path):
    lines = Path(path).read_text().splitlines()
    n_vert, n_face = (int(x) for x in lines[0].split())
    vertices = np.array([[float(x) for x in ln.split()[1:4]] for ln in lines[1:1 + n_vert]])
    faces = np.array([[int(x) - 1 for x in ln.split()[1:4]] for ln in lines[1 + n_vert:1 + n_vert + n_face]])
    if len(vertices) != n_vert or len(faces) != n_face:
        sys.exit(f'{path} does not hold the {n_vert} vertices and {n_face} plates of its header')
    return vertices, faces


def coarsen(vertices, q_full, q):
    """Vertices and plates of the model at grid size q from the table at q_full.

    The table holds the 6 (q_full + 1)^2 grid vertices of the six cube faces in the order
    face, j, i. Vertices on the edges of the cube appear more than once and are merged by
    position. Every quadrilateral is cut along its shorter diagonal.
    """
    step = q_full // q
    grid = vertices.reshape(6, q_full + 1, q_full + 1, 3)[:, ::step, ::step]
    flat = grid.reshape(-1, 3)
    # The copies of an edge vertex differ by up to 3 m, so merge vertices closer than 10 m.
    cell, tol = 0.02, 0.01
    seen, ids_flat, points = {}, np.empty(len(flat), dtype=int), []
    for n, p in enumerate(flat):
        key = tuple(np.floor(p / cell).astype(int))
        found = None
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for m in seen.get((key[0] + dx, key[1] + dy, key[2] + dz), ()):
                        if np.linalg.norm(points[m] - p) < tol:
                            found = m
        if found is None:
            found = len(points)
            points.append(p)
            seen.setdefault(key, []).append(found)
        ids_flat[n] = found
    ids = ids_flat.reshape(6, q + 1, q + 1)
    points = np.array(points)
    faces = []
    for f in range(6):
        for j in range(q):
            for i in range(q):
                a, b, c, d = ids[f, j, i], ids[f, j + 1, i], ids[f, j + 1, i + 1], ids[f, j, i + 1]
                pa, pb, pc, pd = points[[a, b, c, d]]
                if np.sum(np.cross(pa, pc) ** 2) <= np.sum(np.cross(pd, pb) ** 2):
                    faces += [(a, c, d), (a, b, c)]
                else:
                    faces += [(a, b, d), (d, b, c)]
    return points, np.array(faces)


def check_mesh(vertices, faces):
    """Drop unused vertices, then require a closed, consistently wound, outward mesh."""
    used = np.unique(faces)
    remap = -np.ones(len(vertices), dtype=int)
    remap[used] = np.arange(len(used))
    vertices, faces = vertices[used], remap[faces]
    edges = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    directed = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    n_directed = len({(int(a), int(b)) for a, b in directed})
    if not np.all(counts == 2) or n_directed != len(directed):
        sys.exit('the mesh is not closed or its plates are not wound consistently')
    a, b, c = (vertices[faces[:, i]] for i in range(3))
    volume = float(np.sum(np.einsum('ij,ij->i', a, np.cross(b, c))) / 6.0)
    if volume <= 0.0:
        faces = faces[:, ::-1]
        volume = -volume
        print('the plates were wound inward, so the winding was reversed')
    return vertices, faces, volume


def spin_from_constants(path):
    text = Path(path).read_text()
    def values(key):
        match = re.search(rf'{key}\s*=\s*\(([^)]*)\)', text)
        return [float(x) for x in match.group(1).split()]
    ra, dec = values('BODY2000433_POLE_RA')[0], values('BODY2000433_POLE_DEC')[0]
    pm0, rate = values('BODY2000433_PM')[:2]
    a, d = math.radians(ra), math.radians(dec)
    pole = np.array([math.cos(d) * math.cos(a), math.cos(d) * math.sin(a), math.sin(d)])
    obliquity = math.radians(23.4392911)                      # the value PTM3D uses for the ecliptic
    seed = np.array([1.0, 0.0, 0.0])
    if abs(pole[0]) > 0.9:                                    # a pole near the equinox
        seed = np.array([0.0, math.cos(obliquity), math.sin(obliquity)])
    e1 = seed - (seed @ pole) * pole                          # the axis PTM3D measures W from
    e1 /= np.linalg.norm(e1)
    node = np.array([-math.sin(a), math.cos(a), 0.0])        # the axis the constants file measures W from
    delta = math.degrees(math.atan2(pole @ np.cross(e1, node), e1 @ node))
    return dict(pole_ra_deg=ra, pole_dec_deg=dec, prime_meridian_deg=(pm0 + delta) % 360.0,
                period_hours=24.0 * 360.0 / rate, delta_deg=delta)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', default='data/asteroids/eros_gaskell')
    ap.add_argument('--q', type=int, default=32, choices=(8, 16, 32, 64),
                    help='grid size of the model: 12 q^2 plates (default 32, 12,288 plates)')
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    table = download('ver64q.tab', out)
    constants = download('eros_alex.tpc', out)
    vertices, faces = read_table(table)
    if args.q != 64:
        vertices, faces = coarsen(vertices, 64, args.q)
    vertices, faces, volume = check_mesh(vertices, faces)
    radius = (3.0 * volume / (4.0 * math.pi)) ** (1.0 / 3.0)
    with open(out / 'shape.obj', 'w') as fh:
        for x, y, z in vertices:
            fh.write(f'v {x:.5f} {y:.5f} {z:.5f}\n')
        for i, j, k in faces:
            fh.write(f'f {i + 1} {j + 1} {k + 1}\n')
    print(f'shape.obj: {len(vertices)} vertices, {len(faces)} plates, closed, outward, '
          f'volume {volume:.1f} km^3, equal-volume radius {radius:.3f} km')
    spin = spin_from_constants(constants)
    print(f"spin: pole RA {spin['pole_ra_deg']:.3f}, Dec {spin['pole_dec_deg']:.3f}, period "
          f"{spin['period_hours']:.8f} h, angle from the PTM3D axis to the ascending node "
          f"{spin['delta_deg']:.4f}")
    print(f"      generic_body.spin: prime_meridian_deg {spin['prime_meridian_deg']:.6f} at "
          f"epoch 'JD 2451545.000000'")


if __name__ == '__main__':
    sys.exit(main())
