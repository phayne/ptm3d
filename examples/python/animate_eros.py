#!/usr/bin/env python3
"""Animate (433) Eros over one rotation: the body, and the temperature at one site.

The animation on the front page of the repository comes from this script. Three panels follow
one rotation in 96 steps. (a) is the shape model colored by surface temperature, seen from a
fixed direction, with the Sun fixed in the picture so that the body turns under it. (b) is the
surface temperature of one face against time, and (c) is the temperature of that face against
depth, with the range that each depth passes through in the rotation. Everything plotted is
read from the output of the run. The script draws the body with its own z-buffer, so a
concave shape such as the saddle of Eros is ordered correctly.

The shape is the Gaskell model of Eros from NEAR images (NASA Planetary Data System,
doi:10.26033/d0gq-9427). Commands, from the repository root, after building ``ptm_thermal_sim``
(INSTALL.md) and fetching the kernels (``python3 scripts/fetch_kernels.py --bodies moon``)::

    PYTHONPATH=src/python python3 examples/python/fetch_asteroid.py 433 --shape sphere \\
        --start 1993-01-01 --stop 1999-01-01 --epoch 1998-12-23T15:41:23 \\
        --outdir data/asteroids/eros_gaskell
    python3 examples/python/fetch_eros_shape.py --outdir data/asteroids/eros_gaskell
    ./build/bin/ptm_thermal_sim --config examples/asteroid_eros_animation.yaml
    PYTHONPATH=src/python python3 examples/python/animate_eros.py output/asteroid_eros_animation \\
        --out docs/figures/eros_diurnal_cycle.gif

The first command downloads the orbit of Eros from JPL Horizons, the second the shape from the
Planetary Data System (4 MB, checked against a recorded checksum). The script needs numpy,
matplotlib and Pillow (the first command also spiceypy). The run takes about 5 minutes and
about 10 GB of memory (500 rotations of spin-up, then three rotations in 288 outputs) and the
animation a few minutes more. The second rotation is drawn. The first is discarded because the
model settles slightly in its first rotation after the spin-up, and the output that opens the
third is used to report how well the second closes on itself. ``--mp4`` and ``--png`` also
write a video and a single frame. Encoding to GIF and MP4 uses ``ffmpeg`` when it is installed
and Pillow otherwise. ``--face`` chooses another site, and ``--list-sites`` prints the faces
that the automatic choice ranks.

The picture uses lighting to show the shape. Temperature is interpolated for display only: each
vertex of the shape gets the area-weighted mean of the temperatures of the plates around it, the
rasterizer interpolates that temperature linearly across each plate, and the color map is applied
to the interpolated temperature at every pixel. Without this the plates, which are two triangles
of slightly different tilt per quad of the shape grid, show as speckle along the terminator.
Panels (b) and (c) are the exact values of the face of the site. A pixel is shaded by the cosine
of the angle between the Sun and the smoothed surface normal (the area-weighted mean of the
normals of the plates around each vertex), multiplied by the share of lit plates around the vertex,
a plate being lit when the run's direct solar flux on it is not zero, and both are interpolated
across the plate. A pixel that is dark at night or inside a shadow of the terrain gets only the
ambient part, and the edge of a shadow is softened by about one plate. The color of a pixel is its
temperature on the color bar times ``ambient + (1 - ambient) * diffuse``, so a pixel that looks
straight into the Sun shows exactly the color of the color bar and any other pixel a darker version
of it. The script prints how far the display is from the face values.
"""
import argparse
import json
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

from ptm import Results, plots
from ptm.results import depth_entries

SOLAR_CONSTANT = 1361.0          # W/m2 at 1 AU, as in the engine
SITE_LAT_RANGE = (20.0, 45.0)    # automatic site: planetocentric latitude between these (either hemisphere)
AMBIENT = 0.60                   # share of the color that no lighting removes
TMIN, TMAX = 50.0, 270.0         # color limits of the picture and the color bar (K)
CAMERA_ELEVATION = 10.0          # degrees above the equator plane
CAMERA_AZIMUTH = -25.0           # degrees from the Sun, about the spin axis
SUN_AZIMUTH = 0.0
SUPERSAMPLE = 3
SITE_COLOR = '#00b7d9'
FIGURE_IN = (10.0, 6.0)          # width, height of the figure in inches
PANEL_A_IN = 5.55                # width of the picture of the body, inches
SUN_BAND_KM = 3.4                # free space above the body for the Sun arrow
FILL = 0.94                      # the widest view of the body fills this share of the picture


# ------------------------------------------------------------------------------ the run

class Cycle:
    """The arrays of one run, for the second of its three rotations."""

    def __init__(self, run_dir, face=None):
        self.run = Results(run_dir, mesh_units='km')
        run = self.run
        mesh = run.mesh
        self.vertices = np.asarray(mesh['vertices'], float)
        self.faces = np.asarray(mesh['faces'], int)
        self.n_faces = len(self.faces)
        a, b, c = (self.vertices[self.faces[:, i]] for i in range(3))
        cross = np.cross(b - a, c - a)
        self.area = 0.5 * np.linalg.norm(cross, axis=1)
        self.normal = cross / (2.0 * self.area[:, None])
        self.center = (a + b + c) / 3.0
        # The engine turns a normal around when it points toward the center of the body
        # (center . normal < 0), which includes some faces on the steep walls of a concavity.
        # The lit, shadow and night classes follow the engine, the picture the true outward normal.
        self.engine_flipped = np.sum(self.normal * self.center, axis=1) < 0
        self.engine_normal = np.where(self.engine_flipped[:, None], -self.normal, self.normal)
        self.lat = np.degrees(np.arcsin(self.center[:, 2] / np.linalg.norm(self.center, axis=1)))
        self.lon = np.degrees(np.arctan2(self.center[:, 1], self.center[:, 0]))

        elapsed_h = run.times().elapsed_s / 3600.0
        total = len(elapsed_h)
        if total % 3:
            raise SystemExit('expected three rotations in a multiple of three outputs')
        self.n = total // 3                           # outputs per rotation
        self.window = slice(self.n, 2 * self.n + 1)   # the second rotation and the output that closes it
        self.period_h = float(elapsed_h[self.n])
        self.t_h = elapsed_h[:self.n]

        sun = np.array([[p['sun_x'], p['sun_y'], p['sun_z']] for p in run.body_positions()], float)
        self.sun_distance_au = float(np.mean(np.linalg.norm(sun, axis=1))) / 1.495978707e8
        self.sun = sun / np.linalg.norm(sun, axis=1)[:, None]     # body-fixed, per output
        cos = self.engine_normal @ self.sun.T                     # (faces, outputs)
        flux = np.array(run.stack('direct_solar')).T              # (faces, outputs) W/m2
        self.temperature = np.array(run.stack('temperature')).T   # (faces, outputs) K
        expected = SOLAR_CONSTANT / self.sun_distance_au ** 2
        lit = flux > 0.05 * expected * np.maximum(cos, 1e-9)
        shadow = (cos > 0.0) & ~lit                               # terrain shadow, Sun up
        night = cos <= 0.0
        w = self.window
        self.sun, self.cos, self.flux, self.temperature = self.sun[w], cos[:, w], flux[:, w], self.temperature[:, w]
        self.lit, self.shadow, self.night = lit[:, w], shadow[:, w], night[:, w]

        stats = Path(run_dir) / 'statistics'
        manifest = json.loads((stats / 'depth_map_z.json').read_text())
        self.skin_depth_m = float(manifest['skin_depth_m'])
        entries = depth_entries(stats)
        names = sorted(entries, key=lambda f: entries[f]['resolved_m'])
        if not names or not all(entries[f].get('on_node') for f in names):
            raise SystemExit('every depth frame must be on a solver node (no interpolation)')
        self.depth_m = np.array([0.0] + [entries[f]['resolved_m'] for f in names])
        self.depth_fields = ['temperature'] + names
        self._neighbors = None
        self.face = self.choose_face() if face is None else int(face)

    def neighbors(self, face):
        """Faces that share a vertex with ``face``."""
        if self._neighbors is None:
            owners = [[] for _ in range(len(self.vertices))]
            for f, tri in enumerate(self.faces):
                for v in tri:
                    owners[v].append(f)
            self._neighbors = owners
        out = {g for v in self.faces[face] for g in self._neighbors[v]}
        out.discard(face)
        return sorted(out)

    def rank_sites(self):
        """Faces with a night, a day and the longest terrain shadow with the Sun up.

        A face counts only if the faces around it are shadowed for nearly as long, so that a
        single small face on a crater wall is not chosen. The score is the smaller of the shadow
        durations of the face and of the median of its neighbors.
        """
        n = self.n
        count = self.shadow[:, :n].sum(axis=1)
        rows = []
        for i in range(self.n_faces):
            if not SITE_LAT_RANGE[0] <= abs(self.lat[i]) <= SITE_LAT_RANGE[1] or self.engine_flipped[i]:
                continue
            if self.night[i, :n].sum() < 10 or self.lit[i, :n].sum() < 10:
                continue
            around = [g for g in self.neighbors(i) if not self.engine_flipped[g]]
            rows.append((int(min(count[i], np.median(count[around]))), int(count[i]), i))
        rows.sort(reverse=True)
        return rows

    def choose_face(self):
        return self.rank_sites()[0][2]

    def site_profiles(self):
        """Temperature (levels, outputs) of the chosen face in the window."""
        rows = [np.asarray(self.run.stack(f, faces=[self.face]))[self.window, 0] for f in self.depth_fields]
        return np.array(rows)


# -------------------------------------------------------------------------- the picture

def camera_basis():
    """Right, up and toward-camera unit vectors in the frame of the picture (z along the spin)."""
    e, a = math.radians(CAMERA_ELEVATION), math.radians(SUN_AZIMUTH + CAMERA_AZIMUTH)
    toward = np.array([math.cos(e) * math.cos(a), math.cos(e) * math.sin(a), math.sin(e)])
    right = np.cross([0.0, 0.0, 1.0], toward)
    right /= np.linalg.norm(right)
    up = np.cross(toward, right)
    return right, up, toward


def spin_rotation(sun_body):
    """Rotation about the spin axis that puts the Sun at the fixed azimuth of the picture."""
    alpha = math.radians(SUN_AZIMUTH) - math.atan2(sun_body[1], sun_body[0])
    c, s = math.cos(alpha), math.sin(alpha)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def view_extents(cycle):
    """Smallest box (x0, x1, y0, y1 in km) that holds the body at every output of the cycle."""
    right, up, _ = camera_basis()
    x0 = y0 = np.inf
    x1 = y1 = -np.inf
    for k in range(cycle.n):
        v = cycle.vertices @ spin_rotation(cycle.sun[k]).T
        x, y = v @ right, v @ up
        x0, x1, y0, y1 = min(x0, x.min()), max(x1, x.max()), min(y0, y.min()), max(y1, y.max())
    return x0, x1, y0, y1


def vertex_normals(vertices, faces, face_normal, face_area):
    """Area weighted vertex normals."""
    out = np.zeros_like(vertices)
    for i in range(3):
        np.add.at(out, faces[:, i], face_normal * face_area[:, None])
    return out / np.linalg.norm(out, axis=1)[:, None]


def rasterize(px, depth, faces, width, height):
    """z-buffer: the nearest face at every pixel. px is (V, 2) in pixels, depth grows toward us."""
    zbuf = np.full((height, width), -np.inf)
    face_id = np.full((height, width), -1, dtype=np.int32)
    x, y = px[:, 0], px[:, 1]
    for f, (i, j, k) in enumerate(faces):
        det = (y[j] - y[k]) * (x[i] - x[k]) + (x[k] - x[j]) * (y[i] - y[k])
        if det >= 0.0:                      # clockwise on the screen: faces away from the camera
            continue
        x0 = max(int(math.floor(min(x[i], x[j], x[k]))), 0)
        x1 = min(int(math.ceil(max(x[i], x[j], x[k]))), width - 1)
        y0 = max(int(math.floor(min(y[i], y[j], y[k]))), 0)
        y1 = min(int(math.ceil(max(y[i], y[j], y[k]))), height - 1)
        if x1 < x0 or y1 < y0:
            continue
        gx, gy = np.meshgrid(np.arange(x0, x1 + 1) + 0.5, np.arange(y0, y1 + 1) + 0.5)
        l1 = ((y[j] - y[k]) * (gx - x[k]) + (x[k] - x[j]) * (gy - y[k])) / det
        l2 = ((y[k] - y[i]) * (gx - x[k]) + (x[i] - x[k]) * (gy - y[k])) / det
        l3 = 1.0 - l1 - l2
        inside = (l1 >= -1e-9) & (l2 >= -1e-9) & (l3 >= -1e-9)
        if not inside.any():
            continue
        z = l1 * depth[i] + l2 * depth[j] + l3 * depth[k]
        sub_z = zbuf[y0:y1 + 1, x0:x1 + 1]
        win = inside & (z > sub_z)
        sub_z[win] = z[win]
        face_id[y0:y1 + 1, x0:x1 + 1][win] = f
    return zbuf, face_id


class Renderer:
    """Draws the shaded, colored body for one output, in a fixed frame."""

    def __init__(self, cycle, width, height, extents):
        from matplotlib import colormaps
        self.c = cycle
        self.w, self.h, self.ss = width, height, SUPERSAMPLE
        self.right, self.up, self.toward = camera_basis()
        self.vnormal = vertex_normals(cycle.vertices, cycle.faces, cycle.normal, cycle.area)
        self.vertex_area = np.zeros(len(cycle.vertices))
        for i in range(3):
            np.add.at(self.vertex_area, cycle.faces[:, i], cycle.area)
        self.cmap = colormaps['inferno']
        x0, x1, y0, y1 = extents
        self.scale = self.ss * width / (x1 - x0)         # pixels per km, fixed for the whole cycle
        self.cx, self.cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)

    def project(self, v_view):
        x = v_view @ self.right
        y = v_view @ self.up
        px = np.stack([self.w * self.ss / 2.0 + (x - self.cx) * self.scale,
                       self.h * self.ss / 2.0 - (y - self.cy) * self.scale], axis=1)
        return px, v_view @ self.toward

    def vertex_mean(self, values):
        """Area-weighted mean of a per-face quantity at every vertex."""
        c = self.c
        acc = np.zeros(len(c.vertices))
        for i in range(3):
            np.add.at(acc, c.faces[:, i], c.area * values)
        return acc / self.vertex_area

    def render(self, k):
        c = self.c
        rot = spin_rotation(c.sun[k])
        v = c.vertices @ rot.T
        sun_view = rot @ c.sun[k]
        px, depth = self.project(v)
        W, H = self.w * self.ss, self.h * self.ss
        zbuf, fid = rasterize(px, depth, c.faces, W, H)
        covered = fid >= 0
        ids = fid[covered]                                          # faces of the covered pixels
        tri = c.faces[ids]
        gy, gx = np.nonzero(covered)
        gx, gy = gx + 0.5, gy + 0.5
        xs, ys = px[tri, 0], px[tri, 1]
        det = (ys[:, 1] - ys[:, 2]) * (xs[:, 0] - xs[:, 2]) + (xs[:, 2] - xs[:, 1]) * (ys[:, 0] - ys[:, 2])
        det = np.where(det == 0.0, 1.0, det)
        l1 = ((ys[:, 1] - ys[:, 2]) * (gx - xs[:, 2]) + (xs[:, 2] - xs[:, 1]) * (gy - ys[:, 2])) / det
        l2 = ((ys[:, 2] - ys[:, 0]) * (gx - xs[:, 2]) + (xs[:, 0] - xs[:, 2]) * (gy - ys[:, 2])) / det
        weights = np.stack([l1, l2, 1.0 - l1 - l2], axis=1)

        def across(vertex_values):                                  # linear interpolation inside the plate
            return np.sum(weights * vertex_values[tri], axis=1)

        t_vertex = self.vertex_mean(c.temperature[:, k].astype(float))
        t_pixel = across(t_vertex)
        diffuse = np.clip(self.vnormal @ c.sun[k], 0.0, 1.0) * self.vertex_mean(c.lit[:, k].astype(float))
        factor = AMBIENT + (1.0 - AMBIENT) * np.clip(across(diffuse), 0.0, 1.0)
        base = self.cmap(np.clip((t_pixel - TMIN) / (TMAX - TMIN), 0, 1))[:, :3]
        rgb = np.ones((H, W, 3))
        rgb[covered] = base * factor[:, None]
        t_image = np.full((H, W), np.nan)
        t_image[covered] = t_pixel
        rgb = rgb.reshape(self.h, self.ss, self.w, self.ss, 3).mean(axis=(1, 3))
        shift = t_pixel - c.temperature[ids, k]                     # displayed minus face value
        info = dict(px=px, depth=depth, zbuf=zbuf, sun_view=sun_view, rot=rot, t_image=t_image,
                    t_vertex=t_vertex, shift=shift)
        return np.clip(rgb, 0.0, 1.0), info

    def site(self, k, info):
        """Screen position of the site in axes pixels, and whether it is in view."""
        c = self.c
        p = (c.center[c.face] + 0.03 * c.normal[c.face]) @ info['rot'].T
        px, depth = self.project(p[None, :])
        ix, iy = int(px[0, 0]), int(px[0, 1])
        W, H = self.w * self.ss, self.h * self.ss
        visible = 0 <= ix < W and 0 <= iy < H and depth[0] >= info['zbuf'][iy, ix] - 0.05
        return px[0, 0] / self.ss, px[0, 1] / self.ss, bool(visible)


def color_fidelity(cycle, renderer):
    """How far the shaded color of a lit face is from its color-bar color, as a temperature.

    For every face that is lit at an output and faces the Sun within 60 degrees, the shaded
    color is mapped back to the nearest temperature on the color bar.
    """
    ramp = np.linspace(TMIN, TMAX, 2201)
    colors = renderer.cmap((ramp - TMIN) / (TMAX - TMIN))[:, :3]
    out, ratio = [], []
    for k in range(0, cycle.n, 8):
        sel = cycle.lit[:, k] & (cycle.cos[:, k] > 0.5)
        base = renderer.cmap(np.clip((cycle.temperature[sel, k] - TMIN) / (TMAX - TMIN), 0, 1))[:, :3]
        cosine = renderer.vnormal[cycle.faces[sel]].mean(axis=1) @ cycle.sun[k]
        factor = AMBIENT + (1.0 - AMBIENT) * np.clip(cosine, 0, 1)
        shown = base * factor[:, None]
        nearest = ramp[np.argmin(((shown[:, None, :] - colors[None]) ** 2).sum(-1), axis=1)]
        out.append(nearest - cycle.temperature[sel, k])
        ratio.append(factor)
    out, ratio = np.concatenate(out), np.concatenate(ratio)
    return dict(n=int(out.size), median=float(np.median(out)), p95=float(np.percentile(np.abs(out), 95)),
                factor_median=float(np.median(ratio)))


def display_check(cycle, renderer, k):
    """How the displayed temperature at output k differs from the face values (the display only)."""
    info, (sx, sy, seen) = renderer.info, renderer.marker
    shift = np.abs(info['shift'])
    face = cycle.face
    t_face = float(cycle.temperature[face, k])
    t_site = float(np.mean(info['t_vertex'][cycle.faces[face]]))       # at the center of the face
    ix, iy = int(sx * renderer.ss), int(sy * renderer.ss)
    h, w = info['t_image'].shape
    t_marker = float(info['t_image'][iy, ix]) if seen and 0 <= ix < w and 0 <= iy < h else float('nan')
    return dict(k=k, lit=bool(cycle.lit[face, k]), seen=seen, max=float(shift.max()), p99=float(np.percentile(shift, 99)),
                median=float(np.median(shift)), face=t_face, center=t_site, pixel=t_marker, n=int(shift.size))


def report_display(cycle, diag):
    """Print the difference between the display and the face values of the site."""
    use = [d for d in diag if d['lit'] and d['seen'] and np.isfinite(d['pixel'])]
    err_c = np.array([abs(d['center'] - d['face']) for d in use])
    err_p = np.array([abs(d['pixel'] - d['face']) for d in use])
    print(f'marker versus plotted temperature of the site, {len(use)} lit outputs with the site in view: '
          f'face center {np.median(err_c):.2f} K median, {err_c.max():.2f} K worst; '
          f'pixel under the marker {np.median(err_p):.2f} K median, {err_p.max():.2f} K worst')
    worst = max(diag, key=lambda d: d['max'])
    print(f'displayed minus face temperature over all pixels and outputs: largest {worst["max"]:.1f} K '
          f'(output {worst["k"]}), median of the per-output 99th percentile {np.median([d["p99"] for d in diag]):.1f} K, '
          f'median {np.median([d["median"] for d in diag]):.2f} K')


# ------------------------------------------------------------------------------- figure

def build_figure(cycle, profiles, dpi):
    plt = plots.styled()
    import matplotlib as mpl
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    from matplotlib.patches import Patch
    mpl.rcParams.update({'font.size': 12, 'axes.labelsize': 12.5, 'xtick.labelsize': 11.5,
                         'ytick.labelsize': 11.5, 'legend.fontsize': 10.5, 'figure.dpi': dpi,
                         'savefig.bbox': 'standard'})
    W_in, H_in = FIGURE_IN
    n, c = cycle.n, cycle
    extents = view_extents(cycle)
    x0, x1, y0, y1 = extents
    pad = (x1 - x0) * (1.0 / FILL - 1.0)
    extents = (x0 - pad / 2, x1 + pad / 2, y0 - pad / 2, y1 + pad / 2 + SUN_BAND_KM)
    aspect = (extents[3] - extents[2]) / (extents[1] - extents[0])
    a_h = PANEL_A_IN * aspect                              # inches
    fig = plt.figure(figsize=(W_in, H_in), dpi=dpi, facecolor='white')

    def axes(left, bottom, width, height):                # inches to figure fractions
        return fig.add_axes([left / W_in, bottom / H_in, width / W_in, height / H_in])

    top = H_in - 0.30                                      # room for the panel label
    ax_a = axes(0.10, top - a_h, PANEL_A_IN, a_h)
    cax = axes(0.90, top - a_h - 0.24, 3.0, 0.11)
    b_top = top - a_h - 0.95
    ax_b = axes(0.98, 0.66, 4.65, b_top - 0.66)
    ax_c = axes(6.95, 0.66, 2.20, H_in - 0.66 - 0.60)

    # (a) the body
    ax_a.set_axis_off()
    renderer = Renderer(c, int(round(PANEL_A_IN * dpi)), int(round(a_h * dpi)), extents)
    image = ax_a.imshow(np.ones((renderer.h, renderer.w, 3)), extent=(0, renderer.w, renderer.h, 0),
                        interpolation='none', zorder=1)
    ax_a.set_xlim(0, renderer.w)
    ax_a.set_ylim(renderer.h, 0)
    marker = ax_a.plot([], [], 'o', ms=9, mfc=SITE_COLOR, mec='white', mew=1.7, zorder=5)[0]
    sun_arrow = ax_a.annotate('', xy=(0, 0), xytext=(0, 0), zorder=6,
                              arrowprops=dict(arrowstyle='-|>', color='0.25', lw=1.6, mutation_scale=13))
    sun_text = ax_a.text(0, 0, 'Sun', color='0.25', fontsize=11.5, ha='left', va='center', zorder=6)
    ax_a.text(0.0, 1.0, '(a)  433 Eros', transform=ax_a.transAxes, ha='left', va='bottom',
              fontsize=12.5, clip_on=False)
    cb = fig.colorbar(ScalarMappable(Normalize(TMIN, TMAX), cmap='inferno'), cax=cax,
                      orientation='horizontal')
    cb.set_ticks(np.arange(50, 271, 50))
    fig.text((0.90 + 3.0 + 0.15) / W_in, (top - a_h - 0.24 + 0.055) / H_in, 'Surface temperature (K)',
             ha='left', va='center', fontsize=12.5)

    # (b) surface temperature of the site
    t = c.t_h
    ts = profiles[0, :n]
    face = c.face
    night, shadow = c.night[face, :n], c.shadow[face, :n]
    step = c.period_h / n
    for mask, style in ((night, dict(facecolor='0.90', edgecolor='none')),
                        (shadow, dict(facecolor='none', edgecolor='0.45', hatch='////', linewidth=0.0))):
        edges = np.flatnonzero(np.diff(np.r_[0, mask.astype(int), 0]))
        for a, b in zip(edges[::2], edges[1::2]):
            ax_b.axvspan(a * step - 0.5 * step, b * step - 0.5 * step, zorder=0.1, **style)
    ax_b.plot(t, ts, color='0.72', lw=1.6, ls='-', zorder=2)
    elapsed, = ax_b.plot([], [], color='black', lw=2.2, ls='-', zorder=3)
    now_b, = ax_b.plot([], [], 'o', ms=8, mfc=SITE_COLOR, mec='white', mew=1.4, zorder=4, clip_on=False)
    ax_b.set_xlim(-0.5 * step, c.period_h - 0.5 * step)
    lo, hi = ts.min(), ts.max()
    ax_b.set_ylim(lo - 0.06 * (hi - lo), hi + 0.40 * (hi - lo))
    ax_b.set_xlabel('Time since local midnight (h)')
    ax_b.set_ylabel('Surface\ntemperature (K)')
    ax_b.set_xticks(range(0, int(c.period_h) + 1))
    ax_b.legend(handles=[Patch(facecolor='0.90', edgecolor='0.8', linewidth=0.8, label='Night'),
                         Patch(facecolor='none', edgecolor='0.45', hatch='////', linewidth=0.0,
                               label='Shadow of terrain')],
                loc='upper center', ncol=2, frameon=True, framealpha=1.0, edgecolor='none', handlelength=1.4,
                columnspacing=1.4, borderaxespad=0.15)
    lat, lon = c.lat[face], c.lon[face] % 360.0
    ax_b.text(0.0, 1.03, f'(b)  Site at {abs(lat):.0f}°{"N" if lat >= 0 else "S"}, {lon:.0f}°E',
              transform=ax_b.transAxes, ha='left', va='bottom', fontsize=12.5)

    # (c) the profile in depth
    z_cm = c.depth_m * 100.0
    delta_cm = c.skin_depth_m * 100.0
    low, high = profiles[:, :n].min(axis=1), profiles[:, :n].max(axis=1)
    ax_c.fill_betweenx(z_cm, low, high, color='0.85', lw=0, zorder=1)
    ax_c.plot(low, z_cm, color='0.72', lw=1.0, ls='-', zorder=2)
    ax_c.plot(high, z_cm, color='0.72', lw=1.0, ls='-', zorder=2)
    profile, = ax_c.plot([], [], marker='o', color='black', lw=2.0, ls='-', ms=5.2, mfc=SITE_COLOR,
                         mec='black', mew=0.8, zorder=4)
    ax_c.set_ylim(z_cm[-1] + 0.35, -0.2)
    span = high.max() - low.min()
    ax_c.set_xlim(low.min() - 0.06 * span, high.max() + 0.24 * span)
    ax_c.set_xticks(np.arange(np.ceil(low.min() / 20) * 20, high.max() + 1, 20))
    ax_c.set_xlabel('Temperature (K)')
    ax_c.set_ylabel('Depth (cm)')
    ax_c.legend(handles=[Patch(facecolor='0.85', edgecolor='none', label='Diurnal range')],
                loc='lower right', frameon=False, handlelength=1.4, borderaxespad=0.2)
    from matplotlib.ticker import AutoMinorLocator, FixedLocator
    # centimeters on the left only, skin depths on the right only, both from the same zero
    ax_c.set_yticks(np.arange(0, int(z_cm[-1]) + 1))
    ax_c.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax_c.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False)
    sec = ax_c.secondary_yaxis('right', functions=(lambda z: z / delta_cm, lambda d: d * delta_cm))
    sec.set_ylabel('Depth (skin depths)')
    sec.yaxis.set_major_locator(FixedLocator(np.arange(0, int(z_cm[-1] / delta_cm) + 1)))
    sec.yaxis.set_minor_locator(AutoMinorLocator(2))
    sec.tick_params(axis='y', which='both', left=False, right=True, labelleft=False, labelright=True)
    ax_c.text(0.0, 1.02, '(c)  Depth profile at the site', transform=ax_c.transAxes,
              ha='left', va='bottom', fontsize=12.5)

    def draw(k):
        rgb, info = renderer.render(k)
        image.set_data(rgb)
        sx, sy, seen = renderer.site(k, info)
        marker.set_data([sx] if seen else [], [sy] if seen else [])
        renderer.marker = (sx, sy, seen)
        s = info['sun_view']
        sdx, sdy = float(s @ renderer.right), -float(s @ renderer.up)
        norm = math.hypot(sdx, sdy) or 1.0
        x0, y0 = 0.04 * renderer.w, 0.20 * SUN_BAND_KM * renderer.scale / renderer.ss  # upper left
        length = 0.05 * renderer.w
        tip = (x0 + length * sdx / norm, y0 + length * sdy / norm)
        sun_arrow.set_position((x0, y0))
        sun_arrow.xy = tip
        sun_text.set_position((tip[0] + 0.015 * renderer.w, tip[1]))
        elapsed.set_data(t[:k + 1], ts[:k + 1])
        now_b.set_data([t[k]], [ts[k]])
        profile.set_data(profiles[:, k], z_cm)
        renderer.info = info
        return rgb

    return fig, draw, renderer


# ------------------------------------------------------------------------------ encoding

def encode_gif(frames, path, delay_ms, ffmpeg):
    """A looping GIF with one palette for all frames, dithered so the color ramp does not band."""
    path = Path(path)
    if ffmpeg:
        pal = path.with_suffix('.palette.png')
        pattern = str(Path(frames[0]).parent / 'frame_%03d.png')
        fps = 1000.0 / delay_ms
        subprocess.run([ffmpeg, '-y', '-v', 'error', '-framerate', f'{fps:g}', '-i', pattern,
                        '-vf', 'palettegen=max_colors=256:stats_mode=full:reserve_transparent=0', str(pal)],
                       check=True)
        subprocess.run([ffmpeg, '-y', '-v', 'error', '-framerate', f'{fps:g}', '-i', pattern, '-i', str(pal),
                        '-lavfi', 'paletteuse=dither=bayer:bayer_scale=4:diff_mode=rectangle', '-loop', '0',
                        '-map_metadata', '-1', str(path)], check=True)
        pal.unlink()
        return 'ffmpeg palettegen/paletteuse'
    from PIL import Image
    images = [Image.open(f).convert('RGB') for f in frames]
    sample = Image.new('RGB', (images[0].width, images[0].height * 8))
    for i, k in enumerate(np.linspace(0, len(images) - 1, 8).astype(int)):
        sample.paste(images[k], (0, i * images[0].height))
    palette = sample.quantize(colors=256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
    quantized = [im.quantize(palette=palette, dither=Image.Dither.FLOYDSTEINBERG) for im in images]
    quantized[0].save(path, save_all=True, append_images=quantized[1:], duration=delay_ms, loop=0,
                      optimize=False, disposal=1)
    return 'Pillow, one shared palette'


def encode_mp4(frames, path, fps, ffmpeg):
    if not ffmpeg:
        print('ffmpeg not found: no MP4 written')
        return
    pattern = str(Path(frames[0]).parent / 'frame_%03d.png')
    subprocess.run([ffmpeg, '-y', '-v', 'error', '-framerate', f'{fps:g}', '-i', pattern, '-vf',
                    'pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white', '-c:v', 'libx264', '-pix_fmt', 'yuv420p',
                    '-crf', '14', '-map_metadata', '-1', str(path)], check=True)


# ---------------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_dir', help='output directory of examples/asteroid_eros_animation.yaml')
    ap.add_argument('--out', default='eros_diurnal_cycle.gif', help='the GIF to write')
    ap.add_argument('--mp4', help='also write an MP4 at --mp4-dpi')
    ap.add_argument('--png', help='also write one frame, chosen by --png-step')
    ap.add_argument('--png-step', type=int, default=None, help='output index of the PNG (default: 55 percent into the rotation)')
    ap.add_argument('--face', type=int, help='face of the site (default: the automatic choice)')
    ap.add_argument('--list-sites', action='store_true', help='print the ranking of the automatic choice and stop')
    ap.add_argument('--dpi', type=float, default=130.0, help='resolution of the GIF frames')
    ap.add_argument('--mp4-dpi', type=float, default=200.0)
    ap.add_argument('--delay-ms', type=int, default=60, help='time between GIF frames')
    ap.add_argument('--frames-dir', help='where to keep the frames (default: a temporary directory)')
    args = ap.parse_args()

    cycle = Cycle(args.run_dir, args.face)
    if args.list_sites:
        for score, count, i in cycle.rank_sites()[:12]:
            print(f'face {i:5d}  lat {cycle.lat[i]:6.1f}  lon {cycle.lon[i] % 360:6.1f}  '
                  f'{count} of {cycle.n} outputs in terrain shadow (neighbors {score})')
        return
    profiles = cycle.site_profiles()
    n = cycle.n
    closure = np.abs(profiles[:, n].astype(float) - profiles[:, 0].astype(float))
    print(f'site: face {cycle.face}, latitude {cycle.lat[cycle.face]:.1f}, '
          f'east longitude {cycle.lon[cycle.face] % 360:.1f}')
    print(f'{cycle.n_faces} faces, {int(cycle.engine_flipped.sum())} of them with a normal that the engine turns around')
    print(f'sun {cycle.sun_distance_au:.3f} AU, period {cycle.period_h:.4f} h, {n} outputs per rotation')
    print('stored depths (cm):', np.round(cycle.depth_m * 100, 2).tolist(), f'skin depth {cycle.skin_depth_m*100:.3f} cm')
    print('loop closure |T(P) - T(0)| by level (K):', np.round(closure, 3).tolist(), f'max {closure.max():.3f}')
    ffmpeg = shutil.which('ffmpeg')
    keep = Path(args.frames_dir) if args.frames_dir else None
    with tempfile.TemporaryDirectory() as scratch:
        frames_dir = keep or Path(scratch)
        frames_dir.mkdir(parents=True, exist_ok=True)
        fig, draw, renderer = build_figure(cycle, profiles, args.dpi)
        fid = color_fidelity(cycle, renderer)
        print(f'color fidelity of lit faces (cos > 0.5), shaded color mapped to temperature: '
              f'median {fid["median"]:+.1f} K, 95th percentile of |error| {fid["p95"]:.1f} K '
              f'(median light factor {fid["factor_median"]:.2f}, ambient {AMBIENT}, smooth normals)')
        frames = []
        diag = []
        for k in range(n):
            draw(k)
            diag.append(display_check(cycle, renderer, k))
            frame = frames_dir / f'frame_{k:03d}.png'
            fig.savefig(frame, dpi=args.dpi, facecolor='white')
            frames.append(frame)
            if args.png and k == (args.png_step if args.png_step is not None else int(0.55 * n)):
                fig.savefig(args.png, dpi=args.dpi * 1.5, facecolor='white')
        report_display(cycle, diag)
        how = encode_gif(frames, args.out, args.delay_ms, ffmpeg)
        size = Path(args.out).stat().st_size / 1e6
        print(f'wrote {args.out}: {n} frames, {args.delay_ms} ms each, {size:.2f} MB ({how})')
        if args.mp4:
            hi_fig, hi_draw, _ = build_figure(cycle, profiles, args.mp4_dpi)
            hi_dir = frames_dir / 'hi'
            hi_dir.mkdir(exist_ok=True)
            hi_frames = []
            for k in range(n):
                hi_draw(k)
                hi_fig.savefig(hi_dir / f'frame_{k:03d}.png', dpi=args.mp4_dpi, facecolor='white')
                hi_frames.append(hi_dir / f'frame_{k:03d}.png')
            encode_mp4(hi_frames, args.mp4, 1000.0 / args.delay_ms, ffmpeg)
            print(f'wrote {args.mp4}')


if __name__ == '__main__':
    sys.exit(main())
