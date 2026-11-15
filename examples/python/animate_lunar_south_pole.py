#!/usr/bin/env python3
"""Animate the lunar south polar region over one solar day, with the Earth in the sky.

The animation of the south polar region on the front page of the repository comes from this
script. The region is the cap poleward of 87 S (radius 91 km), seen obliquely from a camera fixed
in the Moon-fixed frame, looking across the pole toward the Earth. The surface is colored by the
temperature of the run, with radiation exchange between faces (scattered sunlight and thermal
emission) included, and lit from the direction of the Sun at every output. The Earth is drawn at its
true position, angular size and phase. Everything shown comes from the run or from SPICE.

Data. The LOLA polar elevation models (data set LRO-L-LOLA-4-GDR-V1.0, NASA Planetary Data System;
Smith et al., 2010, doi:10.1029/2010GL043751) of the region poleward of 85 S at 40 m/pixel (for the
simulated region) and of the region poleward of 80 S at 80 m/pixel (for the terrain from 91 to
300 km that casts shadows into it), 115 MB each::

    B=https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/polar/img
    mkdir -p data/dem/lola data/dem/south_pole
    for f in LDEM_85S_40M LDEM_80S_80M; do
        curl -L -o data/dem/lola/$f.IMG $B/$f.IMG
        curl -L -o data/dem/lola/$f.LBL $B/$f.LBL
    done
    shasum -a 256 data/dem/lola/*.IMG
    #  1609e9b860748e0b9952a7239c51b37d3b2a799d35e9fc3c258526711cc1e51d  LDEM_85S_40M.IMG
    #  e6591eeea35402ca5e124ca00f5bf0774381a2d4fa919144b455d82b8ff65fdc  LDEM_80S_80M.IMG

The meshes (``dtm_to_ply_delaunay`` needs GDAL, see tutorial 5, and the ring builder the GDAL
Python bindings and scipy). The simulated region is sampled every 18 pixels (720 m), which gives
99,887 faces, just below the 100,000 faces above which the model computes the view factors in
compressed (H-matrix) form::

    ./build/bin/dtm_to_ply_delaunay data/dem/lola/LDEM_85S_40M.LBL data/dem/south_pole/south_pole_cap.ply \\
        --body-radius 1737400 --downsample 18 --roi-center -90 0 --roi-radius 91
    python3 scripts/build_lola_ring.py --lola data/dem/lola/LDEM_80S_80M.LBL --site-lat -90 --site-lon 0 \\
        --r-core-km 91 --r-ext-km 300 --max-points 150000 \\
        --core-ply data/dem/south_pole/south_pole_cap.ply --out data/dem/south_pole/south_pole_ring.ply

The run and the animation, from the repository root, after building ``ptm_thermal_sim`` and
fetching the kernels (``python3 scripts/fetch_kernels.py --bodies moon``)::

    OMP_NUM_THREADS=8 ./build/bin/ptm_thermal_sim --config examples/lunar_south_pole_animation.yaml \\
        --thermal-threads 8
    PYTHONPATH=src/python python3 examples/python/animate_lunar_south_pole.py \\
        output/lunar_south_pole_animation --ring data/dem/south_pole/south_pole_ring.ply \\
        --out docs/figures/lunar_south_pole_solar_day.gif

Measured on a 2019 Mac Pro with 8 threads, CPU only: the meshes take about 15 s (the ring builder
needs about 5 GB of memory), the run 18 minutes with a peak memory of 9.1 GB (8 lunations of
spin-up, then 301 outputs, three iterations of the radiation exchange), and its view factors,
35 million nonzero entries, take 31 s and are saved in output/lunar_south_pole.vf (425 MB). The
run directory is 0.8 GB. The animation takes about 2 minutes and 2 GB, or 11 minutes and 6 GB
with ``--mp4``. The script needs numpy, matplotlib, Pillow and spiceypy.

The window. The Sun circles the horizon of the pole once per solar day, which is 29.5147 days from
the new Moon of 2029-09-08T10:31:36 UTC (the time for the sub-solar longitude to return to its
starting value). A window from one new Moon to the next, 29.37 days here, does not close: the
interval between new Moons depends on the libration of the Moon, the solar day does not. The run
covers two solar days in 301 outputs. The first lets the model settle at the output cadence after
the spin-up, whose forcing is sampled 24 times per lunation, and is not drawn. The second is drawn.
The surface temperature at its end differs from that at its start by +0.99 K (median over the
faces), by less than 5.4 K on 99% of the faces, and by up to 80 K on six faces where a sliver of
the Sun clears the horizon at the start of the day but not at its end. The lunar month is not
exactly periodic: the sub-solar latitude moves between -1.59 and -1.45 degrees, and the floors
that the Sun never reaches warm by 0.8 K over the drawn day (median; 1.3 K over the day before it),
which the closure includes.

The picture. The terrain is rasterized once, with its own z-buffer and perspective, since the
camera does not move. Heights are drawn at true scale. Temperature is interpolated for
display only: each vertex takes the area-weighted mean of the temperatures of the faces around it,
interpolated across each face. The color of a pixel is its temperature on the color bar (a power
law in temperature, so that the cold shadowed floors and the sunlit slopes are both resolved) times
``ambient + (1 - ambient) * min(1, gain * diffuse)``. The diffuse term is the cosine between the
smoothed surface normal and the direction of the Sun, multiplied by the share of the solar disk that
the model sees from each face (its direct flux over the flux it would receive unshadowed), so the
shadows are those of the model. Terrain outside the simulated region, and the rest of the Moon out
to its limb, are drawn in gray with a fixed shading that does not follow the Sun, since they are
not simulated. The Earth and the Sun
are drawn at their true angular sizes where they are above the terrain, from SPICE positions in
the frame of the model (IAU_MOON) moved from the center of the Moon to the camera. The Earth is lit
from the direction of the Sun as seen from the Earth, which gives its phase. The script prints its
position and phase against independent SPICE numbers (sub-Earth point, elevation from the pole,
phase angle, illuminated fraction).

``--mp4`` and ``--png`` also write a video and one frame. The GIF is encoded with ``ffmpeg`` when
it is installed and with Pillow otherwise. The video needs ``ffmpeg``. The GIF draws every second output (75 frames, 9.4 hours
apart) at 1000 x 500 pixels, which keeps it under 8 MB with a full palette of 256 colors, and the
video draws every output.
"""
import argparse
import math
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from animate_eros import encode_gif, encode_mp4, rasterize, vertex_normals  # noqa: E402
from ptm import Results, plots  # noqa: E402
from ptm.results import read_ply_mesh_minimal  # noqa: E402

R_MOON = 1737400.0                # m, the reference sphere of the LOLA heights and of the mesh
R_EARTH = 6371.0                  # km, mean radius, for the angular size of the Earth
R_SUN = 695700.0                  # km
SOLAR_CONSTANT = 1361.0           # W/m2 at 1 AU, as in the engine
AU_KM = 1.495978707e8

# The camera, in the frame of the pole. An azimuth is the east longitude of a horizontal direction
# at the pole, an elevation is measured from the horizontal plane of the pole (its zenith is -z).
VIEW_AZIMUTH = 0.0                # the camera looks toward longitude 0, the side of the Earth
CAMERA_DISTANCE_KM = 190.0        # horizontal distance of the camera from the pole, on the far side
CAMERA_HEIGHT_KM = 40.0           # height of the camera above the horizontal plane of the pole
CAMERA_PITCH = 6.5                # degrees below the horizontal of the pole at the center of the image
FIELD_OF_VIEW = 55.0              # degrees, horizontal

AMBIENT = 0.65                    # share of the color that no lighting removes
GAIN = 3.0                        # the Sun term reaches full color at an incidence cosine of 1/GAIN
TMIN, TMAX, GAMMA = 20.0, 330.0, 0.6   # color scale: ((T - TMIN) / (TMAX - TMIN)) ** GAMMA
CBAR_TICKS = (20, 40, 60, 100, 150, 200, 300)
SURROUND = 0.13                   # gray of the terrain outside the simulated region
EARTH_DAY = np.array([0.80, 0.87, 0.97])
SUN_COLOR = np.array([1.0, 1.0, 1.0])
SUPERSAMPLE = 3
FIGURE_IN = (10.0, 5.0)
STYLE = dict(ambient=AMBIENT, gain=GAIN, tmin=TMIN, tmax=TMAX, gamma=GAMMA)


# ---------------------------------------------------------------------------- geometry

def pole_frame():
    """Unit vectors toward longitude 0, toward longitude 90 E, and the zenith, at the south pole."""
    return np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, -1.0])


def direction(azimuth, elevation):
    x, y, z = pole_frame()
    a, e = math.radians(azimuth), math.radians(elevation)
    return math.cos(e) * (math.cos(a) * x + math.sin(a) * y) + math.sin(e) * z


def azimuth_elevation(v):
    """Azimuth (east longitude) and elevation in degrees of a direction, in the frame of the pole."""
    x, y, z = pole_frame()
    v = np.asarray(v, float)
    v = v / np.linalg.norm(v, axis=-1, keepdims=True)
    return np.degrees(np.arctan2(v @ y, v @ x)), np.degrees(np.arcsin(np.clip(v @ z, -1, 1)))


POLE = np.array([0.0, 0.0, -R_MOON])


class Camera:
    """A pinhole camera fixed in the Moon-fixed frame."""

    def __init__(self, width, height, azimuth=VIEW_AZIMUTH, distance_km=CAMERA_DISTANCE_KM,
                 height_km=CAMERA_HEIGHT_KM, pitch=CAMERA_PITCH, fov=FIELD_OF_VIEW):
        self.w, self.h = width, height
        zenith = pole_frame()[2]
        self.position = POLE - 1000.0 * distance_km * direction(azimuth, 0.0) + 1000.0 * height_km * zenith
        self.forward = direction(azimuth, -pitch)
        right = np.cross(self.forward, zenith)
        self.right = right / np.linalg.norm(right)
        self.up = np.cross(self.right, self.forward)
        self.focal = 0.5 * width / math.tan(math.radians(0.5 * fov))   # pixels
        self.fov, self.vfov = fov, 2.0 * math.degrees(math.atan(0.5 * height / self.focal))

    def project(self, points):
        """Pixel coordinates, inverse depth (grows toward the camera) and depth (m) of points."""
        d = points - self.position
        z = d @ self.forward
        zs = np.where(z > 1.0, z, 1.0)
        px = np.stack([0.5 * self.w + self.focal * (d @ self.right) / zs,
                       0.5 * self.h - self.focal * (d @ self.up) / zs], axis=1)
        return px, 1.0 / zs, z

    def rays(self):
        """Unit direction of every pixel center, (H, W, 3)."""
        gx, gy = np.meshgrid(np.arange(self.w) + 0.5, np.arange(self.h) + 0.5)
        d = ((gx - 0.5 * self.w)[..., None] * self.right - (gy - 0.5 * self.h)[..., None] * self.up
             + self.focal * self.forward)
        return d / np.linalg.norm(d, axis=-1, keepdims=True)


def outward(vertices, faces):
    """The faces wound so that their normals point away from the center of the Moon."""
    a, b, c = (vertices[faces[:, i]] for i in range(3))
    flip = np.sum(np.cross(b - a, c - a) * (a + b + c), axis=1) < 0
    faces = faces.copy()
    faces[flip] = faces[flip][:, [0, 2, 1]]
    return faces


def face_geometry(vertices, faces):
    a, b, c = (vertices[faces[:, i]] for i in range(3))
    cross = np.cross(b - a, c - a)
    area = 0.5 * np.linalg.norm(cross, axis=1)
    return area, cross / (2.0 * area[:, None]), (a + b + c) / 3.0


def read_ply(path):
    m = read_ply_mesh_minimal(path)
    return np.asarray(m['vertices'], float), np.asarray(m['faces'], np.int64)


# ------------------------------------------------------------------------------- the run

class SolarDay:
    """The arrays of the last solar day of the run: ``frames`` outputs and the output that closes the loop.

    The run covers two solar days. The first lets the model settle at the output cadence after
    the spin-up, whose forcing is sampled more coarsely, and is not drawn.
    """

    def __init__(self, run_dir, frames):
        self.run = run = Results(run_dir, mesh_units='m')
        self.vertices = np.asarray(run.mesh['vertices'], float)
        self.faces = outward(self.vertices, np.asarray(run.mesh['faces'], np.int64))
        self.area, self.normal, self.center = face_geometry(self.vertices, self.faces)
        times = run.times()
        total = len(times.et)
        if total < frames + 1:
            raise SystemExit(f'the run has {total} outputs, fewer than {frames} frames and the closing output')
        steps = list(range(total - frames - 1, total))
        self.first_output = steps[0]
        self.et = np.asarray(times.et, float)[steps]
        self.utc = [str(u) for u in np.asarray(times.utc)[steps]]
        self.n = frames                                    # output n of the window closes the loop
        pos = run.body_positions()
        self.sun_km = np.array([[pos[i]['sun_x'], pos[i]['sun_y'], pos[i]['sun_z']] for i in steps], float)
        self.earth_km = np.array([[pos[i]['primary_x'], pos[i]['primary_y'], pos[i]['primary_z']] for i in steps], float)
        self.sun = self.sun_km / np.linalg.norm(self.sun_km, axis=1)[:, None]
        self.temperature = np.array(run.stack('temperature', steps=steps)).T          # (faces, outputs) K
        self.flux = np.array(run.stack('direct_solar', steps=steps)).T                # (faces, outputs) W/m2
        dist_au = np.linalg.norm(self.sun_km, axis=1) / AU_KM
        unshadowed = SOLAR_CONSTANT / dist_au[None, :] ** 2 * np.clip(self.normal @ self.sun.T, 0.0, None)
        # The share of the solar disk that each face sees in the model: 0 in shadow, 1 fully lit.
        self.lit = np.where((unshadowed > 1e-6) & (self.flux > 0.0),
                            np.clip(self.flux / np.maximum(unshadowed, 1e-12), 0.0, 1.0), 0.0)

    def day(self, k):
        return (self.et[k] - self.et[0]) / 86400.0


# ----------------------------------------------------------------------- the Earth (SPICE)

def sky_geometry(et, observer_m, kernels):
    """The Earth and the Sun seen from a point, from SPICE.

    Positions in IAU_MOON (the frame of the model) with light time and stellar aberration for an
    observer at the center of the Moon, then moved to the observer: the parallax of a point on the
    surface against the center of the Moon is up to 0.26 degrees for the Earth.
    """
    import spiceypy as sp
    sp.furnsh(str(kernels))
    obs_km = np.asarray(observer_m, float) / 1000.0
    out = []
    for t in et:
        e = np.asarray(sp.spkpos('EARTH', t, 'IAU_MOON', 'LT+S', 'MOON')[0]) - obs_km
        s = np.asarray(sp.spkpos('SUN', t, 'IAU_MOON', 'LT+S', 'MOON')[0]) - obs_km
        d, ds = np.linalg.norm(e), np.linalg.norm(s)
        out.append(dict(dir=e / d, radius=math.asin(R_EARTH / d), dist_km=d,
                        sun_at_earth=(s - e) / np.linalg.norm(s - e),
                        sun_dir=s / ds, sun_radius=math.asin(R_SUN / ds)))
    sp.kclear()
    return out


def spice_reference(et, kernels):
    """Numbers from separate SPICE routines, to check the picture against.

    The sub-Earth point (subpnt), the elevation of the Earth seen from the pole (azlcpo, with its
    own light-time and aberration corrections), the phase angle at the Earth (phaseq) and the
    illuminated fraction of its disk, (1 + cos(phase)) / 2.
    """
    import spiceypy as sp
    sp.furnsh(str(kernels))
    rows = []
    for t in et:
        spoint = sp.subpnt('INTERCEPT/ELLIPSOID', 'MOON', t, 'IAU_MOON', 'LT+S', 'EARTH')[0]
        _, lon, lat = sp.reclat(spoint)
        azel = sp.azlcpo('ELLIPSOID', 'EARTH', t, 'LT+S', False, True, POLE / 1000.0, 'MOON', 'IAU_MOON')[0]
        phase = sp.phaseq(t, 'EARTH', 'SUN', 'MOON', 'LT+S')
        rows.append((math.degrees(lat), math.degrees(lon), math.degrees(azel[2]), math.degrees(phase),
                     0.5 * (1.0 + math.cos(phase))))
    sp.kclear()
    return np.array(rows)


# ------------------------------------------------------------------------------ picture

def norm_value(t):
    x = np.clip((np.asarray(t, float) - STYLE['tmin']) / (STYLE['tmax'] - STYLE['tmin']), 0.0, 1.0)
    return x ** STYLE['gamma']


def colors_of(t):
    from matplotlib import colormaps
    return colormaps[plots.CMAP_TEMPERATURE](norm_value(t))[..., :3]


class Scene:
    """The terrain rasterized once (the camera is fixed), and shaded anew at every output."""

    def __init__(self, sd, ring, camera):
        self.sd, self.cam = sd, camera
        cv = sd.vertices
        self.vnormal = vertex_normals(cv, sd.faces, sd.normal, sd.area)
        self.vertex_area = np.zeros(len(cv))
        for i in range(3):
            np.add.at(self.vertex_area, sd.faces[:, i], sd.area)
        rv, rf = ring
        rf = outward(rv, rf)
        self.n_core = len(sd.faces)
        verts = np.vstack([cv, rv])
        faces = np.vstack([sd.faces, rf + len(cv)])
        px, inv_z, z = camera.project(verts)
        W, H = camera.w, camera.h
        x, y = px[faces, 0], px[faces, 1]
        keep = np.flatnonzero((z[faces] > 1000.0).all(axis=1) & (x.max(1) >= 0) & (x.min(1) < W)
                              & (y.max(1) >= 0) & (y.min(1) < H))
        _, fid = rasterize(px, inv_z, faces[keep], W, H)
        fid = np.where(fid >= 0, keep[np.maximum(fid, 0)], -1)
        self.core = (fid >= 0) & (fid < self.n_core)
        self.surround = fid >= self.n_core
        # perspective-correct barycentric weights at the pixels of the simulated region
        tri = faces[fid[self.core]]
        gy, gx = np.nonzero(self.core)
        gx, gy = gx + 0.5, gy + 0.5
        xs, ys = px[tri, 0], px[tri, 1]
        det = (ys[:, 1] - ys[:, 2]) * (xs[:, 0] - xs[:, 2]) + (xs[:, 2] - xs[:, 1]) * (ys[:, 0] - ys[:, 2])
        det = np.where(det == 0.0, 1.0, det)
        l1 = ((ys[:, 1] - ys[:, 2]) * (gx - xs[:, 2]) + (xs[:, 2] - xs[:, 1]) * (gy - ys[:, 2])) / det
        l2 = ((ys[:, 2] - ys[:, 0]) * (gx - xs[:, 2]) + (xs[:, 0] - xs[:, 2]) * (gy - ys[:, 2])) / det
        w = np.clip(np.stack([l1, l2, 1.0 - l1 - l2], axis=1), 0.0, None) * inv_z[tri]
        self.weights = w / w.sum(axis=1, keepdims=True)
        self.tri = tri
        self.core_face = fid[self.core]
        # the rest of the Moon out to its limb, where neither mesh reaches
        self.rays = camera.rays()
        o = camera.position
        b = self.rays @ o
        disc = b * b - (o @ o - R_MOON ** 2)
        hit = (disc > 0) & (-b - np.sqrt(np.maximum(disc, 0)) > 0)
        self.sphere = hit & (fid < 0)
        self.sky = ~hit & (fid < 0)
        # The surroundings are not simulated: a neutral gray with a fixed relief shading, not the Sun.
        _, rn, _ = face_geometry(verts, faces)
        sf = fid[self.surround]
        self.surround_gray = SURROUND * (0.7 + 0.3 * np.clip(-(rn[sf] @ camera.forward), 0.0, 1.0))

    def vertex_mean(self, values):
        f = self.sd.faces
        acc = np.zeros(len(self.vertex_area))
        for i in range(3):
            np.add.at(acc, f[:, i], self.sd.area * values)
        return acc / self.vertex_area

    def render(self, k, sky):
        sd, cam = self.sd, self.cam
        rgb = np.zeros((cam.h, cam.w, 3))
        t_vertex = self.vertex_mean(sd.temperature[:, k].astype(float))
        diffuse = np.clip(self.vnormal @ sd.sun[k], 0.0, 1.0) * self.vertex_mean(sd.lit[:, k])
        w = self.weights
        t_pixel = np.sum(w * t_vertex[self.tri], axis=1)
        d_pixel = np.clip(STYLE['gain'] * np.sum(w * diffuse[self.tri], axis=1), 0.0, 1.0)
        rgb[self.core] = colors_of(t_pixel) * (STYLE['ambient'] + (1.0 - STYLE['ambient']) * d_pixel)[:, None]
        rgb[self.surround] = self.surround_gray[:, None]
        rgb[self.sphere] = SURROUND * 0.7
        g = sky[k]
        info = dict(shift=t_pixel - sd.temperature[self.core_face, k])
        rays = self.rays[self.sky]
        sub = rgb[self.sky]
        # the Sun at its true size, where it is above the terrain
        sun = np.arccos(np.clip(rays @ g['sun_dir'], -1, 1)) < g['sun_radius']
        sub[sun] = SUN_COLOR
        info['sun_px'] = int(sun.sum())
        # the Earth at its true size and phase, where it is above the terrain
        on = np.arccos(np.clip(rays @ g['dir'], -1, 1)) < g['radius']
        info['earth_disk_px'] = int((np.arccos(np.clip(self.rays @ g['dir'], -1, 1)) < g['radius']).sum())
        info['earth_px'] = int(on.sum())
        info['earth_lit_share'] = info['earth_lit_centroid'] = None
        if on.any():
            r = rays[on]
            perp = r - np.outer(r @ g['dir'], g['dir'])
            s = np.clip(np.linalg.norm(perp, axis=1) / math.sin(g['radius']), 0.0, 1.0)
            normal = perp / np.maximum(np.linalg.norm(perp, axis=1), 1e-15)[:, None] * s[:, None] \
                - np.sqrt(1.0 - s * s)[:, None] * g['dir']                  # the hemisphere facing us
            light = np.clip(normal @ g['sun_at_earth'], 0.0, 1.0)
            sub[on] = EARTH_DAY[None, :] * (light ** 0.6)[:, None]
            lit = light > 0
            info['earth_lit_share'] = float(lit.mean())
            info['earth_lit_centroid'] = r[lit].mean(axis=0) if lit.any() else None
        rgb[self.sky] = sub
        return rgb, info


def downsample(rgb, ss):
    H, W, _ = rgb.shape
    return rgb.reshape(H // ss, ss, W // ss, ss, 3).mean(axis=(1, 3))


# ------------------------------------------------------------------------------- figure

def build_figure(sd, scene, sky, dpi, ss):
    """The full-bleed picture with a color bar and the time drawn in the sky."""
    plt = plots.styled()
    import matplotlib as mpl
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import FuncNorm
    mpl.rcParams.update({'font.size': 12, 'figure.dpi': dpi, 'savefig.bbox': 'standard',
                         'text.color': 'white', 'axes.labelcolor': 'white', 'xtick.color': 'white',
                         'ytick.color': 'white', 'axes.edgecolor': 'white'})
    W_in, H_in = FIGURE_IN
    fig = plt.figure(figsize=(W_in, H_in), dpi=dpi, facecolor='black')
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    w, h = scene.cam.w // ss, scene.cam.h // ss
    image = ax.imshow(np.zeros((h, w, 3)), extent=(0, w, h, 0), interpolation='none')
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    cax = fig.add_axes([7.0 / W_in, 1.0 - 0.52 / H_in, 2.6 / W_in, 0.11 / H_in])
    tmin, tmax, gamma = STYLE['tmin'], STYLE['tmax'], STYLE['gamma']
    # The color bar uses the same mapping as the picture (norm_value), over the same limits.
    norm = FuncNorm((norm_value, lambda x: tmin + np.asarray(x, float) ** (1 / gamma) * (tmax - tmin)),
                    vmin=tmin, vmax=tmax)
    mappable = ScalarMappable(norm, cmap=plots.CMAP_TEMPERATURE)
    ticks = np.array(CBAR_TICKS, float)
    if not np.allclose(mappable.to_rgba(ticks)[:, :3], colors_of(ticks)):
        raise SystemExit('the color bar does not match the colors of the picture')
    cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')
    cb.set_ticks(CBAR_TICKS)
    cb.ax.minorticks_off()
    cb.ax.tick_params(labelsize=11, length=3, color='white')
    cb.outline.set_edgecolor('white')
    cb.outline.set_linewidth(0.6)
    fig.text(cax.get_position().x0, cax.get_position().y1 + 0.03, 'Surface temperature (K)',
             ha='left', va='bottom', fontsize=12)
    when = fig.text(0.35 / W_in, 1.0 - 0.30 / H_in, '', ha='left', va='top', fontsize=12.5, linespacing=1.5)
    span = sd.day(sd.n)

    def draw(k):
        rgb, info = scene.render(k, sky)
        image.set_data(np.clip(downsample(rgb, ss), 0.0, 1.0))
        utc = sd.utc[k]
        when.set_text(f'{utc[:10]}  {utc[11:16]} UTC\nDay {sd.day(k):4.1f} of {span:.1f}')
        return info

    return fig, draw


# --------------------------------------------------------------------------------- checks

def check_earth(sd, sky, scene, kernels, frames):
    """Print the Earth of the picture against independent SPICE numbers."""
    ref = spice_reference(sd.et[:sd.n], kernels)
    pole_sky = sky_geometry(sd.et[:sd.n], POLE, kernels)
    az, el = azimuth_elevation(np.array([g['dir'] for g in pole_sky]))
    print(f'Earth seen from the pole: elevation {el.min():+.2f} to {el.max():+.2f} deg, azimuth '
          f'{az.min():+.2f} to {az.max():+.2f} deg (east longitude), above the horizontal at '
          f'{int((el > 0).sum())} of {sd.n} frames; elevation minus SPICE azlcpo: max '
          f'{np.abs(el - ref[:, 2]).max():.4f} deg')
    print(f'  sub-Earth latitude {ref[:, 0].min():+.2f} to {ref[:, 0].max():+.2f} deg, longitude '
          f'{ref[:, 1].min():+.2f} to {ref[:, 1].max():+.2f} deg; angular diameter '
          f'{2 * math.degrees(min(g["radius"] for g in pole_sky)):.3f} to '
          f'{2 * math.degrees(max(g["radius"] for g in pole_sky)):.3f} deg')
    center = sky_geometry(sd.et[:sd.n], np.zeros(3), kernels)
    eng_earth = sd.earth_km[:sd.n] / np.linalg.norm(sd.earth_km[:sd.n], axis=1)[:, None]

    def angle(a, b):
        return np.degrees(np.arccos(np.clip(np.sum(a * b, axis=1), -1, 1)))
    print(f'  directions written by the model vs SPICE here: Sun max {angle(sd.sun[:sd.n], np.array([g["sun_dir"] for g in center])).max():.5f} deg, '
          f'Earth max {angle(eng_earth, np.array([g["dir"] for g in center])).max():.5f} deg')
    in_view = 0
    for k in frames:
        _, info = scene.render(k, sky)
        if not info['earth_px'] or info['earth_px'] < info['earth_disk_px']:
            continue
        in_view += 1
        g = sky[k]
        toward = g['sun_dir'] - (g['sun_dir'] @ g['dir']) * g['dir']      # the bright limb faces the Sun
        toward /= np.linalg.norm(toward)
        c = info['earth_lit_centroid']
        pa = float('nan')
        if c is not None:
            off = c - (c @ g['dir']) * g['dir']
            pa = math.degrees(math.acos(np.clip(off @ toward / np.linalg.norm(off), -1, 1)))
        print(f'  frame {k:3d}: phase angle {ref[k, 3]:6.1f} deg, illuminated fraction SPICE {ref[k, 4]:.3f}, '
              f'drawn {info["earth_lit_share"]:.3f}; lit side toward the Sun within {pa:.1f} deg')
    print(f'  ({in_view} checked frames have the whole disk of the Earth in view)')


# ----------------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_dir', help='output directory of examples/lunar_south_pole_animation.yaml')
    ap.add_argument('--ring', required=True, help='the distant-terrain ring of the run (PLY)')
    ap.add_argument('--kernels', default='data/kernels/moon.tm', help='SPICE metakernel (run from the repository root)')
    ap.add_argument('--out', default='lunar_south_pole.gif', help='the GIF to write')
    ap.add_argument('--mp4', help='also write an MP4 at --mp4-dpi')
    ap.add_argument('--png', help='also write one frame, chosen by --png-step (at --mp4-dpi when --mp4 is given)')
    ap.add_argument('--png-step', type=int, default=40)
    ap.add_argument('--frames', type=int, default=150, help='outputs in the solar day that is drawn')
    ap.add_argument('--preview', help='only write these outputs (comma separated) as PNG to --frames-dir')
    ap.add_argument('--frames-dir', help='where to keep the frames (default: a temporary directory)')
    ap.add_argument('--dpi', type=float, default=100.0, help='resolution of the GIF frames (100 gives 1000 x 500 pixels)')
    ap.add_argument('--mp4-dpi', type=float, default=200.0)
    ap.add_argument('--delay-ms', type=int, default=60, help='time between frames per output (the GIF multiplies it by --gif-step)')
    ap.add_argument('--gif-step', type=int, default=2, help='draw every this many outputs in the GIF (the MP4 draws every output)')
    ap.add_argument('--ss', type=int, default=SUPERSAMPLE, help='supersampling of the picture')
    ap.add_argument('--azimuth', type=float, default=VIEW_AZIMUTH)
    ap.add_argument('--distance-km', type=float, default=CAMERA_DISTANCE_KM)
    ap.add_argument('--height-km', type=float, default=CAMERA_HEIGHT_KM)
    ap.add_argument('--pitch', type=float, default=CAMERA_PITCH)
    ap.add_argument('--fov', type=float, default=FIELD_OF_VIEW)
    for key, val in STYLE.items():
        ap.add_argument('--' + key, type=float, default=val)
    args = ap.parse_args()
    for key in STYLE:
        STYLE[key] = getattr(args, key)

    if args.frames % args.gif_step:
        raise SystemExit(f'--gif-step {args.gif_step} must divide the {args.frames} frames, or the loop of the GIF '
                         f'would close on a different interval')
    if args.png and not args.mp4 and args.png_step % args.gif_step:
        raise SystemExit(f'--png-step {args.png_step} is not one of the outputs of the GIF (every {args.gif_step}); '
                         f'choose a multiple of --gif-step or add --mp4')
    ffmpeg = shutil.which('ffmpeg')
    if args.mp4 and not ffmpeg:
        raise SystemExit('--mp4 needs ffmpeg, which was not found')
    sd = SolarDay(args.run_dir, args.frames)
    n = sd.n
    closure = sd.temperature[:, n] - sd.temperature[:, 0]
    print(f'{len(sd.faces)} faces, {n} frames over {sd.day(n):.4f} days from {sd.utc[0]} UTC (output {sd.first_output})')
    print(f'surface temperature: 1st-99th percentile {np.percentile(sd.temperature[:, :n], 1):.1f}-'
          f'{np.percentile(sd.temperature[:, :n], 99):.1f} K, median {np.median(sd.temperature[:, :n]):.1f} K')
    print(f'loop closure T(end) - T(start): median {np.median(closure):+.3f} K, median |.| '
          f'{np.median(np.abs(closure)):.3f} K, 99th percentile |.| {np.percentile(np.abs(closure), 99):.3f} K, '
          f'largest |.| {np.abs(closure).max():.2f} K')

    def make(dpi):
        cam = Camera(int(round(FIGURE_IN[0] * dpi)) * args.ss, int(round(FIGURE_IN[1] * dpi)) * args.ss,
                     args.azimuth, args.distance_km, args.height_km, args.pitch, args.fov)
        scene = Scene(sd, read_ply(args.ring), cam)
        return cam, scene, sky_geometry(sd.et, cam.position, args.kernels)

    cam, scene, sky = make(args.dpi)
    print(f'camera {args.distance_km:g} km from the pole toward longitude {(args.azimuth + 180) % 360:g}, '
          f'{args.height_km:g} km above the horizontal plane of the pole, looking toward longitude '
          f'{args.azimuth:g} and {args.pitch:g} deg down; field of view {cam.fov:.1f} x {cam.vfov:.1f} deg; '
          f'Earth {2 * math.degrees(sky[0]["radius"]):.2f} deg = '
          f'{2 * math.tan(sky[0]["radius"]) * cam.focal / args.ss:.0f} pixels across')
    fig, draw = build_figure(sd, scene, sky, args.dpi, args.ss)
    if args.preview:
        out = Path(args.frames_dir or '.')
        out.mkdir(parents=True, exist_ok=True)
        for k in [int(s) for s in args.preview.split(',')]:
            info = draw(k)
            fig.savefig(out / f'preview_{k:03d}.png', dpi=args.dpi, facecolor='black')
            print(f'output {k}: Earth {info["earth_px"]} of {info["earth_disk_px"]} pixels in view, Sun {info["sun_px"]}')
        return
    check_earth(sd, sky, scene, args.kernels, range(0, n, 5))
    with tempfile.TemporaryDirectory() as tmp:
        frames_dir = Path(args.frames_dir) if args.frames_dir else Path(tmp)
        frames_dir.mkdir(parents=True, exist_ok=True)
        old = sorted(frames_dir.glob('frame_*.png')) + sorted(frames_dir.glob('hi/frame_*.png'))
        if old:
            print(f'removing {len(old)} frames of an earlier run from {frames_dir}')
            for f in old:
                f.unlink()
        frames, shifts = [], []
        gif_outputs = list(range(0, n, args.gif_step))
        for i, k in enumerate(gif_outputs):
            info = draw(k)
            shifts.append(np.percentile(np.abs(info['shift']), [50, 99]))
            frame = frames_dir / f'frame_{i:03d}.png'
            fig.savefig(frame, dpi=args.dpi, facecolor='black')
            frames.append(frame)
            if args.png and not args.mp4 and k == args.png_step:
                fig.savefig(args.png, dpi=args.dpi, facecolor='black')
        shifts = np.array(shifts)
        print(f'displayed minus face temperature (vertex interpolation): median {np.median(shifts[:, 0]):.2f} K, '
              f'median of the per-frame 99th percentile {np.median(shifts[:, 1]):.1f} K')
        delay = args.delay_ms * args.gif_step
        how = encode_gif(frames, args.out, delay, ffmpeg)
        print(f'wrote {args.out}: {len(frames)} frames (every {args.gif_step} output(s)), {delay} ms each, '
              f'{Path(args.out).stat().st_size / 1e6:.2f} MB ({how})')
        if args.mp4:
            _, hi_scene, hi_sky = make(args.mp4_dpi)
            hi_fig, hi_draw = build_figure(sd, hi_scene, hi_sky, args.mp4_dpi, args.ss)
            hi_dir = frames_dir / 'hi'
            hi_dir.mkdir(exist_ok=True)
            hi_frames = []
            for k in range(n):
                hi_draw(k)
                hi_fig.savefig(hi_dir / f'frame_{k:03d}.png', dpi=args.mp4_dpi, facecolor='black')
                hi_frames.append(hi_dir / f'frame_{k:03d}.png')
                if args.png and k == args.png_step:
                    shutil.copyfile(hi_frames[-1], args.png)
            encode_mp4(hi_frames, args.mp4, 1000.0 / args.delay_ms, ffmpeg)
            print(f'wrote {args.mp4}')


if __name__ == '__main__':
    sys.exit(main())
