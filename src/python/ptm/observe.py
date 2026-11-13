"""Thermal emission of a simulated body as seen by a distant observer.

Given per-face temperatures from a PTM3D run (or any per-face temperature
array on the run's mesh), an observer, and an observing frequency, this module
computes the flux density integrated over the visible disk and the disk
brightness temperature, for comparison with a telescope or spacecraft
measurement.

Three steps::

    from ptm import Results
    from ptm import observe

    run = Results("output/quickstart_europa_diurnal")
    geom = observe.observer_geometry(run, "2025-01-01T12:00:00", observer="EARTH")
    T = observe.temperature_at(run, geom)
    em = observe.disk_emission(run, T, geom, frequency_ghz=233.0, emissivity=0.9)
    print(em.flux_density_jy, em.tb_planck_k, em.tb_rj_k)

Geometry
--------
:func:`observer_geometry` takes the time at which the observer *receives* the
light. It corrects for light time and stellar aberration (SPICE ``CN+S``), so
the body orientation is that of the emission epoch, and
:func:`temperature_at` interpolates the surface temperatures linearly in time
between the two output steps that bracket that epoch. A time whose emission
epoch falls outside the run's output window is an error. With
``time_match="nearest"`` the observation is instead moved by up to half an
output step so that the emission epoch falls on a step, and the shift is
reported in seconds and as the angle the body rotates in that time.

:func:`explicit_geometry` takes a sub-observer latitude, longitude and
distance instead and needs no SPICE kernels, so it works for any body,
including generic bodies whose rotation is not a SPICE frame.

Emission
--------
:func:`disk_emission` sums over the mesh faces that face the observer::

    F_nu = sum_i  e_i * B_nu(T_i) * A_i * mu_i * c_i / r_i**2

where ``A_i`` is the face area, ``r_i`` the face-to-observer distance,
``mu_i`` the cosine of the emission angle and ``c_i`` the cosine of the angle
between the ray and the line from the observer to the body center (the flux
density is the irradiance on a plane facing the body; ``c_i`` differs from 1
only for an observer within a few body radii). Each face emits as a Lambertian
gray body with the observing-band emissivity ``e_i``. On a convex closed mesh
every face with ``mu_i > 0`` is fully visible, so the sum is exact for the
faceted surface. Meshes that are not convex (shape models with concavities,
terrain patches) are refused with :class:`NonConvexMeshError`, because a face
can then be partly hidden by other terrain.

The temperature is an argument, not something read here: the surface
temperature of an output step today, or any other per-face field (a
depth-weighted or at-depth temperature) without a change to this interface.

Brightness temperatures use the nominal solid angle ``pi R**2 / D**2`` (``R``
the mesh's median vertex radius, ``D`` the distance to the body center), so an
isothermal sphere of temperature ``T`` and emissivity 1 has a Planck brightness
temperature of exactly ``T``. Two are given: the Planck brightness temperature
(the inverse of the Planck function) and the Rayleigh-Jeans brightness
temperature (``c**2 I / (2 k nu**2)``). At millimeter wavelengths and Europa
temperatures they differ by about ``h nu / 2k`` (5.6 K at 233 GHz). The
projected solid angle of the mesh itself is also reported.

What is not included
--------------------
The emission is computed from the per-face temperature you pass, as surface
emission with a single emissivity. At millimeter wavelengths the measured
radiation originates below the surface, over a depth set by the absorption
and scattering of the ice. Used with surface temperatures, the emissivity is
therefore an effective value that also absorbs the difference between the
surface and the emitting layer, and it trades off against thermal inertia.
Roughness beaming, polarization and the instrument beam are not included.
"""

from __future__ import annotations

import math
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

__all__ = [
    "PLANCK_H", "BOLTZMANN_K", "SPEED_C", "JANSKY",
    "planck_nu", "planck_nu_temperature", "rayleigh_jeans_temperature",
    "frequency_hz",
    "Geometry", "explicit_geometry", "observer_geometry", "temperature_at",
    "DiskMesh", "disk_mesh", "DiskEmission", "disk_emission",
    "NonConvexMeshError", "ObservationTimeError",
]

#: Exact SI values (2019 redefinition).
PLANCK_H = 6.62607015e-34       # J s
BOLTZMANN_K = 1.380649e-23      # J / K
SPEED_C = 299792458.0           # m / s
#: One jansky in W m^-2 Hz^-1.
JANSKY = 1.0e-26


class NonConvexMeshError(ValueError):
    """The mesh is open, not connected, inconsistently oriented, or concave,
    so faces facing the observer may be partly hidden."""


class ObservationTimeError(ValueError):
    """The requested time cannot be matched to an output step of the run."""


# ---------------------------------------------------------------------------
# Radiometry
# ---------------------------------------------------------------------------

def frequency_hz(frequency_ghz=None, wavelength_um=None):
    """Observing frequency in Hz from exactly one of GHz or micrometers."""
    if (frequency_ghz is None) == (wavelength_um is None):
        raise ValueError("give exactly one of frequency_ghz or wavelength_um")
    if frequency_ghz is not None:
        nu = float(frequency_ghz) * 1.0e9
    else:
        nu = SPEED_C / (float(wavelength_um) * 1.0e-6)
    if not (nu > 0 and math.isfinite(nu)):
        raise ValueError(f"observing frequency must be positive, got {nu} Hz")
    return nu


def planck_nu(temperature_k, nu_hz):
    """Planck spectral radiance B_nu(T) in W m^-2 Hz^-1 sr^-1."""
    T = np.asarray(temperature_k, dtype=np.float64)
    x = PLANCK_H * nu_hz / (BOLTZMANN_K * T)
    return 2.0 * PLANCK_H * nu_hz ** 3 / SPEED_C ** 2 / np.expm1(x)


def planck_nu_temperature(radiance, nu_hz):
    """Inverse of :func:`planck_nu`: the temperature whose B_nu is ``radiance``."""
    I = np.asarray(radiance, dtype=np.float64)
    return (PLANCK_H * nu_hz / BOLTZMANN_K) / np.log1p(
        2.0 * PLANCK_H * nu_hz ** 3 / (SPEED_C ** 2 * I))


def rayleigh_jeans_temperature(radiance, nu_hz):
    """Rayleigh-Jeans brightness temperature c^2 I / (2 k nu^2)."""
    I = np.asarray(radiance, dtype=np.float64)
    return SPEED_C ** 2 * I / (2.0 * BOLTZMANN_K * nu_hz ** 2)


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

@dataclass
class Geometry:
    """Where the observer is, in the body-fixed frame of the mesh.

    ``observer_km`` points from the body center to the observer. The SPICE
    fields are ``None`` for an explicit geometry.
    """
    observer_km: np.ndarray
    sub_observer_lat: float
    sub_observer_lon: float
    distance_km: float
    step: Optional[int] = None
    step_utc: Optional[str] = None
    sun_km: Optional[np.ndarray] = None
    phase_angle_deg: Optional[float] = None
    sub_solar_lat: Optional[float] = None
    sub_solar_lon: Optional[float] = None
    # SPICE observation only
    observer: Optional[str] = None
    frame: Optional[str] = None
    utc_requested: Optional[str] = None
    utc_observed: Optional[str] = None
    light_time_s: Optional[float] = None
    step_offset_s: Optional[float] = None
    rotation_offset_deg: Optional[float] = None
    sun_direction_check_deg: Optional[float] = None
    #: For a satellite seen by a SPICE observer: the planet it orbits and
    #: whether that planet hides it ("visible", "partially occulted",
    #: "totally occulted", or "in transit" when the body is in front).
    occulter: Optional[str] = None
    occultation: Optional[str] = None
    #: "interpolate" or "nearest" (see observer_geometry), and for
    #: "interpolate" the bracketing output steps and their time weights.
    time_match: Optional[str] = None
    steps: Optional[tuple] = None
    weights: Optional[tuple] = None
    utc_emission: Optional[str] = None
    notes: list = field(default_factory=list)


def _latlon(vec):
    v = np.asarray(vec, dtype=np.float64)
    n = np.linalg.norm(v)
    return (float(np.degrees(np.arcsin(np.clip(v[2] / n, -1.0, 1.0)))),
            float(np.degrees(np.arctan2(v[1], v[0]))))


def _angle_deg(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    c = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    # atan2 form is accurate for small angles, where acos is not.
    s = np.linalg.norm(np.cross(a, b)) / (np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.degrees(np.arctan2(s, c)))


def _attach_sun(geom, run, step):
    """Add the run's own Sun vector at ``step`` (phase angle, sub-solar point)."""
    if run is None or step is None:
        return
    path = Path(run.path) / "body_positions.csv"
    if not path.is_file():
        geom.notes.append("body_positions.csv not found: no phase angle")
        return
    from . import readers
    pos = readers.read_body_positions(path)
    if len(pos) != len(run.times().et):
        geom.notes.append(
            f"body_positions.csv has {len(pos)} rows for "
            f"{len(run.times().et)} output steps: no phase angle")
        return
    row = pos[int(step)]
    sun = np.array([row["sun_x"], row["sun_y"], row["sun_z"]], dtype=np.float64)
    geom.sun_km = sun
    geom.phase_angle_deg = _angle_deg(sun, geom.observer_km)
    geom.sub_solar_lat, geom.sub_solar_lon = _latlon(sun)


def explicit_geometry(sub_observer_lat, sub_observer_lon, distance_km,
                      run=None, step=None):
    """An observer at a given planetocentric sub-observer point and distance.

    Latitude and east-positive longitude are in degrees, in the body-fixed
    frame of the mesh; ``distance_km`` is measured from the body center. If a
    run and step are given, the phase angle and sub-solar point at that step
    are added from the run's ``body_positions.csv``.
    """
    lat = math.radians(float(sub_observer_lat))
    lon = math.radians(float(sub_observer_lon))
    d = float(distance_km)
    if not d > 0:
        raise ValueError(f"distance_km must be positive, got {distance_km}")
    vec = d * np.array([math.cos(lat) * math.cos(lon),
                        math.cos(lat) * math.sin(lon),
                        math.sin(lat)])
    geom = Geometry(observer_km=vec, sub_observer_lat=float(sub_observer_lat),
                    sub_observer_lon=float(sub_observer_lon), distance_km=d,
                    step=None if step is None else int(step))
    if run is not None and step is not None:
        geom.step_utc = run.times().utc[int(step)]
        _attach_sun(geom, run, int(step))
    return geom


# -- SPICE -------------------------------------------------------------------

def _kernel_directory(run):
    """``$PTM_DATA_ROOT/kernels`` or ``<repo>/data/kernels``, as the engine does."""
    import os
    from . import paths
    candidates = []
    data_root = os.environ.get("PTM_DATA_ROOT")
    if data_root:
        candidates.append(Path(data_root).expanduser() / "kernels")
    for start in (Path(run.path), Path.cwd(), Path(__file__)):
        root = paths.find_repo_root(start)
        if root is not None:
            candidates.append(root / "data" / "kernels")
    for cand in candidates:
        if (cand / "lsk").is_dir():
            return cand
    raise FileNotFoundError(
        "no SPICE kernel directory found (looked for "
        + ", ".join(str(c) for c in candidates)
        + "); run python3 scripts/fetch_kernels.py, or set PTM_DATA_ROOT")


def _metakernel_files(metakernel, kernel_dir):
    """Absolute paths of the kernels a PTM3D metakernel lists.

    The shipped metakernels use repository-relative paths, so they are read
    here and their entries resolved against ``kernel_dir`` rather than passed
    to ``furnsh`` from an arbitrary working directory.
    """
    text = Path(metakernel).read_text()
    data = "\n".join(re.findall(r"\\begindata(.*?)(?:\\begintext|$)", text,
                                flags=re.S))
    block = re.search(r"KERNELS_TO_LOAD\s*=\s*\((.*?)\)", data, flags=re.S)
    if block is None:
        raise ValueError(f"{metakernel}: no KERNELS_TO_LOAD")
    files = []
    for entry in re.findall(r"'([^']*)'", block.group(1)):
        entry = entry.replace("$KERNELS/", "")
        path = (kernel_dir / entry) if not Path(entry).is_absolute() else Path(entry)
        if not path.exists():
            raise FileNotFoundError(
                f"{metakernel} lists {entry}, not found at {path}; run "
                f"python3 scripts/fetch_kernels.py")
        files.append(path)
    return files


class _Kernels:
    """Load a set of kernels for the duration of a ``with`` block."""

    def __init__(self, files):
        self.files = [str(f) for f in files]

    def __enter__(self):
        import spiceypy
        for f in self.files:
            spiceypy.furnsh(f)
        return spiceypy

    def __exit__(self, *exc):
        import spiceypy
        for f in reversed(self.files):
            spiceypy.unload(f)
        return False


def _run_body(run):
    prov = getattr(run, "engine_provenance", {}) or {}
    body = prov.get("body")
    if not body:
        raise ValueError(
            f"{run.path}: provenance.txt does not name the body; use "
            f"explicit_geometry() instead")
    if prov.get("generic_body", "false").lower() == "true":
        raise ValueError(
            f"{run.path}: '{body}' is a generic body whose rotation is the "
            f"run's own spin model, not a SPICE frame, so a SPICE observer "
            f"cannot be placed in its mesh frame. Use explicit_geometry() with "
            f"the sub-observer latitude and longitude instead.")
    return body


_OCCULT_CODES = {0: "visible", 1: "in transit", 2: "in transit",
                 3: "in transit", -1: "partially occulted",
                 -2: "partially occulted", -3: "totally occulted"}


def _planet_occultation(sp, body_u, frame, observer, et_obs):
    """Whether the planet a satellite orbits hides it from the observer.

    Satellites have NAIF ids 101-998 that do not end in 99; their planet is
    the id ending in 99 of the same system. Anything else (a planet, the Sun,
    an asteroid) returns ``(None, None)``.
    """
    try:
        code = int(sp.bodn2c(body_u))
    except Exception:
        return None, None
    if not (100 < code < 1000) or code % 100 == 99:
        return None, None
    planet_code = (code // 100) * 100 + 99
    planet = sp.bodc2n(planet_code)
    try:
        if int(sp.bodn2c(str(observer).upper())) == planet_code:
            return None, None
    except Exception:
        pass
    result = sp.occult(body_u, "ELLIPSOID", frame, planet, "ELLIPSOID",
                       f"IAU_{planet}", "CN", str(observer), et_obs)
    return planet, _OCCULT_CODES[int(result)]


def observer_geometry(run, utc, observer="EARTH", kernels=(), body=None,
                      frame=None, metakernel=None, time_match="interpolate",
                      warn_rotation_deg=1.0):
    """Observer geometry at a received time.

    Parameters
    ----------
    run : ptm.Results
        A finished run of a built-in (SPICE) body.
    utc : str
        Time at which the observer receives the light, any format SPICE
        ``str2et`` accepts.
    observer : str
        SPICE name or NAIF id of the observer (``"EARTH"``, a spacecraft).
    kernels : sequence of paths
        Extra kernels, e.g. a spacecraft SPK, loaded with the run body's
        metakernel.
    body, frame, metakernel :
        Overrides. By default the body comes from the run's ``provenance.txt``,
        the frame is ``IAU_<BODY>`` (the frame the engine uses) and the
        metakernel is ``<kernels>/<body>.tm``.
    time_match : ``"interpolate"`` or ``"nearest"``
        ``"interpolate"`` (default) keeps the requested time; the geometry is
        that of the emission epoch, and :func:`temperature_at` interpolates
        the temperatures linearly in time between the two output steps that
        bracket it. The emission epoch must lie inside the output window.
        ``"nearest"`` instead moves the observation by up to half an output
        step so that the emission epoch falls exactly on the nearest step,
        and reports the shift and the rotation of the body it implies.
    warn_rotation_deg : float
        With ``"nearest"``, warn when the shift rotates the body by more than
        this angle.
    """
    if time_match not in ("interpolate", "nearest"):
        raise ValueError(f"time_match must be 'interpolate' or 'nearest', "
                         f"got {time_match!r}")
    if body is None:
        body = _run_body(run)
    body_u = body.upper()
    frame = frame or f"IAU_{body_u}"
    kdir = _kernel_directory(run)
    mk = Path(metakernel) if metakernel else kdir / f"{body.lower()}.tm"
    if not mk.is_file():
        raise FileNotFoundError(f"metakernel {mk} not found")
    files = _metakernel_files(mk, kdir) + [Path(k).expanduser() for k in kernels]

    et_steps = run.times().et
    utc_steps = run.times().utc
    spacing = float(np.median(np.diff(et_steps))) if len(et_steps) > 1 else 0.0

    with _Kernels(files) as sp:
        et_rx = sp.str2et(utc)

        def solve(et):
            spoint, trgepc, srfvec = sp.subpnt(
                "INTERCEPT/ELLIPSOID", body_u, et, frame, "CN+S", observer)
            return np.asarray(spoint) - np.asarray(srfvec), trgepc

        obs_vec, et_emit = solve(et_rx)
        if time_match == "interpolate":
            lo, hi = et_steps[0], et_steps[-1]
        else:
            lo = et_steps[0] - 0.5 * spacing
            hi = et_steps[-1] + 0.5 * spacing
        if not (lo <= et_emit <= hi):
            raise ObservationTimeError(
                f"{utc}: the light {observer} receives then left {body} at "
                f"{sp.et2utc(et_emit, 'ISOC', 1)} (light time "
                f"{et_rx - et_emit:.1f} s), outside this run's output window "
                f"{utc_steps[0]} .. {utc_steps[-1]}")
        step = int(np.argmin(np.abs(et_steps - et_emit)))
        offset = float(et_steps[step] - et_emit)
        bracket, weights, rotation = None, None, None

        if time_match == "nearest":
            # Shift the observation so the emission epoch is exactly the step.
            et_obs = et_rx
            for _ in range(8):
                et_obs = et_obs + (et_steps[step] - et_emit)
                obs_vec, et_emit = solve(et_obs)
                if abs(et_emit - et_steps[step]) < 1.0e-3:
                    break
            else:
                raise RuntimeError("light-time iteration did not converge")
            m_req = np.asarray(sp.pxform("J2000", frame, et_steps[step] - offset))
            m_use = np.asarray(sp.pxform("J2000", frame, et_steps[step]))
            _, rot = sp.raxisa(m_use @ m_req.T)
            rotation = float(np.degrees(rot))
        else:
            et_obs = et_rx
            k1 = int(np.searchsorted(et_steps, et_emit))
            k1 = min(max(k1, 1), len(et_steps) - 1) if len(et_steps) > 1 else 0
            k0 = max(k1 - 1, 0)
            if k1 == k0:
                weights = (1.0, 0.0)
            else:
                w1 = float((et_emit - et_steps[k0]) / (et_steps[k1] - et_steps[k0]))
                weights = (1.0 - w1, w1)
            bracket = (k0, k1)

        # The run's own Sun vector must agree with SPICE in this frame, at
        # the step it was recorded for.
        sun_check, _ = sp.spkpos("SUN", et_steps[step], frame, "LT+S", body_u)
        # The Sun as seen from the body at the emission epoch.
        sun_emit, _ = sp.spkpos("SUN", et_emit, frame, "LT+S", body_u)
        utc_obs = sp.et2utc(et_obs, "ISOC", 3)
        utc_emit = sp.et2utc(et_emit, "ISOC", 3)

        # A satellite can be hidden by its planet. The disk sum does not know
        # that, so look and say.
        occulter, occultation = _planet_occultation(sp, body_u, frame,
                                                    observer, et_obs)

    lat, lon = _latlon(obs_vec)
    geom = Geometry(
        observer_km=np.asarray(obs_vec, dtype=np.float64),
        sub_observer_lat=lat, sub_observer_lon=lon,
        distance_km=float(np.linalg.norm(obs_vec)),
        step=step, step_utc=utc_steps[step], observer=observer, frame=frame,
        utc_requested=str(utc), utc_observed=utc_obs, utc_emission=utc_emit,
        light_time_s=float(et_obs - et_emit),
        step_offset_s=offset if time_match == "nearest" else None,
        rotation_offset_deg=rotation, occulter=occulter,
        occultation=occultation, steps=bracket, weights=weights,
        time_match=time_match)
    if occultation not in (None, "visible", "in transit"):
        msg = (f"{utc}: {body} is {occultation} by {occulter} as seen from "
               f"{observer}. The flux computed here is that of the whole "
               f"visible disk and does not include the occultation.")
        geom.notes.append(msg)
        warnings.warn(msg, stacklevel=2)
    _attach_sun(geom, run, step)
    if geom.sun_km is not None:
        geom.sun_direction_check_deg = _angle_deg(geom.sun_km, sun_check)
        if geom.sun_direction_check_deg > 0.05:
            msg = (f"the run's Sun direction at step {step} differs from SPICE "
                   f"({frame}) by {geom.sun_direction_check_deg:.3f} deg: the "
                   f"run and these kernels may not share a body-fixed frame")
            geom.notes.append(msg)
            warnings.warn(msg, stacklevel=2)
    if time_match == "interpolate":
        geom.sun_km = np.asarray(sun_emit, dtype=np.float64)
        geom.phase_angle_deg = _angle_deg(geom.sun_km, geom.observer_km)
        geom.sub_solar_lat, geom.sub_solar_lon = _latlon(geom.sun_km)
    if rotation is not None and rotation > warn_rotation_deg:
        msg = (f"{utc}: nearest output step {utc_steps[step]} is "
               f"{offset / 60.0:+.1f} min from the emission epoch; the body "
               f"rotates {rotation:.2f} deg in that time. Observation "
               f"evaluated at {utc_obs} instead. Use time_match="
               f"'interpolate', a finer output cadence, or start the run so "
               f"that the epoch falls on a step.")
        geom.notes.append(msg)
        warnings.warn(msg, stacklevel=2)
    return geom


def temperature_at(run, geometry, field="temperature"):
    """Per-face ``field`` at the geometry's emission epoch.

    With an interpolated geometry, the two bracketing output steps are
    blended linearly in time; otherwise the geometry's step is read. The
    result is float64, one value per face.
    """
    if geometry.weights is not None:
        k0, k1 = geometry.steps
        w0, w1 = geometry.weights
        a = np.asarray(run.field(field, k0), dtype=np.float64)
        if w1 == 0.0:
            return a
        b = np.asarray(run.field(field, k1), dtype=np.float64)
        return w0 * a + w1 * b
    if geometry.step is None:
        raise ValueError("this geometry names no output step; build it with "
                         "explicit_geometry(..., run=run, step=N)")
    return np.asarray(run.field(field, geometry.step), dtype=np.float64)


# ---------------------------------------------------------------------------
# Mesh
# ---------------------------------------------------------------------------

@dataclass
class DiskMesh:
    """Face geometry in meters with outward normals, checked to be convex."""
    centers_m: np.ndarray
    normals: np.ndarray
    areas_m2: np.ndarray
    radius_m: float
    n_faces: int
    normals_flipped: bool


def _check_closed_convex(vertices, faces, normals, centers):
    """Raise NonConvexMeshError unless the mesh is a closed, connected,
    consistently oriented, locally convex surface (hence convex)."""
    F = len(faces)
    a = faces
    directed = np.concatenate([a[:, [0, 1]], a[:, [1, 2]], a[:, [2, 0]]])
    opposite = np.concatenate([a[:, 2], a[:, 0], a[:, 1]])
    owner = np.tile(np.arange(F), 3)

    key = np.sort(directed, axis=1)
    order = np.lexsort((key[:, 1], key[:, 0]))
    key_s = key[order]
    _, counts = np.unique(key_s, axis=0, return_counts=True)
    if np.any(counts != 2):
        n_open = int(np.sum(counts == 1))
        n_multi = int(np.sum(counts > 2))
        raise NonConvexMeshError(
            f"the mesh is not a closed surface ({n_open} boundary edges, "
            f"{n_multi} edges shared by more than two faces)")
    first, second = order[0::2], order[1::2]
    # Consistent winding: the two faces traverse a shared edge in opposite
    # directions.
    if np.any(np.all(directed[first] == directed[second], axis=1)):
        raise NonConvexMeshError("the mesh's faces are not consistently oriented")

    # Connected: propagate the smallest face label across shared edges.
    label = np.arange(F)
    f1, f2 = owner[first], owner[second]
    for _ in range(F):
        m = np.minimum(label[f1], label[f2])
        new = label.copy()
        np.minimum.at(new, f1, m)
        np.minimum.at(new, f2, m)
        new = new[new]
        if np.array_equal(new, label):
            break
        label = new
    if len(np.unique(label)) != 1:
        raise NonConvexMeshError(
            f"the mesh has {len(np.unique(label))} separate pieces; one piece "
            f"can hide another")

    # Locally convex: across every edge, the far vertex of the neighbor lies
    # on or below the face's plane.
    edge_len = np.linalg.norm(vertices[directed[first, 0]]
                              - vertices[directed[first, 1]], axis=1)
    tol = 1.0e-6 * edge_len
    for fa, fb in ((first, second), (second, first)):
        h = np.einsum("ij,ij->i",
                      vertices[opposite[fb]] - centers[owner[fa]],
                      normals[owner[fa]])
        bad = h > tol
        if np.any(bad):
            raise NonConvexMeshError(
                f"the mesh is concave at {int(np.sum(bad))} edges (largest "
                f"height above a neighboring face "
                f"{float(np.max(h[bad] / edge_len[bad])):.2e} edge lengths); "
                f"faces turned toward the observer can be hidden by others")


def disk_mesh(source, length_unit=None):
    """Prepare face geometry for :func:`disk_emission`.

    ``source`` is a :class:`ptm.Results` (its ``mesh.ply`` and length unit are
    used) or a dict with ``vertices`` and ``faces`` arrays, in which case
    ``length_unit`` (``"m"`` or ``"km"``) is required. The mesh must be
    centered on the body and convex; see :class:`NonConvexMeshError`.
    """
    if isinstance(source, DiskMesh):
        return source
    if hasattr(source, "mesh") and hasattr(source, "mesh_units"):
        cached = getattr(source, "_disk_mesh", None)
        if cached is not None:
            return cached
        mesh = source.mesh
        unit = source.mesh_units
    else:
        mesh = source
        unit = length_unit
    if unit not in ("m", "km"):
        raise ValueError(f"length_unit must be 'm' or 'km', got {unit!r}")
    scale = 1000.0 if unit == "km" else 1.0
    vertices = np.asarray(mesh["vertices"], dtype=np.float64) * scale
    faces = np.asarray(mesh["faces"], dtype=np.int64)
    v0, v1, v2 = vertices[faces[:, 0]], vertices[faces[:, 1]], vertices[faces[:, 2]]
    centers = (v0 + v1 + v2) / 3.0
    cross = np.cross(v1 - v0, v2 - v0)
    twice_area = np.linalg.norm(cross, axis=1)
    if np.any(twice_area <= 0):
        raise ValueError(f"{int(np.sum(twice_area <= 0))} faces have zero area")
    normals = cross / twice_area[:, None]
    areas = 0.5 * twice_area

    radii = np.linalg.norm(vertices, axis=1)
    if np.min(np.linalg.norm(centers, axis=1)) < 0.05 * np.ptp(vertices, axis=0).max():
        raise NonConvexMeshError(
            "the mesh is not centered on the body (a face lies near the "
            "origin); a disk integral needs a whole-body mesh")
    # Outward orientation from the signed volume.
    volume = np.sum(np.einsum("ij,ij->i", centers, normals) * areas) / 3.0
    flipped = volume < 0
    if flipped:
        normals = -normals
    _check_closed_convex(vertices, faces, normals, centers)
    dm = DiskMesh(centers_m=centers, normals=normals, areas_m2=areas,
                  radius_m=float(np.median(radii)), n_faces=len(faces),
                  normals_flipped=bool(flipped))
    if hasattr(source, "mesh") and hasattr(source, "mesh_units"):
        try:
            source._disk_mesh = dm
        except AttributeError:
            pass
    return dm


# ---------------------------------------------------------------------------
# Disk integral
# ---------------------------------------------------------------------------

@dataclass
class DiskEmission:
    """Disk-integrated emission at one observing frequency.

    ``face_flux_jy`` is each face's contribution (zero for faces turned away).
    """
    frequency_hz: float
    flux_density_jy: float
    tb_planck_k: float
    tb_rj_k: float
    omega_nominal_sr: float
    omega_projected_sr: float
    radius_km: float
    distance_km: float
    visible_faces: int
    face_flux_jy: np.ndarray
    geometry: Geometry

    def as_row(self):
        """A flat dict for a table: geometry and results, with units in the keys."""
        g = self.geometry
        return {
            "utc_requested": g.utc_requested,
            "utc_observed": g.utc_observed,
            "step": g.step,
            "step_utc": g.step_utc,
            "utc_emission": g.utc_emission,
            "time_match": g.time_match,
            "steps": (f"{g.steps[0]}-{g.steps[1]}" if g.steps is not None
                      else None if g.step is None else f"{g.step}"),
            "step_weights": None if g.weights is None else f"{g.weights[0]:.4f},{g.weights[1]:.4f}",
            "step_offset_s": g.step_offset_s,
            "rotation_offset_deg": g.rotation_offset_deg,
            "light_time_s": g.light_time_s,
            "sub_observer_lat_deg": g.sub_observer_lat,
            "sub_observer_lon_deg": g.sub_observer_lon,
            "sub_solar_lat_deg": g.sub_solar_lat,
            "sub_solar_lon_deg": g.sub_solar_lon,
            "phase_angle_deg": g.phase_angle_deg,
            "occultation": g.occultation,
            "distance_km": self.distance_km,
            "frequency_ghz": self.frequency_hz / 1.0e9,
            "flux_density_jy": self.flux_density_jy,
            "tb_planck_k": self.tb_planck_k,
            "tb_rayleigh_jeans_k": self.tb_rj_k,
            "omega_nominal_sr": self.omega_nominal_sr,
            "omega_projected_sr": self.omega_projected_sr,
            "visible_faces": self.visible_faces,
        }


def disk_emission(mesh_source, temperature, geometry, *, emissivity,
                  frequency_ghz=None, wavelength_um=None, length_unit=None):
    """Flux density and disk brightness temperature seen from ``geometry``.

    Parameters
    ----------
    mesh_source : ptm.Results, dict, or DiskMesh
        The mesh (see :func:`disk_mesh`).
    temperature : array, one value per face, kelvin
        The temperature that emits: an output step's surface temperature, or
        any other per-face field.
    geometry : Geometry
        From :func:`observer_geometry` or :func:`explicit_geometry`.
    emissivity : float or array, one value per face
        Observing-band emissivity. This is not the thermal-infrared emissivity
        of the energy balance that produced the temperatures.
    frequency_ghz, wavelength_um : float
        The observing frequency, given one way or the other.
    """
    nu = frequency_hz(frequency_ghz, wavelength_um)
    dm = disk_mesh(mesh_source, length_unit=length_unit)
    T = np.asarray(temperature, dtype=np.float64).ravel()
    if T.size != dm.n_faces:
        raise ValueError(f"temperature has {T.size} values for {dm.n_faces} faces")
    if not np.all(np.isfinite(T)) or np.any(T <= 0):
        bad = int(np.sum(~np.isfinite(T) | (T <= 0)))
        raise ValueError(f"temperature has {bad} non-finite or non-positive values")
    eps = np.broadcast_to(np.asarray(emissivity, dtype=np.float64), T.shape)
    if np.any(eps <= 0) or np.any(eps > 1) or not np.all(np.isfinite(eps)):
        raise ValueError("emissivity must lie in (0, 1]")

    P = np.asarray(geometry.observer_km, dtype=np.float64) * 1000.0
    D = float(np.linalg.norm(P))
    if not D > dm.radius_m:
        raise ValueError("the observer is inside the body")
    r = P[None, :] - dm.centers_m
    dist = np.linalg.norm(r, axis=1)
    mu = np.einsum("ij,ij->i", dm.normals, r) / dist
    visible = mu > 0
    # Flux density is the irradiance on a plane facing the body center, so each
    # ray is weighted by the cosine of its angle to that line. The factor is
    # 1 - O((R/D)^2): negligible from Earth, 3 percent at three body radii.
    cos_los = (r @ P) / (dist * D)
    domega = np.where(visible, dm.areas_m2 * mu * cos_los / dist ** 2, 0.0)
    face_flux = eps * planck_nu(T, nu) * domega          # W m^-2 Hz^-1
    flux = float(np.sum(face_flux))
    omega_nom = math.pi * dm.radius_m ** 2 / D ** 2
    intensity = flux / omega_nom
    return DiskEmission(
        frequency_hz=nu, flux_density_jy=flux / JANSKY,
        tb_planck_k=float(planck_nu_temperature(intensity, nu)),
        tb_rj_k=float(rayleigh_jeans_temperature(intensity, nu)),
        omega_nominal_sr=omega_nom, omega_projected_sr=float(np.sum(domega)),
        radius_km=dm.radius_m / 1000.0, distance_km=D / 1000.0,
        visible_faces=int(np.sum(visible)), face_flux_jy=face_flux / JANSKY,
        geometry=geometry)
