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

    F_nu = sum_i  I_i * A_i * mu_i * c_i / r_i**2

where ``I_i`` is the intensity a face emits toward the observer, ``A_i`` the
face area, ``r_i`` the face-to-observer distance, ``mu_i`` the cosine of the
emission angle and ``c_i`` the cosine of the angle between the ray and the
line from the observer to the body center. The flux density is the
irradiance on a plane facing the body, and ``c_i`` differs from 1 only for an
observer within a few body radii. On a convex closed mesh every face with
``mu_i > 0`` is fully visible, so the sum is exact for the faceted surface.
Meshes that are not convex (shape models with concavities, terrain patches)
raise :class:`NonConvexMeshError`, because a face can then be partly hidden
by other terrain.

There are two ways to compute ``I_i``.

surface
    ``I_i = e B_nu(T_i)`` from the surface temperature, with an emissivity
    ``e`` for the observing band. This is right when the emission comes from
    a layer much thinner than the diurnal thermal skin depth, as in the
    thermal infrared.
subsurface
    ``I_i = t(theta_i) * integral of B_nu(T(z)) (alpha / mu_t) exp(-alpha z
    / mu_t) dz``, along the ray refracted into the surface at the face's own
    emission angle, with the Fresnel transmissivity ``t`` of the near-surface
    refractive index and the absorption coefficient ``alpha = 1 / L_e``. It
    needs temperatures at depth (:func:`depth_temperatures_at`) and is
    needed when the electrical skin depth ``L_e`` is comparable to or larger
    than the thermal skin depth, as for ice at millimeter wavelengths. The
    ratio ``beta = L_e / delta_thermal`` is reported to tell the two cases
    apart. An emissivity fitted with the surface mode absorbs the difference
    between the surface and the emitting layer, so it must not be used here.

Brightness temperatures use the nominal solid angle ``pi R**2 / D**2`` (``R``
the mesh's median vertex radius, ``D`` the distance to the body center), so
an isothermal sphere of temperature ``T`` and emissivity 1 has a Planck
brightness temperature of ``T``. A faceted sphere is slightly smaller than
its circumscribed sphere, and an isothermal 100 K icosphere of 5,120 faces
has a brightness temperature 0.12 K low. Two are given: the Planck
brightness temperature (the inverse of the Planck function) and the
Rayleigh-Jeans brightness temperature (``c**2 I / (2 k nu**2)``). At
millimeter wavelengths and Europa temperatures they differ by about
``h nu / 2k`` (5.6 K at 233 GHz). The projected solid angle of the mesh
itself is also reported.

Thermal-infrared bands
----------------------
:func:`band_emission` integrates over a wavelength band instead of a single
frequency. Each visible face emits the band radiance ``e * integral R(lam)
B_lam(T) dlam`` from its surface temperature, with a relative spectral
response ``R`` (:class:`Band`, a top hat between two wavelengths or a
tabulated response), and a gray Lambertian emissivity ``e``. It also
reflects the fraction ``1 - e`` of the sunlight it receives in the band,
from the run's direct sunlight on each face and a 5772 K blackbody Sun. The
result has the band radiance and band brightness temperature of every face
and of the whole visible disk::

    band = observe.Band.top_hat(8.0, 14.0)
    em = observe.band_emission(run, observe.temperature_at(run, geom), geom,
                               band=band, emissivity=0.9,
                               direct_solar=observe.temperature_at(run, geom, "direct_solar"))
    print(em.flux_w_m2, em.tb_k, em.reflected_fraction_disk)

Temperatures at depth
---------------------
:func:`temperature_at_depths` returns the temperature of every face at
given depths below the surface at the emission epoch, from the frames a run
writes with ``output.map_depths_per_step``, interpolated linearly in depth
between the stored levels (:func:`stored_depths` lists them).
:func:`disk_mean_temperature` averages any per-face temperature over the
visible disk::

    dt = observe.temperature_at_depths(run, geom, [0.0, 0.005, 0.02])
    print(observe.disk_mean_temperature(run, dt.temperature_k, geom))

Limitations
-----------
Scattering by grains and voids in the ice, layering of the electrical
properties, surface roughness and beaming, polarization, the instrument
beam, planet light reflected at the observing wavelength, sunlight
reflected at millimeter wavelengths, and the cosmic background are not
included. Resolved images are not yet available from this module.
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
    "geometry_row",
    "DiskMesh", "disk_mesh", "DiskEmission", "disk_emission",
    "Subsurface", "depth_temperatures_at", "thermal_skin_depth_m",
    "StoredDepths", "stored_depths", "DepthTemperatures",
    "temperature_at_depths", "face_weights", "disk_mean_temperature",
    "SUN_TEMPERATURE_K", "STEFAN_BOLTZMANN", "planck_lambda",
    "Band", "BandEmission", "band_emission",
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
    frame of the mesh, and ``distance_km`` is measured from the body center. If a
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
        + "). Run python3 scripts/fetch_kernels.py, or set PTM_DATA_ROOT")


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
                f"{metakernel} lists {entry}, not found at {path}. Run "
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
            f"{run.path}: provenance.txt does not name the body. Use "
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

    Satellites have NAIF ids 101-998 that do not end in 99, and their planet is
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
        ``"interpolate"`` (default) keeps the requested time. The geometry is
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
               f"{offset / 60.0:+.1f} min from the emission epoch, and the body "
               f"rotates {rotation:.2f} deg in that time. Observation "
               f"evaluated at {utc_obs} instead. Use time_match="
               f"'interpolate', a finer output cadence, or start the run so "
               f"that the epoch falls on a step.")
        geom.notes.append(msg)
        warnings.warn(msg, stacklevel=2)
    return geom


def geometry_row(g):
    """The geometry of an observation as a flat dict for a table, with units
    in the keys (the first columns of ``as_row()``)."""
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
    }


def temperature_at(run, geometry, field="temperature"):
    """Per-face ``field`` at the geometry's emission epoch.

    With an interpolated geometry, the two bracketing output steps are
    blended linearly in time. Otherwise the geometry's step is read. The
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
        raise ValueError("this geometry names no output step. Build it with "
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
            f"the mesh has {len(np.unique(label))} separate pieces, and one piece "
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
                f"{float(np.max(h[bad] / edge_len[bad])):.2e} edge lengths), so "
                f"faces turned toward the observer can be hidden by others")


def disk_mesh(source, length_unit=None):
    """Prepare face geometry for :func:`disk_emission`.

    ``source`` is a :class:`ptm.Results` (its ``mesh.ply`` and length unit are
    used) or a dict with ``vertices`` and ``faces`` arrays, in which case
    ``length_unit`` (``"m"`` or ``"km"``) is required. The mesh must be
    centered on the body and convex (see :class:`NonConvexMeshError`).
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
            "origin). A disk integral needs a whole-body mesh")
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
    #: Emergent intensity of each face toward the observer, W m^-2 Hz^-1
    #: sr^-1 (zero for faces turned away).
    face_intensity: Optional[np.ndarray] = None
    #: "surface" or "subsurface". The fields below are for the subsurface mode only.
    mode: str = "surface"
    skin_depth_m: Optional[float] = None
    refractive_index: Optional[float] = None
    #: L_e / delta_thermal at nadir. Much less than 1, the surface mode is
    #: adequate. Of order 1 or more, the emission samples the diurnal wave.
    beta_nadir: Optional[float] = None
    #: Largest weight any visible face gives to the half-space below the
    #: deepest frame (warned above exp(-7)).
    bottom_weight_max: Optional[float] = None

    def as_row(self):
        """A flat dict for a table: geometry and results, with units in the keys."""
        return {
            **geometry_row(self.geometry),
            "distance_km": self.distance_km,
            "frequency_ghz": self.frequency_hz / 1.0e9,
            "flux_density_jy": self.flux_density_jy,
            "tb_planck_k": self.tb_planck_k,
            "tb_rayleigh_jeans_k": self.tb_rj_k,
            "omega_nominal_sr": self.omega_nominal_sr,
            "omega_projected_sr": self.omega_projected_sr,
            "visible_faces": self.visible_faces,
            "mode": self.mode,
            "skin_depth_m": self.skin_depth_m,
            "refractive_index": self.refractive_index,
            "beta_nadir": self.beta_nadir,
        }


def _emission_physics():
    """The public emission scripts (scripts/emission), which hold the one
    implementation of the dielectric model, Fresnel transmissivity and the
    depth integral."""
    import sys
    d = Path(__file__).resolve().parents[3] / "scripts" / "emission"
    if not (d / "depth_weighted_emission.py").is_file():
        raise ImportError(
            f"subsurface emission needs scripts/emission/ from the PTM3D "
            f"repository (looked in {d}). Install the ptm package from a "
            f"checkout with `pip install -e .`")
    if str(d) not in sys.path:
        sys.path.insert(0, str(d))
    import depth_weighted_emission
    import dielectric_mixing
    return depth_weighted_emission, dielectric_mixing


@dataclass
class Subsurface:
    """Electrical properties of the near surface for subsurface emission.

    Give EITHER an electrical skin depth and a refractive index (both constant
    with depth)::

        Subsurface(skin_depth_m=0.07, refractive_index=1.37)

    OR a porosity for the porous-ice dielectric model of
    scripts/emission/dielectric_mixing.py, evaluated at a reference
    temperature (constant with depth)::

        Subsurface(porosity=0.5, reference_temperature_k=100.0)

    or at each node's own temperature, as the per-site script does
    (``reference_temperature_k=None``). With a porosity, ``skin_depth_m`` may
    also be given to replace the model's absorption while keeping its
    refractive index.

    ``skin_depth_m`` is the electrical skin depth L_e = 1/alpha at nadir, the
    depth over which the emitted intensity falls by 1/e. The surface
    transmissivity t(theta) follows from the refractive index by the Fresnel
    equations: it is not a free parameter here.
    """
    skin_depth_m: Optional[float] = None
    refractive_index: Optional[float] = None
    porosity: Optional[float] = None
    reference_temperature_k: Optional[float] = 100.0
    mixing: str = "maxwell_garnett"

    def resolve(self, nu_hz):
        """``(kappa, n_real, temperature_dependent)`` at frequency ``nu_hz``.
        kappa and n_real are None when they depend on temperature."""
        if self.porosity is None:
            if self.skin_depth_m is None or self.refractive_index is None:
                raise ValueError("Subsurface needs skin_depth_m and "
                                 "refractive_index, or a porosity")
            if not (self.skin_depth_m > 0 and self.refractive_index >= 1):
                raise ValueError("skin_depth_m must be positive and "
                                 "refractive_index at least 1")
            return 1.0 / self.skin_depth_m, float(self.refractive_index), False
        if self.refractive_index is not None:
            raise ValueError("give a porosity or a refractive_index, not both")
        _, dm = _emission_physics()
        if self.reference_temperature_k is None:
            if self.skin_depth_m is not None:
                raise ValueError("skin_depth_m replaces the model's absorption "
                                 "and needs a reference_temperature_k")
            return None, None, True
        f_ghz = nu_hz / 1e9
        eps = complex(dm.porous_ice_permittivity(
            f_ghz, self.reference_temperature_k, self.porosity, self.mixing))
        n_real = float(np.real(np.sqrt(eps)))
        kappa = (1.0 / self.skin_depth_m if self.skin_depth_m is not None
                 else float(dm.absorption_coefficient(
                     f_ghz, self.reference_temperature_k, self.porosity,
                     self.mixing)))
        return kappa, n_real, False


def depth_temperatures_at(run, geometry):
    """Per-step temperatures at depth at the geometry's emission epoch.

    Reads the frames written with ``output.map_depths_per_step`` (and the
    depths they were resolved to from statistics/depth_map_z.json), blended
    in time as :func:`temperature_at` does. Returns ``(depths_m, T)`` with T
    of shape (faces, depths), depths increasing from 0.
    """
    import json
    manifest = Path(run.path) / "statistics" / "depth_map_z.json"
    if not manifest.is_file():
        raise FileNotFoundError(
            f"{manifest} not found: run with output.map_depths_m and "
            f"output.map_depths_per_step: true")
    entries = sorted(json.loads(manifest.read_text())["depths"],
                     key=lambda e: e["resolved_m"])
    present = set(run.step_fields)
    missing = [e["field"] for e in entries if e["field"] not in present]
    if missing:
        raise FileNotFoundError(
            f"{run.path}: no per-step frames for {missing}. Run with "
            f"output.map_depths_per_step: true")
    depths = np.array([e["resolved_m"] for e in entries], dtype=np.float64)
    if depths[0] != 0.0:
        raise ValueError(
            f"the shallowest depth frame is at {depths[0]} m. Include 0.0 in "
            f"output.map_depths_m so the column starts at the surface")
    T = np.stack([temperature_at(run, geometry, e["field"]) for e in entries],
                 axis=1)
    return depths, T


def thermal_skin_depth_m(run):
    """The run's diurnal thermal skin depth at the surface (from
    statistics/depth_map_z.json), or None."""
    import json
    manifest = Path(run.path) / "statistics" / "depth_map_z.json"
    if not manifest.is_file():
        return None
    value = json.loads(manifest.read_text()).get("skin_depth_m")
    return float(value) if value else None


# ---------------------------------------------------------------------------
# Temperatures at depth
# ---------------------------------------------------------------------------
#
# A run with output.map_depths_per_step writes, at every output step, the
# temperature of every face at the depths listed in output.map_depths_m. The
# surface temperature frame is the level at depth 0 (the engine's z = 0 frame
# is equal to it bit for bit). A requested depth that is one of the stored
# levels returns that level's values. Any other depth is interpolated linearly
# in depth between the two stored levels around it. Depths below the deepest
# stored level are refused, never extrapolated.


@dataclass
class StoredDepths:
    """The depth levels a run stored at every output step.

    ``depths_m`` starts at 0, the surface temperature frame, and increases.
    ``fields`` names the per-step field of each level. ``requested_m`` and
    ``on_node`` are the engine's record of each level (NaN and None for the
    surface), and ``native_z_m`` is the solver's own depth grid.
    """
    depths_m: np.ndarray
    fields: tuple
    requested_m: np.ndarray
    on_node: tuple
    native_z_m: Optional[np.ndarray] = None
    skin_depth_m: Optional[float] = None


def _depth_config_advice(depths_m):
    listed = ", ".join(f"{z:g}" for z in sorted(set(float(z) for z in depths_m) - {0.0}))
    return (f"Add to the run's configuration\n"
            f"  output:\n"
            f"    map_depths_m: [{listed}]\n"
            f"    map_depth_snap_to_node: false   # store these depths, not the "
            f"nearest grid nodes\n"
            f"    map_depths_per_step: true\n"
            f"and run it again. Depth 0, the surface, needs no depth frames.")


def stored_depths(run):
    """The depth levels ``run`` stored at every output step (:class:`StoredDepths`).

    Raises FileNotFoundError when the run stored no per-step depth frames.
    """
    import json
    from .results import DEPTH_MANIFEST_NAME, depth_of_field
    if "temperature" not in run.step_fields:
        raise FileNotFoundError(f"{run.path}: no per-step temperature frames")
    present = set(run.step_fields)
    manifest_path = Path(run.path) / "statistics" / DEPTH_MANIFEST_NAME
    levels = {}
    native, skin = None, None
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text())
        for e in manifest.get("depths") or []:
            if e.get("field") in present:
                levels[e["field"]] = (float(e["resolved_m"]),
                                      float(e.get("requested_m", math.nan)),
                                      e.get("on_node"))
        if manifest.get("native_z_m"):
            native = np.asarray(manifest["native_z_m"], dtype=np.float64)
        skin = manifest.get("skin_depth_m")
    else:
        names = [f for f in present if depth_of_field(f) is not None]
        if names:
            warnings.warn(
                f"{manifest_path} not found: the depths of the frames "
                f"{sorted(names)} are taken from their names, which round "
                f"them to the nearest millimeter", stacklevel=2)
        for f in names:
            levels[f] = (depth_of_field(f), math.nan, None)
    if not levels:
        raise FileNotFoundError(
            f"{run.path} has no temperature frames at depth (the run did not "
            f"set output.map_depths_per_step).")
    # The surface frame is the level at depth 0. A stored z = 0 frame is the
    # same values and is not listed twice.
    rows = [(0.0, "temperature", math.nan, True)]
    for name, (z, req, node) in sorted(levels.items(), key=lambda kv: kv[1][0]):
        if z == 0.0:
            continue
        if z < 0.0 or not math.isfinite(z):
            raise ValueError(f"{run.path}: frame {name} has depth {z}")
        rows.append((z, name, req, node))
    z = np.array([r[0] for r in rows])
    if np.any(np.diff(z) <= 0):
        raise ValueError(f"{run.path}: two depth frames have the same depth: "
                         f"{[r[1] for r in rows]}")
    return StoredDepths(depths_m=z, fields=tuple(r[1] for r in rows),
                        requested_m=np.array([r[2] for r in rows]),
                        on_node=tuple(r[3] for r in rows), native_z_m=native,
                        skin_depth_m=float(skin) if skin else None)


@dataclass
class DepthTemperatures:
    """Per-face temperatures at requested depths at one observation.

    ``temperature_k`` has shape (depths, faces), in the order of
    ``depths_m``. For each requested depth, ``lower_m`` and ``upper_m`` are the
    stored levels it lies between (equal when the depth is a stored level),
    and ``interpolated`` says whether it was interpolated in depth.
    """
    depths_m: np.ndarray
    temperature_k: np.ndarray
    lower_m: np.ndarray
    upper_m: np.ndarray
    interpolated: np.ndarray
    stored_depths_m: np.ndarray
    geometry: Geometry = field(repr=False)


def _check_depths(depths_m):
    z = np.atleast_1d(np.asarray(depths_m, dtype=np.float64))
    if z.ndim != 1 or z.size == 0:
        raise ValueError("give one or more depths")
    if not np.all(np.isfinite(z)) or np.any(z < 0):
        raise ValueError(f"depths must be finite and not negative (meters "
                         f"below the surface), got {z.tolist()}")
    if len(set(z.tolist())) != z.size:
        raise ValueError(f"depths repeat: {z.tolist()}")
    return z


def temperature_at_depths(run, geometry, depths_m, stored=None):
    """Per-face temperatures at ``depths_m`` (meters, 0 = surface) at the
    geometry's emission epoch (:class:`DepthTemperatures`).

    Each stored level is read and blended in time exactly as
    :func:`temperature_at` does. A depth equal to a stored level returns that
    level. Any other depth is interpolated linearly in depth between the two
    stored levels around it. Depth 0 is the surface temperature frame and
    needs no depth frames. A depth below the deepest stored level, or any
    depth below the surface in a run without depth frames, is an error.
    """
    z = _check_depths(depths_m)
    if np.all(z == 0.0):
        levels = np.array([0.0])
        fields = ("temperature",)
    else:
        try:
            sd = stored if stored is not None else stored_depths(run)
        except FileNotFoundError as exc:
            raise FileNotFoundError(f"{exc} {_depth_config_advice(z)}") from None
        levels, fields = sd.depths_m, sd.fields
        if z.max() > levels[-1]:
            raise ValueError(
                f"depth {z.max():g} m is below the deepest stored level, "
                f"{levels[-1]:g} m. This run stores temperatures from 0 to "
                f"{levels[-1]:g} m (levels "
                f"{', '.join(f'{v:g}' for v in levels)} m) and they are not "
                f"extrapolated. Add deeper depths to output.map_depths_m (the "
                f"solver's column must reach them).")
    cache = {}

    def level(k):
        if k not in cache:
            cache[k] = temperature_at(run, geometry, fields[k])
        return cache[k]

    lower = np.empty(z.size)
    upper = np.empty(z.size)
    interp = np.zeros(z.size, dtype=bool)
    rows = []
    for i, zi in enumerate(z):
        k = int(np.searchsorted(levels, zi, side="right")) - 1
        if levels[k] == zi:
            rows.append(level(k))
            lower[i] = upper[i] = zi
            continue
        w = (zi - levels[k]) / (levels[k + 1] - levels[k])
        rows.append((1.0 - w) * level(k) + w * level(k + 1))
        lower[i], upper[i] = levels[k], levels[k + 1]
        interp[i] = True
    return DepthTemperatures(depths_m=z, temperature_k=np.stack(rows, axis=0),
                             lower_m=lower,
                             upper_m=upper, interpolated=interp,
                             stored_depths_m=np.asarray(levels),
                             geometry=geometry)


def face_weights(mesh_source, geometry, length_unit=None):
    """``(mu, solid_angle_sr)`` of every face seen from ``geometry``: the
    cosine of the emission angle (not positive for faces turned away) and the
    solid angle the face subtends at the observer (zero for faces turned
    away). These are the weights of the disk sums in this module."""
    dm = disk_mesh(mesh_source, length_unit=length_unit)
    mu, _, domega, _ = _disk_weights(dm, geometry)
    return mu, domega


def disk_mean_temperature(mesh_source, temperature, geometry, length_unit=None):
    """Mean of a per-face temperature over the visible disk, weighted by the
    solid angle (projected area over distance squared) of each visible face,
    as ``mean_temperature_k`` of :func:`band_emission`. ``temperature`` has
    shape (faces,) or (n, faces). Returns a float or an array of n."""
    T = np.asarray(temperature, dtype=np.float64)
    _, domega = face_weights(mesh_source, geometry, length_unit=length_unit)
    if T.shape[-1] != domega.size:
        raise ValueError(f"temperature has {T.shape[-1]} values per row for "
                         f"{domega.size} faces")
    omega = float(np.sum(domega))
    if not omega > 0:
        raise ValueError("no face of the mesh faces the observer")
    if T.ndim == 1:
        return float(np.sum(T * domega) / omega)
    return np.array([float(np.sum(row * domega) / omega) for row in T])


def _disk_weights(dm, geometry):
    """Cosine of the emission angle, visibility, and the solid angle each
    visible face subtends at the observer, in sr (zero for faces turned
    away). Also returns the distance to the body center, m."""
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
    return mu, visible, domega, D


def disk_emission(mesh_source, temperature, geometry, *, emissivity=None,
                  frequency_ghz=None, wavelength_um=None, length_unit=None,
                  subsurface=None, depths_m=None, thermal_skin_depth_m=None):
    """Flux density and disk brightness temperature seen from ``geometry``.

    Two modes:

    surface
        ``temperature`` has one value per face (the surface temperature of
        an output step, from :func:`temperature_at`), and each face emits
        ``emissivity * B_nu(T)``. The emissivity is an effective value for
        the observing band.
    subsurface
        ``temperature`` has shape (faces, depths) on ``depths_m`` (from
        :func:`depth_temperatures_at`) and ``subsurface`` gives the
        electrical properties. Each visible face emits the intensity
        t(theta) * integral of B_nu(T(z)) along the ray refracted into the
        surface, at that face's own emission angle theta, where t is the Fresnel
        transmissivity of the refractive index. No emissivity is given: an
        emissivity fitted with the surface mode already absorbs the
        subsurface effect and must not be used here.

    Parameters
    ----------
    mesh_source : ptm.Results, dict, or DiskMesh
        The mesh (see :func:`disk_mesh`).
    temperature : array
        Per face (surface mode) or per face and depth (subsurface mode), K.
    geometry : Geometry
        From :func:`observer_geometry` or :func:`explicit_geometry`.
    emissivity : float or array, one value per face
        Surface mode only: the observing-band emissivity. This is not the
        thermal-infrared emissivity of the energy balance that produced the
        temperatures.
    frequency_ghz, wavelength_um : float
        The observing frequency, given one way or the other.
    subsurface : Subsurface
        Subsurface mode only: the electrical properties.
    depths_m : array
        Subsurface mode only: the depths of ``temperature``'s columns,
        starting at 0. The accuracy of the depth integral is set by this
        list. A warning is issued when the first spacing exceeds 0.16 times
        the shorter of the electrical skin depth along the most oblique
        refracted ray and the thermal skin depth, and when the column reaches
        an optical depth of less than 7 along any ray.
    thermal_skin_depth_m : float
        Subsurface mode, optional: the diurnal thermal skin depth used for the
        diagnostic beta = L_e / delta_thermal (by default read from the run).
    """
    nu = frequency_hz(frequency_ghz, wavelength_um)
    dm = disk_mesh(mesh_source, length_unit=length_unit)
    T = np.asarray(temperature, dtype=np.float64)
    if subsurface is None:
        if emissivity is None:
            raise ValueError("surface mode needs an emissivity (or pass "
                             "subsurface= for the subsurface mode)")
        T = T.ravel()
        if T.size != dm.n_faces:
            raise ValueError(f"temperature has {T.size} values for {dm.n_faces} faces")
    else:
        if emissivity is not None:
            raise ValueError(
                "the subsurface mode takes no emissivity: the surface "
                "transmissivity follows from the refractive index, and an "
                "emissivity fitted with the surface mode already absorbs the "
                "subsurface effect")
        if depths_m is None or T.ndim != 2 or T.shape[0] != dm.n_faces:
            raise ValueError(
                f"the subsurface mode needs temperature of shape (faces, depths) "
                f"= ({dm.n_faces}, n) and depths_m, not {T.shape}")
        depths_m = np.asarray(depths_m, dtype=np.float64)
        if depths_m.shape != (T.shape[1],) or depths_m[0] != 0.0 or np.any(np.diff(depths_m) <= 0):
            raise ValueError("depths_m must start at 0, increase, and match "
                             "temperature's second axis")
    if not np.all(np.isfinite(T)) or np.any(T <= 0):
        bad = int(np.sum(~np.isfinite(T) | (T <= 0)))
        raise ValueError(f"temperature has {bad} non-finite or non-positive values")
    if subsurface is None:
        eps = np.broadcast_to(np.asarray(emissivity, dtype=np.float64), T.shape)
        if np.any(eps <= 0) or np.any(eps > 1) or not np.all(np.isfinite(eps)):
            raise ValueError("emissivity must lie in (0, 1]")

    mu, visible, domega, D = _disk_weights(dm, geometry)
    extra = {}
    if subsurface is None:
        face_intensity = eps * planck_nu(T, nu)
    else:
        face_intensity, extra = _subsurface_intensity(
            T, depths_m, nu, mu, visible, subsurface,
            thermal_skin_depth_m if thermal_skin_depth_m is not None
            else (thermal_skin_depth_m_of(mesh_source)))
    face_flux = face_intensity * domega                  # W m^-2 Hz^-1
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
        geometry=geometry, face_intensity=np.where(visible, face_intensity, 0.0),
        mode="surface" if subsurface is None else "subsurface", **extra)


def thermal_skin_depth_m_of(source):
    """The run's thermal skin depth when ``source`` is a run, else None."""
    if hasattr(source, "path") and hasattr(source, "step_fields"):
        return thermal_skin_depth_m(source)
    return None


#: The depth integral treats temperature as linear between frames. Its error
#: grows as (dz / scale)^2. For a 30 K diurnal wave, spacing at 0.16 times
#: the shorter scale keeps it below 0.1 K (test_observe_subsurface.py).
FIRST_SPACING_LIMIT = 0.16


def _warn_if_coarse(depths_m, skin_depth_m, n_real, mu_visible, delta_thermal):
    """Warn when the first depth spacing is too coarse for the shorter of the
    two scales the integral must resolve: the electrical skin depth along the
    most oblique refracted ray, and the thermal skin depth."""
    if mu_visible.size == 0 or len(depths_m) < 2:
        return
    sin2 = 1.0 - float(np.min(mu_visible)) ** 2
    mu_t = math.sqrt(max(0.0, 1.0 - sin2 / n_real ** 2))
    scales = [skin_depth_m * mu_t]
    if delta_thermal:
        scales.append(float(delta_thermal))
    scale = min(scales)
    dz = float(depths_m[1] - depths_m[0])
    if dz > FIRST_SPACING_LIMIT * scale:
        which = ("electrical skin depth along the most oblique refracted ray"
                 if scale == scales[0] else "thermal skin depth")
        warnings.warn(
            f"the first depth spacing, {dz:.4g} m, exceeds "
            f"{FIRST_SPACING_LIMIT} times the {which} ({scale:.4g} m), so "
            f"the depth integral may be off by more than 0.1 K. Add shallower "
            f"depths to output.map_depths_m.", stacklevel=4)


def _subsurface_intensity(T, depths_m, nu, mu, visible, sub, delta_thermal):
    """Per-face emergent intensity (zero for faces turned away)."""
    dwe, dmix = _emission_physics()
    kappa, n_real, tdep = sub.resolve(nu)
    intensity = np.zeros(T.shape[0])
    w_bottom_max, tau_min = 0.0, np.inf
    theta = np.degrees(np.arccos(np.clip(mu, -1.0, 1.0)))
    idx = np.nonzero(visible)[0]
    cols = T[idx]
    if tdep:
        f_ghz = nu / 1e9
        kap = np.asarray(dmix.absorption_coefficient(
            f_ghz, cols, sub.porosity, sub.mixing), dtype=float)
        nr = np.real(np.sqrt(np.asarray(dmix.porous_ice_permittivity(
            f_ghz, cols, sub.porosity, sub.mixing))))
        t = dmix.fresnel_transmissivity(nr[:, 0], theta[idx])
    else:
        kap, nr = kappa, n_real
        t = dmix.fresnel_transmissivity(n_real, theta[idx])
    I, wb, tau = dwe.emergent_intensity(depths_m, cols, nu, kap, nr,
                                        theta[idx], t, warn_shallow=False)
    intensity[idx] = I
    if idx.size:
        w_bottom_max = float(np.max(wb))
        tau_min = float(np.min(tau))
    if w_bottom_max > math.exp(-7.0):
        warnings.warn(
            f"the depth frames reach an optical depth of only {tau_min:.2f} "
            f"along some rays. The deepest frame's temperature stands in for "
            f"everything below it and carries up to {100 * w_bottom_max:.2f} "
            f"percent of a face's intensity. Add deeper depths to "
            f"output.map_depths_m.", stacklevel=3)
    if tdep:
        f_ghz = nu / 1e9
        t_ref = float(np.median(T[visible, 0]))
        L_e = 1.0 / float(dmix.absorption_coefficient(f_ghz, t_ref, sub.porosity, sub.mixing))
        n_out = float(np.real(np.sqrt(dmix.porous_ice_permittivity(
            f_ghz, t_ref, sub.porosity, sub.mixing))))
    else:
        L_e, n_out = 1.0 / kappa, n_real
    _warn_if_coarse(depths_m, L_e, n_out, mu[visible], delta_thermal)
    beta = L_e / delta_thermal if delta_thermal else None
    return intensity, dict(skin_depth_m=L_e, refractive_index=n_out,
                           beta_nadir=beta, bottom_weight_max=w_bottom_max)


# ---------------------------------------------------------------------------
# Thermal-infrared bands
# ---------------------------------------------------------------------------
#
# In the thermal infrared the surface is opaque, and each face emits from its
# surface temperature. An instrument with a relative spectral response R(lam)
# records the band radiance
#
#     L = e * integral R(lam) B_lam(T) dlam,
#
# and the band brightness temperature T_b is the temperature of a blackbody
# with the same band radiance. The emissivity e is gray (the same at every
# wavelength, and so the same as in the energy balance that produced T) and
# Lambertian. By Kirchhoff's law the surface then reflects a fraction 1 - e of
# the sunlight it receives in the band, as a Lambertian reflector.

#: Effective temperature of the Sun (IAU 2015 Resolution B3), K. A blackbody
#: of this temperature and the nominal solar radius gives 1361 W m^-2 at 1 au.
SUN_TEMPERATURE_K = 5772.0
#: Stefan-Boltzmann constant, from the exact SI h, k and c.
STEFAN_BOLTZMANN = (2.0 * math.pi ** 5 * BOLTZMANN_K ** 4
                    / (15.0 * PLANCK_H ** 3 * SPEED_C ** 2))

#: Gauss-Legendre nodes per panel, and the largest ratio of a panel's long to
#: short wavelength. Planck's function varies by at most a factor exp(x/10)
#: across a panel (x = hc / lam k T), so 16 nodes integrate it to rounding
#: error for any temperature above a few kelvin (test_observe_band.py).
_GL_NODES = 16
_PANEL_RATIO = 1.1


def planck_lambda(temperature_k, wavelength_m):
    """Planck spectral radiance B_lam(T) in W m^-2 sr^-1 per meter of
    wavelength."""
    T = np.asarray(temperature_k, dtype=np.float64)
    lam = np.asarray(wavelength_m, dtype=np.float64)
    x = PLANCK_H * SPEED_C / (lam * BOLTZMANN_K * T)
    # Far down the Wien tail expm1 overflows to inf and B to 0, its limit.
    with np.errstate(over="ignore"):
        return 2.0 * PLANCK_H * SPEED_C ** 2 / lam ** 5 / np.expm1(x)


def _panels(lo_m, hi_m):
    """Panel edges between two wavelengths, each panel at most _PANEL_RATIO
    long in wavelength ratio."""
    n = max(1, int(math.ceil(math.log(hi_m / lo_m) / math.log(_PANEL_RATIO))))
    return lo_m * (hi_m / lo_m) ** (np.arange(n + 1) / n)


@dataclass(frozen=True)
class Band:
    """A wavelength band and the quadrature that integrates over it.

    Build one with :meth:`top_hat` (a response of 1 between two wavelengths)
    or :meth:`from_file` (a tabulated relative spectral response). The
    response is relative: brightness temperatures do not depend on its scale,
    and band radiances are in W m^-2 sr^-1 per unit of response.
    """
    lo_um: float
    hi_um: float
    #: Quadrature nodes (m) and weights (m, response included).
    nodes_m: np.ndarray = field(repr=False)
    weights_m: np.ndarray = field(repr=False)
    response_file: Optional[str] = None

    @classmethod
    def top_hat(cls, lo_um, hi_um):
        """Response 1 from ``lo_um`` to ``hi_um`` micrometers, 0 elsewhere."""
        lo, hi = float(lo_um), float(hi_um)
        if not (0 < lo < hi and math.isfinite(hi)):
            raise ValueError(f"band limits must satisfy 0 < lo < hi, got "
                             f"{lo_um}, {hi_um} um")
        return cls._build(lo, hi, np.array([lo, hi]), np.array([1.0, 1.0]), None)

    @classmethod
    def from_file(cls, path):
        """A relative spectral response read from a text file of two columns,
        wavelength in micrometers and response, separated by spaces, tabs or
        commas. Lines starting with ``#`` are comments. The response is
        linear between the tabulated wavelengths and zero outside them."""
        rows = []
        for k, line in enumerate(Path(path).read_text().splitlines(), 1):
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.replace(",", " ").split()
            if len(parts) != 2:
                raise ValueError(f"{path}:{k}: expected two columns, "
                                 f"wavelength (um) and response, got {line!r}")
            rows.append((float(parts[0]), float(parts[1])))
        if len(rows) < 2:
            raise ValueError(f"{path}: a response needs at least two rows")
        lam, resp = np.array(rows).T
        if np.any(lam <= 0) or np.any(np.diff(lam) <= 0):
            raise ValueError(f"{path}: wavelengths must be positive and increasing")
        if np.any(resp < 0) or not np.any(resp > 0) or not np.all(np.isfinite(resp)):
            raise ValueError(f"{path}: the response must be finite, not negative, "
                             f"and positive somewhere")
        nz = np.nonzero(resp > 0)[0]
        lo_i, hi_i = max(nz[0] - 1, 0), min(nz[-1] + 1, len(lam) - 1)
        return cls._build(float(lam[lo_i]), float(lam[hi_i]), lam[lo_i:hi_i + 1],
                          resp[lo_i:hi_i + 1], str(path))

    @classmethod
    def _build(cls, lo, hi, lam_um, resp, source):
        # Break points: the tabulated wavelengths (where the response has
        # corners) and geometric panels, so the response is linear and the
        # Planck function smooth on every panel.
        edges = np.union1d(lam_um * 1e-6, _panels(lo * 1e-6, hi * 1e-6))
        x, w = np.polynomial.legendre.leggauss(_GL_NODES)
        a, b = edges[:-1, None], edges[1:, None]
        nodes = (0.5 * (b - a) * x[None, :] + 0.5 * (b + a)).ravel()
        weights = (0.5 * (b - a) * w[None, :]).ravel()
        weights = weights * np.interp(nodes, lam_um * 1e-6, resp)
        return cls(lo_um=lo, hi_um=hi, nodes_m=nodes, weights_m=weights,
                   response_file=source)

    @property
    def label(self):
        if self.response_file:
            return f"response {Path(self.response_file).name} ({self.lo_um:g}-{self.hi_um:g} um)"
        return f"{self.lo_um:g}-{self.hi_um:g} um"

    def radiance(self, temperature_k):
        """Band radiance of a blackbody, integral of R B_lam(T) dlam, in
        W m^-2 sr^-1, with the shape of ``temperature_k``."""
        T = np.asarray(temperature_k, dtype=np.float64)
        B = planck_lambda(T[..., None], self.nodes_m)
        return B @ self.weights_m

    def _dradiance_dT(self, T):
        x = PLANCK_H * SPEED_C / (self.nodes_m * BOLTZMANN_K * T[..., None])
        B = planck_lambda(T[..., None], self.nodes_m)
        # dB/dT = B x e^x / ((e^x - 1) T), written to avoid overflow.
        return (B * x / (-np.expm1(-x)) / T[..., None]) @ self.weights_m

    def brightness_temperature(self, radiance):
        """The temperature whose blackbody band radiance is ``radiance``
        (W m^-2 sr^-1). Zero, negative and non-finite radiances give NaN."""
        L = np.asarray(radiance, dtype=np.float64)
        out = np.full(L.shape, np.nan)
        ok = np.isfinite(L) & (L > 0)
        if not np.any(ok):
            return out if out.ndim else float(out)
        # Start from a table in log radiance, then Newton steps in log
        # radiance, which is close to linear in 1/T on the Wien side.
        t_tab = np.geomspace(1.0, 1.0e5, 2001)
        with np.errstate(divide="ignore", under="ignore"):
            ln_tab = np.log(self.radiance(t_tab))
        good = np.isfinite(ln_tab)
        ln_tab, t_tab = ln_tab[good], t_tab[good]
        target = np.log(L[ok])
        if np.any(target < ln_tab[0]) or np.any(target > ln_tab[-1]):
            raise ValueError(
                f"band radiance outside the range of blackbodies from "
                f"{t_tab[0]:g} to {t_tab[-1]:g} K in the {self.label} band")
        T = np.exp(np.interp(target, ln_tab, np.log(t_tab)))
        for _ in range(50):
            f = np.log(self.radiance(T)) - target
            step = f / (self._dradiance_dT(T) / self.radiance(T))
            T = T - step
            if np.all(np.abs(step) <= 1e-13 * T):
                break
        else:
            raise RuntimeError("band brightness temperature did not converge")
        out[ok] = T
        return out if out.ndim else float(out)

    def solar_fraction(self, sun_temperature_k=SUN_TEMPERATURE_K):
        """Fraction of the Sun's irradiance that falls in the band (weighted
        by the response), for a blackbody Sun: pi * integral R B_lam(T_sun)
        dlam / (sigma T_sun^4)."""
        return float(math.pi * self.radiance(sun_temperature_k)
                     / (STEFAN_BOLTZMANN * sun_temperature_k ** 4))


@dataclass
class BandEmission:
    """Disk-integrated thermal-infrared emission in one band.

    Fluxes are irradiances at the observer on a plane facing the body, in
    W m^-2 (per unit of the response). Radiances are in W m^-2 sr^-1. Face
    arrays have one value per face, zero (radiances) or NaN (temperatures)
    for faces turned away from the observer. ``reflected_fraction_max`` (and
    ``reflected_fraction_disk``) is NaN when a visible sunlit face (or the
    whole disk) has zero thermal band radiance, which happens only below a
    few kelvin, and a warning is issued.
    """
    band: Band
    emissivity: object
    flux_w_m2: float
    thermal_flux_w_m2: float
    reflected_flux_w_m2: float
    mean_radiance: float
    tb_k: float
    tb_thermal_k: float
    mean_temperature_k: float
    omega_projected_sr: float
    omega_nominal_sr: float
    radius_km: float
    distance_km: float
    visible_faces: int
    reflected: bool
    reflected_fraction_disk: float
    reflected_fraction_max: float
    solar_fraction: float
    face_mu: np.ndarray = field(repr=False)
    face_solid_angle_sr: np.ndarray = field(repr=False)
    face_radiance_thermal: np.ndarray = field(repr=False)
    face_radiance_reflected: np.ndarray = field(repr=False)
    face_tb_k: np.ndarray = field(repr=False)
    geometry: Geometry = field(repr=False)

    @property
    def face_radiance(self):
        return self.face_radiance_thermal + self.face_radiance_reflected

    def as_row(self):
        """A flat dict for a table: geometry and results, with units in the keys."""
        return {
            **geometry_row(self.geometry),
            "distance_km": self.distance_km,
            "band_lo_um": self.band.lo_um,
            "band_hi_um": self.band.hi_um,
            "response_file": self.band.response_file,
            "flux_w_m2": self.flux_w_m2,
            "thermal_flux_w_m2": self.thermal_flux_w_m2,
            "reflected_flux_w_m2": self.reflected_flux_w_m2,
            "mean_radiance_w_m2_sr": self.mean_radiance,
            "tb_k": self.tb_k,
            "tb_thermal_k": self.tb_thermal_k,
            "mean_temperature_k": self.mean_temperature_k,
            "reflected_fraction_disk": self.reflected_fraction_disk,
            "reflected_fraction_max_face": self.reflected_fraction_max,
            "omega_projected_sr": self.omega_projected_sr,
            "omega_nominal_sr": self.omega_nominal_sr,
            "visible_faces": self.visible_faces,
        }


def band_emission(mesh_source, temperature, geometry, *, band, emissivity,
                  direct_solar=None, reflected=True, length_unit=None):
    """Thermal-infrared emission in a band, per face and over the visible disk.

    Parameters
    ----------
    mesh_source : ptm.Results, dict, or DiskMesh
        The mesh (see :func:`disk_mesh`).
    temperature : array
        Surface temperature of each face, K (from :func:`temperature_at`).
    geometry : Geometry
        From :func:`observer_geometry` or :func:`explicit_geometry`.
    band : Band
        The wavelength band and response.
    emissivity : float or array, one value per face
        Gray Lambertian emissivity. The in-band reflectance is 1 - emissivity.
    direct_solar : array
        Direct sunlight incident on each face at the emission epoch, W m^-2,
        as the run records it (``temperature_at(run, geom, "direct_solar")``),
        which includes the incidence angle and any eclipse. Required unless
        ``reflected=False``.
    reflected : bool
        Add the sunlight reflected in the band. The in-band solar irradiance
        of a face is its direct sunlight times the fraction of a 5772 K
        blackbody's irradiance that falls in the band.

    The disk-integrated flux is F = sum_j L_j A_j mu_j / r_j^2 over the
    visible faces (with the same small correction for a nearby observer as
    :func:`disk_emission`), and the mean radiance is F divided by the solid
    angle of the visible faces, so an isothermal blackbody mesh of any
    resolution has a band brightness temperature equal to its temperature.
    """
    dm = disk_mesh(mesh_source, length_unit=length_unit)
    if not isinstance(band, Band):
        raise TypeError("band must be an observe.Band (Band.top_hat or Band.from_file)")
    T = np.asarray(temperature, dtype=np.float64).ravel()
    if T.size != dm.n_faces:
        raise ValueError(f"temperature has {T.size} values for {dm.n_faces} faces")
    if not np.all(np.isfinite(T)) or np.any(T <= 0):
        bad = int(np.sum(~np.isfinite(T) | (T <= 0)))
        raise ValueError(f"temperature has {bad} non-finite or non-positive values")
    eps = np.broadcast_to(np.asarray(emissivity, dtype=np.float64), T.shape)
    if np.any(eps <= 0) or np.any(eps > 1) or not np.all(np.isfinite(eps)):
        raise ValueError("emissivity must lie in (0, 1]")
    f_sun = band.solar_fraction()
    if reflected:
        if direct_solar is None:
            raise ValueError(
                "reflected sunlight needs direct_solar, the direct sunlight on "
                "each face (temperature_at(run, geom, 'direct_solar')), or "
                "reflected=False to leave it out")
        E = np.asarray(direct_solar, dtype=np.float64).ravel()
        if E.size != dm.n_faces or not np.all(np.isfinite(E)) or np.any(E < 0):
            raise ValueError("direct_solar must have one finite, non-negative "
                             "value per face")

    mu, visible, domega, D = _disk_weights(dm, geometry)
    L_th = np.where(visible, eps * band.radiance(T), 0.0)
    if reflected:
        L_r = np.where(visible, (1.0 - eps) * E * f_sun / math.pi, 0.0)
    else:
        L_r = np.zeros_like(L_th)
    omega = float(np.sum(domega))
    if not omega > 0:
        raise ValueError("no face of the mesh faces the observer")
    F_th = float(np.sum(L_th * domega))
    F_r = float(np.sum(L_r * domega))
    F = F_th + F_r
    face_tb = np.full(T.shape, np.nan)
    face_tb[visible] = band.brightness_temperature(L_th[visible] + L_r[visible])
    # Reflected over thermal radiance. Below a few kelvin the thermal band
    # radiance underflows to zero, and the ratio of a sunlit face is then
    # undefined: NaN, with a warning, never inf. An unlit face at zero
    # thermal radiance reflects nothing and has ratio 0.
    dark = visible & (L_th <= 0.0)
    undefined = dark & (L_r > 0.0)
    ratio = np.zeros_like(L_th)
    lit_th = visible & ~dark
    ratio[lit_th] = L_r[lit_th] / L_th[lit_th]
    ratio[undefined] = np.nan
    if np.any(undefined):
        warnings.warn(
            f"{int(np.sum(undefined))} visible sunlit faces have zero thermal "
            f"band radiance (temperatures down to {T[undefined].min():.3g} K): "
            f"their reflected/thermal ratio, and the maximum over faces, are NaN",
            stacklevel=2)
    if F_th > 0:
        fraction_disk = F_r / F_th
    else:
        fraction_disk = math.nan if F_r > 0 else 0.0
    return BandEmission(
        band=band, emissivity=emissivity, flux_w_m2=F, thermal_flux_w_m2=F_th,
        reflected_flux_w_m2=F_r, mean_radiance=F / omega,
        tb_k=float(band.brightness_temperature(F / omega)),
        tb_thermal_k=float(band.brightness_temperature(F_th / omega)),
        mean_temperature_k=float(np.sum(T * domega) / omega),
        omega_projected_sr=omega,
        omega_nominal_sr=math.pi * dm.radius_m ** 2 / D ** 2,
        radius_km=dm.radius_m / 1000.0, distance_km=D / 1000.0,
        visible_faces=int(np.sum(visible)), reflected=bool(reflected),
        reflected_fraction_disk=fraction_disk,
        reflected_fraction_max=float(np.max(ratio[visible])),
        solar_fraction=f_sun, face_mu=np.where(visible, mu, 0.0),
        face_solid_angle_sr=domega, face_radiance_thermal=L_th,
        face_radiance_reflected=L_r, face_tb_k=face_tb, geometry=geometry)
