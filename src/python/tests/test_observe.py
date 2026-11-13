#!/usr/bin/env python3
"""Checks of ptm.observe: disk-integrated thermal emission at an observer.

    python3 src/python/tests/test_observe.py [run_dir]

Every check is against a closed form, a hand calculation, or a direct SPICE
call, never against the module itself:

1. isothermal sphere: the projected solid angle of an icosphere converges to
   pi R^2 / D^2 at second order, the brightness temperature of an emissivity-1
   sphere converges to its temperature, and the flux of a sphere at a finite
   distance converges to pi R^2 e B / D^2;
2. a latitude temperature law seen from the equator and from the pole agrees
   with one-dimensional integrals evaluated by Gauss-Legendre quadrature, with
   second-order convergence;
3. (SPICE) the sub-Earth point on Europa at 2025-01-01T00:00:00 UTC matches a
   direct subpnt call and the values recorded when this module was designed,
   and the Sun direction a Europa run recorded matches SPICE; on a second
   date, 2025-01-08, the occultation flag agrees with an independent
   classification from the Earth-Europa and Earth-Jupiter vectors and the
   bodies' radii, before, during and after Jupiter hides Europa;
4. units: the Planck function and a Europa-sized flux density match a hand
   calculation, and the frequency and wavelength routes agree;
5. time matching: nearest step, reported offset, refusal outside the run;
6. meshes that are open, dented, or in two pieces are refused.

Checks 3 and 5 need the Europa kernels (python3 scripts/fetch_kernels.py
--bodies europa); 3b and 5 also need a run of
examples/quickstart_europa_diurnal.yaml, by default in
output/quickstart_europa_diurnal. Missing prerequisites are reported as SKIP
and make the script exit 2, so a skip is never mistaken for a pass.
"""

import math
import os
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "src" / "python"))

from ptm import observe  # noqa: E402
from ptm.results import Results  # noqa: E402

FAILS = []
SKIPS = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def skip(name, why):
    print(f"  SKIP  {name}  ({why})")
    SKIPS.append(name)


# ---------------------------------------------------------------------------

def icosphere(level, radius=1.0):
    """Icosphere with outward winding, vertices on the sphere."""
    t = (1.0 + 5 ** 0.5) / 2.0
    v = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0),
         (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t),
         (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    f = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
         (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
         (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
         (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    verts = [np.array(p, float) / np.linalg.norm(p) for p in v]
    faces = f
    for _ in range(level):
        cache = {}

        def mid(a, b):
            key = (min(a, b), max(a, b))
            if key not in cache:
                m = verts[a] + verts[b]
                verts.append(m / np.linalg.norm(m))
                cache[key] = len(verts) - 1
            return cache[key]
        new = []
        for a, b, c in faces:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            new += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        faces = new
    return {"vertices": np.array(verts) * radius, "faces": np.array(faces)}


def face_latitudes(mesh):
    c = mesh["vertices"][mesh["faces"]].mean(axis=1)
    return np.arcsin(c[:, 2] / np.linalg.norm(c, axis=1))


def orders(errors):
    e = np.abs(np.asarray(errors))
    return np.log2(e[:-1] / e[1:])


R_KM = 1560.8
NU_GHZ = 233.0
LEVELS = [2, 3, 4, 5, 6]


def test_isothermal():
    print("1. isothermal sphere")
    T, eps = 100.0, 0.8
    nu = NU_GHZ * 1e9
    far = observe.explicit_geometry(10.0, 33.0, 1.0e6 * R_KM)
    dev, tb_err, tb_pred = [], [], []
    for lv in LEVELS:
        mesh = icosphere(lv, R_KM)
        em = observe.disk_emission(mesh, np.full(len(mesh["faces"]), T), far,
                                   emissivity=eps, frequency_ghz=NU_GHZ,
                                   length_unit="km")
        # Flux over the mesh's own solid angle is e B(T) exactly.
        I_proj = em.flux_density_jy * observe.JANSKY / em.omega_projected_sr
        exact = observe.planck_nu(T, nu) * eps
        if lv == LEVELS[-1]:
            check("flux / projected solid angle = e B(T)",
                  abs(I_proj / exact - 1) < 1e-12, f"rel {I_proj / exact - 1:.1e}")
        dev.append(em.omega_projected_sr / em.omega_nominal_sr - 1)
        em1 = observe.disk_emission(mesh, np.full(len(mesh["faces"]), T), far,
                                    emissivity=1.0, frequency_ghz=NU_GHZ,
                                    length_unit="km")
        tb_err.append(em1.tb_planck_k - T)
        # For an isothermal body the only error is the mesh's solid-angle
        # deficit: B(Tb) = B(T) * Omega_projected / Omega_nominal exactly.
        pred = float(observe.planck_nu_temperature(
            observe.planck_nu(T, nu) * em1.omega_projected_sr / em1.omega_nominal_sr, nu))
        tb_pred.append(abs(em1.tb_planck_k - pred))
    p = orders(dev)
    print("     levels", LEVELS, "solid-angle deficit", ["%.2e" % d for d in dev],
          "orders", ["%.2f" % o for o in p])
    check("projected solid angle -> pi R^2/D^2 at second order",
          np.all(np.abs(p - 2) < 0.1), f"orders {np.round(p, 3)}")
    check("emissivity-1 Planck brightness temperature = B^-1(B(T) x solid-angle ratio)",
          max(tb_pred) < 1e-9, f"max {max(tb_pred):.1e} K")
    check("emissivity-1 Planck brightness temperature -> T at second order",
          np.all(np.abs(orders(tb_err) - 2) < 0.1),
          f"errors {['%.1e' % e for e in tb_err]} K")
    # A finite distance: the flux of an isothermal sphere is still
    # pi R^2 e B / D^2 exactly (it radiates like a point source).
    D = 3.0 * R_KM
    mesh = icosphere(6, R_KM)
    near = observe.explicit_geometry(0.0, 0.0, D)
    em = observe.disk_emission(mesh, np.full(len(mesh["faces"]), T), near,
                               emissivity=eps, frequency_ghz=NU_GHZ,
                               length_unit="km")
    ref = math.pi * R_KM ** 2 / D ** 2 * eps * observe.planck_nu(T, nu) / observe.JANSKY
    check("finite distance D = 3R: flux = pi R^2 e B / D^2",
          abs(em.flux_density_jy / ref - 1) < 1e-4,
          f"rel {em.flux_density_jy / ref - 1:.2e}")


def test_latitude_law():
    print("2. latitude law T = 60 + 70 cos(lat), seen from the equator and the pole")
    nu = NU_GHZ * 1e9
    eps = 0.9
    D = 1.0e8 * R_KM          # parallel-ray limit to ~1e-8
    x, w = np.polynomial.legendre.leggauss(400)

    def B(lat):
        return eps * observe.planck_nu(60.0 + 70.0 * np.cos(lat), nu)
    lat = 0.5 * math.pi * x                                 # [-pi/2, pi/2]
    eq = (2 * R_KM ** 2 / D ** 2) * 0.5 * math.pi * np.sum(w * B(lat) * np.cos(lat) ** 2)
    lat_p = 0.25 * math.pi * (x + 1)                        # [0, pi/2]
    pole = (2 * math.pi * R_KM ** 2 / D ** 2) * 0.25 * math.pi * np.sum(
        w * B(lat_p) * np.sin(lat_p) * np.cos(lat_p))
    for name, geom, ref in (
            ("equator (sub-observer 0, 0)", observe.explicit_geometry(0, 0, D), eq),
            ("pole (sub-observer 90 N)", observe.explicit_geometry(90, 0, D), pole)):
        errs = []
        for lv in LEVELS:
            mesh = icosphere(lv, R_KM)
            T = 60.0 + 70.0 * np.cos(face_latitudes(mesh))
            em = observe.disk_emission(mesh, T, geom, emissivity=eps,
                                       frequency_ghz=NU_GHZ, length_unit="km")
            errs.append(em.flux_density_jy * observe.JANSKY / ref - 1)
        p = orders(errs)
        print(f"     {name}: rel errors", ["%.2e" % e for e in errs],
              "orders", ["%.2f" % o for o in p])
        check(f"{name}: converges at second order to the quadrature",
              np.all(np.abs(p[1:] - 2) < 0.25) and abs(errs[-1]) < 2e-4,
              f"level 6 rel error {errs[-1]:.2e}")


def test_units():
    print("4. units by hand")
    # Hand calculation, constants typed independently to four digits:
    # h nu / k T = 6.626e-34 * 2.33e11 / (1.381e-23 * 100) = 0.11179
    # 2 h nu^3 / c^2 = 2 * 6.626e-34 * (2.33e11)^3 / (2.998e8)^2 = 1.8636e-16
    # B = 1.8636e-16 / (exp(0.11179) - 1) = 1.5750e-15 W m^-2 Hz^-1 sr^-1
    B = float(observe.planck_nu(100.0, 233e9))
    check("B_nu(100 K, 233 GHz) = 1.575e-15 W m^-2 Hz^-1 sr^-1",
          abs(B / 1.5750e-15 - 1) < 1e-3, f"{B:.5e}")
    # Omega = pi (1560.8 / 6.27e8)^2 = 1.9467e-11 sr; F = B Omega = 3.066e-26 = 3.066 Jy
    mesh = icosphere(6, R_KM)
    g = observe.explicit_geometry(0, 0, 6.27e8)
    em = observe.disk_emission(mesh, np.full(len(mesh["faces"]), 100.0), g,
                               emissivity=1.0, frequency_ghz=233.0, length_unit="km")
    check("Europa-sized 100 K sphere at 6.27e8 km: 3.066 Jy",
          abs(em.flux_density_jy / 3.066 - 1) < 2e-3, f"{em.flux_density_jy:.4f} Jy")
    em_m = observe.disk_emission(
        {"vertices": mesh["vertices"] * 1000.0, "faces": mesh["faces"]},
        np.full(len(mesh["faces"]), 100.0), g, emissivity=1.0,
        frequency_ghz=233.0, length_unit="m")
    check("same mesh in meters gives the same flux",
          abs(em_m.flux_density_jy / em.flux_density_jy - 1) < 1e-12)
    em_w = observe.disk_emission(mesh, np.full(len(mesh["faces"]), 100.0), g,
                                 emissivity=1.0,
                                 wavelength_um=observe.SPEED_C / 233e9 * 1e6,
                                 length_unit="km")
    check("wavelength and frequency routes agree",
          abs(em_w.flux_density_jy / em.flux_density_jy - 1) < 1e-12)
    # T_Planck - T_RJ = h nu / 2k - (h nu / k)^2 / (12 T) + O(T^-3):
    # 5.591 - 0.106 = 5.485 K here.
    x = observe.PLANCK_H * 233e9 / observe.BOLTZMANN_K
    series = x / 2 - x ** 2 / (12 * em.tb_planck_k)
    check("Planck - Rayleigh-Jeans brightness temperature = h nu/2k - (h nu/k)^2/12T",
          abs((em.tb_planck_k - em.tb_rj_k) - series) < 2e-3,
          f"{em.tb_planck_k - em.tb_rj_k:.4f} K vs {series:.4f} K")


def test_guards():
    print("6. meshes that can hide faces are refused")
    mesh = icosphere(3, R_KM)
    g = observe.explicit_geometry(0, 0, 1e6)
    T = np.full(len(mesh["faces"]), 100.0)

    def refused(m, t):
        try:
            observe.disk_emission(m, t, g, emissivity=1.0, frequency_ghz=233,
                                  length_unit="km")
        except observe.NonConvexMeshError as exc:
            return True, str(exc)
        return False, "accepted"
    ok, why = refused({"vertices": mesh["vertices"], "faces": mesh["faces"][1:]}, T[1:])
    check("open mesh refused", ok, why[:60])
    dent = mesh["vertices"].copy()
    dent[0] *= 0.97
    ok, why = refused({"vertices": dent, "faces": mesh["faces"]}, T)
    check("dented mesh refused", ok, why[:60])
    two = {"vertices": np.vstack([mesh["vertices"], mesh["vertices"] * 0.3
                                  + [3 * R_KM, 0, 0]]),
           "faces": np.vstack([mesh["faces"], mesh["faces"] + len(mesh["vertices"])])}
    ok, why = refused(two, np.concatenate([T, T]))
    check("two-piece mesh refused", ok, why[:60])
    inward = {"vertices": mesh["vertices"], "faces": mesh["faces"][:, ::-1]}
    em = observe.disk_emission(inward, T, g, emissivity=1.0, frequency_ghz=233,
                               length_unit="km")
    ref = observe.disk_emission(mesh, T, g, emissivity=1.0, frequency_ghz=233,
                                length_unit="km")
    check("inward winding is detected and corrected",
          em.flux_density_jy == ref.flux_density_jy)


# ---------------------------------------------------------------------------

def kernels_present():
    k = REPO / "data" / "kernels"
    need = ["lsk/naif0012.tls", "pck/pck00011.tpc", "spk/de440s.bsp", "spk/jup365.bsp"]
    return all((k / n).exists() for n in need), k, need


def test_spice(run_dir):
    print("3. SPICE geometry")
    ok, kdir, need = kernels_present()
    if not ok:
        skip("sub-Earth point on Europa", "Europa kernels absent; "
             "python3 scripts/fetch_kernels.py --bodies europa")
        return
    import spiceypy as sp
    files = [str(kdir / n) for n in need]
    for f in files:
        sp.furnsh(f)
    try:
        et_rx = sp.str2et("2025-01-01T00:00:00")
        spoint, trgepc, srfvec = sp.subpnt("INTERCEPT/ELLIPSOID", "EUROPA", et_rx,
                                           "IAU_EUROPA", "CN+S", "EARTH")
        utc_emit = sp.et2utc(trgepc, "ISOC", 3)
    finally:
        for f in files:
            sp.unload(f)
    _, lon_d, lat_d = sp.reclat(spoint)
    with tempfile.TemporaryDirectory() as tmp:
        run = Path(tmp)
        (run / "metadata.txt").write_text("faces=1\ntime_steps=3\n")
        (run / "provenance.txt").write_text("body=Europa\ngeneric_body=false\n")
        rows = [(trgepc + d, utc_emit) for d in (-1600.0, 0.0, 1600.0)]
        (run / "time_steps.txt").write_text(
            "# ET UTC\n" + "".join(f"{e:.6f} {u}\n" for e, u in rows))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            g = observe.observer_geometry(Results(run), "2025-01-01T00:00:00")
            gn = observe.observer_geometry(Results(run), "2025-01-01T00:00:00",
                                           time_match="nearest")
    weight_on_1 = dict(zip(g.steps, g.weights)).get(1, 0.0)
    check("emission epoch on a step: all the time weight on that step",
          abs(weight_on_1 - 1.0) < 1e-9, f"steps {g.steps} weights {g.weights}")
    check("nearest-step mode needs no shift there", gn.step == 1 and
          abs(gn.step_offset_s) < 1e-3, f"offset {gn.step_offset_s:.2e} s")
    check("sub-observer point equals a direct subpnt call",
          abs(g.sub_observer_lat - math.degrees(lat_d)) < 1e-9
          and abs(g.sub_observer_lon - math.degrees(lon_d)) < 1e-9)
    check("matches the design-time values 2.7552 N, 77.6484 E, light time 2091.64 s",
          round(g.sub_observer_lat, 4) == 2.7552
          and round(g.sub_observer_lon, 4) == 77.6484
          and round(g.light_time_s, 2) == 2091.64,
          f"{g.sub_observer_lat:.4f} N {g.sub_observer_lon:.4f} E "
          f"{g.light_time_s:.2f} s")

    test_occultation_second_date(files)

    if run_dir is None:
        skip("run's Sun direction vs SPICE", "no Europa run")
        return
    run = Results(run_dir)
    worst = 0.0
    for utc in ("2025-01-01T12:00:00", "2025-01-02T12:00:00", "2025-01-04T06:00:00"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            g = observe.observer_geometry(run, utc)
        worst = max(worst, g.sun_direction_check_deg)
    check("the run's Sun direction agrees with SPICE in IAU_EUROPA",
          worst < 1e-3, f"worst {worst:.1e} deg")


def _occultation_by_angles(sp, et_rx):
    """Classify Europa against Jupiter from Earth with no call to occult:
    light-time-corrected vectors from spkpos and the bodies' radii from
    bodvrd. "clear" and "total" are decided with the planet's largest and
    smallest radius, so they hold for the ellipsoid. Anything between is
    "edge"."""
    e, _ = sp.spkpos("EUROPA", et_rx, "J2000", "CN", "EARTH")
    j, _ = sp.spkpos("JUPITER", et_rx, "J2000", "CN", "EARTH")
    de, dj = float(np.linalg.norm(e)), float(np.linalg.norm(j))
    r_j = sp.bodvrd("JUPITER", "RADII", 3)[1]
    r_e = max(sp.bodvrd("EUROPA", "RADII", 3)[1])
    sep = math.acos(float(np.dot(e, j)) / (de * dj))
    a_e = math.asin(r_e / de)
    if de < dj:
        return "in front"
    if sep + a_e < math.asin(min(r_j) / dj):
        return "total"
    if sep - a_e > math.asin(max(r_j) / dj):
        return "clear"
    return "edge"


def test_occultation_second_date(files):
    """Jupiter hides Europa from Earth from about 19:31 to 22:11 UTC
    (received) on 2025-01-08, totally from 19:42 to 22:01 by the angular
    classification below."""
    import spiceypy as sp
    expect = {"clear": "visible", "total": "totally occulted",
              "edge": "partially occulted"}
    epochs = ["2025-01-08T19:00:00", "2025-01-08T19:36:00", "2025-01-08T20:50:00",
              "2025-01-08T22:06:00", "2025-01-08T22:40:00"]
    for f in files:
        sp.furnsh(f)
    try:
        et0 = sp.str2et("2025-01-08T15:00:00")
        classes = [_occultation_by_angles(sp, sp.str2et(t)) for t in epochs]
    finally:
        for f in files:
            sp.unload(f)
    with tempfile.TemporaryDirectory() as tmp:
        run = Path(tmp)
        n = 25
        (run / "metadata.txt").write_text(f"faces=1\ntime_steps={n}\n")
        (run / "provenance.txt").write_text("body=Europa\ngeneric_body=false\n")
        (run / "time_steps.txt").write_text("# ET UTC\n" + "".join(
            f"{et0 + 1600.0 * i:.6f} step{i}\n" for i in range(n)))
        flags = []
        for t in epochs:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                flags.append(observe.observer_geometry(Results(run), t).occultation)
    check("the angular classification spans clear, edge and total on 2025-01-08",
          classes == ["clear", "edge", "total", "edge", "clear"], str(classes))
    check("occultation flag agrees with the angular classification, 2025-01-08",
          flags == [expect[c] for c in classes],
          ", ".join(f"{t[11:16]} {f}" for t, f in zip(epochs, flags)))


def test_time(run_dir):
    print("5. matching a time to the output steps")
    if not kernels_present()[0] or run_dir is None:
        skip("time matching", "needs the Europa kernels and a Europa run")
        return
    run = Results(run_dir)
    import spiceypy as sp
    _, kdir, need = kernels_present()
    files = [str(kdir / n) for n in need]
    for f in files:
        sp.furnsh(f)
    try:
        _, et_emit, _ = sp.subpnt("INTERCEPT/ELLIPSOID", "EUROPA",
                                  sp.str2et("2025-01-01T12:00:00"),
                                  "IAU_EUROPA", "CN+S", "EARTH")
    finally:
        for f in files:
            sp.unload(f)
    et = run.times().et

    # Default: the requested time is kept; temperatures are interpolated.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        g = observe.observer_geometry(run, "2025-01-01T12:00:00")
    k0, k1 = g.steps
    w1 = (et_emit - et[k0]) / (et[k1] - et[k0])
    check("bracketing steps and weights match the SPICE emission epoch",
          (k0, k1) == (25, 26) and abs(g.weights[1] - w1) < 1e-12
          and abs(sum(g.weights) - 1.0) < 1e-15,
          f"steps {g.steps}, weights {g.weights[0]:.4f} {g.weights[1]:.4f}")
    blend = (g.weights[0] * run.field("temperature", k0).astype(float)
             + g.weights[1] * run.field("temperature", k1).astype(float))
    check("temperature_at blends the two steps",
          np.array_equal(observe.temperature_at(run, g), blend))
    check("no warning for a visible, in-window time", not caught)

    # Explicit option: move to the nearest step and report the shift.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        gn = observe.observer_geometry(run, "2025-01-01T12:00:00",
                                       time_match="nearest")
    rate = 360.0 / (3.551181 * 86400.0)           # Europa, deg per second
    check("nearest step and its offset are reported",
          gn.step == 26 and abs(gn.step_offset_s - 430.1) < 1.0,
          f"step {gn.step}, {gn.step_offset_s:.1f} s")
    check("rotation offset = Europa's rotation rate x offset",
          abs(gn.rotation_offset_deg / (abs(gn.step_offset_s) * rate) - 1) < 1e-3,
          f"{gn.rotation_offset_deg:.4f} deg")
    check("no warning below the rotation threshold", not caught)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        observe.observer_geometry(run, "2025-01-01T12:00:00",
                                  time_match="nearest", warn_rotation_deg=0.1)
    check("warning above the rotation threshold", len(caught) == 1)
    again = observe.observer_geometry(run, gn.utc_observed, time_match="nearest")
    check("the reported observation time falls exactly on the step",
          again.step == gn.step and abs(again.step_offset_s) < 0.01
          and abs(again.sub_observer_lon - gn.sub_observer_lon) < 1e-6,
          f"{again.step_offset_s:.4f} s")
    # Jupiter hides Europa from Earth near 19:10 UTC on 2025-01-01 (receive
    # time) and not at 20:40; checked independently with SPICE occult.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        hidden = observe.observer_geometry(run, "2025-01-01T19:10:00")
    check("an occulted epoch is flagged and warned",
          hidden.occulter == "JUPITER" and hidden.occultation == "totally occulted"
          and any("occulted" in str(w.message) for w in caught),
          f"{hidden.occultation}")
    seen = observe.observer_geometry(run, "2025-01-01T20:40:00")
    check("a visible epoch is not flagged", seen.occultation == "visible",
          f"{seen.occultation}")
    for utc in ("2025-01-01T00:00:00", "2025-01-10T00:00:00"):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                observe.observer_geometry(run, utc)
            check(f"{utc} outside the run is refused", False)
        except observe.ObservationTimeError:
            check(f"{utc} outside the run is refused", True)


def main():
    run_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        REPO / "output" / "quickstart_europa_diurnal")
    if not (run_dir / "temperature").is_dir():
        run_dir = None
    test_isothermal()
    test_latitude_law()
    test_units()
    test_guards()
    test_spice(run_dir)
    test_time(run_dir)
    print()
    if FAILS:
        print(f"FAILED: {len(FAILS)}: {FAILS}")
        return 1
    if SKIPS:
        print(f"PASS with {len(SKIPS)} SKIPPED: {SKIPS}")
        return 2
    print("PASS: all checks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
