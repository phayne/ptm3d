#!/usr/bin/env python3
"""Checks of ptm.observe's subsurface mode against closed forms.

    python3 src/python/tests/test_observe_subsurface.py

No engine run and no SPICE kernels: synthetic icospheres and columns.

1. Fresnel transmissivity: nadir 1 - ((n-1)/(n+1))^2; at Brewster's angle
   R_p = 0; at 50 deg the independent sine/tangent forms of R_s and R_p.
2. One column at an emission angle, in the Rayleigh-Jeans limit (1 kHz,
   where B_nu is linear in T to 1e-9), with constant properties:
   isothermal T_b = t T; linear profile on a finite column
   t (T0 + g L mu_t (1 - exp(-Z/(L mu_t)))); periodic conduction profile with
   beta = L mu_t / l (amplitude 1/sqrt((1+b)^2+b^2), lag atan(b/(1+b)));
   a 1 um skin depth recovers t T(0) on an exponential profile.
3. Disk: an isothermal sphere seen from far away has flux
   e_bar pi R^2 B(T) / D^2 with e_bar = 2 int_0^1 t(mu) mu dmu, converging
   at second order with the mesh; an isothermal column gives exactly the
   surface mode with emissivity t(theta) face by face; a vanishing skin
   depth gives the surface mode with that emissivity and T(0).
4. Per site: the disk path's intensity at one face equals the per-site script
   (brightness_temperature, temperature-dependent porous ice, Fresnel) at
   that face and angle.
5. Refusals: an emissivity in the subsurface mode; a column not starting at
   the surface.
"""

import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "src" / "python"))
sys.path.insert(0, str(REPO / "src" / "python" / "tests"))
sys.path.insert(0, str(REPO / "scripts" / "emission"))

from ptm import observe  # noqa: E402
from test_observe import icosphere, orders  # noqa: E402
import depth_weighted_emission as dwe  # noqa: E402
import dielectric_mixing as dmix  # noqa: E402

FAILS = []
R_KM = 1560.8
K, C = observe.BOLTZMANN_K, observe.SPEED_C


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def rj(intensity, nu):
    return C ** 2 * intensity / (2 * K * nu ** 2)


def test_fresnel():
    print("1. Fresnel transmissivity")
    n = 1.374
    t0 = float(dmix.fresnel_transmissivity(n, 0.0))
    check("nadir t = 1 - ((n-1)/(n+1))^2", abs(t0 - (1 - ((n - 1) / (n + 1)) ** 2)) < 1e-15,
          f"{t0:.6f}")
    for angle in (math.degrees(math.atan(n)), 50.0):
        ti = math.radians(angle)
        tt = math.asin(math.sin(ti) / n)
        rs = (math.sin(ti - tt) / math.sin(ti + tt)) ** 2
        rp = (math.tan(ti - tt) / math.tan(ti + tt)) ** 2
        t = float(dmix.fresnel_transmissivity(n, angle))
        check(f"{angle:.2f} deg: t = 1 - (R_s + R_p)/2 by the sine/tangent forms"
              + (" (Brewster, R_p = 0)" if angle != 50.0 else ""),
              abs(t - (1 - 0.5 * (rs + rp))) < 1e-14 and (angle == 50.0 or rp < 1e-30),
              f"t {t:.6f}, R_p {rp:.1e}")


def test_column():
    print("2. one column at 50 deg, Rayleigh-Jeans limit, constant properties")
    nu = 1.0e3
    n, theta = 1.374, 50.0
    mu_t = math.sqrt(1 - math.sin(math.radians(theta)) ** 2 / n ** 2)
    t = float(dmix.fresnel_transmissivity(n, theta))
    L = 0.068

    def tb(z, T, L_e):
        I, _, _ = dwe.emergent_intensity(z, T, nu, 1.0 / L_e, n, theta, t,
                                         warn_shallow=False)
        return rj(I, nu)
    z = np.concatenate([[0.0], np.geomspace(1e-5, 2.0, 400)])
    got = tb(z, np.full_like(z, 100.0), L)
    check("isothermal: T_b = t T", abs(got / (t * 100.0) - 1) < 1e-9, f"{got:.9f}")
    for Z in (0.1, 2.0):
        zz = np.concatenate([[0.0], np.geomspace(1e-5, Z, 300)])
        got = tb(zz, 90.0 + 40.0 * zz, L)
        want = t * (90.0 + 40.0 * L * mu_t * -math.expm1(-Z / (L * mu_t)))
        check(f"linear profile, column to {Z} m", abs(got / want - 1) < 1e-9,
              f"{got:.6f} vs {want:.6f}")
    ell, A = 0.035, 20.0
    zz = np.linspace(0.0, 1.5, 30001)
    phase = np.linspace(0, 2 * np.pi, 73)[:-1]
    resp = np.array([tb(zz, 100.0 + A * np.exp(-zz / ell) * np.cos(p - zz / ell), L)
                     for p in phase]) - t * 100.0
    c = np.sum(resp * np.exp(-1j * phase)) * 2 / len(phase)
    b = L * mu_t / ell
    amp, lag = abs(c) / (t * A), -np.angle(c)
    check(f"periodic profile, beta = L mu_t / l = {b:.3f}: amplitude and lag",
          abs(amp - 1 / math.hypot(1 + b, b)) < 1e-6
          and abs(lag - math.atan(b / (1 + b))) < 1e-6,
          f"amp {amp:.6f} vs {1 / math.hypot(1 + b, b):.6f}, lag {math.degrees(lag):.3f} deg")
    zz = np.concatenate([[0.0], np.geomspace(1e-10, 1.0, 400)])
    Ls, d = 1e-6, 0.05
    ks = 1.0 / (Ls * mu_t)
    got = tb(zz, 100.0 + 30.0 * np.exp(-zz / d), Ls)
    want = t * (100.0 + 30.0 * ks * d / (ks * d + 1))
    check("1 um skin depth recovers t T(0) (exact exponential form)",
          abs(got - want) < 1e-6, f"{got:.7f} vs {want:.7f} (t T(0) = {t * 130:.4f})")


def test_disk():
    print("3. disk integrals")
    n, L, T0 = 1.374, 0.068, 100.0
    sub = observe.Subsurface(skin_depth_m=L, refractive_index=n)
    D = 1.0e8 * R_KM
    g = observe.explicit_geometry(10.0, 20.0, D)
    nu = 233e9
    x, w = np.polynomial.legendre.leggauss(400)
    mu = 0.5 * (x + 1)
    theta = np.degrees(np.arccos(mu))
    e_bar = 2 * 0.5 * np.sum(w * dmix.fresnel_transmissivity(n, theta) * mu)
    ref = e_bar * math.pi * R_KM ** 2 / D ** 2 * observe.planck_nu(T0, nu) / observe.JANSKY
    errs = []
    z = np.array([0.0, 0.3, 1.0])
    for lv in (3, 4, 5, 6):
        mesh = icosphere(lv, R_KM)
        nf = len(mesh["faces"])
        em = observe.disk_emission(mesh, np.full((nf, 3), T0), g, frequency_ghz=233.0,
                                   subsurface=sub, depths_m=z, length_unit="km")
        errs.append(em.flux_density_jy / ref - 1)
    p = orders(errs)
    print(f"     disk emissivity e_bar = {e_bar:.6f}; rel errors",
          ["%.2e" % e for e in errs], "orders", ["%.2f" % o for o in p])
    check("isothermal sphere with Fresnel t: flux -> e_bar pi R^2 B / D^2, "
          "second order", np.all(np.abs(p - 2) < 0.25) and abs(errs[-1]) < 2e-4,
          f"ico6 rel error {errs[-1]:.1e}")

    mesh = icosphere(4, R_KM)
    nf = len(mesh["faces"])
    rng = np.random.default_rng(1)
    Ts = 80.0 + 40.0 * rng.random(nf)
    sub_em = observe.disk_emission(mesh, np.repeat(Ts[:, None], 3, axis=1), g,
                                   frequency_ghz=233.0, subsurface=sub, depths_m=z,
                                   length_unit="km")
    dm = observe.disk_mesh(mesh, length_unit="km")
    P = g.observer_km * 1000.0
    r = P - dm.centers_m
    mu_f = np.einsum("ij,ij->i", dm.normals, r) / np.linalg.norm(r, axis=1)
    t_face = np.where(mu_f > 0, dmix.fresnel_transmissivity(
        n, np.degrees(np.arccos(np.clip(mu_f, -1, 1)))), 1.0)
    surf_em = observe.disk_emission(mesh, Ts, g, frequency_ghz=233.0,
                                    emissivity=t_face, length_unit="km")
    check("isothermal columns = surface mode with emissivity t(theta), face by face",
          abs(sub_em.flux_density_jy / surf_em.flux_density_jy - 1) < 1e-12)

    zz = np.concatenate([[0.0], np.geomspace(1e-10, 1.0, 200)])
    cols = Ts[:, None] - 30.0 * (1 - np.exp(-zz[None, :] / 0.05))
    thin = observe.disk_emission(mesh, cols, g, frequency_ghz=233.0,
                                 subsurface=observe.Subsurface(skin_depth_m=1e-8,
                                                               refractive_index=n),
                                 depths_m=zz, length_unit="km")
    check("vanishing skin depth = surface mode with t(theta) and T(0)",
          abs(thin.flux_density_jy / surf_em.flux_density_jy - 1) < 1e-6,
          f"rel {thin.flux_density_jy / surf_em.flux_density_jy - 1:.1e}")
    check("beta = L_e / delta_thermal reported",
          abs(observe.disk_emission(mesh, np.repeat(Ts[:, None], 3, axis=1), g,
                                    frequency_ghz=233.0, subsurface=sub, depths_m=z,
                                    length_unit="km",
                                    thermal_skin_depth_m=0.034).beta_nadir
              - L / 0.034) < 1e-12)


def test_per_site():
    print("4. one face of the disk = the per-site script")
    mesh = icosphere(3, R_KM)
    nf = len(mesh["faces"])
    g = observe.explicit_geometry(5.0, 40.0, 6.3e8)
    z = np.concatenate([[0.0], np.geomspace(2e-3, 1.0, 30)])
    lon = np.linspace(0, 2 * np.pi, nf)
    cols = 95.0 + 25.0 * np.cos(lon)[:, None] * np.exp(-z[None, :] / 0.04)
    sub = observe.Subsurface(porosity=0.5, reference_temperature_k=None)
    em = observe.disk_emission(mesh, cols, g, frequency_ghz=233.0, subsurface=sub,
                               depths_m=z, length_unit="km")
    dm = observe.disk_mesh(mesh, length_unit="km")
    P = g.observer_km * 1000.0
    worst = 0.0
    for i in np.nonzero(em.face_intensity > 0)[0][::37]:
        r = P - dm.centers_m[i]
        theta = math.degrees(math.acos(np.dot(dm.normals[i], r) / np.linalg.norm(r)))
        site = dwe.brightness_temperature(z, cols[i], 233.0, 0.5, theta, "fresnel")
        worst = max(worst, abs(em.face_intensity[i] / site.intensity - 1))
    check("face intensity equals brightness_temperature(...).intensity",
          worst < 1e-12, f"worst rel {worst:.1e}")


def test_refusals():
    print("5. refusals")
    mesh = icosphere(2, R_KM)
    nf = len(mesh["faces"])
    g = observe.explicit_geometry(0, 0, 1e6)
    sub = observe.Subsurface(skin_depth_m=0.05, refractive_index=1.37)
    for name, kwargs in (
            ("emissivity in the subsurface mode",
             dict(subsurface=sub, depths_m=[0.0, 1.0], emissivity=0.85)),
            ("column not starting at the surface",
             dict(subsurface=sub, depths_m=[0.01, 1.0]))):
        try:
            observe.disk_emission(mesh, np.full((nf, 2), 100.0), g,
                                  frequency_ghz=233.0, length_unit="km", **kwargs)
            check(f"refused: {name}", False)
        except ValueError:
            check(f"refused: {name}", True)


def main():
    test_fresnel()
    test_column()
    test_disk()
    test_per_site()
    test_refusals()
    print()
    if FAILS:
        print(f"FAILED: {FAILS}")
        return 1
    print("PASS: all checks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
