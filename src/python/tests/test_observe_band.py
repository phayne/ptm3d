#!/usr/bin/env python3
"""Checks of the thermal-infrared band mode of ptm.observe.

    python3 src/python/tests/test_observe_band.py

Every check is against a closed form or an independent quadrature, never
against the module itself:

1. band radiance: the Gauss-Legendre quadrature over a top-hat band agrees
   with the closed-form series for the integral of the Planck function, and
   a tabulated triangular response agrees with a fine trapezoid sum, both to
   better than 1e-6 relative from 20 K to the temperature of the Sun;
2. band brightness temperature: inverting the band radiance returns the
   temperature, and a narrow band reproduces the monochromatic Planck
   brightness temperature, converging at second order in the band width;
3. isothermal sphere: the band brightness temperature of an emissivity-1
   icosphere equals its temperature at every resolution, and its flux
   converges at second order to e pi (R/D)^2 times the band radiance;
4. sphere in instantaneous radiative equilibrium at zero phase
   (T = T_ss cos(i)^(1/4) on the day side): the thermal and reflected
   fluxes converge with mesh resolution to one-dimensional integrals over
   the emission cosine (thermal: Gauss-Legendre; reflected: the Lambert
   sphere's 2/3);
5. the existing single-frequency calculation is unchanged: its face fluxes
   equal, bit for bit, those of the original expressions written out here;
6. inputs that cannot be used are refused.

No SPICE kernels or model runs are needed.
"""

import math
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "src" / "python"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from ptm import observe  # noqa: E402
from test_observe import icosphere, orders  # noqa: E402

H, K, C = observe.PLANCK_H, observe.BOLTZMANN_K, observe.SPEED_C
FAILS = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def series_tail(x, terms=4000):
    """Integral from x to infinity of u^3 / (e^u - 1) du, by the series
    sum_n e^(-n x) (x^3/n + 3 x^2/n^2 + 6 x/n^3 + 6/n^4)."""
    n = np.arange(1, terms + 1, dtype=np.float64)
    return float(np.sum(np.exp(-n * x) * (x ** 3 / n + 3 * x ** 2 / n ** 2
                                          + 6 * x / n ** 3 + 6 / n ** 4)))


def series_band(T, lo_um, hi_um):
    """Integral of B_lam(T) from lo to hi micrometers, W m^-2 sr^-1."""
    x_lo = H * C / (lo_um * 1e-6 * K * T)       # the larger x
    x_hi = H * C / (hi_um * 1e-6 * K * T)
    return 2 * K ** 4 * T ** 4 / (H ** 3 * C ** 2) * (series_tail(x_hi) - series_tail(x_lo))


TEMPS = [20.0, 50.0, 80.0, 100.0, 130.0, 200.0, 300.0, 1000.0, 5772.0]


def test_band_radiance():
    print("1. band radiance against closed forms")
    worst = 0.0
    for lo, hi in ((8.0, 14.0), (10.5, 11.5), (3.0, 50.0)):
        band = observe.Band.top_hat(lo, hi)
        for T in TEMPS:
            ref = series_band(T, lo, hi)
            worst = max(worst, abs(band.radiance(T) / ref - 1))
    check("top hat vs series, 3 bands x 20-5772 K", worst < 1e-6, f"max rel {worst:.1e}")

    # Stefan-Boltzmann: a band covering essentially all of the emission.
    T = 300.0
    full = observe.Band.top_hat(0.5, 1.0e5).radiance(T) * math.pi
    sb = observe.STEFAN_BOLTZMANN * T ** 4
    check("0.5 um to 10 cm at 300 K = sigma T^4 / pi", abs(full / sb - 1) < 1e-6,
          f"rel {full / sb - 1:.1e}")
    check("5772 K blackbody Sun gives 1361 W m^-2 at 1 au",
          abs(observe.STEFAN_BOLTZMANN * 5772.0 ** 4 * (6.957e8 / 1.495978707e11) ** 2
              - 1361.0) < 1.0)

    # A triangular response in a file, against a fine trapezoid sum.
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "triangle.txt"
        path.write_text("# wavelength_um, response\n7.0, 0.0\n9.0, 1.0\n"
                        "12.0, 0.5\n14.0, 0.0\n")
        band = observe.Band.from_file(path)
    lam = np.linspace(7.0, 14.0, 700001)
    resp = np.interp(lam, [7, 9, 12, 14], [0, 1, 0.5, 0])
    worst = 0.0
    for T in TEMPS[:-1]:
        y = resp * observe.planck_lambda(T, lam * 1e-6)
        ref = float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(lam * 1e-6)))
        worst = max(worst, abs(band.radiance(T) / ref - 1))
    check("tabulated response vs fine trapezoid", worst < 1e-6, f"max rel {worst:.1e}")


def test_brightness_temperature():
    print("2. band brightness temperature")
    band = observe.Band.top_hat(8.0, 14.0)
    T = np.array([25.0, 60.0, 90.0, 110.0, 130.0, 250.0, 400.0])
    back = band.brightness_temperature(band.radiance(T))
    err = float(np.max(np.abs(back - T)))
    check("inverse of the band radiance", err < 1e-9, f"max {err:.1e} K")
    check("zero radiance gives NaN", np.isnan(band.brightness_temperature(0.0)))

    # A narrow band around lam0: e * B_lam(T) over the band has the Planck
    # brightness temperature of e * B_nu(T) at c / lam0, up to O(width^2).
    lam0, T0, e = 10.0, 110.0, 0.9
    nu = C / (lam0 * 1e-6)
    mono = float(observe.planck_nu_temperature(e * observe.planck_nu(T0, nu), nu))
    widths = [0.8, 0.4, 0.2, 0.1, 0.05]
    errs = []
    for w in widths:
        b = observe.Band.top_hat(lam0 - w / 2, lam0 + w / 2)
        errs.append(b.brightness_temperature(e * b.radiance(T0)) - mono)
    p = orders(errs)
    print("     widths (um)", widths, "error (K)", ["%.2e" % x for x in errs],
          "orders", ["%.2f" % o for o in p])
    check("narrow band -> monochromatic Planck brightness temperature",
          abs(errs[-1]) < 1e-4 and np.all(np.abs(p - 2) < 0.05),
          f"{abs(errs[-1]):.1e} K at 0.05 um, order {p[-1]:.2f}")


R_KM, D_KM = 1560.8, 6.3e8
LEVELS = [2, 3, 4, 5, 6]


def test_isothermal():
    print("3. isothermal sphere")
    band = observe.Band.top_hat(8.0, 14.0)
    T, e = 120.0, 0.9
    geom = observe.explicit_geometry(10.0, 33.0, D_KM)
    exact = e * math.pi * (R_KM / D_KM) ** 2 * band.radiance(T)
    tb_err, f_err = [], []
    for lv in LEVELS:
        mesh = icosphere(lv, R_KM)
        n = len(mesh["faces"])
        em1 = observe.band_emission(mesh, np.full(n, T), geom, band=band,
                                    emissivity=1.0, reflected=False, length_unit="km")
        tb_err.append(abs(em1.tb_k - T))
        em = observe.band_emission(mesh, np.full(n, T), geom, band=band,
                                   emissivity=e, reflected=False, length_unit="km")
        f_err.append(em.flux_w_m2 / exact - 1)
    p = orders(f_err)
    print("     levels", LEVELS, "flux rel. error", ["%.2e" % x for x in f_err],
          "orders", ["%.2f" % o for o in p])
    check("emissivity 1: band brightness temperature = T at every level",
          max(tb_err) < 1e-9, f"max {max(tb_err):.1e} K")
    check("flux -> e pi (R/D)^2 integral R B dlam at second order",
          abs(f_err[-1]) < 1e-4 and np.all(np.abs(p - 2) < 0.1),
          f"{f_err[-1]:.1e} at level {LEVELS[-1]}, order {p[-1]:.2f}")


def test_radiative_equilibrium():
    print("4. sphere in radiative equilibrium at zero phase")
    band = observe.Band.top_hat(8.0, 14.0)
    e, S, T_ss = 0.9, 50.0, 130.0
    # Observer and Sun both on the +x axis (zero phase).
    geom = observe.explicit_geometry(0.0, 0.0, D_KM)
    # One-dimensional references. Thermal: 2 pi R^2/D^2 * integral over the
    # emission cosine mu of e L(T_ss mu^(1/4)) mu dmu, by Gauss-Legendre.
    # Spectral integral by the closed-form series, independent of Band.
    x, w = np.polynomial.legendre.leggauss(200)
    mu = 0.5 * (x + 1)
    L = np.array([series_band(T_ss * m ** 0.25, 8.0, 14.0) for m in mu])
    omega_scale = 2 * math.pi * (R_KM / D_KM) ** 2
    F_th_ref = omega_scale * e * float(np.sum(0.5 * w * L * mu))
    # Reflected: a Lambert sphere at zero phase reflects (2/3) of a disk of
    # the same reflectance facing the Sun.
    f_sun = math.pi * series_band(observe.SUN_TEMPERATURE_K, 8.0, 14.0) / (
        observe.STEFAN_BOLTZMANN * observe.SUN_TEMPERATURE_K ** 4)
    F_r_ref = (2.0 / 3.0) * (1 - e) * f_sun * S * (R_KM / D_KM) ** 2
    th_err, r_err = [], []
    for lv in LEVELS:
        mesh = icosphere(lv, R_KM)
        c = mesh["vertices"][mesh["faces"]].mean(axis=1)
        v0, v1, v2 = (mesh["vertices"][mesh["faces"][:, k]] for k in range(3))
        nrm = np.cross(v1 - v0, v2 - v0)
        nrm /= np.linalg.norm(nrm, axis=1)[:, None]
        nrm *= np.sign(np.einsum("ij,ij->i", nrm, c))[:, None]
        cos_i = nrm[:, 0]
        T = np.where(cos_i > 0, T_ss * np.clip(cos_i, 0, None) ** 0.25, 1.0)
        T = np.maximum(T, 1.0)
        em = observe.band_emission(mesh, T, geom, band=band, emissivity=e,
                                   direct_solar=S * np.clip(cos_i, 0, None),
                                   length_unit="km")
        th_err.append(em.thermal_flux_w_m2 / F_th_ref - 1)
        r_err.append(em.reflected_flux_w_m2 / F_r_ref - 1)
    pt, pr = orders(th_err), orders(r_err)
    print("     levels", LEVELS, "thermal rel. error", ["%.2e" % x for x in th_err],
          "orders", ["%.2f" % o for o in pt])
    print("     levels", LEVELS, "reflected rel. error", ["%.2e" % x for x in r_err],
          "orders", ["%.2f" % o for o in pr])
    check("thermal flux converges to the 1-D integral",
          abs(th_err[-1]) < 1e-3 and np.all(np.abs(th_err[1:]) < np.abs(th_err[:-1])),
          f"{th_err[-1]:.1e} at level {LEVELS[-1]}, last order {pt[-1]:.2f}")
    check("reflected flux converges to the Lambert sphere",
          abs(r_err[-1]) < 1e-3 and np.all(np.abs(r_err[1:]) < np.abs(r_err[:-1])),
          f"{r_err[-1]:.1e} at level {LEVELS[-1]}, last order {pr[-1]:.2f}")


def test_single_frequency_unchanged():
    print("5. single-frequency surface mode unchanged")
    mesh = icosphere(4, R_KM)
    rng = np.random.default_rng(7)
    T = 80.0 + 50.0 * rng.random(len(mesh["faces"]))
    worst = 0
    for lat, lon, d in ((3.0, 27.0, D_KM), (-40.0, 200.0, 4.0 * R_KM)):
        geom = observe.explicit_geometry(lat, lon, d)
        em = observe.disk_emission(mesh, T, geom, emissivity=0.9,
                                   frequency_ghz=233.0, length_unit="km")
        # The original expressions, written out.
        dm = observe.disk_mesh(mesh, length_unit="km")
        P = np.asarray(geom.observer_km, dtype=np.float64) * 1000.0
        D = float(np.linalg.norm(P))
        r = P[None, :] - dm.centers_m
        dist = np.linalg.norm(r, axis=1)
        mu = np.einsum("ij,ij->i", dm.normals, r) / dist
        visible = mu > 0
        cos_los = (r @ P) / (dist * D)
        domega = np.where(visible, dm.areas_m2 * mu * cos_los / dist ** 2, 0.0)
        face_flux = 0.9 * observe.planck_nu(T, 233.0e9) * domega
        worst += int(np.sum(face_flux / observe.JANSKY != em.face_flux_jy))
        worst += int(float(np.sum(domega)) != em.omega_projected_sr)
    check("face fluxes and solid angle bit for bit", worst == 0,
          f"{worst} values differ")


def test_refusals():
    print("6. refusals")
    mesh = icosphere(2, R_KM)
    n = len(mesh["faces"])
    geom = observe.explicit_geometry(0.0, 0.0, D_KM)
    band = observe.Band.top_hat(8.0, 14.0)

    def refused(fn):
        try:
            fn()
        except (ValueError, TypeError):
            return True
        return False
    check("reflected sunlight without direct_solar", refused(
        lambda: observe.band_emission(mesh, np.full(n, 100.0), geom, band=band,
                                      emissivity=0.9, length_unit="km")))
    check("emissivity above 1", refused(
        lambda: observe.band_emission(mesh, np.full(n, 100.0), geom, band=band,
                                      emissivity=1.2, reflected=False,
                                      length_unit="km")))
    check("band limits reversed", refused(lambda: observe.Band.top_hat(14, 8)))
    with tempfile.TemporaryDirectory() as tmp:
        bad = Path(tmp) / "bad.txt"
        bad.write_text("10 1\n9 1\n")
        check("response wavelengths not increasing",
              refused(lambda: observe.Band.from_file(bad)))
        bad.write_text("8 0\n14 0\n")
        check("response zero everywhere", refused(lambda: observe.Band.from_file(bad)))


def main():
    test_band_radiance()
    test_brightness_temperature()
    test_isothermal()
    test_radiative_equilibrium()
    test_single_frequency_unchanged()
    test_refusals()
    if FAILS:
        print(f"FAIL: {len(FAILS)} check(s): {', '.join(FAILS)}")
        return 1
    print("PASS: all checks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
