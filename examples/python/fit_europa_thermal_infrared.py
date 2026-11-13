#!/usr/bin/env python3
"""Fit thermal inertia and emissivity to disk-integrated 8-14 um observations.

The observations here are SYNTHETIC. They are made from one PTM3D simulation
of Europa with a known thermal inertia and emissivity, plus noise, and then
fitted with a grid of other simulations, so you can see how well the two
parameters are recovered and how they trade off. No telescope data are used.

    export OMP_NUM_THREADS=8
    python3 examples/python/fit_europa_thermal_infrared.py           # about 18 minutes
    python3 examples/python/fit_europa_thermal_infrared.py --no-run  # refit finished runs

What it does
------------
1. Runs examples/quickstart_europa_diurnal.yaml (one Europa orbit starting
   2025-01-01, ico4 mesh) on a grid of thermal inertia and emissivity, with
   a thermal inertia that is constant with depth and temperature
   (``constant_thermal_inertia``) and the Bond albedo held at 0.5. Each point
   also has a five-hour run with one-minute output steps through the Jupiter
   eclipse of 2025-01-01 (as examples/quickstart_europa_eclipse.yaml).
2. Places Earth at a set of observation times (when the light arrives) and
   computes the 8-14 um flux of the whole visible disk with ``ptm.observe``,
   including the sunlight the surface reflects in the band. Half of the
   times sample Europa's rotation and use the one-orbit runs. The other half
   follow Europa out of Jupiter's shadow and use the one-minute runs.
3. Makes synthetic observations from a run at thermal inertia 70 and
   emissivity 0.92 (not grid values) with 5 percent noise, and fits thermal
   inertia and emissivity by chi-square on a fine grid, with and without the
   eclipse times. The model flux at each time is a bicubic spline in the
   logarithm of thermal inertia and in emissivity.
4. Fits a second set of synthetic observations, made with a Bond albedo of
   0.55, with the same albedo-0.5 grid, to show what a wrong assumed albedo
   does to the fitted emissivity.

The emissivity is gray: the same value sets the energy balance of the
simulation, the band emission, and the in-band reflectance 1 - e. Changing
it changes the temperatures, so every emissivity needs its own simulations.

Needs: a built engine, the Europa SPICE kernels
(``python3 scripts/fetch_kernels.py --bodies europa``), scipy and matplotlib.
The runs take about 2.1 GB of disk.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
import warnings
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src" / "python"))

from ptm import Sweep, observe  # noqa: E402

BASE = REPO / "examples" / "quickstart_europa_diurnal.yaml"
BAND = (8.0, 14.0)                 # micrometers, top hat
ALBEDO = 0.5                       # held fixed in the fit
TRUE_TI = 70.0
TRUE_EMISSIVITY = 0.92
WRONG_ALBEDO = 0.55                # for the second set of observations
NOISE = 0.05                       # fractional, per observation
GRID_TI = [20, 35, 50, 85, 120, 170, 240]        # excludes the true value
GRID_EMISSIVITY = [0.80, 0.85, 0.90, 0.95, 1.00]  # excludes the true value
SEED = 20261116

#: Observation (receive) times, UTC. Europa is visible from Earth at each.
ROTATION_TIMES = ["2025-01-01T08:00:00", "2025-01-02T04:00:00",
                  "2025-01-02T16:00:00", "2025-01-03T04:00:00",
                  "2025-01-03T16:00:00", "2025-01-04T04:00:00"]
ECLIPSE_TIMES = ["2025-01-01T20:10:00", "2025-01-01T20:30:00",
                 "2025-01-01T20:50:00", "2025-01-01T21:10:00",
                 "2025-01-01T21:40:00", "2025-01-01T22:20:00"]

#: The eclipse runs, as examples/quickstart_europa_eclipse.yaml.
ECLIPSE_RUN = {"simulation.start_time": "2025-01-01T17:45:00",
               "simulation.duration_days": 300 / 1440,
               "simulation.timesteps": 300,
               "thermal.equilibration_samples_per_diurnal_cycle": 192}


def make_sweeps(outdir):
    """The model grid, the true run, and the true run at the wrong albedo,
    for the whole orbit and for the eclipse at one-minute steps."""
    sweeps = {}
    for kind, extra in (("orbit", {}), ("eclipse", ECLIPSE_RUN)):
        fixed = dict(extra, **{"thermophysical.albedo_override": ALBEDO})
        sweeps[kind] = {
            "grid": Sweep(BASE, grid={"constant_thermal_inertia": GRID_TI,
                                      "thermophysical.emissivity_override": GRID_EMISSIVITY},
                          outdir=outdir / f"grid_{kind}", overrides=fixed),
            "truth": Sweep(BASE, grid={"constant_thermal_inertia": [TRUE_TI],
                                       "thermophysical.emissivity_override": [TRUE_EMISSIVITY]},
                           outdir=outdir / f"truth_{kind}", overrides=fixed),
            "wrong_albedo": Sweep(
                BASE, grid={"constant_thermal_inertia": [TRUE_TI],
                            "thermophysical.emissivity_override": [TRUE_EMISSIVITY]},
                outdir=outdir / f"truth_albedo{WRONG_ALBEDO:g}_{kind}",
                overrides=dict(extra, **{"thermophysical.albedo_override": WRONG_ALBEDO})),
        }
    return sweeps


def geometries(run, times):
    """Observer geometry at each receive time. A time at which Jupiter hides
    Europa is an error, never a silent data point."""
    out = []
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        for t in times:
            out.append(observe.observer_geometry(run, t, observer="EARTH"))
    return out


def band_fluxes(runs, emissivities, geoms, band):
    """In-band flux (W m^-2), thermal emission plus reflected sunlight, one
    row per run and one column per geometry, and the solid angle of the
    visible disk at each geometry (sr). Each run is observed with its own
    (gray) emissivity. All runs must share one time axis and one mesh."""
    flux = np.empty((len(runs), len(geoms)))
    omega = np.empty(len(geoms))
    for i, (run, e) in enumerate(zip(runs, emissivities)):
        for j, g in enumerate(geoms):
            em = observe.band_emission(
                run, observe.temperature_at(run, g), g, band=band, emissivity=e,
                direct_solar=observe.temperature_at(run, g, "direct_solar"))
            flux[i, j] = em.flux_w_m2
            omega[j] = em.omega_projected_sr
    return flux, omega


def chi2_surface(log_ti, eps, flux_grid, obs, sigma, ti_fine, eps_fine):
    """Chi-square on (emissivity, thermal inertia). flux_grid has shape
    (TI, emissivity, time), and each time's flux is a bicubic spline in
    (log TI, emissivity)."""
    from scipy.interpolate import RectBivariateSpline
    model = np.stack([RectBivariateSpline(log_ti, eps, flux_grid[:, :, j])(
        np.log(ti_fine), eps_fine) for j in range(flux_grid.shape[2])])  # (time, TI, eps)
    resid = (model - obs[:, None, None]) / sigma[:, None, None]
    return np.sum(resid ** 2, axis=0).T, model                 # chi2 (eps, TI)


def interval(values, chi2_profile):
    """Range where the profiled chi-square is within 1 of its minimum."""
    ok = values[chi2_profile <= chi2_profile.min() + 1.0]
    return ok.min(), ok.max()


def fit_report(name, n, c2, ti_fine, eps_fine):
    e_best, t_best = np.unravel_index(np.argmin(c2), c2.shape)
    ti_lo, ti_hi = interval(ti_fine, c2.min(axis=0))
    ep_lo, ep_hi = interval(eps_fine, c2.min(axis=1))
    print(f"{name:34s} ({n:2d} times): thermal inertia {ti_fine[t_best]:6.1f} "
          f"[{ti_lo:.1f}, {ti_hi:.1f}]   emissivity {eps_fine[e_best]:.3f} "
          f"[{ep_lo:.3f}, {ep_hi:.3f}]   min chi-square {c2.min():.1f}")
    return ti_fine[t_best], eps_fine[e_best]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--outdir", default=str(REPO / "output" / "europa_tir_fit"))
    parser.add_argument("--no-run", action="store_true",
                        help="use finished runs, start no simulations")
    args = parser.parse_args()
    outdir = Path(args.outdir)
    band = observe.Band.top_hat(*BAND)

    sweeps = make_sweeps(outdir)
    if not args.no_run:
        t0 = time.time()
        for kind in sweeps.values():
            for sw in kind.values():
                sw.run()
        print(f"simulations: {(time.time() - t0) / 60:.1f} minutes")

    times = ROTATION_TIMES + ECLIPSE_TIMES
    n_rot = len(ROTATION_TIMES)
    flux, true_flux, wrong_flux, geoms, omega = {}, {}, {}, {}, {}
    for kind, tlist in (("orbit", ROTATION_TIMES), ("eclipse", ECLIPSE_TIMES)):
        grid = sweeps[kind]["grid"]
        runs = grid.results()
        if len(runs) != len(GRID_TI) * len(GRID_EMISSIVITY):
            raise SystemExit(f"{grid.outdir}: {len(runs)} finished runs, expected "
                             f"{len(GRID_TI) * len(GRID_EMISSIVITY)}")
        ti = grid.values("constant_thermal_inertia")
        ep = grid.values("thermophysical.emissivity_override")
        # The swept parameters must actually have changed the runs.
        t_first = [r.field("temperature", 0) for r in runs]
        if any(np.array_equal(t_first[0], t) for t in t_first[1:]):
            raise SystemExit("two grid points have identical outputs, so a swept "
                             "key had no effect")
        geoms[kind] = geometries(runs[0], tlist)
        f, omega[kind] = band_fluxes(runs, ep, geoms[kind], band)
        cube = np.empty((len(GRID_TI), len(GRID_EMISSIVITY), len(tlist)))
        for row, a, b in zip(f, ti, ep):
            cube[GRID_TI.index(a), GRID_EMISSIVITY.index(b)] = row
        flux[kind] = cube
        true_flux[kind] = band_fluxes(sweeps[kind]["truth"].results(),
                                      [TRUE_EMISSIVITY], geoms[kind], band)[0][0]
        wrong_flux[kind] = band_fluxes(sweeps[kind]["wrong_albedo"].results(),
                                       [TRUE_EMISSIVITY], geoms[kind], band)[0][0]

    flux_grid = np.concatenate([flux["orbit"], flux["eclipse"]], axis=2)
    clean = np.concatenate([true_flux["orbit"], true_flux["eclipse"]])
    clean_wrong = np.concatenate([wrong_flux["orbit"], wrong_flux["eclipse"]])
    omega_all = np.concatenate([omega["orbit"], omega["eclipse"]])
    geoms_all = geoms["orbit"] + geoms["eclipse"]

    rng = np.random.default_rng(SEED)
    sigma = NOISE * clean
    noise = rng.standard_normal(len(times))
    obs = clean + sigma * noise
    obs_wrong = clean_wrong + NOISE * clean_wrong * noise

    def tb(f):
        return band.brightness_temperature(f / omega_all)
    tb_obs = tb(obs)
    tb_err = 0.5 * (tb(obs + sigma) - tb(obs - sigma))

    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "synthetic_observations.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([f"# SYNTHETIC: {BAND[0]:g}-{BAND[1]:g} um, constant thermal "
                    f"inertia {TRUE_TI:g}, emissivity {TRUE_EMISSIVITY}, Bond albedo "
                    f"{ALBEDO}, {100 * NOISE:.0f} percent noise, seed {SEED}"])
        w.writerow(["utc_observed", "set", "sub_observer_lon_deg", "phase_angle_deg",
                    "flux_w_m2", "sigma_w_m2", "tb_k", "sigma_tb_k"])
        for k, g in enumerate(geoms_all):
            w.writerow([g.utc_observed, "rotation" if k < n_rot else "eclipse",
                        f"{g.sub_observer_lon:.3f}", f"{g.phase_angle_deg:.3f}",
                        f"{obs[k]:.6e}", f"{sigma[k]:.6e}", f"{tb_obs[k]:.3f}",
                        f"{tb_err[k]:.3f}"])

    log_ti = np.log(np.array(GRID_TI, dtype=float))
    eps = np.array(GRID_EMISSIVITY)
    ti_fine = np.geomspace(GRID_TI[0], GRID_TI[-1], 400)
    eps_fine = np.linspace(GRID_EMISSIVITY[0], GRID_EMISSIVITY[-1], 401)

    # How well the spline through the grid reproduces the true run, which
    # is not a grid point.
    _, m_true = chi2_surface(log_ti, eps, flux_grid, obs, sigma,
                             np.array([TRUE_TI]), np.array([TRUE_EMISSIVITY]))
    interp_err = np.max(np.abs(m_true[:, 0, 0] / clean - 1))
    print(f"\nSYNTHETIC observations at {BAND[0]:g}-{BAND[1]:g} um: true thermal "
          f"inertia {TRUE_TI:g}, true emissivity {TRUE_EMISSIVITY}, Bond albedo "
          f"{ALBEDO}, {100 * NOISE:.0f} percent noise")
    print(f"The model grid reproduces the true run to {100 * interp_err:.2f} "
          f"percent in flux (the noise is {100 * NOISE:.0f} percent).\n")

    sets = {"rotation only": np.arange(n_rot),
            "rotation and eclipse": np.arange(len(times))}
    surfaces, best = {}, {}
    for name, idx in sets.items():
        c2, model = chi2_surface(log_ti, eps, flux_grid[:, :, idx], obs[idx],
                                 sigma[idx], ti_fine, eps_fine)
        surfaces[name] = c2
        best[name] = fit_report(name, len(idx), c2, ti_fine, eps_fine)
    print("\nWithout noise (expected widths):")
    for name, idx in sets.items():
        c2, _ = chi2_surface(log_ti, eps, flux_grid[:, :, idx], clean[idx],
                             sigma[idx], ti_fine, eps_fine)
        fit_report(name, len(idx), c2, ti_fine, eps_fine)
    print(f"\nObservations made with Bond albedo {WRONG_ALBEDO}, fitted assuming {ALBEDO}:")
    c2w, _ = chi2_surface(log_ti, eps, flux_grid, obs_wrong, NOISE * clean_wrong,
                          ti_fine, eps_fine)
    surfaces["wrong albedo"] = c2w
    fit_report("rotation and eclipse", len(times), c2w, ti_fine, eps_fine)
    print("\nIntervals: profiled chi-square within 1 of its minimum.")

    # Sensitivity of the rotation-time fluxes at the true values: thermal
    # inertia and emissivity from the grid, emissivity at fixed temperatures
    # from the true run, albedo from the run at the second albedo.
    rot = np.arange(n_rot)

    def model_at(ti_v, ep_v):
        _, m = chi2_surface(log_ti, eps, flux_grid[:, :, rot], obs[rot], sigma[rot],
                            np.array([ti_v]), np.array([ep_v]))
        return m[:, 0, 0]
    d_ti = np.log(model_at(TRUE_TI * 1.01, TRUE_EMISSIVITY)
                  / model_at(TRUE_TI / 1.01, TRUE_EMISSIVITY)) / (2 * np.log(1.01))
    d_ep = np.log(model_at(TRUE_TI, TRUE_EMISSIVITY * 1.01)
                  / model_at(TRUE_TI, TRUE_EMISSIVITY / 1.01)) / (2 * np.log(1.01))
    true_orbit = sweeps["orbit"]["truth"].results()[0]
    hi, lo = (band_fluxes([true_orbit], [TRUE_EMISSIVITY * f], geoms["orbit"], band)[0][0]
              for f in (1.01, 1 / 1.01))
    d_ep_fixed = np.log(hi / lo) / (2 * np.log(1.01))
    d_alb = np.log(clean_wrong[rot] / clean[rot]) / np.log(WRONG_ALBEDO / ALBEDO)
    print("\nChange of the in-band flux at the rotation times for a 1 percent increase:")
    for label, d in (("thermal inertia", d_ti),
                     ("emissivity, temperatures held fixed", d_ep_fixed),
                     ("emissivity, through the energy balance", d_ep),
                     (f"Bond albedo ({ALBEDO} to {WRONG_ALBEDO})", d_alb)):
        print(f"  {label:40s} {d.mean():+.2f} percent")

    # The best-fit model at the observation times, for the figure.
    ti_b, ep_b = best["rotation and eclipse"]
    _, m_best = chi2_surface(log_ti, eps, flux_grid, obs, sigma,
                             np.array([ti_b]), np.array([ep_b]))
    tb_best = tb(m_best[:, 0, 0])

    # The true model every five minutes through the egress and every hour
    # through the orbit, for the figure.
    truth = {k: sweeps[k]["truth"].results()[0] for k in ("orbit", "eclipse")}
    dense = {
        "eclipse": [f"2025-01-01T{h:02d}:{m:02d}:00" for h in (20, 21, 22)
                    for m in range(0, 60, 5) if (20, 0) <= (h, m) <= (22, 30)],
        "orbit": [str(np.datetime64("2025-01-01T01:00:00") + np.timedelta64(h, "h"))
                  for h in range(0, 84)]}
    curves = {}
    for kind, tlist in dense.items():
        hours, values = [], []
        for t in tlist:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                g = observe.observer_geometry(truth[kind], t, observer="EARTH")
            hours.append(_hours([g])[0])
            if g.occultation != "visible":       # a gap while Jupiter hides Europa
                values.append(np.nan)
                continue
            em = observe.band_emission(
                truth[kind], observe.temperature_at(truth[kind], g), g, band=band,
                emissivity=TRUE_EMISSIVITY,
                direct_solar=observe.temperature_at(truth[kind], g, "direct_solar"))
            values.append(em.tb_k)
        curves[kind] = (np.array(hours), np.array(values))

    # The same egress from the one-orbit run, at 26.6-minute output steps.
    coarse = []
    for t in dense["eclipse"]:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            g = observe.observer_geometry(truth["orbit"], t, observer="EARTH")
        if g.occultation != "visible":
            coarse.append(np.nan)
            continue
        coarse.append(observe.band_emission(
            truth["orbit"], observe.temperature_at(truth["orbit"], g), g, band=band,
            emissivity=TRUE_EMISSIVITY,
            direct_solar=observe.temperature_at(truth["orbit"], g, "direct_solar")).tb_k)
    gap = np.abs(np.array(coarse) - curves["eclipse"][1])
    print(f"Egress: the one-orbit run (26.6-minute steps) differs from the "
          f"one-minute run by up to {np.nanmax(gap):.1f} K in band brightness "
          f"temperature.")

    np.savez(outdir / "fit_results.npz", ti_fine=ti_fine, eps_fine=eps_fine,
             grid_ti=np.array(GRID_TI, float), grid_emissivity=eps, flux_grid=flux_grid,
             obs_flux=obs, sigma_flux=sigma,
             chi2_rotation=surfaces["rotation only"],
             chi2_all=surfaces["rotation and eclipse"],
             chi2_wrong_albedo=surfaces["wrong albedo"],
             hours=_hours(geoms_all), n_rotation=n_rot, tb_obs=tb_obs,
             tb_err=tb_err, tb_best=tb_best, best_ti=ti_b, best_emissivity=ep_b,
             true_ti=TRUE_TI, true_emissivity=TRUE_EMISSIVITY,
             curve_orbit_hours=curves["orbit"][0], curve_orbit_tb=curves["orbit"][1],
             curve_eclipse_hours=curves["eclipse"][0],
             curve_eclipse_tb=curves["eclipse"][1])
    plot(outdir)
    print(f"Wrote {_shown(outdir)}/synthetic_observations.csv, fit_results.npz, "
          f"chi2_surface.png, model_vs_observations.png")


def _shown(path):
    """A path as printed: relative to the repository when it lies inside it,
    so the printed output does not depend on where the repository is."""
    path = Path(path).resolve()
    try:
        return path.relative_to(REPO.resolve()).as_posix()
    except ValueError:
        return str(path)


def _hours(geoms):
    return np.array([(np.datetime64(g.utc_observed[:19])
                      - np.datetime64("2025-01-01T00:00:00"))
                     / np.timedelta64(1, "h") for g in geoms])


def plot(outdir):
    """The two figures, from fit_results.npz (notebook 05 draws the same)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter, ScalarFormatter
    r = np.load(Path(outdir) / "fit_results.npz")

    fig, ax = plt.subplots(figsize=(6.0, 4.5))
    for key, label, color, style in (
            ("chi2_rotation", "rotation only", "#2a78d6", "--"),
            ("chi2_all", "rotation and eclipse", "#eb6834", "-")):
        c2 = r[key]
        ax.contour(r["ti_fine"], r["eps_fine"], c2 - c2.min(), levels=[2.30, 6.18],
                   colors=color, linewidths=[2.0, 1.0], linestyles=style)
        ax.plot([], [], style, color=color, lw=2, label=label)
    ax.plot([r["true_ti"]], [r["true_emissivity"]], marker="+", ms=12, mew=2,
            color="#0b0b0b", ls="none", label="true value")
    ax.set_xscale("log")
    ax.set_xticks([20, 30, 50, 70, 100, 150, 240])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r"Thermal inertia (J m$^{-2}$ K$^{-1}$ s$^{-1/2}$)")
    ax.set_ylabel("Emissivity")
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(Path(outdir) / "chi2_surface.png", dpi=200)
    plt.close(fig)

    n = int(r["n_rotation"])
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8), sharey=False,
                             gridspec_kw={"width_ratios": [1.6, 1.0]})
    for ax, sl, key in ((axes[0], slice(0, n), "orbit"), (axes[1], slice(n, None), "eclipse")):
        ax.plot(r[f"curve_{key}_hours"], r[f"curve_{key}_tb"], "-", color="#2a78d6",
                lw=1.5, label="true model")
        ax.errorbar(r["hours"][sl], r["tb_obs"][sl], yerr=r["tb_err"][sl], fmt="o",
                    ms=5, color="#eb6834", label="synthetic observations")
        ax.plot(r["hours"][sl], r["tb_best"][sl], "s", ms=6, mfc="none", mew=1.5,
                color="#0b0b0b", label="best fit")
    fig.supxlabel("Hours after 2025-01-01 00:00 UTC (received at Earth)", fontsize=10)
    axes[0].set_ylabel("Band brightness temperature, 8 to 14 µm (K)")
    axes[0].legend(frameon=False, loc="lower right", fontsize=9)
    fig.tight_layout()
    fig.savefig(Path(outdir) / "model_vs_observations.png", dpi=200)
    plt.close(fig)
    (Path(outdir) / "captions.md").write_text(
        "chi2_surface.png: regions of 68 and 95 percent confidence (two "
        "parameters, thick and thin contours) for thermal inertia and gray "
        "emissivity, fitted to synthetic disk-integrated 8-14 um observations "
        "of Europa with the Bond albedo held at 0.5, with the rotation times "
        "only and with the eclipse times added. The cross marks the values "
        "used to make the observations.\n\n"
        "model_vs_observations.png: band brightness temperature of the model "
        "used to make the synthetic observations (line), the observations with "
        "5 percent flux errors (circles), and the best-fitting model "
        "(squares), through one orbit (left) and as Europa leaves Jupiter's "
        "shadow (right).\n")


if __name__ == "__main__":
    main()
