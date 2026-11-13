#!/usr/bin/env python3
"""Fit thermal inertia and emissivity to disk-integrated brightness temperatures.

The observations here are SYNTHETIC. They are made from one PTM3D simulation
of Europa with a known thermal inertia and emissivity, plus noise, and then
fitted with a grid of other simulations, so you can see how well the two
parameters are recovered and how they trade off. No telescope data are used.

    python3 examples/python/fit_europa_disk_tb.py            # about 2 minutes on a CPU
    python3 examples/python/fit_europa_disk_tb.py --no-run   # refit finished runs

What it does
------------
1. Runs examples/quickstart_europa_diurnal.yaml (one Europa orbit starting
   2025-01-01, ico4 mesh) at several values of a thermal inertia that is
   constant with depth and temperature (``constant_thermal_inertia``; this
   mode also turns off the temperature dependence of conductivity), plus one
   run at the "true" value used to make the observations.
2. Places Earth at a set of observation times (when the light arrives) and
   computes the flux density of the whole visible disk at 233 GHz for every
   run with ``ptm.observe``. Half of the times sample Europa's rotation; the
   other half follow Europa out of Jupiter's shadow on 2025-01-01, while it
   is visible from Earth.
3. Makes synthetic observations from the true run with an emissivity of
   0.85 and 1 percent noise, and fits thermal inertia and emissivity by
   chi-square on a fine grid, with and without the eclipse times.

In this surface mode the emission comes from the surface temperature and
the emissivity is an effective value for the observing band. It scales the
flux exactly, so it needs no extra simulations. It is not the
thermal-infrared emissivity of the surface energy balance
(``thermophysical.emissivity_override``), which changes the temperatures
and needs a new simulation for every value.

Needs: a built engine, the Europa SPICE kernels
(``python3 scripts/fetch_kernels.py --bodies europa``), scipy and matplotlib.
"""

from __future__ import annotations

import argparse
import csv
import sys
import warnings
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src" / "python"))

from ptm import Results, Sweep, observe  # noqa: E402

BASE = REPO / "examples" / "quickstart_europa_diurnal.yaml"
FREQ_GHZ = 233.0
TRUE_TI = 70.0
TRUE_EMISSIVITY = 0.85
NOISE = 0.01                      # fractional, per observation
GRID_TI = [20, 35, 50, 85, 120, 170, 240]   # excludes the true value
SEED = 20261116

#: Observation (receive) times, UTC. Europa is visible from Earth at each.
ROTATION_TIMES = ["2025-01-01T08:00:00", "2025-01-02T04:00:00",
                  "2025-01-02T16:00:00", "2025-01-03T04:00:00",
                  "2025-01-03T16:00:00", "2025-01-04T04:00:00"]
ECLIPSE_TIMES = ["2025-01-01T20:10:00", "2025-01-01T20:30:00",
                 "2025-01-01T20:50:00", "2025-01-01T21:10:00",
                 "2025-01-01T21:40:00", "2025-01-01T22:20:00"]


def run_simulations(outdir, threads, run):
    """The model grid and the true run, each a constant_thermal_inertia sweep."""
    common = {"overrides": {"thermal.threads": threads}} if threads else {}
    grid = Sweep(BASE, grid={"constant_thermal_inertia": GRID_TI},
                 outdir=outdir / "grid", **common)
    truth = Sweep(BASE, grid={"constant_thermal_inertia": [TRUE_TI]},
                  outdir=outdir / "truth", **common)
    if run:
        grid.run()
        truth.run()
    return grid, truth


def disk_fluxes(runs, times):
    """Flux density (Jy) with emissivity 1, one row per run, one column per time."""
    geoms = []
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)   # never fit a hidden epoch
        for t in times:
            geoms.append(observe.observer_geometry(runs[0], t, observer="EARTH"))
    flux = np.array([[observe.disk_emission(r, observe.temperature_at(r, g), g,
                                            emissivity=1.0,
                                            frequency_ghz=FREQ_GHZ).flux_density_jy
                      for g in geoms] for r in runs])
    return flux, geoms


def chi2_surface(log_ti_grid, flux_grid, obs, sigma, ti_fine, eps_fine):
    """Chi-square on (thermal inertia, emissivity). Each time's emissivity-1
    flux is a cubic spline in log thermal inertia; emissivity scales it."""
    from scipy.interpolate import CubicSpline
    model1 = np.array([CubicSpline(log_ti_grid, flux_grid[:, j])(np.log(ti_fine))
                       for j in range(flux_grid.shape[1])])      # (time, TI)
    resid = (eps_fine[:, None, None] * model1[None, :, :]
             - obs[None, :, None]) / sigma[None, :, None]
    return np.sum(resid ** 2, axis=1)                           # (eps, TI)


def interval(values, chi2_profile):
    """Range where the profiled chi-square is within 1 of its minimum."""
    ok = values[chi2_profile <= chi2_profile.min() + 1.0]
    return ok.min(), ok.max()


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--outdir", default=str(REPO / "output" / "europa_disk_fit"))
    parser.add_argument("--threads", type=int, default=0,
                        help="engine threads (0: all cores)")
    parser.add_argument("--no-run", action="store_true",
                        help="use finished runs, start no simulations")
    args = parser.parse_args()
    outdir = Path(args.outdir)

    grid, truth = run_simulations(outdir, args.threads, not args.no_run)
    model_runs = grid.results()
    true_run = truth.results()[0]
    ti = np.array(grid.values("constant_thermal_inertia"), dtype=float)
    order = np.argsort(ti)
    ti = ti[order]
    model_runs = [model_runs[i] for i in order]

    # The swept parameter must actually have changed the runs.
    t0 = [r.field("temperature", 0) for r in model_runs]
    if any(np.array_equal(t0[0], t) for t in t0[1:]):
        raise SystemExit("two thermal-inertia points have identical outputs; "
                         "the swept key had no effect")

    times = ROTATION_TIMES + ECLIPSE_TIMES
    flux_grid, geoms = disk_fluxes(model_runs, times)
    flux_true, _ = disk_fluxes([true_run], times)

    rng = np.random.default_rng(SEED)
    clean = TRUE_EMISSIVITY * flux_true[0]
    sigma = NOISE * clean
    obs = clean + sigma * rng.standard_normal(len(times))

    nu = FREQ_GHZ * 1e9
    omega = np.array([observe.disk_emission(true_run, observe.temperature_at(true_run, g),
                                            g, emissivity=1.0, frequency_ghz=FREQ_GHZ
                                            ).omega_nominal_sr for g in geoms])
    tb_obs = observe.planck_nu_temperature(obs * observe.JANSKY / omega, nu)

    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "synthetic_observations.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["# SYNTHETIC: constant thermal inertia "
                    f"{TRUE_TI}, emissivity {TRUE_EMISSIVITY}, "
                    f"{100 * NOISE:.0f} percent noise, seed {SEED}"])
        w.writerow(["utc_observed", "set", "sub_observer_lon_deg",
                    "phase_angle_deg", "flux_density_jy", "sigma_jy",
                    "tb_planck_k"])
        for k, g in enumerate(geoms):
            w.writerow([g.utc_observed, "eclipse" if times[k] in ECLIPSE_TIMES
                        else "rotation", f"{g.sub_observer_lon:.3f}",
                        f"{g.phase_angle_deg:.3f}", f"{obs[k]:.6f}",
                        f"{sigma[k]:.6f}", f"{tb_obs[k]:.3f}"])

    ti_fine = np.geomspace(ti[0], ti[-1], 400)
    eps_fine = np.linspace(0.60, 1.00, 401)
    sets = {"rotation only": np.arange(len(ROTATION_TIMES)),
            "rotation and eclipse": np.arange(len(times))}
    print("\nSYNTHETIC observations: true thermal inertia "
          f"{TRUE_TI:g}, true emissivity {TRUE_EMISSIVITY}, "
          f"{100 * NOISE:.0f} percent noise\n")
    surfaces = {}
    for name, idx in sets.items():
        c2 = chi2_surface(np.log(ti), flux_grid[:, idx], obs[idx], sigma[idx],
                          ti_fine, eps_fine)
        surfaces[name] = c2
        e_best, t_best = np.unravel_index(np.argmin(c2), c2.shape)
        ti_lo, ti_hi = interval(ti_fine, c2.min(axis=0))
        ep_lo, ep_hi = interval(eps_fine, c2.min(axis=1))
        print(f"{name:22s} ({len(idx)} times): thermal inertia "
              f"{ti_fine[t_best]:6.1f}  [{ti_lo:.1f}, {ti_hi:.1f}]   "
              f"emissivity {eps_fine[e_best]:.3f}  [{ep_lo:.3f}, {ep_hi:.3f}]   "
              f"min chi-square {c2.min():.1f}")
    print("\nIntervals: profiled chi-square within 1 of its minimum.")

    # The true model every 5 minutes through the egress, for the figure.
    dense = [f"2025-01-01T{h:02d}:{m:02d}:00" for h in (20, 21, 22)
             for m in range(0, 60, 5) if (20, 10) <= (h, m) <= (22, 25)]
    flux_dense, geoms_dense = disk_fluxes([true_run], dense)
    omega_dense = np.array([observe.disk_emission(
        true_run, observe.temperature_at(true_run, g), g, emissivity=1.0,
        frequency_ghz=FREQ_GHZ).omega_nominal_sr for g in geoms_dense])
    tb_dense = observe.planck_nu_temperature(
        TRUE_EMISSIVITY * flux_dense[0] * observe.JANSKY / omega_dense, nu)

    plot(outdir, ti_fine, eps_fine, surfaces, geoms, times, tb_obs,
         geoms_dense, tb_dense)
    print(f"Wrote {outdir}/synthetic_observations.csv, chi2_surface.png, "
          f"eclipse_tb.png")


def _hours(geoms):
    return np.array([(np.datetime64(g.utc_observed[:19])
                      - np.datetime64("2025-01-01T00:00:00"))
                     / np.timedelta64(1, "h") for g in geoms])


def plot(outdir, ti_fine, eps_fine, surfaces, geoms, times, tb_obs,
         geoms_dense, tb_dense):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = {"rotation only": "#2a78d6", "rotation and eclipse": "#eb6834"}

    fig, ax = plt.subplots(figsize=(6.0, 4.5))
    for name, c2 in surfaces.items():
        ax.contour(ti_fine, eps_fine, c2 - c2.min(), levels=[2.30, 6.18],
                   colors=colors[name], linewidths=[2.0, 1.0])
        ax.plot([], [], color=colors[name], lw=2, label=name)
    ax.plot([TRUE_TI], [TRUE_EMISSIVITY], marker="+", ms=12, mew=2,
            color="#0b0b0b", ls="none", label="true value")
    ax.set_xscale("log")
    ax.set_ylim(0.75, 0.95)
    ax.set_xlabel(r"Thermal inertia (J m$^{-2}$ K$^{-1}$ s$^{-1/2}$)")
    ax.set_ylabel("Emissivity at 233 GHz")
    ax.set_title("Fit to synthetic disk-integrated observations")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "chi2_surface.png", dpi=200)
    plt.close(fig)

    ecl = [k for k, t in enumerate(times) if t in ECLIPSE_TIMES]
    hours = _hours([geoms[k] for k in ecl])
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    ax.plot(_hours(geoms_dense), tb_dense, color="#2a78d6", lw=2,
            label="true model")
    ax.errorbar(hours, tb_obs[ecl],
                yerr=tb_obs[ecl] * NOISE, fmt="o", ms=6, color="#eb6834",
                label="synthetic observations")
    ax.set_xlabel("Hours after 2025-01-01 00:00 UTC (received at Earth)")
    ax.set_ylabel("Disk brightness temperature (K)")
    ax.set_title("Europa leaving Jupiter's shadow, 233 GHz")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "eclipse_tb.png", dpi=200)
    plt.close(fig)
    (outdir / "captions.md").write_text(
        "chi2_surface.png: regions of 68 and 95 percent confidence (two "
        "parameters) for thermal inertia and 233 GHz emissivity, fitted to "
        "synthetic disk-integrated observations, with and without the times "
        "after the Jupiter eclipse. The cross marks the values used to make "
        "the observations.\n\n"
        "eclipse_tb.png: disk brightness temperature (Planck) of the model "
        "used to make the synthetic observations, and the observations with "
        "1 percent errors, as Europa leaves Jupiter's shadow.\n")


if __name__ == "__main__":
    main()
