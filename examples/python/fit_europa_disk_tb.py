#!/usr/bin/env python3
"""Fit thermal inertia and emissivity to disk-integrated brightness temperatures.

The observations here are SYNTHETIC. They are made from one PTM3D simulation
of Europa with a known thermal inertia and emissivity, plus noise, and then
fitted with a grid of other simulations, so you can see how well the two
parameters are recovered and how they trade off. No telescope data are used.

    python3 examples/python/fit_europa_disk_tb.py            # about 6 minutes on a CPU
    python3 examples/python/fit_europa_disk_tb.py --no-run   # refit finished runs

What it does
------------
1. Runs examples/quickstart_europa_diurnal.yaml (one Europa orbit starting
   2025-01-01, ico4 mesh) at several values of a thermal inertia that is
   constant with depth and temperature (``constant_thermal_inertia``, which
   also turns off the temperature dependence of conductivity), plus one run
   at the "true" value used to make the observations. Each value also has a
   second, five-hour run with one-minute output steps
   that starts before Europa enters Jupiter's shadow on 2025-01-01 and ends
   after it leaves.
2. Places Earth at a set of observation times (when the light arrives) and
   computes the flux density of the whole visible disk at 233 GHz for every
   run with ``ptm.observe``. Half of the times sample Europa's rotation and
   use the one-orbit runs. The other half follow Europa out of Jupiter's
   shadow, while it is visible from Earth, and use the one-minute runs.
3. Makes synthetic observations from the true run with an emissivity of
   0.85 and 1 percent noise, and fits thermal inertia and emissivity by
   chi-square on a fine grid, with and without the eclipse times.

In this surface mode the emission comes from the surface temperature and
the emissivity is an effective value for the observing band. It scales the
flux exactly, so it needs no extra simulations. It is not the
thermal-infrared emissivity of the surface energy balance
(``thermophysical.emissivity_override``), which changes the temperatures
and needs a new simulation for every value.

Output cadence near an eclipse
------------------------------
Europa's surface cools and warms quickly at the edges of the shadow, and the
disk brightness temperature rises by about 2 K per minute at egress. A model
used to fit eclipse data needs output steps of 1 minute or less within 10
minutes of the shadow contacts, and 5 minutes or less through the rest of the
eclipse.

The one-orbit runs step every 26.6 minutes. Near egress their disk
brightness temperature is more than 3 K from that of the one-minute runs,
several times the noise. Halving the step to 30 seconds changes it by less
than 0.04 K.

The one-minute runs sample their spin-up as the one-orbit runs do
(``thermal.equilibration_samples_per_diurnal_cycle: 192``). Without that
setting the spin-up also takes one-minute steps and the run is about thirty
times longer.

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

#: The eclipse runs: one-minute output steps from 21 minutes before Europa
#: enters the shadow (18:06 UTC at Europa) to 22:45, past the last eclipse
#: observation (21:45 at Europa, 22:20 at Earth). The spin-up is sampled as
#: in the one-orbit runs.
ECLIPSE_RUN = {"simulation.start_time": "2025-01-01T17:45:00",
               "simulation.duration_days": 300 / 1440,
               "simulation.timesteps": 300,
               "thermal.equilibration_samples_per_diurnal_cycle": 192}


def run_simulations(outdir, threads, run):
    """The model grid and the true run, each a constant_thermal_inertia sweep,
    for the whole orbit and for the eclipse at one-minute steps."""
    base = {"thermal.threads": threads} if threads else {}
    sweeps = {}
    for kind, extra in (("orbit", {}), ("eclipse", ECLIPSE_RUN)):
        overrides = dict(base, **extra)
        sweeps[kind] = (
            Sweep(BASE, grid={"constant_thermal_inertia": GRID_TI},
                  outdir=outdir / f"grid_{kind}", overrides=overrides),
            Sweep(BASE, grid={"constant_thermal_inertia": [TRUE_TI]},
                  outdir=outdir / f"truth_{kind}", overrides=overrides))
    if run:
        for grid, truth in sweeps.values():
            grid.run()
            truth.run()
    return sweeps


def disk_fluxes(runs, times):
    """Flux density (Jy) with emissivity 1, one row per run and one column per
    time, with the geometries and the disk's nominal solid angle (sr) at each
    time. All runs must share one time axis."""
    geoms = []
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)   # never fit a hidden epoch
        for t in times:
            geoms.append(observe.observer_geometry(runs[0], t, observer="EARTH"))
    emission = [[observe.disk_emission(r, observe.temperature_at(r, g), g,
                                       emissivity=1.0, frequency_ghz=FREQ_GHZ)
                 for g in geoms] for r in runs]
    flux = np.array([[e.flux_density_jy for e in row] for row in emission])
    omega = np.array([e.omega_nominal_sr for e in emission[0]])
    return flux, geoms, omega


def chi2_surface(log_ti_grid, flux_grid, obs, sigma, ti_fine, eps_fine):
    """Chi-square on (thermal inertia, emissivity). Each time's emissivity-1
    flux is a cubic spline in log thermal inertia, and emissivity scales it."""
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

    sweeps = run_simulations(outdir, args.threads, not args.no_run)
    model, true_run = {}, {}
    for kind, (grid, truth) in sweeps.items():
        runs = grid.results()
        ti = np.array(grid.values("constant_thermal_inertia"), dtype=float)
        order = np.argsort(ti)
        ti = ti[order]
        model[kind] = [runs[i] for i in order]
        true_run[kind] = truth.results()[0]
        # The swept parameter must actually have changed the runs.
        t0 = [r.field("temperature", 0) for r in model[kind]]
        if any(np.array_equal(t0[0], t) for t in t0[1:]):
            raise SystemExit("two thermal-inertia points have identical outputs, "
                             "so the swept key had no effect")

    times = ROTATION_TIMES + ECLIPSE_TIMES
    flux_rot, geoms_rot, omega_rot = disk_fluxes(model["orbit"], ROTATION_TIMES)
    flux_ecl, geoms_ecl, omega_ecl = disk_fluxes(model["eclipse"], ECLIPSE_TIMES)
    flux_grid = np.hstack([flux_rot, flux_ecl])
    geoms = geoms_rot + geoms_ecl
    omega = np.concatenate([omega_rot, omega_ecl])
    flux_true = np.hstack([disk_fluxes([true_run["orbit"]], ROTATION_TIMES)[0],
                           disk_fluxes([true_run["eclipse"]], ECLIPSE_TIMES)[0]])

    rng = np.random.default_rng(SEED)
    clean = TRUE_EMISSIVITY * flux_true[0]
    sigma = NOISE * clean
    obs = clean + sigma * rng.standard_normal(len(times))

    nu = FREQ_GHZ * 1e9
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

    # The true model every minute through the egress, for the figure.
    dense = [f"2025-01-01T{h:02d}:{m:02d}:00" for h in (20, 21, 22)
             for m in range(60) if (20, 0) <= (h, m) <= (22, 30)]
    tb_dense = {}
    for kind in ("eclipse", "orbit"):
        flux_dense, geoms_dense, omega_dense = disk_fluxes([true_run[kind]], dense)
        tb_dense[kind] = observe.planck_nu_temperature(
            TRUE_EMISSIVITY * flux_dense[0] * observe.JANSKY / omega_dense, nu)
    worst = int(np.argmax(np.abs(tb_dense["orbit"] - tb_dense["eclipse"])))
    print(f"Egress: the one-orbit runs (26.6-minute steps) differ from the "
          f"one-minute runs by up to "
          f"{abs(tb_dense['orbit'][worst] - tb_dense['eclipse'][worst]):.2f} K "
          f"(at {dense[worst]} UTC).")

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
    ax.plot(_hours(geoms_dense), tb_dense["eclipse"], color="#2a78d6", lw=2,
            label="true model, 1-minute steps")
    ax.plot(_hours(geoms_dense), tb_dense["orbit"], color="#2a78d6", lw=1.2,
            ls="--", label="true model, 26.6-minute steps")
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
        "1 percent errors, as Europa leaves Jupiter's shadow. The solid curve "
        "uses one-minute output steps, as the fit does. The dashed curve is "
        "the same model with the 26.6-minute steps of the one-orbit runs, "
        "linearly interpolated in time.\n")


if __name__ == "__main__":
    main()
