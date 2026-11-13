#!/usr/bin/env python3
"""Fit thermal inertia and electrical skin depth to disk brightness temperatures.

The observations here are SYNTHETIC. They are made from one PTM3D simulation
of Europa with a known thermal inertia and electrical skin depth, plus noise,
and then fitted with a grid of other simulations. No telescope data are used.

    python3 examples/python/fit_europa_disk_tb_subsurface.py            # about 7 minutes on a CPU
    python3 examples/python/fit_europa_disk_tb_subsurface.py --no-run   # refit finished runs

Why the subsurface mode
-----------------------
At millimeter wavelengths ice is partly transparent. The emission that
leaves the surface comes from a layer about one electrical skin depth L_e
thick, and on Europa L_e is comparable to the diurnal thermal skin depth.
The brightness temperature then samples the diurnal temperature wave below
the surface, with a smaller amplitude and a later peak than the surface
temperature. A single effective emissivity cannot describe this, so this
example integrates the emission through the temperature profile of every
face at its own emission angle (``ptm.observe`` subsurface mode).

What it does
------------
1. Runs examples/quickstart_europa_diurnal.yaml at five values of a thermal
   inertia that is constant with depth and temperature, plus the "true"
   value. Each run writes the temperature at 28 depths, from the surface to
   4 m, at every output step (``output.map_depths_per_step``). As in
   fit_europa_disk_tb.py, each value has a one-orbit run for the rotation
   and a five-hour run with one-minute steps for the eclipse.
2. Computes the flux density of the visible disk at 233 and 97.5 GHz for a
   range of electrical skin depths. The refractive index and the Fresnel
   transmissivity of the surface come from the porous-ice dielectric model
   in scripts/emission (porosity 0.5, 100 K) and are not fitted. The skin
   depth is fitted at 233 GHz, and the one at 97.5 GHz follows from the
   dielectric model's ratio of the two. The loss tangent of ice rises with
   frequency, so the skin depth falls faster than the inverse of the
   frequency and the ratio is larger than the ratio of the frequencies.
3. Makes synthetic observations from the true run with the model's own skin
   depths and 1 percent noise, and fits thermal inertia and the 233 GHz skin
   depth by chi-square. The fit is repeated with four sets of observations
   to show what each adds.

Depth list and column
---------------------
The first depth spacing (2 mm) is under 0.16 times the shorter of the
electrical skin depth along the most oblique ray and the thermal skin depth,
for every value in the search. At 97.5 GHz the electrical skin depth is
about 40 cm, and the frames must reach about 7 skin depths along every ray
for the emission from below them to be negligible. The frames therefore go
to 4 m, in columns about 4.5 m deep set through ``thermal.thermal_b``. For
the longest skin depths in the search, the deepest frame still stands in for
about 1 percent of the emission, and the script reports how much.

Below the diurnal wave the temperature rises with depth at the rate set by
the internal heat flow, about 3 K per meter here. The 97.5 GHz emission
samples this gradient.

Needs: a built engine, the Europa SPICE kernels
(``python3 scripts/fetch_kernels.py --bodies europa``), scipy and matplotlib.
About 2.3 GB of disk for the runs.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import time
import warnings
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src" / "python"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from ptm import Sweep, observe  # noqa: E402
from fit_europa_disk_tb import (BASE, ECLIPSE_RUN, ECLIPSE_TIMES,  # noqa: E402
                                ROTATION_TIMES)

FREQS_GHZ = (233.0, 97.5)
POROSITY = 0.5
REFERENCE_T_K = 100.0
TRUE_TI = 70.0
GRID_TI = [30, 45, 60, 90, 130]          # excludes the true value
NOISE = 0.01                             # fractional, per observation
SEED = 20261116
DRAWS = 40                               # other noise draws, for the coverage

#: Depths of the temperature frames: every 2 mm to 3 cm, then 12 more,
#: spaced geometrically, to 4 m.
DEPTHS_M = [round(float(z), 4) for z in np.concatenate(
    [np.arange(0.0, 0.0301, 0.002), np.geomspace(0.03, 4.0, 13)[1:]])]

#: Column depth. With constant_thermal_inertia the diurnal skin depth of these
#: runs is ks sqrt(P / pi) / TI = 3.0 m / TI (ks = 0.0096 W m^-1 K^-1, P = one
#: Europa orbit), and the column is thermal_b skin depths deep.
COLUMN_M = 4.5
SKIN_DEPTH_TIMES_TI = 3.0

#: The 233 GHz electrical skin depths the fit searches, m.
SKIN_DEPTHS_233 = np.geomspace(0.02, 0.16, 57)

COLORS = {"233 GHz, rotation": "#2a78d6",
          "233 GHz, rotation and eclipse": "#eb6834",
          "233 and 97.5 GHz, rotation": "#1f9e74",
          "233 and 97.5 GHz, rotation and eclipse": "#0b0b0b"}


def model_skin_depths():
    """Electrical skin depth (m) and refractive index of porous ice at each
    frequency, from the dielectric model."""
    sub = observe.Subsurface(porosity=POROSITY, reference_temperature_k=REFERENCE_T_K)
    out = {}
    for f in FREQS_GHZ:
        kappa, n_real, _ = sub.resolve(f * 1e9)
        out[f] = (1.0 / kappa, n_real)
    return out


def run_simulations(outdir, threads, run):
    """Grid and true runs, for the orbit and for the eclipse."""
    base = {"thermal.threads": threads} if threads else {}
    base.update({"output.map_depths_m": DEPTHS_M,
                 "output.map_depth_snap_to_node": False,
                 "output.map_depths_per_step": True})
    derived = {"thermal.thermal_b": lambda p: int(math.ceil(
        COLUMN_M * p["constant_thermal_inertia"] / SKIN_DEPTH_TIMES_TI))}
    sweeps = {}
    for kind, extra in (("orbit", {}), ("eclipse", ECLIPSE_RUN)):
        overrides = dict(base, **extra)
        sweeps[kind] = (
            Sweep(BASE, grid={"constant_thermal_inertia": GRID_TI},
                  outdir=outdir / f"grid_{kind}", overrides=overrides, derived=derived),
            Sweep(BASE, grid={"constant_thermal_inertia": [TRUE_TI]},
                  outdir=outdir / f"truth_{kind}", overrides=overrides, derived=derived))
    if run:
        for grid, truth in sweeps.values():
            grid.run()
            truth.run()
    return sweeps


def disk_fluxes(run, times, skin_233, ratio):
    """Flux density (Jy) of one run, shape (frequency, time, skin depth),
    with the geometries. ``skin_233`` are 233 GHz skin depths, and each
    frequency's skin depth is ``ratio[f]`` times that."""
    mesh = observe.disk_mesh(run)
    delta = observe.thermal_skin_depth_m(run)
    flux = np.zeros((len(FREQS_GHZ), len(times), len(skin_233)))
    geoms, short = [], []
    for j, t in enumerate(times):
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)   # never fit a hidden epoch
            g = observe.observer_geometry(run, t, observer="EARTH")
        geoms.append(g)
        z, T = observe.depth_temperatures_at(run, g)
        for i, f in enumerate(FREQS_GHZ):
            for k, L in enumerate(skin_233):
                sub = observe.Subsurface(porosity=POROSITY,
                                         reference_temperature_k=REFERENCE_T_K,
                                         skin_depth_m=L * ratio[f])
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    em = observe.disk_emission(mesh, T, g, frequency_ghz=f,
                                               subsurface=sub, depths_m=z,
                                               thermal_skin_depth_m=delta)
                if caught:
                    short.append((f, L * ratio[f], em.bottom_weight_max,
                                  str(caught[0].message)))
                flux[i, j, k] = em.flux_density_jy
    return flux, geoms, short


def chi2_grid(log_ti, model, obs, sigma, ti_fine):
    """Chi-square on (thermal inertia, skin depth). Each observation's flux is
    a cubic spline in log thermal inertia at every skin depth."""
    from scipy.interpolate import CubicSpline
    m = CubicSpline(log_ti, model, axis=0)(np.log(ti_fine))   # (TI, obs, skin)
    return np.sum(((m - obs[None, :, None]) / sigma[None, :, None]) ** 2, axis=1)


def interval(values, profile):
    """Range where the profiled chi-square is within 1 of its minimum."""
    ok = values[profile <= profile.min() + 1.0]
    return ok.min(), ok.max()


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--outdir", default=str(REPO / "output" / "europa_disk_fit_subsurface"))
    parser.add_argument("--threads", type=int, default=0,
                        help="engine threads (0: all cores)")
    parser.add_argument("--no-run", action="store_true",
                        help="use finished runs, start no simulations")
    args = parser.parse_args()
    outdir = Path(args.outdir)
    start = time.time()

    sweeps = run_simulations(outdir, args.threads, not args.no_run)
    engine_s = time.time() - start
    model, truth = {}, {}
    for kind, (grid, true_sweep) in sweeps.items():
        runs = grid.results()
        ti = np.array(grid.values("constant_thermal_inertia"), dtype=float)
        order = np.argsort(ti)
        ti = ti[order]
        model[kind] = [runs[i] for i in order]
        truth[kind] = true_sweep.results()[0]
        for r in model[kind] + [truth[kind]]:
            depths = r.path / "statistics" / "depth_map_z.json"
            if not depths.is_file():
                raise SystemExit(f"{r.path}: no depth frames")

    props = model_skin_depths()
    L_true = props[FREQS_GHZ[0]][0]
    ratio = {f: props[f][0] / L_true for f in FREQS_GHZ}
    print("\nPorous ice, porosity 0.5, 100 K (dielectric model):")
    for f in FREQS_GHZ:
        print(f"  {f:6.1f} GHz: electrical skin depth {100 * props[f][0]:.2f} cm, "
              f"refractive index {props[f][1]:.3f}")

    sets_times = {"orbit": ROTATION_TIMES, "eclipse": ECLIPSE_TIMES}
    grid_flux, true_flux, geoms, short = [], [], [], []
    for kind, times in sets_times.items():
        per_run = [disk_fluxes(r, times, SKIN_DEPTHS_233, ratio) for r in model[kind]]
        grid_flux.append(np.stack([p[0] for p in per_run]))      # (TI, freq, time, skin)
        short += [s for p in per_run for s in p[2]]
        tf, g, s = disk_fluxes(truth[kind], times, np.array([L_true]), ratio)
        true_flux.append(tf[:, :, 0])
        geoms += g
        short += s
    grid_flux = np.concatenate(grid_flux, axis=2)
    true_flux = np.concatenate(true_flux, axis=1)                 # (freq, time)
    n_times = len(geoms)
    if short:
        worst = max(short, key=lambda s: s[2])
        print(f"\n{len(short)} model evaluations warned about the depth list. "
              f"Largest weight below the column {100 * worst[2]:.3f} percent, at "
              f"{worst[0]} GHz with a skin depth of {100 * worst[1]:.1f} cm:\n  {worst[3]}")
    else:
        print("\nNo model evaluation warned about the depth list or column depth.")

    rng = np.random.default_rng(SEED)
    sigma = NOISE * true_flux
    obs = true_flux + sigma * rng.standard_normal(true_flux.shape)

    # Observations are indexed (frequency, time), flattened in that order.
    model_flat = grid_flux.reshape(len(ti), -1, len(SKIN_DEPTHS_233))
    obs_flat, sigma_flat = obs.ravel(), sigma.ravel()
    n_rot = len(ROTATION_TIMES)
    is_rot = np.tile(np.arange(n_times) < n_rot, len(FREQS_GHZ))
    is_233 = np.repeat(np.arange(len(FREQS_GHZ)) == 0, n_times)
    sets = {"233 GHz, rotation": is_233 & is_rot,
            "233 GHz, rotation and eclipse": is_233,
            "233 and 97.5 GHz, rotation": is_rot,
            "233 and 97.5 GHz, rotation and eclipse": np.ones_like(is_rot)}

    # How well a spline through the grid stands in for a run at the true value.
    from scipy.interpolate import CubicSpline
    k_true = int(np.argmin(np.abs(SKIN_DEPTHS_233 - L_true)))
    spline_true = CubicSpline(np.log(ti), model_flat[:, :, k_true], axis=0)(math.log(TRUE_TI))
    exact_true = []
    for kind, times in sets_times.items():
        exact_true.append(disk_fluxes(truth[kind], times,
                                      np.array([SKIN_DEPTHS_233[k_true]]), ratio)[0][:, :, 0])
    exact_true = np.concatenate(exact_true, axis=1).ravel()
    interp_err = float(np.max(np.abs(spline_true - exact_true) / sigma_flat))

    ti_fine = np.geomspace(ti[0], ti[-1], 300)

    def fit(observed, sel):
        c2 = chi2_grid(np.log(ti), model_flat[:, sel, :], observed[sel],
                       sigma_flat[sel], ti_fine)
        t_best, k_best = np.unravel_index(np.argmin(c2), c2.shape)
        return c2, (ti_fine[t_best], *interval(ti_fine, c2.min(axis=1)),
                    100 * SKIN_DEPTHS_233[k_best],
                    *(100 * v for v in interval(SKIN_DEPTHS_233, c2.min(axis=0))),
                    c2.min())

    print(f"\nSYNTHETIC observations: true thermal inertia {TRUE_TI:g}, true "
          f"233 GHz skin depth {100 * L_true:.2f} cm, {100 * NOISE:.0f} percent noise")
    surfaces, rows = {}, []
    for label, observed in (("with noise", obs_flat),
                            ("without noise (expected widths)", true_flux.ravel())):
        print(f"\nFit {label}:")
        for name, sel in sets.items():
            c2, r = fit(observed, sel)
            if label == "with noise":
                surfaces[name] = c2
            rows.append([name, label, int(sel.sum()), *r])
            print(f"  {name:40s} ({sel.sum():2d} obs): thermal inertia {r[0]:6.1f} "
                  f"[{r[1]:.1f}, {r[2]:.1f}]   L_e(233) {r[3]:5.2f} cm "
                  f"[{r[4]:.2f}, {r[5]:.2f}]   min chi-square {r[6]:.1f}")
    print("\nIntervals: profiled chi-square within 1 of its minimum. Interval ends at "
          f"{ti[0]:g} or {ti[-1]:g} (thermal inertia) or {100 * SKIN_DEPTHS_233[0]:g} or "
          f"{100 * SKIN_DEPTHS_233[-1]:g} cm (skin depth) are the edges of the search.")
    rng_draws = np.random.default_rng(SEED + 1)
    hits_ti = hits_le = 0
    for _ in range(DRAWS):
        drawn = (true_flux + sigma * rng_draws.standard_normal(true_flux.shape)).ravel()
        _, r = fit(drawn, sets["233 and 97.5 GHz, rotation and eclipse"])
        hits_ti += r[1] <= TRUE_TI <= r[2]
        hits_le += r[4] <= 100 * L_true <= r[5]
    print(f"With all observations and {DRAWS} other noise draws, the interval contains "
          f"the true thermal inertia {hits_ti} times and the true skin depth "
          f"{hits_le} times (about 68 percent expected).")
    print(f"Spline in thermal inertia vs the true run, at its skin depth: worst "
          f"{interp_err:.3f} sigma.")

    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "fit_table.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["# SYNTHETIC observations, true thermal inertia "
                    f"{TRUE_TI}, true 233 GHz skin depth {100 * L_true:.2f} cm"])
        w.writerow(["set", "observations", "n_obs", "ti_best", "ti_lo", "ti_hi",
                    "le233_best_cm", "le233_lo_cm", "le233_hi_cm", "chi2_min"])
        for r in rows:
            w.writerow([r[0], r[1], r[2]] + [f"{v:.3f}" for v in r[3:]])
    with open(outdir / "synthetic_observations.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["# SYNTHETIC: constant thermal inertia "
                    f"{TRUE_TI}, porous ice porosity {POROSITY}, skin depths "
                    + ", ".join(f"{100 * props[f][0]:.2f} cm at {f} GHz" for f in FREQS_GHZ)
                    + f", {100 * NOISE:.0f} percent noise, seed {SEED}"])
        w.writerow(["utc_observed", "set", "frequency_ghz", "sub_observer_lon_deg",
                    "flux_density_jy", "sigma_jy", "tb_planck_k"])
        for i, f in enumerate(FREQS_GHZ):
            for j, g in enumerate(geoms):
                omega = math.pi * (observe.disk_mesh(truth["orbit"]).radius_m
                                   / (1000.0 * g.distance_km)) ** 2
                tb = observe.planck_nu_temperature(obs[i, j] * observe.JANSKY / omega,
                                                   f * 1e9)
                w.writerow([g.utc_observed, "rotation" if j < n_rot else "eclipse", f,
                            f"{g.sub_observer_lon:.3f}", f"{obs[i, j]:.6f}",
                            f"{sigma[i, j]:.6f}", f"{tb:.3f}"])

    plot_valley(outdir, ti_fine, surfaces, L_true)
    plot_eclipse(outdir, truth["eclipse"], ratio, L_true, geoms[n_rot:], obs[:, n_rot:],
                 sigma[:, n_rot:])
    print(f"\nWrote {outdir}/fit_table.csv, synthetic_observations.csv, "
          f"valley.png, eclipse_tb.png, captions.md")
    print(f"Wall time: engine runs {engine_s:.0f} s, total {time.time() - start:.0f} s.")


def plot_valley(outdir, ti_fine, surfaces, L_true):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    for name, c2 in surfaces.items():
        ax.contour(ti_fine, 100 * SKIN_DEPTHS_233, (c2 - c2.min()).T, levels=[2.30],
                   colors=COLORS[name], linewidths=2.0)
        ax.plot([], [], color=COLORS[name], lw=2, label=name)
    ax.plot([TRUE_TI], [100 * L_true], marker="+", ms=12, mew=2, color="#d62728",
            ls="none", label="true value")
    from matplotlib.ticker import NullFormatter, ScalarFormatter
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([30, 40, 60, 80, 100, 130])
    ax.set_yticks([2, 3, 4, 6, 8, 10, 16])
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(ScalarFormatter())
        axis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r"Thermal inertia (J m$^{-2}$ K$^{-1}$ s$^{-1/2}$)")
    ax.set_ylabel("Electrical skin depth at 233 GHz (cm)")
    ax.set_title("Fit to synthetic observations (68 percent regions)")
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(outdir / "valley.png", dpi=200)
    plt.close(fig)


def plot_eclipse(outdir, run, ratio, L_true, geoms_obs, obs, sigma):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dense = [f"2025-01-01T{h:02d}:{m:02d}:00" for h in (20, 21, 22)
             for m in range(60) if (20, 0) <= (h, m) <= (22, 30)]
    flux, geoms, _ = disk_fluxes(run, dense, np.array([L_true]), ratio)
    radius_m = observe.disk_mesh(run).radius_m

    def tb(f_jy, g, f_ghz):
        omega = math.pi * (radius_m / (1000.0 * g.distance_km)) ** 2
        return observe.planck_nu_temperature(f_jy * observe.JANSKY / omega, f_ghz * 1e9)

    def hours(gs):
        return np.array([(np.datetime64(g.utc_observed[:19])
                          - np.datetime64("2025-01-01T00:00:00"))
                         / np.timedelta64(1, "h") for g in gs])
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    styles = (("#2a78d6", "-", "o"), ("#eb6834", "--", "s"))
    for i, (f, (color, ls, marker)) in enumerate(zip(FREQS_GHZ, styles)):
        ax.plot(hours(geoms), [tb(flux[i, j, 0], g, f) for j, g in enumerate(geoms)],
                color=color, ls=ls, lw=2, label=f"true model, {f:g} GHz")
        t_obs = np.array([tb(obs[i, j], g, f) for j, g in enumerate(geoms_obs)])
        ax.errorbar(hours(geoms_obs), t_obs, yerr=t_obs * sigma[i] / obs[i], fmt=marker,
                    ms=5, color=color, label=f"synthetic observations, {f:g} GHz")
    ax.set_xlabel("Hours after 2025-01-01 00:00 UTC (received at Earth)")
    ax.set_ylabel("Disk brightness temperature (K)")
    ax.set_title("Europa leaving Jupiter's shadow, subsurface emission")
    ax.legend(frameon=False, fontsize=8, ncol=2, loc="upper center",
              bbox_to_anchor=(0.5, -0.17))
    fig.tight_layout()
    fig.savefig(outdir / "eclipse_tb.png", dpi=200)
    plt.close(fig)
    (outdir / "captions.md").write_text(
        "valley.png: regions of 68 percent confidence (two parameters) for thermal "
        "inertia and the electrical skin depth at 233 GHz, fitted to synthetic "
        "disk-integrated observations. With one frequency and the rotation alone "
        "the region is a long valley, because a higher thermal inertia and a "
        "longer skin depth both reduce the diurnal variation. The eclipse times "
        "and the second frequency each narrow it. The cross marks the values used "
        "to make the observations.\n\n"
        "eclipse_tb.png: disk brightness temperature (Planck) at 233 and 97.5 GHz "
        "of the model used to make the synthetic observations, and the "
        "observations with 1 percent errors, as Europa leaves Jupiter's shadow. "
        "The emission at 97.5 GHz comes from deeper ice, which the eclipse "
        "barely cools.\n")


if __name__ == "__main__":
    main()
