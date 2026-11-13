"""The ``ptm`` command: look at a finished run without writing any Python.

Eight subcommands::

    ptm info   <run_dir>                  what the run is, and what it does not record
    ptm plot   <run_dir> <what> [--save]  diurnal | map | globes | hist
    ptm run    <config.yaml>              start the engine on a config
    ptm sweep  <config.yaml> --vary K=V,V run a grid of them
    ptm advise <config.yaml>              is lateral heat conduction needed for this mesh?
    ptm cadence <run_dir>                 did the output step miss sunlit flashes?
    ptm gpu-selfcheck                     is GPU illumination right on this GPU?
    ptm observe <run_dir> --time UTC ...  disk-integrated emission seen by an observer

``ptm info`` and ``ptm plot`` close the gap the install audit measured: the
documented quickstart ends at a stdlib snippet printing three floats, and every
route from there to a picture needed a second environment (the GUI's twelve
unpinned dependencies, Jupyter -- which is named in no requirements file -- or
ParaView).

``ptm run`` is a single serial run and nothing more. ``ptm sweep`` is the same
engine driver over a grid, with the resume, provenance and collection that the
nineteen bespoke sweep scripts in this tree each re-implemented.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

__all__ = ["main"]

_PLOT_KINDS = ("diurnal", "map", "globes", "hist")


def _add_common(parser):
    parser.add_argument("run_dir", help="a finished PTM3D output directory")


def _build_parser():
    parser = argparse.ArgumentParser(
        prog="ptm",
        description="Inspect and plot PTM3D runs, and start the engine.")
    parser.add_argument("--version", action="store_true",
                        help="print the ptm package version and exit")
    sub = parser.add_subparsers(dest="command")

    info = sub.add_parser(
        "info", help="print what a run directory records (and what it does not)")
    _add_common(info)
    info.add_argument("--json", action="store_true",
                      help="emit provenance() as JSON instead of a report")

    plot = sub.add_parser("plot", help="draw one of the canned figures")
    _add_common(plot)
    plot.add_argument("what", choices=_PLOT_KINDS,
                      help="diurnal (one face over time), map (regridded "
                           "lat/lon), globes (orthographic set), hist "
                           "(distribution)")
    plot.add_argument("--field", default="temperature",
                      help="field name (default: temperature)")
    plot.add_argument("--moment", default="max", choices=("avg", "max", "min"),
                      help="which statistic to draw (default: max); ignored "
                           "with --step")
    plot.add_argument("--step", type=int, default=None,
                      help="draw one output step instead of a statistic")
    plot.add_argument("--lat", type=float, default=None,
                      help="latitude for 'diurnal' (degrees)")
    plot.add_argument("--lon", type=float, default=None,
                      help="longitude for 'diurnal' (degrees)")
    plot.add_argument("--face", type=int, default=None,
                      help="face index for 'diurnal', instead of --lat/--lon")
    plot.add_argument("--resolution", type=float, default=None,
                      help="map grid resolution in degrees (default: estimated "
                           "from the mesh)")
    plot.add_argument("--units", default="ascii", choices=("ascii", "unicode"),
                      help="unit spelling in labels: W/m^2 or W/m2 with a "
                           "superscript (default: ascii)")
    plot.add_argument("--mesh-units", default=None, choices=("m", "km"),
                      help="override the mesh length unit if it cannot be "
                           "inferred")
    plot.add_argument("--save", default=None, metavar="FILE",
                      help="write the figure here (extension picks the format)")
    plot.add_argument("--show", action="store_true",
                      help="open an interactive window")
    plot.add_argument("--dpi", type=int, default=200, help="raster DPI (200)")

    run = sub.add_parser("run", help="run the engine on one config")
    run.add_argument("config", help="a PTM3D YAML config")
    run.add_argument("--exe", default=None,
                     help="engine binary (default: found via PTM_BUILD_DIR, "
                          "then <repo>/build/bin)")
    run.add_argument("--workdir", default=None,
                     help="working directory (default: the repository root, "
                          "which is where the shipped metakernels resolve)")

    advise = sub.add_parser(
        "advise",
        help="say whether lateral heat conduction matters for a config's mesh",
        description="Compare the mesh's facet sizes with the lateral-conduction "
                    "length D*(T) = 8k/(pi h) of the configured surface material "
                    "and say whether thermal.lateral_conduction is recommended. "
                    "Reads the config and the mesh; never runs the engine. The "
                    "engine prints the same advisory at the start of a run.")
    advise.add_argument("config", help="a PTM3D YAML config")
    advise.add_argument("--mesh-units", default=None, choices=("m", "km"),
                        help="override the mesh length unit (default: inferred as "
                             "the engine infers it)")
    advise.add_argument("--json", action="store_true",
                        help="emit the numbers and verdict as JSON")

    selfcheck = sub.add_parser(
        "gpu-selfcheck",
        help="check the GPU shadow-ray tracer on this machine's GPU(s)",
        description="Run build/bin/ptm_gpu_selfcheck: trace shadow rays on each "
                    "OpenCL GPU over two synthetic crater meshes and compare every "
                    "answer with a float64 brute-force ground truth. GPU "
                    "illumination (illumination.use_gpu) was validated on an AMD "
                    "Radeon Pro W6800X; run this before using it on other "
                    "hardware. Takes seconds; needs no data files.",
        epilog="exit status: 0 PASS, 1 FAIL, 2 SKIP (no GPU build or device, or "
               "the self-check binary is missing), 3 INVALID (the check could not "
               "be trusted)")
    selfcheck.add_argument("--exe", default=None,
                           help="self-check binary (default: found via "
                                "PTM_BUILD_DIR, then <repo>/build/bin)")
    selfcheck.add_argument("--device", default="all",
                           help="'all' (default) or one OpenCL GPU index")
    selfcheck.add_argument("--pause-ms", type=int, default=None,
                           help="host sleep between dispatches (default: the "
                                "binary's, 500)")

    cadence = sub.add_parser(
        "cadence",
        help="post-run output-cadence indicator for a finished run",
        description="Area of faces lit in exactly one output frame, overall and among faces with "
                    "frame-sampled T_max < 110 K and < 54 K: the post-run half of the engine's "
                    "[cadence-advisor] block. The total bounds the area whose sunlit events the "
                    "frames missed entirely; the cold-trap split is an indicator, not a bound.")
    cadence.add_argument("run_dir", help="a finished PTM3D output directory with direct_solar/ frames")
    cadence.add_argument("--json", action="store_true", help="emit the numbers as JSON")

    sweep = sub.add_parser(
        "sweep", help="run a grid of simulations from one base config",
        description="Run a parameter grid. Serial by default: there is one "
                    "GPU. Resumable: re-running skips points that finished. "
                    "Reading a finished sweep needs no engine.",
        epilog="thermal inertia is TWO parameters, and they are different physics:\n"
               "  thermal_inertia            ks = TI^2/(rhos*cp(273 K)), kd scaled\n"
               "                             with it; depth contrast preserved.\n"
               "                             Needs an engine with --thermal-inertia.\n"
               "  constant_thermal_inertia   cp = TI^2/(ks*rhos), chi=0, kd=ks,\n"
               "                             rhod=rhos; TI exactly constant in\n"
               "                             temperature AND depth. Any engine.\n"
               "  They give measurably different diurnal amplitudes on the same\n"
               "  config. Pick one deliberately; you cannot sweep both at once,\n"
               "  and neither is a default.\n"
               "\n"
               "examples:\n"
               "  ptm sweep examples/quickstart_asteroid_diurnal.yaml \\\n"
               "      --vary thermal_inertia=50,100,200,400\n"
               "  ptm sweep output/asteroid_sweep --no-run --table sweep.csv "
               "--plot curves --save family.png\n",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sweep.add_argument(
        "target",
        help="a base YAML config, or an existing sweep directory to reopen")
    sweep.add_argument(
        "--vary", action="append", default=[], metavar="NAME=V1,V2,...",
        help="one grid axis. NAME is 'thermal_inertia', "
             "'constant_thermal_inertia', 'albedo', or any dotted config key "
             "(thermal.equilibration_orbits). Repeat for a cartesian product. "
             "See the note below on the two thermal-inertia conventions.")
    sweep.add_argument("--outdir", default=None,
                       help="where the points go (default: the base config's "
                            "output.directory with _sweep appended)")
    sweep.add_argument("--exe", default=None, help="engine binary")
    sweep.add_argument("--workdir", default=None,
                       help="working directory for the engine (default: the "
                            "repository root)")
    sweep.add_argument("--jobs", type=int, default=1,
                       help="concurrent engine processes (default 1). CPU-only: "
                            "refuses to start if any point asks for the GPU.")
    sweep.add_argument("--timeout", type=float, default=None,
                       help="wall-clock ceiling per point, in seconds")
    sweep.add_argument("--no-resume", action="store_true",
                       help="re-run every point, including finished ones")
    sweep.add_argument("--no-run", action="store_true",
                       help="do not launch anything; only --table / --plot the "
                            "points already on disk")
    sweep.add_argument("--dry-run", action="store_true",
                       help="print the engine command for each point, run none")
    sweep.add_argument("--list", action="store_true",
                       help="print the grid and each point's status, run none")
    sweep.add_argument("--use-flags", default="auto",
                       choices=("auto", "yes", "no"),
                       help="deliver the named parameters through the engine "
                            "CLI flags (auto: use them if this binary has "
                            "them). 'no' cannot serve thermal_inertia, which "
                            "has no config-file equivalent.")
    sweep.add_argument("--table", default=None, metavar="FILE.csv",
                       help="write a tidy table of the finished points")
    sweep.add_argument("--plot", default=None, choices=("curves",),
                       help="draw the family-of-curves figure")
    sweep.add_argument("--field", default="temperature",
                       help="field for --table and --plot (default temperature)")
    sweep.add_argument("--save", default=None, metavar="FILE",
                       help="write the figure here")
    sweep.add_argument("--dpi", type=int, default=200, help="raster DPI (200)")

    obs = sub.add_parser(
        "observe",
        help="disk-integrated thermal emission seen by an observer",
        description="Flux density and disk brightness temperature of a run's "
                    "body as seen by a distant observer at given times. With "
                    "--observer, each time is when the observer receives the "
                    "light. It is corrected for light time, and the "
                    "temperatures are interpolated in time between the two "
                    "output steps around the emission epoch. With "
                    "--sub-observer, the geometry is given directly and "
                    "--step chooses the output step. Give --emissivity for "
                    "the surface mode, or the electrical properties for the "
                    "subsurface mode. With --band-um or --response the "
                    "emission is integrated over a thermal-infrared band, "
                    "with the sunlight the surface reflects in that band. "
                    "See docs/observing.md.",
        epilog="examples:\n"
               "  ptm observe output/quickstart_europa_diurnal \\\n"
               "      --time 2025-01-01T12:00:00 --freq-ghz 233 --emissivity 0.9\n"
               "  ptm observe output/quickstart_europa_diurnal --step 40 \\\n"
               "      --sub-observer 0 90 6.3e8 --freq-ghz 233 --emissivity 0.9\n"
               "  ptm observe output/quickstart_europa_observing \\\n"
               "      --time 2025-01-01T12:00:00 --freq-ghz 233 --porosity 0.5\n"
               "  ptm observe output/quickstart_europa_diurnal \\\n"
               "      --time 2025-01-01T12:00:00 --band-um 8 14 --emissivity 0.9\n",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    _add_common(obs)
    obs.add_argument("--time", action="append", default=[], metavar="UTC",
                     help="observation (receive) time, repeat for several")
    obs.add_argument("--observer", default="EARTH",
                     help="SPICE name or NAIF id of the observer (EARTH)")
    obs.add_argument("--kernel", action="append", default=[], metavar="FILE",
                     help="extra SPICE kernel such as a spacecraft SPK, repeatable")
    obs.add_argument("--sub-observer", nargs=3, type=float, default=None,
                     metavar=("LAT", "LON", "DIST_KM"),
                     help="explicit geometry instead of SPICE: planetocentric "
                          "latitude, east longitude (deg), distance from the "
                          "body center (km)")
    obs.add_argument("--step", type=int, action="append", default=[],
                     help="output step(s) for --sub-observer")
    obs.add_argument("--nearest-step", action="store_true",
                     help="move each observation by up to half an output step "
                          "so its emission epoch falls on a step, instead of "
                          "interpolating the temperatures in time (default)")
    band = obs.add_mutually_exclusive_group(required=True)
    band.add_argument("--freq-ghz", type=float, help="observing frequency, GHz")
    band.add_argument("--wavelength-um", type=float,
                      help="observing wavelength, micrometers")
    band.add_argument("--band-um", nargs=2, type=float, default=None,
                      metavar=("LO", "HI"),
                      help="thermal-infrared band with a response of 1 from LO "
                           "to HI micrometers")
    band.add_argument("--response", default=None, metavar="FILE",
                      help="thermal-infrared band with the relative spectral "
                           "response in FILE (two columns: wavelength in "
                           "micrometers, response)")
    obs.add_argument("--emissivity", type=float, default=None,
                     help="surface mode: observing-band emissivity (not the "
                          "thermal-infrared emissivity of the energy balance). "
                          "With a band: the gray emissivity, whose complement "
                          "is the in-band reflectance")
    tir = obs.add_argument_group("thermal-infrared band (--band-um or --response)")
    tir.add_argument("--no-reflected", action="store_true",
                     help="leave out the sunlight reflected in the band")
    tir.add_argument("--faces-csv", default=None, metavar="FILE",
                     help="also write every visible face at every time: "
                          "latitude, longitude, emission cosine, temperature, "
                          "thermal and reflected band radiance, band "
                          "brightness temperature")
    sub = obs.add_argument_group(
        "subsurface mode (needs a run with output.map_depths_per_step)")
    sub.add_argument("--skin-depth-m", type=float, default=None,
                     help="electrical skin depth L_e at nadir, m (with "
                          "--refractive-index, or replacing the porosity "
                          "model's absorption)")
    sub.add_argument("--refractive-index", type=float, default=None,
                     help="real refractive index of the near surface")
    sub.add_argument("--porosity", type=float, default=None,
                     help="porous-ice dielectric model (scripts/emission)")
    sub.add_argument("--reference-temperature-k", default="100",
                     help="temperature at which the porosity model is "
                          "evaluated (constant with depth), or 'node' to use "
                          "each node's temperature (default 100)")
    obs.add_argument("--csv", default=None, metavar="FILE",
                     help="also write the table as CSV")
    return parser


# --- observe ----------------------------------------------------------------

_OBS_COLUMNS = (
    ("utc_emission", "emitted UTC (at body)", "{}"),
    ("steps", "steps", "{}"),
    ("rotation_offset_deg", "shift deg", "{:.3f}"),
    ("sub_observer_lat_deg", "obs lat", "{:.3f}"),
    ("sub_observer_lon_deg", "obs lon", "{:.3f}"),
    ("phase_angle_deg", "phase deg", "{:.2f}"),
    ("distance_km", "dist km", "{:.4e}"),
    ("flux_density_jy", "F_nu Jy", "{:.5g}"),
    ("tb_planck_k", "Tb Planck K", "{:.3f}"),
    ("tb_rayleigh_jeans_k", "Tb RJ K", "{:.3f}"),
    ("occultation", "view", "{}"),
)


_BAND_COLUMNS = (
    ("utc_emission", "emitted UTC (at body)", "{}"),
    ("steps", "steps", "{}"),
    ("rotation_offset_deg", "shift deg", "{:.3f}"),
    ("sub_observer_lat_deg", "obs lat", "{:.3f}"),
    ("sub_observer_lon_deg", "obs lon", "{:.3f}"),
    ("phase_angle_deg", "phase deg", "{:.2f}"),
    ("distance_km", "dist km", "{:.4e}"),
    ("flux_w_m2", "F W/m2", "{:.5g}"),
    ("tb_k", "Tb K", "{:.3f}"),
    ("tb_thermal_k", "Tb thermal K", "{:.3f}"),
    ("mean_temperature_k", "mean T K", "{:.3f}"),
    ("reflected_fraction_disk", "refl/therm", "{:.4f}"),
    ("reflected_fraction_max_face", "max face", "{:.4f}"),
    ("occultation", "view", "{}"),
)


def _do_observe(args, run):
    import csv
    import warnings
    from . import observe

    if args.sub_observer is None and not args.time:
        raise SystemExit("ptm observe: give --time (with --observer) or "
                         "--sub-observer with --step")
    if args.sub_observer is not None and args.time:
        raise SystemExit("ptm observe: --time goes with --observer. With "
                         "--sub-observer choose output steps with --step")
    if args.sub_observer is not None and not args.step:
        raise SystemExit("ptm observe: --sub-observer needs --step")

    geoms = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            if args.sub_observer is not None:
                lat, lon, dist = args.sub_observer
                geoms = [observe.explicit_geometry(lat, lon, dist, run=run, step=s)
                         for s in args.step]
            else:
                match = "nearest" if args.nearest_step else "interpolate"
                geoms = [observe.observer_geometry(run, t, observer=args.observer,
                                                   kernels=args.kernel,
                                                   time_match=match)
                         for t in args.time]
        except (observe.ObservationTimeError, FileNotFoundError, ValueError) as exc:
            raise SystemExit(f"ptm observe: {exc}")
    for w in caught:
        print(f"WARNING: {w.message}", file=sys.stderr)

    if args.band_um is not None or args.response is not None:
        return _do_observe_band(args, run, geoms)
    if args.no_reflected or args.faces_csv:
        raise SystemExit("ptm observe: --no-reflected and --faces-csv go with "
                         "a thermal-infrared band (--band-um or --response)")

    subsurface = None
    if (args.skin_depth_m is not None or args.porosity is not None
            or args.refractive_index is not None):
        if args.emissivity is not None:
            raise SystemExit("ptm observe: the subsurface mode takes no "
                             "--emissivity (the surface transmissivity follows "
                             "from the refractive index)")
        tref = (None if str(args.reference_temperature_k).lower() == "node"
                else float(args.reference_temperature_k))
        subsurface = observe.Subsurface(
            skin_depth_m=args.skin_depth_m, refractive_index=args.refractive_index,
            porosity=args.porosity, reference_temperature_k=tref)
    elif args.emissivity is None:
        raise SystemExit("ptm observe: give --emissivity (surface mode) or the "
                         "subsurface options (--skin-depth-m with "
                         "--refractive-index, or --porosity)")

    rows = []
    for g in geoms:
        try:
            if subsurface is None:
                em = observe.disk_emission(run, observe.temperature_at(run, g), g,
                                           emissivity=args.emissivity,
                                           frequency_ghz=args.freq_ghz,
                                           wavelength_um=args.wavelength_um)
            else:
                depths, T = observe.depth_temperatures_at(run, g)
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    em = observe.disk_emission(run, T, g, subsurface=subsurface,
                                               depths_m=depths,
                                               frequency_ghz=args.freq_ghz,
                                               wavelength_um=args.wavelength_um)
                for w in caught:
                    print(f"WARNING: {w.message}", file=sys.stderr)
        except (observe.NonConvexMeshError, FileNotFoundError, ValueError) as exc:
            raise SystemExit(f"ptm observe: {exc}")
        rows.append(em.as_row())

    em0 = em
    if subsurface is None:
        print(f"{run.path.name}: surface mode, {em0.frequency_hz / 1e9:.3f} GHz, "
              f"band emissivity {args.emissivity}, radius {em0.radius_km:.1f} km")
    else:
        beta = ("unknown (no thermal skin depth recorded)"
                if em0.beta_nadir is None else f"{em0.beta_nadir:.2f}")
        print(f"{run.path.name}: subsurface mode, {em0.frequency_hz / 1e9:.3f} GHz, "
              f"electrical skin depth {100 * em0.skin_depth_m:.2f} cm, "
              f"refractive index {em0.refractive_index:.3f}, "
              f"radius {em0.radius_km:.1f} km")
        print(f"beta = L_e / thermal skin depth = {beta}. Much less than 1, "
              f"the surface mode would serve. Near 1 or more, the emission "
              f"samples the diurnal temperature wave.")
    if args.sub_observer is None:
        if args.nearest_step:
            print(f"observer {args.observer}, times corrected for light time "
                  f"(CN+S) and moved to the nearest output step")
        else:
            print(f"observer {args.observer}, times corrected for light time "
                  f"(CN+S), temperatures interpolated linearly in time "
                  f"between the output steps shown")
    print("brightness temperatures use the solid angle pi R^2/D^2")
    head = [h for _, h, _ in _OBS_COLUMNS]
    if args.sub_observer is None:
        head = ["observed UTC"] + head
    print("  ".join(head))
    for r in rows:
        cells = [fmt.format(r[k]) if r[k] is not None else "-"
                 for k, _, fmt in _OBS_COLUMNS]
        if args.sub_observer is None:
            cells = [r["utc_observed"]] + cells
        print("  ".join(cells))
    if args.csv:
        with open(args.csv, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"wrote {args.csv}")
    return 0


def _do_observe_band(args, run, geoms):
    import csv
    import numpy as np
    from . import observe

    if (args.skin_depth_m is not None or args.porosity is not None
            or args.refractive_index is not None):
        raise SystemExit("ptm observe: the subsurface options go with "
                         "--freq-ghz or --wavelength-um, not with a band")
    if args.emissivity is None:
        raise SystemExit("ptm observe: a band needs --emissivity (the gray "
                         "emissivity of the surface)")
    try:
        band = (observe.Band.top_hat(*args.band_um) if args.band_um is not None
                else observe.Band.from_file(args.response))
    except (OSError, ValueError) as exc:
        raise SystemExit(f"ptm observe: {exc}")
    reflected = not args.no_reflected
    if reflected and "direct_solar" not in run.step_fields:
        raise SystemExit(
            f"ptm observe: {run.path} has no direct_solar frames, which the "
            f"reflected sunlight needs. Use --no-reflected to leave it out.")

    rows, faces = [], []
    lat, lon = run.lat, run.lon
    for k, g in enumerate(geoms):
        try:
            T = observe.temperature_at(run, g)
            E = observe.temperature_at(run, g, "direct_solar") if reflected else None
            em = observe.band_emission(run, T, g, band=band,
                                       emissivity=args.emissivity,
                                       direct_solar=E, reflected=reflected)
        except (observe.NonConvexMeshError, FileNotFoundError, ValueError) as exc:
            raise SystemExit(f"ptm observe: {exc}")
        rows.append(em.as_row())
        if args.faces_csv:
            vis = np.nonzero(em.face_mu > 0)[0]
            for j in vis:
                faces.append((k, rows[-1]["utc_observed"] or rows[-1]["step_utc"],
                              int(j), lat[j], lon[j], em.face_mu[j], T[j],
                              em.face_radiance_thermal[j],
                              em.face_radiance_reflected[j], em.face_tb_k[j]))

    em0 = em
    print(f"{run.path.name}: band {band.label}, gray emissivity "
          f"{args.emissivity}, reflected sunlight "
          f"{'included' if reflected else 'left out'}, radius "
          f"{em0.radius_km:.1f} km")
    print(f"fraction of a 5772 K blackbody Sun in the band: "
          f"{em0.solar_fraction:.4e}")
    if args.sub_observer is None:
        if args.nearest_step:
            print(f"observer {args.observer}, times corrected for light time "
                  f"(CN+S) and moved to the nearest output step")
        else:
            print(f"observer {args.observer}, times corrected for light time "
                  f"(CN+S), temperatures interpolated linearly in time "
                  f"between the output steps shown")
    print("band brightness temperatures use the solid angle of the visible "
          "faces. refl/therm: reflected over thermal flux of the disk. max "
          "face: the largest ratio on any visible face")
    head = [h for _, h, _ in _BAND_COLUMNS]
    if args.sub_observer is None:
        head = ["observed UTC"] + head
    print("  ".join(head))
    for r in rows:
        cells = [fmt.format(r[k]) if r[k] is not None else "-"
                 for k, _, fmt in _BAND_COLUMNS]
        if args.sub_observer is None:
            cells = [r["utc_observed"]] + cells
        print("  ".join(cells))
    if args.csv:
        with open(args.csv, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"wrote {args.csv}")
    if args.faces_csv:
        with open(args.faces_csv, "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["time_index", "utc", "face", "lat_deg", "lon_deg",
                             "mu", "temperature_k", "radiance_thermal_w_m2_sr",
                             "radiance_reflected_w_m2_sr", "tb_k"])
            for row in faces:
                writer.writerow([row[0], row[1], row[2], f"{row[3]:.4f}",
                                 f"{row[4]:.4f}", f"{row[5]:.6f}",
                                 f"{row[6]:.4f}", f"{row[7]:.6e}",
                                 f"{row[8]:.6e}", f"{row[9]:.4f}"])
        print(f"wrote {args.faces_csv} ({len(faces)} rows)")
    return 0


# --- info -------------------------------------------------------------------

def _print_info(run):
    prov = run.provenance()
    print(f"run          {prov['run_dir']}")
    print(f"faces        {prov['n_faces']}")

    time = prov["time"]
    if time:
        cadence = ("" if time["cadence_hours"] is None
                   else f", {time['cadence_hours']:.3g} h cadence")
        print(f"steps        {time['n_steps']} over "
              f"{time['elapsed_hours']:.4g} h "
              f"({time['elapsed_hours'] / 24.0:.4g} d){cadence}")
        print(f"window       {time['first_utc']} .. {time['last_utc']}")
    meta = prov["metadata"]
    if meta.get("failed_timesteps"):
        print(f"WARNING      {meta['failed_timesteps']} timestep(s) failed to "
              f"export")

    unit = prov["mesh_units"]
    corroborated = prov["mesh_units_corroborated_by"]
    how = (f"consistent with the {corroborated} radius under this unit"
           if corroborated
           else "from face-center magnitude alone -- NOT corroborated")
    print(f"mesh units   {unit} ({how})")
    print(f"mesh radius  {prov['mesh_median_radius_m'] / 1e3:.4g} km (median "
          f"|face center|)")

    print(f"per-step     {', '.join(prov['step_fields']) or '(none)'}")
    stats = prov["statistics"]
    if stats:
        grouped = {}
        for field, moment in stats:
            grouped.setdefault(field, []).append(moment or "(scalar)")
        print("statistics   " + "; ".join(
            f"{f} [{', '.join(m)}]" for f, m in sorted(grouped.items())))
    else:
        print("statistics   (none)")

    print(f"units style  {prov['units_style']} "
          f"(flux spelled {run.units_of('illumination')})")
    if prov["config_file"]:
        print(f"config       {prov['config_file']}")
        print(f"  sha256     {prov['config_sha256']}")
    engine = prov["engine"]
    if engine:
        # The point of provenance.txt is that these are answerable at all; a
        # run that has them should not make the reader open the file.
        print(f"engine       {prov['engine_version']} @ "
              f"{prov['engine_git_sha']} "
              f"({engine.get('engine_build_type', '?')}, GPU "
              f"{'on' if engine.get('engine_gpu_enabled') == 'true' else 'off'})")
        print(f"  binary     {engine.get('engine_binary', '?')}")
        print(f"  md5        {prov['engine_binary_md5']}")
        if prov["command_line"]:
            print(f"  invoked    {prov['command_line']}")
    if prov["started_utc"]:
        source = "provenance.txt" if engine.get("run_started_utc") else "heartbeat.txt"
        print(f"started      {prov['started_utc']} ({source})")

    print("\nnot recorded by the run:")
    for item in prov["not_recorded"]:
        print(f"  - {item}")


# --- plot -------------------------------------------------------------------

def _do_plot(args, run):
    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    from . import plots

    # save= is deliberately not forwarded: the writing happens once, below, so
    # that --dpi applies to the file that is actually kept.
    common = dict(show=args.show)
    if args.what == "diurnal":
        if args.face is None and (args.lat is None or args.lon is None):
            raise SystemExit(
                "ptm plot diurnal needs --face, or both --lat and --lon")
        fig = plots.plot_diurnal(run, field=args.field, face=args.face,
                                 lat=args.lat, lon=args.lon, **common)
    elif args.what == "map":
        fig = plots.plot_map(run, args.field, moment=args.moment,
                             step=args.step, resolution=args.resolution,
                             **common)
    elif args.what == "globes":
        fig = plots.plot_globes(run, args.field, moment=args.moment,
                                step=args.step,
                                mark_subsolar=args.step is not None, **common)
    else:
        fig = plots.plot_histogram(run, args.field, moment=args.moment,
                                   step=args.step, **common)

    if args.save:
        fig.savefig(args.save, dpi=args.dpi)
        print(f"wrote {args.save}")
    elif not args.show:
        print("nothing written: pass --save FILE or --show", file=sys.stderr)
        return 1
    return 0


# --- run --------------------------------------------------------------------

def _do_run(args):
    from . import paths, run as run_mod

    config = Path(args.config).expanduser().resolve()
    if not config.is_file():
        raise SystemExit(f"config not found: {config}")

    root = paths.find_repo_root(config.parent) or paths.find_repo_root()
    executable = (Path(args.exe).expanduser() if args.exe
                  else paths.find_engine_executable(
                      search_roots=[root] if root else []))
    if executable is None or not Path(executable).exists():
        raise SystemExit(
            "ptm_thermal_sim not found. Build it, or set PTM_BUILD_DIR, or "
            "pass --exe. See INSTALL.md.")

    workdir = args.workdir or (str(root) if root else None)
    config_dict = {}
    try:
        import yaml
        config_dict = yaml.safe_load(config.read_text()) or {}
    except Exception as exc:                                   # noqa: BLE001
        print(f"[ptm] could not parse {config} for CLI flags ({exc}); "
              f"passing it to the engine unread", file=sys.stderr)

    command = run_mod.build_command(executable, config, config_dict)
    print(f"[ptm] {' '.join(command)}")
    if workdir:
        print(f"[ptm] working directory: {workdir}")
    for msg in run_mod.long_run_advisory(config_dict):
        print(msg, file=sys.stderr)

    process = run_mod.EngineProcess(
        command, working_dir=workdir,
        on_output=lambda line: print(line),
        on_warning=lambda line: print(f"WARNING {line}", file=sys.stderr),
        on_error=lambda line: print(f"ERROR   {line}", file=sys.stderr),
        on_heartbeat=lambda msg, stalled: print(f"[ptm] {msg}", file=sys.stderr),
    )
    status = process.run()
    print(f"[ptm] engine exited {status}")
    return status


# --- advise -------------------------------------------------------------------

def _do_advise(args):
    from .lateral_advisor import advise_config

    config = Path(args.config).expanduser()
    if not config.is_file():
        raise SystemExit(f"config not found: {config}")
    from . import paths
    root = paths.find_repo_root(config.resolve().parent) or paths.find_repo_root()
    try:
        advice = advise_config(config, mesh_units=args.mesh_units, repo_root=root)
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(f"ptm advise: {exc}")
    if args.json:
        import json
        print(json.dumps(advice.as_dict(), indent=2))
    else:
        print("\n".join(advice.lines))
    return 0


# --- gpu-selfcheck -------------------------------------------------------------

def _do_gpu_selfcheck(args):
    import subprocess
    from . import paths

    root = paths.find_repo_root()
    executable = (Path(args.exe).expanduser() if args.exe
                  else paths.find_engine_executable(
                      "ptm_gpu_selfcheck", search_roots=[root] if root else []))
    if executable is None or not Path(executable).exists():
        print("ptm gpu-selfcheck: SKIP -- ptm_gpu_selfcheck not found. It is built "
              "only in a GPU build (-DPTM_ENABLE_GPU=ON -DGPU_BACKEND=OPENCL, target "
              "ptm_gpu_selfcheck); set PTM_BUILD_DIR or pass --exe. Nothing was "
              "checked.", file=sys.stderr)
        return 2
    command = [str(executable), "--device", str(args.device)]
    if args.pause_ms is not None:
        command += ["--pause-ms", str(args.pause_ms)]
    print(f"[ptm] {' '.join(command)}", file=sys.stderr)
    return subprocess.run(command).returncode


# --- cadence ------------------------------------------------------------------

def _do_cadence(args):
    from .cadence_advisor import format_postrun, postrun
    try:
        r = postrun(args.run_dir)
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(f"ptm cadence: {exc}")
    if args.json:
        import json
        print(json.dumps(r, indent=2))
    else:
        print(format_postrun(r))
    return 0


# --- sweep --------------------------------------------------------------------

def _parse_scalar(text):
    """``'50'`` -> 50, ``'0.05'`` -> 0.05, ``'true'`` -> True, else the string."""
    lowered = text.strip().lower()
    if lowered in ("true", "false"):
        return lowered == "true"
    for cast in (int, float):
        try:
            return cast(text)
        except ValueError:
            pass
    return text


def _parse_vary(specs):
    """``['ti=1,2', ...]`` -> an ordered ``{name: [values]}``."""
    grid = {}
    for spec in specs:
        if "=" not in spec:
            raise SystemExit(
                f"--vary needs NAME=V1,V2,...; got {spec!r}")
        name, _, values = spec.partition("=")
        name = name.strip()
        items = [v for v in values.split(",") if v.strip() != ""]
        if not items:
            raise SystemExit(f"--vary {name}: no values given")
        if name in grid:
            raise SystemExit(f"--vary {name} given twice")
        grid[name] = [_parse_scalar(v) for v in items]
    return grid


def _print_sweep_list(sweep):
    print(f"base   {sweep.base_config_path}")
    print(f"outdir {sweep.outdir}")
    print(f"engine {sweep.exe} (md5 {sweep.exe_md5})")
    for name, values in sweep.axes.items():
        route = sweep.route.get(name, "config")
        print(f"axis   {name} = {values}  [{route}]")
    print(f"points {len(sweep)} ({'x'.join(str(n) for n in sweep.shape)})")
    for point in sweep:
        complete, why = point.completeness()
        mark = "done" if complete else "todo"
        print(f"  {mark}  {point.name:<40s} {why}")


def _do_sweep(args):
    from . import sweep as sweep_mod

    target = Path(args.target).expanduser()
    if target.is_dir():
        if args.vary:
            raise SystemExit(
                f"{target} is an existing sweep directory; its grid comes from "
                f"sweep.json, so --vary cannot be combined with it")
        sweep = sweep_mod.load_sweep(target, exe=args.exe) if args.exe \
            else sweep_mod.load_sweep(target)
    else:
        if not target.is_file():
            raise SystemExit(f"not a config file or sweep directory: {target}")
        if not args.vary:
            raise SystemExit("ptm sweep needs at least one --vary NAME=V1,V2")
        use_flags = {"auto": "auto", "yes": True, "no": False}[args.use_flags]
        try:
            sweep = sweep_mod.Sweep(
                target, grid=_parse_vary(args.vary), outdir=args.outdir,
                exe=args.exe, workdir=args.workdir, use_flags=use_flags,
                timeout=args.timeout)
        except ValueError as exc:
            raise SystemExit(f"ptm sweep: {exc}")

    if args.list:
        try:
            _print_sweep_list(sweep)
        except (FileNotFoundError, RuntimeError) as exc:
            raise SystemExit(f"ptm sweep: {exc}")
        return 0

    status = 0
    if not args.no_run:
        # A missing engine and an axis this engine cannot deliver are both
        # ordinary operator errors with a paragraph of explanation attached;
        # they belong on stderr as a message, not as a traceback.
        try:
            summary = sweep.run(resume=not args.no_resume, jobs=args.jobs,
                                dry_run=args.dry_run)
        except (FileNotFoundError, RuntimeError, ValueError) as exc:
            raise SystemExit(f"ptm sweep: {exc}")
        if summary["failed"]:
            status = 1
        if args.dry_run:
            return status

    if args.table:
        import numpy as np
        field = args.field
        rows = sweep.write_table(
            args.table,
            **{f"{field}_max": lambda r: float(np.max(r.stat(field, "max"))),
               f"{field}_min": lambda r: float(np.min(r.stat(field, "min"))),
               f"{field}_mean": lambda r: float(np.mean(r.stat(field, "avg")))})
        print(f"wrote {args.table} ({len(rows)} rows)")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        fig = sweep.plot_curves(field=args.field)
        if not args.save:
            raise SystemExit("--plot needs --save FILE")
        fig.savefig(args.save, dpi=args.dpi)
        print(f"wrote {args.save}")
    return status


# --- entry point ------------------------------------------------------------

def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.version:
        from . import __version__
        print(__version__)
        return 0
    if args.command is None:
        parser.print_help()
        return 1

    if args.command == "run":
        return _do_run(args)
    if args.command == "sweep":
        return _do_sweep(args)
    if args.command == "advise":
        return _do_advise(args)
    if args.command == "cadence":
        return _do_cadence(args)
    if args.command == "gpu-selfcheck":
        return _do_gpu_selfcheck(args)

    from .results import Results
    units = getattr(args, "units", "ascii")
    mesh_units = getattr(args, "mesh_units", None)
    try:
        run = Results(args.run_dir, units_style=units, mesh_units=mesh_units)
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(f"ptm: {exc}")

    if args.command == "info":
        if args.json:
            import json
            print(json.dumps(run.provenance(), indent=2, default=str))
            return 0
        _print_info(run)
        return 0
    if args.command == "observe":
        return _do_observe(args, run)
    return _do_plot(args, run)


if __name__ == "__main__":
    raise SystemExit(main())
