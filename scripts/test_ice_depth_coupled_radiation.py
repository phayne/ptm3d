#!/usr/bin/env python3
"""The ice-stability depth must be computed from the radiation-coupled temperatures.

``thermal.compute_ice_stability`` accumulates the time mean <E(z)> of the
sublimation rate at every node of every column, and solves the stability depth
from it. In a run with radiation exchange, every face that receives indirect
flux (scattered sunlight and thermal emission from the surrounding terrain) is
re-solved by the coupling re-thermal, and its reported temperatures come from
that re-solve. Until 2026-11 the accumulator was filled only by the whole-mesh
solve, before the coupling, so ice_depth.bin ignored the indirect flux on every
face (PTM-629).

The case is the public synthetic bowl crater `data/meshes/lunar_crater_80N.ply`
(1 km across, 200 m deep, at 80N), one lunation on the CPU. The shadowed floor
is warmed almost only by the crater walls. Two runs share every setting except
radiation exchange:

  off      radiation exchange off;
  coupled  radiation exchange on, coupling re-thermal of every face carrying
           indirect flux.

Both pin the spin-up (``seed_anchor: warm``, ``seed_indirect: off``) and the
coupled run sets ``eq_indirect_source: off``, so the whole-mesh solve is the
same integration in both runs and any difference comes from the coupling alone.
Both set ``checkpoint_lookahead_steps: 1`` so that every sample of the
accumulation window is a recorded frame of the floor site's depth series.

Gates:
  G1  the coupled temperatures differ from the uncoupled ones (the case is not
      vacuous): the floor site's mean surface temperature rises by > 10 K;
  G2  the coupled ice product differs from the uncoupled one, in the direction
      physics requires: indirect flux only adds energy, so no <E(z)> sum
      decreases (relative 1e-5), no stability depth becomes shallower, and at
      least 100 faces become deeper. With the defect the two products are
      byte-identical;
  G3  the accumulator describes the integration the temperatures come from: at
      the floor site (a re-solved face), <E(z_j)> from ice_depth_sum_E.bin
      equals the mean of E(T(z_j, t)) over the recorded depth series of the
      same window, at every node below the surface, to 1e-4 relative (float32
      sums, six-decimal CSV). With the defect this fails by orders of magnitude;
  G4  the same reconstruction holds at a sunlit site on the plain in both runs
      (the floor is below the 30 K cutoff of E with radiation off, so the
      reconstruction is checked where it compares nonzero numbers);
  G5  provenance.txt records which temperatures the product was built from;
  G6  the SURFACE node, which decides the surface class: <E(0)> from the
      accumulator equals the mean E of the recorded surface temperature to 1%
      at the floor (coupled) and the plain. Not 1e-4: the accumulator samples
      the column's solved surface node, the series the re-evaluated output
      surface (measured 0.1-0.2%, about 0.005 K); a lost, zeroed or
      wrong-window node 0 is off by orders of magnitude.

Needs the Moon SPICE kernels (python3 scripts/fetch_kernels.py --bodies moon)
and build/bin/ptm_thermal_sim. Standard library only. About 25 s with 4 threads.

  python3 scripts/test_ice_depth_coupled_radiation.py [--threads N] [--keep DIR]
"""
import argparse
import json
import math
import os
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
MESH = REPO / "data/meshes/lunar_crater_80N.ply"
STEPS = 118
LOCATIONS = "name,lat_deg,lon_deg\nfloor,80.0,0.0\nplain,80.023085,0.0\n"


def config(work, label, threads, radiation, depths_file=None):
    cfg = {
        "simulation": {"start_time": "2026-04-05T00:00:00", "duration_days": 29.53,
                       "timesteps": STEPS, "body": "Moon"},
        "mesh": {"import_file": str(MESH), "scale_factor": 1.0, "use_topographic": False},
        "thermal": {"run_thermal_mesh": True, "threads": threads, "equilibration_orbits": 1.0,
                    "equilibration_steps_per_day": -1.0, "use_gpu": False,
                    "compute_statistics": True, "enable_shadow_flux": False,
                    "crater_dd_ratio": 8.0, "geothermal_heat_flow": -1.0,
                    "compute_ice_stability": True, "seed_anchor": "warm",
                    "seed_indirect": "off", "eq_indirect_source": "off",
                    "checkpoint_lookahead_steps": 1,
                    "final_state_file": str(work / (label + "_final.ptmt"))},
        "thermophysical": {"use_default_albedo": True, "incidence_angle_albedo": True},
        "radiation": {"enabled": False},
        "illumination": {"include_primary_reflection": True, "use_batch_processing": True,
                         "use_albedo_map": False},
        "output": {"directory": str(work / label), "format": "binary",
                   "write_raw_stats": True, "locations_file": str(work / "locs.csv")},
        "advanced": {"verbose": False},
    }
    if depths_file:
        cfg["output"]["depths_file"] = str(depths_file)
    if radiation:
        cfg["radiation"] = {"enabled": True, "use_gpu": False, "use_hmatrix": True,
                            "use_aca": True, "aca_tolerance": 0.001, "hmatrix_threads": threads,
                            "view_factor_file": str(work / "crater80N.vf"),
                            "max_coupling_iterations": 5, "coupling_temp_threshold": 2.0,
                            "coupling_flux_threshold": 0.1, "coupling_face_fraction": 0.01,
                            "band_separated": True}
    path = work / (label + ".yaml")
    path.write_text(json.dumps(cfg, indent=1))           # JSON is valid YAML
    return path


def run(exe, cfg, log):
    with open(log, "w") as fh:
        r = subprocess.run([exe, "--config", str(cfg)], cwd=REPO, stdout=fh,
                           stderr=subprocess.STDOUT)
    if r.returncode != 0:
        sys.stdout.write(Path(log).read_text()[-4000:])
        raise SystemExit("[test_ice_depth_coupled_radiation] FAIL: engine exited %d (log %s; "
                         "the work directory is kept)" % (r.returncode, log))
    return Path(log).read_text()


def read_bin(path):
    b = Path(path).read_bytes()
    n = struct.unpack_from("<I", b, 0)[0]
    return struct.unpack_from("<%df" % n, b, 4)


def read_raw(rundir):
    raw = Path(rundir) / "statistics" / "_raw"
    b = (raw / "ice_depth_sum_E.bin").read_bytes()
    nl, nf = struct.unpack_from("<II", b, 0)
    sum_e = struct.unpack_from("<%df" % (nl * nf), b, 8)
    n_samples = int((raw / "ice_depth_n_samples.txt").read_text().split()[0])
    g = (raw / "depth_grid.bin").read_bytes()
    z = struct.unpack_from("<%dd" % struct.unpack_from("<I", g, 0)[0], g, 4)
    return nl, nf, sum_e, n_samples, z


def sublimation_e(t):
    """The engine's E(T) (iceAccumulateBatchColumn), in double precision."""
    if t <= 30.0:
        return 0.0
    return 6.21e23 * 611.657 * math.exp(-6139.06 * (1.0 / t - 1.0 / 273.16)) / math.sqrt(t)


def depth_series(rundir, site):
    rows = [ln for ln in (Path(rundir) / "locations" / (site + "_depths.csv")).read_text()
            .splitlines() if ln and not ln.startswith("#")]
    return [[float(x) for x in r.split(",")[1:]] for r in rows[1:]]


def site_face(log, site):
    for ln in log.splitlines():
        if ln.startswith("[location-columns] site '%s' -> face " % site):
            return int(ln.rsplit(" ", 1)[1])
    raise SystemExit("[test_ice_depth_coupled_radiation] FAIL: the log does not name the "
                     "%s site's face (the depth series was not recorded from the mesh "
                     "column)" % site)


def worst_reconstruction(rundir, site, face, n_eq):
    """Largest relative gap between <E(z_j)> from the accumulator and from the
    recorded series, over nodes j >= 1, and the number of nodes compared (a node
    where both are 0, i.e. colder than 30 K all window, compares nothing). Node 0
    is skipped: under the default column step the accumulator holds the solved
    surface node, the series the re-evaluated output surface."""
    nl, nf, sum_e, n_samples, _ = read_raw(rundir)
    series = depth_series(rundir, site)
    # Interval t ends at frame t + 1; the window is [n_eq, n_eq + n_samples).
    frames = series[n_eq + 1:n_eq + 1 + n_samples]
    if len(frames) != n_samples:
        raise SystemExit("[test_ice_depth_coupled_radiation] FAIL: %d recorded frames for %d "
                         "samples" % (len(frames), n_samples))
    worst, n_cmp = 0.0, 0
    for j in range(1, nl):
        want = sum_e[j * nf + face] / n_samples
        got = sum(sublimation_e(fr[j]) for fr in frames) / n_samples
        if want > 0.0 or got > 0.0:
            n_cmp += 1
            worst = max(worst, abs(got - want) / max(abs(want), abs(got)))
    return worst, n_cmp


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--engine", default=str(REPO / "build/bin/ptm_thermal_sim"))
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--keep", default=None,
                    help="work directory to keep (default: a temp dir, removed on success)")
    a = ap.parse_args()

    if not os.access(a.engine, os.X_OK):
        raise SystemExit("[test_ice_depth_coupled_radiation] FAIL: %s is not built" % a.engine)
    if not (REPO / "data/kernels/spk/de440s.bsp").exists():
        raise SystemExit("[test_ice_depth_coupled_radiation] FAIL: Moon SPICE kernels missing; "
                         "run python3 scripts/fetch_kernels.py --bodies moon")
    work = Path(a.keep).resolve() if a.keep else Path(tempfile.mkdtemp(prefix="ptm_ice_rad_"))
    work.mkdir(parents=True, exist_ok=True)
    (work / "locs.csv").write_text(LOCATIONS)
    os.environ["OMP_NUM_THREADS"] = str(a.threads)

    # A first radiation-off run without a depth list gives the node depths, so
    # that the depth series of the gated runs is recorded exactly at the nodes.
    run(a.engine, config(work, "grid", a.threads, False), work / "grid.log")
    _, _, _, _, z = read_raw(work / "grid")
    depths = work / "depths.txt"
    depths.write_text("".join(repr(v) + "\n" for v in z))

    log_off = run(a.engine, config(work, "off", a.threads, False, depths), work / "off.log")
    log_on = run(a.engine, config(work, "coupled", a.threads, True, depths), work / "coupled.log")
    face = site_face(log_on, "floor")
    face_plain = site_face(log_on, "plain")
    if site_face(log_off, "floor") != face or site_face(log_off, "plain") != face_plain:
        raise SystemExit("[test_ice_depth_coupled_radiation] FAIL: the sites map to different "
                         "faces in the two runs")
    n_eq = None
    for ln in log_on.splitlines():
        if ln.startswith("[Ice stability] CPU SIMD-batch accumulator armed: <E(z)> over output "
                         "intervals ["):
            n_eq = int(ln.split("[", 2)[2].split(",")[0])
    if n_eq is None:
        raise SystemExit("[test_ice_depth_coupled_radiation] FAIL: no accumulator window in the log")

    fails = []

    def gate(name, ok, detail):
        print("  %-3s %-74s %s" % (name, detail, "PASS" if ok else "FAIL"))
        if not ok:
            fails.append(name)

    t_off = read_bin(work / "off/statistics/temperature_avg.bin")
    t_on = read_bin(work / "coupled/statistics/temperature_avg.bin")
    gate("G1", t_on[face] - t_off[face] > 10.0,
         "floor face %d mean surface T: %.2f K off, %.2f K coupled"
         % (face, t_off[face], t_on[face]))

    z_off = read_bin(work / "off/statistics/ice_depth.bin")
    z_on = read_bin(work / "coupled/statistics/ice_depth.bin")
    _, _, e_off, _, _ = read_raw(work / "off")
    _, _, e_on, _, _ = read_raw(work / "coupled")
    n_lower = sum(1 for x, y in zip(e_off, e_on) if y < x * (1.0 - 1e-5))
    n_shallow = sum(1 for x, y in zip(z_off, z_on) if y < x)
    n_deep = sum(1 for x, y in zip(z_off, z_on) if y > x)
    n_surf_off = sum(1 for v in z_off if v == 0.0)
    n_surf_on = sum(1 for v in z_on if v == 0.0)
    gate("G2", n_lower == 0 and n_shallow == 0 and n_deep >= 100,
         "ice: %d faces deeper, %d shallower, %d <E> sums lower; surface ice %d -> %d"
         % (n_deep, n_shallow, n_lower, n_surf_off, n_surf_on))

    w_on, n_on = worst_reconstruction(work / "coupled", "floor", face, n_eq)
    gate("G3", w_on <= 1e-4 and n_on >= 10,
         "coupled floor: <E(z)> accumulator vs series, %d nodes, worst rel %.2e" % (n_on, w_on))
    w_off, n_off = worst_reconstruction(work / "off", "plain", face_plain, n_eq)
    w_pl, n_pl = worst_reconstruction(work / "coupled", "plain", face_plain, n_eq)
    gate("G4", max(w_off, w_pl) <= 1e-4 and min(n_off, n_pl) >= 10,
         "plain, off / coupled: %d / %d nodes, worst rel %.2e / %.2e"
         % (n_off, n_pl, w_off, w_pl))

    def provenance(label):
        for ln in (work / label / "provenance.txt").read_text().splitlines():
            if ln.startswith("ice_depth_temperatures="):
                return ln.split("=", 1)[1]
        return ""
    p_off, p_on = provenance("off"), provenance("coupled")
    gate("G5", p_off.startswith("whole_mesh_solve") and p_on.startswith("radiation_coupled"),
         "provenance: off '%s', coupled '%s'" % (p_off[:16], p_on[:17]))

    def node0_gap(rundir, site, face):
        nl, nf, sum_e, n_samples, _ = read_raw(rundir)
        frames = depth_series(rundir, site)[n_eq + 1:n_eq + 1 + n_samples]
        want = sum_e[face] / n_samples
        got = sum(sublimation_e(fr[0]) for fr in frames) / n_samples
        return abs(got - want) / max(abs(want), abs(got), 1e-300), want
    g_fl, e_fl = node0_gap(work / "coupled", "floor", face)
    g_pl, e_pl = node0_gap(work / "coupled", "plain", face_plain)
    gate("G6", g_fl <= 1e-2 and g_pl <= 1e-2 and e_fl > 0 and e_pl > 0,
         "surface node <E(0)>: floor rel %.1e, plain rel %.1e" % (g_fl, g_pl))

    if fails:
        print("\n[test_ice_depth_coupled_radiation] FAIL: %s (inputs, logs and products kept "
              "in %s)" % (", ".join(fails), work))
        return 1
    if not a.keep:
        shutil.rmtree(work, ignore_errors=True)
    print("\n[test_ice_depth_coupled_radiation] PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
