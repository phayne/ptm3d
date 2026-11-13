#!/usr/bin/env python3
"""output.map_depths_per_step: per-step temperature frames at depth.

    python3 src/python/tests/test_depth_frames.py [path/to/ptm_thermal_sim]

Runs the Europa quickstart and a short 80 N crater run (CPU) with depth
frames on and checks them against products the engine already writes from
the same column:

1. every requested depth has one frame per output step, readable with
   ``Results.field``;
2. the z = 0 frame equals the surface temperature frame bit for bit at every
   step;
3. the frames reduce to statistics/<field>_{avg,max,min}.bin (max and min
   exactly, the mean to float32 rounding);
4. at two named sites, the frames equal locations/<site>_depths.csv at every
   output step and depth, to the CSV's printed precision plus float32
   rounding;
5. the CPU scalar solver path writes complete frames too, with z = 0 equal to
   its own surface frames;
6. the run refuses map_depths_per_step without map_depths_m and with a batch
   checkpoint;
7. with no spin-up the window starts at step 0, whose frame is stored where
   the initial column is set rather than after a time step; checks 1-4 hold
   there too, including the first and last frames at the two sites;
8. the time-outer driver (used with lateral conduction) gives frames equal bit
   for bit to the default driver with lateral conduction off, and with
   lateral conduction on it writes complete frames whose z = 0 frame equals
   the surface frame and that reduce to the depth statistics. Check 8 uses
   the 80 N crater (a mesh in meters, radiation off) over two days.

Every solver site that stores a frame is reached: frames at a window step
can only come from the site that solves that step, and the engine refuses to
write frames with a missing (step, face). The default runs (1-4) store window
step 0 after a time step, check 7 stores it at the initial column, check 5
uses the scalar fold and check 8 the time-outer driver.

Needs the built engine and the Europa kernels (python3
scripts/fetch_kernels.py --bodies europa). Exit 2 (SKIP) if either is absent.
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import yaml

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "src" / "python"))

from ptm.results import Results  # noqa: E402

DEPTHS = [0.0, 0.004, 0.013, 0.03, 0.07, 0.15, 0.3]
SITES = {"eq_east": (0.0, 90.0), "mid_north": (45.0, -30.0)}
FAILS = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def run(exe, tmp, tag, **output):
    example = output.pop("_example", "quickstart_europa_diurnal.yaml")
    cfg = yaml.safe_load(open(REPO / "examples" / example))
    cfg["output"]["directory"] = str(tmp / tag)
    thermal = output.pop("_thermal", None)
    simulation = output.pop("_simulation", None)
    env = dict(os.environ, **output.pop("_env", {}))
    cfg["output"].update(output)
    if thermal:
        cfg["thermal"].update(thermal)
    if simulation:
        cfg["simulation"].update(simulation)
    path = tmp / f"{tag}.yaml"
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    proc = subprocess.run([str(exe), "--config", str(path)], cwd=REPO,
                          capture_output=True, text=True, env=env)
    (tmp / f"{tag}.log").write_text(proc.stdout + proc.stderr)
    return proc.returncode, proc.stdout + proc.stderr, tmp / tag


def field_name(depth):
    return f"temperature_z{depth:.3f}m".replace(".", "p", 1)


def check_run(r, log, fields, label, sites=True):
    """Checks 1-4 on one run: complete frames, z = 0 equal to the surface,
    reduction to the depth statistics, and (with sites) the site depth series
    at every step, reported separately for the first and last frames."""
    steps = r.steps("temperature")
    n = len(steps)
    present = all(r.steps(f) == steps for f in fields)
    check(f"{label}: one frame per output step at every depth", present,
          f"{len(fields)} depths x {n} steps")
    same = all(np.array_equal(r.field(fields[0], k), r.field("temperature", k))
               for k in steps)
    check(f"{label}: z = 0 frame equals the surface temperature frame, every step",
          same)
    worst_avg, exact_extrema = 0.0, True
    for f in fields:
        stack = np.stack([r.field(f, k) for k in steps]).astype(np.float64)
        worst_avg = max(worst_avg, float(np.max(np.abs(
            stack.mean(axis=0) - r.stat(f, "avg")))))
        exact_extrema &= np.array_equal(stack.max(axis=0).astype(np.float32),
                                        r.stat(f, "max"))
        exact_extrema &= np.array_equal(stack.min(axis=0).astype(np.float32),
                                        r.stat(f, "min"))
    check(f"{label}: frames reduce to statistics max and min exactly", exact_extrema)
    check(f"{label}: frames average to statistics avg within float32 rounding",
          worst_avg < 1e-4, f"worst {worst_avg:.1e} K")
    if not sites:
        return
    worst = {"all": 0.0, "first": 0.0, "last": 0.0}
    for name in SITES:
        face = None
        for line in log.splitlines():
            if f"site '{name}' -> face" in line:
                face = int(line.rsplit(" ", 1)[1])
        series = r.depth_series(name)
        frames = np.array([[r.field(f, k)[face] for f in fields] for k in steps])
        if series.temperature.shape != frames.shape:
            check(f"{label}: site {name}: series and frames have the same shape",
                  False, f"{series.temperature.shape} vs {frames.shape}")
            return
        diff = np.abs(series.temperature - frames)
        worst["all"] = max(worst["all"], float(np.max(diff)))
        worst["first"] = max(worst["first"], float(np.max(diff[0])))
        worst["last"] = max(worst["last"], float(np.max(diff[-1])))
    check(f"{label}: frames at each site's face equal its depth series",
          worst["all"] < 2e-5, f"worst {worst['all']:.1e} K over "
          f"{len(SITES)} sites x {n} steps x {len(fields)} depths")
    check(f"{label}: first and last frames equal the site series",
          max(worst["first"], worst["last"]) < 2e-5,
          f"first {worst['first']:.1e} K, last {worst['last']:.1e} K")


def main():
    exe = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "build" / "bin" / "ptm_thermal_sim"
    kernels = REPO / "data" / "kernels" / "spk" / "jup365.bsp"
    if not exe.is_file() or not kernels.exists():
        print(f"SKIP: needs {exe} and the Europa kernels")
        return 2

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        sites = tmp / "sites.csv"
        sites.write_text("name,lat,lon\n" + "".join(
            f"{n},{la},{lo}\n" for n, (la, lo) in SITES.items()))
        depths_file = tmp / "depths.txt"
        depths_file.write_text("\n".join(repr(d) for d in DEPTHS) + "\n")
        common = dict(map_depths_m=DEPTHS, map_depth_snap_to_node=False,
                      map_depths_per_step=True)

        print("1-4. CPU SIMD-batch path (the default)")
        rc, log, out = run(exe, tmp, "simd", locations_file=str(sites),
                           depths_file=str(depths_file), **common)
        check("run completes", rc == 0, f"rc {rc}")
        if rc != 0:
            print(log[-2000:])
            return 1
        fields = [field_name(d) for d in DEPTHS]
        z0 = fields[0]
        check_run(Results(out), log, fields, "default run")

        print("5. CPU scalar path")
        rc, log, out = run(exe, tmp, "scalar", _thermal={"use_simd_batch": False},
                           **common)
        check("scalar run completes", rc == 0, f"rc {rc}")
        if rc == 0:
            r = Results(out)
            ok = all(r.steps(f) == r.steps("temperature") for f in fields)
            ok &= all(np.array_equal(r.field(z0, k), r.field("temperature", k))
                      for k in r.steps("temperature"))
            check("scalar path: complete frames, z = 0 equals surface", ok)

        print("6. refusals")
        rc, log, _ = run(exe, tmp, "nodepths", map_depths_per_step=True)
        check("refused without map_depths_m",
              rc != 0 and "map_depths_per_step is true but output.map_depths_m is empty" in log,
              f"rc {rc}")
        rc, log, _ = run(exe, tmp, "ckpt", _thermal={"batch_checkpoint": True}, **common)
        check("refused with a batch checkpoint",
              rc != 0 and "map_depths_per_step is set with thermal.batch_checkpoint" in log,
              f"rc {rc}")

        print("7. no spin-up: the window starts at step 0 (Europa)")
        rc, log, out = run(exe, tmp, "nospinup", locations_file=str(sites),
                           depths_file=str(depths_file),
                           _thermal={"equilibration_orbits": 0.0}, **common)
        check("no spin-up run completes", rc == 0, f"rc {rc}")
        if rc == 0:
            check("the statistics window starts at step 0",
                  "abs steps [0," in log)
            check_run(Results(out), log, fields, "no spin-up")

        print("8. time-outer driver and lateral conduction (80 N crater, 2 days)")
        crater = dict(_example="tutorial04_crater80N_radoff.yaml",
                      _simulation={"duration_days": 2.0, "timesteps": 16}, **common)
        short = {"equilibration_orbits": 1.0}
        rc_b, log_b, out_b = run(exe, tmp, "crater_batch", _thermal=short, **crater)
        rc_t, log_t, out_t = run(exe, tmp, "crater_timeouter", _thermal=short,
                                 _env={"PTM_THERMAL_TIME_OUTER": "1"}, **crater)
        check("both crater runs complete", rc_b == 0 and rc_t == 0,
              f"rc {rc_b}, {rc_t}")
        if rc_b == 0 and rc_t == 0:
            check("the second run used the time-outer driver",
                  "TIME-OUTER driver" in log_t and "TIME-OUTER driver" not in log_b)
            b, t = Results(out_b), Results(out_t)
            same = all(np.array_equal(b.field(f, k), t.field(f, k))
                       for f in fields + ["temperature"]
                       for k in b.steps("temperature"))
            check("time-outer frames equal default-driver frames bit for bit",
                  same and b.steps("temperature") == t.steps("temperature"),
                  f"{len(fields)} depths x {len(b.steps('temperature'))} steps")
        # lateral_activation_c far above its default makes edges on this
        # meter-scale mesh active, so the coupled network is not empty.
        rc, log, out = run(exe, tmp, "crater_lateral",
                           _thermal=dict(short, lateral_conduction=True,
                                         lateral_activation_c=1.0e4), **crater)
        check("lateral conduction run completes", rc == 0, f"rc {rc}")
        if rc == 0:
            active = [line for line in log.splitlines() if "[lateral] activation:" in line]
            n_active = int(active[0].split("activation:")[1].split()[0]) if active else 0
            check("lateral edges are active", n_active > 0, f"{n_active} edges")
            check_run(Results(out), log, fields, "lateral", sites=False)
            if rc_b == 0:
                lat = Results(out)
                moved = max(float(np.max(np.abs(lat.field(f, k) - b.field(f, k))))
                            for f in fields for k in lat.steps("temperature"))
                check("lateral frames differ from the uncoupled run", moved > 0,
                      f"max {moved:.3f} K")

    print()
    if FAILS:
        print(f"FAILED: {FAILS}")
        return 1
    print("PASS: all checks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
