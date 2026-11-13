#!/usr/bin/env python3
"""output.map_depths_per_step: per-step temperature frames at depth.

    python3 src/python/tests/test_depth_frames.py [path/to/ptm_thermal_sim]

Runs the Europa quickstart (CPU, about 15 s per run) with depth frames on and
checks them against products the engine already writes from the same column:

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
   checkpoint.

Needs the built engine and the Europa kernels (python3
scripts/fetch_kernels.py --bodies europa). Exit 2 (SKIP) if either is absent.
"""

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
    cfg = yaml.safe_load(open(REPO / "examples" / "quickstart_europa_diurnal.yaml"))
    cfg["output"]["directory"] = str(tmp / tag)
    thermal = output.pop("_thermal", None)
    cfg["output"].update(output)
    if thermal:
        cfg["thermal"].update(thermal)
    path = tmp / f"{tag}.yaml"
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    proc = subprocess.run([str(exe), "--config", str(path)], cwd=REPO,
                          capture_output=True, text=True)
    (tmp / f"{tag}.log").write_text(proc.stdout + proc.stderr)
    return proc.returncode, proc.stdout + proc.stderr, tmp / tag


def field_name(depth):
    return f"temperature_z{depth:.3f}m".replace(".", "p", 1)


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
        r = Results(out)
        steps = r.steps("temperature")
        n = len(steps)
        fields = [field_name(d) for d in DEPTHS]
        present = all(r.steps(f) == steps for f in fields)
        check("one frame per output step at every depth", present,
              f"{len(fields)} depths x {n} steps")
        z0 = fields[0]
        same = all(np.array_equal(r.field(z0, k), r.field("temperature", k))
                   for k in steps)
        check("z = 0 frame equals the surface temperature frame, every step", same)
        worst_avg, exact_extrema = 0.0, True
        for f in fields:
            stack = np.stack([r.field(f, k) for k in steps]).astype(np.float64)
            worst_avg = max(worst_avg, float(np.max(np.abs(
                stack.mean(axis=0) - r.stat(f, "avg")))))
            exact_extrema &= np.array_equal(stack.max(axis=0).astype(np.float32),
                                            r.stat(f, "max"))
            exact_extrema &= np.array_equal(stack.min(axis=0).astype(np.float32),
                                            r.stat(f, "min"))
        check("frames reduce to statistics max and min exactly", exact_extrema)
        check("frames average to statistics avg within float32 rounding",
              worst_avg < 1e-4, f"worst {worst_avg:.1e} K")
        worst_site = 0.0
        for name in SITES:
            face = None
            for line in log.splitlines():
                if f"site '{name}' -> face" in line:
                    face = int(line.rsplit(" ", 1)[1])
            series = r.depth_series(name)
            frames = np.array([[r.field(f, k)[face] for f in fields] for k in steps])
            if series.temperature.shape != frames.shape:
                check(f"site {name}: series and frames have the same shape", False,
                      f"{series.temperature.shape} vs {frames.shape}")
                continue
            worst_site = max(worst_site, float(np.max(np.abs(series.temperature - frames))))
        check("frames at each site's face equal its depth series",
              worst_site < 2e-5, f"worst {worst_site:.1e} K over "
              f"{len(SITES)} sites x {n} steps x {len(DEPTHS)} depths")

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

    print()
    if FAILS:
        print(f"FAILED: {FAILS}")
        return 1
    print("PASS: all checks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
