#!/usr/bin/env python3
"""The GPU illumination path must repeat the CPU path's lit/dark decisions (PTM-794 drift gate).

The GPU direct-solar path (`calculateDirectSolarBatchWithFiniteDiskGPU`) is a
hand-written mirror of the CPU (Embree) path: ray origins, near-distance
cutoff, center-ray fallback to disk sampling, disk pattern, rise-epsilon. Any
edit to one side that is not made on the other reopens PTM-794. This gate runs
the engine twice on a small core+ring case built from the shipped synthetic
crater (`data/meshes/lunar_crater_80N.ply`) plus a generated occluder ring that
rises away from the crater and shadows it at low Sun:

  CPU   illumination.use_gpu: false                      (CPU Embree)
  HOST  illumination.use_gpu: true, PTM_BVH_HOST_KERNEL=1, PTM_OPENCL_DEVICE_TYPE=cpu
        (the GPU path's host logic AND its kernel text, executed on the CPU;
        no GPU is used)

in two illumination environments (the validated fine-mesh recipe, and the
default), and requires:

  P1  0 lit/dark disagreements in direct_solar between CPU and HOST, every
      frame, every face, in both environments;
  P2  the case is not vacuous: lit face-frames exist, and the ring shadows
      core face-frames that a ring-less CPU run lights (control run NORING);
  P3  every run exits 0 and writes all its frames; HOST says the GPU tracer
      ran (not a silent CPU fallback), and provenance.txt records 0 fallbacks.

Needs a GPU (OpenCL) build of ptm_thermal_sim and the lunar SPICE kernels; it
never touches a GPU. A CPU-only build exits 2 (SKIP), loudly.

Usage: python3 scripts/test_gpu_host_parity.py [--binary build/bin/ptm_thermal_sim] [--keep DIR]
"""
import argparse
import math
import os
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CORE = REPO / "data/meshes/lunar_crater_80N.ply"
RECIPE = {"PTM_ILLUM_TRUE_NORMALS": "1", "PTM_ILLUM_RISE_EPS": "1", "PTM_ILLUM_RISE_KEPS": "0.2",
          "PTM_ILLUM_RISE_HABS": "0", "PTM_ILLUM_SEAL_HOLES": "1"}
FRAMES = 36


def read_core_vertices():
    lines = CORE.read_text().split("end_header\n", 1)[1].split("\n")
    head = CORE.read_text().split("end_header\n", 1)[0]
    nv = int(head.split("element vertex ")[1].split()[0])
    return [tuple(float(x) for x in lines[i].split()[:3]) for i in range(nv)]


def write_ring(path, verts):
    """A coarse annulus 1.5-8 km from the crater center that rises 450 m outward."""
    n = len(verts)
    c = [sum(v[k] for v in verts) / n for k in range(3)]
    R = math.sqrt(sum(x * x for x in c))
    up = [x / R for x in c]
    e = [-up[1], up[0], 0.0]
    el = math.sqrt(e[0] ** 2 + e[1] ** 2)
    e = [x / el for x in e]
    nn = [up[1] * e[2] - up[2] * e[1], up[2] * e[0] - up[0] * e[2], up[0] * e[1] - up[1] * e[0]]
    nr, na = 16, 64
    pts = []
    for i in range(nr + 1):
        r = 1500.0 + (8000.0 - 1500.0) * (i / nr) ** 1.5
        h = 20.0 + 450.0 * (r - 1500.0) / 6500.0 - r * r / (2.0 * R)
        for j in range(na):
            a = 2.0 * math.pi * (j + 0.5 * (i % 2)) / na
            x, y = r * math.cos(a), r * math.sin(a)
            pts.append([c[k] + e[k] * x + nn[k] * y + up[k] * h for k in range(3)])
    faces = []
    for i in range(nr):
        for j in range(na):
            a, b = i * na + j, i * na + (j + 1) % na
            cc, d = (i + 1) * na + j, (i + 1) * na + (j + 1) % na
            faces += [(a, b, d), (a, d, cc)]
    with open(path, "w") as f:
        f.write("ply\nformat ascii 1.0\ncomment synthetic occluder ring for test_gpu_host_parity\n"
                f"element vertex {len(pts)}\nproperty double x\nproperty double y\nproperty double z\n"
                f"element face {len(faces)}\nproperty list uchar int vertex_indices\nend_header\n")
        for p in pts:
            f.write(f"{p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
        for t in faces:
            f.write(f"3 {t[0]} {t[1]} {t[2]}\n")


def config(out, ring, use_gpu):
    occ = f"  occluder_file: {ring}\n" if ring else ""
    return ("simulation:\n  start_time: '2026-04-05T00:00:00'\n  duration_days: 29.53\n"
            f"  timesteps: {FRAMES}\n  body: Moon\n"
            f"mesh:\n  import_file: {CORE}\n{occ}  scale_factor: 1.0\n  use_topographic: false\n"
            "thermal:\n  run_thermal_mesh: false\n  use_gpu: false\n"
            "radiation:\n  enabled: false\n"
            "illumination:\n  include_primary_reflection: false\n  use_batch_processing: true\n"
            f"  use_gpu: {'true' if use_gpu else 'false'}\n  gpu_devices: 1\n  bvh_cache_file: 'off'\n"
            f"output:\n  directory: {out}\n  format: binary\n")


def frames(run):
    out = {}
    for p in sorted((run / "direct_solar").glob("*.bin")):
        b = p.read_bytes()
        n = struct.unpack("<I", b[:4])[0]
        out[p.name] = struct.unpack(f"<{n}f", b[4:4 + 4 * n])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--binary", default=str(REPO / "build/bin/ptm_thermal_sim"))
    ap.add_argument("--keep", default=None, help="keep the runs in this directory")
    a = ap.parse_args()
    binary = Path(a.binary).resolve()

    probe = subprocess.run([str(binary), "--list-gpu-devices"], cwd=REPO, capture_output=True,
                           text=True, timeout=120, env={**os.environ, "PTM_OPENCL_DEVICE_TYPE": "cpu"})
    if "without OpenCL support" in (probe.stdout + probe.stderr) or probe.returncode != 0:
        print("[test_gpu_host_parity] SKIP: this binary has no OpenCL GPU support; nothing was "
              "checked. Build with -DPTM_ENABLE_GPU=ON.")
        sys.exit(2)

    fails = 0

    def check(name, ok, detail=""):
        nonlocal fails
        print(f"  [{name}] {'PASS' if ok else 'FAIL'}" + (f"  {detail}" if detail and not ok else ""))
        fails += not ok

    work = Path(a.keep) if a.keep else Path(tempfile.mkdtemp(prefix="ptm_gpu_host_parity_"))
    work.mkdir(parents=True, exist_ok=True)
    try:
        ring = work / "ring.ply"
        write_ring(ring, read_core_vertices())

        def run(name, env_extra, use_gpu, with_ring=True):
            out = work / name
            cfg = work / f"{name}.yaml"
            cfg.write_text(config(out, ring if with_ring else None, use_gpu))
            env = {k: v for k, v in os.environ.items() if not k.startswith("PTM_ILLUM_")
                   and k not in ("PTM_BVH_HOST_KERNEL", "PTM_OPENCL_DEVICE_TYPE")}
            env.update(env_extra)
            p = subprocess.run([str(binary), "--config", str(cfg)], cwd=REPO, env=env,
                               capture_output=True, text=True, timeout=1800)
            (work / f"{name}.log").write_text(p.stdout + p.stderr)
            return p.returncode, p.stdout + p.stderr, out

        host_env = {"PTM_BVH_HOST_KERNEL": "1", "PTM_OPENCL_DEVICE_TYPE": "cpu"}
        for label, recipe in (("recipe", RECIPE), ("default", {})):
            print(f"[{label} environment]")
            rc_c, log_c, cpu = run(f"cpu_{label}", recipe, False)
            rc_h, log_h, host = run(f"host_{label}", {**recipe, **host_env}, True)
            check(f"P3 {label}: both runs exit 0", rc_c == 0 and rc_h == 0, f"rc {rc_c}/{rc_h}")
            check(f"P3 {label}: HOST ran the GPU path on the host kernel",
                  "GPU ray tracer enabled" in log_h and "PTM_BVH_HOST_KERNEL=1" in log_h
                  and "falling back" not in log_h and "Disabling GPU ray tracing" not in log_h)
            prov = (host / "provenance.txt").read_text() if (host / "provenance.txt").exists() else ""
            check(f"P3 {label}: provenance records the GPU path with no fallback",
                  "illumination_gpu_fallback_pass=none" in prov, prov[-300:])
            fc, fh = frames(cpu), frames(host)
            check(f"P3 {label}: all {FRAMES} frames written", len(fc) == FRAMES and len(fh) == FRAMES,
                  f"{len(fc)}/{len(fh)}")
            lit = flips = diff = 0
            dmax = 0.0
            for k in fc:
                for x, y in zip(fc[k], fh.get(k, ())):
                    lit += x > 0
                    flips += (x > 0) != (y > 0)
                    if x != y:
                        diff += 1
                        dmax = max(dmax, abs(x - y))
            print(f"    {lit:,} CPU-lit face-frames; lit/dark disagreements {flips}; values differing "
                  f"{diff} (max |dF| {dmax:.3g} W/m2)")
            check(f"P1 {label}: 0 lit/dark disagreements CPU vs HOST", flips == 0 and lit > 0,
                  f"{flips} disagreements")
            if label == "recipe":
                rc_n, _, noring = run("noring_recipe", recipe, False, with_ring=False)
                fn = frames(noring)
                shadowed = sum((n_ > 0) and (x <= 0) for k in fc for x, n_ in zip(fc[k], fn.get(k, ())))
                print(f"    ring-shadowed core face-frames (lit without the ring): {shadowed:,}")
                check("P2 the ring shadows the core (non-vacuous)", rc_n == 0 and shadowed >= 100,
                      f"{shadowed}")
    finally:
        if not a.keep:
            shutil.rmtree(work, ignore_errors=True)

    print(f"[test_gpu_host_parity] {'PASS' if not fails else 'FAIL'}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
