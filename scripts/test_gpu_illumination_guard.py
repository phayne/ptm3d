#!/usr/bin/env python3
"""GPU illumination runs with a NOTICE, not a refusal (PTM-794 fixed).

The GPU shadow-ray tracer was validated against CPU Embree on an AMD Radeon Pro
W6800X for lunar meshes. A run that asks for it is accepted and told what was
validated and how to check other hardware (`ptm gpu-selfcheck`). The opt-in
that the unreleased guard required (`illumination.gpu_experimental`,
`--gpu-illumination-experimental`) is retired; a config or command line that
still carries it is told so by the ordinary unknown-key/argument warning.

Runs `ptm_thermal_sim --validate-config` (parse only: no mesh, no SPICE, no
simulation) on small generated configs and checks:

  G1  illumination.use_gpu: true is ACCEPTED (rc 0, "ACCEPTED"), with no
      refusal, and prints the NOTICE naming the validated GPU, PTM-794 and
      `ptm gpu-selfcheck`;
  G2  the CLI flag --use-gpu-illumination alone is accepted the same way;
  G3  the default (CPU) config prints no NOTICE and no refusal;
  G4  the retired key illumination.gpu_experimental is reported as an unknown
      key (not silently dropped), and the run is still accepted;
  G5  the retired flag --gpu-illumination-experimental is reported as an
      unrecognized argument, not silently accepted.

The NOTICE is printed in every build (CPU-only as well as OpenCL), so this test
needs no GPU.

Usage: python3 scripts/test_gpu_illumination_guard.py [--binary build/bin/ptm_thermal_sim]
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
BASE = "simulation:\n  body: Moon\nthermal:\n  run_thermal_mesh: true\n"


def validate(binary, cfg, extra=()):
    p = subprocess.run([str(binary), "--config", str(cfg), "--validate-config", *extra],
                       cwd=REPO, capture_output=True, text=True, timeout=120)
    return p.returncode, p.stdout + p.stderr


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--binary", default=str(REPO / "build/bin/ptm_thermal_sim"))
    a = ap.parse_args()
    fails = 0

    def check(name, ok, detail=""):
        nonlocal fails
        print(f"  [{name}] {'PASS' if ok else 'FAIL'}" + (f"  {detail}" if detail and not ok else ""))
        fails += not ok

    notice = "[illumination] NOTICE: GPU illumination requested"
    refused = "FATAL [illumination]"

    def notice_complete(out):
        return (notice in out and "W6800X" in out and "PTM-794" in out
                and "ptm gpu-selfcheck" in out)

    with tempfile.TemporaryDirectory(prefix="ptm_gpu_illum_guard_") as d:
        d = Path(d)
        gpu = d / "gpu.yaml"
        gpu.write_text(BASE + "illumination:\n  use_gpu: true\n")
        rc, out = validate(a.binary, gpu)
        check("G1 use_gpu: true accepted", rc == 0 and "ACCEPTED" in out and refused not in out,
              out[-600:])
        check("G1 NOTICE names the validated GPU, PTM-794 and the self-check",
              notice_complete(out), out[-600:])

        plain = d / "plain.yaml"
        plain.write_text(BASE)
        rc, out = validate(a.binary, plain, ["--use-gpu-illumination"])
        check("G2 CLI flag accepted with the NOTICE",
              rc == 0 and refused not in out and notice_complete(out), out[-600:])

        rc, out = validate(a.binary, plain)
        check("G3 CPU default: no NOTICE, no refusal",
              rc == 0 and notice not in out and refused not in out, out[-600:])

        retired = d / "retired.yaml"
        retired.write_text(BASE + "illumination:\n  use_gpu: true\n  gpu_experimental: true\n")
        rc, out = validate(a.binary, retired)
        check("G4 retired key reported, run still accepted",
              rc == 0 and "unknown key: illumination.gpu_experimental" in out, out[-600:])

        rc, out = validate(a.binary, plain,
                           ["--use-gpu-illumination", "--gpu-illumination-experimental"])
        check("G5 retired CLI flag reported",
              "unrecognized command-line argument '--gpu-illumination-experimental'" in out,
              out[-600:])

    print(f"[test_gpu_illumination_guard] {'PASS' if not fails else 'FAIL'}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
