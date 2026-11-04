#!/usr/bin/env python3
"""GPU illumination is opt-in experimental (PTM-794): refused without the opt-in, loud with it.

Runs `ptm_thermal_sim --validate-config` (parse only: no mesh, no SPICE, no
simulation) on small generated configs and checks:

  G1  illumination.use_gpu: true WITHOUT illumination.gpu_experimental is REFUSED
      (nonzero exit) and the refusal names PTM-794 and the opt-in key;
  G2  the CLI flag --use-gpu-illumination alone is refused the same way;
  G3  use_gpu: true WITH gpu_experimental: true is accepted and prints the
      EXPERIMENTAL warning;
  G4  the CLI pair --use-gpu-illumination --gpu-illumination-experimental is
      accepted with the warning, and the opt-in flag is not an unknown argument;
  G5  the default (CPU) config prints neither the refusal nor the warning, and
      gpu_experimental is not reported as an unknown key.

The guard is checked in every build (CPU-only as well as OpenCL), so this test
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

    refused = "FATAL [illumination]: illumination.use_gpu is true"
    warned = "WARNING: EXPERIMENTAL GPU illumination enabled"

    with tempfile.TemporaryDirectory(prefix="ptm_gpu_illum_guard_") as d:
        d = Path(d)
        no_opt = d / "no_opt.yaml"
        no_opt.write_text(BASE + "illumination:\n  use_gpu: true\n")
        rc, out = validate(a.binary, no_opt)
        check("G1 config without opt-in refused", rc != 0 and refused in out, out[-600:])
        check("G1 refusal names PTM-794 and the key",
              "PTM-794" in out and "gpu_experimental: true" in out)
        check("G1 refused before validate-config accepts", "ACCEPTED" not in out)

        plain = d / "plain.yaml"
        plain.write_text(BASE)
        rc, out = validate(a.binary, plain, ["--use-gpu-illumination"])
        check("G2 CLI flag without opt-in refused", rc != 0 and refused in out, out[-600:])

        opt = d / "opt.yaml"
        opt.write_text(BASE + "illumination:\n  use_gpu: true\n  gpu_experimental: true\n")
        rc, out = validate(a.binary, opt)
        check("G3 opt-in accepted", rc == 0 and "ACCEPTED" in out and refused not in out,
              out[-600:])
        check("G3 opt-in warns loudly", warned in out and "PTM-794" in out)
        check("G3 opt-in key is a known key", "unknown key: illumination.gpu_experimental" not in out)

        rc, out = validate(a.binary, plain,
                           ["--use-gpu-illumination", "--gpu-illumination-experimental"])
        check("G4 CLI opt-in accepted with warning",
              rc == 0 and warned in out and refused not in out, out[-600:])
        check("G4 CLI opt-in flag is known",
              "unrecognized command-line argument '--gpu-illumination-experimental'" not in out)

        rc, out = validate(a.binary, plain)
        check("G5 CPU default: no refusal, no warning",
              rc == 0 and refused not in out and warned not in out, out[-600:])

    print(f"[test_gpu_illumination_guard] {'PASS' if not fails else 'FAIL'}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
