#!/usr/bin/env python3
"""Test: a chunked run names its unmerged flux-cube pieces (PTM-598).

Under scripts/run_chunked_simulation.py each chunk runs its own flux-cube pass
and writes its own .fcb, and ptm_merge_chunks merges statistics only. Before
this check existed a 2-chunk window ended with "[chunked] DONE." and the
per-chunk products sat in chunk_NNN/flux_cubes/ with no mention anywhere.

The arms exercise the two messages the runner prints:
  * the launch-time warning, when flux_cubes.enabled and more than one chunk;
  * the end-of-run listing of every .fcb and per-cube CSV directory, with the
    output-step range each covers, as a WARNING for more than one chunk and as
    a plain locator for one chunk; nothing at all when no chunk wrote cubes.

No engine is run. Needs PyYAML (the runner imports it).
  python3 scripts/test_chunked_flux_cube_notice.py
"""
import importlib.util
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load_runner():
    spec = importlib.util.spec_from_file_location(
        "run_chunked_simulation", HERE / "run_chunked_simulation.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def make_plan(root: Path, steps):
    plan = []
    for i, ts in enumerate(steps):
        plan.append({"idx": i, "ts": ts, "chunk_dir": root / f"chunk_{i:03d}"})
    return plan


def touch(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"FCB1")


def main():
    rc = load_runner()
    fails = 0

    def check(ok, what):
        nonlocal fails
        print(("PASS " if ok else "FAIL ") + what)
        if not ok:
            fails += 1

    check(rc.flux_cubes_enabled({"flux_cubes": {"enabled": True}}),
          "flux_cubes.enabled: true is detected")
    check(not rc.flux_cubes_enabled({"flux_cubes": {"enabled": False}})
          and not rc.flux_cubes_enabled({}) and not rc.flux_cubes_enabled({"flux_cubes": None}),
          "disabled, absent and null flux_cubes are not")

    w = rc.flux_cube_split_warning(3)
    check("WARNING" in w and "3 chunks" in w and "no tool merges them" in w
          and "--chunks 1" in w,
          "launch warning names the chunk count, the missing merge and the remedy")

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)

        plan = make_plan(root / "two", [59, 59])
        check(rc.flux_cube_pieces_report(plan) is None,
              "no flux_cubes directory in any chunk -> no report")

        for p in plan:
            touch(p["chunk_dir"] / "flux_cubes" / "flux_cubes_site_h0p35.fcb")
        (plan[1]["chunk_dir"] / "flux_cubes" / "csv_flux_cubes_site_h0p35").mkdir()
        r = rc.flux_cube_pieces_report(plan) or ""
        print("    " + r.replace("\n", "\n    "))
        check(r.startswith("[chunked] WARNING") and "NOT merged" in r,
              "two chunks -> a WARNING that the product was not merged")
        check("output steps [0, 59)" in r and "output steps [59, 118)" in r,
              "each piece carries the output-step range it covers")
        check(r.count("flux_cubes_site_h0p35.fcb") == 2
              and "csv_flux_cubes_site_h0p35" in r,
              "every .fcb and the per-cube CSV directory are named")

        one = make_plan(root / "one", [118])
        touch(one[0]["chunk_dir"] / "flux_cubes" / "flux_cubes.fcb")
        r1 = rc.flux_cube_pieces_report(one) or ""
        print("    " + r1.replace("\n", "\n    "))
        check("WARNING" not in r1 and "whole window" in r1
              and "output steps [0, 118)" in r1 and "flux_cubes.fcb" in r1,
              "one chunk -> no warning, and the location of the product is given")

    print("\n[test_chunked_flux_cube_notice] %s (%d failure(s))"
          % ("PASS" if fails == 0 else "FAIL", fails))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
