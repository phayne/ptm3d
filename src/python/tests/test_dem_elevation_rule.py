#!/usr/bin/env python3
"""DEM tools agree on how a stored raster value becomes a height in meters.

The rule lives in scripts/dem_elevation.py and, line for line, in
src/cpp/tools/dem_elevation.h. Three groups of checks:

* **The rule** (always runs, no GDAL). Every case in scripts/dem_elevation_cases.tsv
  through the Python rule, and through build/bin/test_dem_elevation (the C++ rule)
  when that is built.

The binaries are taken from build/bin, or from $PTM_DEM_TOOLS_BIN_DIR (to run
these checks against another build). Missing prerequisites SKIP, naming what is
missing.

    python3 src/python/tests/test_dem_elevation_rule.py
"""
from __future__ import annotations

import atexit
import os
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import dem_elevation  # noqa: E402
import test_dem_tools as tdt  # noqa: E402  (synthetic-DTM WKT and PLY reader)

_TMP = []
atexit.register(lambda: [shutil.rmtree(d, ignore_errors=True) for d in _TMP])
BIN_DIR = Path(os.environ.get("PTM_DEM_TOOLS_BIN_DIR", REPO_ROOT / "build" / "bin"))
CASES = REPO_ROOT / "scripts" / "dem_elevation_cases.tsv"
R = 1737400.0


class Skip(Exception):
    pass


def _bin(name):
    p = BIN_DIR / name
    if not p.is_file() or not os.access(p, os.X_OK):
        raise Skip(f"missing binary: {p}. Build it with 'cmake --build build --target {name}'.")
    return p


def _gdal():
    try:
        from osgeo import gdal, osr
        import numpy as np
    except ImportError as exc:
        raise Skip(f"GDAL Python bindings or numpy not importable: {exc}")
    gdal.UseExceptions()
    return gdal, osr, np


def _cases():
    out = []
    for line in CASES.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        c = line.split("\t")
        assert len(c) == 11, f"malformed case line: {line!r}"
        out.append(c)
    return out


# ---------------------------------------------------------------- the rule
def test_python_rule_matches_the_case_table():
    o = (lambda s: None if s == "-" else float(s))
    bad = []
    for c in _cases():
        try:
            t = dem_elevation.resolve(o(c[1]), o(c[2]), "" if c[3] == "-" else c[3], c[4] == "1",
                                      o(c[5]), o(c[6]), float(c[7]))
            ok = (c[8] == "OK" and abs(t.scale - float(c[9])) <= 1e-9 * max(1, abs(float(c[9])))
                  and abs(t.offset - float(c[10])) <= 1e-9 * max(1, abs(float(c[10]))))
            got = f"OK {t.scale:g} {t.offset:g}"
        except dem_elevation.DemElevationError:
            ok, got = c[8] == "ERROR", "ERROR"
        if not ok:
            bad.append(f"{c[0]}: got {got}, expected {c[8]} {c[9]} {c[10]}")
    assert not bad, "Python rule disagrees with the case table:\n  " + "\n  ".join(bad)


def test_cpp_rule_matches_the_case_table():
    exe = _bin("test_dem_elevation")
    proc = subprocess.run([str(exe), str(CASES)], capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0, f"C++ rule disagrees with the case table:\n{proc.stdout}{proc.stderr}"
    assert f"PASS ({len(_cases())} cases)" in proc.stdout, proc.stdout


TESTS = [
    test_python_rule_matches_the_case_table,
    test_cpp_rule_matches_the_case_table,
]


def _pytest_wrap():
    try:
        import pytest
    except ImportError:
        return

    def make(fn):
        def wrapper():
            try:
                fn()
            except Skip as exc:
                pytest.skip(str(exc))
        wrapper.__name__ = fn.__name__
        wrapper.__doc__ = fn.__doc__
        return wrapper
    for fn in TESTS:
        globals()[fn.__name__] = make(fn)


_pytest_wrap()


def main():
    fails = skips = 0
    for t in TESTS:
        try:
            t()
            print(f"ok:   {t.__name__}")
        except Skip as e:
            skips += 1
            print(f"SKIP: {t.__name__}: {e}")
        except AssertionError as e:
            fails += 1
            print(f"FAIL: {t.__name__}: {e}")
    print(f"{len(TESTS) - fails - skips} passed, {fails} failed, {skips} skipped")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
