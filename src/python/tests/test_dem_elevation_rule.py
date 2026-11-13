#!/usr/bin/env python3
"""DEM tools agree on how a stored raster value becomes a height in meters.

The rule lives in scripts/dem_elevation.py and, line for line, in
src/cpp/tools/dem_elevation.h. Three groups of checks:

* **The rule** (always runs, no GDAL). Every case in scripts/dem_elevation_cases.tsv
  through the Python rule, and through build/bin/test_dem_elevation (the C++ rule)
  when that is built.

* **dtm_to_ply_delaunay end to end** (needs the binary and GDAL's Python bindings).
  Synthetic south-polar GeoTIFFs: an Int16 LOLA-style raster with scale 0.5 and a
  reference radius 600 m above the body radius must give |XYZ| = R + DN*0.5 + 600
  at every vertex; the same with the radius equal to R gives R + DN*0.5; a Float32
  raster in meters is unchanged; a raster in kilometers that does not say so stops
  with an error instead of writing a mesh.

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


# ------------------------------------------------------ dtm_to_ply_delaunay
N, PX = 120, 20.0      # 2.4 km square at the south pole


def _dn_grid(np):
    xx, yy = np.meshgrid(np.arange(N), np.arange(N))
    return (400 + 3 * xx - 2 * yy + 50 * ((xx * yy) % 7)).astype(np.int16)


def _write_polar(path, array, gdal_type, scale=None, offset=None, nodata=None):
    gdal, osr, np = _gdal()
    ds = gdal.GetDriverByName("GTiff").Create(str(path), N, N, 1, gdal_type)
    ds.SetGeoTransform((-N / 2 * PX, PX, 0.0, N / 2 * PX, 0.0, -PX))
    srs = osr.SpatialReference()
    srs.ImportFromWkt(tdt._POLAR_STEREO_WKT)
    ds.SetProjection(srs.ExportToWkt())
    b = ds.GetRasterBand(1)
    if scale is not None:
        b.SetScale(scale)
    if offset is not None:
        b.SetOffset(offset)
    if nodata is not None:
        b.SetNoDataValue(nodata)
    b.WriteArray(array)
    ds.FlushCache()
    del ds


def _convert(src, out, extra=()):
    exe = _bin("dtm_to_ply_delaunay")
    return subprocess.run([str(exe), str(src), str(out), "--body-radius", str(R), "--no-qc"]
                          + list(extra), capture_output=True, text=True, timeout=600)


def _check_radii(mesh, label):
    bad = [i for i in range(mesh.n_vertices)
           if abs(mesh.radius[i] - (R + mesh.elev[i])) > 1e-3]
    assert not bad, f"{label}: |XYZ| != R + elevation at {len(bad)} vertices"


def _dtm_case(tag, array, gdal_type, scale, offset, expected_shift, extra=()):
    gdal, _, np = _gdal()
    with tempfile.TemporaryDirectory() as td:
        src, out = Path(td) / f"{tag}.tif", Path(td) / f"{tag}.ply"
        _write_polar(src, array, gdal_type, scale, offset)
        proc = _convert(src, out, extra)
        assert proc.returncode == 0, f"{tag}: exit {proc.returncode}\n{proc.stdout[-1500:]}{proc.stderr[-1500:]}"
        mesh = tdt._read_ply(out)
    _check_radii(mesh, tag)
    stored = np.asarray(array, dtype=np.float64)
    expect = (stored * (scale if scale is not None else 1.0)) + expected_shift
    # The mesh keeps every grid point at --downsample 1 (no ROI): compare the sets.
    got = np.sort(np.asarray(mesh.elev))
    want = np.sort(expect.ravel())
    assert got.size == want.size, f"{tag}: {got.size} vertices, expected {want.size}"
    err = np.abs(got - want).max()
    assert err < 1e-2, (f"{tag}: vertex heights differ from stored*scale + {expected_shift:g} "
                        f"by up to {err:.3f} m (median shift {np.median(got - want):.3f} m)")
    return proc


def test_dtm_reference_radius_above_the_body_radius_is_kept():
    """PTM-912: an offset 600 m above R gives heights 600 m higher, not dropped."""
    gdal, _, np = _gdal()
    proc = _dtm_case("lola_r_plus_600", _dn_grid(np), gdal.GDT_Int16, 0.5, R + 600.0, 600.0)
    assert "Elevation: height (m) = stored value * 0.5 [file metadata]" in proc.stdout, \
        "the log must say which scale and offset were used and where they came from"


def test_dtm_lola_style_label_is_unchanged():
    gdal, _, np = _gdal()
    _dtm_case("lola_r", _dn_grid(np), gdal.GDT_Int16, 0.5, R, 0.0)


def test_dtm_float_meters_without_metadata_is_unchanged():
    gdal, _, np = _gdal()
    _dtm_case("float_m", _dn_grid(np).astype(np.float32) * 0.5, gdal.GDT_Float32, None, None, 0.0)


def test_dtm_kilometers_without_a_unit_stop_with_an_error():
    gdal, _, np = _gdal()
    with tempfile.TemporaryDirectory() as td:
        src, out = Path(td) / "km.tif", Path(td) / "km.ply"
        _write_polar(src, _dn_grid(np).astype(np.float32) / 1000.0, gdal.GDT_Float32, 1.0, 1737.4)
        proc = _convert(src, out)
        assert proc.returncode != 0 and not out.exists(), (
            f"a kilometer raster without a unit was converted (exit {proc.returncode}); "
            f"its heights would be 1000x too small\n{proc.stdout[-800:]}")
        assert "kilometers" in proc.stderr, proc.stderr
        # ... and converts correctly with an explicit scale and offset
        proc = _convert(src, out, ["--elev-scale", "1000", "--elev-offset", "0"])
        assert proc.returncode == 0, proc.stderr[-800:]
        mesh = tdt._read_ply(out)
        want = np.sort((_dn_grid(np).astype(np.float32) / 1000.0).astype(np.float64).ravel() * 1000)
        assert np.abs(np.sort(mesh.elev) - want).max() < 1e-2


TESTS = [
    test_python_rule_matches_the_case_table,
    test_cpp_rule_matches_the_case_table,
    test_dtm_reference_radius_above_the_body_radius_is_kept,
    test_dtm_lola_style_label_is_unchanged,
    test_dtm_float_meters_without_metadata_is_unchanged,
    test_dtm_kilometers_without_a_unit_stop_with_an_error,
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
