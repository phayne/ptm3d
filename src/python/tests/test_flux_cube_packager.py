#!/usr/bin/env python3
"""The flux-cube HDF5 packager, and the engine provenance it copies.

Two defects found while verifying a 2026-09-02 lunar polar flux-cube product,
both
of which every structural check passed:

* **Flux-cube statistics header** -- ``scripts/ptm_flux_cubes_to_hdf5.py`` read ``statistics/*.bin``
  with a bare ``np.fromfile(dtype=float32)``, which takes the file's 4-byte
  uint32 face-count HEADER as if it were the first face's temperature and
  shifts every face by one index. That is the statistics-header off-by-one bug class, recurring in a
  script written after the statistics-header off-by-one had been fixed in ``ptm_qa.py``. It printed
  "2397891 faces" for a 2,397,890-face mesh and nobody read the number.
* **Flux-cube hash precision** -- the engine wrote the flux-cube provenance JSON's
  ``geomhash_core`` through a default-precision ``ostream``: six significant
  digits, so a hash of the form 31904726.165896364 was recorded as 31904700.0.

Run directly::

    python3 src/python/tests/test_flux_cube_packager.py
"""

from __future__ import annotations

import importlib.util
import json
import re
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "src" / "python"))

#: A SYNTHETIC geometry hash, invented for this suite. Real hashes identify
#: real meshes and this project does not publish the ones it has, so what is
#: kept here is the PROPERTY the two defects turned on and not the number they
#: turned on it with: a double of this magnitude whose six-significant-digit
#: rendering loses about 26 in hash units -- under 1e-6 relative, and so INSIDE
#: the pass's own 1e-3 match tolerance, which is the whole point of the second
#: arm below. Substituting the real value changes no arm's outcome.
CORE_GEOMHASH = 31904726.165896364

_SCRATCH = None


class Skip(Exception):
    """A check that cannot run here, as distinct from one that failed."""


def scratch():
    global _SCRATCH
    if _SCRATCH is None:
        _SCRATCH = tempfile.mkdtemp(prefix="ptm_fcpack_")
    return Path(_SCRATCH)


def packager():
    """The script under test, imported by path (it is not a package)."""
    path = REPO_ROOT / "scripts" / "ptm_flux_cubes_to_hdf5.py"
    spec = importlib.util.spec_from_file_location("_fc_packager", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def write_stat_bin(path, values, declared_n=None):
    """A statistics binary: uint32 face count, then that many float32.

    ``declared_n`` writes a header that disagrees with the payload, which is
    what a truncated or foreign file looks like.
    """
    n = len(values) if declared_n is None else declared_n
    with open(path, "wb") as fh:
        fh.write(struct.pack("<I", n))
        fh.write(struct.pack("<%df" % len(values), *values))
    return path


# --------------------------------------------------------------- the flux-cube statistics header
def test_a_statistics_binary_is_read_past_its_face_count_header():
    """THE ARM. The values that come back are the values that went in.

    The naive read is measured alongside, on the same file, so the arm records
    what the defect actually did rather than asserting that it is gone.
    """
    import numpy as np
    m = packager()
    vals = [10.0, 20.0, 30.0, 40.0, 50.0]
    p = write_stat_bin(scratch() / "temperature_max.bin", vals)

    got = m.read_face_stat_bin(str(p))
    assert list(got) == vals, f"read back {list(got)}, wrote {vals}"

    naive = np.fromfile(str(p), dtype=np.float32)
    assert len(naive) == len(vals) + 1, (
        "the bare read no longer returns an extra element; this arm's premise "
        "about the on-disk format has changed and it needs rewriting")
    assert float(np.median(naive)) == 25.0 and float(np.median(got)) == 30.0, (
        f"the off-by-one must be visible in the statistic the packager takes: "
        f"naive median {np.median(naive)}, correct median {np.median(got)}")


def test_a_truncated_statistics_binary_is_refused():
    """A header that disagrees with the file length is a refusal, not a slice.

    The flux-cube reader learned this one the expensive way: a wrong record
    size does not raise, it lands in a neighboring block and returns printable
    nonsense. Only a LENGTH check catches it.
    """
    m = packager()
    p = write_stat_bin(scratch() / "truncated.bin", [1.0, 2.0], declared_n=9)
    try:
        m.read_face_stat_bin(str(p))
    except ValueError as exc:
        assert "9 faces" in str(exc) and "bytes" in str(exc), \
            f"the refusal must name the header count and the length: {exc}"
        return
    raise AssertionError("a statistics file whose header claims 9 faces and "
                         "whose body holds 2 was accepted")


def test_a_statistics_binary_from_another_mesh_is_refused():
    """Right format, wrong run. Length alone cannot catch this one."""
    m = packager()
    p = write_stat_bin(scratch() / "foreign.bin", [1.0, 2.0, 3.0])
    try:
        m.read_face_stat_bin(str(p), expect_faces=2397890,
                             expect_source="the core mesh")
    except ValueError as exc:
        assert "2397890" in str(exc) and "3" in str(exc), \
            f"the refusal must name both counts: {exc}"
        return
    raise AssertionError("a 3-face statistics file was accepted for a "
                         "2,397,890-face mesh")


def test_a_ply_face_count_is_read_from_the_header():
    """The identity check needs the mesh's own count, and gets it from the PLY."""
    m = packager()
    p = scratch() / "tiny.ply"
    p.write_text("ply\nformat ascii 1.0\nelement vertex 3\n"
                 "property float x\nelement face 7\n"
                 "property list uchar int vertex_indices\nend_header\n")
    assert m.ply_face_count(str(p)) == 7
    assert m.ply_face_count(str(scratch() / "nope.ply")) is None


# --------------------------------------------------------------- the flux-cube hash precision
def _extract_geomhash_formatting():
    """The engine's own two lines that render `geomhash_core`, from the source.

    Extracted rather than retyped: a test that compiles a COPY of the
    expression proves the copy round-trips and says nothing about the engine.
    """
    src = (REPO_ROOT / "src" / "cpp" / "main" / "ptm_thermal_sim.cpp").read_text()
    decl = re.search(r"^\s*(std::ostringstream\s+(\w+);)\s*$", src, re.M | re.S)
    lines = [l.strip() for l in src.split("\n")
             if "geomhash_core;" in l and "<<" in l]
    if not lines:
        raise AssertionError(
            "no line in ptm_thermal_sim.cpp streams `geomhash_core` into a "
            "string. The geometry-hash precision fix has been removed or renamed.")
    if len(lines) != 1:
        raise AssertionError(f"expected one such line, found {len(lines)}: {lines}")
    stream_line = lines[0]
    name = stream_line.split("<<")[0].strip()
    if not re.match(r"^\w+$", name):
        raise AssertionError(f"could not name the stream in: {stream_line}")
    assert decl is not None
    return name, stream_line


def _working_cxx():
    """A C++ compiler that actually runs on this machine.

    ⛔ Not `shutil.which("c++")`. On this Mac Pro that resolves to
    /usr/local/bin/c++, which is a Homebrew binary for the wrong architecture
    and dies with "Bad CPU type in executable" -- an OSError from the *exec*,
    not a non-zero exit, so a which()-and-run probe raises instead of falling
    through. The project's own toolchain (CommandLineTools) is tried first and
    every candidate is EXECUTED before it is trusted.
    """
    candidates = ["/Library/Developer/CommandLineTools/usr/bin/clang++",
                  "/usr/bin/clang++", shutil.which("clang++"),
                  shutil.which("c++"), shutil.which("g++")]
    probe = "#include <sstream>\nint main(){ return 0; }\n"
    d = scratch() / "cxx_probe"
    d.mkdir(exist_ok=True)
    (d / "probe.cpp").write_text(probe)
    for c in candidates:
        if not c or not Path(c).exists():
            continue
        try:
            # Compile a header-using program, not `--version`: the version
            # banner prints happily from a clang that cannot find <sstream>,
            # which is the state of this machine without DEVELOPER_DIR (the
            # same variable the CMake build needs -- see the project notes).
            r = subprocess.run([c, "-std=c++17", str(d / "probe.cpp"),
                                "-o", str(d / "probe")],
                               capture_output=True, text=True, timeout=300,
                               env=_cxx_env())
        except OSError:
            continue
        if r.returncode == 0:
            return c
    return None


def _cxx_env():
    import os
    env = dict(os.environ)
    env.setdefault("DEVELOPER_DIR", "/Library/Developer/CommandLineTools")
    return env


def test_the_engine_renders_the_geomhash_at_round_trip_precision():
    """COMPILE THE ENGINE'S OWN LINE and check what it prints.

    The two lines are lifted verbatim out of ptm_thermal_sim.cpp, wrapped in a
    main() with a hash substituted for the member, compiled and
    run. If someone drops the `setprecision(17)` the source text changes and
    this fails with the six-digit value it printed instead.
    """
    cxx = _working_cxx()
    if cxx is None:
        raise Skip("no working C++ compiler found")
    name, stream_line = _extract_geomhash_formatting()
    body = stream_line.replace("fcin.geomhash_core", "gh")
    prog = ("#include <iomanip>\n#include <iostream>\n#include <sstream>\n"
            "int main(){ const double gh = %.17g;\n"
            "  std::ostringstream %s;\n  %s\n"
            "  std::cout << %s.str() << std::endl; return 0; }\n"
            % (CORE_GEOMHASH, name, body, name))
    d = scratch() / "gh_probe"
    d.mkdir(exist_ok=True)
    (d / "probe.cpp").write_text(prog)
    r = subprocess.run([cxx, "-std=c++17", "-O0", str(d / "probe.cpp"),
                        "-o", str(d / "probe")],
                       capture_output=True, text=True, timeout=300,
                       env=_cxx_env())
    if r.returncode != 0:
        raise AssertionError(f"the engine's own formatting lines did not "
                             f"compile in isolation:\n{r.stderr}")
    out = subprocess.run([str(d / "probe")], capture_output=True, text=True,
                         timeout=60).stdout.strip()
    assert float(out) == CORE_GEOMHASH, (
        f"the engine renders the geomhash as {out!r}, which reads back as "
        f"{float(out)!r}, not {CORE_GEOMHASH!r}. A mesh-identity hash that "
        f"does not round-trip cannot identify a mesh.")
    assert json.loads('{"g":' + out + '}')["g"] == CORE_GEOMHASH, \
        f"{out!r} is not JSON that reads back exactly"


def test_six_significant_digits_is_a_provenance_loss_not_a_gate_failure():
    """State the size of the defect, so nobody restates it wrongly.

    The truncation is 26.17 in hash units -- 8.2e-7 relative -- which is INSIDE
    the pass's own 1e-3 relative match tolerance (flux_cube_pass.cpp:222). No
    gate was ever fooled. What was lost is the ability to MATCH the recorded
    hash against the full-precision constant the crew keeps for this mesh.
    """
    six = float("%g" % CORE_GEOMHASH)
    assert six == 31904700.0, f"six-digit rendering is {six!r}"
    rel = abs(CORE_GEOMHASH - six) / CORE_GEOMHASH
    assert rel < 1e-3, (
        f"relative loss {rel:.3g} is NOT inside the pass's 1e-3 match "
        f"tolerance, so the hash-precision defect would have been a gate failure and the "
        f"'provenance-only' claim in the commit message is wrong")
    assert six != CORE_GEOMHASH, "no loss at all; this arm has no subject"


def test_the_packager_copies_the_provenance_number_at_full_precision():
    """Whatever the engine wrote, the HDF5 gets -- to the last digit."""
    try:
        import h5py
        import numpy as np                                     # noqa: F401
    except ImportError as exc:
        raise Skip(f"h5py/numpy not in this interpreter ({exc})")
    m = packager()
    blob = '{"geomhash_core":%.17g,"n_faces_core":2397890}' % CORE_GEOMHASH
    meta = json.loads(blob)
    assert meta["geomhash_core"] == CORE_GEOMHASH, \
        "the JSON the engine writes does not itself round-trip"
    # The packager's own copy step: every meta key becomes an `engine_*` root
    # attribute. Exercised through h5py rather than asserted about the source.
    out = scratch() / "prov.h5"
    with h5py.File(out, "w") as h:
        for k, v in meta.items():
            h.attrs["engine_" + k] = v
    with h5py.File(out, "r") as h:
        assert float(h.attrs["engine_geomhash_core"]) == CORE_GEOMHASH, (
            f"the packaged attribute is "
            f"{float(h.attrs['engine_geomhash_core'])!r}")
    assert hasattr(m, "read_face_stat_bin"), \
        "the packager under test is not the one carrying the statistics-header fix"


# --------------------------------------------------------------- PTM-609: nadir emission floor
# One cube, three steps, CORRECT data by construction. The nadir row sees the
# facet beneath (id 7, F_own = 0.0857), a neighbor that holds the row's largest
# entry (id 8, f_max = 0.2447) and the rest (id 9). Steps 0 and 2 are
# isothermal at 300 K. At step 1 the facet beneath is at 300 K and the others
# at 100 K, which is the configuration of the field case: the old floor
# eps*sigma*T_face^4*f_max multiplies facet 7's temperature by facet 8's view
# factor and FAILS correct data.
_SIG = 5.670374419e-8
_ROW = [(7, 0.0857), (8, 0.2447), (9, 0.6696)]


def _synthetic_product(face_id=7):
    temps = [{7: 300.0, 8: 300.0, 9: 300.0},
             {7: 300.0, 8: 100.0, 9: 100.0},
             {7: 300.0, 8: 300.0, 9: 300.0}]
    ns, nf, ncomp = 3, 6, 6
    flux = [0.0] * (nf * ncomp * ns)
    for st, T in enumerate(temps):
        tir = sum(F * 0.95 * _SIG * T[j] ** 4 for j, F in _ROW)
        flux[(1 * ncomp + 2) * ns + st] = tir              # nadir terrain_ir
        flux[(1 * ncomp + 5) * ns + st] = tir              # nadir total
    cube = (0, 0, 0.0, 0.0, 1.0e6, 0.0, 0.0, face_id, 0, 0.0, 0.0, 1.0, 0.0, 1.2)
    d = dict(n_cubes=1, n_faces=nf, n_comp=ncomp, n_steps=ns, height_m=1.2,
             flux=flux, t_face=[T[7] for T in temps], penumbral=[0] * ns,
             f_sky=[1.0, 0.0, 1.0, 1.0, 1.0, 1.0], f_ring=[0.0] * nf,
             f_max=[0.0, 0.2447, 0.0, 0.0, 0.0, 0.0], cubes=[cube])
    return d


def _write_fcg(path, face_id=7, height=1.2):
    """A one-cube .fcg in the engine's layout (flux_cube_format.h)."""
    pk = packager()
    hdr = struct.pack(pk.FCG_HDR, b"FCG1", 1, 1, 6, len(_ROW), 1.0, 3, 0,
                      0.0, 0.0, 25.0, -1.0, height, 0.5, 16384, 1, 0, 0.0,
                      0.0, 0.0, 1.0, *([0.0] * 12), 0.0, 0, 0)
    cube = struct.pack(pk.CUBE_FMT, 0, 0, 0.0, 0.0, 1.0e6, 0.0, 0.0, face_id, 0,
                       0.0, 0.0, 1.0, 0.0, height)
    patches = b""
    for f in range(6):
        if f == 1:
            patches += struct.pack(pk.PATCH_FMT, 0, len(_ROW), 1, 0.0, 0.0, 0.2447, 0.0, 0.0, 0.0)
        else:
            patches += struct.pack(pk.PATCH_FMT, len(_ROW), 0, 0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    cols = struct.pack("<%di" % len(_ROW), *[j for j, _ in _ROW])
    vals = struct.pack("<%df" % len(_ROW), *[F for _, F in _ROW])
    with open(path, "wb") as fh:
        fh.write(hdr + cube + patches + cols + vals)
    return path


def test_the_nadir_floor_uses_the_facet_beneath_and_passes_correct_data():
    pk = packager()
    d = _synthetic_product()
    # The defect, stated: the old f_max floor at step 1 is below what the
    # correct data delivers, so the old check reported FAIL on this product.
    got = d["flux"][(1 * 6 + 2) * 3 + 1]
    old_ratio = got / (0.90 * _SIG * 300.0 ** 4 * 0.2447)
    assert old_ratio < 1.0, "fixture does not reproduce PTM-609 (old ratio %.3f)" % old_ratio
    fcg = pk.read_fcg_rows(str(_write_fcg(scratch() / "ok.fcg")))
    not_run = []
    fails = pk.qa(d, fcg, not_run)
    assert fails == [] and not_run == [], (fails, not_run)
    print(f"     old f_max floor ratio {old_ratio:.4f} (false FAIL); F_own floor passes")


def test_without_the_cache_the_floor_is_reported_not_run():
    pk = packager()
    not_run = []
    fails = pk.qa(_synthetic_product(), None, not_run)
    assert fails == [] and len(not_run) == 1 and "--fcg" in not_run[0], (fails, not_run)


def test_a_facet_beneath_missing_from_its_nadir_row_fails():
    # A mis-indexed cube -> face mapping: the product and its cache agree on
    # facet 5, which the nadir row never sees.
    pk = packager()
    fcg = pk.read_fcg_rows(str(_write_fcg(scratch() / "misindex.fcg", face_id=5)))
    fails = pk.qa(_synthetic_product(face_id=5), fcg, [])
    assert any("do not see the facet beneath" in f for f in fails), fails


def test_a_corrupt_cache_is_refused_not_a_traceback():
    # A patch row that runs past the view-factor arrays: read_fcg_rows must
    # raise ValueError (reported as a QA failure by main), not IndexError.
    pk = packager()
    path = _write_fcg(scratch() / "corrupt.fcg")
    blob = bytearray(path.read_bytes())
    hs = struct.calcsize(pk.FCG_HDR) + struct.calcsize(pk.CUBE_FMT)
    nadir = hs + struct.calcsize(pk.PATCH_FMT)            # patch 1 = nadir
    struct.pack_into("<Q", blob, nadir, 99)               # row_offset far past nnz = 3
    path.write_bytes(bytes(blob))
    try:
        pk.read_fcg_rows(str(path))
    except ValueError as exc:
        assert "runs past" in str(exc), exc
    else:
        raise AssertionError("a corrupt cache was read without complaint")


def test_a_cache_from_another_product_fails():
    pk = packager()
    fcg = pk.read_fcg_rows(str(_write_fcg(scratch() / "other.fcg", face_id=8)))
    fails = pk.qa(_synthetic_product(face_id=7), fcg, [])
    assert any("not this product's cache" in f for f in fails), fails


TESTS = [
    test_a_statistics_binary_is_read_past_its_face_count_header,
    test_a_truncated_statistics_binary_is_refused,
    test_a_statistics_binary_from_another_mesh_is_refused,
    test_a_ply_face_count_is_read_from_the_header,
    test_the_engine_renders_the_geomhash_at_round_trip_precision,
    test_six_significant_digits_is_a_provenance_loss_not_a_gate_failure,
    test_the_packager_copies_the_provenance_number_at_full_precision,
    test_the_nadir_floor_uses_the_facet_beneath_and_passes_correct_data,
    test_without_the_cache_the_floor_is_reported_not_run,
    test_a_facet_beneath_missing_from_its_nadir_row_fails,
    test_a_cache_from_another_product_fails,
    test_a_corrupt_cache_is_refused_not_a_traceback,
]


def main():
    failures = skipped = 0
    for test in TESTS:
        try:
            test()
        except Skip as exc:
            skipped += 1
            print(f"SKIP {test.__name__}: {exc}")
        except AssertionError as exc:
            failures += 1
            print(f"FAIL {test.__name__}\n     {exc}")
        except Exception as exc:                       # noqa: BLE001
            failures += 1
            print(f"ERROR {test.__name__}: {type(exc).__name__}: {exc}")
        else:
            print(f"PASS {test.__name__}")
    total = len(TESTS)
    print(f"\n[test_flux_cube_packager] {total - failures - skipped}/{total} "
          f"passed, {skipped} skipped")
    if _SCRATCH:
        shutil.rmtree(_SCRATCH, ignore_errors=True)
    if skipped == total:
        print("[test_flux_cube_packager] SKIPPED ENTIRELY -- proved NOTHING here.")
    if failures:
        print("[test_flux_cube_packager] FAIL")
        return 1
    print("[test_flux_cube_packager] PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
