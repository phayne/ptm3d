#!/usr/bin/env python3
"""What an ice file in statistics/ is: verified from its record, never assumed.

``statistics/ice_depth.bin`` was the flux-threshold depth before the retention
criterion and is the retention depth since. A run into a reused directory, or
a merge into one, could leave one next to the other's record, and every reader
presented it as the retention depth. ``ptm.ice_products`` decides from
``ice_retention.txt`` and the bytes on disk. Both directions are asserted:

1. a retention product whose files match the record's checksums is accepted
   as such by ``ptm.ice_products``, ``Results.stat`` (no warning), the site
   lookup script and the viewer label;
2. an ``ice_depth.bin`` with no record (an engine before the retention
   criterion) is the legacy flux depth: ``Results.stat`` warns and says so,
   the lookup and the label call it the flux-threshold depth;
3. a file that does not match the record (a stale flux file next to a newer
   record), a record without checksums, and a companion without a record are
   refused: ``Results.stat`` raises, the lookup withholds the value, the label
   says "Criterion Unverified".

No engine, no GPU. Run directly::

    python3 src/python/tests/test_ice_products.py
"""

from __future__ import annotations

import shutil
import struct
import subprocess
import sys
import tempfile
import warnings
import zlib
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "src" / "python"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np                                              # noqa: E402

from ptm import ice_products                                    # noqa: E402
from ptm.results import Results                                 # noqa: E402
from test_results_api import _icosphere, make_run               # noqa: E402

LOOKUP = REPO_ROOT / "scripts" / "property_maps" / "lookup_face_at_site.py"


def _bin_bytes(values):
    values = np.asarray(values, dtype="<f4")
    return struct.pack("<I", values.size) + values.tobytes()


def write_retention(stats, n, depth=0.05, checksums=True):
    """The six retention files and their record, as the engine writes them."""
    stats.mkdir(parents=True, exist_ok=True)
    files = {
        "ice_depth.bin": np.full(n, depth),
        "ice_table_layer.bin": np.full(n, 5.0),
        "ice_retention_class.bin": np.full(n, 1.0),
        "ice_toplayer_lifetime.bin": np.full(n, 0.2),
        "ice_column_lifetime.bin": np.full(n, 0.9),
        "ice_exposed_loss_tau.bin": np.full(n, 1e-3),
    }
    lines = ["criterion: ice retention after an age (retention-v2 (2026-11-18))",
             "age_gyr: 1", "files: ..."]
    for name, vals in files.items():
        b = _bin_bytes(vals)
        (stats / name).write_bytes(b)
        if checksums:
            lines.append("file_crc32: %s %08x" % (name, zlib.crc32(b)))
    (stats / "ice_retention.txt").write_text("\n".join(lines) + "\n")


def _plain_ply(path):
    """The input-mesh layout the lookup script reads: float32 xyz, faces with
    no extra properties (it cannot read the engine's mesh.ply, whose faces
    carry latitude and longitude)."""
    v, f = _icosphere(2)
    v = v * 1737.4e3
    with open(path, "wb") as h:
        h.write(b"ply\nformat binary_little_endian 1.0\nelement vertex %d\n"
                b"property float x\nproperty float y\nproperty float z\n"
                b"element face %d\nproperty list uchar int vertex_indices\nend_header\n"
                % (len(v), len(f)))
        h.write(v.astype("<f4").tobytes())
        for a, b, c in f:
            h.write(struct.pack("<B3i", 3, a, b, c))
    return path


def lookup(run):
    mesh = _plain_ply(run / "plain.ply")
    r = subprocess.run([sys.executable, str(LOOKUP), "--mesh", str(mesh),
                        "--stats-dir", str(run / "statistics"), "--lat", "0", "--lon", "0"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    return next(ln for ln in r.stdout.splitlines() if ln.strip().startswith("ice_depth"))


def _gui_label():
    """visualization_panel.ice_field_label, or None without PyQt6."""
    try:
        sys.path.insert(0, str(REPO_ROOT / "gui"))
        from gui_panels.visualization_panel import ice_field_label
        return ice_field_label
    except Exception:                                       # noqa: BLE001
        return None


def test_verified_retention_is_accepted(tmp):
    run = make_run(tmp / "verified")
    res = Results(run)
    write_retention(run / "statistics", res.n_faces)
    for name in ice_products.RETENTION_FILES:
        p = ice_products.provenance(run / "statistics", name)
        assert p.criterion == ice_products.RETENTION, (name, p)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        v = res.stat("ice_depth", None)
    assert np.allclose(v, 0.05)
    assert res.ice_provenance().verified
    line = lookup(run)
    assert "ice table at depth after the retention age" in line and "0.05000 m" in line, line
    assert ice_products.display_label(run / "statistics", "ice_depth.bin", "Ice Table Depth",
                                      " (m)") == "Ice Table Depth (m)"
    lab = _gui_label()
    if lab is not None:
        assert lab(run, "ice_depth", "Ice Table Depth (m)") == "Ice Table Depth (m)"


def test_no_record_is_the_legacy_flux_depth(tmp):
    run = make_run(tmp / "legacy")
    res = Results(run)
    (run / "statistics" / "ice_depth.bin").write_bytes(_bin_bytes(np.full(res.n_faces, 0.02)))
    p = ice_products.provenance(run / "statistics")
    assert p.criterion == ice_products.LEGACY_FLUX, p
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        v = res.stat("ice_depth", None)
    assert np.allclose(v, 0.02)
    assert any("NOT the ice-table depth" in str(x.message) for x in w), [str(x.message) for x in w]
    line = lookup(run)
    assert "flux-threshold depth" in line and "0.02000 m" in line, line
    label = ice_products.display_label(run / "statistics", "ice_depth.bin", "Ice Table Depth", " (m)")
    assert label == "Ice Depth, Flux Criterion, Legacy Run (m)", label
    lab = _gui_label()
    if lab is not None:
        assert lab(run, "ice_depth", "Ice Table Depth (m)") == label


def test_stale_file_next_to_a_record_is_refused(tmp):
    run = make_run(tmp / "stale")
    res = Results(run)
    write_retention(run / "statistics", res.n_faces)
    # A later run with an engine from before the retention criterion writes
    # its flux ice_depth.bin over the retention one and leaves the record.
    (run / "statistics" / "ice_depth.bin").write_bytes(_bin_bytes(np.full(res.n_faces, 0.3)))
    p = ice_products.provenance(run / "statistics")
    assert p.criterion == ice_products.UNVERIFIED and "checksum" in p.reason, p
    try:
        res.stat("ice_depth", None)
    except ice_products.IceProvenanceError as exc:
        assert "unverified" in str(exc)
    else:
        raise AssertionError("Results.stat returned an ice_depth.bin that does not match its record")
    try:
        ice_products.require_retention(run / "statistics")
    except ice_products.IceProvenanceError:
        pass
    else:
        raise AssertionError("require_retention accepted a stale ice_depth.bin")
    line = lookup(run)
    assert "UNVERIFIED" in line and "withheld" in line and "0.30000" not in line, line
    label = ice_products.display_label(run / "statistics", "ice_depth.bin", "Ice Table Depth", " (m)")
    assert label == "Ice Table Depth (m) [Criterion Unverified]", label
    # The companions still match: they are accepted on their own.
    assert ice_products.provenance(run / "statistics", "ice_retention_class.bin").verified


def test_record_without_checksums_is_unverified(tmp):
    run = make_run(tmp / "nochecksum")
    res = Results(run)
    write_retention(run / "statistics", res.n_faces, checksums=False)
    for name in ice_products.RETENTION_FILES:
        assert ice_products.provenance(run / "statistics", name).criterion == ice_products.UNVERIFIED
    try:
        res.stat("ice_retention_class", None)
    except ice_products.IceProvenanceError:
        pass
    else:
        raise AssertionError("a record without checksums was trusted")


def test_companion_without_record_is_unverified(tmp):
    run = make_run(tmp / "orphan")
    res = Results(run)
    (run / "statistics" / "ice_retention_class.bin").write_bytes(_bin_bytes(np.ones(res.n_faces)))
    p = ice_products.provenance(run / "statistics", "ice_retention_class.bin")
    assert p.criterion == ice_products.UNVERIFIED, p
    assert ice_products.provenance(run / "statistics") is None    # no ice_depth.bin
    line = lookup(run)
    assert "missing" in line, line


def test_non_ice_fields_are_untouched(tmp):
    run = make_run(tmp / "plain")
    res = Results(run)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        res.stat("temperature", "max")
    assert ice_products.display_label(run / "statistics", "temperature_max.bin", "T", " (K)") == "T (K)"


TESTS = [
    test_verified_retention_is_accepted,
    test_no_record_is_the_legacy_flux_depth,
    test_stale_file_next_to_a_record_is_refused,
    test_record_without_checksums_is_unverified,
    test_companion_without_record_is_unverified,
    test_non_ice_fields_are_untouched,
]


def main():
    tmp = Path(tempfile.mkdtemp(prefix="ptm_ice_products_"))
    failures = 0
    try:
        for test in TESTS:
            try:
                test(tmp)
            except AssertionError as exc:
                failures += 1
                print(f"FAIL {test.__name__}\n     {exc}")
            except Exception as exc:                   # noqa: BLE001
                failures += 1
                print(f"ERROR {test.__name__}: {type(exc).__name__}: {exc}")
            else:
                print(f"PASS {test.__name__}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print(f"\n[test_ice_products] {len(TESTS) - failures}/{len(TESTS)} passed"
          f"{'' if _gui_label() else ' (viewer label checked through ptm.ice_products only: no PyQt6)'}")
    if failures:
        print("[test_ice_products] FAIL")
        return 1
    print("[test_ice_products] PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
