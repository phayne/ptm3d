#!/usr/bin/env python3
"""scripts/fetch_data.py against file:// URLs, and the Europa terrain manifest entry.

No network and no downloads: each arm builds a tiny archive in a temporary
directory, writes a manifest that points at it with a file:// URL, and runs the
real script as a subprocess. What is pinned:

* a dataset published in its OWN record (a per-dataset "base_url") is fetched
  from there, verified and unpacked where its "unpack_to" says;
* --base-url (a staging copy) still overrides it;
* a wrong checksum is refused and nothing is unpacked;
* a placeholder record is refused with a message, never a silent skip;
* the shipped "europa-terrain" entry is well formed and unpacks to exactly the
  directory every configs/europa_*.yaml reads its mesh and view factors from.

Optional arm: with PTM_EUROPA_TERRAIN_STAGE=<dir holding the staged
europa_terrain_socet_v1.tar.gz>, the real archive is fetched through file:// and
every file the configurations name must arrive. Without it the arm SKIPS and
says so.

Run directly::

    python3 src/python/tests/test_fetch_data.py
"""

from __future__ import annotations

import hashlib
import io
import json
import os
import re
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
FETCH = REPO_ROOT / "scripts" / "fetch_data.py"
MANIFEST = REPO_ROOT / "data" / "data_manifest.json"


class Skip(Exception):
    """An arm that cannot run here, as distinct from one that failed."""


def _archive(path, files):
    """A gzip tar of {name: bytes}; returns (sha256, size)."""
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tf:
        for name, data in files.items():
            ti = tarfile.TarInfo(name)
            ti.size = len(data)
            tf.addfile(ti, io.BytesIO(data))
    path.write_bytes(buf.getvalue())
    return hashlib.sha256(buf.getvalue()).hexdigest(), len(buf.getvalue())


def _manifest(path, dataset, manifest_base="file:///nonexistent/manifest-level"):
    path.write_text(json.dumps({"manifest_version": 1, "base_url": manifest_base,
                                "datasets": {"demo": dataset}}, indent=2))


def _run(*args):
    r = subprocess.run([sys.executable, str(FETCH), *args], capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def _setup(tmp, sha=None, base_url=None):
    tmp = Path(tmp)
    (tmp / "record").mkdir()
    real_sha, size = _archive(tmp / "record" / "demo.tar.gz",
                              {"a.ply": b"mesh bytes", "a.vf": b"view factor bytes"})
    ds = {"description": "test", "tier": "optional", "license": "CC-BY-4.0",
          "base_url": base_url if base_url is not None else (tmp / "record").as_uri(),
          "files": [{"name": "demo.tar.gz", "sha256": sha or real_sha, "size": size,
                     "archive": True, "unpack_to": "data/external/demo_set/"}]}
    _manifest(tmp / "manifest.json", ds)
    return tmp


def test_a_dataset_with_its_own_record_is_fetched_from_it():
    with tempfile.TemporaryDirectory() as t:
        tmp = _setup(t)
        rc, out = _run("--manifest", str(tmp / "manifest.json"), "--dataset", "demo",
                       "--dest", str(tmp / "data"))
        assert rc == 0, out
        target = tmp / "data" / "external" / "demo_set"
        assert (target / "a.ply").read_bytes() == b"mesh bytes", out
        assert (target / "a.vf").read_bytes() == b"view factor bytes", out
        rc, out = _run("--manifest", str(tmp / "manifest.json"), "--dataset", "demo",
                       "--dest", str(tmp / "data"), "--verify")
        assert rc == 0 and "unpacked content OK" in out, out


def test_base_url_on_the_command_line_overrides_the_dataset_record():
    with tempfile.TemporaryDirectory() as t:
        tmp = _setup(t, base_url="https://zenodo.org/records/ZENODO_PLACEHOLDER_X/files")
        rc, out = _run("--manifest", str(tmp / "manifest.json"), "--dataset", "demo",
                       "--dest", str(tmp / "data"), "--base-url", (tmp / "record").as_uri())
        assert rc == 0, out
        assert (tmp / "data" / "external" / "demo_set" / "a.ply").exists(), out


def test_a_wrong_checksum_is_refused_and_nothing_is_unpacked():
    with tempfile.TemporaryDirectory() as t:
        tmp = _setup(t, sha="0" * 64)
        rc, out = _run("--manifest", str(tmp / "manifest.json"), "--dataset", "demo",
                       "--dest", str(tmp / "data"))
        assert rc != 0, "a checksum mismatch was accepted:\n" + out
        assert not (tmp / "data" / "external" / "demo_set" / "a.ply").exists(), (
            "files were unpacked from an archive that failed its checksum")


def test_a_placeholder_record_is_refused_out_loud():
    with tempfile.TemporaryDirectory() as t:
        tmp = _setup(t, base_url="https://zenodo.org/records/ZENODO_PLACEHOLDER_X/files")
        rc, out = _run("--manifest", str(tmp / "manifest.json"), "--dataset", "demo",
                       "--dest", str(tmp / "data"))
        assert rc != 0, "a placeholder record was fetched:\n" + out
        assert "not published yet" in out, out


def _europa_needs():
    """{file name} every configs/europa_*.yaml reads under data/external/sutton_DEMs/."""
    need = set()
    for p in (REPO_ROOT / "configs").glob("europa_*.yaml"):
        need |= set(re.findall(r"data/external/sutton_DEMs/([\w.]+)", p.read_text()))
    return need


def test_the_europa_terrain_entry_is_well_formed():
    m = json.loads(MANIFEST.read_text())
    ds = m["datasets"].get("europa-terrain")
    assert ds, "data_manifest.json has no europa-terrain entry"
    assert ds["license"] == "CC-BY-4.0", ds["license"]
    assert "16904632" in ds["attribution"], "the source dataset is not credited"
    (fe,) = ds["files"]
    assert re.fullmatch(r"[0-9a-f]{64}", fe["sha256"]), fe["sha256"]
    assert fe["size"] > 100_000_000 and fe["archive"] is True, fe
    assert fe["unpack_to"] == "data/external/sutton_DEMs/", fe["unpack_to"]
    need = _europa_needs()
    assert len(need) == 10, f"expected 5 meshes + 5 view factors, configs name {sorted(need)}"
    # The record id appears in both the DOI and the base URL, and they must agree.
    rec_doi = ds["doi"].rsplit("zenodo.", 1)[1]
    rec_url = ds["base_url"].split("/records/")[1].split("/")[0]
    assert rec_doi == rec_url, f"DOI record {rec_doi!r} != base_url record {rec_url!r}"


def test_the_staged_europa_archive_delivers_every_configured_file():
    stage = os.environ.get("PTM_EUROPA_TERRAIN_STAGE", "").strip()
    if not stage:
        raise Skip("set PTM_EUROPA_TERRAIN_STAGE to the directory holding the staged "
                   "europa_terrain_socet_v1.tar.gz to run this arm")
    with tempfile.TemporaryDirectory() as t:
        rc, out = _run("--dataset", "europa-terrain", "--dest", t,
                       "--base-url", Path(stage).resolve().as_uri())
        assert rc == 0, out
        have = set(os.listdir(Path(t) / "external" / "sutton_DEMs"))
        missing = sorted(_europa_needs() - have)
        assert not missing, f"the archive lacks files the configs name: {missing}"


TESTS = [
    test_a_dataset_with_its_own_record_is_fetched_from_it,
    test_base_url_on_the_command_line_overrides_the_dataset_record,
    test_a_wrong_checksum_is_refused_and_nothing_is_unpacked,
    test_a_placeholder_record_is_refused_out_loud,
    test_the_europa_terrain_entry_is_well_formed,
    test_the_staged_europa_archive_delivers_every_configured_file,
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
    ran = len(TESTS) - skipped
    print(f"\n[test_fetch_data] {ran - failures}/{ran} passed, {skipped} skipped")
    if failures:
        print("[test_fetch_data] FAIL")
        return 1
    print("[test_fetch_data] PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
