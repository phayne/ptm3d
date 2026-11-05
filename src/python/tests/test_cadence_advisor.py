#!/usr/bin/env python3
"""The Python output-cadence advisor (ptm.cadence_advisor) against the engine's header and a synthetic run.

1. Constants equal those in src/cpp/thermal/cadence_advisor.h (re-read from source).
2. theta/omega for the Moon is 1.05 h (Fermi 2026-11-04) and scales as 1/d; the recommended step is
   0.4 theta/omega; the missed-flash dT reproduces Fermi's Europa estimate (~16 K at 35.5 min).
3. postrun() on a synthetic run directory: R_ct = (area of ALL faces lit in 1 frame, or 1-2 frames when
   provenance.txt says the step is below theta/omega) / (area with T_max < thr); format_postrun() gives
   the engine's line format; the no-T_max form; `ptm cadence`.

Run directly (no pytest needed)::

    python3 src/python/tests/test_cadence_advisor.py
"""
from __future__ import annotations

import re
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "src" / "python"))

import numpy as np  # noqa: E402

from ptm import cadence_advisor as ca  # noqa: E402


def test_constants_match_header():
    text = (REPO_ROOT / "src/cpp/thermal/cadence_advisor.h").read_text()
    def c(name):
        m = re.search(r"constexpr\s+double\s+" + name + r"\s*=\s*([-+0-9.eE]+)\s*;", text)
        assert m, name
        return float(m.group(1))
    assert c("kRecC") == ca.REC_C
    assert c("kPeak5W") == ca.PEAK5_W
    assert c("kLitThresholdWm2") == ca.LIT_THRESHOLD_WM2


def test_theta_over_omega():
    th, om, to = ca.theta_over_omega(body="Moon")
    assert abs(th - 0.533) < 0.001 and abs(om - 0.508) < 0.001 and abs(to - 1.05) < 0.01, (th, om, to)
    _, _, to2 = ca.theta_over_omega(distance_au=2.0, solar_day_days=29.530588)
    assert abs(to2 / to - 0.5) < 1e-3   # theta ~ 1/d
    a, b = ca.recommended_step(th, om)
    assert abs(a - 0.4 * to) < 1e-12 and abs(b - 0.3 * to) < 1e-12
    F = 1361.0 * 0.4 / 5.2 ** 2
    assert abs(ca.missed_flash_delta_t(F, 65.0, 0.9, 110.0, 35.5 * 60) - 16.0) < 1.5


def _write(path, arr):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(np.uint32(len(arr)).tobytes())
        fh.write(np.asarray(arr, dtype="<f4").tobytes())


def _run_dir(tmp, with_tmax=True):
    # 4 unit-area triangles (two squares); face k lit in k frames over 3 frames
    V = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [2, 0, 0], [2, 1, 0]], float) * np.sqrt(2.0)
    F = np.array([[0, 1, 2], [0, 2, 3], [1, 4, 5], [1, 5, 2]])
    with open(tmp / "mesh.ply", "w") as fh:
        fh.write("ply\nformat ascii 1.0\nelement vertex 6\nproperty float x\nproperty float y\nproperty float z\n"
                 "element face 4\nproperty list uchar int vertex_indices\nend_header\n")
        for v in V:
            fh.write("%.9g %.9g %.9g\n" % tuple(v))
        for f in F:
            fh.write("3 %d %d %d\n" % tuple(f))
    lit = np.array([[0, 1, 1, 0], [0, 0, 1, 1], [0, 0, 1, 1]])   # frames x faces: counts 0,1,3,2
    for t, row in enumerate(lit):
        _write(tmp / "direct_solar" / f"direct_solar_{t:06d}.bin", row * 100.0)
    if with_tmax:
        _write(tmp / "statistics" / "temperature_max.bin", [40.0, 50.0, 300.0, 100.0])


def test_postrun_synthetic():
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        _run_dir(tmp)          # faces lit in 0, 1, 3, 2 frames; T_max 40, 50, 300, 100 K; equal areas
        r = ca.postrun(tmp)
        assert r["n_frames"] == 3 and r["max_frames_lit"] == 1, r
        assert abs(r["indicator_area_frac"] - 0.25) < 1e-12, r                 # face 1
        assert abs(r["ct_area_frac_110K"] - 0.75) < 1e-12                       # faces 0, 1, 3
        assert abs(r["ct_overstatement_bound_110K"] - 0.25 / 0.75) < 1e-12
        assert abs(r["ct_overstatement_bound_54K"] - 0.25 / 0.5) < 1e-12        # faces 0, 1
        line = ca.format_postrun(r)
        assert line.startswith("[cadence-advisor] post-run: estimated upper bound on the fractional overstatement "
                               "of the cold-trap area = (area of all faces lit in exactly one output frame, 25% of "
                               "the mesh) / (area with frame-sampled T_max below the threshold): 33.3% at 110 K "
                               "(cold traps 75% of the mesh), 50% at 54 K (50%)."), line
        # a step below theta/omega in provenance.txt switches the count to one or two frames (faces 1 and 3)
        (tmp / "provenance.txt").write_text("cadence_dt_out_s=1500\ncadence_sun_disk_deg=0.533\n"
                                            "cadence_sun_rate_deg_per_h=0.508\n")
        r2 = ca.postrun(tmp)
        assert r2["max_frames_lit"] == 2 and abs(r2["indicator_area_frac"] - 0.5) < 1e-12, r2
        assert "one or two output frames" in ca.format_postrun(r2)
        from ptm import cli
        assert cli.main(["cadence", str(tmp)]) == 0
        assert cli.main(["cadence", str(tmp), "--json"]) == 0
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        _run_dir(tmp, with_tmax=False)
        line = ca.format_postrun(ca.postrun(tmp))
        assert line.endswith("25% of the mesh area (no per-face T_max for this run, so no cold-trap estimate)."), line


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS {t.__name__}")
        except Exception as exc:  # noqa: BLE001
            failed += 1
            print(f"FAIL {t.__name__}: {type(exc).__name__}: {exc}")
    print(f"{len(tests) - failed}/{len(tests)} passed" + ("" if not failed else " -- FAIL"))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
