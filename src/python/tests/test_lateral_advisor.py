#!/usr/bin/env python3
"""The Python lateral-conduction advisor (ptm.lateral_advisor) against the paper and the engine.

Checks, in the order they would break:

1. D*(T) reproduces the Moon3D paper's values for lunar regolith (Hayne et al.,
   in revision, Section 2) and the same pinned numbers as the engine's
   ``test_lateral_advisor`` (0.8650 mm, 6.5289 mm, 2.0493 cm, 4.3531 cm).
2. The thresholds, representative temperatures and Martinez coefficients equal
   the constants in the C++ headers (re-read from source, so the two advisors
   cannot drift apart silently).
3. Verdicts on synthetic meshes: millimeter -> RECOMMENDED, 250 m -> not
   needed, 1 m -> not needed (rock caveat), a 0.5%-area fine patch ->
   BORDERLINE, empty -> UNMEASURED.
4. advise_config end to end on a PLY written to a temporary directory, a
   procedural icosphere config, lateral-ON wording, constant thermal inertia,
   and an unrecognized conductivity law (a NOTE, never silence).

Run directly (no pytest needed)::

    python3 src/python/tests/test_lateral_advisor.py
"""

from __future__ import annotations

import math
import re
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "src" / "python"))

import numpy as np  # noqa: E402

from ptm import lateral_advisor as la  # noqa: E402
from ptm.body_table import BODY_CAPABILITIES  # noqa: E402


def _moon():
    p = BODY_CAPABILITIES["Moon"]
    return la.SurfaceMaterial(ks=p["ks"], rhos=p["rhos"], cp0=p["cp0"],
                              cp_coeff=tuple(p["cpCoeff"]), emissivity=p["emissivity"],
                              chi=2.7, day_s=p["solar_day_days"] * 86400.0)


def _grid(L, s, x0=0.0):
    n = int(round(L / s))
    xs = x0 + np.arange(n + 1) * s
    X, Y = np.meshgrid(xs, np.arange(n + 1) * s)
    V = np.column_stack([X.ravel(), Y.ravel(), np.zeros(X.size)])
    i, j = np.meshgrid(np.arange(n), np.arange(n))
    a = (j * (n + 1) + i).ravel()
    F = np.vstack([np.column_stack([a, a + 1, a + n + 2]),
                   np.column_stack([a, a + n + 2, a + n + 1])])
    return V, F


def _verdict(V, F, m):
    return la.evaluate(*la.facet_sizes(V, F), m)


def test_paper_values():
    m = _moon()
    d = {T: la.critical_length(m, T) for T in (300.0, 110.0, 70.0, 50.0)}
    assert round(d[300.0] * 1e4) == 9, d          # 0.9 mm
    assert round(d[110.0] * 1e4) == 65, d         # 6.5 mm
    assert round(d[70.0] * 1e2) == 2, d           # 2 cm
    assert round(d[50.0] * 1e2) == 4, d           # 4 cm
    for T, pin in ((300.0, 8.650237e-4), (110.0, 6.528925e-3),
                   (70.0, 2.049285e-2), (50.0, 4.353098e-2)):
        assert abs(d[T] / pin - 1.0) < 1e-6, (T, d[T], pin)


def _header_constants(path, names):
    text = (REPO_ROOT / path).read_text()
    out = {}
    for name in names:
        hit = re.search(r"constexpr\s+double\s+" + name + r"\s*=\s*([-+0-9.eE]+)\s*;", text)
        assert hit, f"{path}: {name} not found"
        out[name] = float(hit.group(1))
    return out


def test_constants_match_engine():
    c = _header_constants("src/cpp/thermal/lateral_advisor.h", (
        "kAdvisorColdT", "kAdvisorIceT", "kAdvisorSunlitT", "kAdvisorRecommendAreaFraction",
        "kAdvisorBorderlineFloor", "kAdvisorBorderlineFactor", "kAdvisorRockCaveatFactor"))
    assert c["kAdvisorColdT"] == la.ADVISOR_COLD_T
    assert c["kAdvisorIceT"] == la.ADVISOR_ICE_T
    assert c["kAdvisorSunlitT"] == la.ADVISOR_SUNLIT_T
    assert c["kAdvisorRecommendAreaFraction"] == la.RECOMMEND_AREA_FRACTION
    assert c["kAdvisorBorderlineFloor"] == la.BORDERLINE_FLOOR
    assert c["kAdvisorBorderlineFactor"] == la.BORDERLINE_FACTOR
    assert c["kAdvisorRockCaveatFactor"] == la.ROCK_CAVEAT_FACTOR
    names9 = ("kMartinezAc", "kMartinezBc", "kMartinezCc", "kMartinezDc", "kMartinezEc",
              "kMartinezFc", "kMartinezGc", "kMartinezHc", "kMartinezIc")
    names10 = ("kMartinezA1", "kMartinezA2", "kMartinezB1", "kMartinezB2")
    m = _header_constants("src/cpp/thermal/conductivity_model.h", names9 + names10 + ("kMartinezTMin",))
    assert tuple(m[n] for n in names9) == la._MARTINEZ_EQ9
    assert tuple(m[n] for n in names10) == la._MARTINEZ_EQ10
    assert m["kMartinezTMin"] == la._MARTINEZ_TMIN


def test_synthetic_verdicts():
    m = _moon()
    assert _verdict(*_grid(0.3, 2.2e-3), m).verdict == "RECOMMENDED"
    assert _verdict(*_grid(50000.0, 220.0), m).verdict == "not needed"
    r = _verdict(*_grid(300.0, 0.88), m)
    assert r.verdict == "not needed"
    assert any("In rock" in line for line in la.format_advice(r, m))
    V1, F1 = _grid(20.0, 0.88)
    V2, F2 = _grid(1.42, 0.0176, x0=25.0)
    V = np.vstack([V1, V2])
    F = np.vstack([F1, F2 + len(V1)])
    rb = _verdict(V, F, m)
    assert rb.verdict == "BORDERLINE", (rb.frac_cold, rb.frac_cold_x2)
    re_ = la.evaluate([], [], m)
    assert re_.verdict == "UNMEASURED"
    assert "did NOT run" in la.format_advice(re_, m)[0]


def _write_ply(path, V, F):
    with open(path, "w") as fh:
        fh.write("ply\nformat ascii 1.0\nelement vertex %d\nproperty float x\nproperty float y\n"
                 "property float z\nelement face %d\nproperty list uchar int vertex_indices\n"
                 "end_header\n" % (len(V), len(F)))
        for v in V:
            fh.write("%.9g %.9g %.9g\n" % tuple(v))
        for f in F:
            fh.write("3 %d %d %d\n" % tuple(f))


def test_advise_config_end_to_end():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        V, F = _grid(0.1, 2.2e-3)
        _write_ply(tmp / "fine.ply", V, F)
        cfg = {"simulation": {"body": "Moon"}, "mesh": {"import_file": "fine.ply"},
               "thermal": {"run_thermal_mesh": True}}
        import yaml
        (tmp / "fine.yaml").write_text(yaml.safe_dump(cfg))
        a = la.advise_config(tmp / "fine.yaml")
        assert a.verdict == "RECOMMENDED", a.lines
        assert a.lines[0].startswith("[lateral-advisor] RECOMMENDED")
        # local patch, magnitude alone says km: read as meters, other reading stated
        assert any("Read as kilometers instead" in line for line in a.lines), a.lines

        cfg["thermal"]["lateral_conduction"] = True
        (tmp / "on.yaml").write_text(yaml.safe_dump(cfg))
        assert "is ON, and appropriate" in la.advise_config(tmp / "on.yaml").lines[0]

        cfg["thermophysical"] = {"conductivity_model": "not_a_law"}
        (tmp / "bad.yaml").write_text(yaml.safe_dump(cfg))
        assert any("is not recognized" in line for line in la.advise_config(tmp / "bad.yaml").lines)

        cfg["thermophysical"] = {"constant_thermal_inertia": 50.0}
        mat, _ = la.surface_material_from_config(cfg)
        assert mat.chi == 0.0 and mat.cp_coeff == ()
        assert abs(mat.cp0 - 50.0 ** 2 / (mat.ks * mat.rhos)) < 1e-9

    ico = la.advise_config({"simulation": {"body": "Moon"}, "mesh": {"resolution": 6}})
    assert ico.verdict == "not needed", ico.lines
    assert "procedural icosphere" in ico.lines[3]


def test_cli():
    from ptm import cli
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        V, F = _grid(10.0, 1.0)
        _write_ply(tmp / "m.ply", V, F)
        (tmp / "c.yaml").write_text("simulation: {body: Moon}\nmesh: {import_file: m.ply}\n")
        assert cli.main(["advise", str(tmp / "c.yaml")]) == 0
        assert cli.main(["advise", str(tmp / "c.yaml"), "--json"]) == 0


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
