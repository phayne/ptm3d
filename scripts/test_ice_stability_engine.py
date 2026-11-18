#!/usr/bin/env python3
"""Ice stability end to end: production defaults, the criterion itself, and a chunked run.

Complements ``test_ice_depth_coupled_radiation.py`` (which pins the spin-up so
that the coupled and uncoupled whole-mesh solves are the same integration) and
the analytic unit gate ``test_ice_stability_criterion``. The case is the public
synthetic bowl crater ``data/meshes/lunar_crater_80N.ply`` (1 km across, 200 m
deep, at 80N), on the CPU, one lunation at 118 steps.

Runs:
  off       radiation exchange off, every other setting at its default;
  coupled   radiation exchange on, every other setting at its default (seed
            anchor, seed indirect and the eq-indirect pass are NOT pinned);
  chunk_000/001  the coupled case as two chunks of 59 steps, built the way
            scripts/run_chunked_simulation.py builds them (start-time offset,
            warm start from the previous chunk's final state, look-ahead 1),
            with thermal.ice_stability.flux_threshold_kgm2s set to the
            literature value 1 kg m^-2 Gyr^-1, then merged by ptm_merge_chunks
            exactly as the runner calls it (no --J0).

Gates:
  D1  defaults, sign: coupled vs off, no face shallower, >= 100 faces deeper,
      and no SURFACE-node <E> sum lower (the eq-indirect spin-up may lower a
      few deep-node sums by < 5%; that count is printed, not gated);
  D2  the criterion, on every face of both runs, checked from the saved sums
      against Schorghofer & Taylor (2007) Eq. 11 with J0 and l read from the
      config: surface-stable exactly when mu <E(0)> <= J0; an at-depth face's
      depth brackets the root of z = (mu l / 2 J0) <E(z)> (log-linear <E>)
      within the solver's bisection precision; a censored face has
      J(z_i) > J0 at every node;
  D3  provenance.txt states the threshold in physical units;
  C1  merged chunks: ice_depth.bin is the D2 solution of the summed sums AT THE
      CONFIGURED threshold (and differs from the 3e-14 solution, so the gate
      is not vacuous), with the summed sample count;
  C2  each chunk records radiation-coupled ice temperatures and its threshold;
  C3  ptm_merge_chunks --J0 3e-14 writes the second layer from the same sums,
      warns that it overrides the chunks' threshold, and matches D2 at 3e-14;
  C4  ptm_merge_chunks refuses chunks recorded with different thresholds.

Needs the Moon SPICE kernels (python3 scripts/fetch_kernels.py --bodies moon),
build/bin/ptm_thermal_sim and build/bin/ptm_merge_chunks. Standard library only.

  python3 scripts/test_ice_stability_engine.py [--threads N] [--keep DIR]
"""
import argparse
import json
import math
import os
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
MESH = REPO / "data/meshes/lunar_crater_80N.ply"
STEPS, DAYS = 118, 29.53
MU, ELL = 2.99e-26, 75e-6
J_CODE = 3.0e-14
J_LIT = 1.0 / (1e9 * 365.25 * 86400.0)        # 1 kg m^-2 Gyr^-1
TAG = "[test_ice_stability_engine]"


def config(work, label, threads, radiation, steps=STEPS, days=DAYS, offset_s=None,
           initial=None, J0=None):
    th = {"run_thermal_mesh": True, "threads": threads, "equilibration_orbits": 1.0,
          "equilibration_steps_per_day": -1.0, "use_gpu": False,
          "compute_statistics": True, "enable_shadow_flux": False,
          "crater_dd_ratio": 8.0, "geothermal_heat_flow": -1.0,
          "compute_ice_stability": True,
          "final_state_file": str(work / (label + "_final.ptmt"))}
    sim = {"start_time": "2026-04-05T00:00:00", "duration_days": days,
           "timesteps": steps, "body": "Moon"}
    if offset_s is not None:                     # run_chunked_simulation.py:469-484
        sim["start_time_offset_seconds"] = float(offset_s)
        th["checkpoint_lookahead_steps"] = 1
    if initial is not None:
        th["initial_state_file"] = str(initial)
    if J0 is not None:
        th["ice_stability"] = {"flux_threshold_kgm2s": J0}
    cfg = {"simulation": sim,
           "mesh": {"import_file": str(MESH), "scale_factor": 1.0, "use_topographic": False},
           "thermal": th,
           "thermophysical": {"use_default_albedo": True, "incidence_angle_albedo": True},
           "radiation": {"enabled": False},
           "illumination": {"include_primary_reflection": True, "use_batch_processing": True,
                            "use_albedo_map": False},
           "output": {"directory": str(work / label), "format": "binary",
                      "write_raw_stats": True},
           "advanced": {"verbose": False}}
    if radiation:
        cfg["radiation"] = {"enabled": True, "use_gpu": False, "use_hmatrix": True,
                            "use_aca": True, "aca_tolerance": 0.001, "hmatrix_threads": threads,
                            "view_factor_file": str(work / "crater80N.vf"),
                            "max_coupling_iterations": 5, "coupling_temp_threshold": 2.0,
                            "coupling_flux_threshold": 0.1, "coupling_face_fraction": 0.01,
                            "band_separated": True}
    path = work / (label + ".yaml")
    path.write_text(json.dumps(cfg, indent=1))           # JSON is valid YAML
    return path


def run(cmd, log, expect_ok=True):
    with open(log, "w") as fh:
        r = subprocess.run([str(c) for c in cmd], cwd=REPO, stdout=fh, stderr=subprocess.STDOUT)
    if expect_ok and r.returncode != 0:
        sys.stdout.write(Path(log).read_text()[-4000:])
        raise SystemExit("%s FAIL: %s exited %d (log %s; the work directory is kept)"
                         % (TAG, Path(cmd[0]).name, r.returncode, log))
    return r.returncode, Path(log).read_text()


def read_bin(path):
    b = Path(path).read_bytes()
    n = struct.unpack_from("<I", b, 0)[0]
    return struct.unpack_from("<%df" % n, b, 4)


def read_raw(rundir):
    raw = Path(rundir) / "statistics" / "_raw"
    b = (raw / "ice_depth_sum_E.bin").read_bytes()
    nl, nf = struct.unpack_from("<II", b, 0)
    sum_e = struct.unpack_from("<%df" % (nl * nf), b, 8)
    n_samples = int((raw / "ice_depth_n_samples.txt").read_text().split()[0])
    g = (raw / "depth_grid.bin").read_bytes()
    z = struct.unpack_from("<%dd" % struct.unpack_from("<I", g, 0)[0], g, 4)
    return nl, nf, sum_e, n_samples, z


def criterion_check(z_eng, nl, nf, sum_e, n, z, J0):
    """Every face against S&T Eq. 11 with the surface rule, independently of the
    engine's bisection. Returns (n_bad, classes {surface, depth, censored}, worst
    depth error in units of the tolerance)."""
    pref = 0.5 * MU * ELL / J0
    bad, cls, worst = 0, [0, 0, 0], 0.0
    zcap = float(struct.unpack("<f", struct.pack("<f", z[-1]))[0])
    for f in range(nf):
        E = [sum_e[j * nf + f] / n for j in range(nl)]
        ze = z_eng[f]
        surface = not (E[0] > 0) or MU * E[0] <= J0
        if surface:
            cls[0] += 1
            bad += ze != 0.0
            continue
        cross = next((i for i in range(1, nl) if z[i] >= (pref * E[i] if E[i] > 0 else 0.0)), None)
        if cross is None:
            cls[2] += 1
            bad += ze != zcap
            continue
        cls[1] += 1
        i = cross
        lo, hi = z[i - 1], z[i]
        l0, l1 = math.log(max(E[i - 1], 1e-300)), math.log(max(E[i], 1e-300))

        def g(x):
            t = (x - lo) / (hi - lo)
            return x - pref * math.exp(l0 + t * (l1 - l0))
        tol = (hi - lo) / 512 * 1.01 + 1e-7
        a, b = max(lo, ze - tol), min(hi, ze + tol)
        ok = g(a) <= 0 <= g(b) or abs(g(ze)) <= 1e-12
        if not ok:
            # distance to the true root, in tolerances (bisection to 1e-12)
            x0, x1 = lo, hi
            for _ in range(80):
                xm = 0.5 * (x0 + x1)
                x0, x1 = (xm, x1) if g(xm) < 0 else (x0, xm)
            worst = max(worst, abs(ze - x0) / tol)
            bad += 1
    return bad, cls, worst


def provenance(rundir, key):
    for ln in (Path(rundir) / "provenance.txt").read_text().splitlines():
        if ln.startswith(key + "="):
            return ln.split("=", 1)[1]
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--engine", default=str(REPO / "build/bin/ptm_thermal_sim"))
    ap.add_argument("--merge", default=str(REPO / "build/bin/ptm_merge_chunks"))
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--keep", default=None)
    a = ap.parse_args()
    for exe in (a.engine, a.merge):
        if not os.access(exe, os.X_OK):
            raise SystemExit("%s FAIL: %s is not built" % (TAG, exe))
    if not (REPO / "data/kernels/spk/de440s.bsp").exists():
        raise SystemExit("%s FAIL: Moon SPICE kernels missing; run python3 "
                         "scripts/fetch_kernels.py --bodies moon" % TAG)
    work = Path(a.keep).resolve() if a.keep else Path(tempfile.mkdtemp(prefix="ptm_ice_eng_"))
    work.mkdir(parents=True, exist_ok=True)
    os.environ["OMP_NUM_THREADS"] = str(a.threads)

    fails = []

    def gate(name, ok, detail):
        print("  %-3s %-78s %s" % (name, detail, "PASS" if ok else "FAIL"), flush=True)
        if not ok:
            fails.append(name)

    # ---- defaults: off vs coupled -------------------------------------------
    run([a.engine, "--config", config(work, "off", a.threads, False)], work / "off.log")
    run([a.engine, "--config", config(work, "coupled", a.threads, True)], work / "coupled.log")
    z_off, z_on = read_bin(work / "off/statistics/ice_depth.bin"), read_bin(work / "coupled/statistics/ice_depth.bin")
    nl, nf, e_off, n_off, zg = read_raw(work / "off")
    _, _, e_on, n_on, _ = read_raw(work / "coupled")
    n_shallow = sum(1 for x, y in zip(z_off, z_on) if y < x)
    n_deep = sum(1 for x, y in zip(z_off, z_on) if y > x)
    surf_lower = sum(1 for f in range(nf) if e_on[f] < e_off[f] * (1 - 1e-5))
    deep_lower = sum(1 for k in range(nf, nl * nf) if e_on[k] < e_off[k] * (1 - 1e-5))
    gate("D1", n_shallow == 0 and n_deep >= 100 and surf_lower == 0,
         "defaults: %d deeper, %d shallower, %d surface sums lower (%d deep-node sums lower)"
         % (n_deep, n_shallow, surf_lower, deep_lower))
    for label, zz, ee, nn in (("off", z_off, e_off, n_off), ("coupled", z_on, e_on, n_on)):
        bad, cls, worst = criterion_check(zz, nl, nf, ee, nn, zg, J_CODE)
        gate("D2", bad == 0 and cls[0] + cls[1] > 0,
             "%s: Eq. 11 + surface rule at 3e-14, %d bad faces (surface/depth/censored %d/%d/%d)"
             % (label, bad, *cls) + ("" if bad == 0 else ", worst %.1f tol" % worst))
    th = provenance(work / "coupled", "ice_depth_threshold")
    gate("D3", "kg m^-2 Gyr^-1" in th and "of ice per Gyr" in th and "bare ice at 114.4 K" in th,
         "provenance: '%s'" % th[:70])

    # ---- chunked, literature threshold --------------------------------------
    half, days_half = STEPS // 2, DAYS * (STEPS // 2) / STEPS
    cadence_s = DAYS * 86400.0 / STEPS
    c0 = config(work, "chunk_000", a.threads, True, steps=half, days=days_half, offset_s=0.0, J0=J_LIT)
    run([a.engine, "--config", c0], work / "chunk_000.log")
    c1 = config(work, "chunk_001", a.threads, True, steps=half, days=days_half,
                offset_s=half * cadence_s, initial=work / "chunk_000_final.ptmt", J0=J_LIT)
    run([a.engine, "--config", c1], work / "chunk_001.log")
    rc, mlog = run([a.merge, "--chunks", work / "chunk_000", work / "chunk_001", "--output",
                    work / "merged"], work / "merge.log")
    z_m = read_bin(work / "merged/statistics/ice_depth.bin")
    raws = [read_raw(work / c) for c in ("chunk_000", "chunk_001")]
    sum_m = [x + y for x, y in zip(raws[0][2], raws[1][2])]
    n_m = raws[0][3] + raws[1][3]
    bad, cls, worst = criterion_check(z_m, nl, nf, sum_m, n_m, zg, J_LIT)
    bad_code, _, _ = criterion_check(z_m, nl, nf, sum_m, n_m, zg, J_CODE)
    thr = (work / "merged/statistics/ice_depth_threshold.txt").read_text()
    gate("C1", bad == 0 and bad_code > 0 and "3.17e-17" in thr and ("total samples = %d" % n_m) in mlog,
         "merged at the configured 1 kg/m2/Gyr: %d bad (%d/%d/%d), %d differ from 3e-14"
         % (bad, *cls, bad_code))
    ok2 = all(provenance(work / c, "ice_depth_temperatures").startswith("radiation_coupled")
              and "flux_threshold_kgm2s 3.16" in (work / c / "statistics/_raw/ice_depth_params.txt").read_text()
              for c in ("chunk_000", "chunk_001"))
    gate("C2", ok2, "chunks: radiation-coupled ice temperatures, threshold recorded")
    rc, olog = run([a.merge, "--chunks", work / "chunk_000", work / "chunk_001", "--output",
                    work / "merged_3e-14", "--J0", "3e-14"], work / "merge_override.log")
    z_o = read_bin(work / "merged_3e-14/statistics/ice_depth.bin")
    bad_o, cls_o, _ = criterion_check(z_o, nl, nf, sum_m, n_m, zg, J_CODE)
    gate("C3", bad_o == 0 and "OVERRIDES" in olog,
         "second layer from the same sums (--J0 3e-14): %d bad (%d/%d/%d), warned"
         % (bad_o, *cls_o))
    mixed = work / "mixed" / "chunk_001"
    if mixed.exists():
        shutil.rmtree(mixed)
    shutil.copytree(work / "chunk_001", mixed)
    (mixed / "statistics/_raw/ice_depth_params.txt").write_text(
        "flux_threshold_kgm2s 3e-14\ngrain_hop_length_m 7.5e-05\n")
    rc, xlog = run([a.merge, "--chunks", work / "chunk_000", mixed, "--output", work / "merged_mixed"],
                   work / "merge_mixed.log", expect_ok=False)
    gate("C4", rc != 0 and "different ice criteria" in xlog,
         "chunks with different thresholds refused (exit %d)" % rc)

    if fails:
        print("\n%s FAIL: %s (kept in %s)" % (TAG, ", ".join(fails), work))
        return 1
    if not a.keep:
        shutil.rmtree(work, ignore_errors=True)
    print("\n%s PASS" % TAG)
    return 0


if __name__ == "__main__":
    sys.exit(main())
