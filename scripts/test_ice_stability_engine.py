#!/usr/bin/env python3
"""Ice stability end to end: production defaults, both criteria, offline re-evaluation, chunks.

Complements ``test_ice_depth_coupled_radiation.py`` (which pins the spin-up so
that the coupled and uncoupled whole-mesh solves are the same integration) and
the analytic unit gates ``test_ice_retention`` (retention criterion) and
``test_ice_stability_criterion`` (flux criterion). The case is the public
synthetic bowl crater ``data/meshes/lunar_crater_80N.ply`` (1 km across, 200 m
deep, at 80N), on the CPU, one lunation at 118 steps.

Runs:
  off       radiation exchange off, every other setting at its default;
  coupled   radiation exchange on, every other setting at its default (seed
            anchor, seed indirect and the eq-indirect pass are NOT pinned);
  flux_only the off case with thermal.ice_stability.criterion: flux;
  chunk_000/001  the coupled case as two chunks of 59 steps, built the way
            scripts/run_chunked_simulation.py builds them (start-time offset,
            warm start from the previous chunk's final state, look-ahead 1),
            with thermal.ice_stability.flux_threshold_kgm2s set to the
            literature value 1 kg m^-2 Gyr^-1 (a config written for the old
            product), then merged by ptm_merge_chunks exactly as the runner
            calls it.
The flux-criterion products of off and coupled are computed from their saved
sums with ptm_merge_chunks --ice-only (the documented offline command).

Gates (retention criterion, the default product ice_depth.bin):
  D1  defaults, sign: coupled vs off, no face shallower, >= 100 faces deeper,
      and no SURFACE-node <E> sum lower;
  R1  every face of off and coupled against an independent port of the
      reference script (scratch retention.py of the specification: trapezoid
      in m, linear root in t): the same class on every face, ice-table depth
      within 1e-4 m, top-layer lifetime within 1e-3;
  R2  nesting: every face retained is stable under J0 = 3e-14 (1 m of exposed
      ice per Gyr), and every face stable under 1 kg m^-2 Gyr^-1 (1 mm per
      Gyr) is retained;
  R3  one implementation: ptm_merge_chunks --ice-only on the coupled run
      reproduces the engine's six retention files byte for byte;
  R4  offline re-evaluation at 3 Gyr: warns that it overrides the run, records
      the override, and no ice table is shallower than at 1 Gyr;
  D3  provenance.txt and statistics/ice_retention.txt state the criterion, age,
      inventory, density profile, transport length and grid.
Gates (flux criterion, legacy):
  D2  the flux product of off and coupled (offline, 3e-14) against Schorghofer
      & Taylor (2007) Eq. 11 with the surface rule, face by face;
  L1  criterion: flux writes ice_depth_flux.bin and no ice_depth.bin, says so,
      and equals the offline flux product from the same sums byte for byte.
Gates (chunks):
  C1  merged chunks: ice_depth.bin is the reference solution of the summed sums
      (class for class, 1e-4 m) and ice_depth_flux.bin the Eq. 11 solution at
      the CONFIGURED threshold (not 3e-14), with the summed sample count;
  C2  each chunk records radiation-coupled ice temperatures, its threshold and
      the retention settings, and its log carries the notice that
      flux_threshold_kgm2s now also writes ice_depth_flux.bin;
  C3  ptm_merge_chunks --J0 3e-14 writes the second flux layer from the same
      sums, warns that it overrides the chunks' threshold, and matches D2;
  C4  ptm_merge_chunks refuses chunks recorded with different thresholds;
  C5  coupled chunked vs direct: the merged retention product against the
      direct coupled run. Not an identity: each chunk's coupling re-thermal
      re-seeds its own spin-up (PTM-862), so classes must agree on >= 99% of
      faces and the numbers are printed;
  C7  chunks written before the retention criterion: the merge warns and still
      writes their flux product; mixed with new chunks and given the density
      profile, it writes both products, identical to the new-format merge;
  C6  chunked = direct where the two ARE the same integration: radiation off,
      both started from the same mid-window checkpoint, one 59-step run against
      29 + 30 steps merged; the same class on every face, depth within 1e-4 m,
      the same sample count.

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
           initial=None, J0=None, criterion=None):
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
    if J0 is not None or criterion is not None:
        th["ice_stability"] = {}
        if J0 is not None:
            th["ice_stability"]["flux_threshold_kgm2s"] = J0
        if criterion is not None:
            th["ice_stability"]["criterion"] = criterion
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


GYR = 1e9 * 365.25 * 86400.0
NSUB = 16


def read_record(rundir):
    rec = {}
    for ln in (Path(rundir) / "statistics/_raw/ice_depth_params.txt").read_text().splitlines():
        k, v = ln.split(None, 1)
        rec[k] = v.strip()
    return rec


def reference_retention(sum_e, n, nl, nf, z, rec, age_gyr=1.0):
    """Independent port of the specification's reference script (retention.py,
    2026-11-18): m and 1/l by the trapezoid rule at the sub-step points, the root
    by linear interpolation of t. Pore-filling inventory from the run's
    recorded density profile. Returns per face (depth, class 0/1/2, t_top Gyr)."""
    rs, rd, H = (float(rec[k]) for k in ("rho_surface_kgm3", "rho_deep_kgm3", "density_scale_m"))
    rg, ri, ell = float(rec["grain_density_kgm3"]), float(rec["ice_density_kgm3"]), float(rec["grain_hop_length_m"])
    zf, iv, sv = [], [], []
    for i in range(nl - 1):
        for k in range(NSUB):
            s = k / NSUB
            zf.append(z[i] + s * (z[i + 1] - z[i])); iv.append(i); sv.append(s)
    zf.append(z[-1]); iv.append(nl - 2); sv.append(1.0)
    m = [ri * (1.0 - (rd + (rs - rd) * math.exp(-x / H)) / rg) for x in zf]
    lam = [0.0]
    for k in range(1, len(zf)):
        lam.append(lam[-1] + (zf[k] - zf[k - 1]) / ell)
    w = [m[k] * (1 + 2 * lam[k]) / MU for k in range(len(zf))]
    zh, T = 0.5 * (z[0] + z[1]), age_gyr * GYR
    k_top = NSUB // 2
    out = []
    for f in range(nf):
        le = [math.log(max(sum_e[j * nf + f] / n, 1e-200)) for j in range(nl)]
        t, prev_g, zt, ttop = 0.0, None, None, None
        for k in range(len(zf)):
            g = w[k] * math.exp(-(le[iv[k]] + sv[k] * (le[iv[k] + 1] - le[iv[k]])))
            if k:
                t0 = t
                t += 0.5 * (g + prev_g) * (zf[k] - zf[k - 1])
                if zt is None and t >= T:
                    zt = zf[k - 1] + (T - t0) / max(t - t0, 1e-300) * (zf[k] - zf[k - 1])
            if k == k_top:
                ttop = t
            prev_g = g
        cls = 2 if zt is None else (0 if zt < zh else 1)
        out.append((zt if zt is not None else z[-1], cls, ttop / GYR))
    return out


def read_vec(path):
    return read_bin(path)


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
         "retention defaults: %d deeper, %d shallower, %d surface sums lower (%d deep-node sums lower)"
         % (n_deep, n_shallow, surf_lower, deep_lower))

    # ---- R1: engine vs the reference script, every face ---------------------
    ret = {}
    for label, ee, nn in (("off", e_off, n_off), ("coupled", e_on, n_on)):
        rec = read_record(work / label)
        refs = reference_retention(ee, nn, nl, nf, zg, rec)
        st = work / label / "statistics"
        z_e, c_e, tt_e = (read_bin(st / x) for x in ("ice_depth.bin", "ice_retention_class.bin",
                                                     "ice_toplayer_lifetime.bin"))
        ret[label] = (z_e, c_e)
        n_cls = sum(1 for f in range(nf) if (0 if c_e[f] == 3 else c_e[f]) != refs[f][1])
        dz = max(abs(z_e[f] - refs[f][0]) for f in range(nf) if refs[f][1] < 2)
        dt = max(abs(tt_e[f] / refs[f][2] - 1) for f in range(nf) if 0 < refs[f][2] < 1e30)
        cnt = [sum(1 for c in c_e if c == k) for k in (0, 1, 2, 3)]
        gate("R1", n_cls == 0 and dz <= 1e-4 and dt <= 1e-3,
             "%s vs reference: %d class diffs, max |dz| %.2g m, t_top %.1e (top/depth/none/uf %d/%d/%d/%d)"
             % (label, n_cls, dz, dt, *cnt))

    # ---- flux criterion offline from the saved sums (3e-14, 1 kg/m2/Gyr) ----
    def offline(label, out, *extra, expect_ok=True):
        return run([a.merge, "--ice-only", "--allow-nonstandard-names", "--chunks", work / label,
                    "--output", work / out, *extra], work / (out + ".log"), expect_ok)
    zf = {}
    for label in ("off", "coupled"):
        for tag, J in (("code", J_CODE), ("lit", J_LIT)):
            offline(label, "%s_flux_%s" % (label, tag), "--criterion", "flux", "--J0", "%.17g" % J)
            zf[label, tag] = read_bin(work / ("%s_flux_%s/statistics/ice_depth_flux.bin" % (label, tag)))
        for k, J in (("code", J_CODE), ("lit", J_LIT)):
            ee, nn = (e_off, n_off) if label == "off" else (e_on, n_on)
            bad, cls, worst = criterion_check(zf[label, k], nl, nf, ee, nn, zg, J)
            gate("D2", bad == 0 and cls[0] + cls[1] > 0,
                 "%s flux product at %.3g (offline): %d bad faces (surface/depth/censored %d/%d/%d)"
                 % (label, J, bad, *cls) + ("" if bad == 0 else ", worst %.1f tol" % worst))

    # ---- R2: nesting ---------------------------------------------------------
    zcap = float(struct.unpack("<f", struct.pack("<f", zg[-1]))[0])
    for label in ("off", "coupled"):
        _, c_e = ret[label]
        kept = [c in (0, 1, 3) for c in c_e]
        st_code = [z < zcap for z in zf[label, "code"]]
        st_lit = [z < zcap for z in zf[label, "lit"]]
        v1 = sum(1 for k, s_ in zip(kept, st_code) if k and not s_)
        v2 = sum(1 for k, s_ in zip(kept, st_lit) if s_ and not k)
        gate("R2", v1 == 0 and v2 == 0,
             "%s nesting: retained %d, stable at 1 m/Gyr %d, at 1 mm/Gyr %d; violations %d, %d"
             % (label, sum(kept), sum(st_code), sum(st_lit), v1, v2))

    # ---- R3: one implementation (offline == engine, byte for byte) ----------
    offline("coupled", "coupled_offline")
    names = ("ice_depth.bin", "ice_table_layer.bin", "ice_retention_class.bin",
             "ice_toplayer_lifetime.bin", "ice_column_lifetime.bin", "ice_exposed_loss_tau.bin")
    same = [n_ for n_ in names if (work / "coupled/statistics" / n_).read_bytes()
            == (work / "coupled_offline/statistics" / n_).read_bytes()]
    gate("R3", len(same) == len(names),
         "ptm_merge_chunks --ice-only reproduces the engine's retention files: %d of %d identical"
         % (len(same), len(names)))

    # ---- R4: re-evaluation at another age -----------------------------------
    rc, rlog = offline("coupled", "coupled_3gyr", "--retention-age-gyr", "3")
    z3 = read_bin(work / "coupled_3gyr/statistics/ice_depth.bin")
    c3 = read_bin(work / "coupled_3gyr/statistics/ice_retention_class.bin")
    shallower = sum(1 for x, y in zip(ret["coupled"][0], z3) if y < x - 1e-7)
    prov3 = (work / "coupled_3gyr/statistics/ice_retention.txt").read_text()
    gate("R4", "OVERRIDES" in rlog and "age_gyr: 3" in prov3 and "--retention-age-gyr" in prov3
         and shallower == 0,
         "3 Gyr offline: warned, recorded; retained %d (1 Gyr: %d), %d ice tables shallower"
         % (sum(1 for c in c3 if c != 2), sum(1 for c in ret["coupled"][1] if c != 2), shallower))

    # ---- D3: provenance --------------------------------------------------------
    pr = provenance(work / "coupled", "ice_retention")
    rt = (work / "coupled/statistics/ice_retention.txt").read_text()
    gate("D3", provenance(work / "coupled", "ice_criterion") == "retention"
         and "age_gyr: 1" in pr and "pore_filling" in pr and "75 um" in pr and "H = 0.06" in pr
         and "z_bot" in rt and "top layer 0 to" in rt,
         "provenance: criterion, age, inventory, density, transport length, grid")

    # ---- L1: the legacy criterion alone ------------------------------------
    rc, flog = run([a.engine, "--config", config(work, "flux_only", a.threads, False, criterion="flux")],
                   work / "flux_only.log")
    fz = work / "flux_only/statistics"
    eq_sums = read_raw(work / "flux_only")[2] == e_off
    gate("L1", not (fz / "ice_depth.bin").exists() and "criterion: flux (legacy" in flog
         and eq_sums and (fz / "ice_depth_flux.bin").read_bytes()
         == (work / "off_flux_code/statistics/ice_depth_flux.bin").read_bytes()
         and provenance(work / "flux_only", "ice_depth_flux_threshold").startswith("J0 = 3e-14"),
         "criterion: flux: ice_depth_flux.bin only, warned, = offline flux product from the same sums")

    # ---- chunked, literature threshold set (a config for the old product) ---
    half, days_half = STEPS // 2, DAYS * (STEPS // 2) / STEPS
    cadence_s = DAYS * 86400.0 / STEPS
    c0 = config(work, "chunk_000", a.threads, True, steps=half, days=days_half, offset_s=0.0, J0=J_LIT)
    _, l0 = run([a.engine, "--config", c0], work / "chunk_000.log")
    c1 = config(work, "chunk_001", a.threads, True, steps=half, days=days_half,
                offset_s=half * cadence_s, initial=work / "chunk_000_final.ptmt", J0=J_LIT)
    _, l1 = run([a.engine, "--config", c1], work / "chunk_001.log")
    rc, mlog = run([a.merge, "--chunks", work / "chunk_000", work / "chunk_001", "--output",
                    work / "merged"], work / "merge.log")
    raws = [read_raw(work / c) for c in ("chunk_000", "chunk_001")]
    sum_m = [x + y for x, y in zip(raws[0][2], raws[1][2])]
    n_m = raws[0][3] + raws[1][3]
    z_mf = read_bin(work / "merged/statistics/ice_depth_flux.bin")
    bad, cls, worst = criterion_check(z_mf, nl, nf, sum_m, n_m, zg, J_LIT)
    bad_code, _, _ = criterion_check(z_mf, nl, nf, sum_m, n_m, zg, J_CODE)
    thr = (work / "merged/statistics/ice_depth_flux_threshold.txt").read_text()
    refs_m = reference_retention(sum_m, n_m, nl, nf, zg, read_record(work / "chunk_000"))
    z_mr = read_bin(work / "merged/statistics/ice_depth.bin")
    c_mr = read_bin(work / "merged/statistics/ice_retention_class.bin")
    n_cls = sum(1 for f in range(nf) if (0 if c_mr[f] == 3 else c_mr[f]) != refs_m[f][1])
    dz = max(abs(z_mr[f] - refs_m[f][0]) for f in range(nf) if refs_m[f][1] < 2)
    gate("C1", bad == 0 and bad_code > 0 and "3.17e-17" in thr and ("total samples = %d" % n_m) in mlog
         and n_cls == 0 and dz <= 1e-4,
         "merged: retention vs reference %d class diffs, |dz| %.2g m; flux at 1 kg/m2/Gyr %d bad, %d differ from 3e-14"
         % (n_cls, dz, bad, bad_code))
    ok2 = all(provenance(work / c, "ice_depth_temperatures").startswith("radiation_coupled")
              and read_record(work / c).get("flux_threshold_kgm2s", "").startswith("3.16")
              and read_record(work / c).get("criterion") == "retention"
              and read_record(work / c).get("rho_surface_kgm3") is not None
              for c in ("chunk_000", "chunk_001"))
    gate("C2", ok2 and "NOTICE: this config sets thermal.ice_stability.flux_threshold_kgm2s" in l0
         and (work / "chunk_000/statistics/ice_depth_flux.bin").exists(),
         "chunks: coupled temperatures, both criteria recorded, notice for the flux key")
    rc, olog = run([a.merge, "--chunks", work / "chunk_000", work / "chunk_001", "--output",
                    work / "merged_3e-14", "--J0", "3e-14"], work / "merge_override.log")
    z_o = read_bin(work / "merged_3e-14/statistics/ice_depth_flux.bin")
    bad_o, cls_o, _ = criterion_check(z_o, nl, nf, sum_m, n_m, zg, J_CODE)
    gate("C3", bad_o == 0 and "OVERRIDES" in olog,
         "second flux layer from the same sums (--J0 3e-14): %d bad (%d/%d/%d), warned"
         % (bad_o, *cls_o))
    mixed = work / "mixed" / "chunk_001"
    if mixed.exists():
        shutil.rmtree(mixed)
    shutil.copytree(work / "chunk_001", mixed)
    prm = mixed / "statistics/_raw/ice_depth_params.txt"
    prm.write_text(prm.read_text().replace(read_record(work / "chunk_001")["flux_threshold_kgm2s"], "3e-14", 1))
    rc, xlog = run([a.merge, "--chunks", work / "chunk_000", mixed, "--output", work / "merged_mixed"],
                   work / "merge_mixed.log", expect_ok=False)
    gate("C4", rc != 0 and "different ice criteria" in xlog,
         "chunks with different thresholds refused (exit %d)" % rc)
    # Chunks written before the retention criterion (two-line record), alone and
    # mixed with new ones: the merge still writes their flux product, and the
    # retention product once the density profile is stated.
    legacy = work / "legacy"
    if legacy.exists():
        shutil.rmtree(legacy)
    for c in ("chunk_000", "chunk_001"):
        shutil.copytree(work / c, legacy / c)
        rr_ = read_record(work / c)
        (legacy / c / "statistics/_raw/ice_depth_params.txt").write_text(
            "flux_threshold_kgm2s %s\ngrain_hop_length_m %s\n"
            % (rr_["flux_threshold_kgm2s"], rr_["grain_hop_length_m"]))
    rc_l, llog = run([a.merge, "--chunks", legacy / "chunk_000", legacy / "chunk_001", "--output",
                      work / "merged_legacy"], work / "merge_legacy.log", expect_ok=False)
    ml = work / "merged_legacy/statistics"
    rc_m, mxlog = run([a.merge, "--chunks", legacy / "chunk_000", work / "chunk_001", "--output",
                       work / "merged_mixed_fmt", "--rho-s", "1100", "--rho-d", "1800", "--H-m", "0.06"],
                      work / "merge_mixed_fmt.log", expect_ok=False)
    mm = work / "merged_mixed_fmt/statistics"
    flux_m = (work / "merged/statistics/ice_depth_flux.bin").read_bytes()
    gate("C7", rc_l == 0 and not (ml / "ice_depth.bin").exists() and "NO ice_depth.bin" in llog
         and (ml / "ice_depth_flux.bin").read_bytes() == flux_m
         and rc_m == 0 and (mm / "ice_depth_flux.bin").read_bytes() == flux_m
         and (mm / "ice_depth.bin").read_bytes() == (work / "merged/statistics/ice_depth.bin").read_bytes(),
         "pre-retention chunks: flux product kept, warned; mixed formats + density: both products identical")
    c_dir = ret["coupled"][1]
    agree = sum(1 for x, y in zip(c_mr, c_dir) if x == y)
    dzc = [abs(x - y) for x, y, c1_, c2_ in zip(z_mr, ret["coupled"][0], c_mr, c_dir) if c1_ == c2_ == 1]
    gate("C5", agree >= 0.99 * nf,
         "coupled chunked vs direct: classes agree on %d of %d faces; at-depth |dz| max %.2g m"
         % (agree, nf, max(dzc) if dzc else 0.0))
    # The identity itself. A chunked window differs from a single-shot window
    # in its spin-up (each chunk-0 equilibrates over its own sub-window,
    # PTM-802; each coupled chunk re-seeds its re-thermal, PTM-862), so the
    # comparison starts BOTH from the same checkpoint: the off-radiation state
    # at mid-window (offchunk_000's final state). Direct: the second half in
    # one run. Chunked: the same 59 steps as 29 + 30, merged.
    o0 = config(work, "offchunk_000", a.threads, False, steps=half, days=days_half, offset_s=0.0)
    run([a.engine, "--config", o0], work / "offchunk_000.log")
    mid = work / "offchunk_000_final.ptmt"
    q1, q2 = 29, half - 29
    dd = config(work, "warm_direct", a.threads, False, steps=half, days=days_half,
                offset_s=half * cadence_s, initial=mid)
    run([a.engine, "--config", dd], work / "warm_direct.log")
    wa = config(work, "warm_000", a.threads, False, steps=q1, days=q1 * cadence_s / 86400.0,
                offset_s=half * cadence_s, initial=mid)
    run([a.engine, "--config", wa], work / "warm_000.log")
    wb = config(work, "warm_001", a.threads, False, steps=q2, days=q2 * cadence_s / 86400.0,
                offset_s=(half + q1) * cadence_s, initial=work / "warm_000_final.ptmt")
    run([a.engine, "--config", wb], work / "warm_001.log")
    run([a.merge, "--chunks", work / "warm_000", work / "warm_001", "--output",
         work / "warm_merged", "--allow-nonstandard-names"], work / "warm_merge.log")
    z_wd = read_bin(work / "warm_direct/statistics/ice_depth.bin")
    c_wd = read_bin(work / "warm_direct/statistics/ice_retention_class.bin")
    z_wm = read_bin(work / "warm_merged/statistics/ice_depth.bin")
    c_wm = read_bin(work / "warm_merged/statistics/ice_retention_class.bin")
    n_cd = sum(1 for x, y in zip(c_wm, c_wd) if x != y)
    dzo = max(abs(x - y) for x, y in zip(z_wm, z_wd))
    dto = max(abs(x - y) for x, y in zip(read_bin(work / "warm_merged/statistics/temperature_avg.bin"),
                                         read_bin(work / "warm_direct/statistics/temperature_avg.bin")))
    n_wd = read_raw(work / "warm_direct")[3]
    n_wm = read_raw(work / "warm_000")[3] + read_raw(work / "warm_001")[3]
    gate("C6", n_cd == 0 and dzo <= 1e-4 and n_wd == n_wm,
         "same start state, 1 run vs 2 chunks: %d class diffs, max |dz| %.2g m, samples %d/%d "
         "(T_avg max |dT| %.2g K)" % (n_cd, dzo, n_wd, n_wm, dto))

    if fails:
        print("\n%s FAIL: %s (kept in %s)" % (TAG, ", ".join(fails), work))
        return 1
    if not a.keep:
        shutil.rmtree(work, ignore_errors=True)
    print("\n%s PASS" % TAG)
    return 0


if __name__ == "__main__":
    sys.exit(main())
