#!/usr/bin/env python3
"""End-to-end analytic gates for the flux-cube product on the synthetic 80N crater.

The unit test (test_flux_cube_pass) checks the pass arithmetic on synthetic
callbacks. This test checks what the engine actually delivers: temperatures,
illumination and the radiosity solution flowing through to cube fluxes, on the
public synthetic bowl crater `data/meshes/lunar_crater_80N.ply` (1 km across,
200 m deep, at 80N 0E, on a ~2 km square plain).

It builds two geometry caches (cube heights 0.35 m and 2.0 m) at five named
locations -- the bowl-floor center and four points on the flat plain 700 m from
the center -- runs one lunation (118 output steps) with two-band radiation on
the CPU, and checks these limits:

  plain (slope 0, sky fraction 1 / 0 / 0.5 on the zenith / nadir / side faces)
    G1  direct on every face = S(t) * max(0, n.s) whenever the Sun is clearly up
        and the cube is not penumbral (one shared visibility per cube);
    G2  nadir terrain_ir = eps*sigma*T_face^4 (an infinite plane sees no terrain
        from below the horizon, so its indirect flux is ~0);
    G3  side / nadir terrain_ir = 0.5 (side faces see half the plane);
    G4  side / nadir terrain_vis = 0.5 at sunlit steps;
    G5  zenith terrain_ir and terrain_vis = 0 exactly (the zenith face sees sky);
  bowl floor (slope ~1 deg, shadowed for much of the lunation)
    G6  zenith terrain_ir = the facet beneath's incident indirect infrared from
        the engine's own view-factor solution (indirect_ir / eps), median ratio.

G1-G5 are closed forms independent of PTM3D. G6 compares two parts of the
same engine (the cube's ray-traced view-factor row against the run's H-matrix
view factors), so it is a consistency check, not a validation.

Needs the Moon SPICE kernels (python3 scripts/fetch_kernels.py --bodies moon),
build/bin/ptm_thermal_sim and build/bin/ptm_flux_cube_geom. Standard library
only. About 30 s with 4 threads.

  python3 scripts/test_flux_cube_analytic.py [--threads N] [--keep DIR]
"""
import argparse
import importlib.util
import math
import os
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SIGMA = 5.670374419e-8
S_1AU = 1361.0                 # constants.cpp SOLAR_CONSTANT_1AU
KM_PER_AU = 149597870.7
EPS = 0.95                     # lunar default emissivity
LOCATIONS = """name,lat_deg,lon_deg
floor,80.0,0.0
plain_n,80.023085,0.0
plain_s,79.976915,0.0
plain_e,80.0,0.132940
plain_w,80.0,-0.132940
"""
PLAIN = ("plain_n", "plain_s", "plain_e", "plain_w")
HEIGHTS = (("h0p35", 0.35), ("h2p00", 2.0))

CONFIG = """simulation:
  start_time: '2026-04-05T00:00:00'
  duration_days: 29.53
  timesteps: 118
  body: Moon
mesh:
  import_file: {mesh}
  scale_factor: 1.0
  use_topographic: false
thermal:
  run_thermal_mesh: true
  threads: {threads}
  equilibration_orbits: 1.0
  equilibration_steps_per_day: -1.0
  use_gpu: false
  compute_statistics: true
  enable_shadow_flux: false
  crater_dd_ratio: 8.0
  geothermal_heat_flow: -1.0
thermophysical:
  use_default_albedo: true
  incidence_angle_albedo: true
radiation:
  enabled: true
  use_gpu: false
  use_hmatrix: true
  use_aca: true
  aca_tolerance: 0.001
  hmatrix_threads: {threads}
  view_factor_file: {work}/crater80N.vf
  max_coupling_iterations: 5
  coupling_temp_threshold: 2.0
  coupling_flux_threshold: 0.1
  coupling_face_fraction: 0.01
  band_separated: true
illumination:
  include_primary_reflection: true
  use_batch_processing: true
  use_albedo_map: false
output:
  directory: {work}/out
  format: binary
advanced:
  verbose: false
flux_cubes:
  enabled: true
  cache_files:
  - {work}/crater_h0p35.fcg
  - {work}/crater_h2p00.fcg
  locations:
  - {{name: floor, lat_deg: 80.0, lon_deg: 0.0}}
  - {{name: plain_n, lat_deg: 80.023085, lon_deg: 0.0}}
  - {{name: plain_s, lat_deg: 79.976915, lon_deg: 0.0}}
  - {{name: plain_e, lat_deg: 80.0, lon_deg: 0.132940}}
  - {{name: plain_w, lat_deg: 80.0, lon_deg: -0.132940}}
"""


def load_packager():
    """The .fcb reader of scripts/ptm_flux_cubes_to_hdf5.py, so the layout is
    spelled in one place. Standard library only."""
    spec = importlib.util.spec_from_file_location(
        "_fc_packager", REPO / "scripts" / "ptm_flux_cubes_to_hdf5.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_fcb(pk, path):
    d = pk.read_fcb(str(path))
    return dict(nf=d["n_faces"], ncomp=d["n_comp"], ns=d["n_steps"], flux=d["flux"],
                tf=d["t_face"], pen=d["penumbral"], fsky=d["f_sky"], cubes=d["cubes"],
                names=d["names"], U=d["frame_up"], N=d["frame_north"], E=d["frame_east"],
                h=d["height_m"])


def F(d, cu, fa, co, st):
    return d["flux"][((cu * d["nf"] + fa) * d["ncomp"] + co) * d["ns"] + st]


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def median(xs):
    s = sorted(xs)
    n = len(s)
    return 0.5 * (s[(n - 1) // 2] + s[n // 2]) if n else float("nan")


def read_frame(path):
    b = Path(path).read_bytes()
    n = struct.unpack_from("<I", b, 0)[0]
    return struct.unpack_from("<%df" % n, b, 4)


def run(cmd, log, cwd):
    with open(log, "w") as fh:
        r = subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        sys.stdout.write(Path(log).read_text()[-4000:])
        raise SystemExit("[test_flux_cube_analytic] FAIL: %s exited %d (log %s; the work "
                         "directory is kept)" % (Path(cmd[0]).name, r.returncode, log))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--engine", default=str(REPO / "build/bin/ptm_thermal_sim"))
    ap.add_argument("--geom", default=str(REPO / "build/bin/ptm_flux_cube_geom"))
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--keep", default=None,
                    help="work directory to keep (default: a temp dir, removed on success)")
    a = ap.parse_args()

    for exe in (a.engine, a.geom):
        if not os.access(exe, os.X_OK):
            raise SystemExit("[test_flux_cube_analytic] FAIL: %s is not built" % exe)
    if not (REPO / "data/kernels/spk/de440s.bsp").exists():
        raise SystemExit("[test_flux_cube_analytic] FAIL: Moon SPICE kernels missing; run "
                         "python3 scripts/fetch_kernels.py --bodies moon")
    work = Path(a.keep) if a.keep else Path(tempfile.mkdtemp(prefix="ptm_fc_analytic_"))
    work.mkdir(parents=True, exist_ok=True)
    mesh = REPO / "data/meshes/lunar_crater_80N.ply"
    (work / "locs.csv").write_text(LOCATIONS)
    env_threads = str(a.threads)
    os.environ["OMP_NUM_THREADS"] = env_threads
    for tag, h in HEIGHTS:
        run([a.geom, "--core", str(mesh), "--center-lat", "80", "--center-lon", "0",
             "--radius", "-1", "--locations", str(work / "locs.csv"), "--height-m", str(h),
             "--edge-m", "0.5", "--threads", env_threads,
             "--out", str(work / ("crater_%s.fcg" % tag))], work / ("geom_%s.log" % tag), REPO)
    cfg = work / "run.yaml"
    cfg.write_text(CONFIG.format(mesh=mesh, threads=a.threads, work=work))
    run([a.engine, "--config", str(cfg)], work / "engine.log", REPO)

    out = work / "out"
    rows = (out / "body_positions.csv").read_text().splitlines()[1:]
    sun = [tuple(float(x) for x in r.split(",")[1:4]) for r in rows]

    fails = []
    pk = load_packager()

    def gate(name, ok, detail):
        print("  %-4s %-62s %s" % (name, detail, "PASS" if ok else "FAIL"))
        if not ok:
            fails.append(name)

    for tag, h in HEIGHTS:
        d = read_fcb(pk, out / "flux_cubes" / ("flux_cubes_crater_%s.fcb" % tag))
        ns = d["ns"]
        if len(sun) != ns:
            raise SystemExit("[test_flux_cube_analytic] FAIL: body_positions.csv has %d rows "
                             "for %d output steps" % (len(sun), ns))
        normals = [d["U"], tuple(-x for x in d["U"]), d["N"], tuple(-x for x in d["N"]),
                   d["E"], tuple(-x for x in d["E"])]
        idx = {n: k for k, n in enumerate(d["names"])}
        print("\n--- height %.2f m: %d output steps ---" % (d["h"], ns))

        # G1: direct = S max(0, n.s) on all six faces when clearly lit.
        worst, n_lit = 0.0, 0
        for nm in PLAIN:
            cu = idx[nm]
            for st in range(ns):
                s = sun[st]
                dist = math.sqrt(dot(s, s))
                sh = tuple(x / dist for x in s)
                if dot(sh, d["U"]) < math.sin(math.radians(0.5)) or d["pen"][cu * ns + st]:
                    continue
                S = S_1AU / (dist / KM_PER_AU) ** 2
                n_lit += 1
                for fa in range(6):
                    want = S * max(0.0, dot(normals[fa], sh))
                    worst = max(worst, abs(F(d, cu, fa, 0, st) - want) / S)
        gate("G1", n_lit >= 40 and worst <= 5e-5,
             "direct = S max(0,n.s): %d lit cube-steps, worst |d|/S %.1e" % (n_lit, worst))

        # G2-G5 on the plain.
        r_nad, r_ir, r_vis, zen_max = [], [], [], 0.0
        for nm in PLAIN:
            cu = idx[nm]
            for st in range(ns):
                T = d["tf"][cu * ns + st]
                nad_ir = F(d, cu, 1, 2, st)
                r_nad.append(nad_ir / (EPS * SIGMA * T ** 4))
                for fa in range(2, 6):
                    r_ir.append(F(d, cu, fa, 2, st) / nad_ir)
                nad_vis = F(d, cu, 1, 1, st)
                if nad_vis > 1.0:
                    for fa in range(2, 6):
                        r_vis.append(F(d, cu, fa, 1, st) / nad_vis)
                zen_max = max(zen_max, abs(F(d, cu, 0, 2, st)), abs(F(d, cu, 0, 1, st)))
        m = median(r_nad)
        gate("G2", abs(m - 1.0) <= 2e-3,
             "nadir terrain_ir / eps sigma T^4: median %.5f (n=%d)" % (m, len(r_nad)))
        m = median(r_ir)
        gate("G3", abs(m - 0.5) <= 0.01,
             "side / nadir terrain_ir: median %.4f (n=%d)" % (m, len(r_ir)))
        m = median(r_vis)
        gate("G4", len(r_vis) >= 100 and abs(m - 0.5) <= 0.01,
             "side / nadir terrain_vis (lit): median %.4f (n=%d)" % (m, len(r_vis)))
        fsky_zen = min(d["fsky"][idx[nm] * 6 + 0] for nm in PLAIN)
        gate("G5", fsky_zen == 1.0 and zen_max == 0.0,
             "zenith terrain = 0: max %.3g W/m^2, zenith f_sky %.4f" % (zen_max, fsky_zen))

        # G6: bowl floor zenith terrain_ir vs the engine's own indirect IR.
        cu = idx["floor"]
        face = d["cubes"][cu][7]
        r = []
        for st in range(ns):
            q = read_frame(out / "indirect_ir" / ("indirect_ir_%06d.bin" % st))[face] / EPS
            if q > 0.05:
                r.append(F(d, cu, 0, 2, st) / q)
        m = median(r)
        gate("G6", len(r) >= 50 and abs(m - 1.0) <= 0.05,
             "floor zenith terrain_ir / (Q_ind_ir/eps): median %.4f (n=%d)" % (m, len(r)))

    if fails:
        print("\n[test_flux_cube_analytic] FAIL: %s (inputs, logs and products kept in %s)"
              % (", ".join(fails), work))
        return 1
    if not a.keep:
        shutil.rmtree(work, ignore_errors=True)
    print("\n[test_flux_cube_analytic] PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
