#!/usr/bin/env python3
"""Checks of the temperatures at depth of ptm.observe (``--depth-m``).

    python3 src/python/tests/test_observe_depth.py

A synthetic run directory is written with the layout the engine uses
(mesh.ply, time_steps.txt, per-step frames, statistics/depth_map_z.json). Its
frames hold the periodic conduction solution

    T(z, t) = T0 + A exp(-z/d) cos(w t - z/d + phase of the face),

so every check is against a closed form, never against the module itself.
No SPICE kernels or engine runs are needed.

1. depth 0 is the surface temperature frame, value for value, at an output
   step and between two steps, and with or without depth frames in the run;
2. a requested depth equal to a stored level returns that level's values
   exactly;
3. at depths and times between the stored ones the analytic profile is
   recovered within the stated bound A h^2 exp(-z/d) / (4 d^2) for level
   spacing h, plus the time-interpolation bound A (w dt)^2 / 8, and the
   largest error falls by a factor of about 4 each time the spacing is halved;
4. several depths in one call equal separate single-depth calls, value for
   value, and the disk mean at depth 0 equals mean_temperature_k of
   band_emission;
5. errors: a depth below the surface in a run without depth frames (naming
   output.map_depths_per_step and output.map_depths_m), a depth below the
   deepest stored level (naming the range), a negative or repeated depth;
6. the command line (--depth-m, --csv, --faces-csv, --depths-npz,
   --list-depths) writes the values the API returns;
7. mutation: an off-by-one in the choice of the stored level below a depth
   makes check 3 fail, so the check can fail.
"""

import contextlib
import csv
import io
import json
import math
import struct
import sys
import tempfile
import textwrap
import inspect
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "src" / "python"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from ptm import observe  # noqa: E402
from ptm.results import Results  # noqa: E402
from ptm import cli  # noqa: E402
from test_observe import icosphere  # noqa: E402

FAILS = []
R_KM = 1560.8
T0, A = 100.0, 20.0           # K
DELTA = 0.02                  # skin depth, m
PERIOD = 3.0 * 3600.0         # s
OMEGA = 2.0 * math.pi / PERIOD
DT = 60.0                     # output step, s
N_STEPS = 40
ET0 = 789000000.0


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def analytic(z, t, phase):
    """T(z, t) for every face (phase per face)."""
    return T0 + A * np.exp(-z / DELTA) * np.cos(OMEGA * t - z / DELTA + phase)


def field_name(z):
    return f"temperature_z{int(z)}p{int(round((z - int(z)) * 1000)):03d}m"


def write_frame(path, values):
    values = np.asarray(values, dtype=np.float32)
    with open(path, "wb") as handle:
        handle.write(struct.pack("<i", values.size))
        handle.write(values.tobytes())


def write_run(root, levels, mesh_level=2):
    """A run directory with depth frames at ``levels`` (meters, may include
    0). Returns the per-face phase."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    mesh = icosphere(mesh_level, R_KM)
    v, f = mesh["vertices"], mesh["faces"]
    centers = v[f].mean(axis=1)
    lat = np.degrees(np.arcsin(centers[:, 2] / np.linalg.norm(centers, axis=1)))
    lon = np.degrees(np.arctan2(centers[:, 1], centers[:, 0]))
    normals = v / np.linalg.norm(v, axis=1)[:, None]
    with open(root / "mesh.ply", "wb") as handle:
        handle.write(textwrap.dedent(f"""\
            ply
            format binary_little_endian 1.0
            element vertex {len(v)}
            property double x
            property double y
            property double z
            property double nx
            property double ny
            property double nz
            element face {len(f)}
            property list uchar int vertex_indices
            property float latitude
            property float longitude
            end_header
            """).encode())
        for p, n in zip(v, normals):
            handle.write(struct.pack("<6d", *p, *n))
        for tri, la, lo in zip(f, lat, lon):
            handle.write(struct.pack("<B3i2f", 3, *[int(i) for i in tri], la, lo))
    (root / "metadata.txt").write_text(
        f"faces={len(f)}\ntime_steps={N_STEPS}\nsuccessful_timesteps={N_STEPS}\n"
        f"failed_timesteps=0\nis_custom_mesh=false\ncomponents=temperature\n")
    import datetime
    t0 = datetime.datetime(2025, 1, 1)
    with open(root / "time_steps.txt", "w") as handle:
        handle.write("# Time steps metadata\n# Column 1: Ephemeris time, seconds "
                     "past J2000 (TDB)\n# Column 2: UTC time\n")
        for k in range(N_STEPS):
            utc = (t0 + datetime.timedelta(seconds=k * DT)).strftime("%Y-%m-%dT%H:%M:%S.000")
            handle.write(f"{ET0 + k * DT:.6f}\t{utc}\n")
    phase = np.radians(lon)
    for name, z in [("temperature", 0.0)] + [(field_name(z), z) for z in levels]:
        d = root / name
        d.mkdir(exist_ok=True)
        for k in range(N_STEPS):
            write_frame(d / f"{name}_{k:06d}.bin", analytic(z, k * DT, phase))
    if levels:
        (root / "statistics").mkdir(exist_ok=True)
        entries = [{"requested_m": z, "resolved_m": z, "residual_m": 0.0,
                    "on_node": False, "field": field_name(z)} for z in levels]
        (root / "statistics" / "depth_map_z.json").write_text(json.dumps(
            {"schema": "ptm3d-depth-map/1", "skin_depth_m": DELTA,
             "native_z_m": [0.0, 0.001, 1.0], "depths": entries}))
    return phase


def geometry_between(k0, w1):
    """A geometry at the time fraction ``w1`` between steps k0 and k0 + 1, as
    observer_geometry builds it with time_match='interpolate'."""
    g = observe.explicit_geometry(10.0, 30.0, 6.3e8)
    g.steps, g.weights, g.step = (k0, k0 + 1), (1.0 - w1, w1), k0
    return g


def stored_frame(run, name, g):
    return observe.temperature_at(run, g, name)


def test_surface_and_levels(run, levels):
    print("1-2. depth 0 and stored levels")
    for g, what in ((observe.explicit_geometry(0, 0, 6.3e8, run=run, step=7), "at a step"),
                    (geometry_between(7, 0.37), "between steps")):
        dt = observe.temperature_at_depths(run, g, [0.0])
        check(f"depth 0 equals the surface frame {what}",
              np.array_equal(dt.temperature_k[0], observe.temperature_at(run, g)))
        ok = True
        for z in levels:
            dz = observe.temperature_at_depths(run, g, [z])
            ok &= (np.array_equal(dz.temperature_k[0], stored_frame(run, field_name(z), g))
                   and not dz.interpolated[0])
        check(f"each of {len(levels)} stored levels returned exactly {what}", ok)


def max_error(run, phase, levels, times):
    """Largest error against the analytic profile at the midpoints between
    stored levels, and the bound there."""
    zs = np.concatenate([[0.0], levels])
    mids = 0.5 * (zs[:-1] + zs[1:])
    err, bound = 0.0, 0.0
    for k0, w1 in times:
        g = geometry_between(k0, w1)
        dt = observe.temperature_at_depths(run, g, mids)
        t = (k0 + w1) * DT
        for i, z in enumerate(mids):
            e = np.max(np.abs(dt.temperature_k[i] - analytic(z, t, phase)))
            h = zs[i + 1] - zs[i]
            b = (A * h ** 2 * math.exp(-zs[i] / DELTA) / (4 * DELTA ** 2)
                 + (A * (OMEGA * DT) ** 2 / 8 if 0 < w1 < 1 else 0.0) + 1e-4)
            err = max(err, e)
            if e > b:
                bound = max(bound, e / b)
    return err, bound


def test_analytic(tmp):
    print("3. analytic profile between stored depths and times")
    errs = []
    for n in (6, 12, 24, 48):
        levels = list(np.linspace(0, 3 * DELTA, n + 1)[1:])
        root = Path(tmp) / f"grid{n}"
        phase = write_run(root, levels, mesh_level=1)
        run = Results(root, mesh_units="km")
        # At output steps (depth error alone) and between them.
        e_step, over_step = max_error(run, phase, levels, [(5, 0.0), (20, 0.0)])
        e_mid, over_mid = max_error(run, phase, levels, [(5, 0.5), (21, 0.3)])
        errs.append(e_step)
        check(f"h = {3 * DELTA / n * 1000:.2f} mm: within the bound at steps "
              f"({e_step:.2e} K) and between steps ({e_mid:.2e} K)",
              over_step == 0.0 and over_mid == 0.0,
              f"worst error/bound {max(over_step, over_mid):.2f}" if max(over_step, over_mid) else "")
    ratios = [errs[i] / errs[i + 1] for i in range(len(errs) - 1)]
    check("error falls about 4x per halving of the spacing (second order)",
          all(3.3 < r < 4.5 for r in ratios), ", ".join(f"{r:.2f}" for r in ratios))
    return errs


def test_multiple_and_mean(run):
    print("4. several depths in one call; disk mean")
    g = geometry_between(11, 0.62)
    z = [0.0, 0.0031, 0.01, 0.017, 0.04]
    together = observe.temperature_at_depths(run, g, z)
    alone = [observe.temperature_at_depths(run, g, [zi]).temperature_k[0] for zi in z]
    check("one call with five depths equals five single-depth calls",
          all(np.array_equal(together.temperature_k[i], alone[i]) for i in range(len(z))))
    flags = [bool(f) for f in together.interpolated]
    check("interpolated flags", flags == [False, True, False, True, False], str(flags))
    band = observe.Band.top_hat(8.0, 14.0)
    em = observe.band_emission(run, observe.temperature_at(run, g), g, band=band,
                               emissivity=0.9, reflected=False)
    means = observe.disk_mean_temperature(run, together.temperature_k, g)
    check("disk mean at depth 0 equals band_emission's mean_temperature_k",
          means[0] == em.mean_temperature_k, f"{float(means[0])!r} vs {em.mean_temperature_k!r}")
    check("disk mean of an isothermal field is that temperature",
          abs(observe.disk_mean_temperature(run, np.full(run.n_faces, 123.0), g) - 123.0) < 1e-12)


def test_errors(tmp, run):
    print("5. errors")
    root = Path(tmp) / "nodepth"
    write_run(root, [], mesh_level=1)
    bare = Results(root, mesh_units="km")
    g = geometry_between(3, 0.5)
    dt = observe.temperature_at_depths(bare, g, [0.0])
    check("depth 0 works without depth frames",
          np.array_equal(dt.temperature_k[0], observe.temperature_at(bare, g)))
    try:
        observe.temperature_at_depths(bare, g, [0.0, 0.01])
        check("a depth below the surface without frames is refused", False)
    except FileNotFoundError as exc:
        msg = str(exc)
        check("a depth below the surface without frames is refused, naming the keys",
              "map_depths_per_step: true" in msg and "map_depths_m: [0.01]" in msg)
    g = geometry_between(3, 0.5)
    try:
        observe.temperature_at_depths(run, g, [0.0, 0.05])
        check("a depth below the deepest level is refused", False)
    except ValueError as exc:
        check("a depth below the deepest level is refused, naming the range",
              "from 0 to 0.04 m" in str(exc), str(exc)[:80])
    for bad, what in (([-0.001], "negative"), ([0.01, 0.01], "repeated"),
                      ([float("nan")], "NaN"), ([], "empty")):
        try:
            observe.temperature_at_depths(run, g, bad)
            check(f"{what} depth refused", False)
        except ValueError:
            check(f"{what} depth refused", True)


def run_cli(argv):
    out, err = io.StringIO(), io.StringIO()
    code = 0
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            code = cli.main(argv)
        except SystemExit as exc:
            code = exc.code
    return code, out.getvalue(), err.getvalue()


def test_cli(tmp, root, run):
    print("6. command line")
    tmp = Path(tmp)
    z = [0.0, 0.004, 0.02]
    code, out, err = run_cli(["observe", str(root), "--sub-observer", "10", "30", "6.3e8",
                              "--step", "4", "--step", "9", "--depth-m"] + [str(v) for v in z]
                             + ["--csv", str(tmp / "t.csv"), "--faces-csv", str(tmp / "f.csv"),
                                "--depths-npz", str(tmp / "d.npz")])
    check("ptm observe --depth-m runs without a band", code == 0, err.strip()[:200])
    if code != 0:
        return
    check("the header names the interpolated depth",
          "0.004 m is not a stored level: interpolated linearly in depth between 0 and 0.005 m" in out)
    geoms = [observe.explicit_geometry(10, 30, 6.3e8, run=run, step=s) for s in (4, 9)]
    api = [observe.temperature_at_depths(run, g, z) for g in geoms]
    means = [observe.disk_mean_temperature(run, a.temperature_k, g) for a, g in zip(api, geoms)]
    with open(tmp / "t.csv") as handle:
        rows = list(csv.DictReader(handle))
    check("--csv holds the disk means", all(
        abs(float(r[f"mean_temperature_k_z{v:g}m"]) - m[i]) < 1e-9
        for r, m in zip(rows, means) for i, v in enumerate(z)))
    with open(tmp / "f.csv") as handle:
        frows = list(csv.DictReader(handle))
    worst = max(abs(float(r[f"temperature_k_z{v:g}m"])
                    - api[int(r["time_index"])].temperature_k[i, int(r["face"])])
                for r in frows for i, v in enumerate(z))
    mu = [observe.face_weights(run, g)[0] for g in geoms]
    check("--faces-csv holds every visible face at every depth",
          len(frows) == sum(int(np.sum(m > 0)) for m in mu) and worst <= 5e-5,
          f"{len(frows)} rows, largest difference {worst:.1e} K")
    d = np.load(tmp / "d.npz")
    check("--depths-npz holds (time, depth, face) equal to the API",
          d["temperature_k"].shape == (2, 3, run.n_faces)
          and np.array_equal(d["temperature_k"], np.stack([a.temperature_k for a in api]))
          and np.array_equal(d["mu"], np.stack(mu)), str(d["temperature_k"].shape))
    code, out, _ = run_cli(["observe", str(root), "--list-depths"])
    check("--list-depths lists the stored levels",
          code == 0 and "0.005" in out and "0.04 " in out and "the surface" in out)
    code, out, err = run_cli(["observe", str(root), "--sub-observer", "10", "30", "6.3e8",
                              "--step", "4"])
    check("without --depth-m a band or frequency is still required",
          code == 2 and "one of the arguments --freq-ghz" in err)


def test_mutation(tmp):
    print("7. mutation: off-by-one in the stored level below a depth")
    src = inspect.getsource(observe.temperature_at_depths)
    original = 'int(np.searchsorted(levels, zi, side="right")) - 1'
    if original not in src:
        check("mutation site found", False)
        return
    mutated = src.replace(original, 'min(int(np.searchsorted(levels, zi, side="right")), len(levels) - 2)')
    namespace = dict(vars(observe))
    exec(compile(mutated, "<mutant>", "exec"), namespace)
    saved = observe.temperature_at_depths
    observe.temperature_at_depths = namespace["temperature_at_depths"]
    n_before = len(FAILS)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            levels = list(np.linspace(0, 3 * DELTA, 13)[1:])
            root = Path(tmp) / "mutant"
            phase = write_run(root, levels, mesh_level=1)
            mrun = Results(root, mesh_units="km")
            _, over = max_error(mrun, phase, levels, [(5, 0.0)])
            g = observe.explicit_geometry(0, 0, 6.3e8, run=mrun, step=5)
            exact = all(np.array_equal(
                observe.temperature_at_depths(mrun, g, [z]).temperature_k[0],
                observe.temperature_at(mrun, g, field_name(z))) for z in levels[:-1])
    finally:
        observe.temperature_at_depths = saved
    del FAILS[n_before:]
    check("the analytic check fails on the mutant", over > 1.0,
          f"worst error/bound {over:.1f}")
    check("the stored-level check fails on the mutant", not exact)


def main():
    with tempfile.TemporaryDirectory() as tmp:
        levels = [0.005, 0.01, 0.02, 0.03, 0.04]
        root = Path(tmp) / "run"
        write_run(root, levels)
        run = Results(root, mesh_units="km")
        test_surface_and_levels(run, levels)
        test_analytic(tmp)
        test_multiple_and_mean(run)
        test_errors(tmp, run)
        test_cli(tmp, root, run)
        test_mutation(tmp)
    if FAILS:
        print(f"FAIL: {len(FAILS)} check(s): {', '.join(FAILS)}")
        return 1
    print("PASS: all checks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
