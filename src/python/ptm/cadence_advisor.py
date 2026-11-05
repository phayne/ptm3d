"""Is the output step short enough to see a run's sunlit flashes? The Python half of the advisor.

The engine prints the full advisory during every run (``src/cpp/thermal/cadence_advisor.h``; silence it
with ``illumination.cadence_advisor: false``): the Sun's angular rate omega and disk theta from the
run's own epochs, the recommended output step

    dt_out <= 0.4 theta / omega          (about 25 min on the Moon),

the 5% peak-flux limit 0.3 theta/omega, and, after the run, an estimated upper bound on how much the
frame-sampled cold-trap area is overstated. Sunlit flashes on rough terrain include grazing partial-disk
events with no lower limit on duration, so no step catches all of them, and the size of the mesh's
occluders does not predict how short they are (Fermi, n4990 cadence ladder, 2026-11-04).

This module provides, without the engine:

* :func:`theta_over_omega` -- the solar-disk crossing time for a body;
* :func:`postrun` -- the post-run estimate from a FINISHED run directory, with the same numbers and the
  same text as the engine's ``[cadence-advisor] post-run`` line:
  R_ct(thr) = (area of ALL faces lit in exactly one output frame, or in one or two frames when the step
  is shorter than theta/omega) / (area with frame-sampled T_max < thr);
* ``ptm cadence <run_dir>`` on the shell.

The step and theta/omega are read from the run's own ``provenance.txt`` (written by the engine's
advisor); without them the one-frame count is used, and the result says so.
"""

from __future__ import annotations

import math
from pathlib import Path

__all__ = ["LIT_THRESHOLD_WM2", "REC_C", "PEAK5_W", "SUN_RADIUS_KM", "theta_over_omega",
           "recommended_step", "missed_flash_delta_t", "postrun", "format_postrun"]

LIT_THRESHOLD_WM2 = 1.0   # cadence_advisor.h kLitThresholdWm2
REC_C = 0.4               # kRecC
PEAK5_W = 0.3             # kPeak5W
SUN_RADIUS_KM = 695700.0
AU_KM = 1.495978707e8
SIGMA = 5.670374419e-8


def theta_over_omega(body=None, distance_au=None, solar_day_days=None, declination_deg=0.0):
    """(theta [deg], omega [deg/h], theta/omega [h]) for a body, at its mean distance by default.

    omega = (360 / P_syn) cos(declination); theta = 2 asin(R_sun / d). For an implemented body the
    distance and solar day come from the engine's own table (``ptm.body_table``).
    """
    if body is not None:
        from .body_table import BODY_CAPABILITIES
        props = BODY_CAPABILITIES[body]
        distance_au = props["rAU"] if distance_au is None else distance_au
        solar_day_days = props["solar_day_days"] if solar_day_days is None else solar_day_days
    theta = math.degrees(2.0 * math.asin(SUN_RADIUS_KM / (distance_au * AU_KM)))
    omega = 360.0 / (solar_day_days * 24.0) * math.cos(math.radians(declination_deg))
    return theta, omega, theta / omega


def recommended_step(theta_deg, omega_deg_per_h):
    """(recommended output step 0.4 theta/omega, 5% peak-flux limit 0.3 theta/omega), in hours."""
    return REC_C * theta_deg / omega_deg_per_h, PEAK5_W * theta_deg / omega_deg_per_h


def missed_flash_delta_t(F_abs, gamma, emissivity, T, dt_s):
    """The most a missed full-sun flash of duration dt can raise the surface [K]:
    min[(2 F/Gamma) sqrt(dt/pi), F/(4 eps sigma T^3)] (cadence_advisor.h missedFlashDeltaT)."""
    if not (F_abs > 0 and gamma > 0 and dt_s > 0):
        return 0.0
    return min(2.0 * F_abs / gamma * math.sqrt(dt_s / math.pi), F_abs / (4.0 * emissivity * SIGMA * T ** 3))


def _read_stat(path):
    import numpy as np
    raw = Path(path).read_bytes()
    n = int(np.frombuffer(raw[:4], "<u4")[0])
    v = np.frombuffer(raw[4:4 + 4 * n], "<f4")
    if v.size != n:
        raise ValueError(f"{path}: header says {n} values, file holds {v.size}")
    return v


def _provenance(run):
    out = {}
    p = run / "provenance.txt"
    if p.exists():
        for line in p.read_text(errors="replace").splitlines():
            k, sep, v = line.partition("=")
            if sep:
                out[k.strip()] = v.strip()
    return out


def postrun(run_dir):
    """The post-run estimate for a finished run. Returns a dict (area fractions of the mesh)."""
    import numpy as np
    from .results import read_ply_mesh_minimal
    run = Path(run_dir)
    frames = sorted((run / "direct_solar").glob("direct_solar_*.bin"))
    if not frames:
        raise FileNotFoundError(f"{run}: no direct_solar/ frames (the estimate needs per-step direct solar)")
    mesh = read_ply_mesh_minimal(run / "mesh.ply")
    V, F = mesh["vertices"], mesh["faces"]
    a, b, c = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    area = np.where(np.isfinite(area), area, 0.0)
    n_lit = np.zeros(len(F), dtype=np.int64)
    for f in frames:
        d = _read_stat(f)
        if d.size != len(F):
            raise ValueError(f"{f}: {d.size} faces, mesh has {len(F)}")
        n_lit += d > LIT_THRESHOLD_WM2
    prov = _provenance(run)
    max_lit, rule = 1, "one frame (no cadence numbers in provenance.txt)"
    try:
        dt = float(prov["cadence_dt_out_s"])
        to = float(prov["cadence_sun_disk_deg"]) / float(prov["cadence_sun_rate_deg_per_h"]) * 3600.0
        max_lit = 2 if 0 < dt < to else 1
        rule = "from provenance.txt"
    except (KeyError, ValueError, ZeroDivisionError):
        pass
    tot = area.sum()
    ind = (n_lit >= 1) & (n_lit <= max_lit)
    out = {"n_frames": len(frames), "max_frames_lit": max_lit, "rule": rule,
           "indicator_area_frac": float(area[ind].sum() / tot)}
    tmax_path = run / "statistics" / "temperature_max.bin"
    if tmax_path.exists():
        T = _read_stat(tmax_path)
        for thr, key in ((110.0, "110K"), (54.0, "54K")):
            ct = area[T < thr].sum()
            out[f"ct_area_frac_{key}"] = float(ct / tot)
            out[f"ct_overstatement_bound_{key}"] = float(area[ind].sum() / ct) if ct > 0 else -1.0
    return out


def _pct(x):
    return "%.3g%%" % (100.0 * x)


def format_postrun(r):
    """The engine's post-run line, character for character (cadence_advisor.h formatPostRun)."""
    phrase = "exactly one output frame" if r["max_frames_lit"] == 1 else "one or two output frames"
    if "ct_area_frac_110K" not in r:
        return (f"[cadence-advisor] post-run: faces lit in {phrase}: {_pct(r['indicator_area_frac'])} of the mesh "
                f"area (no per-face T_max for this run, so no cold-trap estimate).")
    est = lambda v: "n/a (no area below the threshold)" if v < 0 else _pct(v)  # noqa: E731
    return ("[cadence-advisor] post-run: estimated upper bound on the fractional overstatement of the cold-trap "
            f"area = (area of all faces lit in {phrase}, {_pct(r['indicator_area_frac'])} of the mesh) / (area with "
            f"frame-sampled T_max below the threshold): {est(r['ct_overstatement_bound_110K'])} at 110 K (cold traps "
            f"{_pct(r['ct_area_frac_110K'])} of the mesh), {est(r['ct_overstatement_bound_54K'])} at 54 K "
            f"({_pct(r['ct_area_frac_54K'])}). On the n4990 cadence ladder (one phase) this estimate was "
            "4.74 / 1.24 / 1.26% at 5.9 h / 2 h / 1 h steps (one frame, one frame, one or two frames), against "
            "true overstatements of 2.5 / 0.71-0.87 / 0.38-0.54%.")
