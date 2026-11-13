"""DAMIT rotation state expressed as PTM3D ``generic_body.spin`` keys.

DAMIT documents the orientation of a model at time t as

    r_ecl = R_z(lambda) R_y(90 deg - beta) R_z(phi0 + 2 pi (t - t0) / P) r_ast

where r_ast is a vector in the frame of the shape file, r_ecl the same vector
in the J2000 ecliptic frame, (lambda, beta) the ecliptic pole, P the sidereal
period, t0 the reference epoch (a Julian date) and phi0 the rotation angle at
t0. R_i(theta) rotates a vector anticlockwise by theta about axis i. The light
curves DAMIT fits are given at light-time-corrected epochs, so t is the time at
the asteroid. DAMIT does not name the time scale of t0. It is used here as a
TDB Julian date, the scale of ``generic_body.spin.epoch``. Models with a YORP
term carry it in the DAMIT database, not in spin.txt, and it is not included.

PTM3D builds the body-fixed frame of a generic body (src/cpp/orbital/
body_fixed_frame.cpp) from the pole P_hat, a reference axis e1 and the angle W:

    e1 = the J2000 vernal equinox (ecliptic +x) projected onto the plane
         normal to the pole, or ecliptic +y when |x . P_hat| > 0.9
    e2 = P_hat x e1
    x_body = cos W e1 + sin W e2,  z_body = P_hat
    W(t) = prime_meridian_deg + 360 (t - spin.epoch) / P

The DAMIT axes at phase phi are the same construction with a different
reference axis. R_z(lambda) R_y(90 deg - beta) carries +z to the pole and +x to

    u = (sin beta cos lambda, sin beta sin lambda, -cos beta),

so x_ast = cos phi u + sin phi (P_hat x u). Both frames are right-handed with
the same third axis, so they coincide when W = phi + delta, where delta is the
angle about the pole from e1 to u. The shape file then needs no rotation, and

    spin.epoch         = t0
    prime_meridian_deg = phi0 + delta
"""
from __future__ import annotations

import math

#: Obliquity the engine uses to rotate pole_ra_deg/pole_dec_deg back to the
#: ecliptic (orbital::SpinState::OBLIQUITY_J2000_RAD, the ECLIPJ2000 value of
#: SPICE). Converting the DAMIT pole with the same value makes the round trip
#: exact.
ENGINE_OBLIQUITY_DEG = 23.4392911


def _unit_pole(lambda_deg: float, beta_deg: float) -> tuple[float, float, float]:
    lam, bet = math.radians(lambda_deg), math.radians(beta_deg)
    return (math.cos(bet) * math.cos(lam), math.cos(bet) * math.sin(lam), math.sin(bet))


def phase_offset_deg(lambda_deg: float, beta_deg: float) -> float:
    """Angle delta (degrees) about the pole from the PTM3D reference axis e1
    to the DAMIT zero-phase axis u, for the ecliptic pole (lambda, beta)."""
    p = _unit_pole(lambda_deg, beta_deg)
    lam, bet = math.radians(lambda_deg), math.radians(beta_deg)
    u = (math.sin(bet) * math.cos(lam), math.sin(bet) * math.sin(lam), -math.cos(bet))
    seed = (1.0, 0.0, 0.0)
    if abs(p[0]) > 0.9:            # the engine's choice for a pole near the equinox
        seed = (0.0, 1.0, 0.0)
    proj = sum(s * q for s, q in zip(seed, p))
    e1 = [s - proj * q for s, q in zip(seed, p)]
    n = math.sqrt(sum(c * c for c in e1))
    e1 = [c / n for c in e1]
    e2 = (p[1] * e1[2] - p[2] * e1[1],
          p[2] * e1[0] - p[0] * e1[2],
          p[0] * e1[1] - p[1] * e1[0])
    return math.degrees(math.atan2(sum(a * b for a, b in zip(u, e2)),
                                   sum(a * b for a, b in zip(u, e1))))


def generic_body_spin(lambda_deg: float, beta_deg: float, period_hours: float,
                      jd0: float | None = None, phi0_deg: float | None = None) -> dict:
    """The ``generic_body.spin`` keys for a DAMIT spin state.

    Returns ``period_hours``, ``pole_ra_deg`` and ``pole_dec_deg``, and, when
    both ``jd0`` and ``phi0_deg`` are given, ``prime_meridian_deg`` and
    ``epoch`` (a string ``'JD <t0>'``). Without a phase the last two are
    omitted. The period is written to the precision DAMIT gives, because the
    phase at a date far from t0 depends on it.
    """
    from horizons.horizons_client import ecliptic_pole_to_equatorial

    ra, dec = ecliptic_pole_to_equatorial(lambda_deg, beta_deg, ENGINE_OBLIQUITY_DEG)
    spin = {
        "period_hours": float(period_hours),
        "pole_ra_deg": round(ra, 9),
        "pole_dec_deg": round(dec, 9),
    }
    if jd0 is not None and phi0_deg is not None:
        w0 = (float(phi0_deg) + phase_offset_deg(lambda_deg, beta_deg)) % 360.0
        spin["prime_meridian_deg"] = round(w0, 9)
        spin["epoch"] = f"JD {float(jd0):.6f}"
    return spin
