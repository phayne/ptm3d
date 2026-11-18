"""Reference implementation of the ice retention criterion, for tests only.

Provenance. This is the closed-form integrator written by the reviewer of the
retention criterion on 2026-11-18 as a check of the engine's numerics. It was
written separately from ``src/cpp/ice/ice_retention.cpp`` and shares no code
with it; this copy is a port from NumPy to the standard library, with the
method unchanged. It checks the numerics of the stated model, not the physics:
no data-anchored validation of ice retention exists.

Model (the engine's): pore-filling ice m(z) = rho_ice (1 - rho(z)/rho_g), with
rho(z) = rho_d + (rho_s - rho_d) exp(-z/H); constant transport length ell, so
Lambda(z) = z/ell; <E> log-linear in depth between nodes; the time to bring the
ice table to z is t(z) = int_0^z m (1 + 2 Lambda) / (mu <E>) dz'.

Method (independent of the engine's choices): ``nsub`` equal sub-steps per node
interval (default 64, the engine uses 16) with m at each sub-step midpoint; on
a sub-step the integrand is (A + B x) exp(-l_a - c x), integrated exactly; the
root t(z_tau) = tau by 60 bisections of the exact partial integral (the engine
uses Newton).

Standard library only, so that the end-to-end test needs no NumPy.
"""
import math

MU = 2.99e-26                       # kg, one water molecule
GYR = 1e9 * 365.25 * 86400.0        # s
RATE_FLOOR = 1e-200


def _partial(ea, A, B, c, x):
    """int_0^x (A + B s) exp(-l_a - c s) ds, ea = exp(-l_a)."""
    u = c * x
    if abs(u) < 1e-4:
        return ea * (A * x * (1 - u / 2 + u * u / 6) + B * x * x * (0.5 - u / 3 + u * u / 8))
    e = math.exp(-u)
    return ea * ((A / c) * (1 - e) + B * ((1 - e) / (c * c) - x * e / c))


def retention_column(E, z, m_of_z, ell, age_gyr=1.0, nsub=64):
    """One column. E: time-mean rate per node (molecules m^-2 s^-1); z: node
    depths (z[0] = 0). Returns (z_tau or None when the column is exhausted,
    t_top in s, t_column in s)."""
    T = age_gyr * GYR
    lE = [math.log(max(e, RATE_FLOOR)) for e in E]
    z_half = 0.5 * z[1]
    t = 0.0
    z_tau = None
    t_top = None
    for i in range(len(z) - 1):
        dz = z[i + 1] - z[i]
        h = dz / nsub
        c = (lE[i + 1] - lE[i]) / dz
        for k in range(nsub):
            za = z[i] + k * h
            m = m_of_z(za + 0.5 * h)
            A = m * (1 + 2 * za / ell) / MU
            B = m * 2 / ell / MU
            ea = math.exp(-(lE[i] + c * k * h))
            dt = _partial(ea, A, B, c, h)
            if z_tau is None and t + dt >= T:
                lo, hi = 0.0, h
                for _ in range(60):
                    x = 0.5 * (lo + hi)
                    if t + _partial(ea, A, B, c, x) >= T:
                        hi = x
                    else:
                        lo = x
                z_tau = za + 0.5 * (lo + hi)
            t += dt
            if i == 0 and k + 1 == nsub // 2:
                t_top = t
    if t_top is None:                   # nsub odd: integrate to z_1/2 directly
        t_top = _top_lifetime(lE, z, m_of_z, ell, nsub)
    return z_tau, t_top, t


def _top_lifetime(lE, z, m_of_z, ell, nsub):
    zh, t = 0.5 * z[1], 0.0
    c = (lE[1] - lE[0]) / z[1]
    h = zh / nsub
    for k in range(nsub):
        za = k * h
        m = m_of_z(za + 0.5 * h)
        t += _partial(math.exp(-(lE[0] + c * za)), m * (1 + 2 * za / ell) / MU, m * 2 / ell / MU, c, h)
    return t


def pore_filling(rho_s, rho_d, H, rho_g, rho_ice=917.0):
    """m(z) of the engine's pore-filling inventory."""
    return lambda zz: rho_ice * (1 - (rho_d + (rho_s - rho_d) * math.exp(-zz / H)) / rho_g)


def retention_faces(sum_e, n, nl, nf, z, rec, age_gyr=1.0, nsub=64):
    """Every face of a run's saved sums (layer-major, as in
    statistics/_raw/ice_depth_sum_E.bin), with the run's recorded settings
    (statistics/_raw/ice_depth_params.txt as a dict). Returns per face
    (depth, class 0 top layer / 1 at depth / 2 none, t_top in Gyr). Class 3
    (surface-rate underflow) is reported as 0."""
    m_of_z = pore_filling(*(float(rec[k]) for k in ("rho_surface_kgm3", "rho_deep_kgm3",
                                                    "density_scale_m", "grain_density_kgm3")),
                          rho_ice=float(rec["ice_density_kgm3"]))
    ell = float(rec["grain_hop_length_m"])
    zh = 0.5 * z[1]
    out = []
    for f in range(nf):
        E = [sum_e[j * nf + f] / n for j in range(nl)]
        zt, ttop, _ = retention_column(E, z, m_of_z, ell, age_gyr, nsub)
        cls = 2 if zt is None else (0 if zt < zh else 1)
        out.append((zt if zt is not None else z[-1], cls, ttop / GYR))
    return out
