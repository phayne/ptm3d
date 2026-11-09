#!/usr/bin/env python3
"""Python (numba) port of the engine's DEFAULT CPU SIMD column step (Ada 2026-11-08, spin-up study).

Mirrors ThermalBatchSoA::advanceTimestep (thermal_model.cpp): per substep
  1. k = kc (1 + R350 T^3), cp = cp0 + T (c1 + T (c2 + T (c3 + T c4)))   [Hayne 2017, Moon defaults]
  2. surface node: quasi-static Newton on  F + k0 (T1 - T0)/dz0 = eps sigma T0^4   (solveSurfaceBC)
  3. Crank-Nicolson on interior nodes 1..n-2 with heat1d's g1/g2 stencil (k_{i-1}, k_i), new T0 in both halves
  4. Thomas; bottom node T_{n-1} = T_{n-2} + Qb dz_{n-2} / k_{n-2}   (PTM-853 link k)
  5. per-node floor T >= T_geo = (Qb/(eps sigma))^0.25
Substeps: per output interval dt_out, n_sub = ceil(dt_out / (P/768)), flux linear in time at substep midpoints.
Grid: the run's own depth grid (statistics/_raw/depth_grid.bin). Properties: Moon defaults, Qb from the run.
"""
import numpy as np, numba as nb

SIG = 5.670374419e-8
EPS = 0.95
KS, KD, RHOS, RHOD, HPAR, CHI = 7.4e-4, 3.4e-3, 1100.0, 1800.0, 0.06, 2.7
CP = (-3.6125, 2.7431, 2.3616e-3, -1.2340e-5, 8.9093e-9)
PDAY = 29.530588 * 86400.0


class Grid:
    def __init__(self, z, Qb):
        self.z = np.asarray(z, float); n = len(self.z)
        self.dz = np.diff(self.z)
        d3z = self.dz[1:] * self.dz[:-1] * (self.dz[1:] + self.dz[:-1])
        self.g1 = 2.0 * self.dz[1:] / d3z; self.g2 = 2.0 * self.dz[:-1] / d3z
        e = np.exp(-self.z / HPAR)
        self.kc = KD + (KS - KD) * e; self.rho = RHOD + (RHOS - RHOD) * e
        self.Qb = float(Qb); self.R350 = CHI / 350.0 ** 3
        self.Tgeo = (self.Qb / (EPS * SIG)) ** 0.25
        self.Rk = np.concatenate(([0.0], np.cumsum(self.dz / self.kc[:-1])))   # integral dz / k_c (link k)


@nb.njit(cache=True, fastmath=False)
def _substep(T, F, dt, z, dz, g1, g2, kc, rho, R350, Qb, Tgeo, a, b, c, d, k, cp, cpp, dpp):
    n = T.shape[0]
    for i in range(n):
        Ti = T[i]
        k[i] = kc[i] * (1.0 + R350 * Ti * Ti * Ti)
        cp[i] = -3.6125 + Ti * (2.7431 + Ti * (2.3616e-3 + Ti * (-1.2340e-5 + Ti * 8.9093e-9)))
    # surface Newton (solveSurfaceBC)
    es = EPS * SIG; T1 = T[1]; k0 = k[0]; dz0 = dz[0]
    Tn = T1
    if F > 0.0:
        Teq = np.sqrt(np.sqrt(F / es))
        if Teq > Tn: Tn = Teq
    for it in range(10):
        T03 = Tn * Tn * Tn
        f = F + k0 * (T1 - Tn) / dz0 - es * T03 * Tn
        df = -k0 / dz0 - 4.0 * es * T03
        hr = 4.0 * es * T03
        if it > 0 and abs(f) <= 1e-4 * hr: break
        dT = f / df
        Tnew = Tn - dT
        if Tnew < 0.5 * Tn: Tnew = 0.5 * Tn
        if Tnew > 4.0 * Tn: Tnew = 4.0 * Tn
        Tn = Tnew
    T[0] = Tn
    # CN interior
    m = n - 2
    for j in range(m):
        i = j + 1
        al = g1[j] * k[i - 1]; be = g2[j] * k[i]
        r = dt / (rho[i] * cp[i])
        a[j] = -0.5 * r * al; b[j] = 1.0 + 0.5 * r * (al + be); c[j] = -0.5 * r * be
        d[j] = T[i] + 0.5 * r * (al * T[i - 1] - (al + be) * T[i] + be * T[i + 1])
    d[0] -= a[0] * T[0]; a[0] = 0.0
    d[m - 1] -= c[m - 1] * T[n - 1]; c[m - 1] = 0.0
    cpp[0] = c[0] / b[0]; dpp[0] = d[0] / b[0]
    for j in range(1, m):
        den = b[j] - a[j] * cpp[j - 1]
        cpp[j] = c[j] / den; dpp[j] = (d[j] - a[j] * dpp[j - 1]) / den
    T[m] = dpp[m - 1]
    for j in range(m - 2, -1, -1):
        T[j + 1] = dpp[j] - cpp[j] * T[j + 2]
    T[n - 1] = T[n - 2] + Qb * dz[n - 2] / k[n - 2]
    for i in range(n):
        if T[i] < Tgeo: T[i] = Tgeo


@nb.njit(parallel=True, cache=True)
def march(T, Fgrid, dts, z, dz, g1, g2, kc, rho, R350, Qb, Tgeo, target):
    """Advance T (nf, nz) through the output grid: Fgrid (nt+1, nf) flux at the interval ends, dts (nt,)."""
    nf, nz = T.shape
    nt = dts.shape[0]
    for f in nb.prange(nf):
        a = np.empty(nz); b = np.empty(nz); c = np.empty(nz); d = np.empty(nz)
        k = np.empty(nz); cp = np.empty(nz); cpp = np.empty(nz); dpp = np.empty(nz)
        Tc = T[f].copy()
        for t in range(nt):
            dto = dts[t]
            ns = int(np.ceil(dto / target))
            if ns < 1: ns = 1
            h = dto / ns
            Fa = Fgrid[t, f]; Fb = Fgrid[t + 1, f]
            for s in range(ns):
                al = (s + 0.5) / ns
                _substep(Tc, (1.0 - al) * Fa + al * Fb, h, z, dz, g1, g2, kc, rho, R350, Qb, Tgeo,
                         a, b, c, d, k, cp, cpp, dpp)
        T[f] = Tc
    return T


def advance(G, T, Fgrid, dts):
    return march(np.ascontiguousarray(T, float), np.ascontiguousarray(Fgrid, float), np.asarray(dts, float),
                 G.z, G.dz, G.g1, G.g2, G.kc, G.rho, G.R350, G.Qb, G.Tgeo, PDAY / 768.0)


# ---------------------------------------------------------------- state files (PTMT v1 body + v2 trailer)
import struct, zlib


def read_state(path):
    raw = open(path, 'rb').read()
    magic, nf, nl = struct.unpack('<III', raw[:12])
    assert magic == 0x50544D54, 'not a PTMT file'
    T = np.frombuffer(raw[12:12 + 4 * nf * nl], '<f4').reshape(nl, nf).T.astype(float)   # (nf, nl)
    return T, raw


def write_state_like(template_raw, T, path):
    """Write T (nf, nl) into a copy of a PTMT v2 file (same header and trailer), CRC32 recomputed."""
    magic, nf, nl = struct.unpack('<III', template_raw[:12])
    body = np.ascontiguousarray(np.asarray(T, float).T, dtype='<f4').tobytes()
    assert len(body) == 4 * nf * nl
    rest = template_raw[12 + 4 * nf * nl:]
    pre = template_raw[:12] + body + rest[:-4]
    crc = zlib.crc32(pre) & 0xFFFFFFFF
    open(path, 'wb').write(pre + struct.pack('<I', crc))
