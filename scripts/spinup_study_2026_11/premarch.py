#!/usr/bin/env python3
"""Rank-1 seasonal PRE-MARCH seed, Python prototype (Ada 2026-11-08; design: Fermi 2026-11-08).

  premarch.py WIN(203009|203102) NFINE [--years Y] [--project | --no-project] [--no-r2] [--tag TAG]

1. Forcing torus F(phi_season, tau_local) per face from the engine's Tier-2 pre-pass table
   (PTM_DUMP_SEASONAL_TABLE: K phases x NS samples over the eclipse year ending at the window start; composed
   absorbed flux, Earthshine as phase means, no indirect -- radiation off here).
   Seasonal phase: linear between adjacent phase centers (cyclic). Local time: the sample with the same ABSOLUTE
   synodic phase inside each phase's own lunation (linear between samples, cyclic).
   Real 1/r^2: each table sample is rescaled from the Earth-Sun distance at ITS epoch to the distance at the
   marching epoch (Keplerian Earth orbit; the Moon's own offset recurs with the synodic phase).
2. March the engine's own column (col_engine.py, parity 3.8e-6 K = float32 state quantization) on a virtual
   grid of P/120 (the engine eq cadence, P/768 substeps) through Y table periods ending at t_inj = W - NFINE*P,
   from a plain conductive seed at the table-mean flux.
3. Slow-mode projection (default on): at the end of period Y-1, Aitken extrapolation of the deep-node series at
   the period boundaries (the seasonal wave is exactly periodic on the torus there; the diurnal wave is ~0 at
   0.47 m), applied as a uniform column shift; then one more period. Faces with a non-geometric sequence are not
   projected (counted).
4. Write the state at t_inj into a copy of the truth dump TRD<N>_<win>.bin (same epoch/trailer; CRC recomputed),
   for the engine run "P": equilibration_orbits = NFINE with PTM_SEED_STATE_FILE.
"""
import sys, time, json, argparse, numpy as np
sys.path.insert(0, '/Users/paha3326/projects/spinup_2026-11-08/tools')
import col_engine as C, score as S

K = '/Users/paha3326/projects/spinup_2026-11-08'
WET = {'203009': 967723269.182627, '203102': 982411269.185171}
PSYN = C.PDAY
# Keplerian Earth orbit (J2000): perihelion 2000-01-03 05:18 UT ~ ET +148,700 s; anomalistic year; eccentricity
T_PERI, Y_ANOM, ECC = 148700.0, 365.259636 * 86400.0, 0.0167086


def inv_r2(et):
    M = 2 * np.pi * ((np.asarray(et) - T_PERI) / Y_ANOM)
    E = M.copy()
    for _ in range(6): E = E - (E - ECC * np.sin(E) - M) / (1 - ECC * np.cos(E))
    return 1.0 / (1 - ECC * np.cos(E)) ** 2


def load_table(path):
    raw = open(path, 'rb').read(); o = 8
    Kp, NS, nf = np.frombuffer(raw[o:o + 12], '<i4'); o += 12
    period, cycle, tau = np.frombuffer(raw[o:o + 24], '<f8'); o += 24
    ets = np.frombuffer(raw[o:o + 8 * Kp * NS], '<f8').reshape(Kp, NS); o += 8 * Kp * NS
    Sin = np.frombuffer(raw[o:o + 8 * Kp * NS], '<f8').reshape(Kp, NS); o += 8 * Kp * NS
    A = np.frombuffer(raw[o:], '<f4').reshape(nf, Kp, NS).astype(float)
    return dict(K=int(Kp), NS=int(NS), nf=int(nf), period=float(period), cycle=float(cycle), ets=ets, S=Sin, A=A)


def torus_flux(tab, t, r2=True):
    """F (len(t), nf) on the torus at epochs t. For each phase k, the sample with the SAME absolute synodic phase
    as t (linear between samples of that phase's lunation) sits at epoch s_k(t), seasonal coordinate
    u_k = (s_k - tc0) mod P; the flux is linear in season between the two u_k that bracket t's own coordinate."""
    Kp, NS, P = tab['K'], tab['NS'], tab['period']; ets, A = tab['ets'], tab['A']; C_ = tab['cycle']
    tc0 = ets[0].mean(); dts = C_ / NS; nt = len(t)
    pos = np.mod(t[:, None] - ets[None, :, 0], C_) / dts                 # (nt, K)
    i0 = np.floor(pos).astype(int) % NS; wi = pos - np.floor(pos); i1 = (i0 + 1) % NS
    seam = i0 == NS - 1                     # never interpolate across the phase lunation's seam (1 lunation of season)
    wi = np.where(seam, 0.0, wi); i1 = np.where(seam, i0, i1)
    sk = ets[None, :, 0] + pos * dts                                       # (nt, K) sample epochs
    uk = sk - tc0                                                          # in [-C/2, P + C/2) roughly, increasing in k
    ut = np.mod(t - tc0, P)
    # bracket: k = last phase with uk <= ut (cyclic wrap through P)
    kk = (uk <= ut[:, None]).sum(1) - 1
    lo = np.where(kk < 0, Kp - 1, kk); hi = np.where(kk < 0, 0, (kk + 1) % Kp)
    r = np.arange(nt)
    ulo = uk[r, lo] - np.where(kk < 0, P, 0.0); uhi = uk[r, hi] + np.where(kk == Kp - 1, P, 0.0)
    w = np.clip((ut - ulo) / (uhi - ulo), 0.0, 1.0)
    corr_now = inv_r2(t) if r2 else np.ones(nt)
    out = np.zeros((nt, tab['nf']))
    for kidx, ww in ((lo, 1.0 - w), (hi, w)):
        a0, a1, wl = i0[r, kidx], i1[r, kidx], wi[r, kidx]
        Fk = (1 - wl)[:, None] * A[:, kidx, a0].T + wl[:, None] * A[:, kidx, a1].T
        if r2:   # rescale from the sample epoch's Earth-Sun distance to the marching epoch's
            ce = (1 - wl) * inv_r2(ets[kidx, a0]) + wl * inv_r2(ets[kidx, a1])
            Fk *= (corr_now / ce)[:, None]
        out += ww[:, None] * Fk
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('win'); ap.add_argument('nfine', type=float)
    ap.add_argument('--years', type=int, default=4); ap.add_argument('--no-project', action='store_true')
    ap.add_argument('--no-r2', action='store_true'); ap.add_argument('--tag', default='P')
    ap.add_argument('--table', default=None)
    a = ap.parse_args()
    tab = load_table(a.table or f'{K}/tables/tab120_{a.win}.bin')
    P = tab['period']; Wet = WET[a.win]; N = int(a.nfine)
    t_inj = Wet - a.nfine * PSYN
    z = S.grid(f'{K}/runs/TR64_{a.win}'); G = C.Grid(z, 0.022)
    # seed: plain conductive profile at the table-mean flux
    Fm = tab['A'].reshape(tab['nf'], -1).mean(1)
    T0s = ((Fm + G.Qb) / (C.EPS * C.SIG)) ** 0.25
    T = T0s[:, None] + G.Qb * G.Rk[None, :]
    nper = int(np.ceil(P / (PSYN / 120.0))); dtp = P / nper     # virtual eq-cadence grid, whole periods
    t0 = time.time(); deep = [T[:, -1].copy()]; nproj = 0; proj_info = {}
    for y in range(a.years):
        ts = t_inj - (a.years - y) * P + dtp * np.arange(nper + 1)
        F = torus_flux(tab, ts, r2=not a.no_r2)
        T = C.advance(G, T, F, np.full(nper, dtp))
        deep.append(T[:, -1].copy())
        if (not a.no_project) and y == a.years - 2 and len(deep) >= 3:
            d1 = deep[-2] - deep[-3]; d2 = deep[-1] - deep[-2]
            with np.errstate(divide='ignore', invalid='ignore'):
                rho = d2 / d1
            ok = np.isfinite(rho) & (rho > 0.0) & (rho < 0.9) & (np.abs(d1) > 1e-6)
            shift = np.where(ok, d2 * rho / (1.0 - rho), 0.0)
            T = T + shift[:, None]; T = np.maximum(T, G.Tgeo)
            nproj = int(ok.sum())
            proj_info = dict(n_projected=nproj, shift_max=float(np.abs(shift).max()),
                             shift_p99=float(np.percentile(np.abs(shift), 99)),
                             rho_med=float(np.median(rho[ok])) if ok.any() else None)
        print(f'  period {y+1}/{a.years}: deep-node change max {np.abs(deep[-1]-deep[-2]).max():.4f} K, '
              f'p99 {np.percentile(np.abs(deep[-1]-deep[-2]), 99):.4f} K  ({time.time()-t0:.0f} s)', flush=True)
    tmpl = open(f'{K}/states/TRD{N}_{a.win}.bin', 'rb').read()
    out = f'{K}/states/{a.tag}{N}_{a.win}.bin'
    C.write_state_like(tmpl, T, out)
    Ttr, _ = C.read_state(f'{K}/states/TRD{N}_{a.win}.bin')
    e = T - Ttr
    info = dict(win=a.win, nfine=a.nfine, years=a.years, project=not a.no_project, r2=not a.no_r2, proj=proj_info,
                wall_s=time.time() - t0, seed_vs_truth_at_tinj=dict(
                    surf_max=float(np.abs(e[:, 0]).max()), surf_p99=float(np.percentile(np.abs(e[:, 0]), 99)),
                    deep_max=float(np.abs(e[:, -1]).max()), deep_p99=float(np.percentile(np.abs(e[:, -1]), 99)),
                    deep_mean=float(e[:, -1].mean()), all_max=float(np.abs(e).max())))
    json.dump(info, open(out + '.json', 'w'), indent=1)
    print(json.dumps(info))


if __name__ == '__main__':
    main()
