#!/usr/bin/env python3
"""Real-year REPLAY pre-march, Python prototype (Ada 2026-11-08). Alternative to the (season x local-time) torus.

  replay.py WIN NFINE [--years-block 1|2] [--replays R] [--no-project] [--tag TAG]

Block = the engine's own composed forcing (PTM_DUMP_ABSORBED of an eq-24 run) over the last B eclipse years
before the window, ON THE RUN'S OWN EQ GRID (so the samples are the truth's samples: no quadrature mismatch).
The block is replayed as a periodic forcing (one seam per block, linear over one step) R times from a plain
conductive seed at the block-mean flux, with an Aitken projection of the whole column at the block boundaries
(exact for periodic forcing: the slow mode decays geometrically from boundary to boundary), then marched from
the block start to t_inj = W - NFINE*P and written into a copy of TRD<N>_<win>.bin for the engine (eq NFINE).
"""
import sys, time, json, argparse, numpy as np
sys.path.insert(0, '/Users/paha3326/projects/spinup_2026-11-08/tools')
import col_engine as C, score as S

K = '/Users/paha3326/projects/spinup_2026-11-08'
PYR = 346.62 * 86400.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('win'); ap.add_argument('nfine', type=int)
    ap.add_argument('--years-block', type=int, default=1); ap.add_argument('--replays', type=int, default=4)
    ap.add_argument('--no-project', action='store_true'); ap.add_argument('--tag', default='R')
    ap.add_argument('--coarse', type=int, default=1, help='keep every m-th sample of the block OLDER than t_inj')
    a = ap.parse_args()
    raw = open(f'{K}/tables/abs_E24_{a.win}.bin', 'rb').read(); ns, nf = np.frombuffer(raw[8:24], '<u8')
    F = np.frombuffer(raw[24:], '<f8').reshape(ns, nf); et = np.fromfile(f'{K}/tables/abs_E24_{a.win}.bin.et', '<f8')
    neq = 2880; dt = et[1] - et[0]
    nb = int(round(a.years_block * PYR / dt))            # block length in eq steps
    i0 = neq - nb
    assert i0 >= 0
    if a.coarse > 1:   # coarse cadence on the part of the block older than t_inj (the fine N lunations stay exact)
        k_inj = nb - 120 * a.nfine
        keep = list(range(0, k_inj, a.coarse)) + list(range(k_inj, nb))
        ts = np.array(keep, float)
        Fb = np.vstack([F[i0 + np.array(keep)], F[i0:i0 + 1]]); dts = np.diff(np.append(ts, nb)) * dt
        print(f'  coarse: {len(keep)} block samples ({k_inj // a.coarse} coarse + {120 * a.nfine} fine) of {nb}')
    else:
        Fb = np.vstack([F[i0:neq], F[i0:i0 + 1]])            # nb intervals, periodic (seam: last sample -> first)
        dts = np.full(nb, dt)
    z = S.grid(f'{K}/runs/TR64_{a.win}'); G = C.Grid(z, 0.022)
    Fm = F[i0:neq].mean(0)
    T = (((Fm + G.Qb) / (C.EPS * C.SIG)) ** 0.25)[:, None] + G.Qb * G.Rk[None, :]
    t0 = time.time(); hist = [T.copy()]; info = {}
    for r in range(a.replays):
        T = C.advance(G, T, Fb, dts); hist.append(T.copy())
        d = hist[-1] - hist[-2]
        print(f'  replay {r+1}/{a.replays}: boundary change max {np.abs(d).max():.4f} p99 {np.percentile(np.abs(d), 99):.4f} K '
              f'({time.time()-t0:.0f} s)', flush=True)
        if (not a.no_project) and r == a.replays - 2 and len(hist) >= 3:
            d1 = (hist[-2] - hist[-3])[:, -1]; d2 = (hist[-1] - hist[-2])[:, -1]   # deep node sets the ratio
            with np.errstate(divide='ignore', invalid='ignore'):
                rho = d2 / d1
            ok = np.isfinite(rho) & (rho > 0) & (rho < 0.9) & (np.abs(d1) > 1e-6)
            ext = np.where(ok, rho / (1 - rho), 0.0)
            T = np.maximum(T + ext[:, None] * (hist[-1] - hist[-2]), G.Tgeo)
            info = dict(n_projected=int(ok.sum()), rho_med=float(np.median(rho[ok])) if ok.any() else None,
                        shift_max=float(np.abs(ext[:, None] * (hist[-1] - hist[-2])).max()))
    # final partial pass: block start -> t_inj (= index neq - 120*N)
    k = (nb - 120 * a.nfine) if a.coarse == 1 else len(range(0, nb - 120 * a.nfine, a.coarse))
    T = C.advance(G, T, Fb[:k + 1], dts[:k])
    tmpl = open(f'{K}/states/TRD{a.nfine}_{a.win}.bin', 'rb').read()
    out = f'{K}/states/{a.tag}{a.nfine}_{a.win}.bin'
    C.write_state_like(tmpl, T, out)
    Ttr, _ = C.read_state(f'{K}/states/TRD{a.nfine}_{a.win}.bin'); e = T - Ttr
    res = dict(win=a.win, nfine=a.nfine, years_block=a.years_block, replays=a.replays, proj=info,
               wall_s=time.time() - t0, seed_vs_truth=dict(surf_max=float(np.abs(e[:, 0]).max()),
               surf_p99=float(np.percentile(np.abs(e[:, 0]), 99)), deep_max=float(np.abs(e[:, -1]).max()),
               deep_p99=float(np.percentile(np.abs(e[:, -1]), 99)), deep_mean=float(e[:, -1].mean())))
    json.dump(res, open(out + '.json', 'w'), indent=1); print(json.dumps(res))


if __name__ == '__main__':
    main()
