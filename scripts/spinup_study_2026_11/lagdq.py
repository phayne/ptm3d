#!/usr/bin/env python3
"""Lagged-vs-converged indirect flux (Ada 2026-11-09; Fermi step 1, no new physics).
usage: lagdq.py WIN
TRUTH = series/TRf64_<win>.ids (eq 64, converged, full-span re-seed): the last 3000 steps are the E24 grid.
LAG   = series/LAG1_E24_<win>.ids (eq 24, ONE coupling sweep: the series a frozen history would carry).
Reports, per face class (from the radiation-off truth TR64_<win>) and per lunation of the 24-lunation history:
  dq/q = <q_lag>/<q_conv> - 1 of the lunation-mean (and history-mean) indirect flux, for faces with <q_conv> > 0.05 W/m2;
  candidate (a): rescale the lag per face by <q_conv>/<q_lag> measured on the last N_fine eq lunations (N_fine = 2, 4, 8);
  split into the visible and IR bands.
Pass guide (Fermi): |dq/q| <~ 2.5% for N_fine = 4, <~ 4% for N_fine = 8."""
import sys, numpy as np
sys.path.insert(0, '/Users/paha3326/projects/spinup_2026-11-08/tools'); import score as S
K = '/Users/paha3326/projects/spinup_2026-11-08'


def load(p):
    raw = open(p, 'rb').read()
    assert raw[:8] == b'PTMIDS01'
    nF, nT, lo = np.frombuffer(raw[8:32], '<u8'); bands = raw[32]
    a = np.frombuffer(raw[33:], '<f4'); n = int(nF * nT)
    out = {'tot': a[:n].reshape(nF, nT), 'lo': int(lo)}
    if bands: out['vis'] = a[n:2 * n].reshape(nF, nT); out['ir'] = a[2 * n:3 * n].reshape(nF, nT)
    return out


w = sys.argv[1]
Tr = load(f'{K}/series/TRf64_{w}.ids'); Lg = load(f'{K}/series/' + (sys.argv[2] if len(sys.argv) > 2 else f'LAG1_E24_{w}') + '.ids')
nE = Lg['tot'].shape[1]; assert nE == 3000, nE
off = Tr['tot'].shape[1] - nE
cls = S.classes(f'{S.R}/TR64_{w}')
P = 120; nl = 24   # eq lunations of E24 (steps 0..2879), then the 120-step window


def lun(a, b0, b1):  # mean over eq lunations [b0, b1)
    return a[:, b0 * P:b1 * P].astype(float).mean(1)


for band in ('tot', 'ir', 'vis'):
    if band not in Lg: continue
    C = Tr[band][:, off:]; L = Lg[band]
    print(f'== {w} band {band}')
    for N in (2, 4, 8):
        hist = slice(0, (nl - N) * P)
        qc = C[:, hist].astype(float).mean(1); ql = L[:, hist].astype(float).mean(1)
        ok = qc > 0.05
        dq = ql / np.where(ok, qc, 1) - 1
        # candidate (a): per-face ratio from the last N fine eq lunations
        rc = lun(C, nl - N, nl); rl = lun(L, nl - N, nl)
        ratio = np.where(rl > 1e-6, rc / np.maximum(rl, 1e-12), 1.0)
        dqa = ql * ratio / np.where(ok, qc, 1) - 1
        s = f'  N_fine {N}: history = first {nl - N} lunations;'
        for nm, m in cls.items():
            mm = m & ok
            if mm.sum() < 5: s += f' {nm} n{mm.sum()};'; continue
            a, b = np.abs(dq[mm]), np.abs(dqa[mm])
            s += (f' {nm} n{mm.sum()} lag med {np.median(dq[mm])*100:+.1f}% |p90| {np.percentile(a,90)*100:.1f}% |p99| '
                  f'{np.percentile(a,99)*100:.1f}% | rescaled med {np.median(dqa[mm])*100:+.1f}% |p90| {np.percentile(b,90)*100:.1f}% '
                  f'|p99| {np.percentile(b,99)*100:.1f}%;')
        print(s)
    # per-lunation (season) mesh median of dq/q, all faces with q > 0.05
    row = []
    for k in range(nl):
        qc = lun(C, k, k + 1); ql = lun(L, k, k + 1); ok = qc > 0.05
        row.append(f'{np.median(ql[ok]/qc[ok]-1)*100:+.0f}' if ok.sum() > 5 else 'na')
    print('  per-lunation median dq/q [%] (lunation 0 = W-24P): ' + ' '.join(row))
