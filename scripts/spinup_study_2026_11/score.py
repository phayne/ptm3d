#!/usr/bin/env python3
"""Spin-up decisive-test scorer (Ada 2026-11-08).

score.py TRUTH ARM [ARM ...] [--L8 L8ARM] [--json out.json]
Per ARM vs TRUTH (same output window):
  surface T_min/T_avg/T_max per-face max |e| and p99 |e|, split by class:
    never-lit   : E_abs == 0 over the truth run's WHOLE history (64 eq lunations + window), from the truth ledger
    lit-earlier : E_abs > 0 in the history but == 0 in the window
    lit-now     : E_abs > 0 in the window
  7 cm T_avg max/p99 (all faces); class flips at 54 / 110 K on surface T_max; ice max |dz| [cm] and node changes;
  never-worse vs L8 (fraction of faces where max_q |e_arm| <= max_q |e_L8| + 0.02 K).
PASS bars (Sagan/Fermi): surface max <= 0.5 and p99 <= 0.25 per class; 7 cm T_avg max <= 0.5; 0 flips;
ice max|dz| <= 0.1 cm.
Stats: '<f4' with a one-word header (dropped).
"""
import sys, os, json, numpy as np

R = '/Users/paha3326/projects/spinup_2026-11-08/runs'


def st(d, q): return np.fromfile(f'{d}/statistics/{q}.bin', '<f4')[1:].astype(float)


def grid(d):
    raw = open(f'{d}/statistics/_raw/depth_grid.bin', 'rb').read(); return np.frombuffer(raw[4:], '<f8')


def classes(truth):
    import csv
    nf = len(st(truth, 'temperature_avg'))
    tot = np.zeros(nf); win = np.zeros(nf)
    L = np.genfromtxt(f'{truth}.ledger.csv', delimiter=',', names=True, usecols=('face', 'cycle', 'E_abs'))
    f = L['face'].astype(int); c = L['cycle']; e = L['E_abs']
    np.add.at(tot, f, e); np.add.at(win, f[c >= 1000], e[c >= 1000])
    return {'never-lit': tot == 0, 'lit-earlier': (tot > 0) & (win == 0), 'lit-now': win > 0}


def err(a, t, q): return st(a, q) - st(t, q)


def score(truth, arm, cls, l8=None):
    out = {'arm': os.path.basename(arm)}
    worst = np.zeros(len(st(truth, 'temperature_avg')))
    for q in ('min', 'avg', 'max'):
        e = err(arm, truth, f'temperature_{q}'); worst = np.maximum(worst, np.abs(e))
        for nm, m in cls.items():
            ax = np.abs(e[m])
            out[f'{nm}:T{q}'] = ((float(ax.max()), float(np.percentile(ax, 99)), float(e[m].mean())) if m.any()
                                 else (0.0, 0.0, 0.0))
    if os.path.exists(f'{arm}/statistics/temperature_z0p071m_avg.bin'):
        e7 = err(arm, truth, 'temperature_z0p071m_avg'); ax = np.abs(e7)
        out['7cm:Tavg'] = (float(ax.max()), float(np.percentile(ax, 99)), float(e7.mean()))
    else:
        out['7cm:Tavg'] = (float('nan'), float('nan'), float('nan'))   # no depth maps on this config
    tA, tT = st(arm, 'temperature_max'), st(truth, 'temperature_max')
    out['flips54'] = int(((tA < 54) != (tT < 54)).sum()); out['flips110'] = int(((tA < 110) != (tT < 110)).sum())
    if not os.path.exists(f'{truth}/statistics/_raw/depth_grid.bin'):
        out['ice_dz_cm'] = float(np.abs(st(arm, 'ice_depth') - st(truth, 'ice_depth')).max() * 100); out['ice_nodechg'] = -1
    else:
      z = grid(truth).astype(np.float32); Zt = st(truth, 'ice_depth'); Za = st(arm, 'ice_depth')
      out['ice_dz_cm'] = float(np.abs(Za - Zt).max() * 100)
      out['ice_nodechg'] = int((np.searchsorted(z, Zt, 'right') != np.searchsorted(z, Za, 'right')).sum())
    if l8 is not None:
        w8 = np.zeros_like(worst)
        for q in ('min', 'avg', 'max'): w8 = np.maximum(w8, np.abs(err(l8, truth, f'temperature_{q}')))
        out['noworse_vs_L8'] = float(np.mean(worst <= w8 + 0.02))
    surf_ok = all(out[f'{nm}:T{q}'][0] <= 0.5 and out[f'{nm}:T{q}'][1] <= 0.25
                  for nm in cls for q in ('min', 'avg', 'max') if cls[nm].any())
    out['PASS'] = bool(surf_ok and not (out['7cm:Tavg'][0] > 0.5) and out['flips54'] == 0 and out['flips110'] == 0
                       and out['ice_dz_cm'] <= 0.1)
    return out


def fmt(o, cls):
    s = f"{o['arm']:16s}"
    for nm in cls:
        mx = max(o[f'{nm}:T{q}'][0] for q in ('min', 'avg', 'max'))
        p9 = max(o[f'{nm}:T{q}'][1] for q in ('min', 'avg', 'max'))
        s += f" | {nm[:6]} max {mx:5.2f} p99 {p9:5.2f}"
    s += (f" | 7cm max {o['7cm:Tavg'][0]:5.2f} p99 {o['7cm:Tavg'][1]:5.2f} | flips {o['flips54']}/{o['flips110']}"
          f" | ice {o['ice_dz_cm']:.3f} cm ({o['ice_nodechg']})")
    if 'noworse_vs_L8' in o: s += f" | nw(L8) {o['noworse_vs_L8']*100:5.1f}%"
    return s + (' PASS' if o['PASS'] else '')


if __name__ == '__main__':
    args = sys.argv[1:]; l8 = None; js = None
    if '--L8' in args: i = args.index('--L8'); l8 = f'{R}/{args[i+1]}'; del args[i:i+2]
    if '--json' in args: i = args.index('--json'); js = args[i+1]; del args[i:i+2]
    truth = f'{R}/{args[0]}'; cls = classes(truth)
    print(f"truth {args[0]}: " + ', '.join(f'{k} {v.sum()}' for k, v in cls.items()))
    res = []
    for a in args[1:]:
        if not os.path.isdir(f'{R}/{a}/statistics'): print(f'{a:16s} MISSING'); continue
        o = score(truth, f'{R}/{a}', cls, l8); res.append(o); print(fmt(o, cls))
    if js: json.dump(res, open(js, 'w'), indent=1)
