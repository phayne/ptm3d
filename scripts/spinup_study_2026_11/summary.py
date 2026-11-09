#!/usr/bin/env python3
"""Markdown summary table of spin-up arms vs truth (Ada 2026-11-08). usage: summary.py TRUTH L8 ARM [ARM ...]
Columns per class (never-lit / lit-earlier / lit-now): max and p99 of |e| over surface T_min, T_avg, T_max
(the worst of the three quantities); 7 cm T_avg max / p99; class flips 54/110 K; ice max |dz|; never-worse vs L8."""
import sys, os
sys.path.insert(0, '/Users/paha3326/projects/spinup_2026-11-08/tools')
import score as S

args = sys.argv[1:]
cls_run = None
if '--cls' in args:   # radiation-on: E_abs includes indirect flux, so classes come from the radiation-off truth
    i = args.index('--cls'); cls_run = args[i + 1]; del args[i:i + 2]
truth, l8 = args[0], args[1]
T = f'{S.R}/{truth}'; cls = S.classes(f'{S.R}/{cls_run}' if cls_run else T)
if cls_run: print(f'(classes from the radiation-off truth `{cls_run}`: direct + Earthshine only)')
print(f'Truth `{truth}`: ' + ', '.join(f'{k} {int(v.sum())}' for k, v in cls.items()))
print()
print('| arm | never-lit max / p99 | lit-earlier max / p99 | lit-now max / p99 | 7 cm T_avg max / p99 | flips 54/110 | ice max dz [cm] | no worse than L8 | PASS |')
print('|---|---|---|---|---|---|---|---|---|')
for a in args[2:]:
    if not os.path.isdir(f'{S.R}/{a}/statistics'):
        print(f'| {a} | missing |'); continue
    o = S.score(T, f'{S.R}/{a}', cls, f'{S.R}/{l8}')
    cells = []
    for nm in cls:
        mx = max(o[f'{nm}:T{q}'][0] for q in ('min', 'avg', 'max')); p9 = max(o[f'{nm}:T{q}'][1] for q in ('min', 'avg', 'max'))
        cells.append(f'{mx:.2f} / {p9:.2f}')
    c7 = o['7cm:Tavg']
    print(f"| {a} | " + ' | '.join(cells) + f" | {c7[0]:.2f} / {c7[1]:.2f} | {o['flips54']}/{o['flips110']} | "
          f"{o['ice_dz_cm']:.3f} | {o['noworse_vs_L8']*100:.1f}% | {'yes' if o['PASS'] else 'no'} |")
