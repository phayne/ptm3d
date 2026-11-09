#!/usr/bin/env python3
"""Wall-time split per arm from PTM_TIMING_JSON (Ada 2026-11-08). usage: timing.py ARM [ARM ...]
radiosity = radiosity_sweep; thermal = thermal_whole_mesh + rethermal; illumination; total; coupling sweeps."""
import sys, json, os
R = '/Users/paha3326/projects/spinup_2026-11-08/runs'
print('| arm | total s | radiosity s | thermal s (whole-mesh + re-thermal) | illumination s | sweeps |')
print('|---|---|---|---|---|---|')
for a in sys.argv[1:]:
    p = f'{R}/{a}.timing.json'
    if not os.path.exists(p): print(f'| {a} | missing |'); continue
    d = json.load(open(p)); iv = d['intervals']
    g = lambda k: iv.get(k, {}).get('wall_s', 0.0)
    tot = d['total'].get('wall_s', 0.0) if isinstance(d.get('total'), dict) else d.get('total', 0.0)
    ci = d.get('coupling_iterations'); n = len(ci) if isinstance(ci, list) else ci
    print(f"| {a} | {tot:.0f} | {g('radiosity_sweep'):.0f} | {g('thermal_whole_mesh') + g('rethermal'):.0f} | "
          f"{g('illumination'):.1f} | {n} |")
