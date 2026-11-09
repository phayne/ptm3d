#!/usr/bin/env python3
"""Torus forcing fidelity (Ada 2026-11-08): torus F vs the engine's true composed forcing over the 8 eq lunations
before the window (abs_PAR8_<win>.bin), per lunation and over all 8. usage: torusfid.py WIN [table]"""
import sys,numpy as np; sys.path.insert(0,'/Users/paha3326/projects/spinup_2026-11-08/tools'); import premarch as PM
K='/Users/paha3326/projects/spinup_2026-11-08'; w=sys.argv[1]
tab=PM.load_table(sys.argv[2] if len(sys.argv)>2 else f'{K}/tables/tab120_{w}.bin'); nf=tab['nf']
te=tab['ets'].ravel(); Fe=PM.torus_flux(tab,te,r2=False); A=tab['A'].reshape(nf,-1).T
print(f"K={tab['K']} NS={tab['NS']} self-reproduction max |d| {np.abs(Fe-A).max():.3g} W/m2")
raw=open(f'{K}/tables/abs_PAR8_{w}.bin','rb').read(); ns,nf=np.frombuffer(raw[8:24],'<u8')
F=np.frombuffer(raw[24:],'<f8').reshape(ns,nf); et=np.fromfile(f'{K}/tables/abs_PAR8_{w}.bin.et','<f8')
Ft=PM.torus_flux(tab,et[:960]); Ftrue=F[:960]
for L in range(8):
    sl=slice(L*120,(L+1)*120); a=Ftrue[sl].mean(0); b=Ft[sl].mean(0); d=b-a
    print(f'lunation {L}: true {a.mean():.3f} torus {b.mean():.3f} |d| max {np.abs(d).max():.2f} p99 {np.percentile(np.abs(d),99):.3f}; rel p99 (lit>1) {np.percentile(np.abs(d[a>1])/a[a>1],99) if (a>1).any() else 0:.3f}')
a=Ftrue.mean(0); b=Ft.mean(0); d=b-a; print(f'8-lunation mean: |d| max {np.abs(d).max():.2f} p99 {np.percentile(np.abs(d),99):.3f} rel p99 (lit>1) {np.percentile(np.abs(d[a>1])/a[a>1],99):.3f}')
