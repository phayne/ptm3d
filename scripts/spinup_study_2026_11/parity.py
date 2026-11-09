#!/usr/bin/env python3
"""Column parity (Ada 2026-11-08): Python port vs engine on the engine's own dumped forcing.
From the truth state at W - 8P (TRD8), march the 960 eq steps of abs_PAR8 and compare with the engine's state at W."""
import sys, time, numpy as np
sys.path.insert(0, '/Users/paha3326/projects/spinup_2026-11-08/tools'); import col_engine as C, score as S
K = '/Users/paha3326/projects/spinup_2026-11-08'
raw = open(f'{K}/tables/abs_PAR8_203009.bin', 'rb').read(); ns, nf = np.frombuffer(raw[8:24], '<u8')
F = np.frombuffer(raw[24:], '<f8').reshape(ns, nf); et = np.fromfile(f'{K}/tables/abs_PAR8_203009.bin.et', '<f8')
z = S.grid(f'{K}/runs/TR64_203009'); G = C.Grid(z, 0.022)
T0, _ = C.read_state(f'{K}/states/TRD8_203009.bin'); Te, _ = C.read_state(f'{K}/states/PAR8_W_203009.bin')
neq = 960; t = time.time()
T = C.advance(G, T0, F[:neq + 1], np.diff(et[:neq + 1]))
print(f'march {neq} steps: {time.time()-t:.1f} s; T0 surface range {T0[:,0].min():.2f}-{T0[:,0].max():.2f}')
d = T - Te
for j in (0, 1, 3, 6, 10, 15, 22):
    print(f'node {j:2d} z={z[j]*100:7.3f} cm: max|py-eng| {np.abs(d[:,j]).max():.2e}  mean {d[:,j].mean():+.2e}')
print('worst face', np.unravel_index(np.abs(d).argmax(), d.shape), np.abs(d).max())
