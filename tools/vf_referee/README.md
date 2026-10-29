# Offline view-factor referee (`tools/vf_referee/`)

An independent Monte Carlo estimate of the view-factor matrix, used to build
production view factors for large meshes and to referee the in-engine builder.

| File | Role |
|---|---|
| `mc_ref_pp.cpp` | Per-partner Monte Carlo sweep (Embree 4); one chunk of faces per process. |
| `ensure_mc_ref_pp.py` | Builds `mc_ref_pp` from the current source (cached by source and compile-command hash) so a stale binary is never used. |
| `run_pp_sweep.sh` | Runs a sweep through `ensure_mc_ref_pp.py`. |
| `assemble_vf.py` | Symmetrizes the sweep into a reciprocal CSR matrix, with a `.provenance.json` sidecar the GUI and engine read. |
| `assemble_vf_ooc.py` | Out-of-core twin of `assemble_vf.py`, byte-identical output in bounded memory. |
| `vf_stats.py` | Sparsity and row-sum statistics of a CSR. |
| `test_*.py` | Regression gates; `test_vf_assembler_gate.py` and `test_vf_cross_agreement.py` run in CI. |

The engine loads the result through `PTM_HMATRIX_VF_REPLACE`; see
`docs/radiation-exchange.md` for the method and `docs/ENVIRONMENT.md` for the
variables.

A CMake build can also produce `build/bin/mc_ref_pp`
(`-DPTM_BUILD_VF_REFEREE=ON`). The gates do not use it unless it is passed
explicitly (`--mc-ref-pp` or `$PTM_MC_REF_PP`).
