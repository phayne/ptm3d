# Configuration Guide

This document lists the YAML configuration keys PTM3D reads, with the defaults the
parser applies when a key is absent (`src/cpp/tools/config_parser.cpp`).

> If you change a key or a default in the parser, regenerate the Python layer's copy of
> this table as well. See [CONTRIBUTING.md](../CONTRIBUTING.md#changing-an-engine-configuration-key).

> Runtime behavior also depends on `PTM_*` environment variables.
> [ENVIRONMENT.md](ENVIRONMENT.md) has the complete table, including the
> `engine_env:` block with which a configuration declares its own environment.

## Structure

A configuration file has up to ten top-level sections: `simulation`, `mesh`,
`thermal`, `illumination`, `radiation`, `output`, `generic_body`,
`thermophysical`, `gpu` and `advanced`, plus the optional `engine_env` block.
Every section is optional, and absent keys take the defaults below.

```yaml
simulation:
  body: "Moon"
  start_time: "2025-01-01T00:00:00"
  duration_days: 29.53
  timesteps: 72

mesh:
  resolution: 5
  import_file: ""

thermal:
  run_thermal_mesh: true
  threads: 0
  equilibration_orbits: 4

illumination:
  include_primary_reflection: true
  use_batch_processing: true

radiation:
  enabled: true
  band_separated: true
  max_coupling_iterations: 16

output:
  directory: "output/simulation"
  format: "binary"
  vtk_mode: "both"
```

Runnable examples are in [`../examples/`](../examples/). Start with
`quickstart_moon_diurnal.yaml`.

## simulation

| Key | Type | Default | Description |
|---|---|---|---|
| `body` | string | `"Europa"` | Planetary body name |
| `start_time` | string | `"2025-01-01T00:00:00"` | Start time (ISO 8601) |
| `start_time_offset_seconds` | float | 0.0 | Ephemeris-time offset added after `start_time` is resolved (see [The output epoch grid](#the-output-epoch-grid)) |
| `duration_days` | float | 3.55 | Simulated duration of the output window |
| `timesteps` | int | 72 | Number of output timesteps |

### The output epoch grid

The output grid is half-open. A window of `timesteps` samples over
`duration_days` has cadence

```
dt = duration_days * 86400 / timesteps      [ephemeris seconds]
```

and samples `start_time + i*dt` for `i` in `[0, timesteps)`. The last sample
is therefore one cadence before the end of the window, and the next window
of the same simulation begins exactly at `start_time + timesteps*dt`.
Consecutive windows therefore tile, with no epoch sampled twice and none
skipped.

An equilibration prefix does not change this grid. `equilibration_orbits` adds
its own samples before `start_time`, so output frame 0 is at `start_time` exactly,
whatever the equilibration length or cadence.

`start_time_offset_seconds` shifts the whole grid by a number of ephemeris
seconds. Use it, rather than a shifted UTC `start_time` string, whenever a
window boundary must fall on an exact multiple of the cadence. The difference
TDB − UTC carries a periodic relativistic term (~1.7 ms amplitude, annual period)
and steps by a whole second at a leap second, so adding to the UTC string is not
the same shift. `scripts/run_chunked_simulation.py` uses it for every chunk.
See [chunked-run-equivalence.md](chunked-run-equivalence.md) for runs divided
into consecutive windows.

## mesh

| Key | Type | Default | Description |
|---|---|---|---|
| `resolution` | int | 5 | Icosphere subdivision level when no mesh is imported |
| `import_file` | string | `""` | Mesh file to load (PLY/OBJ/VTK) |
| `occluder_file` | string | `""` | Shadow-only occluder mesh, added to the ray-tracing scene but never solved or written. See [extended-occluder-architecture.md](extended-occluder-architecture.md) |
| `scale_factor` | float | 1.0 | Vertex scale factor (see [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md) before using it) |
| `export_basic` | bool | false | Write the generated mesh before solving |
| `use_topographic` | bool | false | Build the mesh from a topographic map |
| `topo_directory` | string | `""` | Directory holding the topography maps |
| `latitude_override` | float | -999.0 | Site latitude for topographic meshing (sentinel = unset) |
| `longitude_override` | float | -999.0 | Site longitude for topographic meshing (sentinel = unset) |

## thermal

| Key | Type | Default | Description |
|---|---|---|---|
| `run_thermal_mesh` | bool | false | Run the thermal solver (illumination only when false) |
| `threads` | int | 0 | OpenMP thread count (0 = auto) |
| `thermal_m` | int | 10 | Layers within the upper skin depth |
| `thermal_n` | int | 15 | Growth parameter $n$ of the fixed geometric depth grid ([thermal-model.md §4.2](thermal-model.md#42-geometric-depth-grid)). Every run uses the adaptive grid, so this value and the `--thermal-n` flag have no effect |
| `thermal_b` | int | 15 | Skin depths to the bottom layer |
| `equilibration_orbits` | float | 1.0 | Number of diurnal cycles of ephemeris-driven spin-up integrated before the output window. The number is the same on the CPU scalar, CPU SIMD and GPU persistent solvers. `equilibration_cycles` changes how it is reached and not the total. Under the legacy seed the CPU solvers also run an analytic pre-conditioner that this key does not count. See [Equilibration length](#equilibration-length) and [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md), item 2 |
| `equilibration_steps_per_day` | float | -1.0 | Equilibration steps per Earth day (86,400 s) on every body. On the Moon, 24 gives about 709 steps per lunation. -1 = unset: `equilibration_samples_per_diurnal_cycle` applies if set, otherwise the output cadence ([thermal-model.md §7.2](thermal-model.md#72-decoupled-equilibration-cadence)) |
| `equilibration_samples_per_diurnal_cycle` | float | -1.0 | Equilibration steps per diurnal cycle, on any body. Used only when `equilibration_steps_per_day` is unset. On a lunar mesh the largest surface-temperature error was about 7 K at 30, 4 K at 120 and 3 K at 240 ([thermal-model.md §7.2](thermal-model.md#72-decoupled-equilibration-cadence)) |
| `equilibration_cycles` | int | -1 | Replays of the equilibration forcing block used to reach `equilibration_orbits`, GPU persistent solver only (-1 = 1 replay). The block is `equilibration_orbits / equilibration_cycles` cycles long, so the integrated total does not change and only the flux RAM and illumination cost fall |
| `equilibration_target_K` | float | 0.1 | Accuracy targeted by the spin-up estimate, and the threshold that the retrospective certificate must meet. See [Retrospective equilibration certificate](#retrospective-equilibration-certificate) |
| `equilibration_stop_on_certificate` | bool | false | End the analytic pre-conditioner when the certificate reaches `equilibration_target_K` instead of at a fixed cycle count. CPU SIMD batch solver only. On the other solvers the run stops at startup with an error when this key is set |
| `eq_extrapolation` | map | off | Pass-stride extrapolation of the equilibration state (see below) |
| `compute_statistics` | bool | true | Compute and write per-face statistics |
| `use_gpu` | bool | false | OpenCL thermal solver (a CPU-only build stops at startup with an error) |
| `gpu_devices` | int | 1 | Devices in the thermal pool (0 = all enumerated) |
| `gpu_device_ids` | list[int] | [] | Explicit device indices. Overrides `gpu_devices` while non-empty. Must be a list. A scalar stops the run with an error ([explanation](#keys-the-engine-reads-as-lists)) |
| `gpu_batch_size` | int | 0 | Faces per GPU batch (0 = automatic cap) |
| `gpu_inter_batch_pause_ms` | int | 0 | Pause between consecutive batches on a device |
| `auto_calibrate_gpu` | bool | true | Probe the host once for a dispatch ceiling when no machine profile exists |
| `thermal_inertia_map` | string | `""` | Per-face `float32` thermal-inertia map |
| `h_parameter_map` | string | `""` | Per-face `float32` H-parameter map (mutually exclusive with the above) |
| `use_default_property_map` | bool | false | Auto-discover the bundled property map for the body |
| `property_maps_directory` | string | `"data/property_maps"` | Where to look for property maps |
| `geothermal_heat_flow` | float | -1.0 | Basal heat flow in W m<sup>-2</sup> (-1 = body default: 0.018 Moon, 0.030 Europa, 0.020 Mercury). Values > 0 are applied by every solver, and the log confirms with `[Thermo override] applied at solver: Qb ...`. A value of `0` selects the body default, not zero heat flow |
| `heat_flow_mask_file` | string | `""` | Face-index file for locally elevated basal flux |
| `heat_flow_mask_value` | float | 1.0 | Basal flux on the masked faces (W m<sup>-2</sup>) |
| `integration_days` | float | -1.0 | Statistics-integration window (-1 = the output window) |
| `integration_steps` | int | -1 | Same window, expressed in steps |
| `integration_offset_days` | float | 0.0 | Offset of the integration window |
| `integration_offset_steps` | int | 0 | Same offset, expressed in steps |
| `initial_state_file` | string | `""` | Load a per-face `T(z)` checkpoint and skip equilibration (see [Thermal-state checkpoints](#thermal-state-checkpoints)) |
| `final_state_file` | string | `""` | Write a per-face `T(z)` checkpoint at the end of the run |
| `checkpoint_lookahead_steps` | int | 0 | Extra output steps integrated past the export window, so `final_state_file` is the state at the next window's first epoch (see [Thermal-state checkpoints](#thermal-state-checkpoints)) |
| `state_restore_allow_mismatch` | list | `[]` | Per-field exemption from the checkpoint identity checks. Lists one or more of `depth_grid`, `mesh`, `config`, `body`, `grid_m`, `grid_n`, `grid_b`, `refgrid_H`, `epoch`. A mismatch in a listed field gives a warning instead of an error. No single setting exempts all fields |
| `checkpoint_every_output_steps` | int | 0 | Mid-run checkpoint cadence (0 = chunk boundaries only). Has no effect on the GPU solver. See [Batch-completion checkpointing](#batch-completion-checkpointing) |
| `batch_checkpoint` | bool | false | Save each completed face batch so that an interrupted run resumes from it instead of re-integrating the window. CPU SIMD batch solver only. The other solvers stop with an error when it is set. See [Batch-completion checkpointing](#batch-completion-checkpointing) |
| `batch_checkpoint_dir` | string | `""` | Directory for the partial files. Empty = `<output.directory>/_thermal_state/partial` |
| `batch_checkpoint_max_gib` | float | 64.0 | The run stops at startup if the projected set of partial files exceeds this size. 0 disables the check |
| `batch_checkpoint_flush_mib` | float | 0.0 | How much completed work to accumulate before writing a partial. 0 derives it as one sixteenth of the projected set, clamped to [1, 64] MiB |
| `eq_indirect_source` | string | `auto` | Source of the indirect flux applied during equilibration. `auto` estimates it from the first sweep of the coupling loop, `converged` from the converged sweep (more accurate, at about 1.9 times the run time), and `off` applies none ([radiation-exchange.md §7.7](radiation-exchange.md)) |
| `eq_indirect_passes` | int | 1 | Number of times that estimate is refreshed and the equilibration repeated. Must be ≥ 1. Prefer `eq_indirect_source: converged` to a second `auto` pass ([radiation-exchange.md §7.7](radiation-exchange.md)) |
| `eq_indirect_flux_file` | string | `""` | Explicit per-face `'QIND'` mean-indirect file. Overrides `eq_indirect_source` |
| `eq_indirect_subcycle` | string | `phase` | The indirect flux applied during equilibration when the output window is shorter than one diurnal cycle (chunk 0 of a chunked run, or an eclipse or flyby window). The mean over such a window is biased by the part of the cycle it covers (+61% on a 968-face crater). `phase` (default) takes the indirect flux of the coupling sweep over the last whole cycle of the equilibration span and applies it to every equilibration step by phase, so that the equilibration carries the time-varying indirect flux as a full-cycle run does. A constant mean spread evenly over the night biases the column warm (+4.8 K on the seed of that crater). `phase` requires that the equilibration span be swept against real temperatures (always on the CPU, and on the GPU under `gpu_rethermal: seed_series`), together with `seed_indirect: series` and `seed_window` auto or eq_span. Otherwise `whole_cycle` applies, with a warning. `whole_cycle` applies a constant whole-cycle mean estimated from the illumination tables and the radiosity operator. `legacy` keeps the window mean. If the equilibration span holds less than one whole cycle, the window mean is kept with a warning. Runs whose output window covers a cycle are unaffected. Recorded as `eq_indirect_subcycle` / `eq_indirect_subcycle_effect` in `provenance.txt` |
| `gpu_flux_forcing` | string | `interpolate` | Surface forcing within an output interval on the GPU solver. `interpolate` interpolates the absorbed flux linearly in time at each substep midpoint, as the CPU solvers do. `hold` holds the starting flux of the interval for the whole interval, which places the forcing half an interval early. On a 4990-face crater at 6 h cadence it gave T_max 1.27 K lower on average, and up to 21.9 K lower, than the CPU solver. It is kept only to reproduce archived GPU products, and it gives a warning. The interval after the last frame of a session interpolates toward the flux of the next epoch where the run has it (the equilibration phase into the output window), wraps to frame 0 on a replayed equilibration pass, and otherwise holds. That interval affects only the end-of-run state. The CPU solvers ignore the key. Also `PTM_GPU_FLUX_FORCING`. Recorded as `gpu_flux_forcing` in `provenance.txt` |
| `gpu_rethermal` | string | `seed_series` | How a GPU radiation-coupling re-thermal pass starts its columns. `seed_series` follows the CPU method (`PTM_CPU_RETHERMAL=seed-series`). Each pass re-seeds the selected faces from the composed flux of that pass (absorbed plus the current indirect series), re-equilibrates on it, then integrates the output window. The whole-mesh GPU solve then keeps the surface temperatures of the equilibration span, so the coupling sweep solves the enclosure there as well. This turns off `rad-skip-eq` and costs a radiosity sweep over the equilibration steps. `resume` restarts each pass from the post-equilibration column of iteration 0. On a 4990-face crater it gave T_min 1.9 K warmer than the CPU solver. It is kept only to reproduce archived GPU products, and it gives a warning. Warm-started chunks restore their checkpoint in either mode. Also `PTM_GPU_RETHERMAL`. Recorded as `gpu_rethermal` in `provenance.txt` |
| `indirect_mean_out_file` | string | `""` | Write the converged per-face output-window mean indirect flux as a `'QIND'` file |
| `compute_ice_stability` | bool | false | Compute the subsurface ice-stability depth |
| `enable_shadow_flux` | bool | false | Sub-grid crater shadowing correction |
| `crater_dd_ratio` | float | 8.0 | Diameter-to-depth ratio for that correction |
| `spinup_seed` | string | `physical` | The initial condition. `physical` seeds every face from its own cycle-mean conductive steady profile, so a permanently shadowed column starts on the geothermal floor rather than 42–47 K above it. `legacy` selects the original seed of each solver and the analytic pre-conditioner. Use it to reproduce archived results. An unrecognized name stops the run at startup |
| `seed_anchor` | string | `lumped` | The surface temperature from which the seed builds its conductive profile. Read only under `spinup_seed: physical`. `lumped` anchors at the cycle mean. `warm` uses the fourth-moment anchor, which by Jensen's inequality is never below the mean (+45 K on the lit faces of the example lunar 80N crater) |
| `seed_window` | string | `auto` | The forcing from which the seed is built. `auto` is `legacy` except in an eq-indirect re-solve whose output window is shorter than one whole cycle (chunk 0 of a chunked run, or an eclipse or flyby window). There the series seed has no effect and the legacy seed would fall back to the warm anchor with the plain profile, so `eq_span` is used instead for the rest of that run, including CPU re-thermals. This removes part of the warm bias of a chunked run relative to a single run, not all of it. The difference in T_min between `--chunks 2` and `--chunks 1` was -0.44 K (CPU, 968 faces) and +1.69 K (CPU, 140k-face crater). With the GPU thermal solver the difference was +2.45 K, and +1.58 K with `eq_indirect_subcycle: whole_cycle`. `eq_span` (opt-in) reads the forcing that the equilibration steps integrate (samples `[0, equilibration_steps)` of the composed absorbed flux, on the equilibration cadence, including the indirect flux those steps carry) and uses the `lumped` anchor on every solve. On a run whose output window covers a year it is less accurate than `legacy`, which then seeds on the exact annual mean. `legacy` keeps the original seed window, for reproducing archived results. `eq_span` and `auto` fall back to `legacy`, with a warning, when the run has less than ~0.9 cycles of equilibration. Also `--seed-window` and `PTM_SEED_WINDOW`. Recorded as `seed_window` / `seed_window_source` / `seed_window_effect` in `provenance.txt` |
| `seed_seasonal` | string | `off` | Two-timescale (seasonal) seed correction. At a polar site the deep column retains the seasonal cycle, which a seed or spin-up built from the run's own window cannot represent. Dark-season windows are then several K too cold on lit faces and up to ~25 K too cold on seasonally shadowed faces, and bright-season windows are too warm. With `on` (any body, with the 346.62-day eclipse year for the Moon and the orbital year for other bodies) or `auto` (the Moon only, and only on an output span under half an eclipse year that is not warm-started, see `seed_seasonal_window_days`), the run illuminates `seed_seasonal_phases` phases of one diurnal cycle each (`seed_seasonal_samples` samples) over the seasonal period ending at the output window, solves the periodic seasonal storage state of each face on its own column (harmonic balance in Kirchhoff G space), and moves the seed to it. Coupled runs add a per-phase indirect estimate from the radiosity operator in the eq-indirect re-solve. The correction is applied on the CPU SIMD and GPU solvers. The CPU scalar solver logs a notice and leaves the seed uncorrected. Warm-started chunks are not corrected. The accuracy depends on a partial cancellation of two biases. In a bright-season crater window about 0.26% of lit faces are up to ~1 K worse than with the default seed. In a dark-season window the correction is too warm on cold traps (median +0.28 K in T_max, with 6% of faces worse by more than 0.1 K, up to 1.3 K). Scope and evidence are in [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md) item 11. Also `--seed-seasonal` and `PTM_SEED_SEASONAL`. Recorded as `seed_seasonal*` in `provenance.txt` |
| `seed_seasonal_phases` | int | 12 | Seasonal phases K for `seed_seasonal` (4-64) |
| `seed_seasonal_samples` | int | 24 | Samples per diurnal cycle for each phase (8-480). The pre-pass illuminates K x samples extra epochs |
| `seed_seasonal_min_correction_K` | float | 0.5 | Per-face threshold for `seed_seasonal`. A face is corrected only where the correction moves its column by more than this value. The number of corrected faces is recorded in `provenance.txt` |
| `seed_seasonal_window_days` | float | 0.0 | The total output window, for the `auto` criterion (the Moon, output span under half an eclipse year, not warm-started). 0 = this run's own `duration_days`. A chunked run with `seed_seasonal` enabled should set it to the whole window, since each chunk's `duration_days` is only its slice |
| `seed_anchor_refuse_two_timescale` | bool | false | Reject the lumped anchor on a face that carries both diurnal and seasonal power, where one conductance cannot represent both. Off by default because a rejected face receives `warm`, the anchor with the +45 K bias. Detection and reporting are always on |
| `analytic_seed_cycles` | float | -1.0 | Diurnal cycles of the CPU analytic clear-sky pre-conditioner: `cos(lat)·cos(h)` insolation with no topography, no shadowing, no IR and no indirect flux. It is a forcing, not a seed. It is not counted by `equilibration_orbits`, and the GPU solver has no equivalent. -1 follows `spinup_seed` (physical gives 0, legacy gives Moon 24, Europa 2, Mercury 4). 0 turns it off. A positive value is used as given and reported in the log |
| `lateral_conduction` | bool | false | Lateral heat conduction between neighboring face columns, on every capacitive node of the depth column. CPU SIMD batch solver only (`use_gpu: false`). See [Criteria for lateral conduction](thermal-model.md#21-criteria-for-lateral-conduction) |
| `lateral_scheme` | string | `monolithic` | `monolithic` solves the vertical Crank-Nicolson rows and a backward-Euler lateral term as one coupled system per substep (the validated scheme). `balanced`, an operator split, is kept only for comparison |
| `lateral_picard` | int | 0 | Extra solves per substep with the material properties re-evaluated at the mid-step temperature (0-4, `monolithic` only). 2 keeps energy conservation errors below 0.1% under large lateral contrasts |
| `lateral_activation_c` | float | 2.0 | Couple an edge only where `min(l_i, l_j) < c max(D*_i, D*_j)`, with `l` the mean facet edge and `D*` evaluated at each face's starting temperature. A value ≤ 0 couples every edge |
| `lateral_cg_rel_tol` | float | 1e-8 | Relative residual of each lateral linear solve. An unconverged solve stops the run |
| `lateral_zero_coupling` | bool | false | Run the lateral solver with every edge decoupled, so that no heat crosses any edge. Both runs of a comparison with and without lateral conduction then use the same solver, and their difference is the lateral term alone. For comparisons only |
| `rethermal_all_faces` | bool | false | During radiation coupling, re-integrate every face instead of the subset whose forcing changed. Implied by `lateral_conduction`. In a comparison with and without lateral conduction, set it on the run without it so that both runs re-integrate the same faces |
| `lateral_advisor` | bool | true | Print the pre-run `[lateral-advisor]` block, which compares the facet sizes with `D*` and reports whether `lateral_conduction` is recommended. It affects the log only, and `false` turns it off |

The checkpoint files are used by the chunked-run workflow
(`scripts/run_chunked_simulation.py` and `ptm_merge_chunks`).


### Thermal-state checkpoints

`thermal.initial_state_file` and `thermal.final_state_file` pass a whole-mesh
per-face `T(z)` column from one run to the next. Every whole-mesh solver (the GPU
persistent solver, the CPU SIMD batch solver, the CPU scalar solver, and the
streamed CPU solver) reads and writes them in the same format
(`uint32` magic `0x50544D54` `'PTMT'`, `uint32 n_faces`, `uint32 n_layers`, then
`n_layers × n_faces` `float32` in layer-major order). A checkpoint written by
either backend can be restored on the other.

A checkpoint has three properties.

- A checkpoint that cannot be read or written stops the run. A missing file, a
  wrong magic number, a truncated body, a face- or layer-count mismatch, or a
  destination that cannot be opened for writing ends the run with a message
  naming the reason. The run never falls back to a cold start, because the
  output of a cold-started chunk would contain the spin-up transient.
- The file is `float32` on every solver. The CPU solvers carry `T` in
  `float64`, so writing a checkpoint rounds each value once, at the write.
  Restoring is exact. A CPU-to-file-to-CPU round trip therefore reproduces the
  file bit for bit and the in-memory column to `float32` precision.
- The write is atomic. The model writes a temporary file beside the
  destination, `fsync`s it and renames it into place, so the destination is
  always either the previous complete checkpoint or the new one.


#### The v2 identity trailer

Checkpoints are written in the PTMT v2 format. The v1 header and body are
unchanged, and an identity trailer is appended after them. Files in the v1
format still load.

The trailer records what the column was integrated with:

| Field | What it is for |
|---|---|
| `z[n_layers]` | the depth grid, so a restore onto a different grid with the same layer count is caught |
| `et_seconds` | the SPICE ephemeris time the column is valid at |
| `output_cadence_s` | the writing run's output cadence |
| `refgrid_H`, `grid_m/n/b` | the grid parameters, named individually so that a mismatch report identifies the parameter |
| `mesh_geomhash` | a face-area hash of the mesh |
| `config_identity_hash` | the resolved thermophysical scalars and grid parameters (see below for what it does not cover) |
| `body_name`, `engine_git_sha` | provenance |
| `flags` | pre- or post-radiation-coupling, whether every face was patched, and whether the writing run had `coupling_select_dt > 0` |
| `crc32` | a checksum over the whole file, so a length-correct but zero-filled cloud-sync placeholder cannot pass as a checkpoint |

`ptm_state_info <file>` prints all of it, verifies the checksum, and, given
`--body`/`--m`/`--n`/`--b`, compares the depth grid recorded in the file
with the grid those settings would build.

The physics hash does not cover everything. Per-face maps (albedo, thermal
inertia, `H`) and a per-face `Qb` mask change the physics face by face and are
outside the hash, as are the view factors. Two runs that differ only in a loaded
map can produce the same hash. `src/cpp/thermal/state_identity.h` lists what is
hashed, and the log states this limit whenever a checkpoint is used.

`W(0)` is reserved and always written as NaN. It is intended as an input to a
warm-start admissibility criterion that is not implemented. Computing it requires
a reduction over the surface time series of the output window, which is available
at the time of the write on some solvers and not on others.
`ptm_state_info` reports it as "not recorded".


#### Identity checks at restore

A restore compares the trailer with the current run. Each check is
independent, and a failed check stops the run with a message that names it.

| Clause | Behavior |
|---|---|
| `depth_grid` | the recorded `z[]` differs elementwise (relative 1e-12): error, naming the first node |
| `mesh` | `mesh_geomhash` differs: error |
| `config` | `config_identity_hash` differs: error |
| `body` | a different body: error |
| `grid_m` / `grid_n` / `grid_b` / `refgrid_H` | error, naming the parameter |
| `epoch` | see below |
| batch-completion partial | error, and this check cannot be exempted |

`thermal.state_restore_allow_mismatch` takes a list of those field names and
turns each named check into a warning. The checks cannot all be turned off with
a single setting, so each exemption is an explicit choice.

The epoch check distinguishes a hand-off between consecutive windows from a
warm start at an unrelated time.

- The checkpoint epoch equals the first integrated epoch of this run. The
  states tile the timeline, as in the chunked workflow, and the run continues
  without a message.
- The checkpoint is behind by up to one output cadence. This is a hand-off
  written without `checkpoint_lookahead_steps: 1`, so the restored column is one
  step old. The run prints a warning with both epochs and the cadence, and
  continues.
- The checkpoint is behind by more than one cadence. The stored column is not
  the state at the epoch of this run, and the run stops with an error.
- The checkpoint is ahead of the first epoch of this run. The two windows
  overlap, and the run would re-integrate a span the column has already
  passed. The run stops with a separate error message.

A warm start from a different season changes the results. In a test at one
lunar south-polar site it moved `temperature_min` by 5.75 K rms and 49.9 K
maximum, and shifted the ice-stability depth product by 19.3 percentage points.
Chaining is admissible only where the states tile the timeline, which is what
`scripts/run_chunked_simulation.py` does within a window. It is not admissible
across windows or seasons.

A v1 file still loads. It carries no trailer, so the restore prints one line
for each quantity it could not verify. Write the state again in the v2 format
to enable the checks.

Setting `initial_state_file` sets `equilibration_orbits` to 0 and turns off the
analytic clear-sky pre-conditioner, because a restored column is a continuation
and a new spin-up would discard it. Both changes are reported in the log.

The solver integrates
`total_time_steps − 1` intervals, so by default `final_state_file` holds the
state at the last exported frame. On a half-open grid the next window begins
one cadence later, so a hand-off needs one extra grid point.
`thermal.checkpoint_lookahead_steps: 1` supplies it. The run integrates one
output step past the export window and writes the checkpoint there. That extra
step is neither written as a frame nor included in `statistics/`. It is also
excluded from the forcing mean of the equilibration seed, so it does not change
the initial condition of the run. `scripts/run_chunked_simulation.py` sets it
for every chunk, and it is not needed outside the chunked workflow.

One case gives a notice in the log and not an error. With radiation coupling on
a CPU solver, the coupled re-run keeps no `T(z)` column, so the checkpoint holds
the uncoupled column on the faces that were re-run.


### Batch-completion checkpointing

`thermal.batch_checkpoint` limits the work lost when a run is interrupted. It is
off by default, and with it off the output is unchanged.

It is not a checkpoint every N steps. The thermal solvers integrate one face
batch through all timesteps before the next batch starts. At no time during the
run do all faces share an output step, and a whole-mesh state written mid-run
would include faces that have not been solved. For this reason
`thermal.checkpoint_every_output_steps` gives a warning and has no effect.

At a batch boundary, however, a well-defined set of faces is finished for the
entire window. The run saves these as a self-describing partial file (`PTMB`)
per flush, which holds the surface time series of each completed batch, its
end-of-window `T(z)` column and its slice of the depth-map accumulator. On
restart the run validates the partial files it finds, restores those batches and
integrates only the rest.

A partial file holds the surface time series, not only the column, because
every statistic and every per-step product is a reduction of it. Per face that
is `n_timesteps × 8 B` plus `n_layers × 4 B`. On a 2,397,890-face mesh with 1440
output steps the partial set is about 26 GiB per window. The log gives the exact
projection at startup, and the run stops if it exceeds
`batch_checkpoint_max_gib`.

Two losses remain. A run that stops inside a batch loses the whole time
integration of that batch. A resumed run also repeats the illumination and
radiation setup, because only the thermal integration is saved. To limit the
loss further, divide the run into chunks with
`scripts/run_chunked_simulation.py`, which gives the same results as a single
run.

Batch-completion checkpointing is implemented on the CPU SIMD batch whole-mesh
solver only. With the key set, the GPU persistent solver, the CPU scalar solver
and any radiation-coupled run stop at startup with an error:

| Combination | Reason for the error |
|---|---|
| `thermal.use_gpu: true` | Not implemented on the GPU solver. Resume equivalence is tested by killing a run abruptly, which cannot be done safely with OpenCL queues in flight |
| `--no-simd-batch` | The depth grid of the CPU scalar solver is per face, so one partial identity cannot describe the set |
| `radiation.enabled: true` | The whole-mesh solve is entered more than once, and a partial set that survived from the first pass would be restored in the second pass in place of its results |

Every partial file carries the depth grid, epoch, output cadence, mesh geometry
hash, configuration identity hash, body and engine SHA of the run, and a CRC32
over the whole file. On resume, any mismatch stops the run with an error that
names the field. The file also carries an input digest, a checksum of the bytes
that the batch integrated (its absorbed-flux rows, its seeded column, its
per-face `Q_b`). The digest detects differences that the configuration hash does
not. A different albedo map, different view factors, or radiation on instead of
off all leave `config_identity_hash` unchanged and all change the digest.

Completeness is computed from the partial files that validated, and not from
`manifest.json`. The manifest is an index. If it names a file that is missing or
does not validate, the resume stops with an error. If the manifest is behind the
files on disk, because a flush completed and the manifest rewrite did not, the
resume reports this and continues. A partial file cannot be loaded as
`thermal.initial_state_file`, because its format and magic number are different.

The crash-resume test is `scripts/batch_checkpoint_crash_gate.py`. It kills the
run at several points, resumes, and requires the finished output files to be
identical to those of an uninterrupted run.


### Equilibration length

`thermal.equilibration_orbits` is the number of diurnal cycles of
ephemeris-driven forcing that the solver integrates before the output window
opens. This holds on all three solvers and at any `equilibration_cycles`:

```
block_cycles = equilibration_orbits / equilibration_cycles
integrated   = equilibration_cycles * block_cycles = equilibration_orbits
```

`equilibration_cycles` controls memory and cost, not the length of the spin-up.
Raising it shortens the assembled forcing block and replays it more times, which
reduces the equilibration flux array and the illumination calculation by that
factor while integrating the same number of cycles. A long equilibration may
not fit in memory otherwise: a 25-cycle equilibration block is about 6.1 GB on a
509,337-face lunar mesh. There is a trade-off, because *N* replays of a
short block are not the same forcing as one long block of real ephemeris. The
block should be a whole number of diurnal cycles, so that each replay restarts
at the point in the cycle where the previous one ended. The log gives a warning
when it is not.

Only the GPU persistent solver replays the block. A CPU run prints a notice and
integrates the same total length as one block of the full size, which uses more
memory and computes more illumination but is the same spin-up.

Every run prints the configured and the integrated number before it starts, so
that they can be compared directly:

```
Equilibration (configured): 20.00 day/night cycles  (thermal.equilibration_orbits)
Equilibration (integrated): 20.00 day/night cycles = 2 pass(es) x 10.00-cycle forcing block
```

`equilibration_orbits` does not count the pre-conditioner. Under `spinup_seed: legacy` the CPU scalar and
SIMD solvers first run an analytic clear-sky pre-conditioner (`NYEARSEQ * round(year/day)` cycles:
24 on the Moon, 2 on Europa, 4 on Mercury) with flat, unshadowed insolation. Under the default
`spinup_seed: physical` it is off (`analytic_seed_cycles` resolves to 0) unless requested. The
banner prints its length in both cases.

### Retrospective equilibration certificate

`equilibration_orbits` sets the length of the spin-up. It does not show whether
that length was sufficient. The certificate tests this with a conservation law,
and not with the rate at which the temperature is still changing.

The quantity is the cycle-mean interface heat flux at a cut depth
`z_c = min(5 delta_ref, z_max/2)`, compared with the geothermal flux that drives
the column, and converted to kelvin:

```
r_c   = qbar_c - Qb                       [W/m^2]
E     = |r_c| / h,   h = 4 eps sigma T^3
```

Below the forced layer an equilibrated column carries exactly `Qb` at every
depth, because nothing there stores or releases net heat over a cycle. `r_c`
therefore depends on the spatial gradient of the state and not on its rate of
change. The distinction is significant in single precision. When float32
rounding freezes a column, it stops changing without having converged, and a
test based on the change between passes then reports convergence. In an
example lunar state that did not change over 100 further lunations, the
certificate reported 0.26 K against a 0.1 K target.

The choice of `h` depends on the product:

| Observable | `h` |
|---|---|
| diurnal maximum | `4 eps sigma T_max^3` |
| diurnal minimum | `4 eps sigma T_min^3` |
| diurnal mean | the harmonic cycle mean of `4 eps sigma T0(t)^3` |
| deep products | `E_deep = |r_c| tau1 / C_below` (unbiased, not conservative: pointwise 0.91–1.09; reported, never used to pass) |

The mean uses the harmonic cycle mean and not the arithmetic one. The surface
displacement caused by a flux residual is largest where the radiative
stiffness is smallest, so the two are anti-correlated and the arithmetic mean
understates the error (by a factor of 8.7 for a 25 %/75 % duty cycle between
350 K and 100 K). The minimum is in the table for the same reason. It is the
instant of the cycle with the smallest radiative stiffness and therefore the
surface observable most sensitive to a flux residual. A lunar spin-up stopped on
the mean alone ended 13.6 cycles early and left `temperature_min` 0.69 K from
its value after the full spin-up. Stopped on all three observables, it was within
0.06 K.

After a run, `build/bin/ptm_eq_certify` evaluates the certificate from a
saved `PTMT` checkpoint and the statistics beside it. It needs no mesh, no view
factors and no GPU, so it can be applied to a finished run, or to the last chunk
of a running one, while another job uses the GPU:

```
build/bin/ptm_eq_certify --config <run>/chunk_000/config.yaml \
    --state <run>/_thermal_state/thermal_state_after_chunk_000.bin \
    --stats <run>/chunk_000/statistics --period-ref annual
```

Exit status 0 = every scored face passed, 2 = at least one FAIL, 3 = at least
one INDETERMINATE, 1 = the tool could not run. The tool does not treat a
snapshot of a forced column as a cycle mean, and it does not substitute the
arithmetic mean of `h` for the harmonic one. Where the data cannot decide, it
reports INDETERMINATE and lists the input that would decide it. `Qb` is read
from the bottom boundary condition of the state itself and not from the config,
so the certificate uses the heat flow with which the column was integrated.

During a run, the certificate is used through `thermal.equilibration_stop_on_certificate` (default false).
When true, the analytic clear-sky pre-conditioner ends when all three surface
observables reach `equilibration_target_K` instead of at its fixed cycle count,
and the achieved `E` is reported in both cases. The key controls only that
phase. The ephemeris-driven equilibration prefix cannot end early, because its
steps are epochs on the output grid and stopping it would move the output
window. The key is implemented on the CPU SIMD batch solver. On every other
solver, and in a run with no pre-conditioner enabled, the run stops at startup
with an error when the key is set.

### Equilibration extrapolation (`thermal.eq_extrapolation`)

The GPU equilibration phase replays one forcing block several times
(`equilibration_cycles`). Because the same block is replayed unchanged, the
temperature field after each pass is a fixed-point iteration, and its limit can
be estimated from a few passes. This option applies reduced-rank extrapolation
to that sequence and then runs the final pass from the extrapolated state. The
final pass refines the extrapolated state, checks it, and provides the state
from which the output window starts.

The default is `off`, which leaves the equilibration unchanged.

```yaml
thermal:
  equilibration_orbits: 50   # integrated diurnal cycles
  equilibration_cycles: 5    # ...as 5 replays of a 10-cycle block
  eq_extrapolation:
    method: rre        # off | rre          (default off)
    order_p: 4         # 3 or 4; needs equilibration_cycles >= order_p + 1
    tol_K: 0.5         # certificate reporting threshold, K
    hold_K: 1.0e-3     # per-pass increment below which a layer is held, K
```

An `order_p` below 3 stops the run with an error, because a second-order fit is
less accurate than no extrapolation. If `equilibration_cycles` is smaller
than `order_p + 1` the run gives a warning and uses plain cycling.

Every run reports a one-line summary of how many faces were extrapolated, how
many the certificate declined, the 95th-percentile error estimate, the fraction
of layers held, and the number of near-degenerate windows.

On a 509,337-face lunar polar mesh at
`order_p: 4` with five passes, the extrapolated state was 55 times closer (sunlit
faces) and 64 times closer (seasonal-terminator faces) to a 200-lunation reference
than four passes of plain cycling, for 1.18 times the run time of four passes. The
benefit is confined to faces that receive flux.

About 70 percent of that mesh receives no flux during equilibration, and those
faces gain nothing from extrapolation, although they have not converged.
On the GPU solver their deep columns stop changing between passes while still
about 0.36 K short of the steady state, and a sequence that has stopped
changing cannot be extrapolated. A column that is unchanged between passes
is not necessarily converged. See
[KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md) sections 1 and 3.

The method has one known limitation. The fit assigns a weight to the first snapshot in its
window, which is the initial guess. For a depth layer that converged early while
shallower layers were still evolving, that weight can pull the layer back toward
the initial guess. The displacement is bounded by

    displacement  <=  |gamma_0| * (initial-guess error at that depth)

with `gamma_0` the weight on the first snapshot, about 0.09. The displacement is
therefore roughly ten percent of the initial-guess error and not an
amplification of it. The `hold_K` test suppresses it by holding any layer whose
last-pass increment is below the threshold, and the run reports how many windows
remain susceptible. On the mesh above, 73 of 509,337 faces retained a
displacement, the largest 0.47 K at 1 m depth, all on faces receiving no flux
during equilibration. Raising `hold_K` holds more layers and exchanges
extrapolation speed for reliability. The layers it holds converge at the
ordinary rate.

## illumination

| Key | Type | Default | Description |
|---|---|---|---|
| `include_primary_reflection` | bool | true | Include reflected sunlight and thermal IR from the primary body |
| `use_batch_processing` | bool | true | Trace illumination rays in vectorized batches (faster, same results) |
| `primary_disk_body_cull` | bool | true | On open patches, test primary-disk (Earthshine, planetshine) probe and sample rays against the body sphere R_occ, as the direct Sun is tested. With `false`, a face tilted toward a primary below the limb receives its IR and reflected light through the body. That setting reproduces archived runs and prints a warning. Recorded in `provenance.txt` (`primary_disk_body_cull`, `primary_disk_body_cull_effect`) |
| `cadence_advisor` | bool | true | Print the `[cadence-advisor]` block: the output step against the recommended 0.4 θ/ω, the 5% peak-flux limit, the largest temperature rise a missed flash can cause, and after the run the estimated upper bound on the overstatement of the cold-trap area. It affects the log and `provenance.txt` only, and `false` turns it off. See [Choosing the output cadence](thermal-model.md#22-choosing-the-output-cadence) |
| `use_albedo_map` | bool | true | Load a spatially varying albedo map when one is available |
| `albedo_map` | string | `""` | Explicit albedo FITS map |
| `albedo_maps_directory` | string | `"data/albedo_maps"` | Where to look for albedo maps |
| `positions_file` | string | `""` | Precomputed body/Sun positions instead of SPICE lookups |
| `use_gpu` | bool | false | GPU BVH shadow-ray tracer instead of CPU Embree. Validated against CPU Embree on an AMD Radeon Pro W6800X for lunar meshes. On other GPU hardware run `ptm gpu-selfcheck` first ([details](KNOWN_LIMITATIONS.md#physics-and-numerics)). The run prints a NOTICE to this effect |
| `gpu_devices` | int | 1 | Devices for the ray tracer (0 = all enumerated) |
| `gpu_device_ids` | list[int] | [] | Explicit device indices. Must be a list, and a scalar stops the run with an error ([explanation](#keys-the-engine-reads-as-lists)) |
| `gpu_inter_batch_pause_ms` | int | 0 | Pause between consecutive batches on a device |
| `bvh_cache_file` | string | `""` | GPU BVH cache path (empty = auto under the output directory, and `off` disables). Reused only if built from exactly this mesh and build settings (content key). Otherwise it is rebuilt, and the log gives the reason |
| `solar_disk_scheme` | string | `"legacy"` | Treatment of the solar disk for direct sunlight. `legacy`: center ray clear -> full disk, else 9 samples. `chord`: a 25-probe screen, then a 3-column edge bisection integrated exactly and checked with 128 rays, 1024-sample fallback, body sphere per ray. `reference`: validation only (~12,900 rays per partly lit face-frame). CPU batch solver only. With a non-legacy scheme the direct solar flux is computed on the CPU instead of by the GPU tracer, with a warning. See [solar-disk-scheme.md](solar-disk-scheme.md) |
| `solar_limb_darkening` | string | `"none"` | `none` (uniform solar disk) or `neckel_labs` (bolometric Neckel & Labs 1994 law, F/I = 0.829). The total flux of the disk is unchanged, while partly visible disks and slivers change. Requires a non-legacy `solar_disk_scheme` (an error with `legacy`) |
| `solar_disk_chord_check_adaptive_n` | int | `0` | Faster chord integration that changes the results, off by default. N > 0: when every screen probe agrees with the fitted edge, check the fit with N points instead of 128. `chord` only (an error otherwise). See [solar-disk-scheme.md](solar-disk-scheme.md) |
| `solar_disk_chord_temporal` | float | `0.0` | Faster chord integration that changes the results, off by default. Brackets each column's edge search by +/-delta disk radii around the edge extrapolated from the face's two previous frames (consecutive equal Sun steps only). Results depend on the frame sequence. `chord` only (an error otherwise) |
| `solar_disk_chord_fallback_n` | int | `1024` | Chord fallback Hammersley samples (8-4096). The results change when it is not 1024. `chord` only (an error otherwise) |

## radiation

Radiation exchange is configured from YAML as well as from the command line, and
the command line takes precedence. The full key list, with the physics behind each, is in
[radiation-exchange.md §9](radiation-exchange.md). The most commonly used keys are:

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | bool | false | Enable inter-face radiation exchange |
| `band_separated` | bool | true | Band-separated visible + thermal-IR radiosity (false = single gray band) |
| `exchange_selection` | string | `all` | Which faces keep the computed indirect flux. `all` (default): every face keeps it, and the result does not depend on how the run was split into chunks. `shadow`: only faces shadowed at least `min_shadow_fraction` of the time keep it. This setting reproduces earlier published runs and is best paired with `coupling_mode: legacy` |
| `min_shadow_fraction` | float | 0.3 | The shadow-fraction cut used by `exchange_selection: shadow`. Read only under that setting, so it has no effect in a default configuration. The shadow fraction is computed over the illumination series of the run itself, so a warm-started chunk (which has no equilibration samples) reclassifies faces and sets their indirect flux to zero |
| `view_factor_file` | string | `""` | Save/load path for the view-factor matrix. A relative path is resolved against the directory of the configuration file |
| `use_hmatrix` | bool | false | Hierarchical-matrix compression for large meshes |
| `use_monte_carlo_vf` | bool | false | Monte Carlo ray-gather view-factor builder |
| `mc_rays_per_face` | int | 100000 | Rays cast per face by that builder. Read only when it is selected |
| `nearfield_quadrature_order` | int | 3 | Order of the symmetric double-area quadrature for near-field pairs. 3 is the validated value. Order 1 is the centroid estimator with no near-field treatment. It is +51.7 % in error on the closed form for a touching pair (order 3 is +18.5 %), and the log identifies it at startup as a diagnostic setting |
| `coupling_mode` | string | `bounded` | Convergence criterion of the thermal–radiation coupling loop. `bounded` (default): stop once the largest per-face temperature change falls below `coupling_temp_threshold`, which bounds the residual on every face. `coupling_flux_threshold` is not part of this test. `fixed`: no early exit. Every window runs exactly `max_coupling_iterations` sweeps, so the result cannot depend on how the window was chunked. `legacy`: stop once only `coupling_face_fraction` of faces are still changing. This mode has no error bound, and the result of a coupled run then depends on the chunking. It is the only mode that reads `coupling_face_fraction`. The value is case-sensitive: `Bounded` is a startup error |
| `max_coupling_iterations` | int | 16 | Maximum number of Picard iterations (see radiation-exchange.md §7.5) |
| `coupling_temp_threshold` | float | 0.5 | Temperature convergence threshold (K), used under `coupling_mode: bounded` |
| `coupling_flux_threshold` | float | 0.1 | Indirect-flux cut used to select faces for another thermal solve, in W/m². It is not a convergence criterion, and it is read only under the non-default `coupling_select_mode: legacy`. The default selection ignores it |
| `beaming` | bool | false | Smith (1967) anisotropic beamed thermal emission. Used for lunar polar site and regional configurations, and off on sub-cm meshes (see below) |
| `rms_slope_deg` | float | 0.0 | RMS sub-facet slope $w$ in degrees, the parameter of the beaming model. No effect unless `beaming: true` |
| `use_gpu` | bool | false | OpenCL radiosity solver |
| `gpu_batch_size` | int or `auto` | `auto` (0) | Timesteps per multi-GPU radiosity pool dispatch (`radiation.gpu_devices` > 1). `auto` sizes it from free host RAM at the start of each sweep ([below](#radiosity-pool-batch-size-gpu_batch_size)). An integer overrides it. Results are identical at any value |
| `gpu_inter_batch_pause_ms` | int | 0 | Idle gap between pool batches on a device. The host work between batches is done during the gap |
| `verbose` | bool | false | Detailed radiation progress |

### Radiosity pool batch size: `gpu_batch_size`

The multi-GPU radiosity pool solves the timesteps of a sweep in batches. Every timestep is an
independent solve, so the batch size changes run time and host memory but not the results. The
host buffers of a batch take 72 bytes per face per timestep (about 0.13 GB per timestep on a
1.8 M-face mesh). GPU memory does not grow with the batch.

`auto` (the default, also `0` or an absent key) chooses the batch at the start of each sweep:

1. The budget is 0.25 x (available host RAM - reserve), where available is free plus inactive
   memory on macOS (`freeram` on Linux) and reserve = max(4 GiB, 10% of physical RAM). Memory the
   run already holds when the sweep starts (view factors, flux stores, thermal state) is not
   counted as available. The factor 0.25 leaves room for the batch buffers that stay allocated
   into the re-thermal phase, other jobs on the host, and memory that is inactive but slow to
   reclaim.
2. The batch is at most the sweep. A sweep that fits whole is split into two even batches when
   each half is at least 150 timesteps, because with a single batch the host post-processing runs
   after the GPUs finish instead of beside them. On a 1.8 M-face mesh with a 3000 ms pause, the
   pause added no run time from 156 timesteps per batch, and larger batches gained nothing more.
3. The sweep takes the fewest batches the budget allows, evened over the sweep and rounded up to
   whole steps per device when that keeps the batch count and stays within the budget. Otherwise
   the batch stays uneven (an extra batch would cost every device a full pause) and the device
   that takes the remainder steps rotates from batch to batch. The floor is min(sweep, 4 x devices). If the budget affords less than the floor,
   the floor is used and a warning gives the memory it needs.

The log prints the chosen batch and every input (available and total RAM, reserve, budget,
bytes per timestep, timesteps affordable, floor), and `provenance.txt` records it as
`rad_pool_batch=`. `PTM_RAD_POOL_BATCH_BUDGET_MB` replaces the RAM-derived budget with a fixed
one. If the host memory query fails, a fixed 4096 MiB budget is used, with a warning.

An explicit integer takes precedence (clamped to the sweep). A value below the floor is raised
to it, with a warning. The CLI equivalent is `--rad-gpu-batch-size N|auto`.

A larger batch means longer GPU work between rests, while the inter-batch pause still
applies between batches. Measured dispatch budgets, batch sizes and pauses are given in
[gpu-wedge-guardrails.md](dev/gpu-wedge-guardrails.md#g7-dispatch-budget-batch-size-and-pauses-measured-2026-11).

### Beamed emission: `beaming` + `rms_slope_deg`

Rough regolith beams thermal IR toward the Sun, because the Sun-facing micro-facets are the hot
ones. In a polar crater this redistributes emission from the sunlit wall onto the shadowed floor
and warms the floor by several kelvin. The physics, the validation and the full cost table are in
[radiation-exchange.md §8](radiation-exchange.md). This section describes the configuration keys.

Both keys are required:

```yaml
radiation:
  beaming: true            # THE SWITCH. Default false; nothing beams without it.
  rms_slope_deg: 32.5      # THE KNOB. RMS sub-facet slope w, degrees.
```

- `beaming: true` with `rms_slope_deg` unset or $\le 0$ is a startup error.
- `beaming: false` with a nonzero `rms_slope_deg` runs Lambertian emission and prints one line
  stating that the slope has no effect. The value may be kept in the file as documentation.
- A config with no `beaming` key and `rms_slope_deg > 0` still beams, and prints a
  `[config] NOTE` asking for the explicit key. New configs should set both keys.

The slope $w$ is a surface property, not a fit parameter. Published lunar values span
20–37°: 32.5° (Jhoti et al. 2024, Table 1, lunar highlands, the validated polar value and the
recommended default), ~20° (Bandfield et al. 2015, Diviner global), ~37° (Rubanenko et al. 2020).

Whether beaming applies depends on the facet scale. Beaming is a sub-facet model.
It represents roughness that the mesh does not resolve, at the scale where surface elements become
thermally isolated from one another, which is millimeters to centimeters on lunar regolith.

Beaming should therefore be included only when the facet scale is much larger than ~cm. At or below
that scale the mesh already resolves the anisothermality, and beaming would count it twice. The
facet scales, $\sqrt{\text{mean face area}}$, of several meshes are:

| Mesh | $\sqrt{\bar{A}}$ | beaming |
|---|---|---|
| `sunlit_radcore_r50` (50 km core, 140k faces) | 237 m | on |
| `haworth_radcore_r32c` (32 km core, 100k faces) | 181 m | on |
| a 3 km lunar polar site core (2.4 M faces) | 3.46 m | on |
| Moon3D testbed scan (509,337 faces) | 1.96 mm | off |

Beaming is used for site and regional lunar polar configurations (meters and
up). It is not enabled in the year-long configurations of the Moon3D paper (sub-cm facets), in `examples/validation/baseline_crater80N_*.yaml` (the
regression baselines) or in the `tutorial04_crater80N_*` tutorial configs. It
should not be carried to Europa or asteroid configs, because $w = 32.5°$ is a value for lunar regolith
validated on lunar polar-crater geometry.

The cost is about +22% run time and 5–6.5 GiB of host RAM on a 140k-face lunar polar mesh.
The memory term scales with the number of nonzeros in the view-factor matrix and reaches an
estimated 22–28 GiB on a 2.4 M-face polar site, so check the memory before enabling it on a large
mesh. Every run with beaming on prints a `[beam-cost]` breakdown.

Two caveats apply. With beaming, the results depend on how a run is divided into chunks,
because the emitter and receiver sets are defined per window (the run gives a warning).
Beaming also does not remove the known winter and margin warm biases.

## flux_cubes

A passive probe that samples the flux field at fixed points in space. The block
is at the top level and not nested under `output:`, so that it is not affected by
the body defaults that override `thermophysical:` on a built-in body. The
physics, the geometry cache, the ray budget and the outputs are described in
[flux-cubes.md](flux-cubes.md). This table lists the keys.

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | bool | false | Run the flux-cube pass. Turning it on does not change any mesh field |
| `cache_file` | string | `""` | Geometry cache built by `ptm_flux_cube_geom`. Required when enabled, unless `cache_files` is given |
| `cache_files` | list[string] | — | Multi-height: one cache per cube height, all consumed by one run. Setting both this key and `cache_file` is an error |
| `locations` | list[map] | — | Named locations (`name`, `lat_deg`, `lon_deg`, optional `height_m`, `edge_m`). The parser's only sequence of maps, so its per-entry fields have no dotted keys |
| `write_binary` | bool | true | Write `flux_cubes/flux_cubes.fcb` |
| `write_csv` | bool | true | Write `flux_cubes/csv/cube_XXXXX.csv` |
| `write_csv_max_cubes` | int | 0 | Maximum number of per-cube CSV files on a long window (0 = every cube). A named location is always written |
| `shadow_probe_rays` | int | 5 | 1 = cube center only, 5 = center plus four rim points |

The geometry keys (`spacing_m`, `height_factor`, `n_rays`, …) are accepted and
reported as having no effect. Their values are stored in the cache header, and the
run uses those. A warning names any key the parser does not read
(`warnUnreadFluxCubeKeys()`).

### Cubes at more than one height

A geometry cache carries one height. `FcgHeader.height_m` is a scalar and
is part of the cache key, so cubes at 0.35 m, 1.0 m and 2.0 m require three
caches. Build them with three `ptm_flux_cube_geom --height-m` runs and list them:

```yaml
flux_cubes:
  enabled: true
  cache_files:
    - site_h0p35.fcg
    - site_h1p00.fcg
    - site_h2p00.fcg
```

The run then makes one thermal solve and three flux-cube passes, and writes three
`.fcb` files. The terrain temperature field does not depend on the height of the
cubes, so one thermal solve serves all heights.

Outputs are named after the file stem of each cache
(`flux_cubes/flux_cubes_site_h0p35.fcb`), and every `.fcb` carries its own
`height_m` in the binary header and in the JSON. Each file is therefore
self-describing and can be passed to `scripts/ptm_flux_cubes_to_hdf5.py`, which
takes an explicit `.fcb` path. Two caches whose file-name stems collide stop the
run with an error before any pass, so no partial product is written.

A single cache (whether written as `cache_file:` or as a one-element
`cache_files:`) writes `flux_cubes/flux_cubes.fcb`, with the same flux data in
both cases.

This section is set in the configuration file only. The GUI does not expose it.

## shadow_depth

The shadow-depth pass computes the maximum height above the local surface at which the Sun is not even
partially visible (the umbral ceiling), per face and per output epoch. Like
`flux_cubes`, the block is at the top level.

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | bool | false | Run the shadow-depth pass. It does not change any other output |
| `h_max_m` | float | 5.0 | Search ceiling, meters. Crew mobility and thermal constraints |
| `march_steps` | int | 16 | March steps from the surface to `h_max_m` |
| `bisect_rounds` | int | 4 | Bisections inside the bracket the march found. Final bracket = `h_max_m / (march_steps * 2^bisect_rounds)` = 0.0195 m at the defaults |
| `solar_disk_samples` | int | 128 | Disk samples when a probe that is blocked on every ray is examined further |
| `write_frames` | bool | true | Per-epoch `shadow_depth/shadow_depth_<index>.bin` (index zero-padded to at least 6 digits) |
| `monotonicity_faces` | int | 10000 | Faces sampled by the overhang audit (0 disables) |
| `monotonicity_epochs` | int | 10 | Epochs sampled by the overhang audit (0 disables) |

Height is measured along the body radial, not along the face normal.
`ptm_flux_cube_geom` places its cubes along the same vertical, so the two
products share a height coordinate.

The Sun is treated as a finite disk, which makes this the umbral ceiling
and not the penumbral one. The disk test
(`src/cpp/illumination/solar_disk_visibility.h`) is shared with the flux-cube pass.
Neglecting the solar angular radius would change the result by
`d * sec^2(eps) * theta_sun` (0.19 m at a 20 m occluder distance and 45 deg
elevation, ten times the search resolution).

The search marches upward and then bisects. Bisection alone assumes that the
transition from blocked below to clear above is single-valued, which an overhang
violates. The march finds the first clear step from the surface upward, which
bounds the contiguous umbral column above the surface.

The output files are:

| Path | Contents |
|---|---|
| `shadow_depth/shadow_depth_<index>.bin` | one frame per output epoch (index zero-padded to at least 6 digits): `uint32` face count then `float32[F]` |
| `statistics/shadow_depth_{avg,max,min}.bin` | over resolved, nonzero epochs only. `NaN` where a face had none |
| `statistics/shadow_lit_frac.bin` | fraction of evaluated epochs with the Sun visible at the surface |
| `statistics/shadow_unresolved_frac.bin` | fraction umbral above `h_max_m` |
| `statistics/shadow_depth.json` | the manifest: search ladder, resolution, sentinel encoding, solar-disk treatment, the resolved ray `tnear`, the epoch window, and the consistency-check counts. Written last, so its presence shows that the maps beside it are complete |

Three sentinel values are used. `0.0` means the Sun is visible at the surface.
`+inf` (float32 `0x7F800000`) means no clear height was found at or below
`h_max_m`, which is the most common value inside a permanently shadowed region.
`NaN` means the face was not evaluated.

In the GUI's 3D viewer, these six products (the per-epoch series, the three
statistics and the two fractions) appear in the Field drop-down like any other
field. Because two of the three encodings are not on the color scale, the
shadow-depth fields are scaled to `[0, h_max_m]` and not autoscaled, and the
sentinels are drawn in separate colors: gray for `NaN` (all six fields), and for
`+inf` a color that the colormap in use does not contain, magenta by default.
Magenta is the top of `cool` and close to a point on `hsv`, both of which the
Colormap dropdown offers, so with those colormaps the viewer uses the first
candidate color that the colormap does not contain and reports its choice. A
pane above the view gives that legend for the selected field, and a warning if
no candidate is distinguishable. The array itself is not modified, so a probed
facet still reports `inf`.

For chunked runs, `scripts/run_chunked_simulation.py` sets
`output.skip_time_series` for every chunk, so it stops with an error on a config
with `write_frames: true`. Set `write_frames: false` (the window statistics and the
manifest still merge correctly across chunks, from `statistics/_raw` sum/count
side files) or run the window unchunked.

## gpu

These keys set the GPU safety policy that applies across solvers: whether GPU
radiosity is allowed at all, and which devices each pool may use. They are
distinct from the per-solver `gpu_devices` / `gpu_device_ids` / `gpu_batch_size`
keys, which tune a GPU solver that is already allowed. The rationale, the
platform table and the cost of each safeguard are in
[gpu-wedge-guardrails.md](dev/gpu-wedge-guardrails.md).

| Key | Type | Default | Description |
|---|---|---|---|
| `min_vram_margin_gib` | float | 2.0 | Minimum VRAM headroom per device that a resident radiosity CSR replica must leave. With less, GPU radiosity is not used (see the platform table below). 0 = fit check only |
| `placement_policy` | string | `auto` | `auto`, `shared` or `disjoint`: whether the radiosity replica and the GPU thermal pool may occupy the same devices |
| `placement_margin_gib` | float | 4.0 | Margin below which `auto` separates the two pools onto disjoint devices, and below which the thin-margin advisory prints |
| `thin_margin_policy` | string | `separate` | Action when the margin is below `min_vram_margin_gib`: `separate` places the two pools on disjoint devices and keeps GPU radiosity, and `refuse` gives up GPU radiosity |
| `hard_min_vram_margin_gib` | float | 0.25 | Under `separate`, the margin below which GPU radiosity is given up in any case |
| `dispatch_target_s` | float | 1.5 | Wall-clock budget [s] against which each GPU thermal dispatch is sized (the adaptive self-chunking of the persistent kernel). Changes dispatch boundaries only, not results. `PTM_GPU_DISPATCH_TARGET_S` overrides it. The valid range is 0.1-10 s, and other values are ignored with a warning. Recorded in `provenance.txt` |

The failure that these keys prevent, a GPU that stops responding and can take
the machine down with it, was observed on macOS with AMD (Navi21) devices, so the policy depends on the platform.

| | macOS + AMD | every other host |
|---|---|---|
| `min_vram_margin_gib` | enforced: with a thin margin GPU radiosity is given up and radiosity runs on the CPU | advisory: prints the same message and proceeds |
| `placement_policy: auto` | separates the pools when the margin is thin | never separates, and prints a note when it would have |
| `placement_policy: disjoint` | separates | separates (an explicit setting applies on every host) |

Setting `min_vram_margin_gib` in the config (or `PTM_GPU_MIN_VRAM_MARGIN_GIB`
in the environment) makes it enforced on every host.

With the 2.0 default, a site whose CSR replica fits by less than 2 GiB
has a thin margin. By default (`thin_margin_policy: separate`) the radiosity
replica and the GPU thermal pool are then placed on disjoint devices and GPU
radiosity is kept, at an estimated +7-10 min per window on a multi-million-face
polar site. GPU radiosity is given up only below `hard_min_vram_margin_gib`
(0.25), or when the pools cannot be separated. The alternative is expensive. On
a lunar polar site that fit by 0.64 GiB, CPU radiosity took 108.9 min per window
and the GPU pool 31.4 min. The margin is needed because the GPU thermal pool
allocates roughly 0.77 GiB per device on the same devices while the replica is
resident, and the fit check does not count it. Set `thin_margin_policy: refuse`
to give up GPU radiosity instead. See
[gpu-wedge-guardrails.md](dev/gpu-wedge-guardrails.md).

`auto` separates the pools only when all of the following hold: the thermal solver is on the GPU, the two
device sets overlap, the margin is below `placement_margin_gib`, the host is
macOS + AMD, and neither side has explicit device ids. With explicit device
ids, `auto` leaves the placement as given and prints a warning. The
single-GPU radiosity solver cannot select its device, so separation is not
available there, and the log reports this.

```yaml
gpu:
  min_vram_margin_gib: 2.0
  placement_policy: auto
  placement_margin_gib: 4.0
```

## output

| Key | Type | Default | Description |
|---|---|---|---|
| `directory` | string | `"output"` | Output directory |
| `format` | string | `"binary"` | `binary` or `vtk` |
| `vtk_mode` | string | `"both"` | `thermal`, `complete` or `both` |
| `include_primary_ir` | bool | false | Include the primary-body IR field in VTK output |
| `thermal_vtk_file` | string | `""` | Explicit thermal VTK path |
| `complete_vtk_file` | string | `""` | Explicit complete VTK path |
| `locations_file` | string | `""` | CSV of point locations to export |
| `series_output` | string | `""` | Time-series output path |
| `location_data` | string | `""` | Per-location output path |
| `depths_file` | string | `""` | Depth list for location exports |
| `skip_time_series` | bool | false | Do not write the per-step files. Statistics are still computed |
| `write_raw_stats` | bool | false | Write merge-friendly raw statistics for the chunked workflow |
| `map_depths_m` | list of float | `[]` | Depths (m) at which to write mesh-wide temperature maps. `[]` = off |
| `map_depth_snap_to_node` | bool | true | Snap each requested depth to the nearest solved grid node |
| `map_depth_stats` | list of string | `[avg, max, min]` | Which statistics to write per depth |
| `skip_mesh` | bool | false | Do not copy `mesh.ply` into the output directory |

### Two-band radiosity and effective-albedo products

Three products are written without a key of their own, by any run that
produces them:

| On disk | Written by | Shape |
|---|---|---|
| `indirect_solar/`, `statistics/indirect_solar_{avg,max,min}.bin` | a two-band radiosity run | per-step series + statistics |
| `indirect_ir/`, `statistics/indirect_ir_{avg,max,min}.bin` | a two-band radiosity run | per-step series + statistics |
| `statistics/effective_albedo.bin` | every run with a per-face flux path | one value per face |

All three appear in the Field drop-down of the GUI's 3D viewer. Two conventions
apply to them.

- The indirect bands are absorbed flux at the facet, while
  `direct_solar`, `primary_reflected` and `primary_ir` (same units, same mesh,
  adjacent in the same menu) are stored as incident flux. The two-band solve
  produces absorbed band fluxes, and `indirect_solar + indirect_ir` equals the
  value coupled into the thermal boundary condition (the run checks this
  identity). The viewer notes the convention in the pane above the view when
  one of them is selected.
- `effective_albedo` is `NaN` wherever a facet was never lit during the
  integration window, because there is no incident flux from which to form the
  ratio. On a 2,397,890-facet lunar polar run at 1 h cadence this was 26,057
  facets, the count the log reports as never lit. The viewer scales the color
  bar over the finite values only and draws those facets gray, with a legend
  entry for gray.

### Subsurface temperature maps (`map_depths_m`)

Mesh-wide temperature statistics at depth are written beside the surface
statistics and reduced over the same integration window:

```
statistics/temperature_z0p050m_{avg,max,min}.bin   uint32 face count + float32 per face
statistics/depth_map_z.json                        specified and sampled depths
statistics/_raw/temperature_z0p050m_{sum,count}.bin  (with write_raw_stats, for the chunk merge)
```

The file name carries the resolved depth to three decimals with `.` written
`p`, so that names sort correctly and are valid on any filesystem. Two requests
that resolve to the same name stop the run with an error.

`[]` turns the product off. Turning it on changes no other output file.

Depth maps do not appear in the Field drop-down of the GUI's 3D viewer. They
can be exported (GeoTIFF/HDF5/COG) and read from Python, and
`ptm.results.detect_available_statistics` reports them. The viewer's field
menu is built from a fixed list of names, which cannot include names that carry
the depth.

Both list keys must be written as YAML lists. `map_depths_m: 0.05` and
`map_depths_m: "0.05, 0.25"` are scalars, and either stops the run at config
load with a message that gives the correct form. `map_depths_m:` with nothing
after it is allowed and means the same as `[]` (off). Values are not split or
converted: a quoted `"a, b"` is one value, not two. The same applies to every
list key. See
[Keys the engine reads as lists](#keys-the-engine-reads-as-lists).

A depth below the deepest solved node stops the run with an error that names
`thermal.thermal_b`, the YAML key that sets the depth of the grid. The internal
name of this parameter is `grid_b`, which is not a configuration key.

Snapping to the nearest node is the default. A depth between two grid nodes is
sampled by linear interpolation, and the interpolation error can be
significant. At a lunar-equatorial diurnal amplitude *A* = 160 K it is ~0.195 K
for an instantaneous value, and it grows with harmonic number because harmonic
*h* propagates with skin depth `delta/sqrt(h)`. A depth that coincides with a
node has no such error, so the default snaps and reports the snapped depth and
the residual. Set `map_depth_snap_to_node: false` to sample the depth as
requested. The run then prints the interpolation error, and `depth_map_z.json`
carries the error for each harmonic.

The avg, max and min statistics do not have the same error. Averaging commutes
with interpolation in *z*, while max and min do not, because adjacent layers
peak at different times. For this reason `max`/`min` are accumulated in the
solver at the requested depth and not interpolated afterwards from per-layer
extrema, which would give an error of `1.229e-3 * A`, largest at layer 0.
`depth_map_z.json` therefore reports two coefficient arrays per depth.
`bias_instantaneous_per_A_by_harmonic` is a rigorous upper bound
(on a time series and, since `|max f - max g| <= max|f - g|`, on the maps too),
while `bias_max_per_A_by_harmonic` is the smaller typical error in the max.

z = 0 is the surface product. The `z=0` map is the same quantity as
`statistics/temperature_{avg,max,min}.bin`, sampled at the same instants from the
same surface node. On the CPU SIMD batch solver all three are identical to
their surface counterparts. On the GPU solver `max` and `min` are identical and
`avg` differs by the float32 running-sum error of the device (7.6e-5 K over a
118-frame window). This identity is a useful check on a new run.

The product is available on the CPU SIMD batch solver (the default),
the GPU persistent solver on any device count, the GPU radiation-coupling
selected-faces re-thermal, and both CPU scalar solvers
(`--no-simd-batch`, and the streamed solver used in its place on a large mesh)
provided no per-face thermal property map is loaded. A radiation-coupled run on
the CPU solver stops with an error. Each error message gives the remedy.

The CPU scalar solver has this condition because it builds a depth grid for
each face. The only per-face input to that grid is `H`, so without a
per-face `H` map every face has the same grid and the product is the same as on
the SIMD batch solver. The run verifies this by comparing the grid of each face
with the reference grid, and stops if any face differs. With
`thermal.h_parameter_map`, `thermal.thermal_inertia_map` or
`thermal.use_default_property_map` set, the grids differ and the run stops with
an error, before the illumination phase when the configuration alone determines
it and at the solver otherwise.

The error follows from the definition of the product. If each face were snapped
to its own nearest node, a map labeled with one depth would describe a surface
whose depth varies from facet to facet. With snapping off, every facet would be
interpolated to the requested depth, but `depth_map_z.json` carries one
`native_z_m`, one `node_j`, one `weight_on_j_plus_1` and one interpolation-error
value per depth, each of which would then be a distribution over facets. A map
on a per-face grid needs a per-face error budget, which is not written.

On the example 12,800-facet lunar crater at depths 0, 0.05 and 0.25 m, the
z = 0 map from the CPU scalar solver is identical to
`statistics/temperature_*.bin`, as on the other solvers. The maps of the CPU
scalar and CPU SIMD batch solvers agree to within the difference between their
surface products: max |ΔT| in `avg` is 0.205 K at the surface, and 0.137 K and
0.118 K at 0.053 m and 0.218 m. The depth sampling adds no error beyond the
difference between the solvers.

For chunked runs with `write_raw_stats: true`, `ptm_merge_chunks` discovers the
depth fields from chunk 0 and stops with an error if any chunk carries a
different depth set or a different statistic set. `depth_map_z.json` is carried
through to the merged directory.

Two products are not implemented: a per-depth time series (5 depths x 1440
frames at 2.4M faces would be 64 GiB) and an export of the per-layer statistics on the native
depth grid of each column.

### Per-location depth series (`locations_file` + `depths_file`)

Setting both keys writes one CSV per named site,
`locations/<name>_depths.csv`. Each file has a `# Location:` comment, a
`# Integration:` comment, a header `time,T(z=0.000 m),T(z=0.050 m),...`, then
one row per step of the whole integration, including the equilibration prefix.
(`statistics/` covers the output window only.)

The CSV is recorded during the whole-mesh solve, so the file and
`statistics/temperature_*.bin` at the same mesh face come from one integration.
Reduced over the statistics window, the CSV matches the statistics file to
within the combined precision of the two files (float32 in the `.bin`, six
decimals in the CSV).

A requested depth below the deepest solved node stops the run before the solve,
with a message that names the deepest depth available and `thermal.thermal_b`.
The example `data/locations/depths.txt` asks for depths down to 1.00 m, and the
default lunar grid (`thermal_b: 15`) reaches 0.4706 m, so this combination
gives the error. Raise `thermal.thermal_b` until the bottom node reaches the
required depth (on the Moon, `45` reaches 1.4119 m), or remove the out-of-range
depths.

The integration that produced a file is recorded in the file. Line 2 of
every CSV is either

```
# Integration: mesh column, recorded during the whole-mesh solve on the CPU SIMD batch path (grid fingerprint ...); the same integration as statistics/temperature_*.bin at this face
```

or

```
# Integration: standalone 1-D column, re-solved for this site (runThermalModelDepthSeries) -- NOT the mesh-wide integration; see the run log
```

Every reader in this repository skips `#` lines.

Three solvers record the column: the CPU SIMD batch solver (the default,
including its streamed variant), the CPU scalar solver, and the CPU
radiation-coupling selected-faces re-thermal. The CSV of a coupled CPU run
therefore carries the re-solved column at the faces affected by coupling, not
the column of iteration zero. The GPU persistent solver keeps no per-frame
column. A GPU run uses the standalone integration instead, reports this in the
log, and records it on line 2 of each CSV. The file is still written and is
valid as a standalone 1-D solution.

The CPU scalar solver with a per-face thermal-inertia / `h_parameter_map`,
which gives every face its own depth grid, also uses the standalone
integration. The recorder compares the node depths of each face with the
reference grid elementwise, not by layer count (two faces can agree on the
number of layers and differ in the top cell), and the log names the comparison
that failed.

The CSV and `map_depths_m` share the solved column, and where both are active
they give the same numbers. Over 4 sites x 4 depths x {avg, max, min}, the CSV
reduced over the statistics window matches the depth map at the face of each
site to within the combined precision of the two files, on the CPU scalar and
CPU SIMD batch solvers and with a per-face `H` map on the SIMD solver. The two
products are sampled by separate routines (the map through
`depthmap::Accumulator::addColumn`, the CSV through `depthmap::sampleColumnAtDepths`),
so the agreement is an independent check.

Where the two products have different requirements, the map is the stricter. On
the CPU scalar solver with a per-face `H` map, `map_depths_m` stops the run at
configuration time. Without `map_depths_m` the same config completes, and the
CSV comes from the standalone integration, as recorded on its line 2. In
neither case is a product written from a grid on which it was not solved.

A per-face H map changes the meaning of the CSV. The recorded column is that of
the mesh solve, and the CPU SIMD batch solver puts every face on one grid at the
reference *H* of the map. A standalone column solved at the *H* of that face on
its own grid is better resolved for the site, but it is not the value that the
mesh-wide field reports there. The CSV is intended to agree with the mesh, and
the run prints both H values at each site.

## thermophysical

On a built-in body (Moon, Mercury, Europa) the parser reads only the keys
below. Plain-named properties such as `ks`, `kd`, `H`, `rhos`, `cp0` or `albedo`
are read only under `generic_body.surface:`. On a named body the run prints one
`[config] WARNING:` line naming every such key it found and does not apply them.
See [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md).

| Key | Type | Default | Description |
|---|---|---|---|
| `use_default_albedo` | bool | true | false marks the albedo as explicitly configured |
| `constant_thermal_inertia` | float | -1.0 | Above 0 selects constant-thermal-inertia mode, which sets the radiative conductivity term to zero |
| `albedo_override` | float | -1.0 | Bond albedo override (sentinel = body default) |
| `emissivity_override` | float | -1.0 | Emissivity override |
| `ks_override` | float | -1.0 | Surface conductivity override |
| `kd_override` | float | -1.0 | Deep conductivity override |
| `H_override` | float | -1.0 | Regolith density/conductivity e-folding depth H, in meters: the scale in `rho(z) = rhod + (rhos-rhod)·exp(-z/H)` and `kc(z) = kd + (ks-kd)·exp(-z/H)`. No effect when `constant_thermal_inertia > 0` (which makes those exponentials multiply zero), and superseded per face by `thermal.h_parameter_map` / `thermal.thermal_inertia_map`. Both cases give a warning |
| `incidence_angle_albedo` | bool | false | Enable the Hapke directional-hemispheric albedo |
| `albedo_model` | string | `""` | `hapke` is equivalent to the key above |
| `albedo_hapke_w` | float | -1.0 | Single-scattering albedo (below 0 = derive from the Bond albedo) |
| `conductivity_model` | string | `"hayne2017"` | `hayne2017` or `martinez2021` (see [martinez-conductivity-model.md](martinez-conductivity-model.md)). An unrecognized value stops the run |

## Keys the engine reads as lists

Seven keys are read as YAML sequences and must be written as lists. A scalar
in place of a list stops the run at config load, with a message that names the
value given and the list form to use.

| Key | What the list sets |
|---|---|
| `output.map_depths_m` | depths of the subsurface temperature maps |
| `output.map_depth_stats` | statistics written per depth |
| `thermal.gpu_device_ids` | explicit device indices |
| `illumination.gpu_device_ids` | explicit device indices |
| `radiation.gpu_device_ids` | explicit device indices |
| `generic_body.shape.semi_axes_km` | body shape |
| `flux_cubes.cache_files` | flux-cube geometry caches |

For `gpu_device_ids` the list form protects the device selection.
`gpu_device_ids` overrides `gpu_devices` only while it is non-empty, and a
scalar would be read as an empty list. The run would then follow
`gpu_devices`: `1` (the default) selects the first enumerated GPU, and `0`/auto
selects every enumerated GPU.

`semi_axes_km` is also checked for length. It must list exactly three axes.

A key written with nothing after it (`gpu_device_ids:`, `map_depths_m:`) is a
null, not a scalar, and means "not set". Values are not split or converted: a
quoted `"0, 1"` is one value, not two.

## generic_body

Bodies without built-in properties are described here, with an ephemeris from
either Keplerian elements or a SPICE metakernel, plus a `surface:` block that
carries the full thermophysical property set. See the asteroid examples in
[`../examples/`](../examples/) and
[docs/tutorials/](tutorials/).

| Key | Type | Default | Description |
|---|---|---|---|
| `shape.semi_axes_km` | list[float] × 3 | body default | Triaxial semi-axes a, b, c in km. Must be exactly three entries. Any other shape or length is an error (see [Keys the engine reads as lists](#keys-the-engine-reads-as-lists)) |
| `shape.type` | string | spherical | `spherical`, `ellipsoidal` or `custom_mesh`. Any other value is an error. `ellipsoidal` without three positive axes, or `custom_mesh` without a mesh, builds a sphere and prints a notice |
| `ephemeris_source` | string | auto | `keplerian`, `spice` or `auto`. Any other value is an error. `auto` resolves to Keplerian (a warning is printed if a `spice_metakernel` is also set) |
| `orbit.epoch` | epoch string | J2000 | Epoch of the osculating elements. `orbit.mean_anomaly_deg` is the mean anomaly at this epoch. See the epoch formats below |
| `spin.prime_meridian_deg` | float | 0.0 | Prime meridian angle W0 at `spin.epoch`, degrees about the pole, measured from the body's reference axis: the J2000 ecliptic +x (the vernal equinox) projected onto the body's equatorial plane, or ecliptic +y if the pole lies within 25.8° of +x (`orbital/body_fixed_frame.cpp`). W(t) = W0 + 360° (t − epoch) / `spin.period_hours` |
| `spin.epoch` | epoch string | J2000 | Epoch of `spin.prime_meridian_deg` |

Both `*.epoch` keys are on the TDB scale, the scale of SPICE ET and of JPL
osculating-element epochs. The accepted formats are `YYYY-MM-DDTHH:MM:SS[.fff]`
(a space may replace the `T`, and the time and seconds are optional), the same with
a ` TDB` suffix, or a Julian date `JD 2451545.0` (the `JD` may be omitted for a
number above 100000). A `UTC`, `Z`, `TT`, `TAI` or `UT` suffix is an error and is
not converted. The value `'2000-01-01T12:00:00'`, used in the example configs,
is therefore J2000 exactly.

A DAMIT spin state (`lambda beta P` / `JD0 phi0`) maps onto these keys by
`spin.epoch: JD <JD0>`, the pole from (lambda, beta) converted to equatorial
RA/Dec, and `prime_meridian_deg = phi0 + delta`. The zero of phi in DAMIT and
the zero of W in PTM3D are different directions in the equatorial plane of the
body, so delta, the angle about the pole between them, is a constant for a given
pole. Compute it once from the two reference vectors, and confirm the result,
including its sign, by comparing the sub-solar point at one epoch with an
independent prediction.

Any key under `generic_body:` (or its `shape:`, `orbit:`, `spin:`, `surface:`
blocks) that is not read is named in a `[config] WARNING` at load
time. The provenance keys written by the HORIZONS tool (`spk_file`, `surface.albedo_*`
other than `albedo`) are named in a `NOTE` instead.

### Resolution of `generic_body.spice_metakernel`

| Value | Resolved to |
|---|---|
| Absolute path (`/…/metakernel.tm`) | Used exactly as given. The kernel directory is not prepended |
| Relative `data/…` path (`data/horizons_cache/…/metakernel.tm`) | Resolved as `mesh.import_file` and every other data key are: under `$PTM_DATA_ROOT`, else under the repository `data/` found from the directory of the config file |
| Any other relative path (`bennu.tm`) | `<kernel dir>/<value>` first, then `<current working directory>/<value>` |

The kernel directory is `$PTM_SPICE_KERNELS`, then
`$PTM_DATA_ROOT/kernels`, then the compiled-in `SPICE_KERNELS_DIR` (CMake cache
variable, default `<repo>/data/kernels`). See
[ENVIRONMENT.md](ENVIRONMENT.md). `ptm_thermal_sim` has no config key or
command-line flag for it (`ptm_camera_sim` has `--kernels-path`). Each choice
replaces the compiled-in directory and does not add to it, so an override
directory must contain the metakernels as well as `lsk/`, `pck/` and `spk/`.

`$PTM_SPICE_KERNELS` is used exactly as given, whether or not the directory is
populated. The other two candidates are tested. A directory with no
`lsk/naif0012.tls` cannot serve any run (not even a Keplerian one, because
`start_time` must still be converted), so it is skipped in favor of the next
candidate. A fresh clone, which contains the `*.tm` files but not `lsk/`, `pck/`
or `spk/`, therefore finds a kernel tree through `$PTM_DATA_ROOT`. A run whose kernel
directory has no `lsk/naif0012.tls`, for example because `$PTM_SPICE_KERNELS` names an empty
directory, stops with the SPICE error `SPICE(NOLEAPSECONDS)`. Download the kernels as described
in [INSTALL.md §5](../INSTALL.md#5-spice-kernel-data), or correct the directory.

Every run prints the kernel directory it resolved and its source (`SPICE kernel
directory: … ($PTM_DATA_ROOT/kernels)`) and the metakernel file it
loaded (`Loaded SPICE meta-kernel: …`). If none of the candidates exists the run
stops and lists every path it tried, both the metakernel candidates and the
kernel-directory candidates, in order, with the reason each was or was
not taken. A different metakernel is never substituted for the one
configured.

A `simulation.body` name other than `Moon`, `Mercury` or `Europa` that has no
`generic_body:` block stops the run at startup, with `rc 1` and a message naming
the three built-in bodies and the `generic_body:` alternative. Use
`generic_body:` for any body other than these three.

Paths inside the metakernel (`PATH_VALUES`) are resolved by CSPICE itself,
relative to the working directory, and are independent of the resolution above. The
`data/kernels/*.tm` files in the repository use `PATH_VALUES = ( 'data/kernels' )`, so a
built-in body must be run from the repository root. A metakernel with
absolute `PATH_VALUES` works from any directory. A `PATH_VALUES` entry longer than
80 characters must be split with the KPL `+` continuation.

## Generating a template

```bash
./build/bin/ptm_thermal_sim --generate-template my_config.yaml
```

## Checking a config without running it

```bash
./build/bin/ptm_thermal_sim --config my_config.yaml --validate-config
```

This option parses the file, applies any command-line overrides, runs every
startup check that depends on the configuration alone, and exits. `rc 0` =
accepted, `rc 1` = rejected, with the error message on stdout/stderr. Nothing is
simulated, no output directory is created, and an existing output directory
(including its `heartbeat.txt`) is not modified.

The check covers YAML syntax and structure, the resolution of data paths under `PTM_DATA_ROOT`,
the enum values (`conductivity_model`, `coupling_mode`, `exchange_selection`,
`gpu.placement_policy`, `gpu.thin_margin_policy`), the numeric ranges
(`nearfield_quadrature_order`, `min_shadow_fraction`, the `gpu.*` margins), the
`beaming`/`rms_slope_deg` pair, and whether this build has the GPU support the
config asks for.

It does not cover anything that needs the mesh, the body or SPICE: mesh
loading and units, generic-body construction, SPICE kernels, property-map face
counts, the `martinez2021` density check, view factors, and the spin-up seed
strings. The output lists these exclusions each time. `rc 0` means the
configuration is well formed, not that the run will complete.

The *Check Configuration* item in the GUI's *Simulation* menu runs this check and shows its
output verbatim.

## Overriding from the command line

Command-line options override YAML settings:

```bash
./build/bin/ptm_thermal_sim --config config.yaml --duration 7.0 --timesteps 168
```

The full flag list is in [USER_GUIDE.md](USER_GUIDE.md) and from
`ptm_thermal_sim --help`.

## Built-in bodies

| Body | Day length | Year length | Distance (AU) |
|---|---|---|---|
| Moon | 29.53 days | 365.25 days | 1.0 |
| Mercury | 176 days | 87.97 days | 0.387 |
| Europa | 3.55 days | 11.86 years | 5.2 |

These three have pre-defined thermophysical property sets. Any other body name
stops the run at startup unless the config has a `generic_body:` block (see
[generic_body](#generic_body)). With that block, the run starts from placeholder
properties, prints a notice on stderr, and replaces them with the values in the block.
