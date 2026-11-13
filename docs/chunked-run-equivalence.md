# Chunked-run equivalence and the output epoch grid

A chunked run divides one long simulation into consecutive windows and merges the
statistics of the windows. The merged result should equal that of a single run over
the whole span. This note describes the construction of the output epoch grid that
ensures this, and quantifies the differences that remain.

## Requirements

**R1, tiling.** The output epochs of the chunks must partition the output epochs of
the single run, so that every epoch is sampled exactly once. A duplicated epoch gives
one instant double weight in every reduction (`temperature_avg`,
`illumination_avg`, the Schorghofer ice accumulator), and a skipped epoch omits an
instant.

**R2, handoff.** The state that chunk *k* passes to chunk *k+1* must be the state at
the first output epoch of chunk *k+1*, not at the last output epoch of chunk *k*.

## Construction of the grid

The grid is half-open. With cadence `dt = duration_days*86400/timesteps`, the
single run samples `t0 + i*dt` for `i` in `[0, N_total)`, and chunk *k* is given the
sub-range `[c_k, c_k + N_k)` of the same grid. R1 then holds by construction.

R2 requires one more grid point than is written to the output, because the solver
integrates `total_time_steps − 1` intervals and its final state is therefore at its
last grid point. Chunk *k* therefore

- integrates `N_k` intervals over `N_k + 1` grid points, `[t_k, t_k + N_k*dt]`
- exports and reduces only the first `N_k`
- writes `final_state_file` at the extra point, which is `t_{k+1}`.

This is the setting `thermal.checkpoint_lookahead_steps: 1`. The extra point is not
written as a frame, does not appear in `time_steps.txt`, is not included in
`statistics/`, and is excluded from the forcing mean of the equilibration seed. Its
only effect is on the epoch of the checkpoint.

Chunk offsets are given as `simulation.start_time_offset_seconds`, an offset in
ephemeris seconds, and not as a shifted UTC `start_time` string. TDB − UTC has a
periodic relativistic term (~1.7 ms amplitude, annual period) and changes by a whole
second at a leap second, so an offset in UTC does not move a boundary by an exact
number of cadence seconds.

## Epoch grid in versions before 2026-08-12

In versions of PTM3D before 2026-08-12,
`IlluminationModel::calculateIlluminationFromMesh` built one uniform grid as

```cpp
et_step = (et_end - et_start) / (time_steps_ - 1);
```

with `time_steps_ = equilibration_steps + output_steps` and `et_end - et_start`
the combined span of equilibration and output. The output cadence was therefore
obtained by dividing the combined span by one fewer interval than there are samples.
This had three consequences, quantified here for an 80 N lunar crater (12,800 faces,
236 output steps over 59.06 days, one equilibration lunation, CPU, radiation off).

1. The cadence was too long by the factor `N_tot/(N_tot − 1)`: 21683.362 s instead of
   21621.966 s, or 0.28 %.
2. Output frame 0 was displaced from the nominal start by
   `equilibration_days*86400/(N_tot − 1)`. A run with a `start_time` of
   `2026-04-05T00:00:00` had its first frame at `2026-04-05T01:59:53.722`.
3. `N_tot` depends on the length of the run, so chunks of the same simulation
   had different cadences and their grids did not tile. The largest epoch
   misalignment was 14,438 s with 2 chunks, 18,061 s with 4, and 18,798 s with an
   unequal division into 5, with one nearly duplicated epoch at each chunk boundary.

The second consequence applies to every run with an equilibration period, and the
first to every run, including single runs. The epoch labels of products generated
with those versions are in error by these amounts. The physics at the shifted epochs
was solved correctly.

A further 0.5 s of quantization arose because the grid was anchored on
`calculateAdjustedStartTime()`, which passes through
`et2utc_c(..., "ISOC", 0, ...)` at whole-second precision. The grid is anchored on
`start_time` directly.

## Comparison of chunked and single runs

The comparisons use the 80 N lunar crater (12,800 faces, 236 output steps over
59.06 days, CPU SIMD, radiation off). The first table gives the differences in the
epoch grid between chunked and single runs.

| chunks | earlier grid: max offset | earlier grid: near-duplicate epochs | half-open grid: max offset | half-open grid: near-duplicate epochs |
|---|---|---|---|---|
| 1 | 0.000 s | 0 | 0.000 s | 0 |
| 2 | 14 438.444 s | 1 | 0.000 s | 0 |
| 4 | 18 060.805 s | 3 | 0.000 s | 0 |
| 5 (48,47,47,47,47) | 18 798.149 s | 4 | 1.1e-6 s | 0 |

`illumination_avg`, `illumination_max`, `direct_solar_avg`, `direct_solar_max`
are identical between the chunked and single runs at every chunk count. They do not
depend on the thermal solver, so they test the epoch grid alone. With the earlier
grid they differed by 0.55 W m<sup>-2</sup> in the mean and 5.7 at most with 2 chunks.

The second table gives the per-frame surface temperature for 4 chunks. All chunks
were warm-started from a shared checkpoint, so that the cold-start seed does not enter
the comparison.

| | earlier grid | half-open grid |
|---|---|---|
| chunk 0, all 59 frames | — | identical |
| chunks 1-3, interior frames | — | 3.05e-5 K max (float32 rounding of the checkpoint) |
| chunks 1-3, frame 0 (the chunk boundary) | — | up to 37.8 K |
| merged `temperature_avg` | 0.2585 mean / 9.798 max K | 0.00081 mean / 0.160 max K |

The mean difference of the merged result is reduced by a factor of 320 and the
maximum by a factor of 61.

## Remaining differences

Two differences remain, and neither arises from the epoch grid.

The first concerns the frame at a chunk boundary. In the comparison above, frame 0 of
every chunk after the first contained the surface temperature restored from the
checkpoint, and not the surface temperature solved at the epoch of that frame.
Because the surface node has no heat capacity, its temperature reflects the epoch at
which it was last evaluated. The difference was up to 37.8 K in that one frame, and
8e-4 K when averaged over 236 samples. In the CPU solvers (SIMD and scalar), frame 0
of a warm-started chunk is solved at its own epoch, so a chunked CPU series equals the
single-run series at every frame, including the chunk boundaries
(`scripts/test_chunk_seam_equivalence.py`). The GPU solver reports the surface node
of the checkpoint at frame 0 and prints a warning when it restores a checkpoint.

The second is that the equilibration seed of a cold-started first chunk depends on
the length of the chunk. The `spinup_seed: physical` seed sets the temperature of
each column from `<F_abs>`, the time mean of the absorbed flux of the face over the
whole forcing table. That table has `equilibration_steps + output_steps` entries, so
chunk 0 of a division into 4 averages over 177 steps where the single run averages
over 354, and the two runs start from different initial conditions. On the crater
above, a chunked run with a cold-started first chunk differs from the single run by
0.717 K in the mean and 13.05 K at most in `temperature_avg` with 4 chunks, compared
with 1.141 and 13.50 with the earlier grid. The chunked and single runs are therefore
not equivalent in this configuration. For the same reason, the initial condition of a
single run changes if more output days are requested at the same equilibration
length.
