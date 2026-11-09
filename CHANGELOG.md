# Changelog

## Unreleased

- **History marcher M1 (part 2): `ptm_history`, the streamed spin-up history job. No simulator change.**
  `ptm_history JOB.yaml [--threads N] [--nice N]` marches ~2 eclipse years of real forcing before a window and
  writes a PTMT v2 handoff state at the window start minus `fine_orbits` lunations (for `PTM_SEED_STATE_FILE`
  / a warm start). Method (Fermi REC, decisions 2026-11-09): seed from 480 samples of the first 12 lunations
  with one IR sweep from radiative-equilibrium walls, periodic replay at P/120 x3 + Aitken; then every epoch
  one canonical width-3 product gives this epoch's IR (from the marching surface) and the visible term with
  TWO bounces (approved: single bounce passed winter by only 0.01 K); the production CPU column, loaded once
  per history step; P/120 substeps on old history, P/768 on the last 10 lunations.
  - **Union epoch grid is mandatory**: the forcing must be a contiguous run of the site's union grid
    (`epoch_grid`), else FATAL (`allow_foreign_epoch_grid: true` runs anyway, loudly).
  - Pausable: `PAUSE` file or SIGUSR1/SIGTERM/SIGINT -> checkpoint, exit 75; resume = rerun (refused while
    `PAUSE` exists); checkpoints every `checkpoint_every_lunations`; identity (forcing, CSR, epoch grid,
    column, params, code) checked on resume. The depth grid must equal `grid_reference_state`'s.
  - MEASURED (n4990 W6h, radiation on): engine handoff vs Fermi's numpy REC marcher <= 7e-13 K; window run
    seeded from it (eq 2, coupled) vs the 64-orbit truth: surface max 0.14 K (Sep 2030), 0.24 K (Feb 2031),
    0 class flips; a truth-seeded control scores 0.00 K. Paused twice and resumed on 3 and 5 threads, the
    handoff is byte-identical to the uninterrupted 8-thread run.
  - Input today is a PTMHF001 forcing file (per-epoch absorbed + reflected direct); the in-process
    illumination stage (classifier + chord + Earth cache) is the next milestone.

- **History marcher M1 (part 1): canonical SpMV entry point, `.hzc` v2, pause/resume. No simulator change.**
  - `SparseVFOperator::matvecInterleavedCanonical(width, X, Y)`, width 1..4, ALWAYS the canonical per-row
    order, on every ISA. The marcher uses it (decision 2026-11-09): it is a new path with no legacy result to
    reproduce, so the reason the canonical kernel is opt-in on x86 (it moves x86 results in the last bits
    relative to the legacy loop) does not apply. MEASURED on the AGT near-final referee CSR (1.29e9 nnz):
    width 2, Mac Pro 0.131-0.149 s canonical vs 0.233 s legacy (the x86 legacy two-input loop costs 1.8x one
    product); MacBook 0.060 s. Width 3 (IR + visible + second visible bounce, one traversal): Mac Pro
    0.186 s, MacBook 0.088 s. Rows are reduced independently in a pinned order, so the result does not depend
    on the thread count. The existing `matvec`/`matvec2` (and `PTM_SPMV_CANONICAL`) are unchanged.
  - `.hzc` file v2 (v1 float32 files are refused with a rebuild message): 16-bit elevation codes rounded
    toward STRADDLE (a round trip can lose a certification, never gain one), MULTI as the (+inf, -inf) pair,
    column runs, streamed writer (bounded RAM), memory-mapped reader with a whole-file checksum.
    `ptm_hzcls` builds tile by tile (`HZ_BUILD_ONLY`, `HZ_TILE`, `HZ_CENSUS_ONLY`); `HZ_ROUNDTRIP` is FATAL
    on any gained certification. MEASURED, one AGT site, union of its 8 windows (14,029 epochs, 1.80M faces):
    6.48 GB, 417 s, peak RSS 1.62 GB (MacBook, 18 threads); the MacBook- and Mac Pro-built files are
    byte-identical on a 20k-face sample.
  - `src/cpp/history/`: `history_checkpoint` (atomic checkpoints with the job identity checked field by
    field, fallback to the previous checkpoint when the latest is corrupt, pause by `PAUSE` file or
    SIGUSR1/SIGTERM/SIGINT, resume refused while `PAUSE` exists) and `test_history_checkpoint`, the gate:
    a toy marcher on the engine's own kernels, paused twice and resumed with different thread counts, must
    end byte-identical to the uninterrupted run.
  - `bench_thermal_batch`: wall per whole-mesh substep of the production CPU column at a chosen face count,
    for the engine's per-substep load/save driver and a per-history-step driver.

- **History marcher M0: horizon-classifier library, `ptm_hzcls`, and `bench_vf_spmv`. No simulator change.**
  New STATIC library `horizon_classifier` (illumination/): the per-face solar horizon envelope that screens the
  history marcher's direct-solar term (CLEAR / BLOCKED / STRADDLE -> exact chord), with the face-plane guard,
  a fall-through census (an unvisited column or a disk outside the built ladder always goes to the exact chord,
  which makes one per-site envelope safe to reuse across windows), a deterministic 0.1% audit selection, and
  the keyed per-site `.hzc` file (scene hash of the traced float32 buffers + ray-predicate hash; never a
  triangle count; float32 values rounded conservatively). Built without `-ffast-math` so its boundary
  comparisons are IEEE-exact. `ptm_hzcls` builds/applies `.hzc` files; `test_horizon_classifier` gates it
  (synthetic geometry, plus `--parity` against the measurement tool: n4990 and an AGT 20k-face sample are
  byte-identical). `bench_vf_spmv` times the production `SparseVFOperator` products on a referee CSR
  (read-only). None of these is linked into `ptm_thermal_sim`.

- **Bottom boundary uses the bottom link's conductivity, k_{N-2} (PTM-853). Default output changes at
  the mK level.** The geothermal-flux boundary set T_{N-1} = T_{N-2} + Qb dz/k with k at the deeper
  node, so the bottom link conducted Qb k_{N-2}/k_{N-1} instead of Qb (about 1e-3 low on the Moon grid).
  Changed on all three paths (`botTemp`, `applyBottomBCBatch`, the GPU persistent kernel). Measured
  default-output change, warmer wherever it moves: Moon G3 4e-5 K; Europa G5 2 mK (ice on 72 faces, <= 1.7e-3 cm);
  Moon W6h production 4 mK; W6h radiation off 4.3 mK at the surface and 4.9 mK at 7 cm, with no ice node
  changes and no 54 K / 110 K class flips. Generic bodies with Qb = 0 are unchanged. On the GPU, sunlit
  and seasonally lit faces move as on the CPU (W6h: median ratio 1.035, worst difference 1.4 mK). Dark-floor
  columns (below ~30 K) do not move: their drift is below one float32 ulp per substep, so float32
  freezes it (PTM-858). The GPU's dark floor was already where the CPU's lands after the fix.
  - **Old thermal states read Qb about 1e-3 low.** A state written before this change carries the
    old bottom gradient, so Qb inferred from it (`eqQbFromState`, `ptm_state_info`) reads
    Qb k_{N-2}/k_{N-1}.
  - **Do not resume a chunked campaign across this change.** Its earlier chunks and the warm-start
    checkpoints were integrated with the old boundary. Rerun the campaign from the start on one
    engine version.
- **Opt-in implicit-midpoint column step, CPU and GPU** (`thermal.column_step: midpoint`; PTM-846, PTM-855).
  The default stays `legacy`, which is byte-identical. The step couples the massless surface into the implicit
  half of the top Crank-Nicolson row, and re-linearizes the surface, k and cp at the step midpoint.
  `thermal.column_step_passes: 2` is the scheme (second order in dt). 0 and 1 run as labeled
  diagnostics; 3 or more are refused, because they are worse at flash faces. Substeps per output interval:
  n_sub = max(8, ceil(dt_out / (P/1024))), with P the solar day, so a substep never exceeds P/1024 on the
  CPU (the `PTM_CN_SUBSTEP_DIV` diagnostic is refused beyond P/384); for the GPU's cap, see below.
  The frame surface and the ice accumulator's surface node are re-solved at the frame epoch with a
  self-consistent k0.
  - **Accuracy.** Moon n4990, radiation off, against a P/49152 reference:
    - 5.9 h cadence: max 0.27 K, p99 0.22 K (legacy: 3.0 K / 2.2 K).
    - 1 h cadence: max 0.24 K (legacy: 7.2 K on the CPU, 1.37 K on the GPU). The legacy GPU is closer
      because its substep is limited by the top cell's Fourier number (~16 substeps per hour here), while
      the legacy CPU takes a fixed P/768 (2 substeps per hour).
  - **Cost.** CPU thermal time x2.5 at 5.9 h and x3.0 at 1 h. GPU kernel time per face-step, against the
    legacy GPU: x1.34 in a 1 h output window, x0.28 at 5.9 h (the cadence of the default equilibration
    prefix). A whole run's ratio depends on its prefix:window split; for example, an 8-orbit prefix at
    5.9 h plus a 29.5 d window at 1 h comes out at x0.47 in thermal wall time. The GPU runs it as a
    separate OpenCL program (`-DPTM_COLUMN_MIDPOINT`), so the legacy kernel is unchanged. The GPU and the CPU agree to 0.011 K at the same substep count.
  - **Watchdog cap on the GPU.** One output step may take at most 200 substeps on the GPU. The rule's
    count reaches that cap at dt_out = 0.195 P. Between 0.195 P and 0.26 P the run continues with 200
    substeps and a warning; beyond 0.26 P, where the capped substep would exceed P/768 (above the 0.5 K
    error bar), it is refused before any illumination or GPU work. For the Moon that is 5.8 d and
    7.7 d. **Europa at a 1-day cadence is refused** (P/710); the longest accepted Europa cadence is ~22 h.
  - **Refused (fatal), with the reason printed:** radiation-coupled runs, lateral conduction and
    time-outer runs, the CPU scalar path, and `PTM_CN_SUBSTEP_DIV` with the GPU.
  - **Warm starts and checkpoints warn.** A midpoint run with `initial_state_file` or `final_state_file`
    prints a warning: chunked and warm-start equivalence is not yet validated on this step.
  - **Validation scope.** Moon n4990 only, with radiation off. Europa, generic bodies, Martinez conductivity
    and chunked runs on the midpoint step are not yet validated.
- **GPU substep watchdog cap is now reported when it binds** (PTM-857). The legacy GPU's Fo-limited
  substep count is floored at max(1 s, dt_out/200). At a 1-day output cadence that floor binds on about
  95% of the faces of a Moon crater mesh; at 5.9 h and 1 h it never binds. The first binding prints a
  warning with its face count, and one summary line at the end of the run gives the totals. Outputs are unchanged.
- **Opt-in energy ledger, off by default** (`PTM_ENERGY_LEDGER=<csv>`; PTM-846). It writes a per-face,
  per-diurnal-cycle budget in double precision: absorbed, emitted, geothermal, T_geo-floor enthalpy,
  top and bottom conducted flux, and the column enthalpy change. Its residuals are normalized to
  max(absorbed, emitted). Products are byte-identical with it on or off. It is implemented on the CPU
  SIMD batch-outer driver only; on the GPU, the time-outer/lateral driver and the scalar path it is
  fatal rather than silently absent.

- **Flux cubes: the body sphere now applies to cube rays on open patches** (PTM-823). The
  flux-cube pass traced its Sun and primary-disk rays through Embree alone, so on an open patch a
  cube ray aimed below the limb that missed the terrain counted as clear sky. Each cube ray now
  first takes the direct Sun's R_occ test, under the same open-mesh gate. Each pass also logs, per
  site, the mesh horizon against the sphere's cap edge, plus the number of cube rays the sphere
  blocked, to the log and to `provenance.txt`. The `.fcb` layout is unchanged. Wherever the mesh
  horizon stays above the cap, the products are identical: on the crater80N test case, all three
  heights match the previous engine outside the engine identity, and every surface output stays
  byte-identical.
- **Chord direct solar about 2x faster, byte-identical** (`solar_disk_scheme: chord`, CPU batch path).
  Every chord ray of a face starts at its centroid on the default rise-epsilon path, so the screen
  and the edge check now go to Embree as 16-ray packets (`rtcOccluded16`, each lane the one-ray
  query's float ray and filter inputs). The edge check traces its points predicted visible and
  predicted blocked as separate packets and still answers with the first wrong-side point in its own
  order. The screen stops after the first packet that has seen both a clear and a blocked probe. The
  per-face class counters are per thread instead of shared atomics (a fully culled 1.8 M-face frame:
  ~46 ms in the timing-kit log, 2.5 ms measured now). The column bisection and the fallback samples stay one ray at a time.
  AGT polar timing-kit window, illumination phase: 1,445 -> 959 s (direct solar ~1,010 s inferred -> 524 s measured).
  Byte-identical to the previous chord build on 21 gates (chord on: G3 x5 env/LD variants, G4, G5,
  AGT 72 h x2, Haworth30; legacy: G3/G3lat/G4/G5, production and KEPS env, AGT 72 h, W6h,
  Haworth30), on the restated-G2 Moon3D set and on two AGT days. `provenance.txt` keeps the
  rule's ray count. Diagnostics: `PTM_SOLAR_CHORD_PACKETS=0|4|8|16`, `PTM_SOLAR_CHORD_PROFILE=1|2`,
  `PTM_SOLAR_CHORD_DUMP=<dir>` (each announced). The library is built with `-ffast-math`: a chunked
  fallback was tried and flipped one grazing sample in 1.3e8 face-frames, so it was reverted (PTM-848).
- **Opt-in, results-changing chord settings, all off by default:**
  `illumination.solar_disk_chord_check_adaptive_n` (fails the restated G2: max dF 42-45 W/m2),
  `solar_disk_chord_fallback_n` (512 passes G2 but saves ~1%), `solar_disk_chord_temporal`
  (-5% direct solar at the AGT site; G2 cannot exercise it, PTM-849). Each prints a NOTICE,
  lands in `provenance.txt`, and is fatal with a scheme other than `chord`.

- **`radiation.gpu_batch_size` defaults to `auto`, sized from free host RAM** (multi-GPU
  radiosity pool). At each sweep the batch is the largest that 0.25 x (available RAM - reserve)
  affords at 72 bytes per face per timestep, capped at the sweep (split in two when each half
  is >= 150 timesteps), evened over the sweep and rounded to the device count when that keeps
  the batch count (otherwise the device taking the remainder steps rotates per batch; next entry), and floored at 4 x devices.
  On the 1.8 M-face polar site it was 48 timesteps under the old fixed 4 GiB budget; the measured knee
  is ~150-192 (radiosity sweep -27..-29% vs batch 49). The batch and its inputs are logged and
  written to `provenance.txt` (`rad_pool_batch=`). Results are byte-identical at any batch size.
  The per-timestep host cost was under-counted as 48 B/face (measured ~63; PTM-833). An explicit
  value below the floor is now raised **with a warning** (PTM-829; it was silent). `auto` is
  accepted in YAML and by `--rad-gpu-batch-size`. Cooling: larger batches mean longer GPU work
  between rests; the pause remains (see `docs/dev/gpu-wedge-guardrails.md`, G7).
- Multi-GPU radiosity sweep: less host work between batches, and the pool's inter-batch
  pause is now visible. The coupling sweep gathers each batch face-major into buffers
  it reuses, moves the per-timestep source albedo into the solver instead of copying it,
  and the solver and manager hand their per-band buffers on by move instead of
  zero-filling and copying them. A batch the host budget caps is evened out to whole
  steps per device when that keeps the batch count (a 312-step sweep: 49 -> 48 on 4
  devices; the device holding the 49th step finished last in every batch), and otherwise
  the device that takes a batch's remainder steps rotates from batch to batch.
  `PTM_LOG_GPU_DISPATCH=1` now prints the radiosity pool's per-device work and pause and
  the sweep's per-batch host phases. On a 1.8 M-face core+ring site
  with `radiation.gpu_inter_batch_pause_ms: 3000`, host time between batches fell from
  ~2.4 s to ~1.1 s and the sweep from 141.0 s to 132.4 s (-6.1%): the sweep is bound by
  that pause, which each device sleeps out after every batch, so most of the host time
  saved is now spent in the pause (12 s -> 158 s slept). Outputs are byte-identical.
- **New, opt-in: solar limb darkening** (`illumination.solar_limb_darkening: none | neckel_labs`,
  default `none`). This is a bolometric law built from Neckel & Labs (1994) Table I, weighted by the
  AM0 disk-center spectral intensity, with F/I = 0.829. Fluxes are normalized by the disk-mean
  intensity, so a fully visible disk still delivers `S cos i`; partly visible disks and slivers
  change. It applies to `solar_disk_scheme` chord and reference, and is fatal with
  legacy. See `docs/solar-disk-scheme.md`.

- **New, opt-in: probe-then-chord solar disk** (`illumination.solar_disk_scheme: legacy | chord |
  reference`, default `legacy`, byte-identical). The legacy rule returns the full disk
  whenever the disk center is visible, drops thin slivers with its 9 samples, and tests the body sphere
  on the center ray only. `chord` screens the disk with 25 probes, bisects the edge in three
  columns (10 halvings), integrates it exactly, checks the fit with 128 rays, falls back to 1024
  samples, and tests the body sphere per ray. CPU batch path; GPU runs fall back to the CPU tracer for direct solar, with a warning;
  eclipse frames keep the legacy heavy path. `provenance.txt` records the scheme and its face-frame
  classes. See `docs/solar-disk-scheme.md`.

- Direct solar on the CPU, rise-epsilon path (`PTM_ILLUM_RISE_EPS=1`, the production illumination
  bundle): each shadow ray is now ONE Embree any-hit query with an argument filter that rejects the
  source triangle and hits rising less than h_eps, instead of a closest-hit march of 1.7-4.6
  `rtcIntersect1` calls. That is the GPU kernel's semantics. Finite-disk samples at or below the face
  plane, which add exactly 0, are no longer traced, and the per-face center ray no longer allocates
  heap vectors or opens a nested OpenMP region. Outputs are byte-identical: G3, G3 lateral, G4 and G5 in
  the default and production environments, G3 with keps 0.2 and with TRUE_NORMALS off, an AGT polar
  72-h window in three environments, W6h and the Haworth R1 30-frame run. Illumination wall time on the
  AGT 72-h window: 148 s -> 94 s. `PTM_ILLUM_RISE_FILTER_AUDIT=1` re-runs the march next to every
  filtered ray and prints the number of disagreements. RISE_EPS2 keeps the march (its OR layer needs
  every hit).

- The GPU thermal dispatch budget is now `gpu.dispatch_target_s` (default 1.5 s,
  was a fixed 2.0 s; env `PTM_GPU_DISPATCH_TARGET_S`; printed at startup and
  written to `provenance.txt`), and the adaptive self-chunking aims at 0.9 of it,
  so the budget is a ceiling: on a 509,337-face Moon3D run with 25 equilibration
  orbits the longest sub-dispatch went from 2,052 ms to 1,367 ms. Dispatch
  boundaries only; outputs are byte-identical.
- Measured GPU settings for the AMD W6800X workstation are in
  `docs/dev/gpu-wedge-guardrails.md` (G7): `thermal.gpu_batch_size: 10240` for
  single-GPU thermal runs (2.8-3.7x faster than 1280, identical results);
  `thermal.gpu_inter_batch_pause_ms` acts only on the multi-GPU pool and the
  radiation re-thermal. GUI tooltips now say what these settings do.
- **Earthshine and planetshine no longer reach faces through the body** (PTM-819;
  `illumination.primary_disk_body_cull`, new, default `true`). On an open patch the primary-disk
  probe and sample rays had no body-sphere test, which the direct Sun has: a face tilted toward a
  primary that had set below the limb received its IR and reflected light. Every primary-disk ray
  now takes the direct Sun's forward ray–sphere test on R_occ, under the same closed/open topology
  gate. Only faces whose disk reaches the body cap are re-evaluated, so every other face is
  byte-identical to before. On the Moon3D 85°S crater over one lunation (40 frames), the fix
  zeroed every face epoch whose Earth disk was wholly behind the body (19% of Earth-lit face
  epochs). It reduced every face epoch whose disk straddled the limb (12%), and it left every other
  face epoch bit-identical. That is checked against an independent float64 test. The removed
  leak is Earth IR plus reflected Earthlight: up to 0.072 W/m² of IR and 0.081 W/m² of reflected
  light per face (incident). Over the lunation (360 steps, radiation off), the 54 K and 110 K
  cold-trap areas did not change. T_max fell by more than 1 K on 2.3 cm², and by up to 5.9 K on 13
  cold faces tilted toward Earth, where the absorbed leak removed was 0.075–0.139 W/m². On the
  Haworth crater over a full year of daily frames, 1.93% of Earth-lit face epochs had the disk
  behind the body and 2.88% straddling. The largest absorbed Earthshine on a face epoch whose disk
  reaches the body (an upper bound on the leak) was 2.86e-2 W/m², outside the ±15 km Diviner
  comparison window; inside it, on the cold-trap faces (T_max < 110 K), at most 4.7e-4 W/m².
  On the Europa Pwyll regional patch (placed on the equator 88° from the sub-Jovian point, where
  Jupiter sits on the horizon), Jupiter's disk straddled the body cap on
  96.6% of the Jupiter-lit face epochs, yet the fix changed 621 of 28,300 faces. The patch-mean
  Jupiter IR plus reflected light fell by 0.03%, by at most 6.4e-4 W/m² absorbed on any face epoch,
  and T_min, T_avg and T_max moved by at most 0.0015 K. Direct solar stayed byte-identical. There
  the re-evaluation pass doubles the primary-disk rays: about 3% more wall clock for the full run,
  44% for an illumination-only run. Unaffected, byte-identical: closed meshes and bodies without
  a primary. Also unaffected: open patches whose primary never sets below the limb. `false`
  reproduces earlier runs and prints a warning. `provenance.txt` records `primary_disk_body_cull` and
  the number of rays it blocked.

- CI: the macOS legs (`python-layer` and `build-and-test` on `macos-latest`) are
  temporarily non-blocking (2026-11-04). GitHub-hosted macOS arm64 runners are
  capacity-constrained and jobs were failing without ever getting a runner. The macOS legs
  still run when a runner is available, but a failed, canceled or unacquired macOS job no
  longer fails the workflow; every run that plans one prints a warning saying so. Linux
  jobs stay blocking. To be reverted once capacity recovers.
- **New, opt-in: seasonal (two-timescale) spin-up seed** (`thermal.seed_seasonal: on | auto`,
  default `off`; `--seed-seasonal`, `PTM_SEED_SEASONAL`). A window-periodic spin-up misses the
  annual thermal wave that the deep column carries at polar sites. The run illuminates 12
  phases of one diurnal cycle each over the seasonal period, solves each face's periodic
  seasonal storage state on its own column, and moves the seed to it; coupled runs add a
  per-phase indirect estimate from the radiosity operator. Measured against continuous
  real-ephemeris references: dark-season polar windows from -2..-7 K (lit) and -16..-25 K
  (seasonally shadowed) T_min error to within ~0.5 K and ~1 K; a bright-season lunar crater
  window from +1.4 K to -0.7 K (lit T_min) and +0.64 K to -0.13 K (cold-trap T_max).
  The default is unchanged (byte-identical) pending review.
- **GPU coupling re-thermals now follow the CPU's method, and chunk 0's equilibration carries a
  time-varying indirect flux** (PTM-812; `thermal.gpu_rethermal`, new, default `seed_series`;
  `thermal.eq_indirect_subcycle` gains `phase`, the new default). The GPU re-thermal resumed
  iteration 0's post-equilibration column; the CPU re-seeds and re-equilibrates on each pass's
  indirect flux. Under `seed_series` the GPU does the same, and its whole-mesh solve keeps the
  equilibration span's temperatures so the coupling sweep solves that span too. On a sub-cycle
  window the equilibration used to carry a constant mean indirect flux, which by Jensen's
  inequality biases the column warm; `phase` folds the sweep's own indirect over the last whole
  equilibration cycle onto the equilibration steps by phase. `resume` and `legacy` reproduce the
  earlier outputs. On a 968-face crater, `--chunks 2` minus `--chunks 1` T_min is now -0.002 K on both
  paths with the same seed path in both arms (was +2.02 K GPU, +0.05 K CPU); on a 4990-face crater the
  GPU agrees with a converged CPU run to 0.05 K (T_min), 0.08 K (T_avg) and 0.09 K (T_max). Cost: each
  GPU coupling re-thermal now carries an equilibration phase, and the coupling sweep covers the
  equilibration span; that run went from 52 s to 171 s, and a 1.8M-face production window took 46 minutes
  against 28 for the earlier build (not a matched A/B). That production window ran with `seed_window` and
  `eq_indirect_subcycle` pinned to `legacy` and only 2 equilibration cycles (the engine warns that its p95
  requirement is 40), on a full-cycle window, so it never exercised the new `phase` default. A warm-started
  chunk restores its checkpoint under either mode, with a one-line note. `provenance.txt` records
  `gpu_rethermal`. Also `PTM_GPU_RETHERMAL`.
- **GPU thermal kernel interpolates the surface flux within each output interval** (PTM-808;
  `thermal.gpu_flux_forcing`, new; default `interpolate`). The persistent GPU kernel held each
  interval's starting flux for every substep, while the CPU solvers interpolate linearly at the
  substep midpoints. On a 4990-face crater at 6 h cadence that put GPU T_max 1.27 K low on
  average and up to 21.9 K low against the CPU. The kernel now interpolates the same way; the
  output frame is still the surface solved with the frame's own flux against the interior at that
  epoch. Against a CPU run with the same re-thermal start mode and converged substeps, GPU T_max
  now agrees to 0.14 K max (0.004 K mean) on that mesh. At production substeps the CPU itself
  sits up to 1.6 K off its own converged T_max there (PTM-810). `hold` reproduces earlier GPU
  outputs byte-for-byte, with a warning. CPU outputs are unchanged. Also `PTM_GPU_FLUX_FORCING`;
  `provenance.txt` records `gpu_flux_forcing`.
- **Chunked runs: chunk 0's equilibration indirect flux is a whole-cycle mean** (`thermal.eq_indirect_subcycle`,
  new; default `whole_cycle`). When the output window is shorter than one cycle, the internal eq-indirect pass used to
  install the output-window mean, a phase-biased sub-cycle mean (3.6 vs 2.3 W/m^2 on a 968-face crater). It now
  installs a whole-cycle estimate built from the equilibration span's illumination tables and the radiosity operator,
  identical on the CPU and GPU thermal paths. On that mesh, `--chunks 2` minus `--chunks 1` T_min went from +2.45 K to
  +1.58 K on the GPU thermal path (-0.44 K to -0.49 K on the CPU). With the seed path forced the same in both arms it went from
  +2.86 K to +2.02 K (GPU) and +0.10 K to +0.05 K (CPU). The remaining GPU split matches a CPU run with
  `PTM_CPU_RETHERMAL=resume` to 0.02 K, so it is the re-thermal start mode (PTM-808), not the indirect flux. Runs whose
  output window covers a cycle are byte-identical. `legacy` reproduces earlier outputs byte-for-byte. `provenance.txt`
  records `eq_indirect_subcycle` and `eq_indirect_subcycle_effect`.
- **Chunked runs: chunk 0 is no longer seeded warm** (`thermal.seed_window`, new; default
  `auto`). The trigger is any run whose output window is shorter than one cycle (chunk 0
  of a chunked run, and also eclipse or flyby windows). In its eq-indirect re-solve the output window is shorter than one cycle,
  the series seed was inert, and the run fell back to the warm anchor with the plain profile.
  Under `auto` that one case now seeds from the equilibration forcing (the forcing the eq
  steps integrate, with the lumped anchor and Kirchhoff profile). Every other run is
  byte-identical to before. On a coupled lunar test mesh, `--chunks 2` minus `--chunks 1`
  T_min went from +5.03 K to -0.44 K mean, and on a 140k-face crater from +4.54 K to +1.69 K;
  the residual is the single-shot run's legacy seed against chunk 0's equilibration seed
  (PTM-802). On the GPU thermal path the 968-face split went from +4.46 K to +2.45 K; the rest
  was chunk 0's sub-cycle mean indirect flux on the equilibration forcing (PTM-806, below). `eq_span` applies the
  equilibration-forcing seed everywhere (opt-in; worse than the default on a run whose output
  window covers a year); `legacy` reproduces earlier outputs byte-for-byte. Also
  `--seed-window` and `PTM_SEED_WINDOW`. `provenance.txt` records `seed_window`,
  `seed_window_source` and `seed_window_effect`.
- New output-cadence advisor. After illumination, the engine measures the Sun's
  angular rate omega and disk theta from the run's own epochs and prints a
  `[cadence-advisor]` block: the output step against the recommended
  0.4 theta/omega (about 25 minutes on the Moon), the 5% peak-flux limit
  0.3 theta/omega, the surface thermal response time, and the largest
  temperature rise a missed flash can cause. After the run it prints an
  estimated upper bound on the fractional overstatement of the cold-trap area at
  110 K and 54 K. A shadow-ray probe decides only whether any face has an
  occluder ("not needed" when none does). Log and `provenance.txt` only;
  silence it with `illumination.cadence_advisor: false`. `ptm cadence <run_dir>`
  recomputes the post-run estimate. Guidance: `docs/thermal-model.md` §2.2,
  "Choosing the output cadence".
- GPU illumination now repeats the CPU path's shadow-ray decisions (PTM-794): ray origins
  and the near-distance cutoff match the CPU, a blocked center ray falls through to
  solar-disk sampling as on the CPU (it used to mean "dark"), the disk sampling pattern is
  the CPU's (the rotation sign was reversed), `PTM_ILLUM_RISE_EPS` is mirrored in the kernel,
  and the horizon cull runs on the GPU path. `PTM_ILLUM_RISE_EPS2` and topological
  self-skip are not mirrored, so with them set direct solar runs on the CPU path, and the run
  says so. Kernel fixes: the traversal stack no longer silently drops subtrees past 30
  entries; the absolute 1e-8 parallel cutoff (which skipped grazing hits on millimeter
  facets) is replaced by a watertight ray-triangle test (Woop, Benthin and Wald 2013, the
  algorithm Embree uses), so grazing rays no longer slip between adjacent triangles; the
  kernel is built without FP contraction and `-cl-fast-relaxed-math` is removed. Errors no longer
  fabricate results: an OpenCL failure used to return every ray occluded, and a failed or
  missing device left rays lit; the run now falls back to the CPU path for the rest of the
  run and says so. The work-group size honors the kernel's own limit. New test-only
  switches: `PTM_BVH_HOST_KERNEL`, `PTM_OPENCL_DEVICE_TYPE`, `PTM_BVH_CL_OPTIONS`
  (`docs/ENVIRONMENT.md`). CPU-path outputs are unchanged.
- GPU illumination is validated on an AMD Radeon Pro W6800X against CPU Embree,
  with the default illumination settings and with the fine-mesh recipe: no
  differences on a 50 km crater core+ring mesh; on a lunar polar core+ring mesh
  (8% of lit face-frames disagreed before the fix) 87 face-frames of 129M differ,
  net energy +0.0000%. Each is one ray flipping at a shadow edge: usually one of
  the 9 disk samples (≤ S·cosθ/9), occasionally the center ray, which moves a
  penumbra face between full Sun and disk sampling (max 174 W/m² on one
  face-frame). `illumination.use_gpu: true` (or
  `--use-gpu-illumination`) prints a NOTICE naming the validated GPU and the new
  self-check; CPU Embree stays the default. The opt-in key
  `illumination.gpu_experimental` and the flag `--gpu-illumination-experimental`,
  added earlier in this release cycle, are retired; if a config or command line
  still has them, the run warns that they are unknown.
- New `ptm gpu-selfcheck` (binary `ptm_gpu_selfcheck`, built in every GPU build):
  checks the GPU shadow-ray tracer on each OpenCL GPU against a float64
  brute-force ground truth on synthetic meshes, in seconds, and reports PASS or
  FAIL. Run it before using GPU illumination on hardware other than the
  W6800X. See `INSTALL.md`.
- The GPU BVH disk cache (`illumination.bvh_cache_file`, on by default) is now
  keyed on its content (vertices, indices, scene offset, build settings); it
  used to be accepted when the triangle COUNT matched, so a different mesh with
  the same face count reused a BVH bounding other triangles (PTM-800). Old
  `BVH1` cache files are rebuilt, and the run says why.
- GPU BVH traversal pads each slab exit distance by 1 + 2γ₃ (Ize 2013), so
  float rounding cannot cull a box a ray enters; no measured change on the
  validation meshes.
- `provenance.txt` records whether GPU illumination served the run, how many
  direct-solar passes it traced, and the first pass that fell back to the CPU
  and why.
- New developer gate `scripts/test_gpu_host_parity.py`: the GPU illumination
  path executed on the host must make the CPU path's lit/dark decisions on a
  small core+ring case (needs a GPU build; uses no GPU).
- `thermal.lateral_scheme` now defaults to `monolithic` (was `balanced`), the
  scheme validated for production by the G15 transient gate. It matters only
  with `thermal.lateral_conduction: true`; set `lateral_scheme: balanced`
  explicitly to keep the split scheme. `lateral_picard` no longer needs the
  scheme named. New checks in `test_config_shape` pin the default.
- Lateral heat conduction is documented as a user option
  (`thermal.lateral_conduction`, CPU SIMD solver only), with every lateral key
  and its default in `docs/configuration.md`. The validated scheme,
  `lateral_scheme: monolithic`, is the default (entry above).
  `docs/KNOWN_LIMITATIONS.md` and `docs/thermal-model.md` no longer say
  lateral conduction is absent.
- New lateral-conduction advisor. Before every thermal run the engine compares
  the mesh's facet sizes with the critical length D*(T) = 8k/(pi h) of the
  configured surface material (Hayne et al., in revision) and prints a
  `[lateral-advisor]` block: RECOMMENDED when facets smaller than D*(50 K)
  cover more than 1% of the surface area, BORDERLINE, or not needed. It never
  changes a setting or an output; `thermal.lateral_advisor: false` silences
  it. The same check runs without the engine as `ptm advise <config>` or
  `ptm.lateral_advisor.advise_config()`. Guidance:
  `docs/thermal-model.md` §2.1, "When to enable lateral conduction".
- A lateral-conduction run now refuses, at startup, a `lateral_scheme` other than
  `balanced` or `monolithic` and `lateral_picard > 0` without `monolithic`
  (previously a misspelled scheme silently ran the balanced split), and refuses
  a mesh in kilometers confirmed by the body radius (the lateral network reads
  vertex coordinates as meters). Valid configurations are unaffected.
- `provenance.txt` records the effective lateral settings (enabled, scheme,
  Picard passes, CG tolerance, activation multiple, zero coupling) and the
  advisor's verdict on every run.
- Europa terrain meshes and view factors are a dataset of their own:
  `python3 scripts/fetch_data.py --dataset europa-terrain` fetches the five
  Galileo SSI SOCET SET terrains (Bland 2022, figshare,
  doi:10.6084/m9.figshare.16904632; CC BY 4.0) into `data/external/sutton_DEMs/`,
  which `configs/europa_*` and the GUI's Europa terrain card read. Its Zenodo
  record id is a placeholder until the record is published; until then the
  fetcher refuses the dataset with a message.
- `fetch_data.py`: a dataset may name its own `base_url` (its own record);
  `--base-url` still overrides it. New test `src/python/tests/test_fetch_data.py`
  (file:// URLs, no network), run in CI.
- `data_manifest.json`: the `moon-topo` attribution named the Moon3D testbed scan;
  it is derived from LRO LOLA (Smith et al. 2010).

## v1.1.1 (2026-10-29): the distribution repository

The same engine as v1.1. This release starts the repository fresh as the PTM3D
distribution, with no prior history; the development repository is archived.

- Removed from the distribution, and kept in a separate research repository:
  campaign orchestration (the polar landing-site, Amundsen and Diviner-validation
  campaigns), study configurations and data, paper-figure and one-off analysis
  scripts, and the Diviner seasonal arrays (59 MB). The Europa terrain grid in
  `configs/`, the cached asteroid shapes and the viewshed tool stay.
- The offline view-factor referee moved from `scripts/ptm12_diag/` to
  `tools/vf_referee/`, and `assemble_vf_haworth.py` is now `assemble_vf.py`.
  Existing CSR sidecars that record the old assembler path still load. A new
  opt-in CMake target, `-DPTM_BUILD_VF_REFEREE=ON`, builds `build/bin/mc_ref_pp`.
- `src/cpp/viewshed/CMakeLists.txt` no longer stops configuration when
  `dark_horse_viewshed_example.cpp` is absent; it warns and skips that target.
- `scripts/run_yearlong_wedge_safe.py` requires `--expect-hole-seal-count`
  (its default was one mesh's value), and `scripts/fetch_horizons_spk.py`
  requires `--spacecraft`, `--start`, `--stop` and `--output`.
- `scripts/ptm_qa.py selftest` drops, and lists, cases whose archived inputs are
  absent, and fails if nothing ran.
- `docs/figures/v2/radiation_exchange_schematic` is re-rendered from the shipped
  synthetic 85N crater (`examples/crater85N_illumination.yaml`).
- The inclusion policy is in `CONTRIBUTING.md`.
- `scripts/test_rethermal_zero_exchange.py` pins `thermal.seed_indirect: off` in both arms. Under
  the v1.1 default (`series`) its two arms started from different spin-up states, and it failed on
  spin-up memory (6.46 mK) rather than on the re-thermal it tests; with the seed pinned it passes at
  1.53e-5 K.

## v1.1: the engine for the Moon3D paper (Hayne et al., in revision)

Results change for radiation-coupled runs; see "Behavior changes".

### Behavior changes

- Spin-up seed for coupled runs. `thermal.seed_indirect` defaults to `series`: in a
  radiation-coupled run the spin-up starts each column from its absorbed flux plus the
  time-varying indirect flux of the first radiation solution, folded over whole diurnal cycles.
  Every face of a coupled run starts closer to its converged state (measured on the shipped
  80°N crater with a 1-cycle spin-up: lit-face mean temperature within 0.36 K of a 40-cycle
  reference instead of 1.79 K). Where it cannot act (radiation off, a warm-started chunk, the
  CPU scalar solver, zero indirect flux, an output window shorter than one diurnal cycle) the run says so and uses
  the previous seed. `provenance.txt` records `seed_indirect`, where it came from, and whether it
  acted. Previous behavior: `thermal.seed_indirect: off` or `--seed-indirect off`.
- CPU coupling passes re-seed with the same series seed (`PTM_CPU_RETHERMAL=seed-series` is
  now the default). On the same crater, never-lit mean temperature moves from 1.75 K to 0.38 K
  from the 40-cycle reference. Previous behavior: `PTM_CPU_RETHERMAL=reseed`.
- Chunked CPU coupled runs write the coupled end-of-window column to the checkpoint, so a
  chunked series matches the single-shot series at the seams.
- View-factor orientation check. A precomputed matrix loaded with `PTM_HMATRIX_VF_REPLACE`
  is refused (exit 3) when it was built for the other side of the surface than the run
  illuminates (`PTM_VF_ORIENT_OVERRIDE=1` proceeds with a warning). The offline referee
  `mc_ref_pp` orients faces by the same rule and refuses a mesh it would have to reverse on more
  than 1% of faces.
- Unknown configuration keys and command-line arguments are reported as warnings instead
  of being ignored silently; `--validate-config` summarizes them ("ACCEPTED WITH WARNINGS").
- `time_steps.txt` labels column 1 correctly: ephemeris seconds past J2000 (the values
  never changed).
- Engine messages no longer carry internal tracking numbers.

### Additions

- `--seed-indirect off|series|mean`; `PTM_TIMING_JSON` (per-stage timing report); view-factor
  statistics sidecars; `assemble_vf_haworth.py --symmetrize-divisor pooled`; an out-of-core view
  factor assembler.
- Apple Silicon (arm64) builds.

### Removed from the repository

- Draft-paper documents, Europa ice-shell research code and documents, internal campaign scripts
  that pointed at private paths, and dated design and audit records.

### Examples and documentation

- Quickstart spin-up raised to 80 cycles; tutorial 2 uses the fine-mesh illumination settings;
  the Moon3D year-long example carries the paper's settings, including Q_b = 0.022 W/m².
- Documentation rewritten for this release; `KNOWN_LIMITATIONS.md` lists only open limitations.

## Fixed in earlier development (after v1.0.3)

Listed so that results from older builds can be read correctly.

- `thermal.geothermal_heat_flow` was ignored on the Moon, Europa and Mercury (all solvers): runs
  used the built-in Q_b whatever the configuration said.
- The CPU SIMD solver used a weaker surface-temperature iteration that could leave a face up to
  tens of kelvin warm at a shadow edge; it now shares the scalar solver's method.
- The CPU solver floored layer thickness at 1 mm, so on fast-rotating, low-inertia bodies its
  top cell was thicker than its own grid.
- The default spin-up seed became the per-face conductive steady profile (`spinup_seed:
  physical`); the 24-cycle clear-sky pre-conditioner of earlier CPU runs is off unless requested.
- Per-step frames are numbered with at least six digits (four before); readers sort by the parsed
  index.
- `equilibration_orbits` means integrated cycles on every solver (earlier GPU runs integrated
  twice the configured value).
- Mercury used its sidereal rotation as the solar day (a factor of three).
- An unrecognized body name silently ran with Europa's properties; it is now refused.
- A 500 K temperature cap applied on one CPU path only; it was removed on all paths.
- The lunar Bond-albedo map was registered 180° off in longitude (the distributed archive is
  corrected); Saturn and Mars primary infrared fluxes were inconsistent with their temperatures.
- The near-field view-factor occlusion test used fixed meter-scale offsets; it is now relative to
  facet size. The GPU radiosity solver now renormalizes view-factor rows above 1 in-kernel.
- Thermal-state checkpoints carry an identity record (grid, mesh, configuration, body) and refuse
  a mismatched restore.
- A `PTM_ILLUM_TRUE_NORMALS` run on a smooth icosphere undid the analytic radial normals; it is
  now declined there.
