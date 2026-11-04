# Changelog

## Unreleased

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
