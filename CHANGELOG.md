# Changelog

## Unreleased

- GPU illumination is now opt-in experimental. `illumination.use_gpu: true`
  (or `--use-gpu-illumination`) is refused at startup unless
  `illumination.gpu_experimental: true` (or `--gpu-illumination-experimental`)
  is also set, and an opted-in run prints an EXPERIMENTAL warning. The GPU
  ray tracer is not yet equivalent to the CPU (Embree) path (PTM-794): on a
  lunar polar mesh the two disagreed on the lit/dark state of 8% of CPU-lit
  face-frames, at grazing Sun and shadow edges. CPU-path runs are unchanged.
  See `docs/KNOWN_LIMITATIONS.md`.
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
