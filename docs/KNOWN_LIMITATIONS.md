# Known limitations

What the model does not do, or does only approximately, in this release. Read this before a
production run. Fixed problems are listed in [CHANGELOG.md](../CHANGELOG.md), not here.

## Physics and numerics

1. **Lateral heat conduction is off by default.** Each face is then an independent 1D column.
   Where facets are smaller than a few centimeters (the Moon3D testbed mesh: 509,337 faces, mean
   facet edge about 2.5 mm, about 1.3 m across), small cold traps are warmed by their neighbors in
   reality but not in the default model; for a pit with depth-to-diameter 0.1 at 85°S, the 110 K
   cold trap disappears below about 3–5 cm diameter (grid-dependent). Lateral conduction is available as an option
   (`thermal.lateral_conduction`, CPU solver only), and the run prints a `[lateral-advisor]` block
   saying whether this mesh needs it; see
   [thermal-model.md §2.1](thermal-model.md#21-when-to-enable-lateral-conduction).

2. **Spin-up length is set by the user.** `thermal.equilibration_orbits` (default 1) sets how many
   diurnal cycles are integrated before the output window. Surface temperatures settle in a few
   cycles; the deep column relaxes over roughly 15 to 23 lunations on the Moon and about 10
   orbits on Europa. The run prints the number of cycles it estimates you need and warns when
   the setting is shorter (the estimate is conservative); for any quantity you report, run N and
   2N cycles and compare.

3. **Radiation exchange during spin-up is approximate.** In a coupled run the spin-up is forced
   with an estimate of the indirect (scattered and re-emitted) flux taken from the first
   radiation solution, folded over whole diurnal cycles. On the shipped 80°N crater a one-cycle spin-up
   then lies within about 0.4 K (median) of a 40-cycle reference, for lit and never-lit faces
   alike. Where this seed cannot act (a warm-started chunk, the CPU scalar solver, an output
   window shorter than one diurnal cycle) the run says so and uses the previous seed, and
   absolute temperatures are uncertain by about 1 to 2 K. The difference between coupled and
   uncoupled runs is better constrained than either absolute value.

4. **The column solver is first order in time.** The surface temperature is solved from the
   previous substep's subsurface state, so the error falls linearly with the substep. At the
   default settings the bias in time-mean temperature is at most about 0.1 K; it is larger for
   peak temperatures at coarse output cadence.

5. **CPU and GPU solvers are not equivalent.** The GPU holds each output interval's absorbed
   flux constant across its substeps, which raises lit-face peak temperatures by about 0.4 K
   (median) and a few kelvin in the tail at daily cadence, and more at coarser cadence; it also
   caps the substep count at 200. Use the CPU solver for accuracy-critical work, and never mix
   backends within one comparison.

6. **No temperature below the geothermal floor.** Every layer is bounded below by
   T = (Q_b/εσ)^(1/4): 24.04 K on the Moon at the default Q_b = 0.018 W/m² (25.28 K at
   0.022 W/m²). A value at the floor is the bound, not a solution. Q_b is constant in time:
   one value per run (`thermal.geothermal_heat_flow`), or one per face with a heat-flow mask.

7. **View factors follow the engine's face-orientation rule.** By default each face normal is
   flipped radially outward; with `PTM_ILLUM_TRUE_NORMALS=1` the triangle winding is kept, with
   one sign per connected surface. A wrong side is not caught by the usual checks: a view-factor
   matrix traced on the underside still passes reciprocity, row-sum and nonzero tests. The
   offline referee applies the same rule and refuses a mesh it would reverse on more than 1% of
   faces. A precomputed matrix loaded with `PTM_HMATRIX_VF_REPLACE` is refused if its recorded
   orientation differs from the run's; a matrix with no orientation record is treated as raw
   winding. Reuse a view-factor file only with the exact mesh it was built from. On closed,
   non-convex bodies the radial rule can choose the wrong side of true overhangs.

8. **Fine-mesh illumination corrections are opt-in.** By default the shadow ray starts about
   1 cm above the face and normals are flipped radially outward, which can light
   millimeter-scale shadowed floors and create false sunlit overhangs. On fine polar meshes,
   enable the full validated set (`PTM_ILLUM_TRUE_NORMALS=1`, `PTM_ILLUM_RISE_EPS=1`,
   `PTM_ILLUM_RISE_KEPS=0.2`, `PTM_ILLUM_RISE_HABS=0`, `PTM_ILLUM_SEAL_HOLES=1`), never
   `PTM_ILLUM_RISE_EPS` alone; see [ENVIRONMENT.md](ENVIRONMENT.md).

9. **Sub-facet roughness (beaming) is a scale choice.** On meter-scale and coarser polar meshes
   beaming is needed to match Diviner (at Haworth it moves the cold-trap summer-maximum bias from
   −10.3 K to −0.8 K median); on meshes whose facets already resolve the roughness (below 5 cm)
   it would count the same physics twice, and the engine refuses it.

10. **GPU illumination is validated on one GPU model (PTM-794).** The GPU shadow-ray path
    repeats the CPU path's decisions (ray origins, near-distance cutoff, center-ray fallback to
    solar-disk sampling, disk sampling pattern, and the `PTM_ILLUM_RISE_EPS` correction of
    item 8; `PTM_ILLUM_RISE_EPS2` and topological self-skip run on the CPU) and uses a
    watertight ray-triangle test like Embree's. On an AMD Radeon Pro W6800X (macOS) it was
    compared with CPU Embree on three lunar meshes. A 50 km crater core+ring mesh: 0
    differences of any kind in 4.2 million face-frames, both with the default illumination
    settings and with the fine-mesh recipe of item 8. The 509,337-face Moon3D mesh: 3 lit/dark
    disagreements in 255,207 lit face-frames, each below 0.07 W/m2. A lunar polar core+ring
    mesh, where the build before the fix disagreed on 8% of lit face-frames: 87 face-frames of
    129M differ, net energy +0.0000%. Each is one ray flipping at a shadow edge: usually one of
    the 9 disk samples (at most S·cosθ/9), occasionally the center ray, which moves a penumbra
    face between full Sun and disk sampling (max 174 W/m² on one face-frame). Traced one by one,
    those center rays pass within 0.07-0.8 mm of a triangle edge, at or below the float32
    spacing of the vertex coordinates there (0.1-1 mm): ties, which the GPU resolves like a
    float64 trace in 15 of 16 cases and CPU Embree (robust mode) counts as hits. Two limits
    follow. The validation covers that GPU and lunar meshes only: another GPU, driver or
    OpenCL compiler can mis-execute the kernel (the macOS OpenCL CPU device does), so run
    `ptm gpu-selfcheck` (or `build/bin/ptm_gpu_selfcheck`) before production use on any other
    hardware; it checks every GPU it finds against a float64 ground truth in seconds and
    reports PASS or FAIL (`scripts/test_gpu_host_parity.py` is the developer gate that keeps the
    GPU path's host logic equal to the CPU path's). And CPU Embree stays the default: on the tested workstation it was
    also faster (70.5 s against 118.8 s of ray tracing on the core+ring mesh), with the
    500 ms wedge-safety pause between dispatches (`illumination.gpu_inter_batch_pause_ms`)
    that the validation runs used; the engine default for that pause is 0. A run with
    `use_gpu: true` prints a NOTICE naming the validated GPU and the self-check.

11. **Spin-up seed and seasonal memory.** Every spin-up of a few cycles before a window, from any
    seed, converges toward the column that window's own forcing would give if it repeated. At
    polar sites the deep column instead remembers the seasons (the annual thermal wave reaches
    0.1-0.2 m), so a single-window run is biased by season: dark-season windows several K too cold
    on lit faces and up to ~25 K too cold on seasonally shadowed faces; bright-season windows
    warm (lunar crater test mesh, against a continuous real-ephemeris reference on the SAME depth
    grid: lit T_min +2.0 K, cold-trap T_max +0.3 K, T_max < 54 K class +0.8 K). Year-long windows
    are not affected (they seed on the annual mean).
    `thermal.seed_seasonal` corrects the seed for this (see configuration.md). What it does and
    does not establish:
    - It is accurate by **partial cancellation of two biases**, not exactly: the per-phase
      quasi-static surface potential from the lumped anchor runs 1.1-2.8 K COLD on lit phases (and
      about 0 on dark phases), while the seasonal annual-mean level it feeds runs WARM by a
      depth-uniform offset (about +0.8 to +3.6 K at the column bottom when the exact per-phase
      potential is used instead). The source of the warm offset is not identified; it is not the
      depth-grid integration, the column linearization, the run protocol after the seed, or
      year-to-year forcing drift.
    - Measured against continuous real-ephemeris references: never worse on 72 of 72
      face-windows at a polar site (window-mean |error| from several K to 0.1-0.4 K). On a lunar
      crater test mesh, bright-season window: about 0.26% of lit faces up to ~1 K worse than the
      old seed (aggregate T_min error 1.6 K -> 0.3 K). Same mesh, DARK-season window (nearly every
      face a cold trap): aggregate error halves (0.62 K -> 0.29 K) but the correction OVERSHOOTS
      warm on cold traps (T_max < 54 K class: median +0.28 K, where the old seed is -0.64 K), and
      6% of faces are worse than the old seed by more than 0.1 K, up to 1.3 K. This is the
      unexplained warm offset below, seen where the cold lit-phase bias has nothing to cancel.
    - Validated scope: the Moon, single windows shorter than half an eclipse year, at about
      85 degrees latitude (every reference so far); CPU SIMD thermal path. Equatorward the
      seasonal bias itself is small (single-window cold-trap bias about 0.13 K at 80 degrees,
      0.04 K at 70, about 0 at 60 on the test mesh), and the per-face gate leaves a face alone
      where the correction is under 0.5 K; that this keeps lower latitudes unchanged is
      expected, not yet shown against a reference.
    - Chunked runs: chunk 0 of a GPU-thermal chunked run still carries a sub-cycle mean indirect
      flux on its equilibration forcing (PTM-806), independent of the seed.

## Configuration

- Very small generic-body properties are lost. A `generic_body.surface` value below 5e-7 is
  rounded to zero and replaced by a placeholder without warning. Keep such values at or above
  5e-7, or set them with the `*_override` keys, which keep full precision.
- Plain property keys are ignored on the Moon, Mercury and Europa (the run warns); use the
  `thermophysical.*_override` keys. A per-face albedo map overrides `albedo_override` for
  absorbed sunlight.
- Unknown keys and command-line flags are reported but do not stop the run. Check the
  warnings, or run `ptm_thermal_sim --config FILE --validate-config` first.
- `mesh.scale_factor` does not work; the engine refuses such a mesh. Rescale the mesh file.

## Platforms

- Of the GPU backends, only OpenCL on macOS with AMD GPUs has been exercised; the HIP and CUDA
  backends are untested. Long GPU runs on macOS with AMD hardware can hang the machine; see
  [dev/gpu-wedge-guardrails.md](dev/gpu-wedge-guardrails.md).
- On Windows, use WSL2 and the Linux instructions.
