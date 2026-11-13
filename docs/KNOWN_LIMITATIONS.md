# Known limitations

This page lists what the model does not do, or does only approximately, in this release. We
recommend reading it before a long simulation. Problems that have been corrected are listed in
[CHANGELOG.md](../CHANGELOG.md).

## Physics and numerics

1. **Lateral heat conduction is off by default.** Each face is then an independent 1D column.
   On meshes with facets smaller than a few centimeters, small cold traps are warmed by
   conduction from their neighbors, and the default model omits this heating. An example is the
   Moon3D testbed mesh, which has 509,337 faces with a mean facet edge of about 2.5 mm and is
   about 1.3 m across. For a pit with a depth-to-diameter ratio of 0.1 at 85°S, the 110 K cold
   trap disappears below a diameter of about 3–5 cm, depending on the grid. Lateral conduction
   is available as an option (`thermal.lateral_conduction`, CPU solver only). The run prints a
   `[lateral-advisor]` block that reports whether the mesh needs it, as described in
   [thermal-model.md §2.1](thermal-model.md#21-criteria-for-lateral-conduction).

2. **Spin-up length is set by the user.** `thermal.equilibration_orbits` (default 1) sets the
   number of diurnal cycles integrated before the output window. Surface temperatures settle
   in a few cycles. The deep column relaxes over roughly 15 to 23 lunations on the Moon and
   about 10 orbits on Europa. The run prints an estimate of the number of cycles needed and
   a warning when the setting is shorter. The estimate is conservative. For any quantity you
   report, run N and 2N cycles and compare.

3. **Radiation exchange during spin-up is approximate.** In a simulation with radiation
   exchange, the spin-up is forced with an estimate of the indirect (scattered and re-emitted)
   flux. The estimate is taken from the first radiation solution and repeated over whole
   diurnal cycles. On the 80°N crater of the examples, a spin-up of one cycle is then within
   about 0.4 K (median) of a 40-cycle reference, for lit and never-lit faces alike. This
   estimate cannot be applied to a warm-started chunk, with the CPU scalar solver, or when the
   output window is shorter than one diurnal cycle. In these cases the log reports that the
   previous estimate is used, and absolute temperatures are uncertain by about 1 to 2 K. The
   difference between runs with and without radiation exchange is better constrained than
   either absolute value.

4. **The column solver is first order in time.** The surface temperature is solved from the
   subsurface state of the previous substep, so the error decreases linearly with the substep.
   At the default settings the bias in time-mean temperature is at most about 0.1 K. It is
   larger for peak temperatures at coarse output cadence (up to 1.6 K at a 6-hour step on a
   lunar crater). A second-order step, `thermal.column_step: midpoint`, is available for runs
   without radiation exchange or lateral conduction
   ([thermal-model.md §6.3](thermal-model.md#63-substeps-within-an-output-interval)).

5. **CPU and GPU solvers are not identical.** Both interpolate the absorbed flux within each
   output interval, but the GPU solver works in single precision and sets its own substep
   count from a stability limit, with at most 200 substeps per output interval. On a lunar
   crater mesh at a 6-hour output step, GPU temperatures agree with a time-converged CPU
   solution to within 0.1 K on nearly every face, and to about 1 K on a few strongly heated
   faces. Use one solver for all runs in a comparison.

6. **No temperature below the geothermal floor.** Every layer is bounded below by
   T = (Q_b/εσ)^(1/4), which is 24.04 K on the Moon at the default Q_b = 0.018 W/m² (25.28 K at
   0.022 W/m²). A temperature at this value is the bound and not a solution of the heat
   equation. Q_b is constant in time. It has one value per run
   (`thermal.geothermal_heat_flow`), or one per face with a heat-flow mask.

7. **View factors depend on the face-orientation criterion.** By default each face normal is
   flipped radially outward. With `PTM_ILLUM_TRUE_NORMALS=1` the triangle winding is kept, with
   one sign per connected surface. The usual checks do not detect a wrong side, because a
   view-factor matrix traced on the underside of a surface still passes the reciprocity,
   row-sum and nonzero tests. The offline view-factor reference tool applies the same criterion
   and stops with an error on a mesh that it would reverse on more than 1% of faces. A
   precomputed matrix loaded with `PTM_HMATRIX_VF_REPLACE` is rejected if its recorded
   orientation differs from that of the run, and a matrix with no orientation record is treated
   as raw winding. Reuse a view-factor file only with the exact mesh from which it was built.
   On closed, non-convex bodies the radial criterion can choose the wrong side of true
   overhangs.

8. **Fine-mesh illumination corrections are opt-in.** By default the shadow ray starts about
   1 cm above the face and normals are flipped radially outward. On a fine mesh this can
   illuminate millimeter-scale shadowed floors and create false sunlit overhangs. On fine polar
   meshes, enable the full validated set (`PTM_ILLUM_TRUE_NORMALS=1`, `PTM_ILLUM_RISE_EPS=1`,
   `PTM_ILLUM_RISE_KEPS=0.2`, `PTM_ILLUM_RISE_HABS=0`, `PTM_ILLUM_SEAL_HOLES=1`). Do not set
   `PTM_ILLUM_RISE_EPS` alone. The variables are described in [ENVIRONMENT.md](ENVIRONMENT.md).

9. **Sub-facet roughness (beaming) depends on the facet scale.** On polar meshes with facets of
   a meter or larger, beaming is needed to match Diviner. At Haworth it changes the median bias
   of the summer maximum temperature in cold traps from −10.3 K to −0.8 K. On meshes whose
   facets already resolve the roughness (below 5 cm), beaming would count the same physics
   twice, and the run stops with an error.

10. **GPU illumination is validated on one GPU model.** The GPU shadow-ray calculation makes
    the same decisions as the CPU calculation. These are the ray origins, the near-distance
    cutoff, the change from a single center ray to solar-disk sampling, the disk sampling
    pattern, and the `PTM_ILLUM_RISE_EPS` correction of item 8. `PTM_ILLUM_RISE_EPS2` and the
    topological self-skip run on the CPU. The ray-triangle test is watertight, like that of
    Embree. On an AMD Radeon Pro W6800X (macOS) the GPU was compared with CPU Embree on three
    lunar meshes.

    - On a 50 km crater mesh with a core and a surrounding ring, there were no differences in
      4.2 million face-frames, both with the default illumination settings and with the
      fine-mesh settings of item 8.
    - On the 509,337-face Moon3D mesh, 3 of 255,207 lit face-frames differed between lit and
      dark, each by less than 0.07 W/m2.
    - On a lunar polar mesh with a core and a ring, 87 of 129M face-frames differ, with a net
      energy difference of +0.0000%.

    Each difference is one ray that changes at a shadow edge. Usually it is one of the 9 disk
    samples, which changes the flux by at most S·cosθ/9. Occasionally it is the center ray,
    which moves a face in the penumbra between full Sun and disk sampling (at most 174 W/m² on
    one face-frame). These center rays pass within 0.07-0.8 mm of a triangle edge, which is at
    or below the float32 spacing of the vertex coordinates there (0.1-1 mm). In such ties the
    GPU agrees with a float64 trace in 15 of 16 cases, and CPU Embree in robust mode counts
    them as hits.

    Two limits follow. First, the validation covers that GPU and lunar meshes only. Another
    GPU, driver or OpenCL compiler can execute the kernel incorrectly, as the macOS OpenCL CPU
    device does. Run `ptm gpu-selfcheck` (or `build/bin/ptm_gpu_selfcheck`) before using any
    other hardware for a long simulation. It checks every GPU it finds against a float64
    reference in seconds and reports PASS or FAIL. For developers,
    `scripts/test_gpu_host_parity.py` tests that the host logic of the GPU calculation equals
    that of the CPU calculation. Second, CPU Embree remains the default. On the tested
    workstation it was also faster, with 70.5 s of ray tracing on the core and ring mesh
    against 118.8 s on the GPU. The GPU time includes the 500 ms pause between dispatches
    (`illumination.gpu_inter_batch_pause_ms`) that the validation runs used. The default for
    that pause is 0. A run with `use_gpu: true` prints a NOTICE that names the validated GPU and
    the self-check.

11. **Spin-up seed and seasonal memory.** A spin-up of a few cycles before a window, from any
    seed, converges toward the column that the forcing of that window would give if it
    repeated. At polar sites the deep column instead retains the seasonal cycle, because the
    annual thermal wave reaches 0.1-0.2 m. A single-window run is therefore biased by season.
    Windows in the dark season are several K too cold on lit faces and up to ~25 K too cold on
    seasonally shadowed faces, and windows in the bright season are too warm. On a lunar crater
    test mesh, compared with a continuous real-ephemeris reference on the same depth grid, the
    bright-season errors are +2.0 K in T_min on lit faces, +0.3 K in T_max in cold traps, and
    +0.8 K in the T_max < 54 K class. Year-long windows are not affected, because their seed is
    the annual mean. `thermal.seed_seasonal` corrects the seed for this effect (see
    configuration.md). Its accuracy and tested range are as follows.
    - The correction is accurate through the partial cancellation of two biases and is not
      exact. The quasi-static surface potential computed for each phase from the lumped anchor
      is 1.1-2.8 K cold on lit phases and about 0 on dark phases. The seasonal annual-mean
      level that it feeds is warm by an offset that is uniform with depth, about +0.8 to
      +3.6 K at the column bottom when the exact per-phase potential is used instead. The
      source of the warm offset is not identified. It is not the depth-grid integration, the
      column linearization, the run protocol after the seed, or year-to-year drift of the
      forcing.
    - In comparisons with continuous real-ephemeris references at a polar site, the corrected
      seed was never worse on 72 of 72 face-windows, and the window-mean |error| fell from
      several K to 0.1-0.4 K. On a lunar crater test mesh in a bright-season window, about
      0.26% of lit faces are up to ~1 K worse than with the uncorrected seed, while the
      aggregate T_min error falls from 1.6 K to 0.3 K. On the same mesh in a dark-season
      window, when nearly every face is a cold trap, the aggregate error halves (0.62 K to
      0.29 K), but the correction is too warm on cold traps. In the T_max < 54 K class the
      median error is +0.28 K, where the uncorrected seed gives -0.64 K, and 6% of faces are
      worse than with the uncorrected seed by more than 0.1 K, up to 1.3 K. This is the
      unexplained warm offset described above, in a case where the cold bias of the lit phases
      does not cancel it.
    - The correction has been tested for the Moon, for single windows shorter than half an
      eclipse year, at about 85 degrees latitude, with the CPU SIMD thermal solver. Equatorward
      the seasonal bias itself is small. On the test mesh the single-window cold-trap bias is
      about 0.13 K at 80 degrees, 0.04 K at 70 and about 0 at 60. The correction is not applied
      to a face where it is under 0.5 K. We expect this to leave lower latitudes unchanged, but
      this has not been shown against a reference.
    - In chunked runs, the equilibration of the first chunk uses by default the indirect flux
      resolved by diurnal phase (`thermal.eq_indirect_subcycle: phase`) with both the CPU and
      GPU solvers. The GPU solver seeds each coupling iteration again
      (`thermal.gpu_rethermal: seed_series`), so the seasonal correction is applied at every
      iteration, as with the CPU solver. On a 4990-face crater window with
      `seed_seasonal: auto`, the mean GPU minus CPU difference in T_min is -0.10 K, compared
      with +2.47 K without this reseeding. Neither setting has been compared together with the
      seasonal seed against a continuous reference.

12. **Peak temperatures sampled at the output times miss brief illumination.**
    `temperature_max` and `direct_solar_max` are taken over the output epochs, so a face lit
    for less than one output step can be recorded as colder than it is. Illumination at grazing
    incidence can be arbitrarily brief, so no step removes the effect. On a decimated Moon3D
    crater mesh over one lunation, this overstated the 110 K cold-trap area by 0.7% (relative)
    at a 2-hour step and by 2.5% at a 5.9-hour step. The run prints a `[cadence-advisor]` block
    with the recommended step (0.4 θ/ω, about 25 minutes on the Moon) and an estimate of the
    overstatement after the run, as described in
    [thermal-model.md §2.2](thermal-model.md#22-choosing-the-output-cadence).

13. **Results computed without the body test for the primary disk include light from below
    the limb.** On an open patch of terrain, rays toward the primary body (Earthshine or
    planetshine) are tested against the body sphere. Without this test, a face tilted toward a
    primary that has set below the limb still receives its infrared and reflected light through
    the body. This occurs only where the primary can be at or below the limb, which is near the
    poles and, on a synchronously rotating satellite, near the longitudes about 90° from the
    sub-primary point at any latitude. Set `illumination.primary_disk_body_cull: false` only to
    reproduce a result computed without the test. The flux-cube pass traces its rays without
    the sphere. The size of the effect in three cases is as follows.
    - On the Moon3D 85°S crater, 19% of the face epochs lit by Earth had the whole disk behind
      the body and 12% had part of it. The excess is at most 0.07 W/m² of infrared and
      0.08 W/m² of reflected Earthlight per face. Over one lunation without radiation exchange,
      the test left the 54 K and 110 K cold-trap areas unchanged. It lowered T_max by more
      than 1 K on 2.3 cm² of the 1.3 m crater, and by up to 5.9 K on 13 faces that were already
      cold (about 30 K) and tilted toward Earth, where the absorbed flux removed was
      0.075–0.139 W/m².
    - On the Haworth crater over a year, 1.93% of the Earth-lit face epochs had the disk behind
      the body and 2.88% had part of it. The largest absorbed Earthshine on such a face epoch,
      which is an upper bound on the excess, was 0.029 W/m², outside the ±15 km window of the
      Diviner comparison. Inside the window, on the cold-trap faces (T_max < 110 K), it was at
      most 4.7e-4 W/m².
    - On the Europa Pwyll regional patch placed on the equator 88° from the sub-Jovian point,
      Jupiter's disk reaches the body on 97% of the Jupiter-lit face epochs, but the change is
      small. It affects 621 of 28,300 faces and 0.03% of the Jupiter light on the patch, and
      changes temperatures by at most 0.0015 K.

    Runs in which the primary stays above the limb, and runs on closed meshes, are unaffected.

## Configuration

- Very small generic-body properties are lost. A `generic_body.surface` value below 5e-7 is
  rounded to zero and replaced by a placeholder without a warning. Keep such values at or above
  5e-7, or set them with the `*_override` keys, which keep full precision.
- Plain property keys are ignored on the Moon, Mercury and Europa, and the run prints a
  warning. Use the `thermophysical.*_override` keys. A per-face albedo map takes precedence
  over `albedo_override` for absorbed sunlight.
- Unknown keys and command-line flags are reported but do not stop the run. Check the
  warnings, or run `ptm_thermal_sim --config FILE --validate-config` first.
- `mesh.scale_factor` is not supported, and the run stops with an error when it is set.
  Rescale the mesh file.

## Platforms

- Of the GPU versions, only OpenCL on macOS with AMD GPUs has been tested. The HIP and CUDA
  versions are untested. Long GPU runs on macOS with AMD hardware can hang the machine. The
  precautions are described in [dev/gpu-wedge-guardrails.md](dev/gpu-wedge-guardrails.md).
- On Windows, use WSL2 and the Linux instructions.
