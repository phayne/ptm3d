# Integration of the direct solar flux over the solar disk

`illumination.solar_disk_scheme` selects how the direct solar flux on a face is integrated over
the finite solar disk. The default is `legacy`. The `chord` scheme is more accurate for a partly
visible disk, and `reference` is a dense sampling used for validation.

## The legacy scheme and its errors

The legacy scheme is used by the CPU calculation (`calculateDirectSolarBatchLegacy`) and by the
GPU kernel. For each face and frame it proceeds as follows.

1. If the whole disk is below the plane of the face, the flux is 0.
2. On an open mesh, if the ray to the center of the disk intersects the body sphere `R_occ`, the
   flux is 0.
3. The ray to the center of the disk is traced. If it is clear, the full disk is taken to be
   visible and the flux is `S cos i`.
4. Otherwise the visible fraction is estimated with 9 rotated Hammersley samples over the disk.

This scheme has three sources of error.

- **Center test.** A disk whose center is visible is never treated as partly hidden. At a
  straight edge this overestimates the flux by 0.105 of `S cos i` on average, and by up to 0.5 at
  each crossing of the edge.
- **Loss of thin slivers.** Nine samples can miss a thin sliver of the disk. A face that sees
  0.1-5% of the disk often receives 0, and this is the range of flux that determines whether a
  face is a cold trap at 54 K or 110 K.
- **Body sphere tested at the center only.** A disk partly above the horizon of the body sphere
  is set to zero whenever its center is below that horizon. On a small open patch the sphere
  horizon is the local horizon, so a whole patch is dark while half of the Sun is above it.

## The chord scheme (`chord`)

The chord scheme locates the edge of the visible part of the disk and integrates the visible
area analytically. It is available in the CPU calculation only. For each face and frame:

1. **Analytic tests without rays.** A disk wholly below the plane of the face gives 0. On an open
   mesh, so does a disk lying entirely inside the occluding cone of the body sphere (a float64
   cone test).
2. **Visibility test with 25 rays.** The rays are placed in a frame aligned with the limb, with
   `y` along the local vertical projected on the sky: one at the center, 16 on the limb
   (r = 0.999) and 8 at r = 0.6. If all are clear, the full disk is returned, clipped by the
   plane of the face (`diskMeanCosAboveHorizon`). If all are blocked, the flux is 0. Sparser
   arrangements are not sufficient (see the validation table below). Nine rays at r = 0.85 miss
   slivers at the limb (1.6% of them), and nine at r = 0.999 leave errors of up to 228 W/m2. An
   arrangement that tests only the center with the top or the bottom of the limb misses
   occluders at the side of the disk, with errors of up to 566 W/m2.
3. **Chord.** Otherwise the disk is partly visible. The edge of the visible region is located by
   bisection in three columns, at `x = 0` and at `x = +/-0.6` of the chord half-width at the edge
   found in the center column, with 10 halvings each. The piecewise-linear edge is then
   integrated exactly with the linear cosine weight (`illumination/solar_disk_chord.h`: piecewise
   Gauss-Legendre in `x = sin(theta)`, exact to ~1e-12).
4. **Verification of the edge.** 128 points spread over the disk are traced, omitting points
   within 0.005 disk radii of the fitted edge and points at which the face does not see the Sun.
   If any traced point lies on the wrong side of the fitted edge, if the center column has no
   edge, or if the side columns are inconsistent with it, the flux is computed from 1024 rotated
   Hammersley samples instead.

A ray is counted as blocked if it intersects the body sphere, if it lies at or below the plane of
the face, or if the terrain blocks it. The terrain test is the same as in the legacy scheme
(rise-epsilon, self-skip or offset). Frames in which the primary overlaps the solar disk
(eclipses) use the dense sampling of the legacy scheme, and are counted in `provenance.txt`.

**`reference`** is intended for validation only. It tests visibility with 57 rays (9 interior
points, a 32-point ring on the limb and a 16-point ring at r = 0.7). If all are clear, it returns
the full disk. Otherwise it uses a 128 x 128 jittered stratified grid, about 12,900 rays inside
the disk. It is the standard against which the chord scheme is tested, and it is too slow for
routine use.

## Ray packets in the chord scheme

The scheme above determines which rays are traced. The way in which they are passed to the ray
tracer affects the computation time but not the result. With the default rise-epsilon ray origin,
every ray of a face starts at the face centroid, so the rays are traced as 16-ray Embree packets
(`rtcOccluded16`, with the same ray and filter inputs in each lane as a single-ray query). The
illumination library is compiled with `-ffast-math`, so identical results are a tested property
of this arrangement and are not guaranteed for any rearrangement of the arithmetic.

- The 25 visibility rays are traced 16 per packet, stopping after the first packet that contains
  both a clear and a blocked ray, since the disk is then partly visible whatever the remaining
  rays give.
- The edge-verification points are traced together, with the points predicted visible and those
  predicted blocked in separate packets. The result is the first wrong-side point in the
  original order of the points.
- The fallback samples and the column bisection are traced one ray at a time. Each halving of
  the bisection depends on the previous one. The fallback takes about 1% of the time, and tracing
  it in packets changes the rounding of its sample directions under `-ffast-math`, which changes
  the result of a grazing sample.

`provenance.txt` records the number of rays that the scheme traces one at a time, independent of
the packet width. The following diagnostics each write a notice to standard error:
`PTM_SOLAR_CHORD_PACKETS=0|4|8|16` (packet width, 0 = one query per ray),
`PTM_SOLAR_CHORD_PROFILE=1` (rays and thread time per stage for each frame, which slows the run)
or `=2` (wall time and class counts per frame only), and `PTM_SOLAR_CHORD_DUMP=<dir>` (the fitted
edge of every partly visible face-frame).

## Optional approximations (off by default)

Three options reduce the cost of the chord scheme and change its results. Each prints a NOTICE
and is recorded in `provenance.txt`. None is a default, and only settings that meet the
validation criteria below should be used.

- `illumination.solar_disk_chord_check_adaptive_n: N` (default 0 = off). The 25 visibility rays
  are compared with the fitted edge, using the same margin and face-plane conditions as the edge
  verification. If every ray agrees, the fit is verified with N points instead of 128, and
  otherwise with the full 128.
- `illumination.solar_disk_chord_fallback_n: N` (default 1024, range 8-4096) sets the number of
  fallback samples.
- `illumination.solar_disk_chord_temporal: delta` (disk radii, default 0 = off). The edge search
  in each column starts from a bracket of +/-delta around the edge extrapolated linearly from the
  fits of the face in the two previous calls. The bracket is used only when the Sun moved by the
  same angle (within 5%, under 2 degrees) in both steps. If the ends of the bracket do not
  straddle an edge of the expected orientation, the full bisection is performed. The edge
  verification is unchanged. The result depends on the sequence of frames, so a chunked or
  restarted run differs from a single-shot run in the first frames after each start.

All three stop the run with an error if the scheme is not `chord`.

The table summarizes tests of these options on a lunar polar site (24 hourly frames of its
sunniest output day) and on the validation set described below.

| option | direct-solar time at the polar site, relative to the default | validation set (Moon3D, 140 frames) | recommendation |
|---|---|---|---|
| `check_adaptive_n` 16 / 32 / 64 | not timed | fails: max dF 45 / 44 / 42 W/m2 (limit 20) | do not use |
| `fallback_n` 512 | fallback is about 1% of the time | passes (max 19.5, sliver loss 4e-5, worst class 0.27%) | no useful gain |
| `fallback_n` 256 | as above | fails: max dF 26.8 | do not use |
| `temporal` 0.04 | -5% (33.0 vs 34.6-35.3 s) | not tested: at the 5.9-h cadence of the set the Sun moves more than 2 degrees per step, so no frame uses the option | at the polar site, max dF 0.84 W/m2 relative to the default, none > 2 |

## Scope and limitations

- **CPU direct solar only.** With `illumination.use_gpu: true`, a scheme other than `legacy`
  moves the direct solar calculation from the GPU to the CPU for the whole run. The log and
  `provenance.txt` record this. The GPU kernel implements the legacy scheme only.
- **Single-face calculation.** The series computed for named locations, and illumination
  computed one face at a time, use the legacy scheme. A warning is issued once.
- **Local-horizon cull.** `setHorizonCull` is a test on the center ray and is not applied with a
  scheme other than `legacy`. A warning is issued once.
- `provenance.txt` records `solar_disk_scheme`, the face-frame counts per class (below-plane,
  sphere-cone, full, zero, chord, both fallbacks, reference grid, eclipse-legacy) and the rays
  traced.

## Limb darkening (`solar_limb_darkening`)

`illumination.solar_limb_darkening: none | neckel_labs` (default `none`, a uniform disk) applies
to the `chord` and `reference` schemes only. Requesting it with `legacy` stops the run with an
error, because the legacy scheme has no limb darkening. Limb darkening redistributes brightness
across the disk. Every flux is divided by the disk-mean intensity, so a fully visible disk above
the plane of the face still gives exactly `S cos i`. Only partly visible disks change, and thin
slivers change most, because they see the dark limb.

`neckel_labs` is a bolometric law constructed from Table I of Neckel & Labs (1994, Sol. Phys.
153, 91). That table gives 30 continuum 5th-order polynomials `P5(mu) = sum A_k mu^k` from 303 to
1099 nm, and the paper gives no wavelength-integrated law. We average the coefficients with the
weight of the disk-center spectral intensity, `I_c ~ F_lambda / (F/I)_lambda`, using the AM0
extraterrestrial spectrum (280-4000 nm) for `F_lambda`:

    A = [0.350358, 1.226714, -1.553505, 1.979649, -1.413788, 0.410569],
    F/I = <I>/I(1) = 0.8293,  I(mu = 0)/I(1) = 0.350

- 26% of the AM0 flux lies outside 303-1099 nm. This flux is assigned the nearest tabulated law,
  that of 303 nm below the range and that of 1099 nm above it. Restricting the average to
  303-1099 nm instead gives F/I = 0.8125 and `I(0)/I(1) = 0.304`. Infrared limb darkening is
  weaker than at 1.1 um, so the adopted law slightly overstates the darkening at the limb. This
  treatment of the infrared is provisional.
- For comparison, the linear law with u = 0.6 has F/I = 0.800 and `I(0)/I(1) = 0.400`.
- A sliver with its edge at y = 0.97 (0.31% of the disk) receives 0.599 of its uniform-disk flux.
- The transcription of Table I satisfies the constraints given in the paper, `sum A_k = 1` and
  `F/I = 2 sum A_k/(k+2)`, for all 30 rows. One digit of A_2 at 365.9 nm was determined from
  these constraints.
- Face-frames in eclipse use the dense sampling of the legacy scheme, which assumes a uniform
  disk. A warning is issued, and `provenance.txt` counts them.

## Validation against the reference scheme

The chord scheme was compared with `reference` on the Moon3D mesh over 140 frames. The criteria
for acceptance are:

- total energy within +/-0.1%
- slivers with F_ref >= 2 W/m2 returned as 0 in at most 0.1% of cases
- max |dF| <= 20 W/m2
- at most 0.5% of the faces in each class see any |dF| > 2 W/m2. The classes are the 54 K and
  110 K cold traps (CT54, CT110) and the faces within +/-10 K of 54 and 110 K, according to
  T_max in the first window of the coupled Moon3D simulation.

The 99th percentile of the equivalent temperature error (keq p99) is reported but is not a
criterion.

| visibility rays / verification points / fallback samples | energy | sliver loss | max dF (W/m2) | worst class fraction | keq p99 | result |
|---|---|---|---|---|---|---|
| center plus top or bottom limb / none / 64 | +7e-5 | 0.28% | 566 | 0.64% | 7.5 K | fails |
| 9 at 0.85 / none / 64 | +1.7e-3 | 1.57% | 228 | 1.15% | 21.5 K | fails |
| 9 at 0.999 / none / 64 | +4e-5 | 0.28% | 228 | 0.64% | 7.3 K | fails |
| 25 / 16 at 0.05 / 256 | +3e-5 | 0.025% | 89 | 0.39% | 3.0 K | fails |
| 25 / 64 at 0.01 / 1024 | +2e-5 | 0.001% | 42 | 0.18% | 1.3 K | fails |
| **25 / 128 at 0.005 / 1024 (default)** | **+1e-5** | **0.001%** | **19.5** | **0.14%** | **1.2 K** | **passes** |

Energy is the sum over the face-frames, unweighted by face area.
