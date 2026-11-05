# Direct solar: the solar-disk scheme

`illumination.solar_disk_scheme` selects how the engine turns the finite solar disk into direct
flux on a face. The default, `legacy`, is the rule every earlier engine used and stays byte-identical.
The design and the error budget are in Fermi's solar-disk spec (2026-11-05); this page describes what
the engine does and what was measured.

## The legacy rule and its errors

For each face and frame the CPU batch path (`calculateDirectSolarBatchLegacy`) and the GPU kernel:

1. return 0 if the whole disk is below the face plane;
2. on an open mesh, return 0 if the ray to the disk **center** hits the body sphere `R_occ`;
3. trace the center ray; if it is clear, return the **full** disk, `S cos i`;
4. otherwise count 9 rotated Hammersley samples over the disk.

Three errors follow from it:

- **Center rule.** A disk whose center is visible is never partial. On a straight edge this
  over-delivers 0.105 of `S cos i` on average and up to 0.5 per edge crossing.
- **Sliver dropout.** 9 samples miss thin slivers. A face that sees 0.1-5% of the disk, which is
  exactly the flux range that decides 54 K and 110 K cold-trap membership, often reads 0.
- **Center-only body sphere.** A disk partly above the sphere horizon is zeroed whenever its center is
  below it. On a small open patch the sphere horizon is the local horizon, so a whole patch goes dark
  while half the Sun is up.

## The chord rule (`chord`, `chord_strict`)

Per face and frame, CPU batch path only:

1. **Analytic pre-culls, no rays.** The whole disk below the face plane gives 0. On an open mesh, so
   does a disk lying entirely inside the body sphere's occluding cone (a float64 cone test).
2. **Probes** in a limb-aligned frame, with `y` along the local vertical projected on the sky:
   - center clear: test the bottom limb (`chord_strict` adds the two lower diagonals). If it is
     clear, return the full disk with the face-plane clip (`diskMeanCosAboveHorizon`);
   - center blocked: test the top limb and the two upper diagonals. If all are blocked, return 0.
3. **Chord.** Otherwise the disk is partly visible. Bisect the edge in three columns, at `x = 0` and
   at `x = +/-0.6` of the chord half-width at the center-column edge, with 10 halvings each. Then
   integrate the piecewise-linear edge **exactly** with the linear cosine weight
   (`illumination/solar_disk_chord.h`: piecewise Gauss-Legendre in `x = sin(theta)`, exact to ~1e-12).
4. **Fallback.** If the center column has no edge, or the side columns disagree with it, use 64
   rotated Hammersley samples.

Every ray's "blocked" means: it hits the body sphere, OR it lies at or below the face plane, OR the
terrain blocks it. The terrain test uses the same Embree predicate as the legacy path (rise-epsilon,
self-skip or offset). Frames in which the primary overlaps the solar disk (eclipses) keep the legacy
heavy path, and are counted as such in `provenance.txt`.

**`reference`** is for validation only. It runs a 57-ray screen (9 probes, a 32-point limb ring and
a 16-point ring at r = 0.7). If the screen is all clear it returns the full disk; otherwise it uses a
128 x 128 jittered stratified grid, about 12,900 rays inside the disk. Use it to gate the chord, not
for production.

## What is and is not covered

- **CPU batch direct solar only.** With `illumination.use_gpu: true`, a non-legacy scheme moves direct
  solar off the GPU tracer for the whole run. The log and `provenance.txt` say so (GPU fallback
  record). The GPU kernel is unchanged.
- **The single-face path** (named-location series, non-batch illumination) still uses the legacy rule,
  and the run warns once.
- **The local-horizon cull** (`setHorizonCull`) is a center-ray test and is not applied under a
  non-legacy scheme; the run warns once.
- `provenance.txt` records `solar_disk_scheme`, plus the face-frame counts per class (below-plane,
  sphere-cone, full, zero, chord, both fallbacks, reference grid, eclipse-legacy) and the rays traced.

## Limb darkening (`solar_limb_darkening`)

`illumination.solar_limb_darkening: none | neckel_labs` (default `none`, a uniform disk). It applies
to the non-legacy schemes only; asking for it with `legacy` is fatal, because the legacy rule would
silently ignore it. Limb darkening redistributes brightness across the disk: every flux divides by
the disk-mean intensity, so a fully visible disk above the face plane still delivers exactly
`S cos i`. Only partly visible disks change, and slivers change most, because they see the dark limb.

`neckel_labs` is a bolometric law built from Neckel & Labs (1994, Sol. Phys. 153, 91), Table I. That
table gives 30 continuum 5th-order polynomials `P5(mu) = sum A_k mu^k`, 303-1099 nm; the paper gives
no wavelength-integrated law. We average the coefficients with the weight of the disk-center
spectral intensity, `I_c ~ F_lambda / (F/I)_lambda`, using the AM0 extraterrestrial spectrum
(280-4000 nm) for `F_lambda`:

    A = [0.350358, 1.226714, -1.553505, 1.979649, -1.413788, 0.410569],
    F/I = <I>/I(1) = 0.8293,  I(mu = 0)/I(1) = 0.350

- 26% of the AM0 flux lies outside 303-1099 nm. It takes the nearest tabulated law: the 303 nm law
  below and the 1099 nm law above. The band-only alternative (303-1099 nm) gives F/I = 0.8125 and
  `I(0)/I(1) = 0.304`. Infrared limb darkening is weaker than at 1.1 um, so the adopted law
  slightly OVERSTATES darkening at the limb. **This choice is pending a physics review.**
- For comparison, the textbook linear law with u = 0.6 has F/I = 0.800 and `I(0)/I(1) = 0.400`.
- A sliver with its edge at y = 0.97 (0.31% of the disk) receives 0.599 of its uniform-disk flux.
- The Table I transcription is checked against the paper's own constraints: `sum A_k = 1` and
  `F/I = 2 sum A_k/(k+2)` hold for all 30 rows (one A_2 digit at 365.9 nm was resolved this way).
- Eclipse face-frames keep the legacy heavy path, which samples a uniform disk. The run warns, and
  `provenance.txt` counts them.
