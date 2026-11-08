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

## The chord rule (`chord`)

Per face and frame, CPU batch path only:

1. **Analytic pre-culls, no rays.** The whole disk below the face plane gives 0. On an open mesh, so
   does a disk lying entirely inside the body sphere's occluding cone (a float64 cone test).
2. **A 25-probe screen** in a limb-aligned frame, with `y` along the local vertical projected on the
   sky: the center, 16 probes on the limb (r = 0.999) and 8 at r = 0.6. If all are clear, return
   the full disk with the face-plane clip (`diskMeanCosAboveHorizon`). If all are blocked, return 0.
   Screens measured against the reference (Moon3D R1 mesh, 140 frames; restated G2, see below):
   - the first version's 2/4 limb ladder (center plus the bottom limb, or center plus the top limb
     and upper diagonals) missed side occluders, up to 566 W/m2 (PTM-838);
   - 9 probes at r = 0.85 miss limb slivers (dropout 1.6%);
   - 9 at 0.999 leave a 228 W/m2 edge-fit tail.
3. **Chord.** Otherwise the disk is partly visible. Bisect the edge in three columns, at `x = 0` and
   at `x = +/-0.6` of the chord half-width at the center-column edge, with 10 halvings each. Then
   integrate the piecewise-linear edge **exactly** with the linear cosine weight
   (`illumination/solar_disk_chord.h`: piecewise Gauss-Legendre in `x = sin(theta)`, exact to ~1e-12).
4. **Edge check and fallback.** Trace 128 points spread over the disk. Skip points within 0.005
   disk radii of the fitted edge, and points where the face does not see the Sun. If any traced
   point lands on the wrong side of the fit, or if the center column has no edge, or if the side
   columns disagree with it, use 1024 rotated Hammersley samples instead.

Every ray's "blocked" means: it hits the body sphere, OR it lies at or below the face plane, OR the
terrain blocks it. The terrain test uses the same Embree predicate as the legacy path (rise-epsilon,
self-skip or offset). Frames in which the primary overlaps the solar disk (eclipses) keep the legacy
heavy path, and are counted as such in `provenance.txt`.

**`reference`** is for validation only. It runs a 57-ray screen (9 probes, a 32-point limb ring and
a 16-point ring at r = 0.7). If the screen is all clear it returns the full disk; otherwise it uses a
128 x 128 jittered stratified grid, about 12,900 rays inside the disk. Use it to gate the chord, not
for production.

## How the chord is traced (cost, results unchanged)

The rule above fixes which rays decide the answer; how they are handed to Embree does not change it.
On the default rise-epsilon path every ray of a face starts at the face centroid, so the engine traces
them as 16-ray Embree packets (`rtcOccluded16`, the same float ray and filter inputs per lane as the
one-ray query). The illumination library is compiled with `-ffast-math`, so "identical" here is a
measured property of this build (gates below), not a guarantee for any restructuring of the arithmetic:

- the screen, 16 probes per packet, stopping after the first packet that has seen both a clear and a
  blocked probe (the face is then partial whatever the rest say);
- the edge check, all its points at once, the points predicted visible and those predicted blocked in
  separate packets; the answer is the first wrong-side point in the check's own order, as before;
- the fallback samples and the column bisection one ray at a time (the bisection because each halving
  depends on the last; the fallback because it is about 1% of the time, and a packet form changed the
  `-ffast-math` code generation of its sample directions, which flipped one grazing sample).

`provenance.txt` keeps the rule's ray count (the rays the one-at-a-time rule traces), so it matches
earlier builds. Diagnostics, each announced on stderr: `PTM_SOLAR_CHORD_PACKETS=0|4|8|16` (packet
width, 0 = one query per ray), `PTM_SOLAR_CHORD_PROFILE=1` (per-frame rays and thread time per stage;
the per-ray timers slow the run) or `=2` (per-frame wall time and class counts only), and
`PTM_SOLAR_CHORD_DUMP=<dir>` (the fitted edge of every partial face-frame).

## Results-changing speedups (off by default)

Each prints a NOTICE and is recorded in `provenance.txt`. Neither is a default: a default change needs
the restated G2 gate below and sign-off.

- `illumination.solar_disk_chord_check_adaptive_n: N` (default 0 = off). The screen probes are
  tested against the fitted edge with the check's own rules (margin, face plane). If every probe agrees,
  the fit is checked with N points instead of 128; otherwise with the full 128.
- `illumination.solar_disk_chord_fallback_n: N` (default 1024, range 8-4096): fallback samples.
- `illumination.solar_disk_chord_temporal: delta` (disk radii, default 0 = off). Each column's edge
  search starts from a bracket of +/-delta around the edge extrapolated linearly from the face's fits in
  the two previous calls, used only when the Sun moved by the same angle (within 5%, under 2 degrees) in
  both steps. If the bracket ends do not straddle an edge of the expected orientation, the full column
  runs. The edge check is unchanged. The result depends on the frame sequence: a chunked or restarted
  run differs from a single-shot run in the first frames after each start.

All three are fatal with a scheme other than `chord`, which would ignore them.

Measured (2026-11-07; AGT polar site, 24 hourly frames of its sunniest output day, and the restated G2
set below):

| option | AGT direct-solar time vs the identical default | restated G2 (Moon3D 140 frames) | verdict |
|---|---|---|---|
| `check_adaptive_n` 16 / 32 / 64 | not timed (fails G2) | FAIL: max dF 45 / 44 / 42 W/m2 (gate 20) | do not use |
| `fallback_n` 512 | fallback is about 1% of the time | PASS (max 19.5, dropout 4e-5, worst class 0.27%) | no useful gain |
| `fallback_n` 256 | as above | FAIL: max dF 26.8 | do not use |
| `temporal` 0.04 | -5% (33.0 vs 34.6-35.3 s) | not exercised: at the set's 5.9-h cadence the Sun steps more than 2 degrees, so no frame uses it (identical to the default by construction) | AGT vs the default: max dF 0.84 W/m2, none > 2 |

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

## Validation (restated G2, Moon3D R1 mesh, first 140 frames, against `reference`)

Gates (Fermi review 2026-11-05):
- energy within +/-0.1%;
- slivers with F_ref >= 2 W/m2 reported as 0 in at most 0.1% of cases;
- max |dF| <= 20 W/m2;
- at most 0.5% of the faces in each class (CT54, CT110, and +/-10 K margins around 54 and 110 K, from
  the R1 coupled initial-window T_max) see any |dF| > 2 W/m2.

The kelvin-equivalent p99 is reported but not gated.

| screen / check / fallback | energy | dropout | max dF (W/m2) | worst class fraction | keq p99 | verdict |
|---|---|---|---|---|---|---|
| 2/4 ladder / none / 64 (first version) | +7e-5 | 0.28% | 566 | 0.64% | 7.5 K | FAIL |
| 9 at 0.85 / none / 64 | +1.7e-3 | 1.57% | 228 | 1.15% | 21.5 K | FAIL |
| 9 at 0.999 / none / 64 | +4e-5 | 0.28% | 228 | 0.64% | 7.3 K | FAIL |
| 25 / 16 at 0.05 / 256 | +3e-5 | 0.025% | 89 | 0.39% | 3.0 K | FAIL |
| 25 / 64 at 0.01 / 1024 | +2e-5 | 0.001% | 42 | 0.18% | 1.3 K | FAIL |
| **25 / 128 at 0.005 / 1024 (default)** | **+1e-5** | **0.001%** | **19.5** | **0.14%** | **1.2 K** | **PASS** |

Energy is the sum over the face-frames, unweighted by face area.
