# Derivation of the primary-body reflected flux ("Jupiter-shine", "Earthshine")

This note derives the irradiance that a surface face receives from sunlight reflected by its
primary body, and states the assumptions and the limits of the result. The expression is
implemented in `src/cpp/illumination/illumination_model.cpp`.

## 1. Geometry and symbols

| symbol | meaning |
| --- | --- |
| `R` | primary's radius (`primary_radius_`, km) |
| `d` | face-to-primary-center distance (km) |
| `α` | primary's angular radius seen from the face, `α = asin(R/d)` |
| `θ` | angle between the face normal and the direction to the primary |
| `a` | phase angle = Sun–primary–face angle (0 = primary "full", π = "new") |
| `F` | solar irradiance at the primary, `S₀/(r_pri/AU)²` |
| `A_B` | primary's Bond albedo (`*_BOND_ALBEDO`) |
| `p` | primary's geometric albedo |
| `Φ(a)` | phase function, normalized `Φ(0) = 1` |
| `q` | phase integral, `q = 2∫₀^π Φ(a) sin a da`; `A_B = p·q` |

The far-field condition `d ≫ R` holds by 4 orders of magnitude for every primary
that PTM3D models. Jupiter seen from Europa is the most demanding case (`α = 5.98°`).

## 2. The result

For a sphere whose surface is Lambertian with hemispherical reflectance `A_B`
(a *Lambert sphere*), the irradiance on a face at distance `d` is

```
E = p · F · sin²α · Φ_L(a) · cos θ                                        (1)

    p     = (2/3) A_B                                                     (2)
    Φ_L(a) = [ sin a + (π − a) cos a ] / π                                (3)
```

`sin²α = (R/d)²`, so (1) is the familiar `E = p F (R/d)² Φ(a)` photometric law
projected onto the face by `cos θ`.

### 2.1 Geometric albedo `p = (2/3) A_B`

`A_B = p·q`, and for a Lambert sphere `q = 3/2` (evaluated numerically in §5), so
`p = A_B/(3/2) = (2/3) A_B`. For Jupiter, with `A_B = 0.343`, this gives `p = 0.2287`.

### 2.2 View factor `sin²α` rather than `Ω/4π`

`sin²α·cos θ` is the view factor from a differential area to a sphere, the
same factor used for the thermal emission of the primary (`Ω/π = 2(1−cos α)`, which
equals `sin²α` to 0.3 % at `α = 6°`). The quantity `Ω/(4π)` is the fraction of the full
sphere of directions that the primary occupies, and it is not an irradiance view
factor.

### 2.3 Derivation of (1)

Consider a surface element of the primary with incidence angle `i` and emission angle
`e`. For a Lambertian surface the emergent radiance is `L = A_B F cos i / π` for
`cos i > 0`. The flux reaching a distant observer is the integral of `L cos e` over
the part of the sphere that is both lit and visible (a lune), divided by `d²`:

```
E_⊥(a) = (1/d²) ∫_{cos i>0, cos e>0} (A_B F cos i / π) cos e dA            (4)
```

With `dA = R² dω` on the unit sphere this gives `E_⊥ = A_B F (R/d)² · (2/3)Φ_L(a)`,
which is (1) with (2) and (3). Section 5 evaluates (4) numerically and agrees with (1)
to 1 part in 10⁷.

The normalization can be checked in two independent ways.

* Zero phase. A Lambert disk of radius `R` facing the beam has radiance
  `F/π` (for `A=1`) and produces the irradiance `(F/π)(πR²)/d² = F(R/d)²`. This is
  the definition of `p`, so `E_⊥(0) = p F (R/d)²`.
* Energy. Integrating (1) over a sphere of radius `d` gives
  `∮E_⊥ d²dΩ = F p R² π q = F πR² (p q) = A_B F πR²`, which is exactly the power that
  the primary intercepts and reflects.

### 2.4 Limits

| phase | `Φ_L` | meaning |
| --- | --- | --- |
| `a = 0` (full) | 1 | maximum; face sees the fully lit hemisphere |
| `a = π/2` (quadrature) | `1/π = 0.3183` | half the disk lit, limb-darkened |
| `a = π` (new) | 0 | face sees only the primary's night side |

`Φ_L` decreases monotonically on `[0, π]`, with `Φ_L(π−ε) → ε³/(3π)`, a cubic
rather than quadratic approach to zero.

For Europa this has a practical consequence. Jupiter's phase angle is tied to local
time on the sub-Jovian hemisphere. Full Jupiter (`a = 0`) occurs when Europa is
between the Sun and Jupiter, which is local midnight on the sub-Jovian hemisphere,
and `a = π` occurs at local noon, which is also the eclipse configuration. The
maximum of `Φ_L` therefore coincides with the time at which the surface is coldest and
`dT/dF = T/(4εσT⁴)` is largest.

## 3. Comparison with an isotropic normalization

A simpler expression for this term, used in earlier versions of PTM3D, scales the
reflected flux by the fraction of the full sphere of directions that the primary
occupies and by the phase function `(1+cos a)/2`:

```cpp
reflected = solar_flux_at_primary * primary_albedo * phase_function
          * (solid_angle / (4.0 * M_PI)) * cos_angle;      // Φ_code=(1+cos a)/2
```

Since `Ω/(4π) = (1−cos α)/2 ≈ sin²α/4`, this is

```
E_code = F · (A_B/4) · sin²α · Φ_code(a) · cos θ,      Φ_code = (1+cos a)/2   (5)
```

Comparing (5), which carries the subscript "code", with (1) and (2):

| factor | Lambert sphere (true) | equation (5) (code) | ratio true/code |
| --- | --- | --- | --- |
| albedo | `p = (2/3)A_B` | `A_B` | 2/3 |
| view factor | `sin²α` | `sin²α/4` | 4 |
| phase | `Φ_L(a)` | `(1+cos a)/2` | `Φ_L/Φ_code` |

so

```
E_true / E_code = (8/3) · Φ_L(a) / Φ_code(a)                               (6)
```

The ratio is 8/3 = 2.667 at full phase, but it is not constant. It falls to 1.698 at
quadrature, equals 1 at `a = 128.7°`, and tends to 0 as `a → π`. Near new phase,
equation (5) overestimates the flux.

### 3.1 Energy conservation

The factor `1/4` in `Ω/(4π)` resembles the `1/4` in the radiosity `J = A_B F/4` of a
sphere that re-radiates isotropically, which might suggest that (5) is correctly
normalized. It is not.

A factor `1/4` does belong to an isotropic re-radiator: an isotropic sphere has
`Φ ≡ 1` and phase integral `q_iso = 4`, hence `p_iso = A_B/4`. Equation (5), however,
multiplies that isotropic normalization by a non-isotropic phase function
`Φ_code = (1+cos a)/2`, whose phase integral is

```
q_code = 2∫₀^π Φ_code sin a da = 2                                          (7)
```

The implied total reflected power is therefore `F πR² · p_iso · q_code
= F πR² · A_B · (2/4) = ½ A_B F πR²`. Equation (5) radiates exactly half of the
sunlight that the primary reflects, so it violates energy conservation by a factor of
2 in addition to having the wrong dependence on phase.

Multiplying (5) by 4 while retaining `Φ_code` does not correct it. That would give
`4 × ½ = 2×` the correct total power and 1.50× the correct flux at full phase. Without
the factor of 4, (5) is 2.667× low at full phase and 2× low in total energy. The
normalization and the phase function must be corrected together. Equation (1) has
`q = 3/2` and closes the energy budget exactly (§2.3).

## 4. Magnitude of Jupiter-shine at Europa

With the constants of the code (`SOLAR_CONSTANT 1361`, `JUPITER_RADIUS 69911 km`,
`JUPITER_BOND_ALBEDO 0.343`, `JUPITER_IR_FLUX 13.8`), Europa at
`d = 671,100 km` and Jupiter at 5.2038 AU, the values at full phase are as follows.
The line marked "fixed code" is equation (1) and the line marked "old code" is
equation (5).

```
F              = 1361/5.2038²      = 50.26  W/m²
α              = asin(69911/671100)= 5.978°     sin²α = 1.0853e-2
p (Lambert)    = (2/3)(0.343)      = 0.2287
E_refl(full)   = 50.26·0.2287·sin²α = 0.1247 W/m²      <-- fixed code
E_refl(pre-fix)= 50.26·0.0858·sin²α = 0.0468 W/m²      <-- old code
E_IR           = 13.8·sin²α         = 0.1498 W/m²      (primary IR, unchanged)
```

The thermal term can be checked through the total power, independently of any view
factor. Jupiter emits `4πR²·13.8 = 8.48e17 W`, which spread over `4πd²` at Europa
gives `0.1498 W/m²`, in agreement.

A value of order 1 W/m² is sometimes quoted for Jupiter-shine at Europa. The values
above are smaller, and the geometry sets a firm upper limit. Setting `p = 1` (a
perfectly white Jupiter) gives

```
E_refl^max = F · 1 · sin²α = 0.545 W/m²
```

and adding the thermal term gives a maximum total of `0.70 W/m² < 1 W/m²` for any
physically possible Jupiter. The correct magnitude is a few tenths of a W/m²:
`0.15 (IR) + 0.12 (reflected) ≈ 0.28 W/m²` at the sub-Jovian point at full phase, or
0.55 % of the 50 W/m² solar flux at Europa. We could not obtain 1 W/m² from any
combination of Jupiter's published radius, albedo and effective temperature, and the
model is not adjusted to match that value.

Using Jupiter's visual geometric albedo `p_V ≈ 0.52` in place of the Lambert value
0.229 would give 0.297 W/m². That combination is not internally consistent. It implies
`q = A_B/p = 0.343/0.52 = 0.66`, which is not Jupiter's phase integral (the
Cassini-based `A_B = 0.503 ± 0.012` of Li et al. 2018 implies `q ≈ 0.97`). In
addition, `p_V` is a V-band quantity and the flux here is bolometric. We therefore use
the Lambert sphere, which is consistent in energy with the `*_BOND_ALBEDO` constants
in `constants.cpp`. A body-specific pair `(p, Φ)` based on observations is not
implemented.

## 5. Numerical verification

The lune integral (4) was evaluated directly on a 4000×8000 surface grid and compared
with the closed form (1)–(3):

## 5. Numerical verification (`/tmp/jup_derive.py`, reproduced in the check script)

Direct evaluation of the lune integral (4) on a 4000×8000 surface grid vs the
closed form (1)–(3):

```
phase[deg]  Phi_L(closed)  numeric/(2/3)   rel.err
   0.0        1.000000       1.000000     +7.7e-08
  30.0        0.880843       0.880843     +9.0e-08
  45.0        0.755409       0.755409     +9.5e-08
  90.0        0.318310       0.318310     +1.0e-07
 135.0        0.048302       0.048302     +3.6e-07
 170.0        0.000562       0.000562     +7.8e-06
 180.0        0.000000       0.000000        exact
phase integral q (Lambert)          = 1.5000   (expect 3/2)
phase integral q ((1+cos a)/2)      = 2.0000   (the pre-fix code, eq. 7)
phase integral q (isotropic)        = 4.0000
```

## 6. Assumptions and limitations

1. Lambert sphere. Real primaries are not Lambertian, and Jupiter is
   backscattering. From §4, a Jupiter with its observed geometric albedo would be
   up to ~2.3× brighter at full phase than the Lambert value used here. This is the
   largest remaining physical uncertainty in the term.
2. Uniform primary. There is no zonal albedo structure and no shadows of Io or
   Ganymede.
3. Far field. `sin²α` is the exact `dA→sphere` view factor. The small-`α`
   step from `Ω/π` to `sin²α` is a 0.3 % effect at Europa and 0.008 % at the Moon.
4. Gray. Reflected light is absorbed by the target with `(1 − A_vis)`
   (`ptm_thermal_sim.cpp`), which is the appropriate spectral band. There is no
   color term.
5. Occlusion weighting. When the disk of the primary is partly hidden (§7),
   the visible fraction is computed as if the disk were uniformly bright.
   At intermediate phase the lit part of the disk is off-center, so a partial
   occultation is weighted slightly incorrectly. The error is bounded by the
   brightness contrast across the disk and applies only to the narrow band of faces
   that are partly occulted.

## 7. Fractional visibility of a finite disk

The primary is not a point source. A single ray to the center of the disk would make
its visibility all-or-nothing, and a face whose normal lies 0.1° below the direction
of the primary would receive exactly zero. The angular radius of the disk is 5.98°
for Jupiter seen from Europa and 0.95° for Earth seen from the Moon.

The generalization of `E = J·(Ω/π)·cos θ` to a partly hidden disk is

```
E = J · (Ω/π) · ⟨cos θ⟩,     ⟨cos θ⟩ ≡ (1/N) Σ_k  vis_k · max(cos θ_k, 0)   (8)
```

where the `N` samples cover the disk and `vis_k` is the terrain visibility of sample
`k`. `⟨cos θ⟩` reduces to `cos θ₀` when the disk is small and fully visible, so (8)
contains the point-source expression as a limit.

### 7.1 Analytic horizon term

On a smooth global mesh (an icosphere, as for Europa and for global models of the
Moon) there is no terrain, and the only effect of the finite disk is at the horizon,
where `cos θ ≤ 0`. This case has a closed form. Writing the altitude of the disk
center above the horizon plane of the face as `h = asin(cos θ₀)` and linearizing
`cos θ ≈ h + x` across the disk (`|x| ≤ α`):

```
⟨cos θ⟩ = α · g(h/α),        g(t) = [ tπ/2 + t²√(1−t²) + t·asin(t)
                                      + (2/3)(1−t²)^{3/2} ] / π           (9)
```

with `g(t≥1) = t`, `g(t≤−1) = 0`, and `g(0) = 2/(3π) = 0.2122`. A face whose normal
points exactly at the horizon therefore still receives 21 % of `α·J·Ω/π`, where a
point source would give 0. `g` is C¹ and monotonic on `[−1,1]`, so the step at the
horizon is replaced by a smooth ramp. The truncation error is `O(α²)` (≲1 % of an
already small partial-band flux at Jupiter/Europa, ≲0.03 % at Earth/Moon).

On an icosphere the band `|h| < α` covers a fraction `≈ sin α` of the faces, which is
10.4 % of a Europa mesh and 1.7 % of a global Moon mesh.

### 7.2 Terrain sampling

The visibility of the disk is first tested with 5 rays, toward the center of the disk
and 4 points on its rim, tracing only directions above the horizon of the face. If all
five are clear, `vis ≡ 1` and (9) is used. If all are blocked, `E = 0`. The rim points
are oriented relative to the horizon, `u = normalize(n̂ − (n̂·d̂)d̂)`, so that the
point `(+1,0)` is the highest point of the disk. With an arbitrary orientation all
five rays can miss a thin crescent above the horizon.

A smooth convex body would always let a ray above the horizon escape. A faceted
icosphere does not. On a Europa ico5 mesh, 1178 of 5.9 M face-steps give a mixed
result, and 697 have all five rays blocked by grazing intersections with neighboring
facets. Evaluating every case with 512 rays finds more such intersections (3857), so
they are real and not an artifact of the five-ray test, but they carry
`⟨cos θ⟩ ≤ 2.3e-6`, or ≤ 3.5e-7 W/m². Global meshes therefore require ~5 rays per
face instead of 1.

When the five rays give a mixed result, (8) is evaluated directly with N = 512 rays.
The sample points are a Hammersley sequence mapped concentrically onto the disk, the
same generator used for the solar disk. The error does not fall as `1/N`. A straight
edge crosses low-discrepancy points in clumps, so the error falls more slowly. The
table gives the largest error as the edge of an occluder is moved across the disk in
steps of 1/500 of the disk radius (`./build/bin/test_primary_disk` §4).

| N | worst \|sampled − analytic\| | 1/N |
| --- | --- | --- |
| 16 | 0.0837 | 0.0625 |
| 32 | 0.0484 | 0.0312 |
| 64 | 0.0381 | 0.0156 |
| 128 | 0.0175 | 0.0078 |
| 256 | 0.0159 | 0.0039 |
| **512** | **0.0082** | 0.0020 |

A power-law fit gives `worst ≈ 0.49·N^(−0.65)`, roughly the inverse square root of N
(an independent fit to the same data gives `N^(−0.68)`). Two other point sets at
N = 512 perform no better: a Fibonacci spiral gives 0.0092 and an equal-area polar
grid 0.0174. N = 512 is the smallest power of two that meets the `< 1 %` criterion.

The five-ray test is an approximation. An occluder that hides part of the interior of
the disk while missing the center and all four rim points is not detected. Setting
`PTM_PRIMARY_DISK_ALWAYS_SAMPLE=1` evaluates (8) with 512 rays in every case, for
validation.

### 7.3 Fully visible and fully hidden disks

The 512-ray evaluation is omitted when both

1. `h ≥ α` (the whole disk is above the horizon plane of the face), and
2. all 5 test rays are clear.

Then `⟨cos θ⟩ = cos θ₀` exactly, which is the point-source expression. Faces with
`h ≤ −α` receive 0, and no rays are traced. The finite-disk treatment therefore
affects only the faces between these two cases.

## 8. Effect on surface temperature

At Europa's sub-Jovian point at full phase (local midnight), the total Jupiter-shine
is `0.1247 + 0.1498 = 0.275 W/m²`, compared with `0.0468 + 0.1498 = 0.197` for
equation (5), a difference of `ΔF_inc = +0.078 W/m²` in incident flux. The surface
gains only the absorbed part. Reflected light is in the solar band and is absorbed
with `(1 − A_vis)` (`ptm_thermal_sim.cpp`), so with Europa's `A_vis = 0.55`,

```
ΔF_abs = (1 − 0.55) × 0.078 = 0.0351 W/m²
ΔT     = T · ΔF_abs / (4 ε σ T⁴),   ε = 0.95
```

| background T | ΔT from the difference in reflected flux (absorbed) |
| --- | --- |
| 130 K (Europa dusk) | +0.07 K |
| 90 K (Europa deep night) | +0.22 K |
| 70 K | +0.47 K |
| 50 K (cold trap) | +1.30 K |

The thermal emission of the primary is in the thermal band, is absorbed with `ε`, and
is not affected.

The horizon ramp of §7 changes the flux in both directions. Faces immediately below the
horizon of a point source gain up to `α·g·J·Ω/π`, and faces immediately above it lose the
difference between `cos θ₀` and (9).
