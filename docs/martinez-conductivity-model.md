# Regolith conductivity models: Hayne et al. (2017) and Martinez & Siegler (2021)

PTM3D can integrate either of two thermal-conductivity laws for planetary
regolith, selected for each run in the YAML configuration:

```yaml
thermophysical:
  conductivity_model: hayne2017     # DEFAULT — unchanged historical behavior
  # conductivity_model: martinez2021
```

An unrecognized value stops the run with an error. The default is not substituted,
because the choice of law affects every result of the run.

## 1. The two laws

### Hayne et al. (2017) (default)

```
k(T, z) = kc(z) * [ 1 + chi * (T/350)^3 ]
kc(z)   = kd - (kd - ks) * exp(-z/H)
rho(z)  = rhod - (rhod - rhos) * exp(-z/H)
```

This is Watson's equation, `k = A + B*T^3`, with a temperature-independent
solid-conduction term `A = kc(z)` and a radiative term whose coefficient is
tied to the solid term through `chi`. Its laboratory basis (Apollo fines,
Cremers & Birkebak) extends down to roughly 150 K. Below that temperature it is
an extrapolation, and the floors of permanently shadowed regions are in that range.

Reference: Hayne, P. O., et al. (2017), *JGR Planets* **122**, 2371,
doi:10.1002/2017JE005387, Appendix A.

### Martinez & Siegler (2021) (optional)

```
Eq. 10:  k(T, rho) = (A1*rho - A2) * k_am(T) + (B1*rho - B2) * T^3
Eq.  9:  k_am(T)   = Ac + Bc*T^-4 + Cc*T^-3 + Dc*T^-2 + Ec*T^-1
                        + Fc*T + Gc*T^2 + Hc*T^3 + Ic*T^4
```

| coefficient | value | | coefficient | value |
|---|---|---|---|---|
| `Ac` | -0.203297      | | `Fc` |  0.01101     |
| `Bc` | -11.472        | | `Gc` | -2.8049e-5   |
| `Cc` |  22.5793       | | `Hc` |  3.35837e-8  |
| `Dc` | -14.3084       | | `Ic` | -1.40021e-11 |
| `Ec` |  3.41742       | | | |

| | value |
|---|---|
| `A1` | 5.0821e-6  |
| `A2` | 5.1e-3     |
| `B1` | 2.0022e-13 |
| `B2` | 1.953e-10  |

`k` in W m^-1 K^-1, `rho` in kg m^-3, `T` in K.

`k_am(T)` is the amorphous-endmember solid conductivity of Woods-Robinson et al.
(2019, their Eq. 30), a "zeroth-order" fit to amorphous SiO2 data from Touloukian
& Ho (1970). The principal physical difference from the Hayne law is that the solid
term falls steeply toward low temperature instead of remaining constant.

Reference: Martinez, A., & Siegler, M. A. (2021), "A Global Thermal Conductivity
Model for Lunar Regolith at Low Temperatures", *JGR Planets* **126**(10),
e2021JE006829, doi:[10.1029/2021JE006829](https://doi.org/10.1029/2021JE006829).

## 2. Scope of the implementation

The implementation uses Eqs. 9 and 10 of the paper as published, with the published
coefficients. The other elements of the paper are not included. Each corresponds to
a separate configuration setting, and including them would confound the comparison
of the two conductivity laws.

| Paper element | Status here |
|---|---|
| Retuned density scale height `H = 5.6 cm` | Not applied. PTM3D keeps the `H` of the body (lunar default 6.0 cm) and prints a note whenever `H != 0.056` with this model. `H` cannot be set for a named body in the configuration, because `thermophysical.H` is not read (`KNOWN_LIMITATIONS.md`). It can be set only through a per-face `thermal.h_parameter_map` or for a generic body. Figures 4-11 of the paper therefore cannot be reproduced exactly from a lunar configuration. |
| Feng et al. (2020) bolometric Bond albedo (their Eq. 4) | Not applied. PTM3D offers Hayne-angular or Hapke. |
| Feng's `ks = 8.0e-4`, `kd = 3.8e-3` | Not relevant: `ks`/`kd` do not enter Eq. 10. |
| `cp(T)` | Identical. The paper uses the `cp` of Hayne et al. (2017). |
| Grid `m = 20`, `n = 30`, 2 m depth, `Qb = 0.018 W/m^2` | Set in the configuration, and not changed by the choice of model. |

Two consequences of this model should be stated with any result that uses it.

1. `ks`, `kd` and `chi` have no effect on conduction. Eq. 10 depends only on the
   local density `rho(z)`. A message at the start of the run states this.
2. On PTM3D's Moon the change is not confined to low temperatures. At the surface
   node (`rho = 1100 kg m^-3`) the Martinez conductivity is 0.63-0.69x the Hayne
   conductivity over the whole range 150-400 K. The two laws cross only near
   `rho ~ 1300-1400 kg m^-3`. This is the behavior of the published model (in
   Figure 3 of the paper, "New 1100" lies below "Hayne 1100" over 0-400 K), and it
   means that selecting this model also changes the daytime thermal inertia.

The ratio `k_Martinez / k_Hayne` with PTM3D's lunar defaults (`ks` 7.4e-4,
`kd` 3.4e-3, `chi` 2.7) is as follows.

| T (K) | rho = 1100 | rho = 1300 | rho = 1800 |
|---|---|---|---|
|  30 | 0.134 | 0.203 | 0.241 |
|  60 | 0.279 | 0.422 | 0.500 |
| 100 | 0.457 | 0.685 | 0.810 |
| 150 | 0.611 | 0.868 | 1.006 |
| 250 | 0.693 | 0.989 | 1.150 |
| 400 | 0.628 | 0.853 | 0.976 |

The surface thermal inertia `Gamma = sqrt(k*rho*cp)` (cp = 555 J kg^-1 K^-1,
rho = 1100) changes from 21.3 to 7.8 at 30 K, from 21.9 to 14.8 at 100 K, and from 29.9 to 24.9 at 250 K.

The figures `docs/figures/conductivity/fig_conductivity_models.{png,pdf}` and
`fig_conductivity_ratio.{png,pdf}` are generated by
`scripts/plot_conductivity_models.py`, which computes the curves from the published
coefficients independently of the C++ implementation.

## 3. Validity

The paper states the following.

* The fit of `k_am` to amorphous SiO2 is accurate to ±3% over 200-500 K, ±8% at
  50 K and 900 K, and ±15% below 10 K.
* The density dependence is based on the NU-LHT-2M simulant data of Zhong et al.
  (2016), 15-205 K, at a single bulk density of ~1700 kg m^-3, rescaled across
  1100-1800 kg m^-3.

The following properties are not stated in the paper. They follow from the
equations and are enforced in the code.

* `(A1*rho - A2)` changes sign at 1003.52 kg m^-3 and `(B1*rho - B2)` at
  975.43 kg m^-3, so Eq. 10 gives a negative conductivity below
  ~1004 kg m^-3. A run in which `min(rhos, rhod)` is below that value stops with
  an error. A run outside the fitted range of 1100-1800 kg m^-3 produces a warning
  and proceeds.
* `k_am(T) -> -inf` as `T -> 0` (the `Bc*T^-4` term). It is positive for
  T > 0.99 K, and `k_am(2.7 K) = 0.0605`. PTM3D limits every node to a minimum of
  `ABSOLUTE_ZERO_K = 2.7 K`, so the polynomial is not evaluated in its divergent
  region. The conductivity routine applies the same 2.7 K minimum itself, so that
  it is also safe for a lower input temperature.
* `k_am(T)` is non-monotonic between ~6.5 K (local max 0.1297) and ~11.9 K
  (local min 0.1240), a 4.4% dip that lies within the ±15% uncertainty quoted in
  the paper below 10 K. It is kept as published. `k` remains positive, so the
  diffusion operator remains well posed.

## 4. Code locations

| Piece | File |
|---|---|
| Law, coefficients, guards, parser | `src/cpp/thermal/conductivity_model.h` |
| Selector (carried to every solver) | `thermal::Config::conductivity_model`, `src/cpp/thermal/thermal_model.h` |
| CPU 1-D profile | `Profile::update_k()`, `thermal_model.cpp` |
| CPU mesh/grid fused refresh | `Profile::update_T()`, `thermal_model.cpp` |
| CPU SIMD batch | `updateMaterialPropertiesBatchSIMD()`, `thermal_model.cpp` |
| Initial geothermal seed | `Profile::init_T()`, `thermal_model.cpp` |
| GPU (live kernel) | `getEmbeddedKernelSource()` in `src/cpp/gpu/thermal_opencl.cpp` (`martinez_kam_at_or_above`, `martinez_k`, `compute_k_model`) |
| GPU host plumbing | `ThermalOpenCLSolver::setMaterialParams(..., int cond_model)` |
| YAML parsing | `src/cpp/tools/config_parser.{h,cpp}` |
| Banner + fatal guards | `src/cpp/main/ptm_thermal_sim.cpp` |
| Anchor tests | `src/cpp/tools/test_conductivity_model.cpp`, built as `build/bin/test_conductivity_model` |

Two features of the implementation should be preserved.

* The selector is carried by `thermal::Config` and not by `thermal::Planet`.
  `Config` is built once in `main` and passed to every solver entry point, whereas
  `Planet` is rebuilt from the body name at every entry point. Every solver entry
  point also reapplies the `ks_override`/`kd_override` settings
  (`applyConfiguredThermo()` in `ptm_thermal_sim.cpp`), and a run in which an
  override was not applied by any solver stops before writing output. Whether the
  selected model is in effect should still be confirmed by comparing the output of
  a run with that of the default model.
* `setMaterialParams` takes `cond_model` with no default argument, so that every
  GPU entry point must state the model explicitly. A default argument would allow
  a GPU entry point to use the Hayne law when the configuration requests Martinez.

`src/cpp/gpu/thermal_kernels.cl` is a reference copy that is not loaded by default
(the OpenCL kernels are built from source embedded in `thermal_opencl.cpp`), and
`src/cpp/gpu/thermal_kernels.cpp` is compiled only for the HIP and CUDA backends.
Neither contains the Martinez law.

## 5. Tests

`build/bin/test_conductivity_model` runs 47 checks, all of which must pass:

* A. `k_am(T)` and `k(T,rho)` agree with 40-digit decimal reference values
  computed from the published coefficients, at `rtol 1e-12`, over 10-400 K and
  rho = 1100 / 1300 / 1800 kg m^-3.
* A3. Domain limits: the density at which the sign changes, `k > 0` over T in
  [2.7, 500] K and rho in [1100, 1800], the behavior at the minimum temperature,
  and the density validity test.
* B. The float32 OpenCL implementation agrees with the double-precision host
  implementation (largest relative difference 1.03e-6). The test contains its own
  copy of the device expression, so a divergence between the two files causes this
  test to fail.
* C. The Hayne law gives results identical to the direct evaluation of its
  expression at 12 `(T, kc)` points.
* D. Config-string parsing accepts the valid spellings and rejects invalid ones.

In addition, any change to this code must keep `test_bowl_crater` passing
(4/4, `F_sky` error 0.492713%) and must leave the results of a run with the default
configuration unchanged.
