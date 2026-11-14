# Observing a simulated body

A telescope or spacecraft radiometer that does not resolve a body measures its thermal emission
integrated over the visible disk, and an imager that does resolve it measures the emission of
each part of the disk. PTM3D computes both from a finished simulation. Given an observer, a time
and a wavelength band, it returns the band radiance and band brightness temperature of every face
on the visible disk, the flux of the whole disk at the observer, and the disk brightness
temperature, for direct comparison with a measurement.

The calculation is available from Python as `ptm.observe` and from the command line as
`ptm observe`. [Tutorial 7](tutorials/07_observing_thermal_infrared.md) applies it to Europa
seen from Earth at 8 to 14 μm, including fits of thermal inertia and emissivity to synthetic
observations and the temperatures below the surface during a Jupiter eclipse.

## Thermal emission in the infrared

In the thermal infrared the surface is opaque. Radiation emitted below the uppermost tens of
micrometers is absorbed before it reaches the surface, and this layer is far thinner than the
depth reached by the diurnal temperature wave. Each face therefore emits at its surface
temperature *T*, and no integration over depth is needed.

An instrument records radiance weighted by its relative spectral response *R*(λ). The band
radiance of a face is

```math
L = \varepsilon \int R(\lambda)\, B_\lambda(T)\, d\lambda, \qquad
B_\lambda(T) = \frac{2 h c^2}{\lambda^5}\, \frac{1}{e^{h c/\lambda k T} - 1}
```

where *B*<sub>λ</sub> is the Planck function and ε the emissivity. We take the emissivity to be
gray and Lambertian. It is the same at every wavelength, and so the same as the emissivity in the
energy balance that produced the temperatures, and the emitted radiance does not depend on
direction. The default response is a top hat, equal to 1 from 8 to 14 μm and 0 elsewhere. Any
band limits can be given, or a tabulated response.

The integral is evaluated by Gauss-Legendre quadrature with 16 nodes on each of a set of panels
whose long and short wavelengths differ by at most 10 percent. A tabulated response is taken as
linear between its wavelengths, which are added to the panel edges. For a top-hat band the
quadrature agrees with the closed-form series for the integral of the Planck function to a few
parts in 10<sup>12</sup> from 20 K to the temperature of the Sun.

### The visible disk

The flux at a distant observer is a sum over the faces turned toward it,

```math
F = \sum_j L_j \, \frac{A_j \mu_j}{\Delta_j^2}
```

where *A*<sub>j</sub> is the area of face *j*, μ<sub>j</sub> the cosine of its emission angle and
Δ<sub>j</sub> its distance from the observer. The terms *A*<sub>j</sub>μ<sub>j</sub>/Δ<sub>j</sub><sup>2</sup>
are the solid angles of the faces, and their sum Ω is the solid angle of the disk. The mean band
radiance of the disk is *F*/Ω, and the band brightness temperature *T*<sub>b</sub> of the disk is
the temperature of a blackbody with that band radiance,

```math
\int R(\lambda)\, B_\lambda(T_\mathrm{b})\, d\lambda = \frac{F}{\Omega}.
```

The band radiance increases monotonically with temperature, and *T*<sub>b</sub> is found by
Newton's method on the logarithm of the band radiance. The same definition applied to each face
gives a map of band brightness temperature over the visible disk, which is what an imager
records. Because Ω is the solid angle of the mesh itself, an isothermal blackbody mesh of any
resolution has a disk brightness temperature equal to its temperature. The flux *F* converges to
that of the smooth body as the mesh is refined, with an error that falls by a factor of four
with each subdivision of an icosphere.

On a convex closed mesh every face turned toward the observer is fully visible, so the sum is
exact for the faceted surface. Meshes that are not convex raise an error (see below).

### Sunlight reflected in the band

A gray surface with emissivity ε reflects the fraction 1 − ε of the radiation it receives at
every wavelength (Kirchhoff's law), and a Lambertian surface reflects it equally in all
directions. The reflected band radiance of a face is

```math
L_\mathrm{r} = \frac{(1 - \varepsilon)\, f_\odot\, E}{\pi}, \qquad
f_\odot = \frac{\pi \int R(\lambda)\, B_\lambda(T_\odot)\, d\lambda}{\sigma T_\odot^4}
```

where *E* is the direct sunlight incident on the face, as the simulation recorded it (it
includes the incidence angle and any eclipse), and *f*<sub>⊙</sub> is the fraction of the Sun's
irradiance that falls in the band for a blackbody Sun at *T*<sub>⊙</sub> = 5772 K. That
blackbody, at the nominal solar radius, gives 1361 W m<sup>−2</sup> at 1 au, the measured total
solar irradiance, so the in-band solar irradiance scales with heliocentric distance exactly as
the simulation's sunlight does. From 8 to 14 μm, *f*<sub>⊙</sub> = 1.108 × 10<sup>−3</sup>.

It is tempting to neglect reflected sunlight in the thermal infrared. At 8 to 14 μm, however,
the thermal emission of a cold surface lies far down the Wien tail of the Planck function, while
the Sun's lies near its peak. For Europa with ε = 0.9 seen from Earth at 12:00 UTC on
2025-01-01, the reflected light adds 1.4 percent to the flux of the disk, and it exceeds 1 percent
of the thermal radiance on three quarters of the visible disk. Shortly after an eclipse by
Jupiter, when the surface is still cold but sunlit again, the reflected light adds 4 percent to
the flux. The term is therefore included by default (`--no-reflected` or `reflected=False`
leaves it out). It raises the disk brightness temperature by 0.2 to 0.4 K in these examples.

On the coldest sunlit faces, near the poles where the Sun is low, the reflected sunlight exceeds
the thermal emission by a factor that grows steeply as the surface cools. In the tutorial's
model these faces are at about 41 K and the factor reaches 2 × 10<sup>4</sup>, but they add
almost nothing to the flux of the disk.

### Dominance of the warmest faces

For surfaces at 90 to 130 K the Planck function peaks at 22 to 32 μm, and the band from 8 to
14 μm lies on its Wien side. There the band radiance rises steeply with temperature, roughly as
*T*<sup>*n*</sup> with *n* ≈ *hc*/λ*kT*. Evaluated for the 8 to 14 μm band, *n* is 12.7 at 90 K
and 9.2 at 130 K. A 1 percent change in temperature therefore changes the band radiance by about
10 percent, ten times the effect of a 1 percent change in emissivity at fixed temperature.

Three consequences follow, and the tutorial measures each of them. The disk brightness
temperature is weighted toward the warmest faces, and it exceeds the mean temperature of the
visible disk by several kelvin even with an emissivity below 1. A property that changes the
peak daytime temperature is amplified about tenfold in the band flux, so on Europa a 1 percent
increase in the Bond albedo lowers the 8 to 14 μm flux by about 3 percent. And the emissivity
acts mainly through the temperature rather than through the emitted radiance. A gray emissivity
also enters the energy balance, where a higher emissivity cools the surface, and on Europa this
effect on the band flux is larger than the direct one and of the opposite sign. Thermal inertia
changes the daytime temperature of Europa only weakly, and its clearest signature is the cooling
and recovery of the surface during an eclipse.

## The observer and the time

Name the observer by its SPICE name or NAIF id. Earth is the default. Any body or spacecraft in
the loaded kernels can be used, and further kernels such as a spacecraft trajectory are added
with `--kernel` (or `kernels=` in Python).

The time you give is the time the observer receives the light, in UTC. The light left the body
one light time earlier, and the body's orientation and temperatures are taken at that emission
epoch, with the light-time and stellar-aberration corrections of SPICE (`CN+S`). For Europa seen
from Earth the light time is about 35 minutes, during which Europa rotates about 2.5°. The
temperatures, and the direct sunlight used for the reflected light, are interpolated linearly
in time between the two output steps around the emission epoch. With `--nearest-step`
(`time_match="nearest"`) the observation is instead moved by up to half an output step, and the
shift is reported. A time whose emission epoch lies outside the output window is an error.

Linear interpolation in time is adequate for a slowly rotating body but not across an eclipse,
where the surface cools and warms within minutes of the shadow contacts. A model of eclipse
observations needs output steps of a minute or less near the contacts.

When the body is a satellite, the geometry also reports whether its planet hides it from the
observer. The flag is `visible`, `partially occulted`, `totally occulted` or `in transit`, the
last meaning in front of the planet. An occultation produces a warning, because the flux
returned is always that of the whole visible disk.

For bodies whose rotation is not a SPICE frame, such as the generic bodies of the asteroid
tutorials, give the geometry directly with `--sub-observer LAT LON DIST_KM` (or
`observe.explicit_geometry`). The latitude is planetocentric, the longitude is east-positive in
the frame of the mesh, and the distance is from the body center. Choose the output step
with `--step`. No SPICE kernels are needed.

## Results and units

For a band (`--band-um` or `--response`, `observe.band_emission`):

| Quantity | Units | Meaning |
|---|---|---|
| `flux_w_m2` | W m<sup>−2</sup> | In-band flux at the observer, thermal and reflected |
| `thermal_flux_w_m2`, `reflected_flux_w_m2` | W m<sup>−2</sup> | The two parts of the flux |
| `mean_radiance_w_m2_sr` | W m<sup>−2</sup> sr<sup>−1</sup> | Flux divided by the solid angle of the visible faces |
| `tb_k` | K | Band brightness temperature of the disk |
| `tb_thermal_k` | K | The same without the reflected sunlight |
| `mean_temperature_k` | K | Mean surface temperature of the visible disk, weighted by projected area |
| `reflected_fraction_disk` | none | Reflected over thermal flux of the disk |
| `reflected_fraction_max_face` | none | Largest ratio of reflected to thermal radiance on any visible face |
| `omega_projected_sr` | sr | Solid angle of the visible faces |
| `omega_nominal_sr` | sr | π*R*<sup>2</sup>/*D*<sup>2</sup> for the median vertex radius *R* and distance *D* |
| `sub_observer_lat_deg`, `sub_observer_lon_deg` | degrees | Sub-observer point at the emission epoch |
| `sub_solar_lat_deg`, `sub_solar_lon_deg`, `phase_angle_deg` | degrees | Illumination at the emission epoch |
| `light_time_s` | s | Light time from the body to the observer |

In Python the result also holds, for every face, the emission cosine, the solid angle, the
thermal and reflected band radiance, and the band brightness temperature (NaN for faces turned
away from the observer). From the command line, `--faces-csv FILE` writes these for every visible
face at every time. The band radiances are per unit of the response, so with a top-hat response
of 1 they are in W m<sup>−2</sup> sr<sup>−1</sup>. Brightness temperatures do not depend on the
scale of the response.

The ratio of reflected to thermal radiance is NaN when a visible face that receives sunlight has a
thermal band radiance of zero, which happens only below a few kelvin, and a warning is issued.

A response file has two columns, the wavelength in micrometers and the relative response,
separated by spaces, tabs or commas. Lines starting with `#` are comments.

## From Python

```python
from ptm import Results, observe

run = Results("output/quickstart_europa_diurnal")
geom = observe.observer_geometry(run, "2025-01-01T12:00:00", observer="EARTH")
band = observe.Band.top_hat(8.0, 14.0)            # or observe.Band.from_file("response.txt")
em = observe.band_emission(run, observe.temperature_at(run, geom), geom, band=band,
                           emissivity=0.9,
                           direct_solar=observe.temperature_at(run, geom, "direct_solar"))
print(f"{em.flux_w_m2:.4e} W m-2, T_b {em.tb_k:.2f} K, "
      f"reflected/thermal {em.reflected_fraction_disk:.4f}")
```

```
1.7935e-12 W m-2, T_b 119.21 K, reflected/thermal 0.0135
```

`band.radiance(T)` and `band.brightness_temperature(L)` convert between temperature and band
radiance for any array. `observe.explicit_geometry(lat, lon, distance_km, run=run, step=n)`
replaces `observer_geometry` when the geometry is given directly. Each result also has
`as_row()`, a dictionary of every quantity in the table above.

## From the command line

```bash
ptm observe output/quickstart_europa_diurnal \
    --time 2025-01-01T12:00:00 --band-um 8 14 --emissivity 0.9
ptm observe output/quickstart_europa_diurnal \
    --time 2025-01-01T12:00:00 --response my_filter.txt --emissivity 0.9 --faces-csv faces.csv
ptm observe output/quickstart_europa_diurnal --step 40 \
    --sub-observer 0 90 6.3e8 --band-um 8 14 --emissivity 0.9
```

Repeat `--time` for several observations, and add `--csv FILE` to write the table with every
column. `ptm observe --help` lists the options.

## Temperatures at depth

The thermal-infrared emission comes from the surface, but the temperature below the surface
shows how much heat the near surface stores and how far a change in sunlight, such as an
eclipse, reaches. `ptm observe` gives the temperature of every visible face at chosen depths, at
the same emission epoch, with the same light-time correction, interpolation in time and
occultation flags as the surface temperature. The run must store these temperatures at every
output step:

```yaml
output:
  map_depths_m: [0.005, 0.01, 0.02, 0.05]
  map_depth_snap_to_node: false   # store these depths, not the nearest grid nodes
  map_depths_per_step: true       # at every output step, not only statistics
```

Depths are in meters below the surface. The surface itself needs no entry, because the surface
temperature is always stored. Each depth adds one value per face per output step, 4 bytes, so
four depths add 28 MB to the 300 one-minute steps of `examples/quickstart_europa_eclipse.yaml`
on 5,120 faces. These frames are written by the CPU solver, without radiation exchange between
faces or `thermal.batch_checkpoint`.

```bash
ptm observe output/quickstart_europa_eclipse --list-depths
ptm observe output/quickstart_europa_eclipse --time 2025-01-01T21:00:00 --depth-m 0 0.0075 0.02
```

```
quickstart_europa_eclipse: temperatures stored at every output step at 5 depths
depth m   requested m   field
0         -             temperature  the surface
0.005     0.005         temperature_z0p005m  between grid nodes
0.01      0.01          temperature_z0p010m  between grid nodes
0.02      0.02          temperature_z0p020m  between grid nodes
0.05      0.05          temperature_z0p050m  between grid nodes
the solver's depth grid: 21 nodes from 0 to 0.5306 m, first spacing 0.00376 m
other depths from 0 to the deepest level are interpolated linearly in depth between these levels
quickstart_europa_eclipse: temperatures at depth, 3 depth(s)
observer EARTH, times corrected for light time (CN+S), temperatures interpolated linearly in time between the output steps shown
temperatures at depth: stored levels 0, 0.005, 0.01, 0.02, 0.05 m (0 is the surface)
  0.0075 m is not a stored level: interpolated linearly in depth between 0.005 and 0.01 m
mean temperature of the visible disk at each depth, weighted by projected area
observed UTC  emitted UTC (at body)  steps  obs lat  obs lon  view  T 0 m K  T 0.0075 m K  T 0.02 m K
2025-01-01T21:00:00.000  2025-01-01T20:25:03.150  160-161  2.751  -11.070  visible  97.839  100.689  102.731
```

`--depth-m` takes one or more depths in meters, with 0 for the surface. It can be given alone,
as here, or with a band or frequency, in which case the table of temperatures at depth follows
the emission table. For each observation it gives the mean temperature of the visible disk at
each depth, weighted by the solid angle of each face as `mean_temperature_k` is. The outputs
are:

| Option | Contents |
|---|---|
| `--csv FILE` | One row per observation, with the column `mean_temperature_k_z<depth>m` for each depth |
| `--faces-csv FILE` | One row per visible face and observation, with the column `temperature_k_z<depth>m` for each depth (after the band columns when a band is given) |
| `--depths-npz FILE` | A NumPy archive with `temperature_k` of shape (time, depth, face) for every face, visible or not. It also holds the emission cosine `mu` and solid angle `solid_angle_sr` of each face at each time (zero or negative for faces turned away), the disk means `mean_temperature_k` (time, depth), the depths with `interpolated`, `lower_m`, `upper_m` and `stored_depths_m`, the times and occultation flags, and the `lat_deg`, `lon_deg`, `area_m2`, `center_m` and `normal` of each face |

From Python, `observe.temperature_at_depths(run, geom, depths_m)` returns the temperatures as an
array of shape (depth, face), `observe.disk_mean_temperature(run, T, geom)` averages them over
the visible disk, and `observe.stored_depths(run)` lists the stored levels.

```python
from ptm import Results, observe

run = Results("output/quickstart_europa_eclipse")
geom = observe.observer_geometry(run, "2025-01-01T21:00:00", observer="EARTH")
dt = observe.temperature_at_depths(run, geom, [0.0, 0.0075, 0.02])
print(dt.temperature_k.shape, dt.interpolated)
print(observe.disk_mean_temperature(run, dt.temperature_k, geom).round(2))
```

```
(3, 5120) [False  True False]
[ 97.84 100.69 102.73]
```

### Interpolation in depth

The stored levels are the surface and the depths in `output.map_depths_m`. A requested depth
equal to a stored level returns that level's temperatures exactly. A depth between two levels
*z*<sub>k</sub> and *z*<sub>k+1</sub> is interpolated linearly in depth, and the output names the
two levels. For a temperature wave of amplitude *A* at the surface and thermal skin depth δ, the
error of linear interpolation between levels a distance *h* apart is at most

```math
\frac{A h^2}{4 \delta^2} \, e^{-z_k/\delta},
```

which is 1.6 percent of the amplitude for *h* = δ/4 and falls as *h*<sup>2</sup>. The error is
largest where the profile is most curved, near the surface and shortly after a sudden change in
sunlight. In the eclipse simulation, the temperature at 1 cm interpolated between the levels at
5 mm and 2 cm differs from the stored level at 1 cm by up to 1.5 K, about 25 minutes after Europa
leaves Jupiter's shadow. Store the depths you need rather than relying on interpolation across
wide gaps. With `map_depth_snap_to_node: false` the engine itself interpolates each stored level
linearly between the nodes of its depth grid, and `statistics/depth_map_z.json` records the
error estimate of each level. With the default, `true`, each requested depth is moved to the
nearest node and the stored depth is the node's.

A depth below the deepest stored level is an error that names the stored range, because the
temperatures are never extrapolated. A depth below the surface in a run that stored no depth
frames is an error that names the keys to add. Depth 0 needs no depth frames.

## Subsurface emission at millimeter wavelengths (experimental)

This mode is under development and has not yet been validated against observations.

At millimeter and centimeter wavelengths ice and regolith are partly transparent, and the
emission comes from a layer about one electrical skin depth *L*<sub>e</sub> thick, where
*L*<sub>e</sub> is the depth over which the emitted intensity falls by a factor of e. When the
ratio β = *L*<sub>e</sub>/δ of the electrical skin depth to the diurnal thermal skin depth δ is
much less than 1 the emission comes from the surface, and when β is near 1 or larger it samples
the temperature wave below the surface. For porous water ice on Europa at 233 GHz, β is
about 2.

`ptm observe` works at a single frequency with `--freq-ghz` or `--wavelength-um` and returns the
flux density in janskys and the disk brightness temperature, in two forms. In the surface form
(`--emissivity`), each face emits *e* B<sub>ν</sub>(*T*<sub>s</sub>) from its surface
temperature, with *e* an effective emissivity for the observing frequency. In the subsurface
form, each face emits

I<sub>ν</sub> = t(θ) ∫ B<sub>ν</sub>(T(z)) (α/μ<sub>t</sub>) exp(−α z/μ<sub>t</sub>) dz,

integrated along the ray refracted into the surface at that face's emission angle θ, with
α = 1/*L*<sub>e</sub>, μ<sub>t</sub> the cosine of the refracted angle and t(θ) the Fresnel
transmissivity averaged over polarization. The electrical properties are given as a porosity of
water ice (`--porosity`, from the mixing model in `scripts/emission/dielectric_mixing.py`) or
directly (`--skin-depth-m` with `--refractive-index`). An emissivity fitted with the surface form
already contains the subsurface effect and must not be used in the subsurface form.

The subsurface form needs the temperature of every face at a set of depths at every output step,
which these keys write.

```yaml
output:
  map_depths_m: [0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.0479, 0.0766,
                 0.1225, 0.1957, 0.3128, 0.5]
  map_depth_snap_to_node: false   # interpolate to the listed depths
  map_depths_per_step: true       # every output step, not only statistics
```

`examples/quickstart_europa_observing.yaml` is the Europa quickstart with these keys. The first
depth spacing should be at most 0.16 times the shorter of the thermal skin depth and
*L*<sub>e</sub>μ<sub>t</sub> for the most oblique ray, and the deepest depth at least 7
*L*<sub>e</sub>. A warning is issued when either is not met. These frames are written by the
CPU solver only, without radiation exchange between faces or `thermal.batch_checkpoint`.

```bash
ptm observe output/quickstart_europa_observing \
    --time 2025-01-01T12:00:00 --freq-ghz 233 --porosity 0.5
```

At a single frequency the brightness temperatures use the nominal solid angle
π*R*<sup>2</sup>/*D*<sup>2</sup> of a sphere of the mesh's median vertex radius, and both the
Planck and the Rayleigh-Jeans brightness temperatures are reported. They differ by about
*h*ν/2*k*, 5.6 K at 233 GHz. On a sphere of 5,120 faces this nominal solid angle differs from the
solid angle of the visible faces, used in the infrared, by about 0.1 percent, and on a body that is
not a sphere it can differ by much more. To compare with a measurement, compare fluxes rather
than brightness temperatures. `examples/python/fit_europa_disk_tb.py` and
`examples/python/fit_europa_disk_tb_subsurface.py` fit thermal inertia with the two forms to
synthetic observations.

## Limitations and caveats

Each face is treated as a smooth Lambertian surface. Roughness at scales smaller than a face
makes a surface appear warmer than a smooth one when it is seen from close to the direction of
the Sun. This effect is called thermal beaming. Earth-based observations of Europa are all made
at phase angles below 12°, where beaming can change the 8 to 14 μm flux significantly. Thermal
properties fitted with a smooth model absorb this effect. The emissivity is taken as gray and
independent of direction. Water ice has a higher emissivity between 8 and 14 μm than the value of
0.9 used for Europa in the tutorial, so the reflected sunlight computed with that value is an
upper estimate. Jupiter's thermal emission reflected by Europa is left out. It is small compared
with the thermal emission of the sunlit disk. The instrument's beam and point-spread function are
not modeled. At millimeter wavelengths the calculation also neglects scattering by grains and
voids, layering of the electrical properties, polarization and reflected sunlight.

The disk sum assumes a convex body, on which every face turned toward the observer is fully
visible. Shape models with concavities and terrain patches stop with an error. Maps are given
face by face, and images resampled onto an instrument's pixels are not yet available.
