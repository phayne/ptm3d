# Observing a simulated body

A telescope or spacecraft radiometer that does not resolve a body measures its thermal emission
integrated over the visible disk. PTM3D computes that quantity from a finished simulation. Given
an observer, a time and a frequency, it returns the flux density of the whole visible disk and
the disk brightness temperature, for direct comparison with a measurement.

The calculation is available from Python as `ptm.observe` and from the command line as
`ptm observe`. [Tutorial 7](tutorials/07_observing_europa.md) applies it to Europa from start to
finish, including fits of thermal inertia to synthetic observations.

## Surface and subsurface emission

Thermal emission leaves a surface from a layer whose thickness depends on the wavelength. In the
thermal infrared the layer is a few micrometers to millimeters thick, far thinner than the
depth reached by the diurnal temperature wave. At millimeter and centimeter wavelengths ice and
regolith are partly transparent, and the emission comes from a layer about one electrical skin
depth *L*<sub>e</sub> thick, where *L*<sub>e</sub> is the depth over which the emitted
intensity falls by a factor of e. The ratio of the electrical skin depth to the diurnal thermal
skin depth δ,

β = *L*<sub>e</sub> / δ,

determines which description applies. When β is much less than 1 the emission comes from the
surface. When β is near 1 or larger, the emission samples the temperature wave below the surface
and its variations are smaller and later than those of the surface. For porous water ice on
Europa at 233 GHz, β is about 2.

PTM3D offers both descriptions.

**Surface mode.** Each face emits *e* B<sub>ν</sub>(*T*<sub>s</sub>), with *T*<sub>s</sub> its
surface temperature and *e* an emissivity for the observing band. Use it in the thermal infrared.
At millimeter wavelengths a fitted band emissivity absorbs the difference between the surface
and the emitting layer. It then trades off against thermal inertia and has no simple physical
meaning.

**Subsurface mode.** Each face emits the intensity

I<sub>ν</sub> = t(θ) ∫ B<sub>ν</sub>(T(z)) (α/μ<sub>t</sub>) exp(−α z/μ<sub>t</sub>) dz,

integrated along the ray refracted into the surface at that face's emission angle θ. Here α =
1/*L*<sub>e</sub> is the absorption coefficient, μ<sub>t</sub> the cosine of the refracted
angle, and t(θ) the Fresnel transmissivity of the surface, averaged over polarization. The
transmissivity follows from the refractive index and is not a free parameter. This mode needs
the temperature at depth at the time of the observation (see below). Do not use an emissivity
fitted with the surface mode here, since that value already contains the subsurface effect.

The emission angle, the refraction, the Fresnel transmissivity and the Planck function are
evaluated for every face separately. Limb darkening and the shape of the emitting layer on the
disk are therefore part of the result.

## The observer and the time

Name the observer by its SPICE name or NAIF id. Earth is the default. Any body or spacecraft in
the loaded kernels can be used, and further kernels such as a spacecraft trajectory are added
with `--kernel` (or `kernels=` in Python).

The time you give is the time the observer receives the light, in UTC. The light left the body
one light time earlier, and the body's orientation and temperatures are taken at that emission
epoch, with the light-time and stellar-aberration corrections of SPICE (`CN+S`). For Europa seen
from Earth the light time is about 35 minutes, during which Europa rotates about 2.5°. The
temperatures are interpolated linearly in time between the two output steps around the
emission epoch. With `--nearest-step` (`time_match="nearest"`) the observation is instead moved
by up to half an output step, and the shift is reported. A time whose emission epoch lies
outside the output window is an error.

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

| Quantity | Units | Meaning |
|---|---|---|
| `flux_density_jy` | Jy | Flux density at the observer, the irradiance on a plane facing the body |
| `tb_planck_k` | K | Disk brightness temperature, the inverse Planck function of the mean intensity |
| `tb_rayleigh_jeans_k` | K | Disk brightness temperature in the Rayleigh-Jeans convention |
| `omega_nominal_sr` | sr | Solid angle π*R*<sup>2</sup>/*D*<sup>2</sup> used for the brightness temperatures |
| `omega_projected_sr` | sr | Solid angle of the visible faces of the mesh |
| `sub_observer_lat_deg`, `sub_observer_lon_deg` | degrees | Sub-observer point at the emission epoch |
| `sub_solar_lat_deg`, `sub_solar_lon_deg`, `phase_angle_deg` | degrees | Illumination at the emission epoch |
| `light_time_s` | s | Light time from the body to the observer |
| `skin_depth_m`, `refractive_index`, `beta_nadir` | m, none, none | Subsurface mode only |

The brightness temperatures use the solid angle of a sphere of the mesh's median vertex radius,
so a uniform sphere at temperature *T* with emissivity 1 has a brightness temperature of *T*. A
faceted sphere is slightly smaller than that sphere. For an isothermal 100 K icosphere of 5,120
faces the brightness temperature is 0.12 K low, and the deficit falls by a factor of four with
each subdivision.

The Planck and Rayleigh-Jeans brightness temperatures differ by about *h*ν/2*k*, 5.6 K at
233 GHz. Radio astronomers usually quote the Rayleigh-Jeans value. State which one a comparison
uses.

## Temperatures at depth for the subsurface mode

The subsurface mode needs the temperature of every face at a set of depths at the time of the
observation. Three keys in `output` write these.

```yaml
output:
  map_depths_m: [0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.0479, 0.0766,
                 0.1225, 0.1957, 0.3128, 0.5]
  map_depth_snap_to_node: false   # interpolate to the listed depths
  map_depths_per_step: true       # every output step, not only statistics
```

Each depth gets a directory of per-step files beside `temperature/`, named for the depth (for
example `temperature_z0p010m/`). They have the same format and numbering as the surface
temperature frames, and `Results.field("temperature_z0p010m", step)` reads them. The list must
start at 0. The files take 4 bytes per face, step and depth, which is 49 MiB for the list above
on the Europa quickstart.

The accuracy of the emission integral is set by this list. Two criteria apply, and a warning
is issued when either is not met.

- The first spacing should be at most 0.16 times the shorter of the thermal skin depth and
  *L*<sub>e</sub>μ<sub>t</sub> for the most oblique visible ray. For a diurnal wave of 30 K this
  keeps the error of the integral below 0.1 K.
- The deepest depth should be at least 7 *L*<sub>e</sub>, so that less than 0.1 percent of the
  emission comes from below it. Below the deepest depth its temperature is assumed to continue.

The deepest depth cannot exceed the bottom of the model column, which is `thermal.thermal_b`
diurnal skin depths. On Europa the default column of 0.53 m is deep enough at 233 GHz. At
97.5 GHz the list must reach about 2.7 m, which needs a larger `thermal_b`.

These frames are written by the CPU solver. They are not available with the GPU solver, with
radiation exchange between faces, or with `thermal.batch_checkpoint`. With any of these the run
stops with an error before the simulation begins. A run split into chunks writes frames in each
chunk directory. Join them with `scripts/stitch_frame_dirs.py --field <name>`, once for each
depth.

## Electrical properties

The subsurface mode takes the electrical properties in one of two forms.

```python
observe.Subsurface(porosity=0.5)                                  # porous ice at 100 K
observe.Subsurface(porosity=0.5, reference_temperature_k=None)    # at each node's temperature
observe.Subsurface(skin_depth_m=0.07, refractive_index=1.37)      # given directly
```

With a porosity, the permittivity of porous water ice comes from the mixing model in
`scripts/emission/dielectric_mixing.py`, at a reference temperature or at each node's own
temperature. A `skin_depth_m` given with a porosity replaces the model's absorption and keeps
its refractive index, which is how a fit varies *L*<sub>e</sub>. The loss tangent of ice rises
with frequency, so *L*<sub>e</sub> falls faster than the inverse of the frequency. Take each
frequency's skin depth from the model rather than scaling one value. The table gives the model's
values for porosity 0.5 at 100 K.

| Frequency (GHz) | *L*<sub>e</sub> (cm) | Refractive index |
|---|---|---|
| 97.5 | 39.05 | 1.374 |
| 233 | 6.84 | 1.374 |

For a single column, such as the depth profile of one site, `brightness_temperature` in
`scripts/emission/depth_weighted_emission.py` evaluates the same integral and returns the Planck
and Rayleigh-Jeans brightness temperatures.

## What is not included

The calculation neglects scattering by grains and voids, changes of the electrical properties
with depth other than through temperature, surface roughness and thermal beaming,
polarization, and the instrument beam. It also leaves out radiation reflected by the surface at
the observing frequency, such as sunlight or the thermal emission of a nearby planet. The cosmic
microwave background has a Rayleigh-Jeans temperature of 0.19 K at 233 GHz, and interferometers
do not detect it.

The disk sum assumes a convex body, on which every face turned toward the observer is fully
visible. Shape models with concavities and terrain patches stop with an error, and resolved
images of the disk are not yet available.

## From Python

```python
from ptm import Results, observe

run = Results("output/quickstart_europa_observing")
geom = observe.observer_geometry(run, "2025-01-01T12:00:00", observer="EARTH")
surface = observe.disk_emission(run, observe.temperature_at(run, geom), geom,
                                frequency_ghz=233.0, emissivity=0.9)
depths, T = observe.depth_temperatures_at(run, geom)
subsurface = observe.disk_emission(run, T, geom, frequency_ghz=233.0,
                                   subsurface=observe.Subsurface(porosity=0.5),
                                   depths_m=depths)
print(f"surface    {surface.flux_density_jy:.4f} Jy  {surface.tb_planck_k:.2f} K")
print(f"subsurface {subsurface.flux_density_jy:.4f} Jy  {subsurface.tb_planck_k:.2f} K  "
      f"beta {subsurface.beta_nadir:.2f}")
```

```
surface    3.1934 Jy  104.19 K
subsurface 3.0011 Jy  98.24 K  beta 1.93
```

`observe.explicit_geometry(lat, lon, distance_km, run=run, step=n)` replaces
`observer_geometry` when the geometry is given directly. Each result also has `as_row()`, a
dictionary of every quantity in the table above.

## From the command line

```bash
ptm observe output/quickstart_europa_diurnal \
    --time 2025-01-01T12:00:00 --freq-ghz 233 --emissivity 0.9
ptm observe output/quickstart_europa_observing \
    --time 2025-01-01T12:00:00 --freq-ghz 233 --porosity 0.5
ptm observe output/quickstart_europa_diurnal --step 40 \
    --sub-observer 0 90 6.3e8 --freq-ghz 233 --emissivity 0.9
```

Repeat `--time` for several observations, and add `--csv FILE` to write the table with every
column. `ptm observe --help` lists the options.
