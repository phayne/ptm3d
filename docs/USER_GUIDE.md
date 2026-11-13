# PTM3D User Guide

PTM3D computes surface and subsurface temperatures of airless bodies on three-dimensional
terrain. This guide describes how a simulation is set up and run, what it writes, and how to
read the results. Installation is covered in [../INSTALL.md](../INSTALL.md), a first run in the
[tutorials](tutorials/README.md), every configuration key in
[configuration.md](configuration.md), and the environment variables in
[ENVIRONMENT.md](ENVIRONMENT.md). Run all commands from the repository root.

---

## How a simulation works

A simulation has five stages.

1. **Mesh.** The surface is a mesh of triangular faces. It can be a global icosphere, an
   icosphere displaced by global topography, or a mesh you import (PLY, OBJ or VTK), for
   example one built from a digital elevation model.
2. **Illumination.** At every time step the model locates the Sun with the NASA SPICE toolkit
   and traces shadow rays from each face toward a finite solar disk. For a moon, it adds the
   sunlight reflected by its planet and the planet's thermal emission (earthshine on the Moon,
   Jupiter-shine on Europa).
3. **Spin-up.** Beneath every face is a one-dimensional column of regolith. Each column is
   integrated under the real forcing for a number of diurnal cycles before the output window,
   so that the subsurface reaches a periodic state.
4. **Output window.** The columns are then integrated through the window you asked for, and
   temperatures and fluxes are recorded at each output time.
5. **Radiation exchange (optional).** On concave terrain, faces exchange scattered sunlight and
   thermal infrared radiation. View factors between faces feed a two-band (visible and
   infrared) radiosity calculation, and the resulting indirect flux is returned to the thermal
   model until the temperatures converge.

The physics and numerics of each stage are described in [thermal-model.md](thermal-model.md),
[illumination-model.md](illumination-model.md) and
[radiation-exchange.md](radiation-exchange.md).

---

## Running a simulation

A simulation is described by a YAML configuration file.

```bash
./build/bin/ptm_thermal_sim --config examples/quickstart_moon_diurnal.yaml
```

Options given on the command line take precedence over the file. To check a configuration
without running it, add `--validate-config`.

```bash
./build/bin/ptm_thermal_sim --config my_config.yaml --validate-config
```

Keys and command-line options that the model does not recognize are listed as warnings and
have no effect. Read these warnings, because a misspelled key leaves the default in place.

The examples in [`../examples/`](../examples/) are written to be run from a clone on a CPU. The
following run in about a minute or less on four cores.

| Example | What it simulates |
|---|---|
| `quickstart_moon_diurnal.yaml` | One lunar day on a 1,280-face global Moon |
| `quickstart_europa_diurnal.yaml` | One Europa day on a 5,120-face global Europa, including a Jupiter eclipse |
| `quickstart_europa_albedo_map.yaml` | The same, with Europa's observed albedo map |
| `tutorial04_crater80N_radoff.yaml`, `tutorial04_crater80N_coupled.yaml` | A 12,800-face bowl crater at 80°N without and with radiation exchange |
| `quickstart_asteroid_diurnal.yaml` | A synthetic 5 km asteroid on a circular orbit |
| `quickstart_asteroid_bennu.yaml` | (101955) Bennu with its OSIRIS-REx shape model and SPICE ephemeris |

---

## The configuration file

A configuration file has sections named `simulation`, `mesh`, `thermal`, `thermophysical`,
`illumination`, `radiation`, `output` and, for bodies without built-in properties,
`generic_body`. Every key has a default, and [configuration.md](configuration.md) lists them
all. This section describes the keys most simulations need.

### Time (`simulation`)

```yaml
simulation:
  body: Moon
  start_time: '2026-04-05T00:00:00'
  duration_days: 29.53
  timesteps: 240
```

`start_time` is in UTC. The output times are `start_time + i * duration_days / timesteps` for
`i = 0, 1, ..., timesteps - 1`, so the first output is at `start_time` and the last is one step
before the end of the window. To obtain a temperature map at a particular time, such as the
time of an observation, choose `start_time` and the cadence so that one output time falls on
it. The spin-up is added before `start_time` and does not change these times.

### Surface (`mesh`)

| Key | Meaning |
|---|---|
| `resolution` | Icosphere level when no mesh is imported. Level *n* has 20 × 4<sup>*n*</sup> faces (level 3 has 1,280, level 5 has 20,480). |
| `use_topographic` | Displace the icosphere by global topography. The Moon topography is in the `moon-topo` dataset (see below). |
| `import_file` | A mesh file to simulate instead of an icosphere. |
| `occluder_file` | A mesh that casts shadows but is not itself simulated. Use it to represent distant terrain around a high-resolution site (see [extended-occluder-architecture.md](extended-occluder-architecture.md)). |

The length unit of an imported mesh is inferred from its size relative to the body radius.
Asteroid shape models are usually in kilometers and DEM meshes in meters. The run log reports
the unit it inferred. To make a mesh from a GeoTIFF elevation model, see
[dem-to-mesh-workflow.md](dem-to-mesh-workflow.md).
[Tutorial 5](tutorials/05_lunar_polar_ice.md) goes through the whole process for a lunar polar
region, from the download of a public DEM to the mesh, the distant terrain and the results.

### Heat conduction (`thermal`)

| Key | Meaning |
|---|---|
| `run_thermal_mesh` | Solve for temperatures. Set it to `true`. With `false` only illumination is computed. |
| `equilibration_orbits` | Diurnal cycles of spin-up before the output window (default 1). |
| `threads` | CPU threads, 0 for all cores. |
| `thermal_m`, `thermal_b` | The depth grid. `thermal_m` is the number of layers in the upper diurnal skin depth (default 10), and `thermal_b` the column depth in diurnal skin depths (default 15). |
| `geothermal_heat_flow` | Basal heat flow in W m<sup>-2</sup>. The default is the body's value (0.018 for the Moon, 0.030 for Europa). |
| `compute_statistics` | Write per-face minimum, mean and maximum over the output window. Leave it on. |
| `compute_ice_stability` | Compute the depth to stable water ice (see [Ice stability](#ice-stability)). |
| `use_gpu` | Use the GPU solver. This needs a GPU build. On a CPU build, set it to `false`. |

The spin-up length determines how well the subsurface has settled. After the spin-up, each run
prints an estimate of the number of cycles each face needs to reach
`equilibration_target_K` (0.1 K by default), and a warning if the configured value is smaller.
The surface temperature settles within a few cycles, but the deeper column and the nighttime
minimum take longer, so the warning is most important for minimum and mean temperatures. The
estimate tends to ask for more cycles than a convergence test shows are needed. To test
convergence directly, repeat a run with a longer spin-up and compare the statistics.

### Surface properties (`thermophysical` and the command line)

The Moon, Mercury and Europa have built-in thermophysical properties. These include regolith
density and conductivity that increase with depth, temperature-dependent conductivity and heat
capacity, albedo, emissivity and basal heat flow. On these bodies a property is changed only
through the override keys.

| Key | Meaning |
|---|---|
| `albedo_override` | Bond albedo at normal incidence |
| `emissivity_override` | Infrared emissivity |
| `ks_override`, `kd_override` | Contact conductivity at the surface and at depth (W m<sup>-1</sup> K<sup>-1</sup>) |
| `H_override` | The e-folding depth of density and conductivity (m) |
| `incidence_angle_albedo` | `true` uses the Hapke directional-hemispheric albedo, `false` the Hayne et al. (2017) empirical albedo |
| `conductivity_model` | `hayne2017` (default) or `martinez2021` |

Plain names such as `albedo`, `ks` or `H` under `thermophysical:` are not used on a built-in
body, and the run lists any it finds in a warning. On a body you define yourself, the
properties belong in `generic_body.surface:` (see [Other bodies](#other-bodies)).

Albedo and thermal inertia can also be set on the command line.

```bash
./build/bin/ptm_thermal_sim --config examples/quickstart_moon_diurnal.yaml --albedo 0.10
./build/bin/ptm_thermal_sim --config examples/quickstart_moon_diurnal.yaml --thermal-inertia 200
```

Thermal inertia Γ = (k ρ c<sub>p</sub>)<sup>1/2</sup> is not a model parameter, so a requested
value of Γ has to be reached by changing *k*, ρ or c<sub>p</sub>. PTM3D offers two ways.

- `--thermal-inertia Γ` scales the surface conductivity to
  k<sub>s</sub> = Γ² / (ρ<sub>s</sub> c<sub>p</sub>(273 K)) and the deep conductivity by the
  same factor. The density profile, the temperature dependence of conductivity and heat
  capacity, and the contrast between surface and depth are kept. Γ is then the nominal value
  at the surface and 273 K, and the effective thermal inertia varies with depth and
  temperature.
- `--constant-thermal-inertia Γ` (or the key `constant_thermal_inertia`) makes the column
  uniform and temperature-independent. It sets k<sub>d</sub> = k<sub>s</sub> and
  ρ<sub>d</sub> = ρ<sub>s</sub>, removes the radiative term in the conductivity, and derives
  c<sub>p</sub> = Γ² / (k<sub>s</sub> ρ<sub>s</sub>). Γ then has the same value everywhere in
  the column, which is the quantity most radiometric fits report.

The two give noticeably different diurnal curves for the same Γ, so state which one you used
when you report a fitted value. They cannot be combined in one run.

When a per-face albedo map is loaded, the map takes precedence over `--albedo` and
`albedo_override`. To vary albedo directly on a body with a map, set
`illumination.use_albedo_map: false`.

### Sunlight and albedo maps (`illumination`)

| Key | Meaning |
|---|---|
| `use_albedo_map` | Use a per-face Bond albedo map when one matches the mesh (default `true`). Maps for the Moon and Europa are in the `albedo-maps` dataset. |
| `albedo_map` | A specific albedo map file. |
| `include_primary_reflection` | Include sunlight reflected by the primary body (Earth or Jupiter). The primary's thermal emission is always included. |

### Self-heating (`radiation`)

| Key | Meaning |
|---|---|
| `enabled` | Exchange scattered sunlight and thermal emission between faces. Needed for craters and rough terrain. It has no effect on a convex body. |
| `view_factor_file` | Where the view-factor matrix is saved after it is computed, and loaded from on later runs. |
| `use_hmatrix` | Store the view factors in compressed form, for meshes above roughly 50,000 faces. |
| `band_separated` | Treat visible and infrared radiation as separate bands (default `true`). |
| `max_coupling_iterations`, `coupling_temp_threshold` | The coupling between radiosity and temperature stops when the largest change in any face's temperature falls below the threshold (default 0.5 K), or after the maximum number of iterations (default 16). |
| `beaming`, `rms_slope_deg` | Directional thermal emission from unresolved roughness (Smith, 1967), with the RMS slope in degrees. |

A relative `view_factor_file` path is taken relative to the directory of the configuration
file, not the directory you run from. If the file is not found, the view factors are computed
again, which can take much longer than the simulation itself on a large mesh. Use an absolute
path when the file is somewhere else.

### What is written (`output`)

| Key | Meaning |
|---|---|
| `directory` | The output directory, relative to the directory you run from |
| `locations_file`, `depths_file` | Sites and depths for per-site time series (see [Temperatures at chosen places and times](#temperatures-at-chosen-places-and-times)) |
| `series_output` | Per-site surface time series. A file name is written in the output directory, an absolute path as given |
| `map_depths_m` | Depths in meters at which to write whole-mesh temperature statistics, as a YAML list such as `[0.05, 0.25]` |
| `skip_time_series` | Write the statistics only, not the per-step files |

---

## Output files

A run writes the following to its output directory.

```
output/quickstart_moon_diurnal/
├── mesh.ply              the mesh, in the face order of every .bin file
├── time_steps.txt        the time of each output step
├── body_positions.csv    Sun and primary-body vectors at each output step
├── metadata.txt          face count, step count, fields written
├── provenance.txt        how the run was made (model version, command line, body, mesh)
├── config_used.yaml      a copy of the configuration file that was read
├── statistics/           per-face statistics over the output window
├── temperature/          surface temperature at each output step
├── illumination/         total absorbed flux at each output step
├── direct_solar/         direct sunlight
├── primary_reflected/    sunlight reflected by the primary body
├── primary_ir/           thermal emission from the primary body
├── indirect_solar/       scattered sunlight (radiation exchange only)
└── indirect_ir/          thermal emission from other faces (radiation exchange only)
```

**Binary files.** Every `.bin` file holds a 4-byte little-endian integer *n*, the number of
faces, followed by *n* 4-byte floats, one per face in mesh order. Temperatures are in kelvin and
fluxes in W m<sup>-2</sup>. The tutorial [The output files](tutorials/03_outputs.md) gives a
reader in a few lines of Python. Reading the file without skipping the 4-byte header shifts
every value by one face.

**Statistics.** For each field, `statistics/<field>_min.bin`, `_avg.bin` and `_max.bin` are the
minimum, mean and maximum at each face over the output window, excluding the spin-up.
`statistics/effective_albedo.bin` is the albedo each face absorbed with, averaged over the
direct sunlight it received. It is not defined (NaN) on faces that were never lit.

**Time steps.** `temperature/temperature_000000.bin` and the files that follow are the surface
temperature at each output step, numbered from zero. Step *i* corresponds to line *i* of
`time_steps.txt`, whose first column is ephemeris time (seconds past J2000, TDB) and whose
second is UTC. The flux files of the same step refer to the same time, so they can be compared
face by face. The index has at least six digits.

**Provenance.** `provenance.txt` records which version of the model produced the run and how it
was invoked. Keep it with any results you publish.

To convert the output for ParaView, NetCDF, HDF5 or GeoTIFF, use the scripts in
[Python tools](#python-tools).

---

## Temperatures at chosen places and times

### Whole-mesh maps at an output time

Read the per-step file of the time you need. The `ptm` Python package does this by step number
and by position.

```python
from ptm import Results

run = Results("output/quickstart_moon_diurnal")
utc = run.times().utc                       # UTC of every output step
T = run.field("temperature", step=120)      # surface temperature of every face
face = run.face_at(lat=0.0, lon=90.0)       # face nearest a given point
t_hours, T_face = run.series("temperature", face=face)
```

Latitudes are planetocentric and longitudes are east-positive. Longitudes between 180° and 360°
are accepted and mean the same as their negative equivalents.

### Time series at named sites

List the sites in a comma-separated file with a header line, or one site per line as
`lat lon name` separated by spaces.

```
name,lat,lon
equator_east,0.0,90.0
high_north,85.0,0.0
```

Name the file in the configuration or on the command line.

```yaml
output:
  locations_file: sites.csv
  series_output: sites_series.csv
```

`series_output` is a surface time series with one row per step and, for each site, the
total, direct solar, primary-reflected and primary infrared fluxes and the surface
temperature. A file name is written inside the output directory, and an absolute path is
used as given. If the file cannot be written, the run ends with an error. The series contains every step of the simulation, beginning with the spin-up, so select the
rows whose time is at or after `start_time` to obtain the output window.

### Temperatures with depth

To follow the subsurface, add a file of depths in meters, one per line.

```yaml
output:
  locations_file: sites.csv
  depths_file: depths.txt
```

Each site then gets `locations/<name>_depths.csv`, with the temperature at every requested depth
at every step. These files also begin with the spin-up. The Python reader returns only the
output window.

```python
profile = run.depth_series("equator_east")      # output window only
profile.depths_m, profile.temperature            # depths and temperatures (step, depth)
```

The deepest depth you can request is the bottom of the model column, which depends on the body
and on `thermal_b`. If you ask for a deeper point, the run stops before the simulation and
reports the deepest available depth. Increase `thermal_b` to reach further.

For maps rather than sites, `map_depths_m` writes the minimum, mean and maximum temperature of
every face over the output window at each requested depth, for example
`statistics/temperature_z0p053m_avg.bin`. Each depth is moved to the nearest node of the depth
grid, and the file name carries the depth actually used.

---

## Ice stability

With `thermal.compute_ice_stability: true`, the run computes for every face the depth below
which water ice is stable against sublimation, and writes it to `statistics/ice_depth.bin` in
meters. Ice is stable at depth *z* when the time-averaged vapor flux escaping through the
overlying regolith falls below a threshold J<sub>0</sub> (Schorghofer & Taylor, 2007). The
threshold and the mean free path of the diffusing molecules are set in
`thermal.ice_stability` as `flux_threshold_kgm2s` (default 3 × 10<sup>-14</sup>
kg m<sup>-2</sup> s<sup>-1</sup>, about 1 mm per billion years) and `grain_hop_length_um`
(default 75 μm).

A value of 0 means ice is stable at the surface. A value equal to the bottom of the model
column means ice is not stable anywhere within the column. The calculation uses the vapor
pressure of water ice only.

The ice depth is available on a CPU with the default solver settings. When radiation exchange is on, the
ice depth is computed from the temperatures before the indirect flux is added, so on faces
heated mostly by their surroundings it can show ice shallower than their temperatures support.
The run prints a notice when this applies. For volatiles other than water, compare the maximum
or mean temperature statistics with the stability temperature of the volatile.
[Tutorial 5](tutorials/05_lunar_polar_ice.md) maps the ice depth and the surface cold traps for
water and carbon dioxide at a lunar polar site.

---

## Other bodies

Any body other than the Moon, Mercury and Europa is described in a `generic_body` block.

```yaml
simulation:
  body: Quickstart Asteroid
generic_body:
  name: Quickstart Asteroid
  ephemeris_source: keplerian     # or spice
  shape:
    type: spherical               # or ellipsoidal, custom_mesh
    radius_km: 5.0
  orbit:
    semi_major_axis_au: 2.7
    eccentricity: 0.0
    inclination_deg: 8.0
    lon_ascending_node_deg: 80.0
    arg_perihelion_deg: 50.0
    mean_anomaly_deg: 0.0
    epoch: '2000-01-01T12:00:00'
  spin:
    period_hours: 5.0
    pole_ra_deg: 0.0
    pole_dec_deg: 90.0
    prime_meridian_deg: 0.0
  surface:
    albedo: 0.05
    emissivity: 0.95
    ks: 0.008333333
    kd: 0.008333333
    rhos: 1500.0
    rhod: 1500.0
    H: 0.02
    cp0: 800.0
    Qb: 0.018
    chi: 0.0
```

The complete example is `examples/quickstart_asteroid_diurnal.yaml`. Every property in
`generic_body.surface` is used. The orbit is either propagated from Keplerian elements or, with
`ephemeris_source: spice`, read from a SPICE kernel named in `spice_metakernel`. Use a SPICE
ephemeris when the results must match real dates and distances, for example during a
spacecraft flyby. `examples/quickstart_asteroid_bennu.yaml` is a complete SPICE example with a
shape model.

The shape is a sphere, a triaxial ellipsoid (`shape.semi_axes_km`), or a facet model
(`shape.type: custom_mesh` with `mesh_file`). The Python packages `horizons` and `damit` in
`src/python/` retrieve orbital elements and physical properties from JPL Horizons and shape
models from the DAMIT database, and write a `generic_body` block. The GUI offers the same
search under **Browse HORIZONS**. A DAMIT model's rotational phase at a given date is not
transferred automatically. The conversion to `spin.prime_meridian_deg` is explained in
[configuration.md](configuration.md#generic_body).

The basal heat flow of a generic body (`surface.Qb`) is zero unless you set it, as in the
example above. Without it, a face that receives no sunlight during the spin-up cools to the
2.7 K floor of the model, which can happen in the polar night of a body with high obliquity.

---

## Long simulations

A season or an orbit on a large mesh can exceed the available memory. The script
`scripts/run_chunked_simulation.py` divides such a simulation into consecutive segments, carries
the subsurface state from each segment to the next, and combines the statistics of all
segments with `ptm_merge_chunks` into one set for the whole simulation.

```bash
python3 scripts/run_chunked_simulation.py --config my_config.yaml \
    --total-days 365 --chunks 4 --output-dir output/my_year
```

The script needs the Python environment described in [../INSTALL.md](../INSTALL.md). For an
imported mesh it can choose the number of segments from the available memory, and for an
icosphere you give it with `--chunks`. Add `--dry-run` to see the segments without running
them, and run the script with `--help` for the other options.

---

## Command-line options

Run `./build/bin/ptm_thermal_sim --help` for the complete list. The most used options are below.

**Simulation**

| Option | Meaning |
|---|---|
| `--config FILE` | Configuration file |
| `--validate-config` | Check the configuration and stop |
| `--body NAME` | Body. Moon, Mercury and Europa have built-in properties |
| `--start TIME` | Start time, ISO 8601 UTC |
| `--duration DAYS` | Length of the output window in days |
| `--timesteps N` | Number of output steps |
| `--resolution N` | Icosphere level |
| `--mesh-file FILE` | Mesh to import |
| `--output DIR` | Output directory |
| `--verbose` | Detailed log |

**Thermal model**

| Option | Meaning |
|---|---|
| `--thermal-mesh` | Solve for temperatures |
| `--thermal-threads N` | CPU threads (0 for all cores) |
| `--equilibration-orbits N` | Diurnal cycles of spin-up |
| `--thermal-m N`, `--thermal-b N` | Depth grid (see `thermal_m`, `thermal_b`) |
| `--geothermal-heat-flow Q` | Basal heat flow (W m<sup>-2</sup>) |
| `--albedo A` | Bond albedo at normal incidence |
| `--thermal-inertia Γ` | Thermal inertia by scaling conductivity |
| `--constant-thermal-inertia Γ` | Uniform, temperature-independent thermal inertia |
| `--seed-indirect MODE` | Initial condition of a coupled run: `series` (default), `mean` or `off` |

**Radiation exchange**

| Option | Meaning |
|---|---|
| `--radiation-exchange` | Turn on radiation exchange between faces |
| `--rad-hmatrix` | Compressed view factors, for large meshes |
| `--rad-view-factor-file FILE` | Save or load the view factors |
| `--rad-coupling-iter N` | Maximum coupling iterations (default 16) |
| `--rad-coupling-temp-tol K` | Temperature convergence threshold (default 0.5 K) |
| `--rad-hmatrix-threads N` | Threads for computing the view factors |

**Illumination and output**

| Option | Meaning |
|---|---|
| `--albedo-map FILE` | Albedo map to use |
| `--no-primary-reflection` | Omit sunlight reflected by the primary body |
| `--locations FILE` | Sites for time series |
| `--depths FILE` | Depths for per-site profiles |
| `--series-output FILE` | Per-site surface series (a name inside the output directory, or an absolute path) |
| `--output-format FORMAT` | `binary` (default) or `vtk` |

**GPU build only**

| Option | Meaning |
|---|---|
| `--use-gpu-thermal` | GPU thermal solver |
| `--rad-use-gpu` | GPU radiosity solver |

---

## Python tools

### The `ptm` package

Install the package from the repository root.

```bash
python3 -m pip install -e .              # reading, plotting, sweeps
python3 -m pip install -e '.[regrid]'    # adds scipy, needed for map plots
```

The `ptm` command summarizes, plots, runs and sweeps simulations.

```bash
ptm info  output/quickstart_moon_diurnal
ptm plot  output/quickstart_moon_diurnal map --save map.png
ptm run   examples/quickstart_moon_diurnal.yaml
ptm sweep examples/quickstart_asteroid_diurnal.yaml --vary thermal_inertia=50,100,200,400
```

In Python, `Results` reads a finished run and `Sweep` runs a grid of simulations from one
configuration file.

```python
from ptm import Results, Sweep

run = Results("output/quickstart_moon_diurnal")
tmax = run.stat("temperature", "max")
run.plot_globes("temperature", "max", save="globes.png")

sweep = Sweep("examples/quickstart_asteroid_diurnal.yaml",
              grid={"thermal_inertia": [50, 100, 200, 400]})
sweep.run()
peak = sweep.collect(lambda r: r.stat("temperature", "max").max())
```

A sweep writes each point to its own directory and never edits the base configuration. A grid
name can be `thermal_inertia`, `constant_thermal_inertia`, `albedo`, or any configuration key
written with dots, such as `thermophysical.emissivity_override`. Before trusting a sweep,
confirm that its points differ. Two points with identical output files mean the swept key had
no effect.

### Brightness temperature

<!-- Placeholder: Ada-9 adds the pointer paragraph for `ptm observe` (disk-integrated values for an observer geometry) and owns the numbers in this section. -->

`scripts/emission/depth_weighted_emission.py` converts a subsurface temperature profile into the
brightness temperature a radiometer would measure at a given frequency, using a porous-ice
dielectric model and Fresnel surface emissivity. Use it with the depth profiles described in
[Temperatures with depth](#temperatures-with-depth).
[`examples/python/fit_thermal_inertia.py`](../examples/python/fit_thermal_inertia.py) is a complete
example. It runs a sweep of Europa simulations over thermal inertia, converts each depth profile
to a 233 GHz brightness-temperature curve, and fits the thermal inertia of a synthetic
observation by interpolating between the sweep points.

### Conversion and plotting scripts (`src/python/`)

| Script | Purpose |
|---|---|
| `vtk_converter.py` | Output to VTK/PVD for ParaView (needs `vtk`) |
| `netcdf_converter.py` | NetCDF (needs `vtk` and `netCDF4`) |
| `ptm_hdf5.py` | HDF5 |
| `ptm_geotiff.py` | GeoTIFF |
| `plot_depth_profiles.py` | Temperature against depth |
| `extract_timeseries.py` | Data at chosen latitudes and longitudes |
| `thermal_vis.py` | Temperature maps |

```bash
python src/python/vtk_converter.py output/quickstart_moon_diurnal output/qs_paraview
```

### Data

Science datasets (lunar topography, property and albedo maps, crater meshes, Europa terrain)
are downloaded with `scripts/fetch_data.py`. Run `python3 scripts/fetch_data.py --list` for the
catalog and [data/README.md](../data/README.md) for sources and licenses. Set `PTM_DATA_ROOT`
to keep the data outside the repository.

The configurations in `configs/` simulate five Europa terrain models from Galileo stereo images
at a range of latitudes. They are written for a GPU build and need the `europa-terrain`
dataset. To run one on a CPU, copy it, set `thermal.use_gpu` and `radiation.use_gpu` to
`false`, and give `radiation.view_factor_file` as the absolute path of the downloaded `.vf`
file so the view factors are loaded rather than recomputed.

---

## Graphical interface

The desktop application in `gui/` builds configurations, runs the model and displays the
results in 3D.

```bash
python3 -m pip install -r gui/requirements.txt
cd gui && python3 run_gui.py
```

[Tutorial 4](tutorials/04_gui.md) and [gui/README.md](../gui/README.md) describe it.

---

## Research modules

**Exosphere and ice evolution.** The program `build/bin/ptm_exosphere_main` uses the surface
temperatures of a PTM3D run to compute water-vapor sublimation, ion sputtering yields, the
transition of amorphous to crystalline ice, the growth of a regolith lag, the feedback of surface
composition on albedo, and the evolution of these over long timescales.
`src/python/vtk_converter_exosphere.py` converts its output to VTK for ParaView.

**Visibility.** The program `build/bin/dark_horse_viewshed` computes line-of-sight visibility
over a terrain mesh.
